#' Shift Lentil Environments to Following Wheat Environments
#'
#' Converts Year 1 lentil environment labels to the matching Year 2 wheat
#' environment labels by replacing a terminal source year with a target year.
#'
#' @param env Character vector of environment names.
#' @param from_year Year suffix in the lentil environment.
#' @param to_year Year suffix in the wheat environment.
#'
#' @return Character vector with shifted environment labels.
#' @export
shift_rotation_environment <- function(env, from_year = "2024", to_year = "2025") {
  sub(paste0(from_year, "$"), to_year, as.character(env))
}

#' Prepare Lentil BLUP Predictors for Legacy Correlations
#'
#' Converts a `model_traits()` BLUP table into the long format used by
#' `plot_legacy_correlations()` and `plot_legacy_driver_correlations()`.
#'
#' @param blups BLUP table with Genotype, Environment, Trait, and Predicted.
#' @param from_year Year suffix for lentil observations.
#' @param to_year Year suffix for following wheat observations.
#' @param trait_prefix Optional prefix added to lentil trait names.
#' @param source_label Label stored in the Source column.
#'
#' @return Long-format data frame with Genotype, Environment, Trait, Predicted.
#' @export
prepare_lentil_blup_predictors <- function(blups,
                                           from_year = "2024",
                                           to_year = "2025",
                                           trait_prefix = "Lentil_",
                                           source_label = "lentil_blup") {
  .lr_check_columns(blups, c("Genotype", "Environment", "Trait", "Predicted"))

  blups |>
    dplyr::filter(grepl(paste0(from_year, "$"), .data[["Environment"]])) |>
    dplyr::transmute(
      Genotype = as.character(.data[["Genotype"]]),
      Environment = shift_rotation_environment(.data[["Environment"]], from_year = from_year, to_year = to_year),
      Trait = paste0(trait_prefix, .data[["Trait"]]),
      Predicted = .data[["Predicted"]],
      Source = source_label,
      N_Observations = NA_integer_
    )
}

#' Prepare Lentil Subsample Predictors for Legacy Correlations
#'
#' Aggregates plot-level lentil subsample traits to one value per lentil
#' genotype and environment, then shifts the environment to the following wheat
#' year so the traits can be correlated with wheat legacy phenotypes.
#'
#' @param data Lentil subsample data frame.
#' @param trait_cols Numeric subsample trait columns. If `NULL`, numeric columns
#'   excluding design and identifier columns are used.
#' @param env_col Environment column.
#' @param gen_col Lentil genotype column.
#' @param type_col Optional type column used to exclude checks.
#' @param check_values Values in `type_col` treated as checks.
#' @param include_checks Logical. Include checks?
#' @param from_year Year suffix for lentil observations.
#' @param to_year Year suffix for following wheat observations.
#' @param filter_from_year Logical. Keep only environments ending in `from_year`
#'   before shifting to `to_year`.
#' @param trait_prefix Prefix added to subsample trait names.
#' @param source_label Label stored in the Source column.
#'
#' @return Long-format data frame with Genotype, Environment, Trait, Predicted.
#' @export
prepare_lentil_subsample_predictors <- function(data,
                                                trait_cols = NULL,
                                                env_col = "ENV",
                                                gen_col = "Lentil",
                                                type_col = "Type",
                                                check_values = "Check",
                                                include_checks = FALSE,
                                                from_year = "2024",
                                                to_year = "2025",
                                                filter_from_year = TRUE,
                                                trait_prefix = "Subsample_",
                                                source_label = "lentil_subsample_mean") {
  .lr_check_columns(data, c(env_col, gen_col))

  out <- data
  if (!include_checks && !is.null(type_col) && type_col %in% names(out)) {
    out <- out[!(out[[type_col]] %in% check_values), ]
  }
  if (filter_from_year) {
    out <- out[grepl(paste0(from_year, "$"), out[[env_col]]), ]
  }
  out <- out[!is.na(out[[gen_col]]) & grepl("^L[0-9]{3}$", out[[gen_col]]), ]

  if (is.null(trait_cols)) {
    excluded_cols <- c(
      env_col, gen_col, type_col,
      "Unique.ID", "Col", "Row", "Block", "BlockRow", "BlockCol",
      "Rep", "Plot", "Wheat", "Combo", "Name", "Name.SNPset"
    )
    trait_cols <- names(out)[vapply(out, is.numeric, logical(1))]
    trait_cols <- setdiff(trait_cols, excluded_cols)
  }
  trait_cols <- trait_cols[trait_cols %in% names(out)]
  if (length(trait_cols) == 0) stop("No subsample trait columns were found.")

  out |>
    dplyr::select(dplyr::all_of(c(env_col, gen_col, trait_cols))) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(trait_cols),
      names_to = "Trait",
      values_to = "Predicted"
    ) |>
    dplyr::filter(!is.na(.data[["Predicted"]])) |>
    dplyr::group_by(.data[[gen_col]], .data[[env_col]], .data[["Trait"]]) |>
    dplyr::summarize(
      Predicted = mean(.data[["Predicted"]], na.rm = TRUE),
      N_Observations = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::transmute(
      Genotype = as.character(.data[[gen_col]]),
      Environment = shift_rotation_environment(.data[[env_col]], from_year = from_year, to_year = to_year),
      Trait = paste0(trait_prefix, .data[["Trait"]]),
      Predicted = .data[["Predicted"]],
      Source = source_label,
      N_Observations = .data[["N_Observations"]]
    )
}

#' Prepare Wheat Legacy Targets for Correlation Analysis
#'
#' Converts one or more `model_predecessor_effect()` outputs into long-format
#' wheat legacy target traits.
#'
#' @param ... Named `model_predecessor_effect()` results. Names are used in
#'   target trait labels, for example `YADJ = avg_predecessor_yadj`.
#' @param prefix Prefix used for target trait labels.
#' @param phenotype_col Phenotype column in each result's `gwas_phenotypes`.
#'
#' @return Long-format data frame with Genotype, Environment, Trait, Predicted.
#' @export
prepare_wheat_legacy_targets <- function(...,
                                         prefix = "Wheat_Legacy_",
                                         phenotype_col = "Predecessor_Phenotype") {
  legacy_list <- list(...)
  if (length(legacy_list) == 0) stop("Provide at least one wheat legacy result.")

  target_names <- names(legacy_list)
  if (is.null(target_names)) target_names <- rep("", length(legacy_list))
  unnamed <- target_names == ""
  if (any(unnamed)) target_names[unnamed] <- paste0("Trait", seq_len(sum(unnamed)))

  out <- vector("list", length(legacy_list))
  for (i in seq_along(legacy_list)) {
    x <- legacy_list[[i]]
    df <- if (is.list(x) && "gwas_phenotypes" %in% names(x)) x$gwas_phenotypes else extract_predecessor_gwas_phenotypes(x)
    .lr_check_columns(df, c("Genotype", "Environment", phenotype_col))
    if (!"N_Plots" %in% names(df)) df$N_Plots <- NA_integer_

    out[[i]] <- df |>
      dplyr::transmute(
        Genotype = as.character(.data[["Genotype"]]),
        Environment = as.character(.data[["Environment"]]),
        Trait = paste0(prefix, target_names[[i]]),
        Predicted = .data[[phenotype_col]],
        Source = "wheat_legacy",
        N_Observations = .data[["N_Plots"]]
      )
  }

  dplyr::bind_rows(out)
}

#' Build Legacy Correlation Input
#'
#' Combines any number of lentil predictor tables with one or more wheat legacy
#' target tables into the long format expected by the legacy correlation plots.
#'
#' @param predictors List of long-format predictor tables.
#' @param targets Long-format wheat legacy target table.
#'
#' @return Combined long-format data frame.
#' @export
build_legacy_correlation_input <- function(predictors, targets) {
  if (!is.list(predictors)) predictors <- list(predictors)
  predictors <- predictors[!vapply(predictors, is.null, logical(1))]

  required_cols <- c("Genotype", "Environment", "Trait", "Predicted")
  lapply(predictors, .lr_check_columns, cols = required_cols)
  .lr_check_columns(targets, required_cols)

  dplyr::bind_rows(c(predictors, list(targets)))
}
