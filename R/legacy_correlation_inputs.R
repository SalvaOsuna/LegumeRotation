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

#' Prepare Lentil Microbiome Predictors for Legacy Correlations
#'
#' Aggregates plot-level Year-1 microbiome traits to one value per lentil
#' genotype and environment, then shifts the environment labels to the Year-2
#' wheat environments so these predictors can be correlated with wheat legacy
#' phenotypes.
#'
#' @param data Microbiome data frame.
#' @param trait_cols Numeric microbiome trait columns. If `NULL`, a curated set
#'   of alpha-diversity, ordination, and symbiont-abundance traits is used when
#'   present; otherwise numeric non-design columns are used.
#' @param env_col Environment column.
#' @param gen_col Lentil genotype column.
#' @param type_col Optional type column used to exclude checks.
#' @param check_values Values in `type_col` treated as checks.
#' @param include_checks Logical. Include checks?
#' @param from_year Year suffix for lentil microbiome observations.
#' @param to_year Year suffix for following wheat observations.
#' @param filter_from_year Logical. Keep only environments ending in `from_year`
#'   before shifting to `to_year`.
#' @param trait_prefix Prefix added to microbiome trait names.
#' @param source_label Label stored in the Source column.
#' @param drop_all_missing Logical. Drop traits with no observed values.
#' @param drop_constant Logical. Drop traits with one observed value across all
#'   retained samples.
#'
#' @return Long-format data frame with Genotype, Environment, Trait, Predicted.
#' @export
prepare_lentil_microbiome_predictors <- function(data,
                                                 trait_cols = NULL,
                                                 env_col = "ENV",
                                                 gen_col = "Lentil",
                                                 type_col = "Type",
                                                 check_values = "Check",
                                                 include_checks = FALSE,
                                                 from_year = "2024",
                                                 to_year = "2025",
                                                 filter_from_year = TRUE,
                                                 trait_prefix = "Microbiome_",
                                                 source_label = "lentil_microbiome_mean",
                                                 drop_all_missing = TRUE,
                                                 drop_constant = TRUE) {
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
    curated_cols <- c(
      "Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson",
      "Pielou_J", "Faith_PD", "Faith_SR", "Aitchison_PC1",
      "Aitchison_PC2", "MicrobiomeLegacyScore",
      "RhizobiaceaeRelAbundance", "BradyrhizobiaceaeRelAbundance",
      "LegumeSymbiontRelAbundance"
    )
    trait_cols <- curated_cols[curated_cols %in% names(out)]

    if (length(trait_cols) == 0) {
      excluded_cols <- c(
        env_col, gen_col, type_col,
        "SampleID", "Sample", "Sample_raw", "primer", "plot", "Unique ID",
        "Unique.ID", "SiteID", "Col", "Row", "Block", "BlockRow",
        "BlockCol", "Rep", "Wheat", "Name", "Name_SNPset", "Name.SNPset",
        "cot_color", "LibrarySizeRaw", "LibrarySizeFiltered",
        "Aitchison_PC1_var", "Aitchison_PC2_var"
      )
      trait_cols <- names(out)[vapply(out, is.numeric, logical(1))]
      trait_cols <- setdiff(trait_cols, excluded_cols)
    }
  }

  trait_cols <- trait_cols[trait_cols %in% names(out)]
  trait_cols <- trait_cols[vapply(out[trait_cols], is.numeric, logical(1))]

  if (drop_all_missing) {
    trait_cols <- trait_cols[vapply(out[trait_cols], function(x) any(!is.na(x)), logical(1))]
  }
  if (drop_constant) {
    trait_cols <- trait_cols[vapply(out[trait_cols], function(x) {
      length(unique(stats::na.omit(x))) > 1
    }, logical(1))]
  }
  if (length(trait_cols) == 0) stop("No microbiome trait columns were found.")

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
      Trait = paste0(trait_prefix, sub("^Microbiome", "", .data[["Trait"]])),
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
#' @param phenotype_col Legacy-value column to use as the wheat target.
#'
#' @return Long-format data frame with Genotype, Environment, Trait, Predicted.
#' @export
prepare_wheat_legacy_targets <- function(...,
                                         prefix = "Wheat_Legacy_",
                                         phenotype_col = "Legacy_Value") {
  legacy_list <- list(...)
  if (length(legacy_list) == 0) stop("Provide at least one wheat legacy result.")

  target_names <- names(legacy_list)
  if (is.null(target_names)) target_names <- rep("", length(legacy_list))
  unnamed <- target_names == ""
  if (any(unnamed)) target_names[unnamed] <- paste0("Trait", seq_len(sum(unnamed)))

  out <- vector("list", length(legacy_list))
  for (i in seq_along(legacy_list)) {
    x <- legacy_list[[i]]
    if (is.list(x) && "design_audit" %in% names(x) && !is.null(x$design_audit) && !isTRUE(x$design_audit$design_ok)) {
      warning(
        "Design audit for wheat legacy target '", target_names[[i]], "' is not fully OK. ",
        "Inspect the predecessor model design_audit before using correlations.",
        call. = FALSE
      )
    }
    df <- if (is.list(x) &&
              "legacy_values" %in% names(x) &&
              phenotype_col %in% names(x$legacy_values)) {
      x$legacy_values
    } else if (is.list(x) &&
               "gwas_phenotypes" %in% names(x) &&
               phenotype_col %in% names(x$gwas_phenotypes)) {
      x$gwas_phenotypes
    } else if (is.data.frame(x) && phenotype_col %in% names(x)) {
      x
    } else {
      extract_predecessor_gwas_phenotypes(x)
    }
    if (!"Genotype" %in% names(df) && "Previous_Genotype" %in% names(df)) {
      df$Genotype <- df$Previous_Genotype
    }
    .lr_check_columns(df, c("Genotype", "Environment", phenotype_col))
    if (!"N_Plots" %in% names(df)) df$N_Plots <- NA_integer_
    if (!"SE" %in% names(df) && "Corrected_Mean_SE" %in% names(df)) df$SE <- df$Corrected_Mean_SE
    if (!"SE" %in% names(df)) df$SE <- NA_real_
    if (!"Reliability" %in% names(df)) df$Reliability <- NA_real_
    df$SE <- suppressWarnings(as.numeric(df$SE))
    df$Reliability <- suppressWarnings(as.numeric(df$Reliability))
    df$Weight <- dplyr::case_when(
      !is.na(df$Reliability) ~ df$Reliability,
      !is.na(df$SE) & df$SE > 0 ~ 1 / (df$SE^2),
      TRUE ~ NA_real_
    )

    out[[i]] <- df |>
      dplyr::transmute(
        Genotype = as.character(.data[["Genotype"]]),
        Environment = as.character(.data[["Environment"]]),
        Trait = paste0(prefix, target_names[[i]]),
        Predicted = .data[[phenotype_col]],
        Source = "wheat_legacy",
        N_Observations = .data[["N_Plots"]],
        SE = .data[["SE"]],
        Reliability = .data[["Reliability"]],
        Weight = .data[["Weight"]]
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
