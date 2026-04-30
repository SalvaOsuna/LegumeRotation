#' Derive Wheat Rotation Traits
#'
#' Calculates standardized biomass, N, C, harvest-index, and C:N traits from a
#' wheat subsample table.
#'
#' @param data Wheat subsample data.
#' @param grain_biomass_col Grain biomass column.
#' @param straw_biomass_col Straw biomass column.
#' @param seed_n_pct_col Seed N percentage column.
#' @param straw_n_pct_col Straw N percentage column.
#' @param seed_c_pct_col Seed C percentage column.
#' @param straw_c_pct_col Straw C percentage column.
#' @param prefix Optional prefix added to derived trait names.
#'
#' @return The original data with derived rotation traits appended.
#' @export
derive_wheat_rotation_traits <- function(data,
                                         grain_biomass_col = "Grain_Biomass_g",
                                         straw_biomass_col = "Straw_Biomass_g",
                                         seed_n_pct_col = "Seed_N_pct",
                                         straw_n_pct_col = "Straw_N_pct",
                                         seed_c_pct_col = "Seed_C_pct",
                                         straw_c_pct_col = "Straw_C_pct",
                                         prefix = "") {
  required_cols <- c(
    grain_biomass_col,
    straw_biomass_col,
    seed_n_pct_col,
    straw_n_pct_col,
    seed_c_pct_col,
    straw_c_pct_col
  )
  .lr_check_columns(data, required_cols)

  total_biomass <- paste0(prefix, "Total_Aboveground_Biomass_g")
  harvest_index <- paste0(prefix, "Harvest_Index")
  grain_n <- paste0(prefix, "Grain_N_Content_g")
  straw_n <- paste0(prefix, "Straw_N_Content_g")
  total_n <- paste0(prefix, "Total_Aboveground_N_g")
  n_hi <- paste0(prefix, "N_Harvest_Index")
  grain_c <- paste0(prefix, "Grain_C_Content_g")
  straw_c <- paste0(prefix, "Straw_C_Content_g")
  total_c <- paste0(prefix, "Total_Aboveground_C_g")
  seed_cn <- paste0(prefix, "Seed_C_to_N")
  straw_cn <- paste0(prefix, "Straw_C_to_N")
  total_cn <- paste0(prefix, "Total_C_to_N")

  data |>
    dplyr::mutate(
      "{total_biomass}" := .data[[grain_biomass_col]] + .data[[straw_biomass_col]],
      "{harvest_index}" := dplyr::if_else(
        .data[[total_biomass]] > 0,
        .data[[grain_biomass_col]] / .data[[total_biomass]],
        NA_real_
      ),
      "{grain_n}" := .data[[grain_biomass_col]] * .data[[seed_n_pct_col]] / 100,
      "{straw_n}" := .data[[straw_biomass_col]] * .data[[straw_n_pct_col]] / 100,
      "{total_n}" := .data[[grain_n]] + .data[[straw_n]],
      "{n_hi}" := dplyr::if_else(
        .data[[total_n]] > 0,
        .data[[grain_n]] / .data[[total_n]],
        NA_real_
      ),
      "{grain_c}" := .data[[grain_biomass_col]] * .data[[seed_c_pct_col]] / 100,
      "{straw_c}" := .data[[straw_biomass_col]] * .data[[straw_c_pct_col]] / 100,
      "{total_c}" := .data[[grain_c]] + .data[[straw_c]],
      "{seed_cn}" := dplyr::if_else(
        .data[[seed_n_pct_col]] > 0,
        .data[[seed_c_pct_col]] / .data[[seed_n_pct_col]],
        NA_real_
      ),
      "{straw_cn}" := dplyr::if_else(
        .data[[straw_n_pct_col]] > 0,
        .data[[straw_c_pct_col]] / .data[[straw_n_pct_col]],
        NA_real_
      ),
      "{total_cn}" := dplyr::if_else(
        .data[[total_n]] > 0,
        .data[[total_c]] / .data[[total_n]],
        NA_real_
      )
    )
}

#' Get Default Wheat Rotation Trait Sets
#'
#' @return A named list of default wheat rotation trait groups.
#' @export
get_wheat_rotation_trait_sets <- function() {
  list(
    biomass = c(
      "Grain_Biomass_g",
      "Straw_Biomass_g",
      "Total_Aboveground_Biomass_g",
      "Harvest_Index"
    ),
    nitrogen = c(
      "Seed_N_pct",
      "Straw_N_pct",
      "Grain_N_Content_g",
      "Straw_N_Content_g",
      "Total_Aboveground_N_g",
      "N_Harvest_Index"
    ),
    carbon = c(
      "Seed_C_pct",
      "Straw_C_pct",
      "Grain_C_Content_g",
      "Straw_C_Content_g",
      "Total_Aboveground_C_g"
    ),
    stoichiometry = c(
      "Seed_C_to_N",
      "Straw_C_to_N",
      "Total_C_to_N"
    ),
    rotation_core = c(
      "Y.ADJ",
      "Total_Aboveground_Biomass_g",
      "Grain_N_Content_g",
      "Total_Aboveground_N_g",
      "N_Harvest_Index",
      "Straw_N_Content_g",
      "Straw_C_to_N"
    )
  )
}

#' Prepare Wheat Subsample Rotation Predictors
#'
#' Checks and optionally derives wheat subsample traits before wheat rotational
#' partner, GxG rotational value, and multi-trait index analyses.
#'
#' @param data Wheat subsample data.
#' @param trait_cols Numeric trait columns to verify.
#' @param env_col Environment column.
#' @param facet_col Facet/network column.
#' @param prev_gen_col Previous lentil genotype column.
#' @param curr_gen_col Current wheat genotype column.
#' @param combo_col Lentil-wheat combination column.
#' @param spatial_cols Row and column position columns.
#' @param type_col Optional check flag column.
#' @param check_values Values in `type_col` treated as checks.
#' @param include_checks Include checks?
#' @param derive_traits Derive standard wheat rotation traits?
#' @param require_spatial Require `spatial_cols`?
#' @param ... Passed to `derive_wheat_rotation_traits()`.
#'
#' @return Checked and optionally derived data.
#' @export
prepare_wheat_subsample_rotation_predictors <- function(data,
                                                        trait_cols = NULL,
                                                        env_col = "ENV",
                                                        facet_col = "Facet",
                                                        prev_gen_col = "Lentil",
                                                        curr_gen_col = "Wheat",
                                                        combo_col = "Combo",
                                                        spatial_cols = c("Row", "Col"),
                                                        type_col = "Type",
                                                        check_values = "Check",
                                                        include_checks = FALSE,
                                                        derive_traits = TRUE,
                                                        require_spatial = TRUE,
                                                        ...) {
  required_cols <- c(env_col, facet_col, prev_gen_col, curr_gen_col, combo_col)
  if (require_spatial) required_cols <- c(required_cols, spatial_cols)
  if (!is.null(type_col) && type_col %in% names(data)) required_cols <- c(required_cols, type_col)
  .lr_check_columns(data, required_cols)

  out <- data
  if (!include_checks && !is.null(type_col) && type_col %in% names(out)) {
    out <- out[!(out[[type_col]] %in% check_values), ]
  }
  out <- out[!is.na(out[[prev_gen_col]]) & !is.na(out[[curr_gen_col]]), ]
  if (derive_traits) out <- derive_wheat_rotation_traits(out, ...)

  if (!is.null(trait_cols)) {
    .lr_check_columns(out, trait_cols)
    non_numeric <- trait_cols[!vapply(out[trait_cols], is.numeric, logical(1))]
    if (length(non_numeric) > 0) {
      stop(
        "Wheat rotation trait columns must be numeric: ",
        paste(non_numeric, collapse = ", "),
        call. = FALSE
      )
    }
  }

  out
}

#' Prepare Wheat Rotation Traits
#'
#' Alias for `prepare_wheat_subsample_rotation_predictors()`.
#'
#' @export
prepare_wheat_rotation_traits <- function(...) {
  prepare_wheat_subsample_rotation_predictors(...)
}

#' Model Wheat Rotational Partner Value
#'
#' Estimates which wheat genotypes are better following crops after lentil
#' within each ENV x Facet network.
#'
#' A positive `Wheat_Rotational_Value` means the wheat genotype performed better
#' than the average wheat genotype in the same ENV x Facet network for the
#' selected trait, after accounting for the lentil predecessor and field
#' structure.
#'
#' @param data Wheat rotation or subsample data.
#' @param trait Wheat trait to analyze.
#' @param env_col Environment/site-year column.
#' @param facet_col Facet/network column.
#' @param prev_gen_col Previous lentil genotype column.
#' @param curr_gen_col Current wheat genotype column.
#' @param spatial_cols Row and column position columns.
#' @param fixed_effect_cols Optional fixed design terms, usually `Rep_combo`.
#' @param type_col Optional column used to remove checks.
#' @param check_values Values in `type_col` to remove when `include_checks = FALSE`.
#' @param include_checks Include check plots?
#' @param method `"SpATS"` or `"lme4"`.
#' @param fit_scope `"env_facet"` or `"env_global"`.
#' @param validate_design Run local factorial design diagnostics?
#'
#' @return A list with `rotational_values`, fitted models, diagnostics, and plots.
#' @export
model_wheat_rotational_value <- function(data,
                                         trait,
                                         env_col = "ENV",
                                         facet_col = "Facet",
                                         prev_gen_col = "Lentil",
                                         curr_gen_col = "Wheat",
                                         spatial_cols = c("Row", "Col"),
                                         fixed_effect_cols = c("Rep_combo"),
                                         type_col = "Type",
                                         check_values = "Check",
                                         include_checks = FALSE,
                                         method = c("SpATS", "lme4"),
                                         fit_scope = c("env_facet", "env_global"),
                                         validate_design = TRUE) {
  method <- match.arg(method)
  fit_scope <- match.arg(fit_scope)
  .lr_warn_rotation_rep_terms(fixed_effect_cols)

  if (method == "SpATS" && !requireNamespace("SpATS", quietly = TRUE)) {
    warning("SpATS is not installed; falling back to method = 'lme4'.")
    method <- "lme4"
  }
  if (method == "lme4" && !requireNamespace("lme4", quietly = TRUE)) {
    stop("Package 'lme4' is required for method = 'lme4'.")
  }

  data <- .lr_prepare_rotation_data(
    data = data,
    trait = trait,
    env_col = env_col,
    prev_gen_col = prev_gen_col,
    curr_gen_col = curr_gen_col,
    combo_col = ".__WheatRotationCombo",
    spatial_cols = spatial_cols,
    baseline_col = facet_col,
    type_col = type_col,
    check_values = check_values,
    include_checks = include_checks
  )

  design_audit <- NULL
  if (validate_design && facet_col %in% names(data)) {
    design_audit <- audit_rotation_design(
      data = data,
      env_col = env_col,
      facet_col = facet_col,
      lentil_col = prev_gen_col,
      wheat_col = curr_gen_col,
      combo_col = ".__WheatRotationCombo"
    )
    if (!isTRUE(design_audit$design_ok)) {
      warning(
        "Rotation design audit found one or more incomplete or unexpected ENV x Facet networks. ",
        "Inspect result$design_audit and result$diagnostics before interpreting Wheat_Rotational_Value.",
        call. = FALSE
      )
    }
  }

  results_list <- list()
  model_list <- list()
  diagnostics_list <- list()

  for (env in unique(data[[env_col]])) {
    env_data <- droplevels(data[data[[env_col]] == env & !is.na(data[[trait]]), ])
    if (nrow(env_data) == 0) next
    groups <- if (fit_scope == "env_facet") unique(as.character(env_data$.Baseline_Group)) else "All"

    for (facet in groups) {
      group_data <- if (fit_scope == "env_facet") {
        droplevels(env_data[as.character(env_data$.Baseline_Group) == facet, ])
      } else {
        env_data
      }
      if (nrow(group_data) == 0) next

      group_diag <- .lr_wheat_rotation_group_diagnostics(
        env_data = group_data,
        trait = trait,
        prev_gen_col = prev_gen_col,
        curr_gen_col = curr_gen_col
      )
      group_diag$ENV <- as.character(env)
      group_diag$Facet <- as.character(facet)

      if (dplyr::n_distinct(group_data[[curr_gen_col]]) < 2 ||
          dplyr::n_distinct(group_data[[prev_gen_col]]) < 2) {
        group_diag$Status <- "skipped_insufficient_genotype_levels"
        diagnostics_list[[paste(env, facet, "diagnostic", sep = "::")]] <- group_diag
        next
      }

      fit <- tryCatch({
        if (method == "SpATS") {
          .lr_fit_wheat_rotation_spats(
            env_data = group_data,
            trait = trait,
            prev_gen_col = prev_gen_col,
            curr_gen_col = curr_gen_col,
            fixed_effect_cols = fixed_effect_cols
          )
        } else {
          .lr_fit_wheat_rotation_lme4(
            env_data = group_data,
            trait = trait,
            prev_gen_col = prev_gen_col,
            curr_gen_col = curr_gen_col,
            fixed_effect_cols = fixed_effect_cols
          )
        }
      }, error = function(e) {
        if (method == "SpATS" && requireNamespace("lme4", quietly = TRUE)) {
          warning(paste("SpATS wheat rotational model failed for", env, facet, "and will fall back to lme4:", e$message))
          tryCatch(
            .lr_fit_wheat_rotation_lme4(
              env_data = group_data,
              trait = trait,
              prev_gen_col = prev_gen_col,
              curr_gen_col = curr_gen_col,
              fixed_effect_cols = fixed_effect_cols
            ),
            error = function(e2) {
              warning(paste("Wheat rotational model failed for", env, facet, ":", e2$message))
              NULL
            }
          )
        } else {
          warning(paste("Wheat rotational model failed for", env, facet, ":", e$message))
          NULL
        }
      })

      if (is.null(fit)) {
        group_diag$Status <- "model_failed"
        diagnostics_list[[paste(env, facet, "diagnostic", sep = "::")]] <- group_diag
        next
      }

      rot_df <- fit$rotational_values
      rot_df$ENV <- as.character(env)
      rot_df$Facet <- if (fit_scope == "env_facet") as.character(facet) else rot_df$Facet
      rot_df$Trait <- trait
      rot_df$Method <- method
      rot_df$Fit_Scope <- fit_scope

      group_diag$Status <- .lr_wheat_rotation_status(group_diag)
      group_diag$Model_Status <- "fit"
      group_diag$Residual_Variance <- .lr_model_residual_variance(fit$model)
      group_diag$SE_Min <- .lr_safe_min(rot_df$SE)
      group_diag$SE_Median <- .lr_safe_median(rot_df$SE)
      group_diag$SE_Max <- .lr_safe_max(rot_df$SE)

      key <- paste(env, facet, sep = "::")
      results_list[[key]] <- rot_df
      model_list[[key]] <- fit$model
      diagnostics_list[[paste(env, facet, "diagnostic", sep = "::")]] <- group_diag
    }
  }

  final_df <- dplyr::bind_rows(results_list)
  diagnostics_df <- dplyr::bind_rows(diagnostics_list)

  if (nrow(final_df) > 0) {
    final_df <- final_df |>
      dplyr::group_by(.data[["ENV"]], .data[["Trait"]]) |>
      dplyr::mutate(
        .Global_Baseline_Mean = .lr_precision_weighted_mean(.data[["Corrected_Wheat_Mean"]], .data[["SE"]]),
        Wheat_Rotational_Value_Global =
          .data[["Corrected_Wheat_Mean"]] - .data[[".Global_Baseline_Mean"]],
        Confidence_Class = .lr_confidence_class(.data[["SE"]])
      ) |>
      dplyr::select(-dplyr::all_of(".Global_Baseline_Mean")) |>
      dplyr::ungroup() |>
      dplyr::arrange(.data[["ENV"]], .data[["Facet"]], dplyr::desc(.data[["Wheat_Rotational_Value"]]))
    final_df$Rank <- ave(
      -final_df$Wheat_Rotational_Value,
      final_df$ENV,
      final_df$Facet,
      final_df$Trait,
      FUN = function(x) rank(x, ties.method = "first")
    )
  }

  ranked_plot <- plot_wheat_rotational_ranked(final_df, trait = trait)
  raw_vs_corrected_plot <- plot_wheat_rotational_raw_vs_corrected(final_df, trait = trait)

  list(
    rotational_values = final_df,
    models = model_list,
    diagnostics = diagnostics_df,
    design_audit = design_audit,
    ranked_plot = ranked_plot,
    raw_vs_corrected_plot = raw_vs_corrected_plot,
    settings = list(
      trait = trait,
      method = method,
      fit_scope = fit_scope,
      include_checks = include_checks,
      question = "Which wheat genotypes are better following crops after lentil within their local ENV x Facet network?"
    )
  )
}

#' Model GxG Rotational Value
#'
#' Estimates pair-specific rotational interaction values for wheat subsample
#' traits. This is a rotation-specific wrapper around
#' `model_pair_compatibility()` with clearer output naming.
#'
#' A positive `GxG_Rotational_Value` means that the observed lentil-wheat pair
#' performed better than expected from the average lentil legacy effect and the
#' average wheat rotational partner effect within the same local 10 x 10
#' network. Values should only be interpreted for observed pairs.
#'
#' `GxG_Rotational_Pct` is computed through `model_pair_compatibility()` using
#' the same denominator guard as `Compatibility_Pct`: if
#' `abs(Expected_Additive_Mean) < pct_min_denominator`, the percentage is set to
#' `NA`. Use `GxG_Rotational_Value` as the primary interpretation for centered,
#' ordinal, near-zero, or negative-denominator traits.
#'
#' @return A list with `gxg_values`, models, diagnostics, and plots.
#' @export
model_gxg_rotational_value <- function(data,
                                       trait,
                                       env_col = "ENV",
                                       facet_col = "Facet",
                                       prev_gen_col = "Lentil",
                                       curr_gen_col = "Wheat",
                                       combo_col = "Combo",
                                       spatial_cols = c("Row", "Col"),
                                       fixed_effect_cols = c("Rep_combo"),
                                       random_effect_cols = c("Block"),
                                       type_col = "Type",
                                       check_values = "Check",
                                       include_checks = FALSE,
                                       compute_pct = TRUE,
                                       pct_min_denominator = 1e-6,
                                       validate_design = TRUE) {
  result <- model_pair_compatibility(
    data = data,
    trait = trait,
    env_col = env_col,
    prev_gen_col = prev_gen_col,
    curr_gen_col = curr_gen_col,
    combo_col = combo_col,
    spatial_cols = spatial_cols,
    baseline_col = facet_col,
    fixed_effect_cols = fixed_effect_cols,
    random_effect_cols = random_effect_cols,
    type_col = type_col,
    check_values = check_values,
    include_checks = include_checks,
    compute_compatibility_pct = compute_pct,
    pct_min_denominator = pct_min_denominator,
    validate_design = validate_design
  )

  gxg_values <- result$compatibility
  if (nrow(gxg_values) > 0) {
    diag_status <- result$diagnostics
    if (nrow(diag_status) > 0 && all(c("Environment", "Baseline_Group", "Status") %in% names(diag_status))) {
      diag_status <- diag_status |>
        dplyr::transmute(
          Environment = .data[["Environment"]],
          Baseline_Group = .data[["Baseline_Group"]],
          Status = gsub("zero_combo_variance", "zero_gxg_variance", .data[["Status"]], fixed = TRUE)
        )
    } else {
      diag_status <- NULL
    }

    gxg_values <- gxg_values |>
      dplyr::mutate(
        ENV = .data[["Environment"]],
        Facet = .data[["Baseline_Group"]],
        GxG_Rotational_Value = .data[["Compatibility_Value"]],
        GxG_Rotational_Pct = .data[["Compatibility_Pct"]],
        N_Reps_Per_Combo = .data[["N_Plots"]],
        Reliability = dplyr::if_else(
          is.na(.data[["SE"]]) | .data[["SE"]] <= 0,
          NA_real_,
          1 / (1 + .data[["SE"]])
        )
      )
    if (!is.null(diag_status)) {
      gxg_values <- gxg_values |>
        dplyr::left_join(diag_status, by = c("Environment", "Baseline_Group"))
    } else {
      gxg_values$Status <- NA_character_
    }
    gxg_values <- gxg_values |>
      dplyr::select(
        dplyr::all_of(c(
          "ENV",
          "Facet",
          "Combo",
          "Previous_Genotype",
          "Current_Genotype",
          "Trait",
          "Raw_Mean",
          "Expected_Additive_Mean",
          "Corrected_Pair_Mean",
          "GxG_Rotational_Value",
          "GxG_Rotational_Pct",
          "SE",
          "Reliability",
          "N_Plots",
          "N_Reps_Per_Combo"
        )),
        dplyr::any_of("Status"),
        dplyr::everything()
      )
  }

  heatmap <- plot_gxg_rotational_heatmap(gxg_values, trait = trait)
  ranked_plot <- plot_gxg_rotational_ranked(gxg_values, trait = trait)
  diagnostics <- result$diagnostics
  if (nrow(diagnostics) > 0 && "Status" %in% names(diagnostics)) {
    diagnostics$Status <- gsub("zero_combo_variance", "zero_gxg_variance", diagnostics$Status, fixed = TRUE)
  }

  list(
    gxg_values = gxg_values,
    models = result$models,
    diagnostics = diagnostics,
    design_audit = result$design_audit,
    heatmap = heatmap,
    ranked_plot = ranked_plot,
    pair_compatibility = result,
    settings = c(
      result$settings,
      list(
        question = "Which observed lentil-wheat combinations are unusually good or bad for rotation-relevant outcomes?",
        warning = "GxG_Rotational_Value applies only to observed combinations within each local ENV x Facet network."
      )
    )
  )
}

#' Plot Ranked Wheat Rotational Values
#'
#' @param x A `model_wheat_rotational_value()` result or rotational value data frame.
#' @param trait Optional trait label.
#' @param n Optional number of strongest absolute values per ENV x Facet.
#' @param show_se Draw SE bars?
#'
#' @return A ggplot object.
#' @export
plot_wheat_rotational_ranked <- function(x, trait = NULL, n = NULL, show_se = FALSE) {
  df <- .lr_get_wheat_rotation_df(x)
  if (nrow(df) == 0) return(NULL)
  .lr_check_columns(df, c("ENV", "Facet", "Current_Genotype", "Wheat_Rotational_Value"))
  if (is.null(trait) && "Trait" %in% names(df)) trait <- paste(unique(as.character(df$Trait)), collapse = ", ")

  plot_df <- df
  if (!is.null(n) && is.finite(n)) {
    plot_df <- plot_df |>
      dplyr::group_by(.data[["ENV"]], .data[["Facet"]]) |>
      dplyr::slice_max(abs(.data[["Wheat_Rotational_Value"]]), n = n, with_ties = FALSE) |>
      dplyr::ungroup()
  }

  plot_df <- plot_df |>
    dplyr::arrange(.data[["ENV"]], .data[["Facet"]], .data[["Wheat_Rotational_Value"]]) |>
    dplyr::mutate(
      Panel = paste(.data[["ENV"]], .data[["Facet"]], sep = " | "),
      Plot_Genotype = paste(.data[["Panel"]], .data[["Current_Genotype"]], sep = "___"),
      Plot_Genotype = factor(.data[["Plot_Genotype"]], levels = unique(.data[["Plot_Genotype"]]))
    )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(
    x = .data[["Plot_Genotype"]],
    y = .data[["Wheat_Rotational_Value"]]
  )) +
    ggplot2::geom_segment(
      ggplot2::aes(xend = .data[["Plot_Genotype"]], y = 0, yend = .data[["Wheat_Rotational_Value"]]),
      color = "gray72",
      linewidth = 0.4
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = .data[["Facet"]]),
      size = 2.2,
      alpha = 0.9
    )

  if (show_se && "SE" %in% names(plot_df)) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = .data[["Wheat_Rotational_Value"]] - .data[["SE"]],
        ymax = .data[["Wheat_Rotational_Value"]] + .data[["SE"]]
      ),
      width = 0.16,
      color = "gray40",
      alpha = 0.55,
      na.rm = TRUE
    )
  }

  p +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~Panel, scales = "free_y", drop = TRUE) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::scale_x_discrete(labels = function(x) sub("^.*___", "", x)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Wheat Rotational Partner Value",
      subtitle = paste("Within ENV x Facet deviation", if (!is.null(trait)) paste("| trait:", trait) else ""),
      caption = "Values are centered within ENV x Facet. Cross-facet rankings compare within-facet deviations, not absolute wheat superiority across a complete 100 x 100 design.",
      x = "Current wheat genotype",
      y = "Wheat_Rotational_Value",
      color = "Facet"
    )
}

#' Plot Raw vs Corrected Wheat Rotational Means
#'
#' @param x A `model_wheat_rotational_value()` result or rotational value data frame.
#' @param trait Optional trait label.
#'
#' @return A ggplot object.
#' @export
plot_wheat_rotational_raw_vs_corrected <- function(x, trait = NULL) {
  df <- .lr_get_wheat_rotation_df(x)
  if (nrow(df) == 0) return(NULL)
  needed <- c("ENV", "Facet", "Current_Genotype", "Raw_Mean", "Corrected_Wheat_Mean", "Baseline_Mean")
  if (!all(needed %in% names(df))) return(NULL)
  if (is.null(trait) && "Trait" %in% names(df)) trait <- paste(unique(as.character(df$Trait)), collapse = ", ")

  plot_df <- df |>
    dplyr::mutate(Panel = paste(.data[["ENV"]], .data[["Facet"]], sep = " | "))

  ggplot2::ggplot(plot_df, ggplot2::aes(
    y = stats::reorder(.data[["Current_Genotype"]], .data[["Wheat_Rotational_Value"]])
  )) +
    ggplot2::geom_segment(
      ggplot2::aes(x = .data[["Raw_Mean"]], xend = .data[["Corrected_Wheat_Mean"]], yend = .data[["Current_Genotype"]]),
      color = "gray60"
    ) +
    ggplot2::geom_vline(
      ggplot2::aes(xintercept = .data[["Baseline_Mean"]]),
      linetype = "dashed",
      color = "gray35",
      alpha = 0.7
    ) +
    ggplot2::geom_point(ggplot2::aes(x = .data[["Raw_Mean"]]), shape = 1, size = 1.8, color = "black") +
    ggplot2::geom_point(
      ggplot2::aes(x = .data[["Corrected_Wheat_Mean"]], color = .data[["Wheat_Rotational_Value"]] > 0),
      size = 2.4
    ) +
    ggplot2::facet_wrap(~Panel, scales = "free_y", drop = TRUE) +
    ggplot2::scale_color_manual(values = c("TRUE" = "forestgreen", "FALSE" = "firebrick"), guide = "none") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Raw vs Model-Corrected Wheat Rotational Value",
      subtitle = paste("Hollow = raw mean, solid = corrected mean", if (!is.null(trait)) paste("| trait:", trait) else ""),
      x = "Wheat trait mean",
      y = "Current wheat genotype"
    )
}

#' Plot Wheat Rotation Trait Heatmap
#'
#' @param x A rotational value data frame with multiple traits.
#' @param value_col Value column to plot.
#' @param standardize Standardize values within ENV x Facet x Trait?
#'
#' @return A ggplot object.
#' @export
plot_wheat_rotation_trait_heatmap <- function(x,
                                              value_col = "Wheat_Rotational_Value",
                                              standardize = TRUE) {
  df <- .lr_get_wheat_rotation_df(x)
  if (nrow(df) == 0) return(NULL)
  .lr_check_columns(df, c("ENV", "Facet", "Current_Genotype", "Trait", value_col))

  plot_df <- df
  fill_col <- value_col
  if (standardize) {
    fill_col <- ".Standardized_Value"
    plot_df <- plot_df |>
      dplyr::group_by(.data[["ENV"]], .data[["Facet"]], .data[["Trait"]]) |>
      dplyr::mutate("{fill_col}" := .lr_zscore(.data[[value_col]])) |>
      dplyr::ungroup()
  }

  plot_df <- plot_df |>
    dplyr::mutate(Panel = paste(.data[["ENV"]], .data[["Facet"]], sep = " | "))

  ggplot2::ggplot(plot_df, ggplot2::aes(
    x = .data[["Trait"]],
    y = .data[["Current_Genotype"]],
    fill = .data[[fill_col]]
  )) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::facet_wrap(~Panel, scales = "free_y", drop = TRUE) +
    ggplot2::scale_fill_gradient2(low = "firebrick", mid = "white", high = "forestgreen", midpoint = 0) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Wheat Rotational Value Across Traits",
      x = "Trait",
      y = "Current wheat genotype",
      fill = if (standardize) "Standardized value" else value_col
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' Plot GxG Rotational Heatmap
#'
#' @param x A `model_gxg_rotational_value()` result or GxG value data frame.
#' @param trait Optional trait label.
#' @param value_col Value column to plot.
#' @param signal_tol Threshold for dropping no-signal panels.
#' @param drop_no_signal Drop panels with no detectable pair signal?
#'
#' @return A ggplot object.
#' @export
plot_gxg_rotational_heatmap <- function(x,
                                        trait = NULL,
                                        value_col = "GxG_Rotational_Value",
                                        signal_tol = 1e-8,
                                        drop_no_signal = TRUE) {
  df <- .lr_get_gxg_df(x)
  if (nrow(df) == 0) return(NULL)
  .lr_check_columns(df, c("Previous_Genotype", "Current_Genotype", "ENV", "Facet", value_col))
  if (is.null(trait) && "Trait" %in% names(df)) trait <- paste(unique(as.character(df$Trait)), collapse = ", ")

  plot_df <- df |>
    dplyr::group_by(.data[["ENV"]], .data[["Facet"]]) |>
    dplyr::mutate(
      Max_Abs_GxG = max(abs(.data[[value_col]]), na.rm = TRUE),
      Pair_Signal = .data[["Max_Abs_GxG"]] > signal_tol
    ) |>
    dplyr::ungroup()

  if (drop_no_signal) plot_df <- dplyr::filter(plot_df, .data[["Pair_Signal"]])
  if (nrow(plot_df) == 0) return(NULL)

  plot_df <- plot_df |>
    dplyr::mutate(Panel = paste(.data[["ENV"]], .data[["Facet"]], sep = " | "))

  ggplot2::ggplot(plot_df, ggplot2::aes(
    x = .data[["Current_Genotype"]],
    y = .data[["Previous_Genotype"]],
    fill = .data[[value_col]]
  )) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::facet_wrap(~Panel, scales = "free", drop = TRUE) +
    ggplot2::scale_fill_gradient2(low = "firebrick", mid = "white", high = "forestgreen", midpoint = 0) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "GxG Rotational Value",
      subtitle = paste("Observed lentil-wheat combinations only", if (!is.null(trait)) paste("| trait:", trait) else ""),
      caption = "GxG_Rotational_Value is pair-specific and should not be extrapolated to unobserved combinations outside the local ENV x Facet network.",
      x = "Current wheat genotype",
      y = "Previous lentil genotype",
      fill = "GxG value"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
      axis.text.y = ggplot2::element_text(size = 6),
      panel.grid = ggplot2::element_blank()
    )
}

#' Plot Ranked GxG Rotational Values
#'
#' @param x A `model_gxg_rotational_value()` result or GxG value data frame.
#' @param trait Optional trait label.
#' @param top_n Number of strongest absolute effects per ENV x Facet.
#' @param drop_no_signal Drop panels with no detectable pair signal?
#' @param signal_tol Threshold for dropping no-signal panels.
#' @param show_se Draw SE bars?
#'
#' @return A ggplot object.
#' @export
plot_gxg_rotational_ranked <- function(x,
                                       trait = NULL,
                                       top_n = 20,
                                       drop_no_signal = TRUE,
                                       signal_tol = 1e-8,
                                       show_se = FALSE) {
  df <- .lr_get_gxg_df(x)
  if (nrow(df) == 0) return(NULL)
  .lr_check_columns(df, c("Combo", "ENV", "Facet", "GxG_Rotational_Value"))
  if (is.null(trait) && "Trait" %in% names(df)) trait <- paste(unique(as.character(df$Trait)), collapse = ", ")

  group_signal <- df |>
    dplyr::group_by(.data[["ENV"]], .data[["Facet"]]) |>
    dplyr::summarize(
      Max_Abs_GxG = max(abs(.data[["GxG_Rotational_Value"]]), na.rm = TRUE),
      Pair_Signal = .data[["Max_Abs_GxG"]] > signal_tol,
      .groups = "drop"
    )

  plot_df <- df |>
    dplyr::inner_join(group_signal, by = c("ENV", "Facet")) |>
    dplyr::filter(!drop_no_signal | .data[["Pair_Signal"]]) |>
    dplyr::group_by(.data[["ENV"]], .data[["Facet"]]) |>
    dplyr::slice_max(abs(.data[["GxG_Rotational_Value"]]), n = top_n, with_ties = FALSE) |>
    dplyr::ungroup()
  if (nrow(plot_df) == 0) return(NULL)

  plot_df <- plot_df |>
    dplyr::arrange(.data[["ENV"]], .data[["Facet"]], .data[["GxG_Rotational_Value"]]) |>
    dplyr::mutate(
      Panel = paste(.data[["ENV"]], .data[["Facet"]], sep = " | "),
      Plot_Combo = paste(.data[["Panel"]], .data[["Combo"]], sep = "___"),
      Plot_Combo = factor(.data[["Plot_Combo"]], levels = unique(.data[["Plot_Combo"]]))
    )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(
    x = .data[["Plot_Combo"]],
    y = .data[["GxG_Rotational_Value"]]
  )) +
    ggplot2::geom_segment(
      ggplot2::aes(xend = .data[["Plot_Combo"]], y = 0, yend = .data[["GxG_Rotational_Value"]]),
      color = "gray72",
      linewidth = 0.4
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = .data[["GxG_Rotational_Value"]] > 0),
      size = 2.1,
      alpha = 0.9
    )

  if (show_se && "SE" %in% names(plot_df)) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = .data[["GxG_Rotational_Value"]] - .data[["SE"]],
        ymax = .data[["GxG_Rotational_Value"]] + .data[["SE"]]
      ),
      width = 0.16,
      color = "gray40",
      alpha = 0.5,
      na.rm = TRUE
    )
  }

  p +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~Panel, scales = "free_y", drop = TRUE) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::scale_x_discrete(labels = function(x) sub("^.*___", "", x)) +
    ggplot2::scale_color_manual(values = c("TRUE" = "forestgreen", "FALSE" = "firebrick"), guide = "none") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Ranked GxG Rotational Values",
      subtitle = paste("Strongest observed pair deviations", if (!is.null(trait)) paste("| trait:", trait) else ""),
      caption = "Rankings are within observed ENV x Facet networks and should not be extrapolated to unobserved pairs.",
      x = "Observed lentil-wheat pair",
      y = "GxG_Rotational_Value"
    )
}

#' Plot Multi-Trait GxG Rotation Heatmap
#'
#' @param x GxG value data frame with multiple traits.
#' @param value_col Value column to plot.
#' @param standardize Standardize within ENV x Facet x Trait?
#'
#' @return A ggplot object.
#' @export
plot_gxg_rotation_trait_heatmap <- function(x,
                                            value_col = "GxG_Rotational_Value",
                                            standardize = TRUE) {
  df <- .lr_get_gxg_df(x)
  if (nrow(df) == 0) return(NULL)
  .lr_check_columns(df, c("ENV", "Facet", "Combo", "Trait", value_col))

  plot_df <- df
  fill_col <- value_col
  if (standardize) {
    fill_col <- ".Standardized_Value"
    plot_df <- plot_df |>
      dplyr::group_by(.data[["ENV"]], .data[["Facet"]], .data[["Trait"]]) |>
      dplyr::mutate("{fill_col}" := .lr_zscore(.data[[value_col]])) |>
      dplyr::ungroup()
  }

  plot_df <- plot_df |>
    dplyr::mutate(Panel = paste(.data[["ENV"]], .data[["Facet"]], sep = " | "))

  ggplot2::ggplot(plot_df, ggplot2::aes(
    x = .data[["Trait"]],
    y = .data[["Combo"]],
    fill = .data[[fill_col]]
  )) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::facet_wrap(~Panel, scales = "free_y", drop = TRUE) +
    ggplot2::scale_fill_gradient2(low = "firebrick", mid = "white", high = "forestgreen", midpoint = 0) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "GxG Rotational Value Across Traits",
      x = "Trait",
      y = "Observed lentil-wheat pair",
      fill = if (standardize) "Standardized value" else value_col
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' Build Wheat Rotation Index
#'
#' Builds a transparent, user-defined decision-support score from selected
#' wheat rotational value traits. This index is not an intrinsic biological
#' property; different breeding goals may require different trait weights.
#'
#' @return A data frame with one wheat rotation index per group.
#' @export
build_wheat_rotation_index <- function(rotational_values,
                                       trait_weights = c(
                                         Y.ADJ = 1,
                                         Total_Aboveground_N_g = 1,
                                         Grain_N_Content_g = 1,
                                         N_Harvest_Index = 1
                                       ),
                                       value_col = "Wheat_Rotational_Value",
                                       group_cols = c("ENV", "Facet", "Current_Genotype"),
                                       standardize_within = c("ENV", "Facet", "Trait")) {
  df <- .lr_get_wheat_rotation_df(rotational_values)
  .lr_check_columns(df, c(group_cols, standardize_within, "Trait", value_col))
  traits <- names(trait_weights)
  if (is.null(traits) || any(!nzchar(traits))) stop("trait_weights must be a named numeric vector.", call. = FALSE)

  df |>
    dplyr::filter(.data[["Trait"]] %in% traits) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(standardize_within))) |>
    dplyr::mutate(.Standardized_Value = .lr_zscore(.data[[value_col]])) |>
    dplyr::ungroup() |>
    dplyr::mutate(.Trait_Weight = unname(trait_weights[as.character(.data[["Trait"]])])) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarize(
      Wheat_Rotation_Index = sum(.data[[".Standardized_Value"]] * .data[[".Trait_Weight"]], na.rm = TRUE),
      N_Traits_Used = dplyr::n_distinct(.data[["Trait"]]),
      Trait_Weights = paste(
        paste(.data[["Trait"]], .data[[".Trait_Weight"]], sep = "="),
        collapse = ";"
      ),
      Index_Interpretation = "User-defined decision-support score; change trait weights for different breeding goals.",
      .groups = "drop"
    ) |>
    dplyr::arrange(.data[["ENV"]], .data[["Facet"]], dplyr::desc(.data[["Wheat_Rotation_Index"]]))
}

#' Build GxG Rotation Index
#'
#' Ranks observed lentil-wheat pairs across several rotation-relevant traits.
#' The index applies only to observed combinations within each local ENV x Facet
#' network.
#'
#' @return A data frame with one GxG rotation index per observed pair.
#' @export
build_gxg_rotation_index <- function(gxg_values,
                                     trait_weights = c(
                                       Y.ADJ = 1,
                                       Total_Aboveground_Biomass_g = 1,
                                       Total_Aboveground_N_g = 1,
                                       Grain_N_Content_g = 1,
                                       Straw_N_Content_g = 1,
                                       N_Harvest_Index = 1
                                     ),
                                     value_col = "GxG_Rotational_Value",
                                     group_cols = c("ENV", "Facet", "Combo", "Previous_Genotype", "Current_Genotype"),
                                     standardize_within = c("ENV", "Facet", "Trait")) {
  df <- .lr_get_gxg_df(gxg_values)
  .lr_check_columns(df, c(group_cols, standardize_within, "Trait", value_col))
  traits <- names(trait_weights)
  if (is.null(traits) || any(!nzchar(traits))) stop("trait_weights must be a named numeric vector.", call. = FALSE)

  df |>
    dplyr::filter(.data[["Trait"]] %in% traits) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(standardize_within))) |>
    dplyr::mutate(.Standardized_Value = .lr_zscore(.data[[value_col]])) |>
    dplyr::ungroup() |>
    dplyr::mutate(.Trait_Weight = unname(trait_weights[as.character(.data[["Trait"]])])) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarize(
      GxG_Rotation_Index = sum(.data[[".Standardized_Value"]] * .data[[".Trait_Weight"]], na.rm = TRUE),
      N_Traits_Used = dplyr::n_distinct(.data[["Trait"]]),
      Trait_Weights = paste(
        paste(.data[["Trait"]], .data[[".Trait_Weight"]], sep = "="),
        collapse = ";"
      ),
      Index_Interpretation = "Observed-pair decision-support score within ENV x Facet only.",
      .groups = "drop"
    ) |>
    dplyr::arrange(.data[["ENV"]], .data[["Facet"]], dplyr::desc(.data[["GxG_Rotation_Index"]]))
}

.lr_get_wheat_rotation_df <- function(x) {
  if (is.data.frame(x)) return(x)
  if (is.list(x) && "rotational_values" %in% names(x)) return(x$rotational_values)
  stop("Expected a wheat rotational-value data frame or a result from model_wheat_rotational_value().")
}

.lr_get_gxg_df <- function(x) {
  if (is.data.frame(x)) return(x)
  if (is.list(x) && "gxg_values" %in% names(x)) return(x$gxg_values)
  stop("Expected a GxG rotational-value data frame or a result from model_gxg_rotational_value().")
}

.lr_zscore <- function(x) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  s <- stats::sd(x, na.rm = TRUE)
  if (is.na(s) || s <= 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

.lr_wheat_rotation_group_diagnostics <- function(env_data,
                                                 trait,
                                                 prev_gen_col,
                                                 curr_gen_col) {
  n_previous <- dplyr::n_distinct(env_data[[prev_gen_col]])
  n_current <- dplyr::n_distinct(env_data[[curr_gen_col]])
  n_combos <- dplyr::n_distinct(paste(env_data[[prev_gen_col]], env_data[[curr_gen_col]], sep = "-"))
  expected_combos <- n_previous * n_current
  reps <- table(paste(env_data[[prev_gen_col]], env_data[[curr_gen_col]], sep = "-"))
  data.frame(
    Trait = trait,
    N_Plots = nrow(env_data),
    N_Previous_Genotypes = n_previous,
    N_Current_Genotypes = n_current,
    N_Combos = n_combos,
    N_Observed_Combos = n_combos,
    Expected_Combos = expected_combos,
    Missing_Combos = expected_combos - n_combos,
    Local_Factorial_Complete = n_combos == expected_combos,
    Expected_10x10 = n_previous == 10 && n_current == 10 && n_combos == 100,
    Min_Reps_Per_Combo = if (length(reps) > 0) min(as.integer(reps), na.rm = TRUE) else NA_integer_,
    Status = "not_fit",
    Model_Status = NA_character_,
    Residual_Variance = NA_real_,
    SE_Min = NA_real_,
    SE_Median = NA_real_,
    SE_Max = NA_real_
  )
}

.lr_wheat_rotation_status <- function(diag) {
  status <- character(0)
  if (!isTRUE(diag$Local_Factorial_Complete)) status <- c(status, "incomplete_local_factorial")
  if (isTRUE(diag$Expected_10x10)) status <- c(status, "complete_10x10")
  if (diag$N_Previous_Genotypes < 2) status <- c(status, "missing_lentils")
  if (diag$N_Current_Genotypes < 2) status <- c(status, "missing_wheats")
  if (!is.na(diag$Missing_Combos) && diag$Missing_Combos > 0) status <- c(status, "missing_combos")
  if (!is.na(diag$Min_Reps_Per_Combo) && diag$Min_Reps_Per_Combo < 2) status <- c(status, "low_replication")
  if (length(status) == 0) status <- "fit"
  paste(unique(status), collapse = ";")
}

.lr_fit_wheat_rotation_lme4 <- function(env_data,
                                        trait,
                                        prev_gen_col,
                                        curr_gen_col,
                                        fixed_effect_cols) {
  fixed_effect_cols <- .lr_existing_cols(env_data, fixed_effect_cols)
  for (col in fixed_effect_cols) env_data[[col]] <- as.factor(env_data[[col]])

  spatial_terms <- "poly(Row_Scaled, 2, raw = TRUE) + poly(Col_Scaled, 2, raw = TRUE) + Row_Scaled:Col_Scaled"
  fixed_terms <- .lr_formula_terms(c(curr_gen_col, fixed_effect_cols, spatial_terms))
  formula_str <- paste(trait, "~", fixed_terms, "+ (1|", prev_gen_col, ")")
  m <- lme4::lmer(as.formula(formula_str), data = env_data)

  env_data$Corrected <- predict(m, re.form = NA)
  .lr_summarize_wheat_rotation_values(
    env_data = env_data,
    trait = trait,
    curr_gen_col = curr_gen_col,
    prev_gen_col = prev_gen_col,
    corrected_col = "Corrected",
    se_df = NULL,
    model = m
  )
}

.lr_fit_wheat_rotation_spats <- function(env_data,
                                         trait,
                                         prev_gen_col,
                                         curr_gen_col,
                                         fixed_effect_cols) {
  fixed_effect_cols <- .lr_existing_cols(env_data, fixed_effect_cols)
  for (col in fixed_effect_cols) env_data[[col]] <- as.factor(env_data[[col]])
  nseg <- .lr_spatial_nseg(env_data)

  fixed_candidates <- list(c(curr_gen_col, fixed_effect_cols))
  if (length(fixed_effect_cols) > 0) {
    fixed_candidates <- c(
      fixed_candidates,
      lapply(fixed_effect_cols, function(x) c(curr_gen_col, x))
    )
  }
  fixed_candidates <- c(fixed_candidates, list(curr_gen_col))

  last_error <- NULL
  m <- NULL
  used_fixed_terms <- NULL
  for (candidate in fixed_candidates) {
    fixed_terms <- .lr_formula_terms(unique(candidate))
    m <- tryCatch({
      SpATS::SpATS(
        response = trait,
        genotype = prev_gen_col,
        genotype.as.random = TRUE,
        fixed = as.formula(paste("~", fixed_terms)),
        spatial = ~ SpATS::PSANOVA(Col_Num, Row_Num, nseg = c(nseg["col"], nseg["row"]), degree = 3, pord = 2),
        data = env_data,
        control = list(tolerance = 1e-03, monitoring = 0)
      )
    }, error = function(e) {
      last_error <<- e
      NULL
    })
    if (!is.null(m)) {
      used_fixed_terms <- fixed_terms
      break
    }
  }
  if (is.null(m)) stop(last_error$message)

  preds <- predict(m, which = curr_gen_col)
  se_df <- preds[, c(curr_gen_col, "standard.errors")]
  colnames(se_df) <- c("Current_Genotype", "SE")

  rot_df <- preds[, c(curr_gen_col, "predicted.values", "standard.errors")]
  colnames(rot_df) <- c("Current_Genotype", "Corrected_Wheat_Mean", "SE")

  raw_stats <- env_data |>
    dplyr::group_by(.data[[curr_gen_col]]) |>
    dplyr::summarize(
      Raw_Mean = mean(.data[[trait]], na.rm = TRUE),
      Facet = .lr_summarize_baseline(.data[[".Baseline_Group"]]),
      N_Plots = dplyr::n(),
      N_Previous_Genotypes = dplyr::n_distinct(.data[[prev_gen_col]]),
      N_Observed_Combos = dplyr::n_distinct(paste(.data[[prev_gen_col]], .data[[curr_gen_col]], sep = "-")),
      .groups = "drop"
    )
  colnames(raw_stats)[1] <- "Current_Genotype"

  baseline_mean <- .lr_precision_weighted_mean(rot_df$Corrected_Wheat_Mean, rot_df$SE)
  rot_df <- rot_df |>
    dplyr::left_join(raw_stats, by = "Current_Genotype") |>
    dplyr::mutate(
      Baseline_Mean = baseline_mean,
      Wheat_Rotational_Value = .data[["Corrected_Wheat_Mean"]] - .data[["Baseline_Mean"]],
      Wheat_Rotational_Value_Global = NA_real_,
      Network_Correction = .data[["Corrected_Wheat_Mean"]] - .data[["Raw_Mean"]],
      Total_Correction = Network_Correction
    )

  out <- list(rotational_values = rot_df, model = m)
  out$fixed_terms <- used_fixed_terms
  out
}

.lr_summarize_wheat_rotation_values <- function(env_data,
                                                trait,
                                                curr_gen_col,
                                                prev_gen_col,
                                                corrected_col,
                                                se_df,
                                                model) {
  rot_df <- env_data |>
    dplyr::group_by(.data[[curr_gen_col]]) |>
    dplyr::summarize(
      Corrected_Wheat_Mean = mean(.data[[corrected_col]], na.rm = TRUE),
      Raw_Mean = mean(.data[[trait]], na.rm = TRUE),
      Facet = .lr_summarize_baseline(.data[[".Baseline_Group"]]),
      N_Plots = dplyr::n(),
      N_Previous_Genotypes = dplyr::n_distinct(.data[[prev_gen_col]]),
      N_Observed_Combos = dplyr::n_distinct(paste(.data[[prev_gen_col]], .data[[curr_gen_col]], sep = "-")),
      .groups = "drop"
    )
  colnames(rot_df)[1] <- "Current_Genotype"

  if (!is.null(se_df)) {
    rot_df <- dplyr::left_join(rot_df, se_df, by = "Current_Genotype")
  } else {
    rot_df$SE <- NA_real_
  }

  baseline_mean <- .lr_precision_weighted_mean(rot_df$Corrected_Wheat_Mean, rot_df$SE)
  rot_df <- rot_df |>
    dplyr::mutate(
      Baseline_Mean = baseline_mean,
      Wheat_Rotational_Value = .data[["Corrected_Wheat_Mean"]] - .data[["Baseline_Mean"]],
      Wheat_Rotational_Value_Global = NA_real_,
      Network_Correction = .data[["Corrected_Wheat_Mean"]] - .data[["Raw_Mean"]],
      Total_Correction = Network_Correction
    )

  list(rotational_values = rot_df, model = model)
}
