#' Model Average Lentil Predecessor Effects on Wheat
#'
#' Estimates which lentil genotypes leave better average conditions for the
#' following wheat crop. Within each environment, the model corrects wheat
#' performance for wheat genotype, optional design terms, and field position.
#'
#' @param data Wheat/rotation data containing previous lentil and current wheat columns.
#' @param trait Wheat trait to analyze, for example `"Y.ADJ"`.
#' @param env_col Environment/site-year column.
#' @param prev_gen_col Previous crop genotype column, usually lentil.
#' @param curr_gen_col Current crop genotype column, usually wheat.
#' @param combo_col Lentil-wheat pair identifier. Created if missing.
#' @param spatial_cols Column names for row and column position.
#' @param baseline_col Network/facet column used to define the comparison
#'   baseline for legacy values.
#' @param fixed_effect_cols Optional additional fixed design/network terms.
#'   For rotation predecessor models this should usually be `"Rep_combo"`,
#'   because the same lentil-wheat combo is replicated within environment.
#'   Avoid `"Rep_gen"`, which indexes genotype-level replication and is used
#'   for independent genotype performance models.
#' @param type_col Optional column used to remove checks.
#' @param check_values Values in `type_col` to remove when `include_checks = FALSE`.
#' @param include_checks Logical. Include check plots in the model?
#' @param method `"SpATS"` or `"lme4"`.
#' @param fit_scope `"env_facet"` fits one model per ENV x Facet network.
#'   `"env_global"` fits one model per environment and centers by facet after
#'   fitting. The default is the local 10 x 10 network.
#' @param validate_design Logical. Run local factorial design diagnostics?
#'
#' @return A list with legacy values, plots, fitted models, and settings.
#' @export
model_predecessor_effect <- function(data,
                                     trait = "Y.ADJ",
                                     env_col = "ENV",
                                     prev_gen_col = "Lentil",
                                     curr_gen_col = "Wheat",
                                     combo_col = "Combo",
                                     spatial_cols = c("Row", "Col"),
                                     baseline_col = "Facet",
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

  data <- .lr_prepare_rotation_data(
    data = data,
    trait = trait,
    env_col = env_col,
    prev_gen_col = prev_gen_col,
    curr_gen_col = curr_gen_col,
    combo_col = combo_col,
    spatial_cols = spatial_cols,
    baseline_col = baseline_col,
    type_col = type_col,
    check_values = check_values,
    include_checks = include_checks
  )

  if (method == "SpATS" && !requireNamespace("SpATS", quietly = TRUE)) {
    warning("SpATS is not installed; falling back to method = 'lme4'.")
    method <- "lme4"
  }
  if (method == "lme4" && !requireNamespace("lme4", quietly = TRUE)) {
    stop("Package 'lme4' is required for method = 'lme4'.")
  }

  results_list <- list()
  model_list <- list()
  diagnostics_list <- list()
  design_audit <- NULL
  if (validate_design && !is.null(baseline_col) && baseline_col %in% names(data) && combo_col %in% names(data)) {
    design_audit <- audit_rotation_design(
      data = data,
      env_col = env_col,
      facet_col = baseline_col,
      lentil_col = prev_gen_col,
      wheat_col = curr_gen_col,
      combo_col = combo_col
    )
    if (!isTRUE(design_audit$design_ok)) {
      warning(
        "Rotation design audit found one or more incomplete or unexpected ENV x Facet networks. ",
        "Inspect result$design_audit before interpreting Legacy_Value rankings.",
        call. = FALSE
      )
    }
  }

  envs <- unique(data[[env_col]])

  for (env in envs) {
    env_data <- droplevels(data[data[[env_col]] == env & !is.na(data[[trait]]), ])
    if (nrow(env_data) == 0) next

    groups <- if (fit_scope == "env_facet") unique(as.character(env_data$.Baseline_Group)) else "All"
    for (baseline_group in groups) {
      group_data <- if (fit_scope == "env_facet") {
        droplevels(env_data[as.character(env_data$.Baseline_Group) == baseline_group, ])
      } else {
        env_data
      }
      if (nrow(group_data) == 0) next

      group_diag <- .lr_predecessor_group_diagnostics(
        env_data = group_data,
        trait = trait,
        prev_gen_col = prev_gen_col,
        curr_gen_col = curr_gen_col,
        combo_col = combo_col
      )
      group_diag$Environment <- env
      group_diag$Baseline_Group <- if (fit_scope == "env_facet") baseline_group else "All"

      fit <- tryCatch({
        if (method == "SpATS") {
          .lr_fit_predecessor_spats(
            env_data = group_data,
            trait = trait,
            prev_gen_col = prev_gen_col,
            curr_gen_col = curr_gen_col,
            combo_col = combo_col,
            spatial_cols = spatial_cols,
            baseline_col = baseline_col,
            fixed_effect_cols = fixed_effect_cols
          )
        } else {
          .lr_fit_predecessor_lme4(
            env_data = group_data,
            trait = trait,
            prev_gen_col = prev_gen_col,
            curr_gen_col = curr_gen_col,
            combo_col = combo_col,
            spatial_cols = spatial_cols,
            baseline_col = baseline_col,
            fixed_effect_cols = fixed_effect_cols
          )
        }
      }, error = function(e) {
        if (method == "SpATS" && requireNamespace("lme4", quietly = TRUE)) {
          warning(paste(
            "SpATS predecessor model failed for", env, baseline_group,
            "and will fall back to lme4:",
            e$message
          ))
          tryCatch({
            .lr_fit_predecessor_lme4(
              env_data = group_data,
              trait = trait,
              prev_gen_col = prev_gen_col,
              curr_gen_col = curr_gen_col,
              combo_col = combo_col,
              spatial_cols = spatial_cols,
              baseline_col = baseline_col,
              fixed_effect_cols = fixed_effect_cols
            )
          }, error = function(e2) {
            warning(paste("Predecessor model failed for", env, baseline_group, ":", e2$message))
            NULL
          })
        } else {
          warning(paste("Predecessor model failed for", env, baseline_group, ":", e$message))
          NULL
        }
      })

      if (is.null(fit)) {
        group_diag$Status <- "model_failed"
        diagnostics_list[[paste(env, baseline_group, "diagnostic", sep = "::")]] <- group_diag
        next
      }

      legacy_df <- fit$legacy_values
      legacy_df$Environment <- env
      legacy_df$Trait <- trait
      legacy_df$Method <- method
      legacy_df$Fit_Scope <- fit_scope
      legacy_df$Corrected_Mean_SE <- legacy_df$SE
      legacy_df$Legacy_Value_SE <- NA_real_
      group_diag$Status <- "fit"
      group_diag$Model_Status <- "fit"
      group_diag$Residual_Variance <- .lr_model_residual_variance(fit$model)
      group_diag$SE_Min <- .lr_safe_min(legacy_df$SE)
      group_diag$SE_Median <- .lr_safe_median(legacy_df$SE)
      group_diag$SE_Max <- .lr_safe_max(legacy_df$SE)
      group_diag$Network_Correction_Min <- .lr_safe_min(legacy_df$Network_Correction)
      group_diag$Network_Correction_Median <- .lr_safe_median(legacy_df$Network_Correction)
      group_diag$Network_Correction_Max <- .lr_safe_max(legacy_df$Network_Correction)

      result_key <- paste(env, baseline_group, sep = "::")
      results_list[[result_key]] <- legacy_df
      model_list[[result_key]] <- fit$model
      diagnostics_list[[paste(env, baseline_group, "diagnostic", sep = "::")]] <- group_diag
    }
  }

  final_df <- dplyr::bind_rows(results_list)
  diagnostics_df <- dplyr::bind_rows(diagnostics_list)
  if (nrow(final_df) > 0) {
    final_df <- final_df |>
      dplyr::group_by(.data[["Environment"]]) |>
      dplyr::mutate(
        .Global_Baseline_Mean = .lr_precision_weighted_mean(.data[["Corrected_Mean"]], .data[["SE"]]),
        Legacy_Value_Global = .data[["Corrected_Mean"]] - .data[[".Global_Baseline_Mean"]],
        Legacy_Pct_Global = (.data[["Legacy_Value_Global"]] / .data[[".Global_Baseline_Mean"]]) * 100
      ) |>
      dplyr::select(-dplyr::all_of(".Global_Baseline_Mean")) |>
      dplyr::ungroup() |>
      dplyr::group_by(.data[["Environment"]], .data[["Trait"]]) |>
      dplyr::mutate(Confidence_Class = .lr_confidence_class(.data[["SE"]])) |>
      dplyr::ungroup() |>
      dplyr::arrange(.data[["Environment"]], dplyr::desc(.data[["Legacy_Value"]]))
    final_df$Rank <- ave(
      -final_df$Legacy_Value,
      final_df$Environment,
      final_df$Trait,
      FUN = function(x) rank(x, ties.method = "first")
    )
  }

  p <- .lr_plot_predecessor(final_df, trait)
  ranked_p <- .lr_plot_predecessor_ranked(final_df, trait)
  correction_p <- .lr_plot_predecessor_correction(final_df, trait)
  gwas_df <- extract_predecessor_gwas_phenotypes(final_df)

  list(
    legacy_values = final_df,
    gwas_phenotypes = gwas_df,
    plot = p,
    ranked_plot = ranked_p,
    gwas_plot = ranked_p,
    correction_plot = correction_p,
    diagnostics = diagnostics_df,
    design_audit = design_audit,
    models = model_list,
    settings = list(
      trait = trait,
      method = method,
      fit_scope = fit_scope,
      include_checks = include_checks,
      baseline_col = baseline_col,
      question = "Which lentil genotypes are associated with higher or lower following-wheat performance relative to the other lentils tested within the same ENV x Facet network?"
    )
  )
}

#' Extract GWAS-Ready Lentil Legacy Values
#'
#' Converts predecessor-effect model output into one legacy value per previous
#' lentil genotype, environment, and trait. The default value is
#' `Legacy_Value`, which is the model-corrected deviation from the relevant
#' ENV x Facet baseline.
#'
#' @param x A result from `model_predecessor_effect()` or a legacy-value data frame.
#' @param phenotype_col Column to use as the GWAS phenotype.
#'
#' @return A data frame with one row per genotype x environment x trait.
#' @export
extract_predecessor_gwas_phenotypes <- function(x, phenotype_col = "Legacy_Value") {
  df <- .lr_get_legacy_df(x)
  if (nrow(df) == 0) return(df)

  .lr_check_columns(df, c("Previous_Genotype", "Environment", phenotype_col))

  if (!"Trait" %in% names(df)) df$Trait <- NA_character_
  if (!"N_Plots" %in% names(df)) df$N_Plots <- 1L
  if (!"N_Wheat_Partners" %in% names(df)) df$N_Wheat_Partners <- NA_integer_
  if (!"Baseline_Group" %in% names(df)) df$Baseline_Group <- "All"
  if (!"Corrected_Mean" %in% names(df)) df$Corrected_Mean <- NA_real_
  if (!"Baseline_Mean" %in% names(df)) df$Baseline_Mean <- NA_real_
  if (!"Legacy_Value_Global" %in% names(df)) df$Legacy_Value_Global <- NA_real_
  if (!"SE" %in% names(df)) df$SE <- NA_real_
  if (!"Method" %in% names(df)) df$Method <- NA_character_

  out <- df |>
    dplyr::group_by(
      Genotype = .data[["Previous_Genotype"]],
      Environment = .data[["Environment"]],
      Trait = .data[["Trait"]]
    ) |>
    dplyr::summarize(
      Legacy_Value = .lr_weighted_mean(.data[[phenotype_col]], .data[["N_Plots"]]),
      Corrected_Mean = .lr_weighted_mean(.data[["Corrected_Mean"]], .data[["N_Plots"]]),
      Facet_Baseline_Mean = .lr_weighted_mean(.data[["Baseline_Mean"]], .data[["N_Plots"]]),
      Global_Legacy_Value = .lr_weighted_mean(.data[["Legacy_Value_Global"]], .data[["N_Plots"]]),
      SE = .lr_weighted_mean(.data[["SE"]], .data[["N_Plots"]]),
      N_Plots = sum(.data[["N_Plots"]], na.rm = TRUE),
      N_Wheat_Partners = sum(unique(.data[["N_Wheat_Partners"]]), na.rm = TRUE),
      Baseline_Group = paste(unique(as.character(.data[["Baseline_Group"]])), collapse = ";"),
      N_Baseline_Groups = dplyr::n_distinct(.data[["Baseline_Group"]]),
      Method = paste(unique(as.character(.data[["Method"]])), collapse = ";"),
      Phenotype_Column = phenotype_col,
      Phenotype_Interpretation = "Model-corrected wheat response deviation from the ENV x Facet baseline; positive values indicate better wheat performance after that lentil genotype.",
      .groups = "drop"
    ) |>
    dplyr::group_by(.data[["Environment"]], .data[["Trait"]]) |>
    dplyr::mutate(
      Rank = rank(-.data[["Legacy_Value"]], ties.method = "first")
    ) |>
    dplyr::ungroup() |>
    dplyr::arrange(.data[["Environment"]], .data[["Trait"]], .data[["Rank"]])

  out
}

#' Plot GWAS-Ready Lentil Legacy Values
#'
#' Shows a single ranked list of lentil genotypes per environment using the
#' facet-corrected `Legacy_Value`. Bars are colored by the wheat-partner
#' facet to preserve design context without splitting the ranking into panels.
#'
#' @param x A result from `model_predecessor_effect()` or a legacy-value data frame.
#' @param trait Optional trait label for the plot subtitle.
#' @param n Optional number of top and bottom genotypes to show per environment.
#' @param show_se Logical. Draw model prediction standard errors?
#'
#' @return A ggplot object.
#' @export
plot_predecessor_gwas_phenotypes <- function(x, trait = NULL, n = NULL, show_se = FALSE) {
  df <- .lr_get_legacy_df(x)
  if (nrow(df) == 0) return(NULL)
  df <- .lr_as_legacy_plot_df(df)
  if (is.null(trait) && "Trait" %in% names(df)) {
    trait <- paste(unique(as.character(df$Trait)), collapse = ", ")
  }
  .lr_plot_predecessor_ranked(df, trait = trait, n = n, show_se = show_se)
}

#' Model Lentil-Wheat Pair Compatibility
#'
#' Estimates which observed lentil-wheat genotype pairs perform better or worse
#' than expected from their average lentil and wheat main effects. The `Combo`
#' random effect is the pair-specific compatibility value.
#'
#' @param data Wheat/rotation data.
#' @param trait Wheat trait to analyze.
#' @param env_col Environment/site-year column.
#' @param prev_gen_col Previous crop genotype column, usually lentil.
#' @param curr_gen_col Current crop genotype column, usually wheat.
#' @param combo_col Lentil-wheat pair identifier. Created if missing.
#' @param spatial_cols Column names for row and column position.
#' @param baseline_col Optional network/facet column added to output and plots.
#' @param fixed_effect_cols Optional additional fixed terms. The default
#'   `"Rep_combo"` adjusts for the two replicated plots of the same observed
#'   lentil-wheat combo within environment. Avoid `"Rep_gen"` in rotation
#'   models because it indexes genotype-level replication, not combo-level
#'   replication.
#' @param random_effect_cols Optional design random effects, for example `"Block"`.
#' @param type_col Optional column used to remove checks.
#' @param check_values Values in `type_col` to remove when `include_checks = FALSE`.
#' @param include_checks Logical. Include check plots in the model?
#' @param min_combo_reps Minimum observed plots per lentil-wheat pair to label a
#'   pair as fully supported in the output.
#' @param compute_compatibility_pct Logical. Compute percent compatibility?
#' @param pct_min_denominator Minimum absolute additive mean used for stable
#'   percent calculations.
#' @param validate_design Logical. Run local factorial design diagnostics?
#'
#' @return A list with pair compatibility values, plots, fitted models, and settings.
#' @export
model_pair_compatibility <- function(data,
                                     trait = "Y.ADJ",
                                     env_col = "ENV",
                                     prev_gen_col = "Lentil",
                                     curr_gen_col = "Wheat",
                                     combo_col = "Combo",
                                     spatial_cols = c("Row", "Col"),
                                     baseline_col = "Facet",
                                     fixed_effect_cols = c("Rep_combo"),
                                     random_effect_cols = c("Block"),
                                     type_col = "Type",
                                     check_values = "Check",
                                     include_checks = FALSE,
                                     min_combo_reps = 2,
                                     compute_compatibility_pct = TRUE,
                                     pct_min_denominator = 1e-6,
                                     validate_design = TRUE) {

  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Package 'lme4' is required for pair compatibility models.")
  }
  .lr_warn_rotation_rep_terms(fixed_effect_cols, random_effect_cols)

  data <- .lr_prepare_rotation_data(
    data = data,
    trait = trait,
    env_col = env_col,
    prev_gen_col = prev_gen_col,
    curr_gen_col = curr_gen_col,
    combo_col = combo_col,
    spatial_cols = spatial_cols,
    baseline_col = baseline_col,
    type_col = type_col,
    check_values = check_values,
    include_checks = include_checks
  )

  results_list <- list()
  model_list <- list()
  diagnostics_list <- list()
  design_audit <- NULL
  if (validate_design && !is.null(baseline_col) && baseline_col %in% names(data) && combo_col %in% names(data)) {
    design_audit <- audit_rotation_design(
      data = data,
      env_col = env_col,
      facet_col = baseline_col,
      lentil_col = prev_gen_col,
      wheat_col = curr_gen_col,
      combo_col = combo_col
    )
    if (!isTRUE(design_audit$design_ok)) {
      warning(
        "Rotation design audit found one or more incomplete or unexpected ENV x Facet networks. ",
        "Inspect result$design_audit and result$diagnostics before interpreting pair compatibility.",
        call. = FALSE
      )
    }
  }

  envs <- unique(data[[env_col]])

  for (env in envs) {
    env_data <- droplevels(data[data[[env_col]] == env & !is.na(data[[trait]]), ])
    if (nrow(env_data) == 0) next

    baseline_groups <- unique(as.character(env_data$.Baseline_Group))

    for (baseline_group in baseline_groups) {
      group_data <- droplevels(env_data[as.character(env_data$.Baseline_Group) == baseline_group, ])
      if (nrow(group_data) == 0) next
      group_diag <- .lr_pair_group_diagnostics(
        env_data = group_data,
        trait = trait,
        prev_gen_col = prev_gen_col,
        curr_gen_col = curr_gen_col,
        combo_col = combo_col
      )
      group_diag$Environment <- env
      group_diag$Baseline_Group <- baseline_group

      if (dplyr::n_distinct(group_data[[prev_gen_col]]) < 2 ||
          dplyr::n_distinct(group_data[[curr_gen_col]]) < 2 ||
          dplyr::n_distinct(group_data[[combo_col]]) < 2) {
        group_diag$Status <- "skipped_insufficient_genotype_levels"
        diagnostics_list[[paste(env, baseline_group, "diagnostic", sep = "::")]] <- group_diag
        next
      }

      fit <- tryCatch({
        .lr_fit_pair_lme4(
          env_data = group_data,
          trait = trait,
          prev_gen_col = prev_gen_col,
          curr_gen_col = curr_gen_col,
          combo_col = combo_col,
          spatial_cols = spatial_cols,
          baseline_col = baseline_col,
          fixed_effect_cols = fixed_effect_cols,
          random_effect_cols = random_effect_cols,
          compute_compatibility_pct = compute_compatibility_pct,
          pct_min_denominator = pct_min_denominator
        )
      }, error = function(e) {
        warning(paste("Pair compatibility model failed for", env, baseline_group, ":", e$message))
        NULL
      })

      if (is.null(fit)) {
        group_diag$Status <- "model_failed"
        diagnostics_list[[paste(env, baseline_group, "diagnostic", sep = "::")]] <- group_diag
        next
      }

      pair_df <- fit$compatibility
      pair_df$Environment <- env
      pair_df$Trait <- trait
      pair_df$Pair_Support <- ifelse(pair_df$N_Plots >= min_combo_reps, "supported", "low_replication")
      pair_df$Model_Singular <- fit$diagnostics$Model_Singular
      pair_df$Combo_Variance <- fit$diagnostics$Combo_Variance
      pair_df$Residual_Variance <- fit$diagnostics$Residual_Variance
      pair_df$Local_Factorial_Complete <- fit$diagnostics$Local_Factorial_Complete
      pair_df$Expected_10x10 <- fit$diagnostics$Expected_10x10
      pair_df <- pair_df[order(pair_df$Compatibility_Value, decreasing = TRUE), ]
      pair_df$Rank <- seq_len(nrow(pair_df))

      result_key <- paste(env, baseline_group, sep = "::")
      results_list[[result_key]] <- pair_df
      model_list[[result_key]] <- fit$model
      fit$diagnostics$Environment <- env
      fit$diagnostics$Baseline_Group <- baseline_group
      diagnostics_list[[paste(env, baseline_group, "diagnostic", sep = "::")]] <- fit$diagnostics
    }
  }

  final_df <- dplyr::bind_rows(results_list)
  diagnostics_df <- dplyr::bind_rows(diagnostics_list)
  p <- plot_pair_compatibility_heatmap(final_df, trait)
  ranked_p <- plot_pair_compatibility_ranked(final_df, trait)

  list(
    compatibility = final_df,
    diagnostics = diagnostics_df,
    plot = p,
    heatmap = p,
    ranked_plot = ranked_p,
    design_audit = design_audit,
    models = model_list,
    settings = list(
      trait = trait,
      include_checks = include_checks,
      baseline_col = baseline_col,
      fixed_effect_cols = fixed_effect_cols,
      min_combo_reps = min_combo_reps,
      compute_compatibility_pct = compute_compatibility_pct,
      pct_min_denominator = pct_min_denominator,
      question = "Which lentil-wheat genotype pairs are specifically compatible?"
    )
  )
}

#' Plot Lentil-Wheat Pair Compatibility as a Heatmap
#'
#' Displays observed lentil-wheat compatibility values as a matrix of previous
#' lentil genotype by current wheat genotype. This is usually easier to read
#' than plotting long combo names on the axis.
#'
#' @param x A result from `model_pair_compatibility()` or its compatibility table.
#' @param trait Optional trait label for the subtitle.
#' @param value_col Compatibility column to plot.
#' @param signal_tol Absolute compatibility threshold below which an ENV x Facet
#'   is dropped from the plot.
#' @param drop_no_signal Logical. Drop ENV x Facet panels with no detectable
#'   pair-specific signal?
#'
#' @return A ggplot object.
#' @export
plot_pair_compatibility_heatmap <- function(x,
                                            trait = NULL,
                                            value_col = "Compatibility_Value",
                                            signal_tol = 1e-8,
                                            drop_no_signal = TRUE) {
  df <- .lr_get_pair_df(x)
  if (nrow(df) == 0) return(NULL)
  .lr_check_columns(df, c("Previous_Genotype", "Current_Genotype", "Environment", "Baseline_Group", value_col))
  if (is.null(trait) && "Trait" %in% names(df)) trait <- paste(unique(as.character(df$Trait)), collapse = ", ")

  plot_df <- df |>
    dplyr::group_by(.data[["Environment"]], .data[["Baseline_Group"]]) |>
    dplyr::mutate(
      Max_Abs_Compatibility = max(abs(.data[[value_col]]), na.rm = TRUE),
      Pair_Signal = Max_Abs_Compatibility > signal_tol
    ) |>
    dplyr::ungroup()

  if (drop_no_signal) {
    plot_df <- plot_df |>
      dplyr::filter(.data[["Pair_Signal"]])
  }
  if (nrow(plot_df) == 0) return(NULL)

  plot_df <- plot_df |>
    dplyr::mutate(
      Panel = paste(.data[["Environment"]], .data[["Baseline_Group"]], sep = " | "),
      Previous_Genotype = factor(as.character(.data[["Previous_Genotype"]])),
      Current_Genotype = factor(as.character(.data[["Current_Genotype"]]))
    )

  ggplot2::ggplot(plot_df, ggplot2::aes(
    x = .data[["Current_Genotype"]],
    y = .data[["Previous_Genotype"]],
    fill = .data[[value_col]]
  )) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::facet_wrap(~Panel, scales = "free", drop = TRUE) +
    ggplot2::scale_fill_gradient2(
      low = "firebrick",
      mid = "white",
      high = "forestgreen",
      midpoint = 0,
      name = "Pair effect"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Lentil-Wheat Pair Compatibility",
      subtitle = paste("Showing only ENV x Facet groups with detectable pair-specific signal", if (!is.null(trait)) paste("| trait:", trait) else ""),
      caption = "Compatibility_Value is a pair-specific deviation within the observed local ENV x Facet network; unobserved global 100 x 100 combinations are not created.",
      x = "Current wheat genotype",
      y = "Previous lentil genotype"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
      axis.text.y = ggplot2::element_text(size = 6),
      panel.grid = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = 8)
    )
}

#' Plot Top Lentil-Wheat Pair Compatibility Values
#'
#' Shows the strongest positive and negative pair effects without drawing all
#' observed pair labels in every facet.
#'
#' @param x A result from `model_pair_compatibility()` or its compatibility table.
#' @param trait Optional trait label for the subtitle.
#' @param n Number of strongest absolute pair effects per ENV x Facet panel.
#' @param show_se Logical. Draw conditional random-effect standard errors?
#' @param signal_tol Absolute compatibility threshold below which an ENV x Facet
#'   is dropped from the plot.
#' @param drop_no_signal Logical. Drop ENV x Facet panels with no detectable
#'   pair-specific signal?
#'
#' @return A ggplot object.
#' @export
plot_pair_compatibility_ranked <- function(x,
                                           trait = NULL,
                                           n = 12,
                                           show_se = FALSE,
                                           signal_tol = 1e-8,
                                           drop_no_signal = TRUE) {
  df <- .lr_get_pair_df(x)
  if (nrow(df) == 0) return(NULL)
  .lr_check_columns(df, c("Combo", "Environment", "Baseline_Group", "Compatibility_Value"))
  if (is.null(trait) && "Trait" %in% names(df)) trait <- paste(unique(as.character(df$Trait)), collapse = ", ")

  group_cols <- c("Environment", "Baseline_Group")
  group_signal <- df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarize(
      Max_Abs_Compatibility = max(abs(.data[["Compatibility_Value"]]), na.rm = TRUE),
      Pair_Signal = Max_Abs_Compatibility > signal_tol,
      .groups = "drop"
    )

  signal_df <- df |>
    dplyr::inner_join(group_signal, by = group_cols) |>
    dplyr::filter(.data[["Pair_Signal"]]) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::slice_max(abs(.data[["Compatibility_Value"]]), n = n, with_ties = FALSE) |>
    dplyr::ungroup()

  if (!drop_no_signal) {
    no_signal_df <- group_signal |>
      dplyr::filter(!.data[["Pair_Signal"]]) |>
      dplyr::transmute(
        Environment = .data[["Environment"]],
        Baseline_Group = .data[["Baseline_Group"]],
        Combo = "No detectable pair-specific signal",
        Compatibility_Value = 0,
        SE = NA_real_,
        Pair_Signal = FALSE,
        Max_Abs_Compatibility = .data[["Max_Abs_Compatibility"]]
      )
    signal_df <- dplyr::bind_rows(signal_df, no_signal_df)
  }

  if (nrow(signal_df) == 0) return(NULL)

  plot_df <- signal_df |>
    dplyr::ungroup() |>
    dplyr::arrange(.data[["Environment"]], .data[["Baseline_Group"]], .data[["Compatibility_Value"]]) |>
    dplyr::mutate(
      Panel = paste(.data[["Environment"]], .data[["Baseline_Group"]], sep = " | "),
      Plot_Combo = paste(.data[["Panel"]], .data[["Combo"]], sep = "___"),
      Plot_Combo = factor(.data[["Plot_Combo"]], levels = unique(.data[["Plot_Combo"]])),
      Point_Class = dplyr::case_when(
        !.data[["Pair_Signal"]] ~ "No detectable pair signal",
        .data[["Compatibility_Value"]] > 0 ~ "Positive",
        TRUE ~ "Negative"
      )
    )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[["Plot_Combo"]], y = .data[["Compatibility_Value"]])) +
    ggplot2::geom_segment(
      ggplot2::aes(xend = .data[["Plot_Combo"]], y = 0, yend = .data[["Compatibility_Value"]]),
      color = "gray72",
      linewidth = 0.4
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = .data[["Point_Class"]]),
      size = 2.1,
      alpha = 0.9
    )

  if (show_se && "SE" %in% names(plot_df)) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data[["Compatibility_Value"]] - .data[["SE"]], ymax = .data[["Compatibility_Value"]] + .data[["SE"]]),
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
    ggplot2::scale_color_manual(
      values = c(
        "Positive" = "forestgreen",
        "Negative" = "firebrick",
        "No detectable pair signal" = "gray35"
      ),
      guide = "none"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Top Lentil-Wheat Pair Compatibility Values",
      subtitle = paste("Showing only ENV x Facet groups with detectable pair-specific signal", if (!is.null(trait)) paste("| trait:", trait) else ""),
      caption = "Rankings are within observed ENV x Facet networks. Compatibility_Value is preferred over Compatibility_Pct for centered, ordinal, or near-zero traits.",
      x = "Observed lentil-wheat pair",
      y = "Pair effect beyond additive lentil and wheat effects"
    )
}

.lr_prepare_rotation_data <- function(data,
                                      trait,
                                      env_col,
                                      prev_gen_col,
                                      curr_gen_col,
                                      combo_col = NULL,
                                      spatial_cols = c("Row", "Col"),
                                      baseline_col = NULL,
                                      type_col = NULL,
                                      check_values = "Check",
                                      include_checks = FALSE) {
  required_cols <- c(trait, env_col, prev_gen_col, curr_gen_col, spatial_cols)
  if (!is.null(baseline_col) && baseline_col %in% names(data)) {
    required_cols <- c(required_cols, baseline_col)
  }
  if (!is.null(type_col) && type_col %in% names(data)) required_cols <- c(required_cols, type_col)
  .lr_check_columns(data, required_cols)

  out <- data
  if (!include_checks && !is.null(type_col) && type_col %in% names(out)) {
    out <- out[!(out[[type_col]] %in% check_values), ]
  }

  out <- out[!is.na(out[[prev_gen_col]]) & !is.na(out[[curr_gen_col]]), ]

  if (!is.null(combo_col)) {
    if (!combo_col %in% names(out)) {
      out[[combo_col]] <- paste(out[[prev_gen_col]], out[[curr_gen_col]], sep = "-")
    }
    out[[combo_col]] <- as.factor(out[[combo_col]])
  }

  out[[env_col]] <- as.factor(out[[env_col]])
  out[[prev_gen_col]] <- as.factor(out[[prev_gen_col]])
  out[[curr_gen_col]] <- as.factor(out[[curr_gen_col]])
  if (!is.null(baseline_col) && baseline_col %in% names(out)) {
    out$.Baseline_Group <- as.factor(out[[baseline_col]])
  } else {
    out$.Baseline_Group <- factor("All")
  }
  out$Row_Num <- as.numeric(as.character(out[[spatial_cols[1]]]))
  out$Col_Num <- as.numeric(as.character(out[[spatial_cols[2]]]))

  out <- out[!is.na(out$Row_Num) & !is.na(out$Col_Num), ]
  out$Row_Scaled <- as.numeric(scale(out$Row_Num))
  out$Col_Scaled <- as.numeric(scale(out$Col_Num))
  droplevels(out)
}

.lr_check_columns <- function(data, cols) {
  missing <- setdiff(unique(cols), names(data))
  if (length(missing) > 0) {
    stop(paste("Missing required columns:", paste(missing, collapse = ", ")))
  }
  invisible(TRUE)
}

.lr_get_legacy_df <- function(x) {
  if (is.data.frame(x)) return(x)
  if (is.list(x) && "legacy_values" %in% names(x)) return(x$legacy_values)
  stop("Expected a legacy-value data frame or a result from model_predecessor_effect().")
}

.lr_get_pair_df <- function(x) {
  if (is.data.frame(x)) return(x)
  if (is.list(x) && "compatibility" %in% names(x)) return(x$compatibility)
  stop("Expected a pair-compatibility data frame or a result from model_pair_compatibility().")
}

.lr_as_legacy_plot_df <- function(df) {
  out <- df
  if (!"Previous_Genotype" %in% names(out) && "Genotype" %in% names(out)) {
    out$Previous_Genotype <- out$Genotype
  }
  if (!"SE" %in% names(out)) out$SE <- NA_real_
  if (!"Baseline_Group" %in% names(out)) out$Baseline_Group <- "All"
  out
}

.lr_weighted_mean <- function(x, w) {
  ok <- !is.na(x)
  if (!any(ok)) return(NA_real_)
  w <- w[ok]
  x <- x[ok]
  if (all(is.na(w)) || sum(w, na.rm = TRUE) <= 0) return(mean(x, na.rm = TRUE))
  stats::weighted.mean(x, w = w, na.rm = TRUE)
}

.lr_precision_weighted_mean <- function(x, se) {
  ok <- !is.na(x)
  if (!any(ok)) return(NA_real_)
  x <- x[ok]
  se <- se[ok]
  if (length(se) != length(x) || all(is.na(se)) || all(se <= 0, na.rm = TRUE)) {
    return(mean(x, na.rm = TRUE))
  }
  w <- ifelse(is.na(se) | se <= 0, NA_real_, 1 / (se^2))
  if (all(is.na(w)) || sum(w, na.rm = TRUE) <= 0) return(mean(x, na.rm = TRUE))
  stats::weighted.mean(x, w = w, na.rm = TRUE)
}

.lr_warn_rotation_rep_terms <- function(fixed_effect_cols = character(0),
                                        random_effect_cols = character(0)) {
  rep_terms <- c(fixed_effect_cols, random_effect_cols)
  if ("Rep_gen" %in% rep_terms) {
    warning(
      "Rep_gen indexes genotype-level replication and is intended for independent genotype performance models. ",
      "Rotation predecessor and pair-compatibility models should usually use Rep_combo, ",
      "because the same lentil-wheat Combo is replicated within ENV.",
      call. = FALSE
    )
  }
}

.lr_existing_cols <- function(data, cols) {
  cols[cols %in% names(data)]
}

.lr_spatial_nseg <- function(env_data) {
  max_r <- max(env_data$Row_Num, na.rm = TRUE)
  max_c <- max(env_data$Col_Num, na.rm = TRUE)
  c(
    row = max(1, min(round(max_r / 4), 35)),
    col = max(1, min(round(max_c / 4), 35))
  )
}

.lr_formula_terms <- function(cols) {
  cols <- cols[!is.na(cols) & nzchar(cols)]
  if (length(cols) == 0) return("1")
  paste(cols, collapse = " + ")
}

.lr_fit_predecessor_spats <- function(env_data,
                                      trait,
                                      prev_gen_col,
                                      curr_gen_col,
                                      combo_col,
                                      spatial_cols,
                                      baseline_col,
                                      fixed_effect_cols) {
  fixed_effect_cols <- .lr_existing_cols(env_data, fixed_effect_cols)
  for (col in fixed_effect_cols) env_data[[col]] <- as.factor(env_data[[col]])

  nseg <- .lr_spatial_nseg(env_data)

  fixed_candidates <- list(c(prev_gen_col, fixed_effect_cols))
  if (length(fixed_effect_cols) > 0) {
    fixed_candidates <- c(
      fixed_candidates,
      lapply(fixed_effect_cols, function(x) c(prev_gen_col, x))
    )
  }
  fixed_candidates <- c(fixed_candidates, list(prev_gen_col))

  last_error <- NULL
  m <- NULL
  used_fixed_terms <- NULL

  for (candidate in fixed_candidates) {
    fixed_terms <- .lr_formula_terms(unique(candidate))

    m <- tryCatch({
      SpATS::SpATS(
        response = trait,
        genotype = curr_gen_col,
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

  dropped_terms <- setdiff(fixed_effect_cols, strsplit(used_fixed_terms, " \\+ ")[[1]])
  if (length(dropped_terms) > 0) {
    warning(paste(
      "SpATS fixed-effect matrix was rank deficient; refit with fixed terms:",
      used_fixed_terms
    ))
  }

  preds <- predict(m, which = prev_gen_col)
  legacy_df <- preds[, c(prev_gen_col, "predicted.values", "standard.errors")]
  colnames(legacy_df) <- c("Previous_Genotype", "Corrected_Mean", "SE")

  raw_stats <- env_data |>
    dplyr::group_by(.data[[prev_gen_col]]) |>
    dplyr::summarize(
      Raw_Mean = mean(.data[[trait]], na.rm = TRUE),
      N_Plots = dplyr::n(),
      N_Wheat_Partners = dplyr::n_distinct(.data[[curr_gen_col]]),
      N_Observed_Combos = if (combo_col %in% names(env_data)) dplyr::n_distinct(.data[[combo_col]]) else NA_integer_,
      Baseline_Group = .lr_summarize_baseline(.data[[".Baseline_Group"]]),
      N_Baseline_Groups = .lr_count_baseline(.data[[".Baseline_Group"]]),
      .groups = "drop"
    )
  colnames(raw_stats)[1] <- "Previous_Genotype"

  grand_mean <- .lr_precision_weighted_mean(legacy_df$Corrected_Mean, legacy_df$SE)
  legacy_df <- legacy_df |>
    dplyr::left_join(raw_stats, by = "Previous_Genotype") |>
    dplyr::group_by(Baseline_Group) |>
    dplyr::mutate(
      Baseline_Mean = .lr_precision_weighted_mean(Corrected_Mean, SE),
      Legacy_Value_Global = Corrected_Mean - grand_mean,
      Legacy_Pct_Global = (Legacy_Value_Global / grand_mean) * 100,
      Legacy_Value = Corrected_Mean - Baseline_Mean,
      Legacy_Pct = (Legacy_Value / Baseline_Mean) * 100,
      Network_Correction = Corrected_Mean - Raw_Mean,
      Total_Correction = Network_Correction
    ) |>
    dplyr::ungroup()

  list(legacy_values = legacy_df, model = m, fixed_terms = used_fixed_terms)
}

.lr_fit_predecessor_lme4 <- function(env_data,
                                     trait,
                                     prev_gen_col,
                                     curr_gen_col,
                                     combo_col,
                                     spatial_cols,
                                     baseline_col,
                                     fixed_effect_cols) {
  fixed_effect_cols <- .lr_existing_cols(env_data, fixed_effect_cols)
  for (col in fixed_effect_cols) env_data[[col]] <- as.factor(env_data[[col]])

  spatial_terms <- "poly(Row_Scaled, 2, raw = TRUE) + poly(Col_Scaled, 2, raw = TRUE) + Row_Scaled:Col_Scaled"
  fixed_terms <- .lr_formula_terms(c(prev_gen_col, fixed_effect_cols, spatial_terms))
  formula_str <- paste(trait, "~", fixed_terms, "+ (1|", curr_gen_col, ")")
  m <- lme4::lmer(as.formula(formula_str), data = env_data)

  env_data$Corrected <- predict(m, re.form = NA)
  legacy_df <- env_data |>
    dplyr::group_by(.data[[prev_gen_col]]) |>
    dplyr::summarize(
      Corrected_Mean = mean(Corrected, na.rm = TRUE),
      Raw_Mean = mean(.data[[trait]], na.rm = TRUE),
      N_Plots = dplyr::n(),
      N_Wheat_Partners = dplyr::n_distinct(.data[[curr_gen_col]]),
      N_Observed_Combos = if (combo_col %in% names(env_data)) dplyr::n_distinct(.data[[combo_col]]) else NA_integer_,
      Baseline_Group = .lr_summarize_baseline(.data[[".Baseline_Group"]]),
      N_Baseline_Groups = .lr_count_baseline(.data[[".Baseline_Group"]]),
      .groups = "drop"
    )
  colnames(legacy_df)[1] <- "Previous_Genotype"

  grand_mean <- .lr_precision_weighted_mean(legacy_df$Corrected_Mean, legacy_df$SE)
  legacy_df <- legacy_df |>
    dplyr::group_by(Baseline_Group) |>
    dplyr::mutate(
      SE = NA_real_,
      Baseline_Mean = .lr_precision_weighted_mean(Corrected_Mean, SE),
      Legacy_Value_Global = Corrected_Mean - grand_mean,
      Legacy_Pct_Global = (Legacy_Value_Global / grand_mean) * 100,
      Legacy_Value = Corrected_Mean - Baseline_Mean,
      Legacy_Pct = (Legacy_Value / Baseline_Mean) * 100,
      Network_Correction = Corrected_Mean - Raw_Mean,
      Total_Correction = Network_Correction
    ) |>
    dplyr::ungroup()

  list(legacy_values = legacy_df, model = m)
}

.lr_fit_pair_lme4 <- function(env_data,
                              trait,
                              prev_gen_col,
                              curr_gen_col,
                              combo_col,
                              spatial_cols,
                              baseline_col,
                              fixed_effect_cols,
                              random_effect_cols,
                              compute_compatibility_pct,
                              pct_min_denominator) {
  fixed_effect_cols <- .lr_existing_cols(env_data, fixed_effect_cols)
  for (col in fixed_effect_cols) env_data[[col]] <- as.factor(env_data[[col]])

  random_effect_cols <- .lr_existing_cols(env_data, random_effect_cols)
  for (col in random_effect_cols) env_data[[col]] <- as.factor(env_data[[col]])

  spatial_terms <- "poly(Row_Scaled, 2, raw = TRUE) + poly(Col_Scaled, 2, raw = TRUE) + Row_Scaled:Col_Scaled"
  random_terms <- c(combo_col, random_effect_cols)
  random_terms <- paste0("(1|", random_terms, ")", collapse = " + ")
  fixed_terms <- .lr_formula_terms(c(prev_gen_col, curr_gen_col, fixed_effect_cols, spatial_terms))
  formula_str <- paste(trait, "~", fixed_terms, "+", random_terms)

  m <- lme4::lmer(as.formula(formula_str), data = env_data)
  model_diag <- .lr_pair_group_diagnostics(
    env_data = env_data,
    trait = trait,
    prev_gen_col = prev_gen_col,
    curr_gen_col = curr_gen_col,
    combo_col = combo_col
  )
  model_diag$Status <- "fit"
  model_diag$Formula <- formula_str
  model_diag$Model_Singular <- lme4::isSingular(m)
  vc <- as.data.frame(lme4::VarCorr(m))
  combo_var <- vc[vc$grp == combo_col, "vcov"]
  if (length(combo_var) == 0) combo_var <- NA_real_
  model_diag$Combo_Variance <- combo_var[[1]]
  model_diag$Residual_Variance <- stats::sigma(m)^2
  model_diag$Status <- .lr_pair_status(model_diag)

  combo_re <- lme4::ranef(m, condVar = TRUE)[[combo_col]]
  post_var <- attr(combo_re, "postVar")
  se <- rep(NA_real_, nrow(combo_re))
  if (!is.null(post_var)) se <- sqrt(post_var[1, 1, ])

  combo_effects <- data.frame(
    Combo = rownames(combo_re),
    Compatibility_Value = combo_re[[1]],
    SE = se,
    row.names = NULL
  )

  env_data$Expected_Additive <- predict(m, re.form = NA)
  env_data$Conditional_Predicted <- predict(m)
  env_data$Compatibility_Value <- combo_effects$Compatibility_Value[match(as.character(env_data[[combo_col]]), combo_effects$Combo)]
  env_data$Marginal_Predicted_With_Pair <- env_data$Expected_Additive + env_data$Compatibility_Value

  pair_stats <- env_data |>
    dplyr::group_by(.data[[combo_col]]) |>
    dplyr::summarize(
      Previous_Genotype = dplyr::first(.data[[prev_gen_col]]),
      Current_Genotype = dplyr::first(.data[[curr_gen_col]]),
      Baseline_Group = .lr_summarize_baseline(.data[[".Baseline_Group"]]),
      Raw_Mean = mean(.data[[trait]], na.rm = TRUE),
      Expected_Additive_Mean = mean(Expected_Additive, na.rm = TRUE),
      Conditional_Predicted_Mean = mean(Conditional_Predicted, na.rm = TRUE),
      Marginal_Corrected_Pair_Mean = mean(Marginal_Predicted_With_Pair, na.rm = TRUE),
      N_Plots = dplyr::n(),
      .groups = "drop"
    )
  colnames(pair_stats)[1] <- "Combo"

  compatibility <- pair_stats |>
    dplyr::left_join(combo_effects, by = "Combo") |>
    dplyr::mutate(
      Corrected_Pair_Mean = .data[["Marginal_Corrected_Pair_Mean"]],
      Compatibility_Pct = if (compute_compatibility_pct) {
        ifelse(
          abs(.data[["Expected_Additive_Mean"]]) < pct_min_denominator,
          NA_real_,
          (.data[["Compatibility_Value"]] / .data[["Expected_Additive_Mean"]]) * 100
        )
      } else {
        NA_real_
      }
    )

  if (compute_compatibility_pct && any(abs(compatibility$Expected_Additive_Mean) < pct_min_denominator, na.rm = TRUE)) {
    warning(
      "Compatibility_Pct can be unstable when Expected_Additive_Mean is near zero. ",
      "Use Compatibility_Value as the primary output for centered, ordinal, or near-zero traits.",
      call. = FALSE
    )
  }

  list(compatibility = compatibility, model = m, diagnostics = model_diag)
}

.lr_pair_group_diagnostics <- function(env_data,
                                       trait,
                                       prev_gen_col,
                                       curr_gen_col,
                                       combo_col) {
  combo_counts <- table(env_data[[combo_col]])
  if (length(combo_counts) == 0) combo_counts <- 0L
  n_previous <- dplyr::n_distinct(env_data[[prev_gen_col]])
  n_current <- dplyr::n_distinct(env_data[[curr_gen_col]])
  n_combos <- dplyr::n_distinct(env_data[[combo_col]])
  expected_combos <- n_previous * n_current
  min_reps <- min(as.integer(combo_counts), na.rm = TRUE)

  data.frame(
    Trait = trait,
    N_Plots = nrow(env_data),
    N_Previous_Genotypes = n_previous,
    N_Current_Genotypes = n_current,
    N_Combos = n_combos,
    Expected_Combos = expected_combos,
    Observed_Combos = n_combos,
    Missing_Combos = expected_combos - n_combos,
    Local_Factorial_Complete = n_combos == expected_combos,
    Expected_10x10 = n_previous == 10 && n_current == 10 && n_combos == 100,
    Min_Reps_Per_Combo = min_reps,
    Median_Reps_Per_Combo = stats::median(as.integer(combo_counts), na.rm = TRUE),
    Max_Reps_Per_Combo = max(as.integer(combo_counts), na.rm = TRUE),
    N_Single_Rep_Combos = sum(as.integer(combo_counts) < 2, na.rm = TRUE),
    Status = "not_fit",
    Formula = NA_character_,
    Model_Singular = NA,
    Combo_Variance = NA_real_,
    Residual_Variance = NA_real_
  )
}

.lr_predecessor_group_diagnostics <- function(env_data,
                                              trait,
                                              prev_gen_col,
                                              curr_gen_col,
                                              combo_col) {
  combo_counts <- if (combo_col %in% names(env_data)) table(env_data[[combo_col]]) else integer(0)
  n_previous <- dplyr::n_distinct(env_data[[prev_gen_col]])
  n_current <- dplyr::n_distinct(env_data[[curr_gen_col]])
  n_combos <- if (combo_col %in% names(env_data)) dplyr::n_distinct(env_data[[combo_col]]) else NA_integer_
  expected_combos <- n_previous * n_current
  min_reps <- if (length(combo_counts) > 0) min(as.integer(combo_counts), na.rm = TRUE) else NA_integer_

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
    Min_Reps_Per_Combo = min_reps,
    Status = "not_fit",
    Model_Status = NA_character_,
    Residual_Variance = NA_real_,
    SE_Min = NA_real_,
    SE_Median = NA_real_,
    SE_Max = NA_real_,
    Network_Correction_Min = NA_real_,
    Network_Correction_Median = NA_real_,
    Network_Correction_Max = NA_real_
  )
}

.lr_pair_status <- function(diag) {
  status <- character(0)
  if (!isTRUE(diag$Local_Factorial_Complete)) status <- c(status, "incomplete_local_factorial")
  if (isTRUE(diag$Expected_10x10)) status <- c(status, "complete_10x10")
  if (!is.na(diag$Min_Reps_Per_Combo) && diag$Min_Reps_Per_Combo < 2) status <- c(status, "low_replication")
  if (isTRUE(diag$Model_Singular)) status <- c(status, "singular_fit")
  if (!is.na(diag$Combo_Variance) && diag$Combo_Variance <= 0) status <- c(status, "zero_combo_variance")
  if (length(status) == 0) status <- "fit"
  paste(unique(status), collapse = ";")
}

.lr_model_residual_variance <- function(model) {
  if (is.null(model)) return(NA_real_)
  if (inherits(model, "merMod")) return(stats::sigma(model)^2)
  if (!is.null(model$psi)) return(model$psi[1])
  NA_real_
}

.lr_confidence_class <- function(se) {
  out <- rep("unknown", length(se))
  ok <- !is.na(se)
  if (!any(ok)) return(out)
  q <- stats::quantile(se[ok], probs = c(0.33, 0.66), na.rm = TRUE, names = FALSE)
  out[ok & se <= q[[1]]] <- "high"
  out[ok & se > q[[1]] & se <= q[[2]]] <- "moderate"
  out[ok & se > q[[2]]] <- "low"
  out
}

.lr_safe_min <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  min(x, na.rm = TRUE)
}

.lr_safe_median <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  stats::median(x, na.rm = TRUE)
}

.lr_safe_max <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}

.lr_plot_predecessor <- function(df, trait, show_se = FALSE) {
  if (nrow(df) == 0) return(NULL)
  facet_formula <- if ("Baseline_Group" %in% names(df)) {
    Baseline_Group ~ Environment
  } else {
    ~Environment
  }

  plot_df <- df |>
    dplyr::mutate(.Plot_Genotype = stats::reorder(.data[["Previous_Genotype"]], .data[["Legacy_Value"]]))

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[[".Plot_Genotype"]], y = .data[["Legacy_Value"]])) +
    ggplot2::geom_segment(
      ggplot2::aes(xend = .data[[".Plot_Genotype"]], y = 0, yend = .data[["Legacy_Value"]]),
      color = "gray72",
      linewidth = 0.4
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = .data[["Legacy_Value"]] > 0),
      size = 2.2,
      alpha = 0.9
    )

  if (show_se && "SE" %in% names(plot_df)) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data[["Legacy_Value"]] - .data[["SE"]], ymax = .data[["Legacy_Value"]] + .data[["SE"]]),
      width = 0.16,
      color = "gray40",
      alpha = 0.55,
      na.rm = TRUE
    )
  }

  p +
    ggplot2::coord_flip() +
    ggplot2::facet_grid(facet_formula, scales = "free_y", space = "free_y") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::scale_color_manual(values = c("TRUE" = "forestgreen", "FALSE" = "firebrick"), guide = "none") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Average Lentil Predecessor Effects",
      subtitle = paste("Wheat trait:", trait, "| within ENV x Facet deviation"),
      caption = "Legacy_Value is centered within ENV x Facet. Cross-facet rankings compare within-facet deviations, not absolute performance in a fully connected 100 x 100 factorial.",
      x = "Previous lentil genotype",
      y = "Effect on wheat performance vs facet mean"
    )
}

.lr_plot_predecessor_ranked <- function(df, trait, n = NULL, show_se = FALSE) {
  if (nrow(df) == 0) return(NULL)
  df <- .lr_as_legacy_plot_df(df)
  .lr_check_columns(df, c("Previous_Genotype", "Environment", "Legacy_Value"))

  plot_df <- df

  if (!is.null(n) && is.finite(n)) {
    plot_df <- plot_df |>
      dplyr::group_by(.data[["Environment"]]) |>
      dplyr::arrange(abs(.data[["Legacy_Value"]]), .by_group = TRUE) |>
      dplyr::slice_tail(n = n) |>
      dplyr::ungroup()
  }

  plot_df <- plot_df |>
    dplyr::arrange(.data[["Environment"]], .data[["Legacy_Value"]]) |>
    dplyr::mutate(
      .Plot_Key = paste(.data[["Environment"]], .data[["Previous_Genotype"]], sep = "___"),
      .Plot_Key = factor(.data[[".Plot_Key"]], levels = unique(.data[[".Plot_Key"]]))
    )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[[".Plot_Key"]], y = .data[["Legacy_Value"]])) +
    ggplot2::geom_segment(
      ggplot2::aes(xend = .data[[".Plot_Key"]], y = 0, yend = .data[["Legacy_Value"]]),
      color = "gray72",
      linewidth = 0.4
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = .data[["Baseline_Group"]]),
      size = 2.2,
      alpha = 0.9
    )

  if (show_se && "SE" %in% names(plot_df)) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data[["Legacy_Value"]] - .data[["SE"]], ymax = .data[["Legacy_Value"]] + .data[["SE"]]),
      width = 0.16,
      color = "gray40",
      alpha = 0.55,
      na.rm = TRUE
    )
  }

  p +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~Environment, scales = "free_y") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::scale_x_discrete(labels = function(x) sub("^.*___", "", x)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "GWAS-Ready Lentil Legacy Values",
      subtitle = paste("Facet-corrected wheat legacy value", if (!is.null(trait)) paste("| trait:", trait) else ""),
      caption = "Legacy_Value is centered within ENV x Facet. Rankings across all lentils compare within-facet deviations, not absolute performance across a fully connected 100 x 100 factorial.",
      x = "Previous lentil genotype",
      y = "Facet-corrected Legacy_Value",
      color = "Wheat-partner facet"
    )
}

.lr_plot_predecessor_correction <- function(df, trait) {
  if (nrow(df) == 0) return(NULL)
  needed <- c("Raw_Mean", "Corrected_Mean", "Legacy_Value")
  if (!all(needed %in% names(df))) return(NULL)

  facet_formula <- if ("Baseline_Group" %in% names(df)) {
    Baseline_Group ~ Environment
  } else {
    ~Environment
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(y = reorder(Previous_Genotype, Legacy_Value))) +
    ggplot2::geom_segment(
      ggplot2::aes(x = Raw_Mean, xend = Corrected_Mean, yend = Previous_Genotype),
      color = "gray60"
    ) +
    ggplot2::geom_point(ggplot2::aes(x = Raw_Mean), shape = 1, size = 1.8, color = "black") +
    ggplot2::geom_point(
      ggplot2::aes(x = Corrected_Mean, color = Legacy_Value > 0),
      size = 2.4
    )

  if ("Baseline_Mean" %in% names(df)) {
    p <- p + ggplot2::geom_vline(
      ggplot2::aes(xintercept = Baseline_Mean),
      linetype = "dashed",
      color = "gray35",
      alpha = 0.7
    )
  }

  p +
    ggplot2::facet_grid(facet_formula, scales = "free_y", space = "free_y") +
    ggplot2::scale_color_manual(values = c("TRUE" = "forestgreen", "FALSE" = "firebrick"), guide = "none") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Raw vs Model-Corrected Lentil Legacy",
      subtitle = paste("Hollow = raw mean, solid = corrected mean | Wheat trait:", trait),
      x = paste("Wheat", trait),
      y = "Previous lentil genotype"
    )
}

.lr_plot_pair_compatibility <- function(df, trait, n = 25) {
  plot_pair_compatibility_ranked(df, trait = trait, n = n)
}

.lr_summarize_baseline <- function(x) {
  if (length(x) == 0) return("All")
  vals <- unique(as.character(x[!is.na(x)]))
  if (length(vals) == 0) return("Missing")
  vals[1]
}

.lr_count_baseline <- function(x) {
  if (length(x) == 0) return(1L)
  length(unique(as.character(x[!is.na(x)])))
}
