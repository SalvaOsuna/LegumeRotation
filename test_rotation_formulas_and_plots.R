# title: Configure the formula and plot validation sandbox.
# This script is the main place to test every modeling and plotting function in /R against the ACTIVATE lentil and wheat phenotype files.
rm(list = ls())

RUN_SPATS_TRAIT_MODELS <- TRUE
RUN_ROTATION_MODELS <- TRUE
SAVE_PLOTS <- FALSE
PLOT_DIR <- file.path("Figures", "formula_test_outputs")

if (SAVE_PLOTS && !dir.exists(PLOT_DIR)) dir.create(PLOT_DIR, recursive = TRUE)

required_packages <- c(
  "dplyr", "tidyr", "ggplot2", "lme4", "SpATS", "reshape2", "ggrepel", "devtools"
)
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Install these packages before running this script: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

source_all_r <- function(path = "R") {
  files <- list.files(path, pattern = "\\.R$", full.names = TRUE)
  invisible(lapply(files, source))
}

safe_run <- function(label, expr) {
  cat("\n---", label, "---\n")
  tryCatch(
    {
      value <- eval.parent(substitute(expr))
      cat("OK:", label, "\n")
      value
    },
    error = function(e) {
      warning(paste(label, "failed:", e$message))
      NULL
    }
  )
}

show_plot <- function(plot_object, name) {
  if (is.null(plot_object)) return(invisible(NULL))
  print(plot_object)
  if (SAVE_PLOTS) {
    ggplot2::ggsave(
      filename = file.path(PLOT_DIR, paste0(name, ".png")),
      plot = plot_object,
      width = 11,
      height = 8,
      dpi = 300
    )
  }
  invisible(plot_object)
}

source_all_r()


# title: Load and lightly standardize the phenotype data.
# The wheat file already contains the previous lentil genotype and wheat facet, so this section creates analysis-ready treatment-only tables.
wheat_pheno <- read.csv(
  file.path("data-raw", "ACTIVATE_wheat_myYraw_rep_gen.csv"),
  na.strings = c("", "NA")
)

lentil_pheno <- read.csv(
  file.path("data-raw", "ACTIVATE_lentil_myYraw_rep_gen.csv"),
  na.strings = c("", "NA")
)

lentil_subsample <- read.delim(
  file.path("data-raw", "ACTIVATE_lentil_subsample.txt"),
  na.strings = c("", "NA")
)

lentil_microbiome <- read.csv(
  file.path("data-raw", "lentil_microbiome_layer_2024.csv"),
  na.strings = c("", "NA")
)

wheat_traits <- c("HD", "HT", "MAT", "LD", "YLD", "Y.ADJ", "TWT", "KWT", "PRO")
lentil_traits <- c("DTE", "DTF", "VegP", "DTM", "RepP", "lodging", "YLD", "PRT", "DS")
lentil_subsample_traits <- c(
  "biomass.g", "straw.g", "seed.g", "n.seeds", "KSW", "HI",
  "LECO.stover.C.pct",  "LECO.stover.N.g.m",
  "LECO.seed.C.pct",   "LECO.seed.N.g.m",
  "NHI.pct", "NHI.rel", "C.N.ratio.stover", "C.N.ratio.seed"
) # traits not included: "NIR.seed.pre", "LECO.stover.P.pct","LECO.seed.P.pct","LECO.stover.N.pct","LECO.seed.N.pct",
lentil_microbiome_traits <- c(
  "Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson",
  "Pielou_J", "Aitchison_PC1", "Aitchison_PC2", "MicrobiomeLegacyScore",
  "RhizobiaceaeRelAbundance", "BradyrhizobiaceaeRelAbundance",
  "LegumeSymbiontRelAbundance"
)

wheat_treatment <- wheat_pheno |>
  dplyr::filter(Type == "Treatment") |>
  dplyr::mutate(
    Previous_Crop_Genotype = Lentil,
    Genotype = Wheat,
    Unique.Row.ID = paste(ENV, Unique.ID, sep = "_")
  )

lentil_treatment <- lentil_pheno |>
  dplyr::filter(grepl("^L[0-9]{3}$", Lentil)) |>
  dplyr::mutate(
    Genotype = Lentil,
    Unique.Row.ID = Unique.ID
  )

#From the package root, the easiest option to activate the functions is:
devtools::load_all()


# title: Audit the structured incomplete factorial rotation design.
# The expected ACTIVATE structure is globally incomplete but locally complete: each ENV x Facet should contain 10 lentils, 10 wheats, and 100 observed combinations.
rotation_design_audit <- safe_run("audit_rotation_design: local 10x10 facets", {
  audit <- audit_rotation_design(wheat_treatment)
  stopifnot(all(audit$design_summary$N_Previous_Genotypes == 10))
  stopifnot(all(audit$design_summary$N_Current_Genotypes == 10))
  stopifnot(all(audit$design_summary$N_Combos == 100))
  stopifnot(all(audit$design_summary$Expected_10x10))
  stopifnot(all(audit$lentil_assignment$Facet_Assignment_OK))
  stopifnot(all(audit$wheat_assignment$Facet_Assignment_OK))
  audit
})
print(rotation_design_audit$design_summary)

# title: Inspect trial structure and trait distributions.
# These checks confirm column names, genotype counts, replicate balance, trait ranges, and basic distribution plots before fitting models.
lentil_met <- safe_run("inspect_met: lentil full trial", {
  inspect_met(
    data = lentil_treatment,
    unique_id_col = "Unique.Row.ID",
    env_col = "ENV",
    genotype_col = "Lentil",
    design_cols = c("Rep_combo","Rep_gen", "Block", "Row", "Col"),
    trait_cols = lentil_traits,
    crop_col = "Crop"
  )
})

wheat_met <- safe_run("inspect_met: wheat rotation trial", {
  inspect_met(
    data = wheat_treatment,
    unique_id_col = "Unique.Row.ID",
    env_col = "ENV",
    genotype_col = "Wheat",
    design_cols = c("Rep_combo","Rep_gen", "Facet", "Block", "Row", "Col"),
    trait_cols = wheat_traits,
    crop_col = "Crop"
  )
})

lentil_trial <- safe_run("inspect_trial: lentil single-trial summary", {
  inspect_trial(
    data = dplyr::filter(lentil_treatment, ENV == unique(lentil_treatment$ENV)[1]),
    genotype_col = "Lentil",
    design_cols = c("Rep_combo","Rep_gen", "Block", "Row", "Col"),
    trait_cols = lentil_traits,
    crop_col = "Crop"
  )
})

lentil_trait_inspection <- safe_run("inspect_traits: lentil traits", {
  inspect_traits(lentil_treatment, trait_cols = lentil_traits, env_col = "ENV", plot = TRUE)
})
show_plot(lentil_trait_inspection$plot, "inspect_lentil_traits")

wheat_trait_inspection <- safe_run("inspect_traits: wheat traits", {
  inspect_traits(wheat_treatment, trait_cols = wheat_traits, env_col = "ENV", plot = TRUE)
})
show_plot(wheat_trait_inspection$plot, "inspect_wheat_traits")


# title: Fit lentil phenotypic trait models with lme4.
# This tests the generic BLUP and variance-output path for lentil traits without spatial smoothing, then checks the trait-correlation plot.
lentil_lme4_models <- safe_run("model_traits lme4: lentil YLD and PRT", {
  model_traits(
    data = lentil_treatment,
    method = "lme4",
    trait_cols = c("YLD", "PRT"),
    gen_col = "Lentil",
    env_col = "ENV",
    fixed_model = ~ Rep_gen + Block,
    random_model = ~ (1 | Lentil)
  )
})

lentil_trait_correlation_plot <- safe_run("plot_trait_correlations: lentil YLD vs PRT", {
  plot_trait_correlations(
    blup_data = lentil_lme4_models$blups,
    trait_x = "YLD",
    trait_y = "PRT"
  )
})
show_plot(lentil_trait_correlation_plot, "lentil_trait_correlation")

# title: Fit wheat phenotypic trait models with lme4.
# This tests the generic BLUP and variance-output path for wheat traits and checks that wheat trait correlations can be plotted.
wheat_lme4_models <- safe_run("model_traits lme4: wheat Y.ADJ and PRO", {
  model_traits(
    data = wheat_treatment,
    method = "lme4",
    trait_cols = c("Y.ADJ", "PRO"),
    gen_col = "Wheat",
    env_col = "ENV",
    fixed_model = ~ Rep_gen + Block,
    random_model = ~ (1 | Wheat)
  )
})

wheat_trait_correlation_plot <- safe_run("plot_trait_correlations: wheat Y.ADJ vs PRO", {
  plot_trait_correlations(
    blup_data = wheat_lme4_models$blups,
    trait_x = "Y.ADJ",
    trait_y = "PRO"
  )
})
show_plot(wheat_trait_correlation_plot, "wheat_trait_correlation")


# title: Fit spatial SpATS trait models and spatial trend plots.
# This tests the spatial formula in model_traits and the heatmap function that visualizes corrected field heterogeneity.
if (RUN_SPATS_TRAIT_MODELS) {
  lentil_spats_models <- safe_run("model_traits SpATS: lentil YLD", {
    model_traits(
      data = lentil_treatment,
      method = "SpATS",
      trait_cols = c("YLD"),
      gen_col = "Lentil",
      env_col = "ENV",
      rep_col = "Rep_gen",
      spatial_cols = c("Row", "Col")
    )
  })

  lentil_spatial_plot <- safe_run("plot_spatial_trends: lentil YLD", {
    plot_spatial_trends(
      spatial_data = lentil_spats_models$spatial_trends,
      traits = "YLD",
      mode = "percentage"
    )
  })
  show_plot(lentil_spatial_plot, "lentil_spatial_trends_yld")

  wheat_spats_models <- safe_run("model_traits SpATS: wheat Y.ADJ", {
    model_traits(
      data = wheat_treatment,
      method = "SpATS",
      trait_cols = c("Y.ADJ"),
      gen_col = "Wheat",
      env_col = "ENV",
      rep_col = "Rep_gen",
      spatial_cols = c("Row", "Col")
    )
  })

  wheat_spatial_plot <- safe_run("plot_spatial_trends: wheat Y.ADJ", {
    plot_spatial_trends(
      spatial_data = wheat_spats_models$spatial_trends,
      traits = "Y.ADJ",
      mode = "percentage"
    )
  })
  show_plot(wheat_spatial_plot, "wheat_spatial_trends_yadj")
}


# title: Test the preferred average predecessor effect formula.
# This answers which lentil genotypes leave better wheat conditions using a facet-aware ENV x Facet baseline and wheat genotype correction.
summarize_environment_legacy_significance <- function(result,
                                                      alpha = 0.05,
                                                      p_adjust_method = "BH") {
  # Approximate environment-level heterogeneity test over facet-centered Legacy_Value estimates.
  # SpATS returns prediction SEs for corrected means, not full contrast SEs for each Legacy_Value.
  if (is.null(result) || !"legacy_values" %in% names(result)) {
    stop("Expected a model_predecessor_effect() result with a legacy_values table.")
  }

  legacy_values <- result$legacy_values
  if (is.null(legacy_values) || nrow(legacy_values) == 0) {
    stop("No legacy_values rows were returned.")
  }

  se_col <- dplyr::case_when(
    "Legacy_Value_SE" %in% names(legacy_values) &&
      any(!is.na(legacy_values[["Legacy_Value_SE"]])) ~ "Legacy_Value_SE",
    "Corrected_Mean_SE" %in% names(legacy_values) ~ "Corrected_Mean_SE",
    "SE" %in% names(legacy_values) ~ "SE",
    TRUE ~ NA_character_
  )

  if (is.na(se_col)) {
    stop("No SE column was available for approximate legacy-effect significance.")
  }

  legacy_values |>
    dplyr::mutate(
      Legacy_Effect_SE = .data[[se_col]]
    ) |>
    dplyr::group_by(.data[["Environment"]], .data[["Trait"]]) |>
    dplyr::group_modify(~{
      df <- .x |>
        dplyr::filter(!is.na(.data[["Legacy_Value"]]))

      n_estimates <- nrow(df)
      if (n_estimates < 2) {
        return(data.frame(
          N_Estimates = n_estimates,
          N_Facets = dplyr::n_distinct(df[["Baseline_Group"]]),
          N_Previous_Genotypes = dplyr::n_distinct(df[["Previous_Genotype"]]),
          Weighted_Mean_Legacy = NA_real_,
          Mean_Abs_Legacy = NA_real_,
          Legacy_SD = NA_real_,
          Max_Abs_Legacy = NA_real_,
          Q_Statistic = NA_real_,
          DF = NA_integer_,
          P_Value = NA_real_,
          Weighting = "insufficient_data"
        ))
      }

      has_se <- is.finite(df[["Legacy_Effect_SE"]]) & df[["Legacy_Effect_SE"]] > 0
      if (all(has_se)) {
        weights <- 1 / (df[["Legacy_Effect_SE"]]^2)
        weighting <- "inverse_SE2"
      } else {
        weights <- rep(1, n_estimates)
        weighting <- "equal"
      }

      weighted_mean <- stats::weighted.mean(df[["Legacy_Value"]], weights, na.rm = TRUE)
      q_stat <- sum(weights * (df[["Legacy_Value"]] - weighted_mean)^2, na.rm = TRUE)
      q_df <- n_estimates - 1L

      data.frame(
        N_Estimates = n_estimates,
        N_Facets = dplyr::n_distinct(df[["Baseline_Group"]]),
        N_Previous_Genotypes = dplyr::n_distinct(df[["Previous_Genotype"]]),
        Weighted_Mean_Legacy = weighted_mean,
        Mean_Abs_Legacy = mean(abs(df[["Legacy_Value"]]), na.rm = TRUE),
        Legacy_SD = stats::sd(df[["Legacy_Value"]], na.rm = TRUE),
        Max_Abs_Legacy = max(abs(df[["Legacy_Value"]]), na.rm = TRUE),
        Q_Statistic = q_stat,
        DF = q_df,
        P_Value = stats::pchisq(q_stat, df = q_df, lower.tail = FALSE),
        Weighting = weighting
      )
    }) |>
    dplyr::ungroup() |>
    dplyr::group_by(.data[["Trait"]]) |>
    dplyr::mutate(
      FDR = stats::p.adjust(.data[["P_Value"]], method = p_adjust_method),
      Significant = !is.na(.data[["FDR"]]) & .data[["FDR"]] < alpha
    ) |>
    dplyr::ungroup() |>
    dplyr::arrange(.data[["Trait"]], .data[["FDR"]], .data[["Environment"]])
}

if (RUN_ROTATION_MODELS) {
  avg_predecessor_yadj <- safe_run("model_predecessor_effect SpATS: Y.ADJ", {
    model_predecessor_effect(
      data = wheat_treatment,
      trait = "Y.ADJ",
      env_col = "ENV",
      prev_gen_col = "Lentil",
      curr_gen_col = "Wheat",
      spatial_cols = c("Row", "Col"),
      baseline_col = "Facet",
      fixed_effect_cols = c("Rep_combo"),
      type_col = "Type",
      include_checks = FALSE,
      method = "SpATS"
    )
  })

  print(head(avg_predecessor_yadj$legacy_values, 20))
  legacy_environment_significance_yadj <- safe_run("environment-level legacy effect significance: Y.ADJ", {
    summarize_environment_legacy_significance(avg_predecessor_yadj)
  })
  print(legacy_environment_significance_yadj)
  show_plot(avg_predecessor_yadj$ranked_plot, "avg_predecessor_yadj_legacy_ranked")
  show_plot(avg_predecessor_yadj$plot, "avg_predecessor_yadj_facet_baseline")
  show_plot(avg_predecessor_yadj$correction_plot, "avg_predecessor_yadj_raw_vs_corrected")

  avg_predecessor_pro <- safe_run("model_predecessor_effect SpATS: PRO", {
    model_predecessor_effect(
      data = wheat_treatment,
      trait = "PRO",
      env_col = "ENV",
      prev_gen_col = "Lentil",
      curr_gen_col = "Wheat",
      spatial_cols = c("Row", "Col"),
      baseline_col = "Facet",
      fixed_effect_cols = c("Rep_combo"),
      type_col = "Type",
      include_checks = FALSE,
      method = "SpATS"
    )
  })

  print(head(avg_predecessor_pro$legacy_values, 20))
  legacy_environment_significance_pro <- safe_run("environment-level legacy effect significance: PRO", {
    summarize_environment_legacy_significance(avg_predecessor_pro)
  })
  print(legacy_environment_significance_pro)
  show_plot(avg_predecessor_pro$ranked_plot, "avg_predecessor_pro_legacy_ranked")
  show_plot(avg_predecessor_pro$plot, "avg_predecessor_pro_facet_baseline")
  show_plot(avg_predecessor_pro$correction_plot, "avg_predecessor_pro_raw_vs_corrected")
}

# title: Export GWAS-ready legacy values.
# This creates one facet-corrected Legacy_Value per genotype, environment, and trait, which can be used as a Year-2 legacy phenotype for lentil GWAS.
if (RUN_ROTATION_MODELS) {
  legacy_values_for_gwas <- safe_run("legacy_values for GWAS: Y.ADJ and PRO", {
    dplyr::bind_rows(
      avg_predecessor_yadj$legacy_values,
      avg_predecessor_pro$legacy_values
    ) |>
      dplyr::select(
        Genotype = Previous_Genotype,
        Environment,
        Trait,
        Legacy_Value,
        Corrected_Mean,
        Baseline_Group,
        Baseline_Mean,
        Legacy_Value_Global,
        SE,
        N_Plots,
        N_Wheat_Partners
      )
  })

  print(head(legacy_values_for_gwas, 20))
}

# title: Test the pair-specific compatibility formula.
# This answers which observed lentil-wheat pairs perform better or worse than expected from additive lentil and wheat effects within each ENV x Facet network.
if (RUN_ROTATION_MODELS) {
  pair_compatibility_yadj <- safe_run("model_pair_compatibility lme4: Y.ADJ", {
    model_pair_compatibility(
      data = wheat_treatment,
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
      include_checks = FALSE
    )
  })

  print(head(pair_compatibility_yadj$compatibility, 20))
  print(pair_compatibility_yadj$diagnostics)
  show_plot(pair_compatibility_yadj$heatmap, "pair_compatibility_yadj_heatmap")
  show_plot(pair_compatibility_yadj$ranked_plot, "pair_compatibility_yadj_ranked")

  pair_yadj_observed_check <- safe_run("pair compatibility Y.ADJ uses only observed combinations", {
    observed <- wheat_treatment |>
      dplyr::distinct(Environment = ENV, Baseline_Group = Facet, Previous_Genotype = Lentil, Current_Genotype = Wheat, Combo)
    output_combos <- pair_compatibility_yadj$compatibility |>
      dplyr::distinct(Environment, Baseline_Group, Previous_Genotype, Current_Genotype, Combo)
    unobserved <- dplyr::anti_join(
      output_combos,
      observed,
      by = c("Environment", "Baseline_Group", "Previous_Genotype", "Current_Genotype", "Combo")
    )
    stopifnot(nrow(unobserved) == 0)
    output_combos
  })

  pair_compatibility_pro <- safe_run("model_pair_compatibility lme4: PRO", {
    model_pair_compatibility(
      data = wheat_treatment,
      trait = "PRO",
      env_col = "ENV",
      prev_gen_col = "Lentil",
      curr_gen_col = "Wheat",
      combo_col = "Combo",
      spatial_cols = c("Row", "Col"),
      baseline_col = "Facet",
      fixed_effect_cols = c("Rep_combo"),
      random_effect_cols = c("Block"),
      type_col = "Type",
      include_checks = FALSE
    )
  })

  print(head(pair_compatibility_pro$compatibility, 20))
  print(pair_compatibility_pro$diagnostics)
  show_plot(pair_compatibility_pro$heatmap, "pair_compatibility_pro_heatmap")
  show_plot(pair_compatibility_pro$ranked_plot, "pair_compatibility_pro_ranked")
}


# title: Build modular legacy-correlation input.
# This combines full-trial lentil BLUPs, lentil subsample chemistry/biomass traits, microbiome traits, and multiple wheat legacy targets so you can ask which Year-1 traits drive each Year-2 legacy phenotype.
legacy_correlation_input <- safe_run("Build correlation input: lentil predictors plus wheat legacy targets", {
  lentil_full_trial_predictors <- prepare_lentil_blup_predictors(
    blups = lentil_lme4_models$blups,
    from_year = "2024",
    to_year = "2025",
    trait_prefix = "FullTrial_"
  )

  lentil_subsample_predictors <- prepare_lentil_subsample_predictors(
    data = lentil_subsample,
    trait_cols = lentil_subsample_traits,
    env_col = "ENV",
    gen_col = "Lentil",
    type_col = "Type",
    check_values = "Check",
    include_checks = FALSE,
    from_year = "2024",
    to_year = "2025",
    filter_from_year = TRUE,
    trait_prefix = "Subsample_"
  )

  lentil_microbiome_predictors <- prepare_lentil_microbiome_predictors(
    data = lentil_microbiome,
    trait_cols = lentil_microbiome_traits,
    env_col = "ENV",
    gen_col = "Lentil",
    from_year = "2024",
    to_year = "2025",
    filter_from_year = TRUE,
    trait_prefix = "Microbiome_"
  )

  wheat_legacy_targets <- prepare_wheat_legacy_targets(
    YADJ = avg_predecessor_yadj,
    PRO = avg_predecessor_pro
  )

  build_legacy_correlation_input(
    predictors = list(
      lentil_full_trial_predictors,
      lentil_subsample_predictors,
      lentil_microbiome_predictors
    ),
    targets = wheat_legacy_targets
  )
})


# title: Test the legacy-correlation heatmap for wheat yield legacy.
# This gives a broad trait-by-trait view of Year-1 lentil predictors and the Year-2 wheat yield legacy phenotype.
legacy_correlation_plot_yadj <- safe_run("plot_legacy_correlations: drivers of Wheat_Legacy_YADJ", {
  plot_legacy_correlations(
    data = legacy_correlation_input,
    target_trait = "Wheat_Legacy_YADJ"
  )
})
show_plot(legacy_correlation_plot_yadj, "legacy_correlation_heatmap_yadj")


# title: Test the legacy-correlation heatmap for wheat protein legacy.
# This repeats the same modular correlation view for Wheat_Legacy_PRO so you can compare whether protein legacy has different lentil drivers.
legacy_correlation_plot_pro <- safe_run("plot_legacy_correlations: drivers of Wheat_Legacy_PRO", {
  plot_legacy_correlations(
    data = legacy_correlation_input,
    target_trait = "Wheat_Legacy_PRO"
  )
})
show_plot(legacy_correlation_plot_pro, "legacy_correlation_heatmap_pro")


# title: Rank candidate lentil drivers of wheat legacy.
# This focused plot is easier than the full heatmap when many subsample traits are included; it reports correlations, R2, p-values, and FDR for the selected wheat legacy target.
legacy_driver_pro <- safe_run("plot_legacy_driver_correlations: drivers of Wheat_Legacy_PRO", {
  plot_legacy_driver_correlations(
    data = legacy_correlation_input,
    target_trait = "Wheat_Legacy_PRO",
    top_n = 20,
    min_pairs = 30
  )
})
print(head(legacy_driver_pro$correlations, 20))
show_plot(legacy_driver_pro$plot, "legacy_driver_correlations_pro")

legacy_driver_yadj <- safe_run("plot_legacy_driver_correlations: drivers of Wheat_Legacy_YADJ", {
  plot_legacy_driver_correlations(
    data = legacy_correlation_input,
    target_trait = "Wheat_Legacy_YADJ",
    top_n = 20,
    min_pairs = 30
  )
})
print(head(legacy_driver_yadj$correlations, 20))
show_plot(legacy_driver_yadj$plot, "legacy_driver_correlations_yadj")


# title: Inspect individual legacy-vs-trait scatter plots.
# Use this helper to check candidate legacy drivers one trait at a time. Each panel fits a linear model within environment and annotates R2 and p-value.
plot_legacy_trait_scatter <- function(data,
                                      target_trait,
                                      predictor_trait,
                                      min_pairs = 3) {
  required_traits <- c(target_trait, predictor_trait)
  missing_traits <- setdiff(required_traits, unique(data$Trait))
  if (length(missing_traits) > 0) {
    stop("Missing trait(s) in legacy_correlation_input: ", paste(missing_traits, collapse = ", "))
  }

  plot_df <- data |>
    dplyr::filter(.data[["Trait"]] %in% required_traits) |>
    dplyr::select(dplyr::all_of(c("Genotype", "Environment", "Trait", "Predicted"))) |>
    tidyr::pivot_wider(
      names_from = "Trait",
      values_from = "Predicted",
      values_fn = list(Predicted = mean)
    ) |>
    dplyr::filter(!is.na(.data[[target_trait]]), !is.na(.data[[predictor_trait]]))

  if (nrow(plot_df) == 0) {
    stop("No complete genotype x environment pairs for ", predictor_trait, " vs ", target_trait)
  }

  stats_df <- plot_df |>
    dplyr::group_by(.data[["Environment"]]) |>
    dplyr::group_modify(~{
      n_pairs <- nrow(.x)
      can_fit <- n_pairs >= min_pairs &&
        length(unique(.x[[target_trait]])) > 1 &&
        length(unique(.x[[predictor_trait]])) > 1

      if (!can_fit) {
        return(data.frame(N_Pairs = n_pairs, R2 = NA_real_, P_Value = NA_real_))
      }

      fit <- stats::lm(stats::reformulate(predictor_trait, target_trait), data = .x)
      fit_summary <- summary(fit)
      data.frame(
        N_Pairs = n_pairs,
        R2 = fit_summary$r.squared,
        P_Value = stats::coef(fit_summary)[2, "Pr(>|t|)"]
      )
    }) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      Label = paste0(
        "R2 = ", ifelse(is.na(.data[["R2"]]), "NA", sprintf("%.3f", .data[["R2"]])),
        "\np = ", ifelse(is.na(.data[["P_Value"]]), "NA", format.pval(.data[["P_Value"]], digits = 3)),
        "\nN = ", .data[["N_Pairs"]]
      )
    )

  ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[[predictor_trait]], y = .data[[target_trait]])) +
    ggplot2::geom_point(alpha = 0.75, size = 2) +
    ggplot2::geom_smooth(method = "lm", se = TRUE, color = "gray30", fill = "gray80", linewidth = 0.7) +
    ggplot2::geom_text(
      data = stats_df,
      ggplot2::aes(x = -Inf, y = Inf, label = .data[["Label"]]),
      inherit.aes = FALSE,
      hjust = -0.05,
      vjust = 1.1,
      size = 3.2
    ) +
    ggplot2::facet_wrap(~Environment, scales = "free") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Legacy Effect vs Candidate Trait",
      subtitle = paste(target_trait, "vs", predictor_trait),
      x = predictor_trait,
      y = target_trait
    )
}

legacy_scatter_yadj_fulltrial_yld <- safe_run("scatter: Wheat_Legacy_YADJ vs FullTrial_YLD", {
  plot_legacy_trait_scatter(
    data = legacy_correlation_input,
    target_trait = "Wheat_Legacy_YADJ",
    predictor_trait = "FullTrial_YLD"
  )
})
show_plot(legacy_scatter_yadj_fulltrial_yld, "legacy_scatter_yadj_fulltrial_yld")

legacy_scatter_pro_subsample_biomass <- safe_run("scatter: Wheat_Legacy_PRO vs Subsample_biomass.g", {
  plot_legacy_trait_scatter(
    data = legacy_correlation_input,
    target_trait = "Wheat_Legacy_PRO",
    predictor_trait = "Subsample_biomass.g"
  )
})
show_plot(legacy_scatter_pro_subsample_biomass, "legacy_scatter_pro_subsample_biomass")


# title: Collect the most important output objects.
# This list keeps the main fitted results together so you can inspect tables, plots, and models after the script finishes.
formula_test_results <- list(
  lentil_met = lentil_met,
  wheat_met = wheat_met,
  lentil_lme4_models = lentil_lme4_models,
  wheat_lme4_models = wheat_lme4_models,
  avg_predecessor_yadj = if (exists("avg_predecessor_yadj")) avg_predecessor_yadj else NULL,
  avg_predecessor_pro = if (exists("avg_predecessor_pro")) avg_predecessor_pro else NULL,
  legacy_environment_significance_yadj = if (exists("legacy_environment_significance_yadj")) legacy_environment_significance_yadj else NULL,
  legacy_environment_significance_pro = if (exists("legacy_environment_significance_pro")) legacy_environment_significance_pro else NULL,
  legacy_values_for_gwas = if (exists("legacy_values_for_gwas")) legacy_values_for_gwas else NULL,
  pair_compatibility_yadj = if (exists("pair_compatibility_yadj")) pair_compatibility_yadj else NULL,
  pair_yadj_observed_check = if (exists("pair_yadj_observed_check")) pair_yadj_observed_check else NULL,
  pair_compatibility_pro = if (exists("pair_compatibility_pro")) pair_compatibility_pro else NULL,
  legacy_correlation_input = legacy_correlation_input,
  lentil_microbiome_predictors = if (exists("lentil_microbiome_predictors")) lentil_microbiome_predictors else NULL,
  legacy_driver_pro = if (exists("legacy_driver_pro")) legacy_driver_pro else NULL,
  legacy_driver_yadj = if (exists("legacy_driver_yadj")) legacy_driver_yadj else NULL,
  legacy_scatter_yadj_fulltrial_yld = if (exists("legacy_scatter_yadj_fulltrial_yld")) legacy_scatter_yadj_fulltrial_yld else NULL,
  legacy_scatter_pro_subsample_biomass = if (exists("legacy_scatter_pro_subsample_biomass")) legacy_scatter_pro_subsample_biomass else NULL
)

cat("\nFormula and plot testing script finished. Inspect formula_test_results for fitted objects.\n")
