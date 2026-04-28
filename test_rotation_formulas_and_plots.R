# title: Configure the formula and plot validation sandbox.
# This script is the main place to test every modeling and plotting function in /R against the ACTIVATE lentil and wheat phenotype files.
rm(list = ls())

RUN_SPATS_TRAIT_MODELS <- TRUE
RUN_ROTATION_MODELS <- TRUE
RUN_LEGACY_WRAPPERS <- TRUE
RUN_FACET_MODEL <- TRUE
SAVE_PLOTS <- FALSE
PLOT_DIR <- file.path("Figures", "formula_test_outputs")

if (SAVE_PLOTS && !dir.exists(PLOT_DIR)) dir.create(PLOT_DIR, recursive = TRUE)

required_packages <- c(
  "dplyr", "tidyr", "ggplot2", "lme4", "SpATS", "reshape2", "ggrepel"
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

wheat_traits <- c("HD", "HT", "MAT", "LD", "YLD", "Y.ADJ", "TWT", "KWT", "PRO")
lentil_traits <- c("DTE", "DTF", "VegP", "DTM", "RepP", "lodging", "YLD", "PRT", "DS")
lentil_subsample_traits <- c(
  "biomass.g", "straw.g", "seed.g", "n.seeds", "KSW", "HI", "NIR.seed.pre",
  "LECO.stover.C.pct", "LECO.stover.N.pct", "LECO.stover.P.pct", "LECO.stover.N.g.m",
  "LECO.seed.C.pct", "LECO.seed.N.pct", "LECO.seed.P.pct", "LECO.seed.N.g.m",
  "NHI.pct", "NHI.rel", "C.N.ratio.stover", "C.N.ratio.seed"
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


# title: Inspect trial structure and trait distributions.
# These checks confirm column names, genotype counts, replicate balance, trait ranges, and basic distribution plots before fitting models.
lentil_met <- safe_run("inspect_met: lentil full trial", {
  inspect_met(
    data = lentil_treatment,
    unique_id_col = "Unique.Row.ID",
    env_col = "ENV",
    genotype_col = "Lentil",
    design_cols = c("Rep_gen", "Block", "Row", "Col"),
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
    design_cols = c("Rep_combo", "Facet", "Block", "Row", "Col"),
    trait_cols = wheat_traits,
    crop_col = "Crop"
  )
})

lentil_trial <- safe_run("inspect_trial: lentil single-trial summary", {
  inspect_trial(
    data = dplyr::filter(lentil_treatment, ENV == unique(lentil_treatment$ENV)[1]),
    genotype_col = "Lentil",
    design_cols = c("Rep_gen", "Block", "Row", "Col"),
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
    fixed_model = ~ 1,
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
    fixed_model = ~ 1,
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
      rep_col = "Rep_combo",
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


# title: Merge or expose rotation links.
# This tests merge_rotation_data while using the wheat file's existing Lentil column as the previous-crop genotype link.
rotation_data <- safe_run("merge_rotation_data: use wheat Lentil column", {
  merge_rotation_data(
    cereal_data = wheat_treatment,
    cereal_prev_col = "Lentil"
  )
})


# title: Test the older simple legacy model.
# This legacy formula is useful as a simple baseline, but the newer facet-aware predecessor model below is the preferred interpretation.
simple_legacy <- safe_run("get_legacy_values: simple predecessor model", {
  get_legacy_values(
    data = rotation_data,
    trait = "Y.ADJ",
    env_col = "ENV",
    rep_col = "Block"
  )
})
show_plot(simple_legacy$plot, "simple_legacy_yadj")


# title: Test the preferred average predecessor effect formula.
# This answers which lentil genotypes leave better wheat conditions using a facet-aware ENV x Facet baseline and wheat genotype correction.
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
  print(head(avg_predecessor_yadj$gwas_phenotypes, 20))
  show_plot(avg_predecessor_yadj$ranked_plot, "avg_predecessor_yadj_gwas_ranked")
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
  print(head(avg_predecessor_pro$gwas_phenotypes, 20))
  show_plot(avg_predecessor_pro$ranked_plot, "avg_predecessor_pro_gwas_ranked")
  show_plot(avg_predecessor_pro$plot, "avg_predecessor_pro_facet_baseline")
  show_plot(avg_predecessor_pro$correction_plot, "avg_predecessor_pro_raw_vs_corrected")
}


# title: Export the GWAS-ready predecessor phenotypes.
# This creates one facet-corrected lentil predecessor value per genotype, environment, and trait, which can be used as a Year-2 legacy phenotype for lentil GWAS.
if (RUN_ROTATION_MODELS) {
  predecessor_gwas_phenotypes <- safe_run("extract_predecessor_gwas_phenotypes: Y.ADJ and PRO", {
    dplyr::bind_rows(
      avg_predecessor_yadj$gwas_phenotypes,
      avg_predecessor_pro$gwas_phenotypes
    )
  })

  print(head(predecessor_gwas_phenotypes, 20))
  show_plot(
    plot_predecessor_gwas_phenotypes(avg_predecessor_yadj, trait = "Y.ADJ"),
    "avg_predecessor_yadj_gwas_ranked_helper"
  )
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
      random_effect_cols = c("Block"),
      type_col = "Type",
      include_checks = FALSE
    )
  })

  print(head(pair_compatibility_yadj$compatibility, 20))
  show_plot(pair_compatibility_yadj$plot, "pair_compatibility_yadj")

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
      random_effect_cols = c("Block"),
      type_col = "Type",
      include_checks = FALSE
    )
  })

  print(head(pair_compatibility_pro$compatibility, 20))
  show_plot(pair_compatibility_pro$plot, "pair_compatibility_pro")
}


# title: Test backward-compatible legacy wrappers.
# These functions call the preferred predecessor-effect engine while preserving the older model_legacy_spats and model_legacy_spats2 interfaces.
if (RUN_LEGACY_WRAPPERS) {
  legacy_spats_wrapper <- safe_run("model_legacy_spats wrapper: Y.ADJ", {
    model_legacy_spats(
      data = rotation_data,
      trait = "Y.ADJ",
      env_col = "ENV",
      prev_gen_col = "Previous_Crop_Genotype",
      curr_gen_col = "Wheat",
      spatial_cols = c("Row", "Col"),
      rep_col = "Rep_combo"
    )
  })
  show_plot(legacy_spats_wrapper$plot, "legacy_spats_wrapper_yadj")
  show_plot(legacy_spats_wrapper$correction_plot, "legacy_spats_wrapper_yadj_correction")

  legacy_spats2_wrapper <- safe_run("model_legacy_spats2 wrapper: Y.ADJ correction plot", {
    model_legacy_spats2(
      data = rotation_data,
      trait = "Y.ADJ",
      env_col = "ENV",
      prev_gen_col = "Previous_Crop_Genotype",
      curr_gen_col = "Wheat",
      spatial_cols = c("Row", "Col"),
      rep_col = "Rep_combo"
    )
  })
  show_plot(legacy_spats2_wrapper$plot, "legacy_spats2_wrapper_yadj")
}


# title: Test the facet-corrected legacy formula.
# This older alternative estimates wheat baselines and lentil observed means separately, then compares lentil performance against its specific wheat partners.
if (RUN_FACET_MODEL) {
  legacy_facet <- safe_run("model_legacy_facet: Y.ADJ", {
    model_legacy_facet(
      data = wheat_treatment,
      trait = "Y.ADJ",
      env_col = "ENV",
      prev_gen_col = "Lentil",
      curr_gen_col = "Wheat",
      spatial_cols = c("Row", "Col"),
      rep_col = "Rep_combo"
    )
  })
  show_plot(legacy_facet$plot, "legacy_facet_yadj")
}


# title: Build modular legacy-correlation input.
# This combines full-trial lentil BLUPs, lentil subsample chemistry/biomass traits, and multiple wheat legacy targets so you can ask which Year-1 traits drive each Year-2 legacy phenotype.
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

  wheat_legacy_targets <- prepare_wheat_legacy_targets(
    YADJ = avg_predecessor_yadj,
    PRO = avg_predecessor_pro
  )

  build_legacy_correlation_input(
    predictors = list(lentil_full_trial_predictors, lentil_subsample_predictors),
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
    min_pairs = 10
  )
})
print(head(legacy_driver_pro$correlations, 20))
show_plot(legacy_driver_pro$plot, "legacy_driver_correlations_pro")

legacy_driver_yadj <- safe_run("plot_legacy_driver_correlations: drivers of Wheat_Legacy_YADJ", {
  plot_legacy_driver_correlations(
    data = legacy_correlation_input,
    target_trait = "Wheat_Legacy_YADJ",
    top_n = 20,
    min_pairs = 10
  )
})
print(head(legacy_driver_yadj$correlations, 20))
show_plot(legacy_driver_yadj$plot, "legacy_driver_correlations_yadj")


# title: Collect the most important output objects.
# This list keeps the main fitted results together so you can inspect tables, plots, and models after the script finishes.
formula_test_results <- list(
  lentil_met = lentil_met,
  wheat_met = wheat_met,
  lentil_lme4_models = lentil_lme4_models,
  wheat_lme4_models = wheat_lme4_models,
  avg_predecessor_yadj = if (exists("avg_predecessor_yadj")) avg_predecessor_yadj else NULL,
  avg_predecessor_pro = if (exists("avg_predecessor_pro")) avg_predecessor_pro else NULL,
  predecessor_gwas_phenotypes = if (exists("predecessor_gwas_phenotypes")) predecessor_gwas_phenotypes else NULL,
  pair_compatibility_yadj = if (exists("pair_compatibility_yadj")) pair_compatibility_yadj else NULL,
  pair_compatibility_pro = if (exists("pair_compatibility_pro")) pair_compatibility_pro else NULL,
  legacy_correlation_input = legacy_correlation_input,
  legacy_driver_pro = if (exists("legacy_driver_pro")) legacy_driver_pro else NULL,
  legacy_driver_yadj = if (exists("legacy_driver_yadj")) legacy_driver_yadj else NULL
)

cat("\nFormula and plot testing script finished. Inspect formula_test_results for fitted objects.\n")
