#' Calculate Robust Legacy Values using SpATS
#'
#' Backward-compatible wrapper for `model_predecessor_effect()`. It answers:
#' which previous lentil genotypes leave better average conditions for wheat?
#'
#' @param data The merged dataframe from merge_rotation_data() or wheat rotation data.
#' @param trait String. The cereal trait to analyze, for example `"Y.ADJ"`.
#' @param env_col String. The cereal environment column.
#' @param prev_gen_col String. Name of the previous crop genotype column.
#' @param curr_gen_col String. Name of the current wheat genotype column.
#' @param spatial_cols Character vector. c("Row", "Col") for the cereal field.
#' @param rep_col String (Optional). Additional fixed design term.
#'
#' @return A list containing the legacy values, plot, and fitted models.
#' @export
model_legacy_spats <- function(data, trait, env_col,
                               prev_gen_col = "Previous_Crop_Genotype",
                               curr_gen_col = "Genotype",
                               spatial_cols = c("Row", "Col"),
                               rep_col = NULL) {
  fixed_effect_cols <- character(0)
  if (!is.null(rep_col)) fixed_effect_cols <- rep_col

  res <- model_predecessor_effect(
    data = data,
    trait = trait,
    env_col = env_col,
    prev_gen_col = prev_gen_col,
    curr_gen_col = curr_gen_col,
    spatial_cols = spatial_cols,
    fixed_effect_cols = fixed_effect_cols,
    method = "SpATS",
    include_checks = FALSE
  )

  list(
    legacy_values = res$legacy_values,
    gwas_phenotypes = res$gwas_phenotypes,
    plot = res$plot,
    ranked_plot = res$ranked_plot,
    gwas_plot = res$gwas_plot,
    correction_plot = res$correction_plot,
    models = res$models
  )
}

#' Calculate Robust Legacy Values using SpATS with Raw-vs-Corrected Output
#'
#' Backward-compatible wrapper for `model_predecessor_effect()`. The returned
#' table includes raw means, corrected means, and network correction values.
#'
#' @param data The merged dataframe from merge_rotation_data() or wheat rotation data.
#' @param trait String. The cereal trait to analyze.
#' @param env_col String. The cereal environment column.
#' @param prev_gen_col String. Name of the previous crop genotype column.
#' @param curr_gen_col String. Name of the current wheat genotype column.
#' @param spatial_cols Character vector. c("Row", "Col") for the cereal field.
#' @param rep_col String (Optional). Additional fixed design term.
#'
#' @return A list containing the legacy values, plot, and fitted models.
#' @export
model_legacy_spats2 <- function(data, trait, env_col,
                                prev_gen_col = "Previous_Crop_Genotype",
                                curr_gen_col = "Genotype",
                                spatial_cols = c("Row", "Col"),
                                rep_col = NULL) {
  res <- model_legacy_spats(
    data = data,
    trait = trait,
    env_col = env_col,
    prev_gen_col = prev_gen_col,
    curr_gen_col = curr_gen_col,
    spatial_cols = spatial_cols,
    rep_col = rep_col
  )

  list(
    legacy_values = res$legacy_values,
    gwas_phenotypes = res$gwas_phenotypes,
    plot = res$correction_plot,
    ranked_plot = res$ranked_plot,
    gwas_plot = res$gwas_plot,
    correction_plot = res$correction_plot,
    models = res$models
  )
}
