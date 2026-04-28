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
#' @param spatial_cols Column names for row and column position.
#' @param baseline_col Network/facet column used to define the comparison
#'   baseline for legacy values.
#' @param fixed_effect_cols Optional additional fixed design/network terms.
#' @param type_col Optional column used to remove checks.
#' @param check_values Values in `type_col` to remove when `include_checks = FALSE`.
#' @param include_checks Logical. Include check plots in the model?
#' @param method `"SpATS"` or `"lme4"`.
#'
#' @return A list with legacy values, plots, fitted models, and settings.
#' @export
model_predecessor_effect <- function(data,
                                     trait = "Y.ADJ",
                                     env_col = "ENV",
                                     prev_gen_col = "Lentil",
                                     curr_gen_col = "Wheat",
                                     spatial_cols = c("Row", "Col"),
                                     baseline_col = "Facet",
                                     fixed_effect_cols = c("Rep_combo"),
                                     type_col = "Type",
                                     check_values = "Check",
                                     include_checks = FALSE,
                                     method = c("SpATS", "lme4")) {

  method <- match.arg(method)
  data <- .lr_prepare_rotation_data(
    data = data,
    trait = trait,
    env_col = env_col,
    prev_gen_col = prev_gen_col,
    curr_gen_col = curr_gen_col,
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
  envs <- unique(data[[env_col]])

  for (env in envs) {
    env_data <- droplevels(data[data[[env_col]] == env & !is.na(data[[trait]]), ])
    if (nrow(env_data) == 0) next

    fit <- tryCatch({
      if (method == "SpATS") {
        .lr_fit_predecessor_spats(
          env_data = env_data,
          trait = trait,
          prev_gen_col = prev_gen_col,
          curr_gen_col = curr_gen_col,
          spatial_cols = spatial_cols,
          baseline_col = baseline_col,
          fixed_effect_cols = fixed_effect_cols
        )
      } else {
        .lr_fit_predecessor_lme4(
          env_data = env_data,
          trait = trait,
          prev_gen_col = prev_gen_col,
          curr_gen_col = curr_gen_col,
          spatial_cols = spatial_cols,
          baseline_col = baseline_col,
          fixed_effect_cols = fixed_effect_cols
        )
      }
    }, error = function(e) {
      if (method == "SpATS" && requireNamespace("lme4", quietly = TRUE)) {
        warning(paste(
          "SpATS predecessor model failed for", env,
          "and will fall back to lme4:",
          e$message
        ))
        tryCatch({
          .lr_fit_predecessor_lme4(
            env_data = env_data,
            trait = trait,
            prev_gen_col = prev_gen_col,
            curr_gen_col = curr_gen_col,
            spatial_cols = spatial_cols,
            baseline_col = baseline_col,
            fixed_effect_cols = fixed_effect_cols
          )
        }, error = function(e2) {
          warning(paste("Predecessor model failed for", env, ":", e2$message))
          NULL
        })
      } else {
        warning(paste("Predecessor model failed for", env, ":", e$message))
        NULL
      }
    })

    if (is.null(fit)) next

    legacy_df <- fit$legacy_values
    legacy_df$Environment <- env
    legacy_df$Trait <- trait
    legacy_df$Method <- method
    legacy_df <- legacy_df[order(legacy_df$Legacy_Value, decreasing = TRUE), ]
    legacy_df$Rank <- seq_len(nrow(legacy_df))

    results_list[[as.character(env)]] <- legacy_df
    model_list[[as.character(env)]] <- fit$model
  }

  final_df <- dplyr::bind_rows(results_list)

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
    models = model_list,
    settings = list(
      trait = trait,
      method = method,
      include_checks = include_checks,
      baseline_col = baseline_col,
      question = "Which lentil genotypes leave better average wheat conditions?"
    )
  )
}

#' Extract GWAS-Ready Lentil Predecessor Phenotypes
#'
#' Converts predecessor-effect model output into one phenotype per previous
#' lentil genotype, environment, and trait. The default phenotype is
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
      Predecessor_Phenotype = .lr_weighted_mean(.data[[phenotype_col]], .data[["N_Plots"]]),
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
      Rank = rank(-.data[["Predecessor_Phenotype"]], ties.method = "first")
    ) |>
    dplyr::ungroup() |>
    dplyr::arrange(.data[["Environment"]], .data[["Trait"]], .data[["Rank"]])

  out
}

#' Plot GWAS-Ready Lentil Predecessor Phenotypes
#'
#' Shows a single ranked list of lentil genotypes per environment using the
#' facet-corrected predecessor phenotype. Bars are colored by the wheat-partner
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
#' @param fixed_effect_cols Optional additional fixed terms. Avoid terms that
#'   are redundant with lentil or wheat genotype.
#' @param random_effect_cols Optional design random effects, for example `"Block"`.
#' @param type_col Optional column used to remove checks.
#' @param check_values Values in `type_col` to remove when `include_checks = FALSE`.
#' @param include_checks Logical. Include check plots in the model?
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
                                     fixed_effect_cols = character(0),
                                     random_effect_cols = c("Block"),
                                     type_col = "Type",
                                     check_values = "Check",
                                     include_checks = FALSE) {

  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Package 'lme4' is required for pair compatibility models.")
  }

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
  envs <- unique(data[[env_col]])

  for (env in envs) {
    env_data <- droplevels(data[data[[env_col]] == env & !is.na(data[[trait]]), ])
    if (nrow(env_data) == 0) next

    baseline_groups <- unique(as.character(env_data$.Baseline_Group))

    for (baseline_group in baseline_groups) {
      group_data <- droplevels(env_data[as.character(env_data$.Baseline_Group) == baseline_group, ])
      if (nrow(group_data) == 0) next
      if (dplyr::n_distinct(group_data[[prev_gen_col]]) < 2 ||
          dplyr::n_distinct(group_data[[curr_gen_col]]) < 2 ||
          dplyr::n_distinct(group_data[[combo_col]]) < 2) {
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
          random_effect_cols = random_effect_cols
        )
      }, error = function(e) {
        warning(paste("Pair compatibility model failed for", env, baseline_group, ":", e$message))
        NULL
      })

      if (is.null(fit)) next

      pair_df <- fit$compatibility
      pair_df$Environment <- env
      pair_df$Trait <- trait
      pair_df <- pair_df[order(pair_df$Compatibility_Value, decreasing = TRUE), ]
      pair_df$Rank <- seq_len(nrow(pair_df))

      result_key <- paste(env, baseline_group, sep = "::")
      results_list[[result_key]] <- pair_df
      model_list[[result_key]] <- fit$model
    }
  }

  final_df <- dplyr::bind_rows(results_list)
  p <- .lr_plot_pair_compatibility(final_df, trait)

  list(
    compatibility = final_df,
    plot = p,
    models = model_list,
    settings = list(
      trait = trait,
      include_checks = include_checks,
      baseline_col = baseline_col,
      question = "Which lentil-wheat genotype pairs are specifically compatible?"
    )
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

.lr_as_legacy_plot_df <- function(df) {
  out <- df
  if (!"Previous_Genotype" %in% names(out) && "Genotype" %in% names(out)) {
    out$Previous_Genotype <- out$Genotype
  }
  if (!"Legacy_Value" %in% names(out) && "Predecessor_Phenotype" %in% names(out)) {
    out$Legacy_Value <- out$Predecessor_Phenotype
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
      Baseline_Group = .lr_summarize_baseline(.data[[".Baseline_Group"]]),
      N_Baseline_Groups = .lr_count_baseline(.data[[".Baseline_Group"]]),
      .groups = "drop"
    )
  colnames(raw_stats)[1] <- "Previous_Genotype"

  grand_mean <- mean(legacy_df$Corrected_Mean, na.rm = TRUE)
  legacy_df <- legacy_df |>
    dplyr::left_join(raw_stats, by = "Previous_Genotype") |>
    dplyr::group_by(Baseline_Group) |>
    dplyr::mutate(
      Baseline_Mean = mean(Corrected_Mean, na.rm = TRUE),
      Legacy_Value_Global = Corrected_Mean - grand_mean,
      Legacy_Pct_Global = (Legacy_Value_Global / grand_mean) * 100,
      Legacy_Value = Corrected_Mean - Baseline_Mean,
      Legacy_Pct = (Legacy_Value / Baseline_Mean) * 100,
      Network_Correction = Corrected_Mean - Raw_Mean
    ) |>
    dplyr::ungroup()

  list(legacy_values = legacy_df, model = m, fixed_terms = used_fixed_terms)
}

.lr_fit_predecessor_lme4 <- function(env_data,
                                     trait,
                                     prev_gen_col,
                                     curr_gen_col,
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
      Baseline_Group = .lr_summarize_baseline(.data[[".Baseline_Group"]]),
      N_Baseline_Groups = .lr_count_baseline(.data[[".Baseline_Group"]]),
      .groups = "drop"
    )
  colnames(legacy_df)[1] <- "Previous_Genotype"

  grand_mean <- mean(legacy_df$Corrected_Mean, na.rm = TRUE)
  legacy_df <- legacy_df |>
    dplyr::group_by(Baseline_Group) |>
    dplyr::mutate(
      SE = NA_real_,
      Baseline_Mean = mean(Corrected_Mean, na.rm = TRUE),
      Legacy_Value_Global = Corrected_Mean - grand_mean,
      Legacy_Pct_Global = (Legacy_Value_Global / grand_mean) * 100,
      Legacy_Value = Corrected_Mean - Baseline_Mean,
      Legacy_Pct = (Legacy_Value / Baseline_Mean) * 100,
      Network_Correction = Corrected_Mean - Raw_Mean
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
                              random_effect_cols) {
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
  env_data$Predicted_With_Pair <- predict(m)

  pair_stats <- env_data |>
    dplyr::group_by(.data[[combo_col]]) |>
    dplyr::summarize(
      Previous_Genotype = dplyr::first(.data[[prev_gen_col]]),
      Current_Genotype = dplyr::first(.data[[curr_gen_col]]),
      Baseline_Group = .lr_summarize_baseline(.data[[".Baseline_Group"]]),
      Raw_Mean = mean(.data[[trait]], na.rm = TRUE),
      Expected_Additive_Mean = mean(Expected_Additive, na.rm = TRUE),
      Corrected_Pair_Mean = mean(Predicted_With_Pair, na.rm = TRUE),
      N_Plots = dplyr::n(),
      .groups = "drop"
    )
  colnames(pair_stats)[1] <- "Combo"

  compatibility <- pair_stats |>
    dplyr::left_join(combo_effects, by = "Combo") |>
    dplyr::mutate(
      Compatibility_Pct = (Compatibility_Value / Expected_Additive_Mean) * 100
    )

  list(compatibility = compatibility, model = m)
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
      subtitle = paste("Wheat trait:", trait, "| deviation from ENV x wheat-facet mean"),
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
      title = "GWAS-Ready Lentil Predecessor Phenotypes",
      subtitle = paste("Facet-corrected wheat response phenotype", if (!is.null(trait)) paste("| trait:", trait) else ""),
      x = "Previous lentil genotype",
      y = "Facet-corrected predecessor phenotype",
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
  if (nrow(df) == 0) return(NULL)

  group_cols <- "Environment"
  if ("Baseline_Group" %in% names(df)) group_cols <- c("Environment", "Baseline_Group")

  plot_df <- df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::slice_max(abs(Compatibility_Value), n = n, with_ties = FALSE) |>
    dplyr::ungroup()

  facet_formula <- if ("Baseline_Group" %in% names(plot_df)) {
    Baseline_Group ~ Environment
  } else {
    ~Environment
  }

  ggplot2::ggplot(plot_df, ggplot2::aes(x = reorder(Combo, Compatibility_Value), y = Compatibility_Value)) +
    ggplot2::geom_col(ggplot2::aes(fill = Compatibility_Value > 0), alpha = 0.85) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = Compatibility_Value - SE, ymax = Compatibility_Value + SE),
      width = 0.2,
      color = "gray40",
      na.rm = TRUE
    ) +
    ggplot2::coord_flip() +
    ggplot2::facet_grid(facet_formula, scales = "free_y", space = "free_y") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::scale_fill_manual(values = c("TRUE" = "forestgreen", "FALSE" = "firebrick"), guide = "none") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Lentil-Wheat Pair Compatibility",
      subtitle = paste("Largest absolute compatibility values within ENV x wheat facet for trait:", trait),
      x = "Lentil-wheat pair",
      y = "Pair effect beyond additive lentil and wheat effects"
    )
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
