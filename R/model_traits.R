#' Model Phenotypic Traits (SpATS or lme4)
#'
#' @param data Clean data.frame from inspect_met().
#' @param trait_cols Character vector of traits to analyze.
#' @param method String. "SpATS" or "lme4".
#' @param gen_col String. Genotype column.
#' @param env_col String. Environment column.
#' @param rep_col String (Optional). Replicate column.
#' @param spatial_cols Character vector (SpATS only). c("Row", "Col").
#' @param check_col String (Optional).
#' @param check_name String (Optional).
#' @param fixed_model Formula (lme4 only).
#' @param random_model Formula (lme4 only).
#'
#' @return A list containing blups, vars, and spatial_trends (if SpATS).
#' @export
#' @importFrom SpATS SpATS predict.SpATS getHeritability PSANOVA obtain.spatialtrend
#' @importFrom lme4 lmer Ranef VarCorr
#' @importFrom dplyr bind_rows mutate select
model_traits <- function(data, method = "SpATS", trait_cols, gen_col, env_col, rep_col = NULL,
                         spatial_cols = c("Row", "Col"),
                         check_col = NULL, check_name = NULL,
                         fixed_model = NULL, random_model = NULL) {

  results_blup <- list()
  results_var <- list()
  results_spatial <- list() # New storage list

  envs <- unique(data[[env_col]])

  for (env in envs) {
    env_data <- data[data[[env_col]] == env, ]

    for (trait in trait_cols) {
      if (all(is.na(env_data[[trait]]))) next

      # ============================
      # METHOD: SpATS
      # ============================
      if (method == "SpATS") {
        tryCatch({
          # 1. Prepare Data Types
          env_data[[gen_col]] <- as.factor(env_data[[gen_col]])

          # Force Numeric Coordinates
          row_col_name <- spatial_cols[1]
          col_col_name <- spatial_cols[2]
          env_data$Row_Num <- as.numeric(as.character(env_data[[row_col_name]]))
          env_data$Col_Num <- as.numeric(as.character(env_data[[col_col_name]]))

          if (!is.null(rep_col)) env_data[[rep_col]] <- as.factor(env_data[[rep_col]])

          # 2. Define Fixed Formula
          fixed_form <- "~ 1"
          if (!is.null(rep_col)) fixed_form <- paste("~", rep_col)

          # 3. Fit Model
          m <- SpATS::SpATS(
            response = trait,
            genotype = gen_col,
            genotype.as.random = TRUE,
            spatial = ~ SpATS::PSANOVA(Col_Num, Row_Num, nseg = c(10, 10), nest.div = 2),
            fixed = as.formula(fixed_form),
            data = env_data,
            control = list(tolerance = 1e-03, monitoring = 0)
          )

          # 4. Extract Standard Results
          h2 <- SpATS::getHeritability(m)

          preds <- predict(m, which = gen_col)
          preds_df <- preds[, c(gen_col, "predicted.values", "standard.errors")]
          colnames(preds_df) <- c("Genotype", "Predicted", "SE")

          results_blup[[paste(env, trait)]] <- preds_df |>
            dplyr::mutate(Environment = env, Trait = trait, Method = "SpATS")

          results_var[[paste(env, trait)]] <- data.frame(
            Environment = env, Trait = trait, H2_Generalized = h2,
            Var_Gen = m$var.comp[gen_col], Var_Res = m$psi[1]
          )

          # 5. NEW: Extract Spatial Trend
          # obtain.spatialtrend returns a list with grid coordinates and the fitted matrix
          st <- SpATS::obtain.spatialtrend(m)

          # Flatten the matrix into a dataframe for ggplot
          # The matrix 'st$fit' has rows = st$row.p and cols = st$col.p
          st_df <- expand.grid(Row = st$row.p, Col = st$col.p)
          st_df$Trend <- as.vector(st$fit) # Flattens by column

          results_spatial[[paste(env, trait)]] <- st_df |>
            dplyr::mutate(Environment = env, Trait = trait)

          cat(paste("✔ Solved:", env, trait, "| H2:", round(h2, 2), "\n"))

        }, error = function(e) {
          warning(paste("✖ SpATS Failed for", env, trait, ":", e$message))
        })
      }

      # ============================
      # METHOD: lme4
      # ============================
      if (method == "lme4") {
        if (is.null(fixed_model) || is.null(random_model)) stop("For lme4, provide formulas.")

        full_formula <- paste(trait, paste(as.character(fixed_model), collapse = " "), "+",
                              paste(as.character(random_model)[2], collapse = " "))

        tryCatch({
          m <- lme4::lmer(as.formula(full_formula), data = env_data)

          vc <- as.data.frame(lme4::VarCorr(m))
          vg <- vc[vc$grp == gen_col, "vcov"]
          ve <- vc[vc$grp == "Residual", "vcov"]
          if(length(vg)==0) vg <- NA # Handle singular fits

          # Simple H2 calculation
          n_reps <- mean(table(env_data[[gen_col]]))
          h2_simple <- if(!is.na(vg)) vg / (vg + (ve / n_reps)) else NA

          # BLUPs
          blups <- lme4::ranef(m)[[gen_col]]
          blups$Genotype <- rownames(blups)
          mu <- lme4::fixef(m)["(Intercept)"]
          blups$Predicted <- blups[[1]] + mu

          results_blup[[paste(env, trait)]] <- blups |>
            dplyr::select(Genotype, Predicted) |>
            dplyr::mutate(Environment = env, Trait = trait, Method = "lme4")

          results_var[[paste(env, trait)]] <- data.frame(
            Environment = env,
            Trait = trait,
            H2_BroadSense = h2_simple,
            Var_Gen = vg,
            Var_Res = ve
          )
          cat(paste("✔ Solved:", env, trait, "| H2:", round(h2_simple, 2), "\n"))

        }, error = function(e) {
          warning(paste("✖ lme4 Failed for", env, trait, ":", e$message))
        })
      }
    }
  }

  return(list(
    blups = dplyr::bind_rows(results_blup),
    vars = dplyr::bind_rows(results_var),
    spatial_trends = dplyr::bind_rows(results_spatial)
  ))
}
