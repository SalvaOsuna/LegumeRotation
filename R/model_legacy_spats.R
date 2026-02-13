#' Calculate Robust Legacy Values using SpATS
#'
#' @param data The merged dataframe from merge_rotation_data().
#' @param trait String. The cereal trait to analyze (e.g., "Yield_kg_ha").
#' @param env_col String. The Cereal Environment column.
#' @param prev_gen_col String. Name of the column for Previous Crop Genotype.
#' @param curr_gen_col String. Name of the column for Current Crop (Wheat) Genotype.
#' @param spatial_cols Character vector. c("Row", "Col") for the Cereal field.
#' @param rep_col String (Optional).
#'
#' @return A list containing the Legacy BLUEs, plot, and models.
#' @export
#' @importFrom SpATS SpATS predict.SpATS PSANOVA
#' @importFrom dplyr mutate arrange desc bind_rows
#' @importFrom ggplot2 ggplot aes geom_col coord_flip theme_minimal labs geom_hline geom_errorbar
model_legacy_spats <- function(data, trait, env_col,
                               prev_gen_col = "Previous_Crop_Genotype",
                               curr_gen_col = "Genotype",
                               spatial_cols = c("Row", "Col"),
                               rep_col = NULL) {

  results_list <- list()
  models_list <- list()

  envs <- unique(data[[env_col]])

  for (env in envs) {
    # 1. Filter and Clean Data for this Env
    env_data <- data[data[[env_col]] == env, ]
    env_data <- env_data[!is.na(env_data[[trait]]), ]

    if(nrow(env_data) == 0) next

    # 2. Prepare Factors
    # Convert Previous Crop (Treatment) and Current Crop (Random) to factors
    env_data[[prev_gen_col]] <- as.factor(env_data[[prev_gen_col]])
    env_data[[curr_gen_col]] <- as.factor(env_data[[curr_gen_col]])

    # 3. CRITICAL FIX: Add Numeric Coordinates TO THE DATAFRAME
    # SpATS requires these to be columns in 'data', not external vectors.
    row_col_name <- spatial_cols[1]
    col_col_name <- spatial_cols[2]

    env_data$Row_Num <- as.numeric(as.character(env_data[[row_col_name]]))
    env_data$Col_Num <- as.numeric(as.character(env_data[[col_col_name]]))

    # 4. Dynamic Segments Calculation
    nseg_r <- max(1, floor(max(env_data$Row_Num, na.rm=TRUE) / 2))
    nseg_c <- max(1, floor(max(env_data$Col_Num, na.rm=TRUE) / 2))

    # 5. Define Formulas
    # Fixed Formula: Must use the exact column name string
    fixed_form_str <- paste("~", prev_gen_col)
    if (!is.null(rep_col)) {
      env_data[[rep_col]] <- as.factor(env_data[[rep_col]])
      fixed_form_str <- paste(fixed_form_str, "+", rep_col)
    }

    cat(paste("Modeling Legacy for:", env, "...\n"))

    tryCatch({
      # 6. Fit SpATS
      # Note: We now use Row_Num and Col_Num inside PSANOVA
      m <- SpATS::SpATS(
        response = trait,
        genotype = curr_gen_col,
        genotype.as.random = TRUE,
        fixed = as.formula(fixed_form_str),
        spatial = ~ SpATS::PSANOVA(Col_Num, Row_Num, nseg = c(nseg_c, nseg_r), degree = 3, pord = 2),
        data = env_data,
        control = list(tolerance = 1e-03, monitoring = 0)
      )

      models_list[[env]] <- m

      # 7. Extract BLUEs (Predicted values for Previous Crop)
      preds <- predict(m, which = prev_gen_col)

      # Clean up the prediction frame
      legacy_df <- preds[, c(prev_gen_col, "predicted.values", "standard.errors")]
      colnames(legacy_df) <- c("Previous_Genotype", "Predicted_Mean", "SE")

      # Calculate Legacy Value (Deviation from Grand Mean of the Env)
      grand_mean <- mean(legacy_df$Predicted_Mean, na.rm=TRUE)

      legacy_df <- legacy_df |>
        dplyr::mutate(
          Environment = env,
          Legacy_Value = Predicted_Mean - grand_mean,
          Legacy_Pct = (Legacy_Value / grand_mean) * 100
        ) |>
        dplyr::arrange(dplyr::desc(Legacy_Value))

      results_list[[env]] <- legacy_df

    }, error = function(e) {
      warning(paste("SpATS Legacy Model Failed for", env, ":", e$message))
    })
  }

  final_df <- dplyr::bind_rows(results_list)

  # 8. Plot
  p <- ggplot2::ggplot(final_df, ggplot2::aes(x = reorder(Previous_Genotype, Legacy_Value), y = Legacy_Value)) +
    ggplot2::geom_col(ggplot2::aes(fill = Legacy_Value > 0), alpha = 0.8) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = Legacy_Value - SE, ymax = Legacy_Value + SE), width = 0.2, color = "gray40") +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~Environment, scales = "free_x") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::scale_fill_manual(values = c("TRUE" = "forestgreen", "FALSE" = "firebrick"), guide = "none") +
    ggplot2::labs(
      title = "Rotational Legacy Values (SpATS Corrected)",
      subtitle = paste("Corrected for Wheat Genotype & Spatial Trend | Trait:", trait),
      y = "Legacy Effect (Deviation from Mean)",
      x = "Previous Lentil Genotype"
    ) +
    ggplot2::theme_minimal()

  return(list(legacy_values = final_df, plot = p, models = models_list))
}




##############################################################
#' Calculate Robust Legacy Values using SpATS (With Network Correction)
#'
#' @param data The merged dataframe from merge_rotation_data().
#' @param trait String. The cereal trait to analyze.
#' @param env_col String.
#' @param prev_gen_col String.
#' @param curr_gen_col String.
#' @param spatial_cols Character vector.
#' @param rep_col String (Optional).
#'
#' @return A list containing the Legacy BLUEs (with Raw vs Corrected), plot, and models.
#' @export
#' @importFrom SpATS SpATS predict.SpATS PSANOVA
#' @importFrom dplyr mutate arrange desc bind_rows group_by summarize left_join
#' @importFrom ggplot2 ggplot aes geom_point geom_segment coord_flip theme_minimal labs geom_hline
model_legacy_spats2 <- function(data, trait, env_col,
                               prev_gen_col = "Previous_Crop_Genotype",
                               curr_gen_col = "Genotype",
                               spatial_cols = c("Row", "Col"),
                               rep_col = NULL) {

  results_list <- list()
  models_list <- list()

  envs <- unique(data[[env_col]])

  for (env in envs) {
    # 1. Filter Data
    env_data <- data[data[[env_col]] == env, ]
    env_data <- env_data[!is.na(env_data[[trait]]), ]
    if(nrow(env_data) == 0) next

    # 2. Prepare Columns (Factors & Numerics)
    env_data[[prev_gen_col]] <- as.factor(env_data[[prev_gen_col]])
    env_data[[curr_gen_col]] <- as.factor(env_data[[curr_gen_col]])

    row_col_name <- spatial_cols[1]
    col_col_name <- spatial_cols[2]
    env_data$Row_Num <- as.numeric(as.character(env_data[[row_col_name]]))
    env_data$Col_Num <- as.numeric(as.character(env_data[[col_col_name]]))

    # 3. Fit SpATS Model
    nseg_r <- max(1, floor(max(env_data$Row_Num, na.rm=TRUE) / 2))
    nseg_c <- max(1, floor(max(env_data$Col_Num, na.rm=TRUE) / 2))

    fixed_form_str <- paste("~", prev_gen_col)
    if (!is.null(rep_col)) {
      env_data[[rep_col]] <- as.factor(env_data[[rep_col]])
      fixed_form_str <- paste(fixed_form_str, "+", rep_col)
    }

    cat(paste("Modeling Legacy for:", env, "...\n"))

    tryCatch({
      m <- SpATS::SpATS(
        response = trait,
        genotype = curr_gen_col,       # Wheat is Random (The Correction Mechanism)
        genotype.as.random = TRUE,
        fixed = as.formula(fixed_form_str),
        spatial = ~ SpATS::PSANOVA(Col_Num, Row_Num, nseg = c(nseg_c, nseg_r), degree = 3, pord = 2),
        data = env_data,
        control = list(tolerance = 1e-03, monitoring = 0)
      )
      models_list[[env]] <- m

      # 4. Extract "Corrected" Values (BLUEs)
      # This predicts how the Lentil would perform with an AVERAGE Wheat
      preds <- predict(m, which = prev_gen_col)
      legacy_df <- preds[, c(prev_gen_col, "predicted.values", "standard.errors")]
      colnames(legacy_df) <- c("Previous_Genotype", "Corrected_Mean", "SE")

      # 5. Calculate "Raw" Means for Comparison
      # This allows us to see the magnitude of the correction
      raw_stats <- env_data |>
        dplyr::group_by(.data[[prev_gen_col]]) |>
        dplyr::summarize(Raw_Mean = mean(.data[[trait]], na.rm=TRUE))

      # 6. Merge and Calculate Legacy Values
      grand_mean_corrected <- mean(legacy_df$Corrected_Mean, na.rm=TRUE)

      legacy_df <- legacy_df |>
        dplyr::left_join(raw_stats, by = setNames(prev_gen_col, "Previous_Genotype")) |>
        dplyr::mutate(
          Environment = env,
          # The Network Correction = How much the model penalized/boosted based on Wheat Partners
          Network_Correction = Corrected_Mean - Raw_Mean,

          # The Final Legacy Value (Deviation from Grand Mean)
          Legacy_Value = Corrected_Mean - grand_mean_corrected
        ) |>
        dplyr::arrange(dplyr::desc(Legacy_Value))

      results_list[[env]] <- legacy_df

    }, error = function(e) {
      warning(paste("Model Failed for", env, ":", e$message))
    })
  }

  final_df <- dplyr::bind_rows(results_list)

  # 7. Visualization: Dumbbell Plot (Raw vs Corrected)
  # This visually PROVES the correction to the user
  p <- ggplot2::ggplot(final_df, ggplot2::aes(y = reorder(Previous_Genotype, Legacy_Value))) +
    # The line connecting Raw to Corrected
    ggplot2::geom_segment(ggplot2::aes(x = Raw_Mean, xend = Corrected_Mean, yend = Previous_Genotype), color = "gray60") +
    # The Raw Point (Hollow)
    ggplot2::geom_point(ggplot2::aes(x = Raw_Mean), color = "black", shape = 1, size = 2) +
    # The Corrected Point (Solid Color based on Legacy)
    ggplot2::geom_point(ggplot2::aes(x = Corrected_Mean, color = Legacy_Value > 0), size = 3) +

    ggplot2::facet_wrap(~Environment, scales = "free_x") +
    ggplot2::scale_color_manual(values = c("TRUE" = "forestgreen", "FALSE" = "firebrick")) +
    ggplot2::labs(
      title = "Raw vs. Network-Corrected Legacy",
      subtitle = "Hollow Circle = Raw Mean (Biased) | Solid Circle = SpATS Corrected (Fair)\nThe distance between points represents the correction for Wheat Genotype + Spatial Field",
      x = paste("Yield (", trait, ")"),
      y = "Previous Lentil Genotype",
      color = "Positive Legacy?"
    ) +
    ggplot2::theme_minimal()

  return(list(legacy_values = final_df, plot = p, models = models_list))
}
