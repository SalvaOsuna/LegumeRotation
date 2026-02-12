#' Calculate Rotational Legacy Values
#'
#' @param data The merged dataframe from merge_rotation_data().
#' @param trait String. The cereal trait to analyze (e.g., "Yield_kg_ha").
#' @param env_col String. The Cereal Environment column.
#' @param rep_col String (Optional). Replicate column.
#'
#' @return A list containing a stats table and a plot.
#' @export
#' @importFrom lme4 lmer fixef
#' @importFrom dplyr group_by summarize mutate arrange desc
#' @importFrom ggplot2 ggplot aes geom_col coord_flip theme_minimal labs geom_hline
get_legacy_values <- function(data, trait, env_col, rep_col = NULL) {

  results_list <- list()
  envs <- unique(data[[env_col]])

  for (env in envs) {
    env_data <- data[data[[env_col]] == env, ]

    # Skip if trait is missing
    if (all(is.na(env_data[[trait]]))) next

    # Define Model: Trait ~ Previous_Crop_Genotype (Fixed) + Rep (Random)
    # We use Fixed effects because we want to compare specific lentil varieties.

    formula_str <- paste(trait, "~ Previous_Crop_Genotype")
    if (!is.null(rep_col)) formula_str <- paste(formula_str, "+ (1|", rep_col, ")")

    tryCatch({
      m <- lme4::lmer(as.formula(formula_str), data = env_data)

      # Extract Fixed Effects (The Legacy Values)
      # Note: lmer uses contrast coding (treatment contrasts by default).
      # To get the actual means, we need to be careful.
      # Simplest robust way for a user: Calculate raw means adjusted for design if possible,
      # but here we will extract the coefficients relative to the intercept.

      # BETTER APPROACH for "Rotational Values":
      # We calculate the deviation of each Previous Genotype from the Grand Mean of that Environment.

      coefs <- lme4::fixef(m)

      # This part requires 'emmeans' for perfect accuracy, but to keep dependencies low
      # we can use aggregate() if the design is balanced, or lme4 predicted values.
      # Let's use the predicted values averaged by group (Robust method).

      env_data$Predicted <- predict(m)

      legacy_table <- env_data |>
        dplyr::group_by(Previous_Crop_Genotype) |>
        dplyr::summarize(
          Mean_Wheat_Yield = mean(Predicted, na.rm=TRUE),
          .groups = "drop"
        ) |>
        dplyr::mutate(
          Environment = env,
          # The "Legacy Value" is the difference from the site average
          Legacy_Value = Mean_Wheat_Yield - mean(Mean_Wheat_Yield),
          Legacy_Pct = (Legacy_Value / mean(Mean_Wheat_Yield)) * 100
        ) |>
        dplyr::arrange(dplyr::desc(Legacy_Value))

      results_list[[env]] <- legacy_table

    }, error = function(e) {
      warning(paste("Model failed for", env, ":", e$message))
    })
  }

  final_df <- dplyr::bind_rows(results_list)

  # Plotting
  p <- ggplot2::ggplot(final_df, ggplot2::aes(x = reorder(Previous_Crop_Genotype, Legacy_Value), y = Legacy_Value, fill = Legacy_Value > 0)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~Environment, scales = "free_x") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::scale_fill_manual(values = c("TRUE" = "forestgreen", "FALSE" = "firebrick"), guide = "none") +
    ggplot2::labs(
      title = "Rotational Legacy Effect",
      subtitle = paste("Effect of Previous Lentil on", trait),
      y = paste("Change in", trait, "(vs Site Mean)"),
      x = "Previous Lentil Genotype"
    ) +
    ggplot2::theme_minimal()

  return(list(values = final_df, plot = p))
}
