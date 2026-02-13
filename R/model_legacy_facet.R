#' Calculate Legacy Values using Facet/Partner Deviation
#'
#' @param data The merged dataframe from merge_rotation_data().
#' @param trait String. The cereal trait to analyze (e.g., "Yield_kg_ha").
#' @param env_col String.
#' @param prev_gen_col String. The Lentil Genotype (Previous Crop).
#' @param curr_gen_col String. The Wheat Genotype (Current Crop).
#'
#' @return A list containing the Facet-Corrected Legacy table and a plot.
#' @export
#' @importFrom dplyr group_by summarize mutate left_join filter select
#' @importFrom ggplot2 ggplot aes geom_col coord_flip theme_minimal labs geom_hline
model_legacy_facet <- function(data, trait, env_col,
                               prev_gen_col = "Previous_Crop_Genotype",
                               curr_gen_col = "Genotype",
                               spatial_cols = c("Row", "Col"),
                               rep_col = NULL) {

  results_list <- list()

  envs <- unique(data[[env_col]])

  for (env in envs) {
    # 1. Isolate Data for this Environment
    env_data <- data |> dplyr::filter(ENV %in% env)
    env_data <- na.omit(env_data)

    if(nrow(env_data) == 0) next

    # 2. Calculate Global Mean per Wheat Genotype (The Baseline)
    # This tells us how Wheat X performs on average across the whole trial (all Reps, all Lentils)
    #wheat_global_means <- env_data |>
    #  dplyr::group_by(.data[[curr_gen_col]]) |>
    #  dplyr::summarize(Global_Wheat_Mean = mean(.data[[trait]], na.rm=TRUE))

    # 2. Calculate Global Mean per Wheat Genotype using SpATS (The Baseline)
    cat(paste("Calculating Wheat Baseline BLUPs for:", env, "...\n"))

    # Ensure Current Crop is a factor
    env_data[[curr_gen_col]] <- as.factor(env_data[[curr_gen_col]])

    # Ensure Spatial columns are numeric
    row_col_name <- spatial_cols[1]
    col_col_name <- spatial_cols[2]
    env_data$Row_Num <- as.numeric(as.character(env_data[[row_col_name]]))
    env_data$Col_Num <- as.numeric(as.character(env_data[[col_col_name]]))

    # Calculate dynamic segments
    nseg_r <- max(1, floor(max(env_data$Row_Num, na.rm=TRUE) / 2))
    nseg_c <- max(1, floor(max(env_data$Col_Num, na.rm=TRUE) / 2))

    # Define Fixed Formula
    fixed_form <- "~ 1"
    if (!is.null(rep_col)) fixed_form <- paste("~", rep_col)

    # Fit SpATS to get the pure Wheat BLUPs (ignoring Lentil effect for the baseline)
    m_wheat <- SpATS::SpATS(
      response = trait,
      genotype = curr_gen_col,
      genotype.as.random = TRUE, # Treat Wheat as random to get BLUPs
      spatial = ~ SpATS::PSANOVA(Col_Num, Row_Num, nseg = c(nseg_c, nseg_r), degree = 3, pord = 2),
      fixed = as.formula(fixed_form),
      data = env_data,
      control = list(tolerance = 1e-03, monitoring = 0)
    )

    # Extract the predicted values (BLUP + Intercept)
    wheat_preds <- predict(m_wheat, which = curr_gen_col)

    # Format to match the expected downstream dataframe
    wheat_global_means <- wheat_preds[, c(curr_gen_col, "predicted.values")]
    colnames(wheat_global_means) <- c(curr_gen_col, "Global_Wheat_Mean")
    # ====================================================================

    # 3. Process Each Lentil (The Facet Analysis)
    lentil_results <- env_data |>
      dplyr::group_by(.data[[prev_gen_col]]) |>
      dplyr::summarize(
        Observed_Yield = mean(.data[[trait]], na.rm=TRUE),

        # KEY STEP: Calculate the "Facet Mean"
        # We look at exactly which Wheats were paired with this Lentil
        # And we take the average of THEIR Global Means.
        Facet_Mean = mean(
          wheat_global_means$Global_Wheat_Mean[
            wheat_global_means[[curr_gen_col]] %in% unique(.data[[curr_gen_col]])
          ],
          na.rm = TRUE
        ),

        # Store the list of partners for reference (optional but useful)
        N_Wheat_Partners = dplyr::n_distinct(.data[[curr_gen_col]])
      ) |>
      dplyr::mutate(
        Environment = env,
        # The Legacy is the specific bump this Lentil gave its partners
        # relative to how those partners performed elsewhere
        Legacy_Value = Observed_Yield - Facet_Mean,
        Legacy_Pct = (Legacy_Value / Facet_Mean) * 100
      ) |>
      dplyr::arrange(dplyr::desc(Legacy_Value))

    results_list[[env]] <- lentil_results
  }

  final_df <- dplyr::bind_rows(results_list)

  # 4. Plotting (Deviation Plot)
  p <- ggplot2::ggplot(final_df, ggplot2::aes(x = reorder(Previous_Crop_Genotype, Legacy_Value), y = Legacy_Value)) +
    ggplot2::geom_col(ggplot2::aes(fill = Legacy_Value > 0)) +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~Environment, scales = "free_x") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::scale_fill_manual(values = c("TRUE" = "forestgreen", "FALSE" = "firebrick"), guide = "none") +
    ggplot2::labs(
      title = "Facet-Corrected Legacy Values",
      subtitle = "Deviation of Observed Yield from the expected mean of the specific Wheat Partners",
      x = "Previous Lentil Genotype",
      y = "Legacy Effect (Observed - Facet Baseline)"
    ) +
    ggplot2::theme_minimal()

  return(list(legacy_values = final_df, plot = p))
}
