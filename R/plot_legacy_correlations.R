#' Plot Trait Correlation Matrix by Environment
#'
#' @param data The combined long-format dataframe (Lentil BLUPs + Wheat Legacy).
#' @param target_trait String (Optional). If provided, orders the plot to put this trait first.
#'
#' @return A ggplot object displaying a faceted correlation heatmap.
#' @export
#' @importFrom dplyr select filter mutate rename bind_rows
#' @importFrom tidyr pivot_wider
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_gradient2 theme_minimal labs theme element_text coord_fixed facet_wrap
plot_legacy_correlations <- function(data, target_trait = "Wheat_Legacy_Benefit") {

  # 1. Pivot to Wide Format (Rows = Genotypes, Cols = Traits)
  wide_data <- data |>
    dplyr::select(Genotype, Environment, Trait, Predicted) |>
    tidyr::pivot_wider(names_from = Trait, values_from = Predicted)

  envs <- unique(wide_data$Environment)
  cor_list <- list()

  # 2. Calculate Correlation Matrix per Environment
  for (env in envs) {
    # Isolate data for this environment and remove metadata columns
    env_data <- wide_data |>
      dplyr::filter(Environment == env) |>
      dplyr::select(-Genotype, -Environment)

    # Drop columns that are entirely NA in this environment
    env_data <- env_data[, colSums(is.na(env_data)) < nrow(env_data)]

    # Skip if less than 2 traits available
    if (ncol(env_data) < 2) next

    # Calculate Pearson correlation (handling missing data pairwise)
    cor_mat <- cor(env_data, use = "pairwise.complete.obs", method = "pearson")

    # Flatten the matrix for ggplot
    cor_df <- as.data.frame(as.table(cor_mat)) |>
      dplyr::rename(Trait_1 = Var1, Trait_2 = Var2, Correlation = Freq) |>
      dplyr::mutate(Environment = env)

    cor_list[[env]] <- cor_df
  }

  if (length(cor_list) == 0) stop("Not enough data to calculate correlations.")

  final_cor <- dplyr::bind_rows(cor_list)

  # 3. Order the factors so the Target Trait is at the top/left for easy reading
  if (target_trait %in% final_cor$Trait_1) {
    trait_levels <- unique(c(target_trait, as.character(unique(final_cor$Trait_1))))
    final_cor$Trait_1 <- factor(final_cor$Trait_1, levels = trait_levels)
    final_cor$Trait_2 <- factor(final_cor$Trait_2, levels = rev(trait_levels)) # Reverse for y-axis
  }

  # 4. Generate the Heatmap
  p <- ggplot2::ggplot(final_cor, ggplot2::aes(x = Trait_1, y = Trait_2, fill = Correlation)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +

    # Add correlation numbers inside the tiles
    ggplot2::geom_text(ggplot2::aes(label = round(Correlation, 2)), size = 3, color = "black") +

    # Red-White-Green color scale (Green = Positive correlation, Red = Negative)
    ggplot2::scale_fill_gradient2(
      low = "firebrick", mid = "white", high = "forestgreen",
      midpoint = 0, limit = c(-1, 1), space = "Lab",
      name = "Pearson (r)"
    ) +

    ggplot2::facet_wrap(~Environment) +
    ggplot2::coord_fixed() + # Keeps the tiles perfectly square

    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Trait Correlations across Environments",
      subtitle = paste("Identifying drivers of", target_trait),
      x = "",
      y = ""
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1, face = "bold"),
      axis.text.y = ggplot2::element_text(face = "bold"),
      panel.grid.major = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = 11)
    )

  return(p)
}
