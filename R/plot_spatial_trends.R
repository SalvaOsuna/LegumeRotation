#' Plot Spatial Trends as Heatmaps
#'
#' @param spatial_data The 'spatial_trends' dataframe from model_traits() output.
#' @param traits Character vector. Specific traits to plot (optional).
#' @param envs Character vector. Specific environments to plot (optional).
#'
#' @return A ggplot object.
#' @export
#' @importFrom dplyr filter
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_viridis_c facet_grid theme_minimal labs theme element_text
plot_spatial_trends <- function(spatial_data, traits = NULL, envs = NULL) {

  # Filter data if requested
  plot_data <- spatial_data
  if (!is.null(traits)) plot_data <- dplyr::filter(plot_data, Trait %in% traits)
  if (!is.null(envs)) plot_data <- dplyr::filter(plot_data, Environment %in% envs)

  if (nrow(plot_data) == 0) stop("No data found for the specified traits/environments.")

  # Create Heatmap
  # We use facet_grid to create a matrix of Traits x Environments
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Col, y = Row, fill = Trend)) +
    ggplot2::geom_tile() +

    # Use Viridis palette (colorblind friendly and standard for scientific heatmaps)
    ggplot2::scale_fill_viridis_c(option = "magma", name = "Trend") +

    # Create the grid: Rows = Traits, Cols = Environments
    ggplot2::facet_grid(Trait ~ Environment, scales = "free") +

    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Spatial Trends (SpATS)",
      subtitle = "Smoothed spatial variation (independent of genotype)",
      x = "Column",
      y = "Row"
    ) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold", size = 10),
      panel.grid = ggplot2::element_blank() # Remove grid lines for cleaner heatmap
    )

  return(p)
}
