#' Plot Spatial Trends as Heatmaps
#'
#' @param spatial_data The 'spatial_trends' dataframe from model_traits().
#' @param traits Character vector. Specific traits to plot.
#' @param envs Character vector. Specific environments to plot.
#' @param mode String. "raw" (default) or "percentage".
#'             "percentage" normalizes the trend as % deviation from the mean.
#'
#' @return A ggplot object.
#' @export
#' @importFrom dplyr filter
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradientn facet_grid scale_y_reverse theme_minimal labs theme element_text
plot_spatial_trends <- function(spatial_data, traits = NULL, envs = NULL, mode = "raw") {

  # Filter data
  plot_data <- spatial_data
  if (!is.null(traits)) plot_data <- dplyr::filter(plot_data, Trait %in% traits)
  if (!is.null(envs)) plot_data <- dplyr::filter(plot_data, Environment %in% envs)

  if (nrow(plot_data) == 0) stop("No data found for the specified traits/environments.")

  # Select Column and Label based on Mode
  if (mode == "percentage") {
    fill_col <- "Trend_Pct"
    legend_lab <- "Deviation (%)"
    subtitle_lab <- "Spatial variation as % of trait mean"
  } else {
    fill_col <- "Trend"
    legend_lab <- "Deviation (Raw)"
    subtitle_lab <- "Spatial variation in original units"
  }

  spatial_palette <- grDevices::colorRampPalette(c("#c44e52", "#f8f9fa", "#2c7a51"))(256)
  fill_range <- range(plot_data[[fill_col]], na.rm = TRUE)
  fill_values <- NULL
  if (all(is.finite(fill_range)) && fill_range[1] < 0 && fill_range[2] > 0) {
    fill_values <- (c(fill_range[1], 0, fill_range[2]) - fill_range[1]) / diff(fill_range)
  }

  # Create Heatmap
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Col, y = Row, fill = .data[[fill_col]])) +
    ggplot2::geom_tile() +

    ggplot2::scale_fill_gradientn(
      colors = spatial_palette,
      values = fill_values,
      name = legend_lab
    ) +

    ggplot2::facet_grid(Trait ~ Environment, scales = "free") +

    ggplot2::scale_y_reverse() +

    ggplot2::theme_minimal() +
    ggplot2::labs(
      #title = paste("Spatial Trends (", tools::toTitleCase(mode), ")", sep=""),
      #subtitle = subtitle_lab,
      x = "Column",
      y = "Row"
    ) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold", size = 10),
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 8)
    )

  return(p)
}
