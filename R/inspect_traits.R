#' Inspect Phenotypic Traits (Stats & Plots)
#'
#' @param data A data.frame. Ideally the 'clean_data' output from inspect_met().
#' @param trait_cols Character vector. Names of the columns to analyze.
#' @param env_col String (Optional). The column defining the Environment. If NULL, ignores groups.
#' @param plot Logical. If TRUE, generates distribution plots.
#'
#' @return A list containing summary stats and a plot.
#' @export
#' @importFrom dplyr group_by summarize across all_of n mean sd median min max
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_histogram geom_violin geom_boxplot facet_wrap theme_minimal labs
inspect_traits <- function(data, trait_cols, env_col = NULL, plot = TRUE) {

  # 1. Check if traits exist
  missing <- setdiff(trait_cols, names(data))
  if (length(missing) > 0) stop(paste("Traits not found:", paste(missing, collapse=", ")))

  # 2. Pivot to Long Format
  data_long <- data |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(trait_cols),
      names_to = "Trait",
      values_to = "Value"
    ) |>
    dplyr::filter(!is.na(Value))

  # 3. Calculate Summary Statistics
  group_vars <- "Trait"
  if (!is.null(env_col)) group_vars <- c(env_col, "Trait")

  summary_stats <- data_long |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::summarize(
      N = dplyr::n(),              # <--- THE FIX: dplyr::n() instead of n()
      Mean = mean(Value, na.rm = TRUE),
      Median = median(Value, na.rm = TRUE),
      Min = min(Value, na.rm = TRUE),
      Max = max(Value, na.rm = TRUE),
      SD = sd(Value, na.rm = TRUE),
      CV_percent = (sd(Value, na.rm = TRUE) / mean(Value, na.rm = TRUE)) * 100,
      .groups = "drop"
    )

  # 4. Generate Plots
  p <- NULL
  if (plot) {
    is_multi_env <- !is.null(env_col) && length(unique(data[[env_col]])) > 1

    if (is_multi_env) {
      # Violin Plot for Multi-Environment
      p <- ggplot2::ggplot(data_long, ggplot2::aes(x = .data[[env_col]], y = Value, fill = .data[[env_col]])) +
        ggplot2::geom_violin(alpha = 0.6, trim = FALSE) +
        ggplot2::geom_boxplot(width = 0.1, fill = "white", outliers = FALSE) +
        ggplot2::facet_wrap(~Trait, scales = "free_y") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
        ggplot2::labs(title = paste0(data$Crop[1],"- Trait Distribution by Environment"), y = "Value")

    } else {
      # Histogram for Single Environment
      p <- ggplot2::ggplot(data_long, ggplot2::aes(x = Value)) +
        ggplot2::geom_histogram(bins = 30, fill = "steelblue", color = "white", alpha = 0.8) +
        ggplot2::facet_wrap(~Trait, scales = "free") +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = paste0(data$Crop[1],"- Trait Distribution (Global)"), y = "Count")
    }
  }

  return(list(stats = summary_stats, plot = p))
}
