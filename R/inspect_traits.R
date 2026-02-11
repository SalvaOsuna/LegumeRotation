#' Inspect Phenotypic Traits (Stats & Plots)
#'
#' @param data A data.frame.Ideally the 'clean_data' output from inspect_met().
#' @param trait_cols Character vector. Names of the columns to analyze.
#' @param env_col String (Optional). The column defining the Environment. If NULL, ignores groups.
#' @param plot Logical. If TRUE, generates distribution plots.
#'
#' @return A list containing:
#'   - summary: A data.frame of summary statistics.
#'   - plot: A ggplot object (if plot=TRUE).
#' @export
#' @importFrom dplyr group_by summarize across all_of n mean sd median min max
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_histogram geom_violin geom_boxplot facet_wrap theme_minimal labs
inspect_traits <- function(data, trait_cols, env_col = NULL, plot = TRUE) {

  # 1. Check if traits exist
  missing <- setdiff(trait_cols, names(data))
  if (length(missing) > 0) stop(paste("Traits not found:", paste(missing, collapse=", ")))

  # 2. Pivot to Long Format (Tidy) for easier processing
  # This stacks all traits into one column "Trait" and values into "Value"
  data_long <- data |>
    tidyr::pivot_longer(
      cols = all_of(trait_cols),
      names_to = "Trait",
      values_to = "Value"
    ) |>
    # Remove NAs immediately to avoid errors in stats
    dplyr::filter(!is.na(Value))

  # 3. Calculate Summary Statistics
  # Define the grouping variables (always by Trait, optionally by Env)
  group_vars <- "Trait"
  if (!is.null(env_col)) group_vars <- c(env_col, "Trait")

  summary_stats <- data_long |>
    dplyr::group_by(across(all_of(group_vars))) |>
    dplyr::summarize(
      N = n(),
      Mean = mean(Value),
      Median = median(Value),
      Min = min(Value),
      Max = max(Value),
      SD = sd(Value),
      CV_percent = (sd(Value) / mean(Value)) * 100,
      .groups = "drop"
    )

  # 4. Generate Plots (Conditional)
  p <- NULL
  if (plot) {
    # Logic: If we have multiple Environments, use Violin Plots
    # If we have only 1 Env (or NULL), use Histograms

    is_multi_env <- !is.null(env_col) && length(unique(data[[env_col]])) > 1

    if (is_multi_env) {
      # VIOLIN PLOT: Value vs Environment (Faceted by Trait)
      p <- ggplot2::ggplot(data_long, ggplot2::aes(x = .data[[env_col]], y = Value, fill = .data[[env_col]])) +
        ggplot2::geom_violin(alpha = 0.6, trim = FALSE) +
        ggplot2::geom_boxplot(width = 0.1, fill = "white", outliers = FALSE) + # Boxplot inside for precision
        ggplot2::facet_wrap(~Trait, scales = "free_y") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
        ggplot2::labs(title = "Trait Distribution by Environment", y = "Value")

    } else {
      # HISTOGRAM: Distribution of values (Faceted by Trait)
      p <- ggplot2::ggplot(data_long, ggplot2::aes(x = Value)) +
        ggplot2::geom_histogram(bins = 30, fill = "steelblue", color = "white", alpha = 0.8) +
        ggplot2::facet_wrap(~Trait, scales = "free") +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = "Trait Distribution (Global)", y = "Count")
    }
  }

  # 5. Return Results
  return(list(
    stats = summary_stats,
    plot = p
  ))
}
