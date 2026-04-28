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
    tidyr::pivot_wider(
      names_from = Trait,
      values_from = Predicted,
      values_fn = list(Predicted = mean)
    )

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

#' Plot Candidate Drivers of a Wheat Legacy Trait
#'
#' Computes correlations between one selected wheat legacy target and all
#' candidate lentil predictor traits, separately by environment.
#'
#' @param data Long-format data frame with Genotype, Environment, Trait, Predicted.
#' @param target_trait Wheat legacy trait to explain, for example `"Wheat_Legacy_PRO"`.
#' @param predictor_traits Optional vector of predictor traits to keep.
#' @param method Correlation method. `"pearson"` by default.
#' @param min_pairs Minimum complete genotype pairs required per correlation.
#' @param top_n Optional number of strongest absolute correlations to display
#'   per environment.
#' @param label_p Logical. Add p-value labels to the plot?
#'
#' @return A list with `correlations` and `plot`.
#' @export
plot_legacy_driver_correlations <- function(data,
                                            target_trait,
                                            predictor_traits = NULL,
                                            method = "pearson",
                                            min_pairs = 5,
                                            top_n = 25,
                                            label_p = TRUE) {
  cor_df <- calculate_legacy_driver_correlations(
    data = data,
    target_trait = target_trait,
    predictor_traits = predictor_traits,
    method = method,
    min_pairs = min_pairs
  )

  if (nrow(cor_df) == 0) stop("No driver correlations could be calculated.")

  plot_df <- cor_df
  if (!is.null(top_n) && is.finite(top_n)) {
    plot_df <- plot_df |>
      dplyr::group_by(.data[["Environment"]]) |>
      dplyr::slice_max(abs(.data[["Correlation"]]), n = top_n, with_ties = FALSE) |>
      dplyr::ungroup()
  }

  plot_df <- plot_df |>
    dplyr::arrange(.data[["Environment"]], .data[["Correlation"]]) |>
    dplyr::mutate(
      Plot_Trait = paste(.data[["Environment"]], .data[["Predictor_Trait"]], sep = "___"),
      Plot_Trait = factor(.data[["Plot_Trait"]], levels = unique(.data[["Plot_Trait"]])),
      P_Label = dplyr::case_when(
        is.na(.data[["P_Value"]]) ~ "p = NA",
        .data[["P_Value"]] < 0.001 ~ "p < 0.001",
        TRUE ~ paste0("p = ", signif(.data[["P_Value"]], 2))
      ),
      H_Just = ifelse(.data[["Correlation"]] >= 0, 0, 1),
      Label_Y = .data[["Correlation"]] + ifelse(.data[["Correlation"]] >= 0, 0.03, -0.03)
    )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[["Plot_Trait"]], y = .data[["Correlation"]])) +
    ggplot2::geom_col(ggplot2::aes(fill = .data[["Correlation"]] > 0), alpha = 0.85) +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~Environment, scales = "free_y") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::scale_x_discrete(labels = function(x) sub("^.*___", "", x)) +
    ggplot2::scale_fill_manual(values = c("TRUE" = "forestgreen", "FALSE" = "firebrick"), guide = "none") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste("Candidate Lentil Drivers of", target_trait),
      subtitle = paste("Correlation method:", method, "| top absolute correlations per environment"),
      x = "",
      y = "Correlation with wheat legacy phenotype"
    )

  if (label_p) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(
        label = .data[["P_Label"]],
        y = .data[["Label_Y"]],
        hjust = .data[["H_Just"]]
      ),
      size = 3,
      color = "gray25"
    )
  }

  list(correlations = cor_df, plot = p)
}

#' Calculate Candidate Driver Correlations
#'
#' @param data Long-format data frame with Genotype, Environment, Trait, Predicted.
#' @param target_trait Wheat legacy trait to explain.
#' @param predictor_traits Optional vector of predictor traits to keep.
#' @param method Correlation method.
#' @param min_pairs Minimum complete genotype pairs required.
#'
#' @return Data frame of environment-specific correlations, R2, p-values, and FDR.
#' @export
calculate_legacy_driver_correlations <- function(data,
                                                 target_trait,
                                                 predictor_traits = NULL,
                                                 method = "pearson",
                                                 min_pairs = 5) {
  .lr_check_columns(data, c("Genotype", "Environment", "Trait", "Predicted"))

  wide_data <- data |>
    dplyr::select(dplyr::all_of(c("Genotype", "Environment", "Trait", "Predicted"))) |>
    tidyr::pivot_wider(
      names_from = "Trait",
      values_from = "Predicted",
      values_fn = list(Predicted = mean)
    )

  if (!target_trait %in% names(wide_data)) {
    stop(paste("Target trait not found:", target_trait))
  }

  candidate_traits <- setdiff(names(wide_data), c("Genotype", "Environment", target_trait))
  if (!is.null(predictor_traits)) candidate_traits <- intersect(candidate_traits, predictor_traits)
  if (length(candidate_traits) == 0) stop("No predictor traits available for correlation.")

  out <- list()
  envs <- unique(wide_data$Environment)

  for (env in envs) {
    env_data <- wide_data[wide_data$Environment == env, , drop = FALSE]

    for (predictor in candidate_traits) {
      pair_df <- env_data[, c(target_trait, predictor), drop = FALSE]
      pair_df <- pair_df[stats::complete.cases(pair_df), , drop = FALSE]
      n_pairs <- nrow(pair_df)
      if (n_pairs < min_pairs) next
      if (length(unique(pair_df[[target_trait]])) < 2 || length(unique(pair_df[[predictor]])) < 2) next

      ct <- tryCatch(
        stats::cor.test(pair_df[[predictor]], pair_df[[target_trait]], method = method),
        error = function(e) NULL
      )
      if (is.null(ct)) next

      cor_value <- unname(ct$estimate)
      out[[paste(env, predictor, sep = "::")]] <- data.frame(
        Environment = env,
        Target_Trait = target_trait,
        Predictor_Trait = predictor,
        Correlation = cor_value,
        R2 = cor_value^2,
        P_Value = ct$p.value,
        N_Pairs = n_pairs,
        Method = method
      )
    }
  }

  final <- dplyr::bind_rows(out)
  if (nrow(final) == 0) return(final)

  final |>
    dplyr::group_by(.data[["Environment"]], .data[["Target_Trait"]]) |>
    dplyr::mutate(FDR = stats::p.adjust(.data[["P_Value"]], method = "BH")) |>
    dplyr::ungroup() |>
    dplyr::arrange(.data[["Environment"]], dplyr::desc(abs(.data[["Correlation"]])))
}
