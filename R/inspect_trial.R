#' Inspect Trial Design and Clean Data
#'
#' @param data A data.frame containing the raw phenotypic data.
#' @param genotype_col String. The name of the column containing genotype/variety IDs.
#' @param design_cols Character vector. Names of columns defining the experimental design (e.g., c("Rep", "Block", "Row", "Col")).
#' @param trait_cols Character vector. Names of columns containing phenotypic traits (e.g., c("Yield", "DTF")).
#' @param crop_col String. The name of the column containing the crop name. Defaults to "Crop".
#'
#' @return A list containing a summary list and the cleaned data.frame.
#' @export
#' @importFrom dplyr select group_by summarize n_distinct filter all_of
inspect_trial <- function(data, genotype_col, design_cols, trait_cols, crop_col = "Crop") {

  # 1. Basic Checks
  required_cols <- c(genotype_col, design_cols, trait_cols, crop_col)
  missing_cols <- setdiff(required_cols, names(data))

  if (length(missing_cols) > 0) {
    stop(paste("Error: The following columns are missing from the dataset:",
               paste(missing_cols, collapse = ", ")))
  }

  # 2. Extract Crop Info
  crop_name <- unique(data[[crop_col]])
  if (length(crop_name) > 1) {
    warning("Multiple crops detected in the 'Crop' column. This might be a mixed dataset.")
  }

  # 3. Check Replicates (Balance)
  # Calculates how many times each genotype appears
  rep_counts <- data |>
    dplyr::group_by(across(all_of(genotype_col))) |>
    dplyr::summarize(Count = n())

  # Check if data is balanced (all genotypes have same number of reps)
  is_balanced <- length(unique(rep_counts$Count)) == 1
  balance_msg <- if (is_balanced) {
    paste("Balanced: All genotypes have", unique(rep_counts$Count), "replicates.")
  } else {
    paste("UNBALANCED: Replicates range from", min(rep_counts$Count), "to", max(rep_counts$Count))
  }

  # 4. Design Summary
  design_summary <- lapply(design_cols, function(col) {
    n_levels <- dplyr::n_distinct(data[[col]])
    paste0(col, ": ", n_levels, " levels")
  })

  # 5. Create Clean Dataset
  # Selects only the columns you asked for
  clean_df <- data |>
    dplyr::select(all_of(c(crop_col, genotype_col, design_cols, trait_cols)))

  # 6. Print Report to Console
  cat("=== Trial Inspection Report ===\n")
  cat("Crop:", paste(crop_name, collapse = ", "), "\n")
  cat("Number of Genotypes:", nrow(rep_counts), "\n")
  cat("Design Structure:\n", paste("  -", unlist(design_summary), collapse = "\n"), "\n")
  cat("Replicate Balance:", balance_msg, "\n")
  cat("Traits Identified:", paste(trait_cols, collapse = ", "), "\n")
  cat("===============================\n")

  # Return the list
  return(list(
    summary = list(
      crop = crop_name,
      n_genotypes = nrow(rep_counts),
      is_balanced = is_balanced,
      rep_stats = rep_counts
    ),
    clean_data = clean_df
  ))
}
