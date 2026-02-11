#' Inspect MET Design and Verify Plot Links
#'
#' @param data A data.frame containing the raw phenotypic data.
#' @param unique_id_col String. The column linking legume and cereal plots (must be unique).
#' @param env_col String. The column identifying the Environment (Site/Year).
#' @param genotype_col String. The name of the column containing genotype/variety IDs.
#' @param design_cols Character vector. Columns defining the design hierarchy (e.g., c("Rep", "Block")).
#' @param trait_cols Character vector. Names of columns containing phenotypic traits.
#' @param crop_col String. The name of the column containing the crop name. Defaults to "Crop".
#'
#' @return A list containing a summary report and the cleaned data.
#' @export
#' @importFrom dplyr group_by summarize n n_distinct filter select all_of
#'
inspect_met <- function(data, unique_id_col, env_col, genotype_col, design_cols, trait_cols, crop_col = "Crop") {

  # 1. Basic Column Check
  required_cols <- c(unique_id_col, env_col, genotype_col, design_cols, trait_cols, crop_col)
  missing_cols <- setdiff(required_cols, names(data))

  if (length(missing_cols) > 0) {
    stop(paste("Error: Missing columns:", paste(missing_cols, collapse = ", ")))
  }

  # 2. Check Unique ID Integrity (The "Golden Thread" Check)
  # A Unique.ID must appear ONLY ONCE in a single dataset (e.g., the Lentil dataset)
  id_counts <- table(data[[unique_id_col]])
  duplicates <- id_counts[id_counts > 1]

  if (length(duplicates) > 0) {
    stop(paste("CRITICAL ERROR: '", unique_id_col, "' is not unique! detected ",
               length(duplicates), " duplicate IDs. Check ID:", names(duplicates)[1]))
  }

  # 3. Check Balance BY Environment
  # We loop through environments because one site might be balanced and another unbalanced
  env_summary <- data |>
    dplyr::group_by(across(all_of(env_col))) |>
    dplyr::summarize(
      n_genotypes = dplyr::n_distinct(.data[[genotype_col]]),
      n_plots = dplyr::n(),
      # Check if all genotypes have same number of replicates in this environment
      is_balanced = (n_plots %% n_genotypes == 0)
    )

  # 4. Design Hierarchy Check (Nested)
  # Example: How many Blocks are in each Rep?
  design_structure <- lapply(design_cols, function(col) {
    counts <- dplyr::n_distinct(data[[col]])
    return(paste0(col, ": ", counts, " levels"))
  })

  # 5. Create Clean Dataset
  clean_df <- data |>
    dplyr::select(all_of(c(unique_id_col, env_col, crop_col, genotype_col, design_cols, trait_cols)))

  # 6. Print Report
  cat("=== MET Inspection Report ===\n")
  cat("Crop Identified:", unique(data[[crop_col]]), "\n")
  cat("Total Environments:", nrow(env_summary), "\n")
  cat("Unique ID Integrity: PASSED (All IDs are unique)\n")
  cat("-----------------------------\n")
  print(env_summary)
  cat("-----------------------------\n")
  cat("Design Structure (Global):\n", paste("  -", unlist(design_structure), collapse = "\n"), "\n")
  cat("=============================\n")

  return(list(
    env_summary = env_summary,
    clean_data = clean_df
  ))
}
