#' Merge Legume and Cereal Datasets
#'
#' @param cereal_data The clean dataframe for Year 2 (Wheat).
#' @param legume_data The clean dataframe for Year 1 (Lentil).
#' @param link_col String. The common column (Unique.ID).
#' @param legume_gen_col String. The column name of the lentil variety in legume_data.
#'
#' @return A merged dataframe ready for rotation analysis.
#' @export
#' @importFrom dplyr select rename left_join all_of
merge_rotation_data <- function(cereal_data, legume_data, link_col = "Unique.ID", legume_gen_col = "Genotype") {

  # 1. Prepare the Legume Data (The Predictor)
  # We only need the ID and the Genotype (and maybe the Environment if needed for matching)
  legume_subset <- legume_data |>
    dplyr::select(dplyr::all_of(c(link_col, legume_gen_col))) |>
    dplyr::rename(Previous_Crop_Genotype = !!legume_gen_col)

  # 2. Join into Cereal Data (The Response)
  merged_df <- cereal_data |>
    dplyr::left_join(legume_subset, by = link_col)

  # 3. Safety Check
  missing_links <- sum(is.na(merged_df$Previous_Crop_Genotype))
  if (missing_links > 0) {
    warning(paste("Warning:", missing_links, "wheat plots could not be linked to a previous lentil genotype. Check your Unique.IDs!"))
  } else {
    message("✔ Successfully linked all rotation plots.")
  }

  return(merged_df)
}
