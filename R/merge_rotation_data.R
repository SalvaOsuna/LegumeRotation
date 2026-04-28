#' Merge Legume and Cereal Datasets
#'
#' @param cereal_data The clean dataframe for Year 2 (Wheat).
#' @param legume_data The clean dataframe for Year 1 (Lentil). Optional when
#'   `cereal_data` already contains the previous lentil genotype.
#' @param link_col String. The common plot-link column.
#' @param legume_gen_col String. The column name of the lentil variety in `legume_data`.
#' @param cereal_prev_col Optional previous lentil genotype column already present
#'   in `cereal_data`, for example `"Lentil"`.
#'
#' @return A merged dataframe ready for rotation analysis.
#' @export
#' @importFrom dplyr select rename left_join all_of
merge_rotation_data <- function(cereal_data, legume_data = NULL, link_col = "Unique.ID",
                                legume_gen_col = "Genotype", cereal_prev_col = NULL) {

  # Many ACTIVATE wheat files already carry the predecessor lentil genotype.
  # In that case, avoid an unnecessary join through plot IDs that may have
  # crop-specific formatting.
  if (!is.null(cereal_prev_col) && cereal_prev_col %in% names(cereal_data)) {
    cereal_data$Previous_Crop_Genotype <- cereal_data[[cereal_prev_col]]
    message("Using '", cereal_prev_col, "' from cereal_data as Previous_Crop_Genotype.")
    return(cereal_data)
  }

  if (is.null(cereal_prev_col) && "Lentil" %in% names(cereal_data)) {
    cereal_data$Previous_Crop_Genotype <- cereal_data$Lentil
    message("Using 'Lentil' from cereal_data as Previous_Crop_Genotype.")
    return(cereal_data)
  }

  if (is.null(legume_data)) {
    stop("Provide legume_data, or set cereal_prev_col to a previous lentil genotype column in cereal_data.")
  }

  if (!link_col %in% names(cereal_data) || !link_col %in% names(legume_data)) {
    stop(paste("The link column", link_col, "must exist in both datasets."))
  }

  if (!legume_gen_col %in% names(legume_data) && "Lentil" %in% names(legume_data)) {
    legume_gen_col <- "Lentil"
  }

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
