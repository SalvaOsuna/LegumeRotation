#' Expected Lentil Facet Number
#'
#' @param lentil_id Lentil genotype ID, for example `"L011"`.
#'
#' @return Integer facet number using the ACTIVATE modulo assignment.
#' @export
expected_lentil_facet <- function(lentil_id) {
  n <- .lr_extract_first_integer(lentil_id)
  ((n - 1) %% 10) + 1
}

#' Expected Wheat Facet Number
#'
#' @param wheat_id Wheat genotype ID, for example `"W011"`.
#'
#' @return Integer facet number using blocks of ten wheat genotypes.
#' @export
expected_wheat_facet <- function(wheat_id) {
  n <- .lr_extract_first_integer(wheat_id)
  ceiling(n / 10)
}

#' Check Rotation Design Completeness
#'
#' Summarizes the observed local factorial structure within each environment
#' and facet.
#'
#' @param data Rotation data.
#' @param env_col Environment column.
#' @param facet_col Facet/network column.
#' @param lentil_col Previous lentil genotype column.
#' @param wheat_col Current wheat genotype column.
#' @param combo_col Lentil-wheat combination column.
#'
#' @return A data frame with local factorial completeness diagnostics.
#' @export
check_rotation_design <- function(data,
                                  env_col = "ENV",
                                  facet_col = "Facet",
                                  lentil_col = "Lentil",
                                  wheat_col = "Wheat",
                                  combo_col = "Combo") {
  required_cols <- c(env_col, facet_col, lentil_col, wheat_col, combo_col)
  .lr_check_columns(data, required_cols)

  data |>
    dplyr::distinct(
      .data[[env_col]],
      .data[[facet_col]],
      .data[[lentil_col]],
      .data[[wheat_col]],
      .data[[combo_col]]
    ) |>
    dplyr::group_by(.data[[env_col]], .data[[facet_col]]) |>
    dplyr::summarise(
      N_Previous_Genotypes = dplyr::n_distinct(.data[[lentil_col]]),
      N_Current_Genotypes = dplyr::n_distinct(.data[[wheat_col]]),
      N_Combos = dplyr::n_distinct(.data[[combo_col]]),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      Expected_Combos = .data[["N_Previous_Genotypes"]] * .data[["N_Current_Genotypes"]],
      Missing_Combos = .data[["Expected_Combos"]] - .data[["N_Combos"]],
      Is_Complete_Local_Factorial = .data[["N_Combos"]] == .data[["Expected_Combos"]],
      Expected_10x10 = .data[["N_Previous_Genotypes"]] == 10 &
        .data[["N_Current_Genotypes"]] == 10 &
        .data[["N_Combos"]] == 100
    )
}

#' Check Lentil-to-Facet Assignment
#'
#' @param data Rotation data.
#' @param lentil_col Previous lentil genotype column.
#' @param facet_col Facet/network column.
#'
#' @return A data frame with expected and observed facet assignment.
#' @export
check_lentil_facet_assignment <- function(data,
                                          lentil_col = "Lentil",
                                          facet_col = "Facet") {
  .lr_check_columns(data, c(lentil_col, facet_col))

  data |>
    dplyr::distinct(
      Lentil = .data[[lentil_col]],
      Facet = .data[[facet_col]]
    ) |>
    dplyr::mutate(
      Expected_Facet_Number = expected_lentil_facet(.data[["Lentil"]]),
      Observed_Facet_Number = .lr_facet_number(.data[["Facet"]]),
      Facet_Assignment_OK = .data[["Expected_Facet_Number"]] == .data[["Observed_Facet_Number"]]
    ) |>
    dplyr::arrange(.data[["Expected_Facet_Number"]], .data[["Lentil"]])
}

#' Check Wheat-to-Facet Assignment
#'
#' @param data Rotation data.
#' @param wheat_col Current wheat genotype column.
#' @param facet_col Facet/network column.
#'
#' @return A data frame with expected and observed facet assignment.
#' @export
check_wheat_facet_assignment <- function(data,
                                         wheat_col = "Wheat",
                                         facet_col = "Facet") {
  .lr_check_columns(data, c(wheat_col, facet_col))

  data |>
    dplyr::distinct(
      Wheat = .data[[wheat_col]],
      Facet = .data[[facet_col]]
    ) |>
    dplyr::mutate(
      Expected_Facet_Number = expected_wheat_facet(.data[["Wheat"]]),
      Observed_Facet_Number = .lr_facet_number(.data[["Facet"]]),
      Facet_Assignment_OK = .data[["Expected_Facet_Number"]] == .data[["Observed_Facet_Number"]]
    ) |>
    dplyr::arrange(.data[["Expected_Facet_Number"]], .data[["Wheat"]])
}

#' Audit ACTIVATE Rotation Design
#'
#' Runs local factorial completeness and genotype-to-facet assignment checks.
#'
#' @param data Rotation data.
#' @param env_col Environment column.
#' @param facet_col Facet/network column.
#' @param lentil_col Previous lentil genotype column.
#' @param wheat_col Current wheat genotype column.
#' @param combo_col Lentil-wheat combination column.
#'
#' @return A list with design summary, assignment checks, and `design_ok`.
#' @export
audit_rotation_design <- function(data,
                                  env_col = "ENV",
                                  facet_col = "Facet",
                                  lentil_col = "Lentil",
                                  wheat_col = "Wheat",
                                  combo_col = "Combo") {
  design_summary <- check_rotation_design(
    data = data,
    env_col = env_col,
    facet_col = facet_col,
    lentil_col = lentil_col,
    wheat_col = wheat_col,
    combo_col = combo_col
  )

  lentil_assignment <- check_lentil_facet_assignment(
    data = data,
    lentil_col = lentil_col,
    facet_col = facet_col
  )

  wheat_assignment <- check_wheat_facet_assignment(
    data = data,
    wheat_col = wheat_col,
    facet_col = facet_col
  )

  list(
    design_summary = design_summary,
    lentil_assignment = lentil_assignment,
    wheat_assignment = wheat_assignment,
    design_ok = all(design_summary$Expected_10x10, na.rm = TRUE) &&
      all(lentil_assignment$Facet_Assignment_OK, na.rm = TRUE) &&
      all(wheat_assignment$Facet_Assignment_OK, na.rm = TRUE)
  )
}

.lr_extract_first_integer <- function(x) {
  x_chr <- as.character(x)
  matches <- regexpr("[0-9]+", x_chr)
  out <- rep(NA_integer_, length(x_chr))
  hit <- matches > 0
  out[hit] <- suppressWarnings(as.integer(substr(
    x_chr[hit],
    matches[hit],
    matches[hit] + attr(matches, "match.length")[hit] - 1
  )))
  out
}

.lr_facet_number <- function(facet) {
  facet_chr <- as.character(facet)
  first <- .lr_extract_first_integer(facet_chr)
  has_range <- grepl("[0-9]+[^0-9]+[0-9]+", facet_chr)
  ifelse(has_range, ceiling(first / 10), first)
}
