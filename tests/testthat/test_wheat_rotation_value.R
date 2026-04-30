make_wheat_rotation_test_data <- function(n_lentil = 3, n_wheat = 3, reps = 2) {
  dat <- expand.grid(
    ENV = "E1",
    Facet = "F1",
    Lentil = paste0("L", seq_len(n_lentil)),
    Wheat = paste0("W", seq_len(n_wheat)),
    Rep_combo = seq_len(reps),
    KEEP.OUT.ATTRS = FALSE
  )
  dat$Row <- seq_len(nrow(dat))
  dat$Col <- rep(seq_len(max(2, n_wheat * reps)), length.out = nrow(dat))
  dat$Block <- rep(seq_len(n_lentil), length.out = nrow(dat))
  dat$Type <- "Treatment"
  dat$Combo <- paste(dat$Lentil, dat$Wheat, sep = "-")
  dat$Total_Aboveground_N_g <-
    as.numeric(factor(dat$Wheat)) +
    as.numeric(factor(dat$Lentil)) / 10 +
    as.numeric(dat$Rep_combo) / 100
  dat
}

test_that("derive_wheat_rotation_traits computes N and C traits correctly", {
  dat <- data.frame(
    Grain_Biomass_g = 10,
    Straw_Biomass_g = 20,
    Seed_N_pct = 3,
    Straw_N_pct = 1,
    Seed_C_pct = 45,
    Straw_C_pct = 40
  )

  out <- derive_wheat_rotation_traits(dat)

  expect_equal(out$Total_Aboveground_Biomass_g, 30)
  expect_equal(out$Grain_N_Content_g, 0.3)
  expect_equal(out$Straw_N_Content_g, 0.2)
  expect_equal(out$Total_Aboveground_N_g, 0.5)
  expect_equal(out$N_Harvest_Index, 0.6)
})

test_that("model_wheat_rotational_value returns expected wheat-level output", {
  skip_if_not_installed("lme4")
  wheat_subsample <- make_wheat_rotation_test_data()

  result <- suppressWarnings(
    model_wheat_rotational_value(
      data = wheat_subsample,
      trait = "Total_Aboveground_N_g",
      env_col = "ENV",
      facet_col = "Facet",
      prev_gen_col = "Lentil",
      curr_gen_col = "Wheat",
      spatial_cols = c("Row", "Col"),
      fixed_effect_cols = c("Rep_combo"),
      include_checks = FALSE,
      method = "lme4",
      validate_design = FALSE
    )
  )

  expect_true("rotational_values" %in% names(result))
  expect_true(all(c(
    "ENV",
    "Facet",
    "Current_Genotype",
    "Wheat_Rotational_Value",
    "Total_Correction"
  ) %in% names(result$rotational_values)))
  expect_equal(nrow(result$rotational_values), dplyr::n_distinct(wheat_subsample$Wheat))
})

test_that("model_gxg_rotational_value returns observed combinations only", {
  skip_if_not_installed("lme4")
  wheat_subsample <- make_wheat_rotation_test_data()

  result <- suppressWarnings(
    model_gxg_rotational_value(
      data = wheat_subsample,
      trait = "Total_Aboveground_N_g",
      env_col = "ENV",
      facet_col = "Facet",
      prev_gen_col = "Lentil",
      curr_gen_col = "Wheat",
      combo_col = "Combo",
      spatial_cols = c("Row", "Col"),
      random_effect_cols = c("Block"),
      include_checks = FALSE,
      validate_design = FALSE
    )
  )

  observed <- wheat_subsample |>
    dplyr::distinct(ENV, Facet, Lentil, Wheat, Combo)

  output <- result$gxg_values |>
    dplyr::distinct(
      ENV,
      Facet,
      Previous_Genotype,
      Current_Genotype,
      Combo
    )

  joined <- dplyr::inner_join(
    output,
    observed,
    by = c(
      "ENV",
      "Facet",
      "Previous_Genotype" = "Lentil",
      "Current_Genotype" = "Wheat",
      "Combo"
    )
  )
  expect_equal(nrow(output), nrow(joined))
  expect_true("Status" %in% names(result$gxg_values))
})

test_that("build_wheat_rotation_index returns one score per wheat genotype", {
  wheat_rotational_values <- data.frame(
    ENV = rep("E1", 6),
    Facet = rep("F1", 6),
    Current_Genotype = rep(c("W1", "W2"), each = 3),
    Trait = rep(c("Total_Aboveground_N_g", "Grain_N_Content_g", "N_Harvest_Index"), times = 2),
    Wheat_Rotational_Value = c(1, 2, 3, -1, -2, -3)
  )

  index <- build_wheat_rotation_index(
    rotational_values = wheat_rotational_values,
    trait_weights = c(
      Total_Aboveground_N_g = 1,
      Grain_N_Content_g = 1,
      N_Harvest_Index = 1
    )
  )

  expect_true("Wheat_Rotation_Index" %in% names(index))
  expect_true("N_Traits_Used" %in% names(index))
  expect_equal(nrow(index), 2)
})
