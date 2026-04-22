dfbase <- icd10am.utils::dfbase
dfdiag <- icd10am.utils::dfdiag
dfwide <- icd10am.utils::dfwide

test_that("classify_charlson returns expected columns in long format", {
  result <- classify_charlson(
    df      = dfbase,
    id      = recordID,
    format  = "long",
    df_diag = dfdiag,
    diag    = "diag",
    diagno  = "diagno"
  )

  expect_true("charlson_score" %in% names(result))
  expect_equal(nrow(result), nrow(dfbase))
  expect_true(all(result$charlson_score >= 0))
})

test_that("classify_charlson wide format produces same scores as long format", {
  # Derive long inputs from dfwide so IDs match
  df_episodes <- dfwide[, c("recordID", "ageyears", "agemonths")]

  diag_cols <- grep("^diag[0-9]+$", names(dfwide), value = TRUE)
  df_diag_long <- tidyr::pivot_longer(
    dfwide[, c("recordID", diag_cols)],
    cols      = dplyr::all_of(diag_cols),
    names_to  = "diagno",
    values_to = "diag"
  ) |>
    dplyr::mutate(diagno = as.integer(gsub("diag", "", diagno))) |>
    dplyr::filter(!is.na(diag), diag != "")

  result_long <- classify_charlson(
    df      = df_episodes,
    id      = recordID,
    format  = "long",
    df_diag = df_diag_long,
    diag    = "diag",
    diagno  = "diagno"
  )

  result_wide <- classify_charlson(
    df          = dfwide,
    id          = recordID,
    format      = "wide",
    diag_prefix = "diag"
  )

  scores_long <- result_long[order(result_long$recordID), "charlson_score"][[1]]
  scores_wide <- result_wide[order(result_wide$recordID), "charlson_score"][[1]]

  expect_equal(scores_long, scores_wide)
})

test_that("classify_charlson gives score 0 for episodes with no matching diagnoses", {
  empty_diag <- data.frame(
    recordID = dfbase$recordID[1],
    diag     = "ZZZZZ",
    diagno   = 1L
  )

  result <- classify_charlson(
    df = dfbase[1, ], id = recordID, format = "long",
    df_diag = empty_diag, diag = "diag", diagno = "diagno"
  )

  expect_equal(result$charlson_score, 0)
})

test_that("classify_charlson produces non-zero scores for known comorbid episodes", {
  result <- classify_charlson(
    df = dfbase, id = recordID, format = "long",
    df_diag = dfdiag, diag = "diag", diagno = "diagno"
  )

  expect_true(any(result$charlson_score > 0))
})

test_that("classify_charlson C01-C17 flag columns are present and binary", {
  result <- classify_charlson(
    df = dfbase, id = recordID, format = "long",
    df_diag = dfdiag, diag = "diag", diagno = "diagno"
  )

  flag_cols <- grep("^C[0-9]{2}$", names(result), value = TRUE)
  expect_gte(length(flag_cols), 1L)
  expect_true(all(unlist(result[, flag_cols]) %in% c(0L, 1L, NA)))
})

test_that("classify_charlson Sundararajan method runs without error", {
  result <- classify_charlson(
    df = dfbase, id = recordID, format = "long",
    df_diag = dfdiag, diag = "diag", diagno = "diagno",
    method = "Sundararajan"
  )

  expect_true("charlson_score" %in% names(result))
  expect_true(all(result$charlson_score >= 0))
})
