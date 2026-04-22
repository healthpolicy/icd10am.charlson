#' Classify episodes using the Charlson Comorbidity Index
#'
#' @description
#' Applies the Charlson Comorbidity Index to admitted patient episodes coded
#' with ICD-10-AM. Each episode receives a total Charlson score and binary
#' flags for each of the 17 Charlson comorbidity categories.
#'
#' Both wide-format (one row per episode, diagnoses in separate columns) and
#' long-format (separate diagnosis data frame) inputs are supported via the
#' `format` argument.
#'
#' @param df A data frame with one row per episode.
#' @param id The episode identifier column (unquoted name).
#' @param format Input format: `"long"` (default) or `"wide"`.
#' @param df_diag Long-format diagnosis data frame with one row per diagnosis.
#'   Required when `format = "long"`. Must contain the episode ID column, a
#'   diagnosis code column, and a diagnosis sequence number column.
#' @param diag Name of the diagnosis code column in `df_diag` (quoted string).
#'   Default: `"diag"`. Used only when `format = "long"`.
#' @param diagno Name of the diagnosis sequence number column in `df_diag`
#'   (quoted string). Default: `"diagno"`. Used only when `format = "long"`.
#' @param diag_prefix Character prefix shared by all diagnosis columns in `df`
#'   (e.g. `"diag"` matches `diag1`, `diag2`, ...). Default: `"diag"`. Used
#'   only when `format = "wide"`.
#' @param method Scoring method. `"Quan"` (default) applies the Quan et al.
#'   (2005) updated weights; `"Sundararajan"` applies the original weights.
#'
#' @return The input `df` with additional columns appended:
#'   \describe{
#'     \item{`charlson_score`}{Total Charlson comorbidity score (numeric).}
#'     \item{`S01` ... `S17`}{Binary flags (0/1) for each Charlson comorbidity
#'       category.}
#'   }
#'   Episodes with no matching diagnoses receive a score of 0.
#'
#' @references
#' Charlson ME, Pompei P, Ales KL, MacKenzie CR (1987). A new method of
#' classifying prognostic comorbidity in longitudinal studies: development and
#' validation. *Journal of Chronic Diseases*, 40(5), 373--383.
#' \doi{10.1016/0021-9681(87)90171-8}
#'
#' Quan H, Sundararajan V, Halfon P, et al. (2005). Coding algorithms for
#' defining comorbidities in ICD-9-CM and ICD-10 administrative data.
#' *Medical Care*, 43(11), 1130--1139.
#' \doi{10.1097/01.mlr.0000182534.19832.83}
#'
#' @export
#'
#' @importFrom dplyr all_of arrange distinct filter first group_by if_else
#'   left_join mutate rename select starts_with summarise across
#' @importFrom tidyr pivot_wider
#' @importFrom rlang .data abort as_name ensym sym
#' @importFrom stats setNames
#' @importFrom icd10am.utils .to_long
#'
#' @examples
#' dfbase <- icd10am.utils::dfbase
#' dfdiag <- icd10am.utils::dfdiag
#'
#' # Long format
#' classify_charlson(
#'   df      = dfbase,
#'   id      = recordID,
#'   format  = "long",
#'   df_diag = dfdiag,
#'   diag    = "diag",
#'   diagno  = "diagno",
#'   method  = "Quan"
#' )
#'
#' # Wide format
#' dfwide <- icd10am.utils::dfwide
#' classify_charlson(
#'   df          = dfwide,
#'   id          = recordID,
#'   format      = "wide",
#'   diag_prefix = "diag",
#'   method      = "Quan"
#' )
classify_charlson <- function(
    df,
    id,
    format      = c("long", "wide"),
    df_diag     = NULL,
    diag        = "diag",
    diagno      = "diagno",
    diag_prefix = "diag",
    method      = c("Quan", "Sundararajan")
) {
  format <- match.arg(format)
  method <- match.arg(method)
  id_str <- rlang::as_name(rlang::ensym(id))

  # --- Normalise to long format ---
  if (format == "wide") {
    converted <- .to_long(df, id = !!rlang::sym(id_str), diag_prefix = diag_prefix)
    df_diag   <- converted$df_diag
    diag      <- "diag"
    diagno    <- "diagno"
  }

  if (is.null(df_diag)) {
    rlang::abort("`df_diag` must be supplied when `format = 'long'`.")
  }

  # --- Rename to standard internal names ---
  df_diag_std <- df_diag |>
    dplyr::rename(
      episode_id = dplyr::all_of(id_str),
      diag       = dplyr::all_of(diag),
      diagno     = dplyr::all_of(diagno)
    )

  # --- Match diagnoses to Charlson codes ---
  # Some ICD-10-AM codes map to multiple Charlson categories; many-to-many is expected.
  df_matched <- df_diag_std |>
    dplyr::left_join(char_elix_codes, by = "diag",
                     relationship = "many-to-many") |>
    dplyr::filter(!is.na(.data$charlson_q), !is.na(.data$charlson_qx)) |>
    dplyr::distinct() |>
    dplyr::left_join(
      charlson_weights |> dplyr::rename(charlson_q = "code"),
      by = "charlson_q"
    ) |>
    dplyr::arrange(.data$charlson_q)

  # --- Compute score ---
  weight_col <- if (method == "Quan") "weight_q" else "weight_s"

  df_score <- df_matched |>
    dplyr::group_by(.data$episode_id) |>
    dplyr::summarise(
      charlson_score = sum(.data[[weight_col]], na.rm = TRUE),
      .groups = "drop"
    )

  # --- Compute condition flags ---
  df_flags <- df_matched |>
    dplyr::select("episode_id", "charlson_q", "charlson_qx") |>
    dplyr::distinct() |>
    tidyr::pivot_wider(
      names_from  = "charlson_q",
      values_from = "charlson_qx",
      values_fn   = list(charlson_qx = dplyr::first)
    ) |>
    dplyr::mutate(
      dplyr::across(dplyr::starts_with("C"), ~ dplyr::if_else(is.na(.x), 0L, 1L))
    )

  # --- Join back to episode data ---
  df |>
    dplyr::left_join(df_score, by = stats::setNames("episode_id", id_str)) |>
    dplyr::left_join(df_flags, by = stats::setNames("episode_id", id_str)) |>
    dplyr::mutate(
      charlson_score = dplyr::if_else(is.na(.data$charlson_score), 0, .data$charlson_score)
    )
}
