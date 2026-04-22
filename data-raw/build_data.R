# data-raw/build_data.R
# Build internal lookup data for icd10am.charlson.
# Run this script interactively from the package root to regenerate R/sysdata.rda.
#
# Source files (edit these when ICD-10-AM editions are updated):
#   charlson_icd10am_codes.xlsx  — ICD-10-AM to Charlson category mapping
#   charlson_weights.xlsx        — Quan and Sundararajan weights per category

library(readxl)

char_elix_codes <- as.data.frame(
  readxl::read_excel("data-raw/charlson_icd10am_codes.xlsx",
                     sheet = "charlson_icd10am_codes",
                     col_types = c("text", "text", "text"))
)

charlson_weights <- as.data.frame(
  readxl::read_excel("data-raw/charlson_weights.xlsx",
                     sheet = "charlson_weights",
                     col_types = c("text", "text", "numeric", "numeric"))
)

usethis::use_data(char_elix_codes, charlson_weights, internal = TRUE, overwrite = TRUE)
