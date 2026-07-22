# Shared controls for the ash-HMM stress tests.
#
# Routine stress suite:
#   testthat::test_dir("tests/testthat")
#
# Extended suite:
#   Sys.setenv(ASH_HMM_EXTENDED_TESTS = "true")
#   testthat::test_dir("tests/testthat")
#
# Exact per-scenario override, useful while debugging:
#   Sys.setenv(ASH_HMM_STRESS_REPS = "5")

ash_hmm_stress_reps <- function(standard, extended) {
  override <- suppressWarnings(as.integer(Sys.getenv("ASH_HMM_STRESS_REPS", "")))
  if (length(override) == 1L && !is.na(override) && override > 0L) {
    return(override)
  }
  extended_requested <- tolower(Sys.getenv(
    "ASH_HMM_EXTENDED_TESTS", "false")) %in% c("1", "true", "yes")
  as.integer(if (extended_requested) extended else standard)
}

ash_hmm_error_summary <- function(errors, maximum = 5L) {
  errors <- errors[nzchar(errors)]
  if (!length(errors)) return("")
  paste(utils::head(errors, maximum), collapse = "\n")
}
