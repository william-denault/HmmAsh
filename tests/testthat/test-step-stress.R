ash_hmm_step_stress_fit <- function(seed, amplitude, noise) {
  set.seed(seed)
  truth <- c(rep(0, 60), rep(amplitude, 60),
             rep(2 * amplitude, 60), rep(0, 60))
  se <- rep(noise, length(truth))
  fit <- fit_ash_hmm(
    stats::rnorm(length(truth), truth, se), se,
    half_grid = 14L,
    grid_shape = 2,
    topology = "full",
    null_state = "pointmass",
    maxiter = 20L,
    tolerance = 1e-6,
    prune_max_loglik_loss = 0.1,
    verbose = FALSE)
  list(fit = fit, truth = truth)
}

testthat::test_that("strong four-step signals recover every step", {
  testthat::skip_on_cran()

  # Increased from 20 to 50 routine replicates; extended mode runs 200.
  n_rep <- ash_hmm_stress_reps(standard = 50L, extended = 200L)
  errors <- rep("", n_rep)
  exact <- valid <- logical(n_rep)

  for (i in seq_len(n_rep)) {
    result <- tryCatch(
      ash_hmm_step_stress_fit(8000L + i, amplitude = 2, noise = 0.4),
      error = function(e) e)
    if (inherits(result, "error")) {
      errors[i] <- sprintf("replicate %d: %s", i, conditionMessage(result))
      next
    }
    fit <- result$fit
    exact[i] <-
      !isTRUE(fit$model_selection$collapsed_to_null) &&
      fit$step_selection$step_count == 4L &&
      fit$step_selection$change_count == 3L
    valid[i] <-
      fit$fitted$mu[1L] == 0 &&
      all(fit$fitted$mu[-1L] > 0) &&
      fit$fitted$effect_support == "nonnegative" &&
      all(abs(rowSums(fit$fitted$transition) - 1) < 1e-8) &&
      all(is.finite(fit$posterior$mean)) &&
      all(fit$posterior$mean >= 0) &&
      all(fit$posterior$probability_ge_zero == 1)
  }

  testthat::expect_false(
    any(nzchar(errors)),
    info = ash_hmm_error_summary(errors))
  testthat::expect_true(
    mean(exact) >= 0.98,
    info = sprintf(paste0("Exact four-step recovery in %d/%d replicates; ",
                         "require >=98%%."), sum(exact), n_rep))
  testthat::expect_true(
    all(valid),
    info = sprintf("Valid constrained fit in %d/%d replicates.",
                   sum(valid), n_rep))
})

testthat::test_that("sub-threshold steps are not over-resolved", {
  testthat::skip_on_cran()

  # Increased from 20 to 50 routine replicates; extended mode runs 200.
  n_rep <- ash_hmm_stress_reps(standard = 50L, extended = 200L)
  errors <- rep("", n_rep)
  collapsed <- exact_null <- logical(n_rep)

  for (i in seq_len(n_rep)) {
    result <- tryCatch(
      ash_hmm_step_stress_fit(9000L + i, amplitude = 0.15, noise = 1),
      error = function(e) e)
    if (inherits(result, "error")) {
      errors[i] <- sprintf("replicate %d: %s", i, conditionMessage(result))
      next
    }
    fit <- result$fit
    collapsed[i] <- isTRUE(fit$model_selection$collapsed_to_null)
    exact_null[i] <-
      length(fit$fitted$mu) == 1L &&
      fit$fitted$mu == 0 &&
      fit$step_selection$step_count == 1L &&
      fit$step_selection$change_count == 0L &&
      all(fit$posterior$mean == 0)
  }

  testthat::expect_false(
    any(nzchar(errors)),
    info = ash_hmm_error_summary(errors))
  testthat::expect_true(
    mean(collapsed) >= 0.98,
    info = sprintf(paste0("Strict null selected in %d/%d weak replicates; ",
                         "require >=98%%."), sum(collapsed), n_rep))
  testthat::expect_true(
    mean(exact_null) >= 0.98,
    info = sprintf(paste0("Exact null output in %d/%d weak replicates; ",
                         "require >=98%%."), sum(exact_null), n_rep))
})

testthat::test_that("recurrent levels are counted as separate contiguous steps", {
  testthat::skip_on_cran()

  # This specifically checks that five steps can use only four occupied levels.
  n_rep <- ash_hmm_stress_reps(standard = 20L, extended = 75L)
  errors <- rep("", n_rep)
  exact <- distinct_counts <- logical(n_rep)

  for (i in seq_len(n_rep)) {
    set.seed(16000L + i)
    truth <- c(rep(0, 80), rep(2, 80), rep(4, 80),
               rep(1.5, 80), rep(0, 80))
    se <- rep(0.35, length(truth))
    fit <- tryCatch(
      fit_ash_hmm(
        stats::rnorm(length(truth), truth, se), se,
        half_grid = 24L,
        grid_shape = 2,
        topology = "full",
        null_state = "pointmass",
        maxiter = 40L,
        tolerance = 1e-7,
        prune_max_loglik_loss = 0.1,
        verbose = FALSE),
      error = function(e) e)

    if (inherits(fit, "error")) {
      errors[i] <- sprintf("replicate %d: %s", i, conditionMessage(fit))
      next
    }
    exact[i] <-
      fit$step_selection$step_count == 5L &&
      fit$step_selection$change_count == 4L
    distinct_counts[i] <-
      fit$step_selection$occupied_state_count <
      fit$step_selection$step_count
    if (any(fit$posterior$mean < 0)) {
      errors[i] <- sprintf("replicate %d returned a negative posterior mean", i)
    }
  }

  testthat::expect_false(
    any(nzchar(errors)),
    info = ash_hmm_error_summary(errors))
  testthat::expect_true(
    mean(exact) >= 0.90,
    info = sprintf(paste0("Exact five-step recovery in %d/%d replicates; ",
                         "require >=90%%."), sum(exact), n_rep))
  testthat::expect_true(
    mean(distinct_counts) >= 0.90,
    info = sprintf(paste0("Distinct state/step counts in %d/%d replicates; ",
                         "require >=90%%."), sum(distinct_counts), n_rep))
})
