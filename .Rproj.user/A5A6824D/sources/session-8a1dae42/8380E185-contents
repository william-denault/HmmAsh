testthat::test_that("Gaussian-null stress fits are exact, positive, and crash-free", {
  testthat::skip_on_cran()

  # Increased from 100 to 150 routine replicates; extended mode runs 500.
  n_rep <- ash_hmm_stress_reps(standard = 150L, extended = 500L)
  errors <- rep("", n_rep)
  collapsed <- exact_null <- positive_centers <- logical(n_rep)

  for (i in seq_len(n_rep)) {
    set.seed(7000L + i)
    n <- if (i %% 2L) 128L else 192L
    se <- if (i %% 2L) rep(1, n) else stats::runif(n, 0.5, 1.5)
    y <- stats::rnorm(n, 0, se)

    fit <- tryCatch(
      fit_ash_hmm(
        y, se,
        half_grid = 10L,
        topology = "full",
        null_state = "pointmass",
        maxiter = 10L,
        tolerance = 1e-6),
      error = function(e) e)

    if (inherits(fit, "error")) {
      errors[i] <- sprintf("replicate %d: %s", i, conditionMessage(fit))
      next
    }

    collapsed[i] <- isTRUE(fit$model_selection$collapsed_to_null)
    positive_centers[i] <- length(fit$fitted$mu) == 1L ||
      all(fit$fitted$mu[-1L] > 0)
    exact_null[i] <-
      identical(as.numeric(fit$fitted$mu), 0) &&
      identical(fit$fitted$null_state, "pointmass") &&
      identical(as.numeric(fit$fitted$transition), 1) &&
      identical(as.numeric(fit$fitted$init_prob), 1) &&
      all(fit$state_probability == 1) &&
      all(fit$posterior$mean == 0) &&
      fit$fitted$mixture_weight[1L, 1L] == 1 &&
      all(fit$fitted$mixture_weight[1L, -1L] == 0) &&
      fit$step_selection$step_count == 1L &&
      fit$step_selection$change_count == 0L
  }

  testthat::expect_false(
    any(nzchar(errors)),
    info = ash_hmm_error_summary(errors))
  testthat::expect_true(
    mean(collapsed) >= 0.99,
    info = sprintf("Strict null selected in %d/%d replicates; require >=99%%.",
                   sum(collapsed), n_rep))
  testthat::expect_true(
    mean(exact_null) >= 0.99,
    info = sprintf(paste0("Exact one-state null returned in %d/%d replicates; ",
                         "require >=99%%."), sum(exact_null), n_rep))
  testthat::expect_true(
    all(positive_centers),
    info = sprintf("Positive-center constraint held in %d/%d replicates.",
                   sum(positive_centers), n_rep))
})

testthat::test_that("adaptive null remains numerically valid on pure noise", {
  testthat::skip_on_cran()

  # This smaller companion test exercises the optional adaptive-null branch.
  n_rep <- ash_hmm_stress_reps(standard = 30L, extended = 100L)
  errors <- rep("", n_rep)
  valid <- logical(n_rep)

  for (i in seq_len(n_rep)) {
    set.seed(12000L + i)
    n <- 128L
    se <- stats::runif(n, 0.5, 1.5)
    fit <- tryCatch(
      fit_ash_hmm(
        stats::rnorm(n, 0, se), se,
        half_grid = 10L,
        topology = "full",
        null_state = "adaptive",
        maxiter = 10L,
        tolerance = 1e-6),
      error = function(e) e)

    if (inherits(fit, "error")) {
      errors[i] <- sprintf("replicate %d: %s", i, conditionMessage(fit))
      next
    }

    valid[i] <-
      all(is.finite(fit$posterior$mean)) &&
      all(is.finite(fit$state_probability)) &&
      all(abs(rowSums(fit$state_probability) - 1) < 1e-8) &&
      all(abs(rowSums(fit$fitted$transition) - 1) < 1e-8) &&
      all(abs(rowSums(fit$fitted$mixture_weight) - 1) < 1e-8) &&
      (length(fit$fitted$mu) == 1L || all(fit$fitted$mu[-1L] > 0))
  }

  testthat::expect_false(
    any(nzchar(errors)),
    info = ash_hmm_error_summary(errors))
  testthat::expect_true(
    all(valid),
    info = sprintf("Valid adaptive-null result in %d/%d replicates.",
                   sum(valid), n_rep))
})
