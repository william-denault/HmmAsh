.ash_hmm_test_internal <- function(name) {
  if (exists(name, mode = "function", inherits = TRUE)) {
    return(get(name, mode = "function", inherits = TRUE))
  }
  getFromNamespace(name, "ashHMM")
}

testthat::test_that("truncated-normal emissions equal direct convolution", {
  y <- c(-2, 0.3, 3)
  se <- c(0.7, 1.2, 0.5)
  mu <- 0.8
  tau <- 1.4

  component_log_emission <- .ash_hmm_test_internal(
    ".ash_hmm_component_log_emission")
  formula <- exp(component_log_emission(
    y, se, mu, tau, "nonnegative"))
  normalizer <- stats::pnorm(mu / tau)
  numerical <- vapply(seq_along(y), function(i) {
    stats::integrate(function(beta) {
      stats::dnorm(y[i], beta, se[i]) *
        stats::dnorm(beta, mu, tau) / normalizer
    }, lower = 0, upper = Inf, rel.tol = 1e-11)$value
  }, numeric(1L))

  testthat::expect_equal(formula, numerical, tolerance = 1e-10)
})

testthat::test_that("one-sided posterior moments equal numerical integration", {
  y <- c(-2, 0.3, 3)
  se <- c(0.7, 1.2, 0.5)
  mu <- 0.8
  tau <- 1.4
  prior_sd <- c(0, tau)
  rho <- matrix(c(0, 1), 1L, 2L)
  hmm_log_emission <- .ash_hmm_test_internal(".ash_hmm_log_emission")
  posterior_summary <- .ash_hmm_test_internal(
    ".ash_hmm_posterior_summary")
  log_emission <- hmm_log_emission(
    y, se, mu, prior_sd, rho, "nonnegative")
  posterior <- posterior_summary(
    y, se, mu, prior_sd, rho, log_emission,
    matrix(1, length(y), 1L), "nonnegative")

  normalizer <- stats::pnorm(mu / tau)
  marginal <- exp(log_emission[, 1L])
  numerical_mean <- numerical_second <- numeric(length(y))
  for (i in seq_along(y)) {
    numerical_mean[i] <- stats::integrate(function(beta) {
      beta * stats::dnorm(y[i], beta, se[i]) *
        stats::dnorm(beta, mu, tau) / normalizer
    }, 0, Inf, rel.tol = 1e-11)$value / marginal[i]
    numerical_second[i] <- stats::integrate(function(beta) {
      beta^2 * stats::dnorm(y[i], beta, se[i]) *
        stats::dnorm(beta, mu, tau) / normalizer
    }, 0, Inf, rel.tol = 1e-11)$value / marginal[i]
  }

  testthat::expect_equal(posterior$mean, numerical_mean, tolerance = 1e-10)
  testthat::expect_equal(
    posterior$sd^2 + posterior$mean^2,
    numerical_second, tolerance = 1e-10)
  testthat::expect_true(all(posterior$mean >= 0))
  testthat::expect_true(all(posterior$probability_ge_zero == 1))
  testthat::expect_true(all(posterior$probability_le_zero == 0))
})

testthat::test_that("main one-sided option guarantees nonnegative posterior effects", {
  set.seed(123)
  truth <- c(rep(0, 60), rep(2, 80), rep(0, 60))
  se <- rep(1, length(truth))
  fit <- fit_ash_hmm(
    stats::rnorm(length(truth), truth, se), se,
    nonnegative_state_means = TRUE,
    half_grid = 12L,
    topology = "full",
    maxiter = 10L,
    null_model_selection = "none",
    verbose = FALSE)

  testthat::expect_identical(fit$fitted$effect_support, "nonnegative")
  testthat::expect_true(all(fit$fitted$mu >= 0))
  testthat::expect_true(all(fit$posterior$mean >= 0))
  testthat::expect_true(all(fit$posterior$probability_ge_zero == 1))
  testthat::expect_equal(
    fit$posterior$lfsr, fit$posterior$probability_zero,
    tolerance = 1e-12)
})

testthat::test_that("support and center options cannot silently conflict", {
  y <- c(-1, 0, 1)
  se <- rep(1, length(y))
  testthat::expect_error(
    fit_ash_hmm(
      y, se,
      nonnegative_state_means = TRUE,
      effect_support = "real",
      verbose = FALSE),
    "conflicts")
  testthat::expect_error(
    fit_ash_hmm(
      y, se,
      nonnegative_state_means = FALSE,
      effect_support = "nonnegative",
      verbose = FALSE),
    "conflicts")
})

testthat::test_that("signed mode retains ordinary real-line effects", {
  set.seed(456)
  truth <- c(rep(0, 40), rep(-2, 80), rep(0, 40))
  se <- rep(0.5, length(truth))
  fit <- fit_ash_hmm(
    stats::rnorm(length(truth), truth, se), se,
    nonnegative_state_means = FALSE,
    half_grid = 12L,
    maxiter = 10L,
    null_model_selection = "none",
    verbose = FALSE)

  testthat::expect_identical(fit$fitted$effect_support, "real")
  testthat::expect_true(any(fit$fitted$mu < 0))
  testthat::expect_true(any(fit$posterior$mean < 0))
})

testthat::test_that("positive_state_means remains a working compatibility alias", {
  set.seed(457)
  y <- stats::rnorm(80)
  se <- rep(1, length(y))
  fit <- fit_ash_hmm(
    y, se,
    positive_state_means = TRUE,
    half_grid = 8L,
    maxiter = 3L,
    verbose = FALSE)

  testthat::expect_true(fit$fitted$nonnegative_state_means)
  testthat::expect_identical(fit$fitted$effect_support, "nonnegative")
  testthat::expect_true(all(fit$posterior$mean >= 0))
})

testthat::test_that("fixed-dimensional one-sided EM objective does not decrease", {
  set.seed(321)
  truth <- c(rep(0, 40), rep(1.5, 60), rep(3, 60), rep(0, 40))
  se <- stats::runif(length(truth), 0.4, 0.8)
  fit <- fit_ash_hmm(
    stats::rnorm(length(truth), truth, se), se,
    nonnegative_state_means = TRUE,
    half_grid = 12L,
    topology = "full",
    maxiter = 12L,
    prune_states = FALSE,
    null_model_selection = "none",
    verbose = FALSE)

  testthat::expect_true(all(diff(fit$history$objective) >= -1e-7))
  testthat::expect_true(all(fit$posterior$mean >= 0))
})
