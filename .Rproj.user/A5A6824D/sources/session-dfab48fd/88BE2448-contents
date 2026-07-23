# Regression and stress tests for the positive-state, null-safe ash-HMM.

source("utils.R")
source("HmmAsh.R")

assert_close <- function(x, y, tolerance = 1e-8, label = "values") {
  error <- max(abs(x - y))
  if (length(x) != length(y) || any(!is.finite(x)) ||
      any(!is.finite(y)) || error > tolerance) {
    stop(sprintf("Mismatch in %s (maximum error %.3g).", label, error),
         call. = FALSE)
  }
}

stopifnot(identical(ash_hmm_version(), "1.4.0-documented-nonnegative"))
stopifnot(is.null(formals(fit_ash_hmm)$mu))
stopifnot(isTRUE(eval(formals(fit_ash_hmm)$nonnegative_state_means)))
stopifnot(is.null(formals(fit_ash_hmm)$positive_state_means))
stopifnot(identical(eval(formals(fit_ash_hmm)$null_state),
                    c("pointmass", "adaptive")))

# 1. The automatic grid is zero followed only by strictly positive candidates.
grid <- ash_mu_grid(c(-2, 0, 3), half_grid = 12L, positive_only = TRUE)
stopifnot(grid[1L] == 0, all(grid[-1L] > 0), !anyDuplicated(grid))

# 2. The constrained center M-step cannot cross zero.
statistics <- list(occupancy = c(10, 10),
                   mean_numerator = c(0, -100),
                   mean_precision = c(10, 10))
update <- .ash_hmm_update_means(
  mu = c(0, 1), anchor = c(0, 1), statistics = statistics,
  learn = TRUE, fixed = 1L, minimum_count = 0,
  eligible = c(FALSE, TRUE), damping = 1,
  minimum = c(0, 1e-8), bounds = "none")
stopifnot(update$mu[1L] == 0, update$mu[2L] >= 1e-8)

# 3. Stress test: independent N(0,1) sequences must not crash. The extended-BIC
# comparison should select the strict null with high frequency, and every
# selected-null fit must be internally exact.
null_selected <- logical(20L)
for (seed in seq_len(20L)) {
  set.seed(seed)
  y0 <- rnorm(128L)
  f0 <- fit_ash_hmm(
    y0, rep(1, 128L), half_grid = 10L, topology = "full",
    maxiter = 12L, tolerance = 1e-6, verbose = FALSE)
  null_selected[seed] <- isTRUE(f0$model_selection$collapsed_to_null)
  if (null_selected[seed]) {
    stopifnot(length(f0$fitted$mu) == 1L,
              f0$fitted$mu == 0,
              all(f0$posterior$mean == 0),
              all(f0$state_probability == 1),
              f0$step_selection$step_count == 1L,
              f0$step_selection$change_count == 0L)
  }
}
if (mean(null_selected) < 0.9) {
  stop(sprintf("Strict null selected in only %.1f%% of Gaussian-null runs.",
               100 * mean(null_selected)), call. = FALSE)
}

# 4. Strong positive piecewise-constant signal: retained centers are strictly
# positive, except the exact null, and the log(T)-penalized decoder recovers the
# five true segments in this reproducible example.
set.seed(4001)
truth <- c(rep(0, 80), rep(2, 80), rep(4, 80), rep(1.5, 80), rep(0, 80))
se <- rep(0.35, length(truth))
y <- rnorm(length(truth), truth, se)
fit_steps <- fit_ash_hmm(
  y, se, half_grid = 24L, grid_shape = 2,
  topology = "full", maxiter = 40L, tolerance = 1e-7,
  prune_max_loglik_loss = 0.1, verbose = FALSE)
stopifnot(!fit_steps$model_selection$collapsed_to_null,
          fit_steps$fitted$mu[1L] == 0,
          all(fit_steps$fitted$mu[-1L] > 0),
          fit_steps$step_selection$step_count == 5L,
          fit_steps$step_selection$change_count == 4L)
assert_close(rowSums(fit_steps$fitted$transition),
             rep(1, length(fit_steps$fitted$mu)),
             label = "transition row sums")

# 5. The null prior is an exact Dirac mass by default, while an adaptive ash
# mixture centered at zero can be requested explicitly.
stopifnot(fit_steps$fitted$null_state == "pointmass",
          1L %in% fit_steps$fitted$fixed_pointmass_states,
          fit_steps$fitted$mixture_weight[1L, 1L] == 1,
          all(fit_steps$fitted$mixture_weight[1L, -1L] == 0))
adaptive_null_fit <- fit_ash_hmm(
  y, se, half_grid = 12L, topology = "full",
  null_state = "adaptive", null_model_selection = "none",
  maxiter = 6L, tolerance = 1e-6, verbose = FALSE)
stopifnot(adaptive_null_fit$fitted$null_state == "adaptive",
          !(1L %in% adaptive_null_fit$fitted$fixed_pointmass_states),
          abs(sum(adaptive_null_fit$fitted$mixture_weight[1L, ]) - 1) < 1e-10)

# 6. Signed means remain available explicitly for backward compatibility.
set.seed(55)
signed_truth <- c(rep(-1.5, 60), rep(0, 60), rep(1.5, 60))
signed_fit <- fit_ash_hmm(
  rnorm(length(signed_truth), signed_truth, 0.4),
  rep(0.4, length(signed_truth)),
  nonnegative_state_means = FALSE,
  fixed_pointmass_states = integer(),
  null_model_selection = "none",
  topology = "full", half_grid = 16L,
  maxiter = 10L, tolerance = 1e-6, verbose = FALSE)
stopifnot(any(signed_fit$fitted$mu < 0), any(signed_fit$fitted$mu > 0))

# The previous argument remains a compatibility alias, with conflict checking.
alias_fit <- fit_ash_hmm(
  rnorm(60L), rep(1, 60L),
  positive_state_means = FALSE,
  null_model_selection = "none",
  half_grid = 6L, maxiter = 1L,
  learn_state_means = FALSE, prune_states = FALSE,
  verbose = FALSE)
stopifnot(!alias_fit$fitted$nonnegative_state_means)
conflict <- try(fit_ash_hmm(
  rnorm(20L), rep(1, 20L),
  nonnegative_state_means = TRUE,
  positive_state_means = FALSE), silent = TRUE)
stopifnot(inherits(conflict, "try-error"))

message(sprintf(
  paste0("All positive-step ash-HMM tests passed; strict null selected in ",
         "%d/20 Gaussian-noise runs, recovered %d steps."),
  sum(null_selected), fit_steps$step_selection$step_count))
