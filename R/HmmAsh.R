# ash_hmm.R
#
# A dependency-free implementation of an adaptive-shrinkage hidden Markov
# model for heteroskedastic normal observations.
#
# Model:
#   y_t | beta_t ~ N(beta_t, se_t^2)
#   Q_t is a finite-state Markov chain with transition matrix A
#   beta_t | Q_t=m ~ sum_l rho[m,l] G_ml
#
# In signed mode, G_ml is N(mu[m], prior_sd[l]^2). In the default one-sided
# mode, continuous G_ml are those normals truncated to [0, Inf); the zero-scale
# component remains a point mass at mu[m].
#
# prior_sd[1] must be zero, so the first component is a point mass at mu[m].
# Candidate means and prior standard deviations may be generated from y and se.
# By default the null state is exactly delta_0; set null_state = "adaptive"
# to estimate an ash mixture centered at zero for that state as well.
# After a short fixed-grid warm-up, supported non-null means can be learned by
# the same EM sufficient statistics used for rho. Voronoi constraints and
# likelihood-checked state pruning prevent dense candidates from collapsing.


#' Fit an adaptive-shrinkage hidden Markov model
#'
#' @description
#' Fits a finite-state hidden Markov model to heteroskedastic normal effect
#' estimates. Each hidden state has an adaptive-shrinkage prior formed from a
#' point mass at its state center and a grid of normal components centered at
#' the same value. Transition probabilities, state-specific mixture weights,
#' and (optionally) supported state centers are learned by generalized
#' EM/Baum--Welch updates.
#'
#' @details
#' If `mu` is `NULL`, a dense grid of candidate state centers is constructed
#' from `y` and screened before the first forward-backward pass. With
#' `nonnegative_state_means = TRUE`, state 1 is fixed at zero, all non-null
#' state centers are constrained to be strictly positive, and every continuous
#' ash component is a normal distribution truncated to `[0, Inf)`. Consequently
#' the latent effects and posterior means are nonnegative as well. Set
#' `nonnegative_state_means = FALSE` to use ordinary Gaussian ash components on
#' the real line and construct a signed symmetric center grid.
#'
#' Low-occupancy or nearly duplicated states can be proposed for pruning.
#' A proposed deletion is accepted only after recomputing the complete HMM
#' marginal likelihood. After fitting, an optional information-criterion gate
#' can replace the fitted model by the exact all-zero HMM, and a separate
#' change-penalized decoder reports contiguous steps.
#'
#' @param y Numeric vector of noisy effect estimates.
#' @param se Numeric vector of known standard errors. It must have the same
#'   length as `y`, contain only finite values, and be strictly positive.
#' @param mu Optional numeric vector of initial state centers. State 1 is the
#'   zero/hub state. When `NULL`, centers are constructed and screened
#'   automatically. Supply `mu` when using state-sized custom initial values.
#' @param prior_sd Optional nondecreasing numeric vector of ash prior standard
#'   deviations. Its first value must be zero, representing a point mass. When
#'   `NULL`, a Stephens-style geometric scale grid is constructed automatically.
#' @param half_grid Number of nonnegative candidates in the automatic mean
#'   grid, including zero. Signed mode appends their nonzero negatives.
#' @param grid_shape Positive power-grid shape parameter. Values above one put
#'   relatively more candidates near the outer part of the range.
#' @param grid_expansion Positive multiplier applied to the empirical maximum
#'   absolute observation when constructing the automatic mean grid.
#' @param grid_max_abs Optional finite positive endpoint used instead of
#'   `max(abs(y))` when constructing the automatic grid.
#' @param nonnegative_state_means Logical; if `TRUE` (default), use one state
#'   centered at zero, constrain every non-null state center to be strictly
#'   positive, and truncate all continuous ash components below at zero. Thus
#'   both prior and posterior effects have support on `[0, Inf)`. If `FALSE`,
#'   allow signed state centers and effects on the real line.
#' @param positive_state_means Deprecated compatibility alias for
#'   `nonnegative_state_means`. Leave as `NULL` in new code. If both arguments
#'   are supplied, they must agree.
#' @param effect_support Either `"auto"` (default), `"nonnegative"`, or
#'   `"real"`. `"auto"` follows `nonnegative_state_means`. An explicit value
#'   must agree with that argument. The nonnegative model uses exact
#'   truncated-normal marginal emissions and posterior moments; it does not
#'   clip an unconstrained posterior.
#' @param positive_mean_floor Strict lower bound for learned non-null centers in
#'   nonnegative mode. The default is a scale-aware numerical value.
#' @param prefilter Logical controlling independent state-grid screening before
#'   HMM fitting. The default is `TRUE` for an automatic grid and `FALSE` for a
#'   supplied grid.
#' @param min_state_count,min_state_fraction A prefilter candidate is retained
#'   when its independent soft count is at least the larger of these absolute
#'   and fractional thresholds. The zero state is always retained.
#' @param screening_prior_sd Nonnegative extra prior standard deviation used
#'   only in the independent prefilter score.
#' @param screening_block_size Positive integer block size for memory-efficient
#'   prefilter calculations.
#' @param sequence_id Vector identifying independent sequences. Adjacent equal
#'   values belong to one sequence; transitions are not counted across changes
#'   in `sequence_id`.
#' @param forward_backward_engine Forward-backward implementation:
#'   `"auto"` (default) uses the fast scaled matrix recursion for dense
#'   transition masks and the sparse log-domain recursion for sparse masks;
#'   `"scaled"` requests the fast recursion with automatic log-domain fallback;
#'   `"log"` always uses the reference log-domain recursion.
#' @param topology Transition topology used when `transition_mask` is `NULL`:
#'   `"full"` allows all transitions; `"hub"` disallows direct transitions
#'   between distinct non-null states.
#' @param transition_mask Optional logical square matrix specifying allowed
#'   transitions between the supplied/retained states.
#' @param init_transition Optional initial row-stochastic transition matrix.
#' @param init_prob Optional initial state-probability vector.
#' @param init_rho Optional initial state-by-scale matrix of ash mixture weights.
#' @param stay_probability Initial self-transition probability used when
#'   constructing a transition matrix automatically. A scalar is recycled by
#'   state.
#' @param null_state Either `"pointmass"` (default), which fixes the null prior
#'   to the Dirac mass at zero, or `"adaptive"`, which learns an ash mixture
#'   centered at zero.
#' @param fixed_pointmass_states Optional integer indices of states whose ash
#'   mixture is fixed entirely on its point component. `NULL` derives the value
#'   from `null_state`; an explicit value is an advanced override.
#' @param shared_mixture Logical; if `TRUE`, all free states share one learned
#'   vector of ash mixture weights. The default learns each state separately.
#' @param estimate_init Logical; estimate initial state probabilities rather
#'   than keeping their initialized values fixed.
#' @param transition_prior Dirichlet parameters for allowed transition
#'   probabilities. Supply a scalar or a state-by-state matrix; values must be
#'   at least one.
#' @param mixture_prior Dirichlet parameters for ash mixture weights. Supply a
#'   scalar, a vector with one entry per scale, or a state-by-scale matrix;
#'   values must be at least one.
#' @param init_prior Dirichlet parameters for the initial distribution, used
#'   only when `estimate_init = TRUE`.
#' @param learn_state_means Logical; learn eligible non-null state centers after
#'   the fixed-grid warm-up.
#' @param fixed_mean_states Integer state indices whose centers are never moved.
#'   State 1 should normally remain fixed.
#' @param mean_update_start First EM iteration at which state-center updates and
#'   center-eligibility decisions are allowed.
#' @param mean_min_effective_count Minimum smoothed state occupancy required for
#'   a center update.
#' @param mean_min_pointmass_weight Minimum fitted weight on the state's point
#'   component required for center learning.
#' @param mean_min_self_transition Minimum self-transition probability required
#'   for center learning. This favors persistent plateau-like states.
#' @param mean_damping Number in `(0, 1]` multiplying each eligible center move.
#' @param mean_bounds Either `"voronoi"` (default), which keeps centers in
#'   nonoverlapping anchor cells, or `"none"`.
#' @param prune_states Logical; enable likelihood-checked dynamic state pruning.
#' @param prune_start First EM iteration eligible for dynamic pruning.
#' @param prune_every Positive number of EM iterations between pruning checks.
#' @param prune_min_state_count,prune_min_state_fraction A state is a
#'   low-occupancy pruning candidate below the larger absolute/fractional cutoff.
#' @param prune_max_fraction Maximum fraction of current states considered in
#'   one deletion batch. It must lie strictly between zero and one.
#' @param prune_max_loglik_loss Maximum allowed decrease in full HMM marginal
#'   log likelihood for an accepted deletion batch.
#' @param merge_distance Distance below which adjacent fitted centers are
#'   considered near duplicates. `NULL` uses `0.05 * median(se)`.
#' @param null_model_selection Either `"bic"` for the strict all-zero safety
#'   comparison or `"none"` to retain the fitted HMM unconditionally.
#' @param null_selection_gamma Nonnegative multiplier for the additional
#'   candidate-grid term in the strict-null information criterion.
#' @param null_bic_margin Nonnegative margin favoring retention of the fitted
#'   HMM over the strict null.
#' @param parameter_tolerance Positive threshold used when counting numerically
#'   active parameters for the strict-null information criterion.
#' @param step_penalty Nonnegative penalty for each state change in the separate
#'   penalized decoder. `NULL` uses `step_penalty_scale * log(length(y))`.
#' @param step_penalty_scale Nonnegative multiplier for the default step penalty.
#'   The default `1.5` accounts for both an added segment level and an unknown
#'   change location, reducing transient one-observation steps.
#' @param maxiter Maximum number of generalized EM iterations.
#' @param tolerance Positive relative convergence tolerance for the penalized
#'   fixed-dimensional objective.
#' @param verbose Logical; print grid selection and per-iteration diagnostics.
#'
#' @return An object of class `ash_hmm_fit`, which is a list containing:
#'   \describe{
#'   \item{call}{The matched function call.}
#'   \item{state_probability}{An observation-by-state matrix of smoothed
#'     probabilities `Pr(Q[t] = m | y)`.}
#'   \item{posterior}{A list containing the marginal posterior `mean`, `sd`,
#'     `probability_ge_zero`, `probability_le_zero`, `probability_zero`, and
#'     local false sign rate `lfsr` for every observation.}
#'   \item{viterbi_state}{The ordinary joint MAP state path under the fitted
#'     transition matrix.}
#'   \item{penalized_state}{The state path from the explicit change-penalized
#'     decoder.}
#'   \item{step_selection}{A list with a segment table, per-sequence counts,
#'     total `step_count`, total `change_count`, `occupied_state_count`, decoded
#'     `state`, and the applied `penalty`.}
#'   \item{boundary_probability}{Posterior probabilities of a state change
#'     between adjacent observations; entries spanning independent sequences
#'     are `NA`.}
#'   \item{fitted}{Fitted centers and scale grid, mixture weights, transition
#'     matrix and mask, initial probabilities, state identifiers and
#'     occupancies, null-state type, effect support, constraints, and
#'     mean-learning settings.}
#'   \item{grid}{Automatic-grid settings, original and retained candidates,
#'     screening statistics, selected scale grid, and pruning history.}
#'   \item{log_likelihood,log_null,log_evidence_ratio}{The fitted HMM marginal
#'     log likelihood, exact all-zero log likelihood, and their difference.}
#'   \item{model_selection}{Details of the optional strict-null comparison,
#'     including criteria, effective dimension, and whether the result was
#'     collapsed to the exact null.}
#'   \item{history}{Per-iteration likelihood, objective, number of states,
#'     moved centers, and pruned-state count.}
#'   \item{mean_history}{Long-form history of state centers and occupancies.}
#'   \item{pruning_history}{One row per removed state, recording its persistent
#'     identifier, original grid index, centers, occupancy, and reason.}
#'   \item{converged}{Logical convergence indicator.}
#'   \item{iterations}{Number of completed generalized EM/model-reduction
#'     iterations recorded after initialization.}
#'   }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' truth <- c(rep(0, 50), rep(2, 50), rep(0, 50))
#' se <- rep(0.5, length(truth))
#' y <- rnorm(length(truth), truth, se)
#'
#' fit <- fit_ash_hmm(y, se, nonnegative_state_means = TRUE)
#' fit$posterior$mean
#' fit$step_selection$segments
#'
#' # Use a signed grid when negative state centers are scientifically possible.
#' signed_fit <- fit_ash_hmm(y - 1, se, nonnegative_state_means = FALSE)
#' }
#'
#' @export
fit_ash_hmm <- function(y, se, mu = NULL, prior_sd = NULL,
                        half_grid = 40L,
                        grid_shape = 1.5,
                        grid_expansion = 1.5,
                        grid_max_abs = NULL,
                        nonnegative_state_means = FALSE,
                        positive_state_means = NULL,
                        effect_support = c("auto", "nonnegative", "real"),
                        positive_mean_floor = NULL,
                        prefilter = NULL,
                        min_state_count = 0.5,
                        min_state_fraction = 1e-5,
                        screening_prior_sd = 0,
                        screening_block_size = 5000L,
                        sequence_id = rep(1L, length(y)),
                        forward_backward_engine = c("auto", "scaled", "log"),
                        topology = c("full", "hub"),
                        transition_mask = NULL,
                        init_transition = NULL,
                        init_prob = NULL,
                        init_rho = NULL,
                        stay_probability = 0.95,
                        null_state = c("pointmass", "adaptive"),
                        fixed_pointmass_states = NULL,
                        shared_mixture = FALSE,
                        estimate_init = FALSE,
                        transition_prior = 1,
                        mixture_prior = 1,
                        init_prior = 1,
                        learn_state_means = TRUE,
                        fixed_mean_states = 1L,
                        mean_update_start = 3L,
                        mean_min_effective_count = 2,
                        mean_min_pointmass_weight = 0,
                        mean_min_self_transition = 0.93,
                        mean_damping = 1,
                        mean_bounds = c("voronoi", "none"),
                        prune_states = TRUE,
                        prune_start = 5L,
                        prune_every = 5L,
                        prune_min_state_count = 1,
                        prune_min_state_fraction = 1e-4,
                        prune_max_fraction = 0.25,
                        prune_max_loglik_loss = 0.05,
                        merge_distance = NULL,
                        null_model_selection = c("bic", "none"),
                        null_selection_gamma = 1,
                        null_bic_margin = 0,
                        parameter_tolerance = 1e-6,
                        step_penalty = NULL,
                        step_penalty_scale = 1.5,
                        maxiter = 10L,
                        tolerance = 1e-5,
                        verbose = TRUE) {
  call <- match.call()
  nonnegative_was_missing <- missing(nonnegative_state_means)
  forward_backward_engine <- match.arg(forward_backward_engine)
  topology <- match.arg(topology)
  mean_bounds <- match.arg(mean_bounds)
  null_model_selection <- match.arg(null_model_selection)
  null_state <- match.arg(null_state)
  effect_support <- match.arg(effect_support)
  .ash_hmm_check_data(y, se, sequence_id)
  if (!is.null(positive_state_means)) {
    if (!is.logical(positive_state_means) ||
        length(positive_state_means) != 1L || is.na(positive_state_means)) {
      stop("'positive_state_means' must be NULL, TRUE, or FALSE.",
           call. = FALSE)
    }
    if (!nonnegative_was_missing &&
        !identical(nonnegative_state_means, positive_state_means)) {
      stop(paste("'nonnegative_state_means' and its compatibility alias",
                 "'positive_state_means' must agree when both are supplied."),
           call. = FALSE)
    }
    nonnegative_state_means <- positive_state_means
  }
  if (!is.logical(nonnegative_state_means) ||
      length(nonnegative_state_means) != 1L || is.na(nonnegative_state_means)) {
    stop("'nonnegative_state_means' must be TRUE or FALSE.", call. = FALSE)
  }
  expected_support <- if (nonnegative_state_means) "nonnegative" else "real"
  if (effect_support == "auto") {
    effect_support <- expected_support
  } else if (effect_support != expected_support) {
    stop(paste0("'effect_support = \"", effect_support,
                "\"' conflicts with 'nonnegative_state_means = ",
                nonnegative_state_means, "'."), call. = FALSE)
  }

  automatic_mu <- is.null(mu)
  if (is.null(prefilter)) prefilter <- automatic_mu
  if (!is.logical(prefilter) || length(prefilter) != 1L || is.na(prefilter)) {
    stop("'prefilter' must be TRUE or FALSE.", call. = FALSE)
  }
  if (automatic_mu && (!is.null(transition_mask) || !is.null(init_transition) ||
                       !is.null(init_prob) || !is.null(init_rho))) {
    stop(paste("Supply 'mu' explicitly when using a custom transition mask,",
               "initial transition matrix, initial probabilities, or ash weights."),
         call. = FALSE)
  }

  if (automatic_mu) {
    if (prefilter) {
      grid_selection <- ash_hmm_select_grid(
        y, se,
        half_grid = half_grid,
        shape = grid_shape,
        expansion = grid_expansion,
        max_abs = grid_max_abs,
        min_state_count = min_state_count,
        min_state_fraction = min_state_fraction,
        screening_prior_sd = screening_prior_sd,
        block_size = screening_block_size,
        positive_only = nonnegative_state_means)
      mu <- grid_selection$selected_mu
    } else {
      mu <- ash_mu_grid(y, half_grid = half_grid, shape = grid_shape,
                        expansion = grid_expansion, max_abs = grid_max_abs,
                        positive_only = nonnegative_state_means)
      grid_selection <- list(
        automatic = TRUE,
        full_mu = mu,
        selected_mu = mu,
        selected_index = seq_along(mu),
        removed_index = integer(),
        effective_count = rep(NA_real_, length(mu)),
        mean_weight = rep(NA_real_, length(mu)),
        cutoff = NA_real_,
        settings = list(prefilter = FALSE))
    }
  } else {
    if (!is.numeric(mu) || !length(mu) || anyNA(mu) || any(!is.finite(mu))) {
      stop("'mu' must be a finite numeric vector.", call. = FALSE)
    }
    if (prefilter) {
      if (!is.null(transition_mask) || !is.null(init_transition) ||
          !is.null(init_prob) || !is.null(init_rho)) {
        stop("Prefilter a supplied grid before providing state-sized initial values.",
             call. = FALSE)
      }
      grid_selection <- .ash_hmm_prefilter_mu(
        y, se, mu,
        min_state_count = min_state_count,
        min_state_fraction = min_state_fraction,
        screening_prior_sd = screening_prior_sd,
        block_size = screening_block_size)
      grid_selection$automatic <- FALSE
      mu <- grid_selection$selected_mu
    } else {
      grid_selection <- list(
        automatic = FALSE,
        full_mu = mu,
        selected_mu = mu,
        selected_index = seq_along(mu),
        removed_index = integer(),
        effective_count = rep(NA_real_, length(mu)),
        mean_weight = rep(NA_real_, length(mu)),
        cutoff = NA_real_,
        settings = list(prefilter = FALSE))
    }
  }
  grid_selection$settings$nonnegative_state_means <-
    nonnegative_state_means
  grid_selection$settings$effect_support <- effect_support

  if (nonnegative_state_means) {
    if (mu[1L] != 0 || (length(mu) > 1L && any(mu[-1L] <= 0))) {
      stop(paste("With 'nonnegative_state_means = TRUE', mu[1] must equal zero",
                 "and every non-null state center must be strictly positive."),
           call. = FALSE)
    }
  } else if (mu[1L] != 0) {
    warning("The hub state is state 1, but mu[1] is not zero.")
  }
  automatic_prior_sd <- is.null(prior_sd)
  if (automatic_prior_sd) prior_sd <- ash_sigma_grid(y, se)
  if (!is.numeric(prior_sd) || !length(prior_sd) || anyNA(prior_sd) ||
      any(!is.finite(prior_sd)) || any(prior_sd < 0) || prior_sd[1L] != 0) {
    stop("'prior_sd' must be finite, nonnegative, and start with zero.", call. = FALSE)
  }
  if (is.unsorted(prior_sd, strictly = FALSE)) {
    stop("'prior_sd' must be sorted in nondecreasing order.", call. = FALSE)
  }

  n <- length(y)
  m <- length(mu)
  l <- length(prior_sd)
  maxiter <- as.integer(maxiter)
  if (maxiter < 1L || tolerance <= 0) {
    stop("Require maxiter >= 1 and tolerance > 0.", call. = FALSE)
  }
  if (!is.logical(learn_state_means) || length(learn_state_means) != 1L ||
      is.na(learn_state_means) || !is.logical(prune_states) ||
      length(prune_states) != 1L || is.na(prune_states)) {
    stop("'learn_state_means' and 'prune_states' must be TRUE or FALSE.",
         call. = FALSE)
  }
  mean_update_start <- as.integer(mean_update_start)
  prune_start <- as.integer(prune_start)
  prune_every <- as.integer(prune_every)
  if (mean_update_start < 1L || prune_start < 1L || prune_every < 1L) {
    stop("Mean/pruning iteration controls must be positive integers.",
         call. = FALSE)
  }
  scalar_nonnegative <- function(x) {
    is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x) && x >= 0
  }
  if (is.null(positive_mean_floor)) {
    positive_mean_floor <- sqrt(.Machine$double.eps) *
      max(1, stats::median(se))
  }
  if (!scalar_nonnegative(mean_min_effective_count) ||
      !scalar_nonnegative(mean_min_pointmass_weight) ||
      mean_min_pointmass_weight > 1 ||
      !scalar_nonnegative(mean_min_self_transition) ||
      mean_min_self_transition > 1 ||
      !is.numeric(mean_damping) || length(mean_damping) != 1L ||
      !is.finite(mean_damping) || mean_damping <= 0 || mean_damping > 1 ||
      !scalar_nonnegative(prune_min_state_count) ||
      !scalar_nonnegative(prune_min_state_fraction) ||
      prune_min_state_fraction > 1 ||
      !is.numeric(prune_max_fraction) || length(prune_max_fraction) != 1L ||
      !is.finite(prune_max_fraction) || prune_max_fraction <= 0 ||
      prune_max_fraction >= 1 ||
      !scalar_nonnegative(prune_max_loglik_loss) ||
      !scalar_nonnegative(positive_mean_floor) ||
      !scalar_nonnegative(null_selection_gamma) ||
      !scalar_nonnegative(null_bic_margin) ||
      !scalar_nonnegative(parameter_tolerance) || parameter_tolerance == 0 ||
      !scalar_nonnegative(step_penalty_scale)) {
    stop("Invalid state-mean learning or pruning control.", call. = FALSE)
  }
  if (nonnegative_state_means && positive_mean_floor <= 0) {
    stop("'positive_mean_floor' must be strictly positive.", call. = FALSE)
  }
  if (is.null(step_penalty)) {
    step_penalty <- step_penalty_scale * log(max(2, n))
  }
  if (!scalar_nonnegative(step_penalty)) {
    stop("'step_penalty' must be NULL or a finite nonnegative scalar.",
         call. = FALSE)
  }
  if (is.null(merge_distance)) merge_distance <- 0.05 * stats::median(se)
  if (!scalar_nonnegative(merge_distance)) {
    stop("'merge_distance' must be NULL or a finite nonnegative scalar.",
         call. = FALSE)
  }
  if (verbose && automatic_mu) {
    message(sprintf("automatic mean grid: retained %d of %d states before HMM fitting",
                    length(mu), length(grid_selection$full_mu)))
  }
  if (is.null(fixed_pointmass_states)) {
    fixed_pointmass_states <- if (null_state == "pointmass") 1L else integer()
  }
  fixed_pointmass_states <- sort(unique(as.integer(fixed_pointmass_states)))
  if (any(fixed_pointmass_states < 1L | fixed_pointmass_states > m)) {
    stop("'fixed_pointmass_states' contains an invalid state index.", call. = FALSE)
  }
  fixed_mean_states <- sort(unique(as.integer(fixed_mean_states)))
  if (any(fixed_mean_states < 1L | fixed_mean_states > m)) {
    stop("'fixed_mean_states' contains an invalid state index.", call. = FALSE)
  }
  if (learn_state_means && mean_bounds == "voronoi" && anyDuplicated(mu)) {
    stop("Initial state means must be distinct with Voronoi mean bounds.",
         call. = FALSE)
  }
  mean_anchor <- mu
  state_id <- seq_len(m)
  grid_index <- grid_selection$selected_index
  pruned_grid_index <- integer()
  pruning_history <- data.frame(
    iteration = integer(), state_id = integer(), grid_index = integer(),
    anchor_mu = numeric(), fitted_mu = numeric(), occupancy = numeric(),
    reason = character(), stringsAsFactors = FALSE)

  if (is.null(transition_mask)) {
    transition_mask <- ash_hmm_transition_mask(m, topology)
  } else {
    if (!is.matrix(transition_mask) || any(dim(transition_mask) != c(m, m))) {
      stop("'transition_mask' must be an M by M matrix.", call. = FALSE)
    }
    transition_mask <- transition_mask != 0
    if (any(rowSums(transition_mask) == 0)) {
      stop("Every state needs at least one allowed outgoing transition.", call. = FALSE)
    }
  }

  if (is.null(init_transition)) {
    A <- .ash_hmm_initial_transition(transition_mask, stay_probability)
  } else {
    A <- .ash_hmm_validate_transition(init_transition, transition_mask)
  }

  if (is.null(init_prob)) init_prob <- .ash_hmm_stationary(A)
  if (!is.numeric(init_prob) || length(init_prob) != m || anyNA(init_prob) ||
      any(!is.finite(init_prob)) || any(init_prob < 0) || sum(init_prob) <= 0) {
    stop("'init_prob' must be a nonnegative length-M vector with positive mass.", call. = FALSE)
  }
  init_prob <- .ash_hmm_normalize(init_prob)

  if (is.null(init_rho)) {
    rho <- matrix(1 / l, m, l)
  } else {
    if (!is.matrix(init_rho) || any(dim(init_rho) != c(m, l)) ||
        anyNA(init_rho) || any(!is.finite(init_rho)) || any(init_rho < 0) ||
        any(rowSums(init_rho) <= 0)) {
      stop("'init_rho' must be a nonnegative M by L matrix with positive row sums.",
           call. = FALSE)
    }
    rho <- init_rho / rowSums(init_rho)
  }
  if (length(fixed_pointmass_states)) {
    rho[fixed_pointmass_states, ] <- 0
    rho[fixed_pointmass_states, 1L] <- 1
  }

  transition_prior <- .ash_hmm_expand_prior(transition_prior, m, m,
                                            "transition_prior")
  mixture_prior <- .ash_hmm_expand_prior(mixture_prior, m, l, "mixture_prior")
  init_prior <- rep(init_prior, length.out = m)
  if (anyNA(init_prior) || any(!is.finite(init_prior)) || any(init_prior < 1)) {
    stop("Every 'init_prior' entry must be finite and at least one.", call. = FALSE)
  }

  mixture_active <- matrix(TRUE, m, l)
  if (length(fixed_pointmass_states)) {
    mixture_active[fixed_pointmass_states, ] <- FALSE
  }

  objective_value <- function(log_likelihood, A, rho, init_prob,
                              transition_prior_current = transition_prior,
                              transition_mask_current = transition_mask,
                              mixture_prior_current = mixture_prior,
                              mixture_active_current = mixture_active,
                              init_prior_current = init_prior) {
    ans <- log_likelihood +
      .ash_hmm_dirichlet_penalty(
        A, transition_prior_current, transition_mask_current) +
      .ash_hmm_dirichlet_penalty(
        rho, mixture_prior_current, mixture_active_current)
    if (estimate_init) {
      use <- init_prior_current > 1
      if (any(use)) {
        if (any(init_prob[use] <= 0)) return(-Inf)
        ans <- ans + sum((init_prior_current[use] - 1) * log(init_prob[use]))
      }
    }
    ans
  }

  log_emission <- .ash_hmm_log_emission(
    y, se, mu, prior_sd, rho, effect_support)
  fb <- .ash_hmm_forward_backward(log_emission, A, init_prob,
                                  transition_mask, sequence_id,
                                  engine = forward_backward_engine)
  objective <- objective_value(fb$log_likelihood, A, rho, init_prob)
  history <- data.frame(iteration = 0L,
                        log_likelihood = fb$log_likelihood,
                        objective = objective,
                        states = m,
                        means_updated = 0L,
                        states_pruned = 0L)
  mean_history <- list(data.frame(
    iteration = 0L, state_id = state_id, grid_index = grid_index,
    anchor_mu = mean_anchor, mu = mu,
    occupancy = colSums(fb$state_probability)))
  mean_state_enabled <- rep(FALSE, m)
  converged <- FALSE

  for (iter in seq_len(maxiter)) {
    free_states <- setdiff(seq_len(m), fixed_pointmass_states)
    statistics <- .ash_hmm_component_statistics(
      y, se, mu, prior_sd, rho, log_emission, fb$state_probability,
      effect_support)
    component_counts <- statistics$counts

    rho_new <- rho
    if (length(free_states)) {
      if (shared_mixture) {
        numerator <- colSums(component_counts[free_states, , drop = FALSE]) +
          colSums(mixture_prior[free_states, , drop = FALSE] - 1)
        common <- .ash_hmm_normalize(
          numerator, fallback = colMeans(rho[free_states, , drop = FALSE]))
        rho_new[free_states, ] <- matrix(
          rep(common, each = length(free_states)), length(free_states), l)
      } else {
        for (state in free_states) {
          numerator <- component_counts[state, ] + mixture_prior[state, ] - 1
          rho_new[state, ] <- .ash_hmm_normalize(
            numerator, fallback = rho[state, ])
        }
      }
    }
    if (length(fixed_pointmass_states)) {
      rho_new[fixed_pointmass_states, ] <- 0
      rho_new[fixed_pointmass_states, 1L] <- 1
    }
    A_new <- matrix(0, m, m)
    for (state in seq_len(m)) {
      allowed <- which(transition_mask[state, ])
      numerator <- fb$transition_counts[state, allowed] +
        transition_prior[state, allowed] - 1
      A_new[state, allowed] <- .ash_hmm_normalize(
        numerator, fallback = A[state, allowed])
    }

    # Long, persistent states are the ones for which estimating a precise
    # plateau center is useful. The self-transition gate avoids moving centers
    # that merely tile a smooth bump or are supported by scattered positions.
    current_mean_gate <- rho_new[, 1L] >= mean_min_pointmass_weight &
      diag(A_new) >= mean_min_self_transition
    if (iter == mean_update_start) {
      mean_state_enabled <- current_mean_gate
    } else if (iter > mean_update_start) {
      # Eligibility is established after the warm-up and can switch off, but
      # not on. This stops a weak state from becoming movable merely because
      # later pruning artificially raises its self-transition probability.
      mean_state_enabled <- mean_state_enabled & current_mean_gate
    }
    persistence_damping <- pmin(1, pmax(
      0, (diag(A_new) - mean_min_self_transition) /
        max(1e-8, 1 - mean_min_self_transition)))
    mean_lower_bound <- if (nonnegative_state_means) {
      c(0, rep(positive_mean_floor, max(0L, m - 1L)))
    } else rep(-Inf, m)
    if (effect_support == "nonnegative") {
      mean_update <- .ash_hmm_update_truncated_means(
        y = y, se = se, mu = mu, anchor = mean_anchor,
        prior_sd = prior_sd, rho = rho, log_emission = log_emission,
        gamma = fb$state_probability, statistics = statistics,
        learn = learn_state_means && iter >= mean_update_start,
        fixed = fixed_mean_states,
        minimum_count = mean_min_effective_count,
        eligible = mean_state_enabled,
        damping = mean_damping * persistence_damping,
        minimum = mean_lower_bound,
        bounds = mean_bounds)
    } else {
      mean_update <- .ash_hmm_update_means(
        mu = mu, anchor = mean_anchor, statistics = statistics,
        learn = learn_state_means && iter >= mean_update_start,
        fixed = fixed_mean_states,
        minimum_count = mean_min_effective_count,
        eligible = mean_state_enabled,
        damping = mean_damping * persistence_damping,
        minimum = mean_lower_bound,
        bounds = mean_bounds)
    }
    mu_new <- mean_update$mu

    init_prob_new <- init_prob
    if (estimate_init) {
      numerator <- fb$initial_counts + init_prior - 1
      init_prob_new <- .ash_hmm_normalize(numerator, fallback = init_prob)
    }

    log_emission_new <- .ash_hmm_log_emission(
      y, se, mu_new, prior_sd, rho_new, effect_support)
    fb_new <- .ash_hmm_forward_backward(
      log_emission_new, A_new, init_prob_new, transition_mask, sequence_id,
      engine = forward_backward_engine)
    objective_new <- objective_value(
      fb_new$log_likelihood, A_new, rho_new, init_prob_new)

    if (objective_new < objective - 1e-7 * (1 + abs(objective))) {
      warning("The fixed-dimension EM objective decreased beyond numerical tolerance.")
    }
    delta_em <- objective_new - objective
    states_pruned_iter <- 0L

    # Dynamic deletion is deliberately outside the fixed-dimensional EM step.
    # Candidate states must be empty or nearly duplicate, and deletion is only
    # accepted when the full HMM marginal log likelihood loses at most the
    # explicitly allowed amount. Groups that fail this check are halved.
    pruning_due <- prune_states && iter >= prune_start &&
      ((iter - prune_start) %% prune_every == 0L) && m > 1L
    if (pruning_due) {
      occupancy_new <- colSums(fb_new$state_probability)
      protected <- sort(unique(c(1L, fixed_pointmass_states,
                                 fixed_mean_states)))
      cutoff <- max(prune_min_state_count, n * prune_min_state_fraction)
      low_candidates <- setdiff(which(occupancy_new < cutoff), protected)
      close_candidates <- .ash_hmm_close_state_candidates(
        mu_new, occupancy_new, merge_distance, protected)
      candidates <- unique(c(low_candidates, close_candidates))
      if (length(candidates)) {
        candidates <- candidates[order(occupancy_new[candidates])]
        maximum_drop <- max(1L, floor(prune_max_fraction * m))
        maximum_drop <- min(maximum_drop, length(candidates),
                            m - length(unique(protected)))
        number_to_try <- maximum_drop
        accepted <- NULL

        while (number_to_try >= 1L && is.null(accepted)) {
          drop <- candidates[seq_len(number_to_try)]
          keep <- setdiff(seq_len(m), drop)
          mask_try <- transition_mask[keep, keep, drop = FALSE]
          A_try <- A_new[keep, keep, drop = FALSE]
          valid_rows <- rowSums(mask_try) > 0 & rowSums(A_try) > 0
          if (all(valid_rows)) {
            A_try <- A_try / rowSums(A_try)
            init_try <- .ash_hmm_normalize(init_prob_new[keep])
            rho_try <- rho_new[keep, , drop = FALSE]
            mu_try <- mu_new[keep]
            transition_prior_try <-
              transition_prior[keep, keep, drop = FALSE]
            mixture_prior_try <- mixture_prior[keep, , drop = FALSE]
            init_prior_try <- init_prior[keep]
            fixed_pointmass_try <- match(fixed_pointmass_states, keep,
                                         nomatch = 0L)
            fixed_pointmass_try <- fixed_pointmass_try[fixed_pointmass_try > 0L]
            fixed_mean_try <- match(fixed_mean_states, keep, nomatch = 0L)
            fixed_mean_try <- fixed_mean_try[fixed_mean_try > 0L]
            mixture_active_try <- matrix(TRUE, length(keep), l)
            if (length(fixed_pointmass_try)) {
              mixture_active_try[fixed_pointmass_try, ] <- FALSE
            }
            log_emission_try <- .ash_hmm_log_emission(
              y, se, mu_try, prior_sd, rho_try, effect_support)
            fb_try <- .ash_hmm_forward_backward(
              log_emission_try, A_try, init_try, mask_try, sequence_id,
              engine = forward_backward_engine)
            loss <- fb_new$log_likelihood - fb_try$log_likelihood
            if (is.finite(loss) && loss <= prune_max_loglik_loss) {
              objective_try <- objective_value(
                fb_try$log_likelihood, A_try, rho_try, init_try,
                transition_prior_try, mask_try, mixture_prior_try,
                mixture_active_try, init_prior_try)
              accepted <- list(
                drop = drop, keep = keep, A = A_try, rho = rho_try,
                mu = mu_try, init_prob = init_try, mask = mask_try,
                transition_prior = transition_prior_try,
                mixture_prior = mixture_prior_try,
                init_prior = init_prior_try,
                mixture_active = mixture_active_try,
                fixed_pointmass_states = fixed_pointmass_try,
                fixed_mean_states = fixed_mean_try,
                log_emission = log_emission_try, fb = fb_try,
                objective = objective_try)
            }
          }
          number_to_try <- floor(number_to_try / 2L)
        }

        if (!is.null(accepted)) {
          removed <- accepted$drop
          reason <- ifelse(
            removed %in% low_candidates & removed %in% close_candidates,
            "low occupancy and near duplicate",
            ifelse(removed %in% close_candidates,
                   "near duplicate", "low occupancy"))
          pruning_history <- rbind(
            pruning_history,
            data.frame(iteration = rep(iter, length(removed)),
                       state_id = state_id[removed],
                       grid_index = grid_index[removed],
                       anchor_mu = mean_anchor[removed],
                       fitted_mu = mu_new[removed],
                       occupancy = occupancy_new[removed],
                       reason = reason, stringsAsFactors = FALSE))
          pruned_grid_index <- c(pruned_grid_index, grid_index[removed])
          state_id <- state_id[accepted$keep]
          grid_index <- grid_index[accepted$keep]
          mean_anchor <- mean_anchor[accepted$keep]
          mean_state_enabled <- mean_state_enabled[accepted$keep]
          A_new <- accepted$A
          rho_new <- accepted$rho
          mu_new <- accepted$mu
          init_prob_new <- accepted$init_prob
          transition_mask <- accepted$mask
          transition_prior <- accepted$transition_prior
          mixture_prior <- accepted$mixture_prior
          init_prior <- accepted$init_prior
          mixture_active <- accepted$mixture_active
          fixed_pointmass_states <- accepted$fixed_pointmass_states
          fixed_mean_states <- accepted$fixed_mean_states
          log_emission_new <- accepted$log_emission
          fb_new <- accepted$fb
          objective_new <- accepted$objective
          states_pruned_iter <- length(removed)
          m <- length(mu_new)
        }
      }
    }

    history <- rbind(
      history,
      data.frame(iteration = iter,
                 log_likelihood = fb_new$log_likelihood,
                 objective = objective_new,
                 states = m,
                 means_updated = sum(mean_update$updated),
                 states_pruned = states_pruned_iter))
    mean_history[[length(mean_history) + 1L]] <- data.frame(
      iteration = iter, state_id = state_id, grid_index = grid_index,
      anchor_mu = mean_anchor, mu = mu_new,
      occupancy = colSums(fb_new$state_probability))

    if (verbose) {
      message(sprintf(
        paste0("iteration %d: logLik = %.10f, objective = %.10f, ",
               "states = %d, means moved = %d, pruned = %d"),
        iter, fb_new$log_likelihood, objective_new, m,
        sum(mean_update$updated), states_pruned_iter))
    }

    A <- A_new
    rho <- rho_new
    mu <- mu_new
    init_prob <- init_prob_new
    log_emission <- log_emission_new
    fb <- fb_new
    objective <- objective_new

    minimum_iterations <- max(
      if (learn_state_means) mean_update_start else 1L,
      if (prune_states) prune_start else 1L)
    pruning_checkpoint_complete <- !prune_states || pruning_due
    if (iter >= minimum_iterations && pruning_checkpoint_complete &&
        states_pruned_iter == 0L &&
        abs(delta_em) <= tolerance * (1 + abs(objective))) {
      converged <- TRUE
      break
    }
  }

  log_null <- sum(stats::dnorm(y, mean = 0, sd = se, log = TRUE))
  initial_selected_mu <- grid_selection$selected_mu
  initial_selected_index <- grid_selection$selected_index
  learned_mean_states <- setdiff(which(mean_state_enabled), fixed_mean_states)
  effective_parameters <- .ash_hmm_effective_parameter_count(
    A, transition_mask, rho, fixed_pointmass_states,
    learned_mean_states, estimate_init, init_prob, parameter_tolerance)
  candidate_count <- max(1L, length(grid_selection$full_mu) - 1L)
  full_bic <- -2 * fb$log_likelihood +
    effective_parameters * log(max(2, n)) +
    2 * null_selection_gamma * log(candidate_count)
  null_bic <- -2 * log_null
  collapse_to_null <- null_model_selection == "bic" &&
    null_bic + null_bic_margin <= full_bic
  model_selection <- list(
    method = null_model_selection,
    selected = if (collapse_to_null) {
      "strict_null"
    } else if (effect_support == "nonnegative") {
      "one_sided_ash_hmm"
    } else {
      "signed_ash_hmm"
    },
    collapsed_to_null = collapse_to_null,
    full_log_likelihood = fb$log_likelihood,
    null_log_likelihood = log_null,
    effective_parameters = effective_parameters,
    candidate_count = candidate_count,
    full_bic = full_bic,
    null_bic = null_bic,
    grid_penalty_gamma = null_selection_gamma,
    margin = null_bic_margin)

  if (collapse_to_null) {
    removed <- if (m > 1L) 2L:m else integer()
    occupancy_before_null <- colSums(fb$state_probability)
    pruning_history <- rbind(
      pruning_history,
      data.frame(iteration = rep(max(history$iteration) + 1L, length(removed)),
                 state_id = state_id[removed],
                 grid_index = grid_index[removed],
                 anchor_mu = mean_anchor[removed],
                 fitted_mu = mu[removed],
                 occupancy = occupancy_before_null[removed],
                 reason = rep("global BIC selected strict null", length(removed)),
                 stringsAsFactors = FALSE))
    pruned_grid_index <- c(pruned_grid_index, grid_index[removed])

    state_id <- state_id[1L]
    grid_index <- grid_index[1L]
    mean_anchor <- 0
    mu <- 0
    m <- 1L
    rho <- matrix(0, 1L, l)
    rho[1L, 1L] <- 1
    A <- matrix(1, 1L, 1L)
    transition_mask <- matrix(TRUE, 1L, 1L)
    init_prob <- 1
    transition_prior <- matrix(transition_prior[1L, 1L], 1L, 1L)
    mixture_prior <- matrix(mixture_prior[1L, ], 1L, l)
    init_prior <- init_prior[1L]
    mixture_active <- matrix(FALSE, 1L, l)
    fixed_pointmass_states <- 1L
    fixed_mean_states <- 1L
    mean_state_enabled <- FALSE
    log_emission <- .ash_hmm_log_emission(
      y, se, mu, prior_sd, rho, effect_support)
    fb <- .ash_hmm_forward_backward(
      log_emission, A, init_prob, transition_mask, sequence_id,
      engine = forward_backward_engine)
    objective <- objective_value(
      fb$log_likelihood, A, rho, init_prob,
      transition_prior, transition_mask, mixture_prior,
      mixture_active, init_prior)
    history <- rbind(
      history,
      data.frame(iteration = max(history$iteration) + 1L,
                 log_likelihood = fb$log_likelihood,
                 objective = objective,
                 states = 1L,
                 means_updated = 0L,
                 states_pruned = length(removed)))
    mean_history[[length(mean_history) + 1L]] <- data.frame(
      iteration = max(history$iteration), state_id = state_id,
      grid_index = grid_index, anchor_mu = 0, mu = 0,
      occupancy = n)
    converged <- TRUE
  }

  posterior <- .ash_hmm_posterior_summary(
    y, se, mu, prior_sd, rho, log_emission, fb$state_probability,
    effect_support)
  viterbi_state <- .ash_hmm_viterbi(log_emission, A, init_prob,
                                    transition_mask, sequence_id)
  penalized_state <- .ash_hmm_penalized_viterbi(
    log_emission, transition_mask, sequence_id, step_penalty)
  step_summary <- .ash_hmm_step_summary(
    penalized_state, sequence_id, mu)
  step_summary$state <- penalized_state
  step_summary$penalty <- step_penalty
  grid_selection$initial_selected_mu <- initial_selected_mu
  grid_selection$initial_selected_index <- initial_selected_index
  grid_selection$selected_mu <- mu
  grid_selection$selected_index <- grid_index
  grid_selection$removed_index <- sort(unique(c(
    grid_selection$removed_index, pruned_grid_index)))
  grid_selection$pruning_history <- pruning_history

  ans <- list(call = call,
              state_probability = fb$state_probability,
              posterior = posterior,
              viterbi_state = viterbi_state,
              penalized_state = penalized_state,
              step_selection = step_summary,
              boundary_probability = fb$boundary_probability,
              fitted = list(mu = mu,
                            prior_sd = prior_sd,
                            mixture_weight = rho,
                            transition = A,
                            transition_mask = transition_mask,
                            forward_backward_engine =
                              forward_backward_engine,
                            init_prob = init_prob,
                            shared_mixture = shared_mixture,
                            null_state = if (1L %in% fixed_pointmass_states) {
                              "pointmass"
                            } else "adaptive",
                            fixed_pointmass_states = fixed_pointmass_states,
                            fixed_mean_states = fixed_mean_states,
                            state_id = state_id,
                            mean_anchor = mean_anchor,
                            state_occupancy = colSums(fb$state_probability),
                            effect_support = effect_support,
                            learn_state_means = learn_state_means,
                            mean_bounds = mean_bounds,
                            mean_state_enabled = mean_state_enabled,
                            mean_update_start = mean_update_start,
                            mean_min_effective_count = mean_min_effective_count,
                            mean_min_pointmass_weight = mean_min_pointmass_weight,
                            mean_min_self_transition = mean_min_self_transition,
                            mean_damping = mean_damping,
                            nonnegative_state_means = nonnegative_state_means,
                            positive_state_means = nonnegative_state_means,
                            positive_mean_floor = positive_mean_floor),
              grid = c(grid_selection,
                       list(automatic_prior_sd = automatic_prior_sd,
                            selected_prior_sd = prior_sd)),
              log_likelihood = fb$log_likelihood,
              log_null = log_null,
              log_evidence_ratio = fb$log_likelihood - log_null,
              model_selection = model_selection,
              history = history,
              mean_history = do.call(rbind, mean_history),
              pruning_history = pruning_history,
              converged = converged,
              iterations = nrow(history) - 1L)
  class(ans) <- "ash_hmm_fit"
  ans
}


#' Backward-compatible ash-HMM entry point
#'
#' Calls [fit_ash_hmm()] using the original argument names `x` and `sd`.
#'
#' @param x Numeric vector of noisy effect estimates.
#' @param sd Numeric vector of known standard errors.
#' @param ... Additional arguments passed to [fit_ash_hmm()].
#' @return An object of class `ash_hmm_fit`; see [fit_ash_hmm()].
#' @export
fit_hmm <- function(x, sd, ...) {
  fit_ash_hmm(y = x, se = sd, ...)
}

#' Fit the exact two-state binary Markov model
#'
#' Fits the point-state special case with a symmetric transition matrix
#' `A(q) = matrix(c(1-q, q, q, 1-q), 2, 2, byrow = TRUE)`. The scalar `q` is
#' selected by likelihood maximization over the requested interval.
#'
#' @param y Numeric vector of noisy observations.
#' @param se Numeric vector of known, strictly positive standard errors.
#' @param state_means Two finite point-state means.
#' @param sequence_id Vector identifying independent contiguous sequences.
#' @param q_interval Increasing length-two search interval contained in `[0,1]`.
#' @param grid_size Number of initial likelihood evaluations over `q_interval`.
#' @param init_prob Nonnegative length-two initial state distribution.
#' @return A list containing the estimated flip probability `q`, transition
#'   matrix, log likelihood, smoothed state probabilities, and boundary
#'   probabilities.
#' @export
fit_binary_markov <- function(y, se, state_means = c(0, 1),
                              sequence_id = rep(1L, length(y)),
                              q_interval = c(0, 0.5),
                              grid_size = 501L,
                              init_prob = c(0.5, 0.5)) {
  .ash_hmm_check_data(y, se, sequence_id)
  if (!is.numeric(state_means) || length(state_means) != 2L ||
      anyNA(state_means) || any(!is.finite(state_means))) {
    stop("'state_means' must contain two finite values.", call. = FALSE)
  }
  if (length(q_interval) != 2L || q_interval[1L] < 0 ||
      q_interval[2L] > 1 || q_interval[1L] >= q_interval[2L]) {
    stop("'q_interval' must be an increasing subinterval of [0, 1].",
         call. = FALSE)
  }
  if (!is.numeric(init_prob) || length(init_prob) != 2L || anyNA(init_prob) ||
      any(!is.finite(init_prob)) || any(init_prob < 0) || sum(init_prob) <= 0) {
    stop("'init_prob' must be a nonnegative length-two vector with positive mass.",
         call. = FALSE)
  }
  grid_size <- as.integer(grid_size)
  if (grid_size < 5L) stop("'grid_size' must be at least five.", call. = FALSE)
  init_prob <- .ash_hmm_normalize(init_prob)
  mask <- matrix(TRUE, 2L, 2L)
  log_emission <- cbind(
    stats::dnorm(y, state_means[1L], se, log = TRUE),
    stats::dnorm(y, state_means[2L], se, log = TRUE))

  evaluate <- function(q, details = FALSE) {
    A <- matrix(c(1 - q, q, q, 1 - q), 2L, 2L, byrow = TRUE)
    if (!details) {
      return(.ash_hmm_forward_log_likelihood(
        log_emission, A, init_prob, sequence_id, mask))
    }
    fb <- .ash_hmm_forward_backward(
      log_emission, A, init_prob, mask, sequence_id)
    list(A = A, fb = fb)
  }

  grid <- seq(q_interval[1L], q_interval[2L], length.out = grid_size)
  values <- vapply(grid, evaluate, numeric(1L))
  local <- which(values[2L:(grid_size - 1L)] >= values[1L:(grid_size - 2L)] &
                   values[2L:(grid_size - 1L)] >= values[3L:grid_size]) + 1L
  candidates <- data.frame(q = c(grid[1L], grid[grid_size]),
                           log_likelihood = c(values[1L], values[grid_size]))
  if (length(local)) {
    for (i in local) {
      refined <- stats::optimize(function(q) -evaluate(q),
                                 interval = c(grid[i - 1L], grid[i + 1L]))
      candidates <- rbind(candidates,
                          data.frame(q = refined$minimum,
                                     log_likelihood = -refined$objective))
    }
  }
  best <- candidates[which.max(candidates$log_likelihood), ]
  details <- evaluate(best$q, details = TRUE)
  viterbi <- .ash_hmm_viterbi(log_emission, details$A, init_prob, mask, sequence_id)

  list(q = best$q,
       expected_run_length = if (best$q > 0) 1 / best$q else Inf,
       state_probability = details$fb$state_probability,
       boundary_probability = details$fb$boundary_probability,
       viterbi_state = viterbi,
       transition = details$A,
       init_prob = init_prob,
       log_likelihood = details$fb$log_likelihood,
       search = list(grid = grid, values = values, candidates = candidates))
}
