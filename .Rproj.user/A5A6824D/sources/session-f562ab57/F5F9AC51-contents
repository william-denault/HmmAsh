# ash_hmm.R
#
# A dependency-free implementation of an adaptive-shrinkage hidden Markov
# model for heteroskedastic normal observations.
#
# Model:
#   y_t | beta_t ~ N(beta_t, se_t^2)
#   Q_t is a finite-state Markov chain with transition matrix A
#   beta_t | Q_t=m ~ sum_l rho[m,l] N(mu[m], prior_sd[l]^2)
#
# prior_sd[1] must be zero, so the first component is a point mass at mu[m].
# The means and prior standard deviations may be generated from y and se, but
# are fixed before EM; A, rho, and optionally the initial distribution are fit.


# Fit the general ash-HMM by a coherent EM/Baum-Welch algorithm.
fit_ash_hmm <- function(y, se, mu = NULL, prior_sd = NULL,
                        half_grid = 100L,
                        grid_shape = 3,
                        grid_expansion = 1.5,
                        grid_max_abs = NULL,
                        prefilter = NULL,
                        min_state_count = 0.5,
                        min_state_fraction = 1e-5,
                        screening_prior_sd = 0,
                        screening_block_size = 5000L,
                        sequence_id = rep(1L, length(y)),
                        topology = c("hub", "full"),
                        transition_mask = NULL,
                        init_transition = NULL,
                        init_prob = NULL,
                        init_rho = NULL,
                        stay_probability = 0.95,
                        fixed_pointmass_states = integer(),
                        shared_mixture = FALSE,
                        estimate_init = FALSE,
                        transition_prior = 1,
                        mixture_prior = 1,
                        init_prior = 1,
                        maxiter = 200L,
                        tolerance = 1e-8,
                        verbose = FALSE) {
  call <- match.call()
  topology <- match.arg(topology)
  .ash_hmm_check_data(y, se, sequence_id)

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
        block_size = screening_block_size)
      mu <- grid_selection$selected_mu
    } else {
      mu <- ash_mu_grid(y, half_grid = half_grid, shape = grid_shape,
                        expansion = grid_expansion, max_abs = grid_max_abs)
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

  if (mu[1L] != 0) {
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
  if (verbose && automatic_mu) {
    message(sprintf("automatic mean grid: retained %d of %d states before HMM fitting",
                    length(mu), length(grid_selection$full_mu)))
  }
  fixed_pointmass_states <- sort(unique(as.integer(fixed_pointmass_states)))
  if (any(fixed_pointmass_states < 1L | fixed_pointmass_states > m)) {
    stop("'fixed_pointmass_states' contains an invalid state index.", call. = FALSE)
  }

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

  objective_value <- function(log_likelihood, A, rho, init_prob) {
    ans <- log_likelihood +
      .ash_hmm_dirichlet_penalty(A, transition_prior, transition_mask) +
      .ash_hmm_dirichlet_penalty(rho, mixture_prior, mixture_active)
    if (estimate_init) {
      use <- init_prior > 1
      if (any(use)) {
        if (any(init_prob[use] <= 0)) return(-Inf)
        ans <- ans + sum((init_prior[use] - 1) * log(init_prob[use]))
      }
    }
    ans
  }

  log_emission <- .ash_hmm_log_emission(y, se, mu, prior_sd, rho)
  fb <- .ash_hmm_forward_backward(log_emission, A, init_prob,
                                  transition_mask, sequence_id)
  objective <- objective_value(fb$log_likelihood, A, rho, init_prob)
  history <- data.frame(iteration = 0L,
                        log_likelihood = fb$log_likelihood,
                        objective = objective)
  converged <- FALSE

  free_states <- setdiff(seq_len(m), fixed_pointmass_states)
  for (iter in seq_len(maxiter)) {
    component_counts <- .ash_hmm_component_counts(
      y, se, mu, prior_sd, rho, log_emission, fb$state_probability)

    rho_new <- rho
    if (length(free_states)) {
      if (shared_mixture) {
        numerator <- colSums(component_counts[free_states, , drop = FALSE]) +
          colSums(mixture_prior[free_states, , drop = FALSE] - 1)
        common <- .ash_hmm_normalize(numerator,
                                     fallback = colMeans(rho[free_states, , drop = FALSE]))
        rho_new[free_states, ] <- matrix(rep(common, each = length(free_states)),
                                         length(free_states), l)
      } else {
        for (state in free_states) {
          numerator <- component_counts[state, ] + mixture_prior[state, ] - 1
          rho_new[state, ] <- .ash_hmm_normalize(numerator, fallback = rho[state, ])
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

    init_prob_new <- init_prob
    if (estimate_init) {
      numerator <- fb$initial_counts + init_prior - 1
      init_prob_new <- .ash_hmm_normalize(numerator, fallback = init_prob)
    }

    log_emission_new <- .ash_hmm_log_emission(y, se, mu, prior_sd, rho_new)
    fb_new <- .ash_hmm_forward_backward(log_emission_new, A_new, init_prob_new,
                                        transition_mask, sequence_id)
    objective_new <- objective_value(fb_new$log_likelihood, A_new, rho_new,
                                     init_prob_new)
    history <- rbind(history,
                     data.frame(iteration = iter,
                                log_likelihood = fb_new$log_likelihood,
                                objective = objective_new))

    if (verbose) {
      message(sprintf("iteration %d: logLik = %.10f, objective = %.10f",
                      iter, fb_new$log_likelihood, objective_new))
    }
    if (objective_new < objective - 1e-7 * (1 + abs(objective))) {
      warning("The EM objective decreased beyond numerical tolerance.")
    }

    delta <- objective_new - objective
    A <- A_new
    rho <- rho_new
    init_prob <- init_prob_new
    log_emission <- log_emission_new
    fb <- fb_new
    objective <- objective_new

    if (abs(delta) <= tolerance * (1 + abs(objective))) {
      converged <- TRUE
      break
    }
  }

  posterior <- .ash_hmm_posterior_summary(
    y, se, mu, prior_sd, rho, log_emission, fb$state_probability)
  viterbi_state <- .ash_hmm_viterbi(log_emission, A, init_prob,
                                    transition_mask, sequence_id)
  log_null <- sum(stats::dnorm(y, mean = 0, sd = se, log = TRUE))

  ans <- list(call = call,
              state_probability = fb$state_probability,
              posterior = posterior,
              viterbi_state = viterbi_state,
              boundary_probability = fb$boundary_probability,
              fitted = list(mu = mu,
                            prior_sd = prior_sd,
                            mixture_weight = rho,
                            transition = A,
                            transition_mask = transition_mask,
                            init_prob = init_prob,
                            shared_mixture = shared_mixture,
                            fixed_pointmass_states = fixed_pointmass_states),
              grid = c(grid_selection,
                       list(automatic_prior_sd = automatic_prior_sd,
                            selected_prior_sd = prior_sd)),
              log_likelihood = fb$log_likelihood,
              log_null = log_null,
              log_evidence_ratio = fb$log_likelihood - log_null,
              history = history,
              converged = converged,
              iterations = nrow(history) - 1L)
  class(ans) <- "ash_hmm_fit"
  ans
}

print.ash_hmm_fit <- function(x, ...) {
  cat("Adaptive-shrinkage HMM fit\n")
  cat("  observations:", nrow(x$state_probability), "\n")
  cat("  states:", ncol(x$state_probability))
  if (isTRUE(x$grid$automatic)) {
    cat(" (retained from", length(x$grid$full_mu), "automatic candidates)")
  }
  cat("\n")
  cat("  iterations:", x$iterations,
      if (x$converged) "(converged)" else "(maximum reached)", "\n")
  cat("  log likelihood:", format(x$log_likelihood, digits = 8), "\n")
  invisible(x)
}

# Backward-compatible two-input entry point matching the original function.
fit_hmm <- function(x, sd, ...) {
  fit_ash_hmm(y = x, se = sd, ...)
}

# Exact two-state reduction used in the supplied binary-Markov paper.
# The state means are point masses and A(q) has one symmetric flip parameter.
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
    fb <- .ash_hmm_forward_backward(log_emission, A, init_prob, mask, sequence_id)
    if (details) return(list(A = A, fb = fb))
    fb$log_likelihood
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
