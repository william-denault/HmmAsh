.ash_hmm_version <- "2.0.0-one-sided-ash"

#' Return the implementation version.
#' @export
ash_hmm_version <- function() .ash_hmm_version

.ash_hmm_logsumexp <- function(x) {
  if (!length(x)) return(-Inf)
  z <- max(x)
  if (!is.finite(z)) return(z)
  z + log(sum(exp(x - z)))
}

.ash_hmm_row_logsumexp <- function(x) {
  if (is.null(dim(x))) return(.ash_hmm_logsumexp(x))
  z <- apply(x, 1L, max)
  ans <- rep(-Inf, nrow(x))
  ok <- is.finite(z)
  if (any(ok)) {
    ans[ok] <- z[ok] + log(rowSums(exp(x[ok, , drop = FALSE] - z[ok])))
  }
  ans
}

.ash_hmm_safe_log <- function(x) {
  ans <- rep(-Inf, length(x))
  ans[x > 0] <- log(x[x > 0])
  dim(ans) <- dim(x)
  ans
}

.ash_hmm_normalize <- function(x, fallback = NULL) {
  s <- sum(x)
  if (is.finite(s) && s > 0) return(x / s)
  if (is.null(fallback)) stop("Cannot normalize a zero-mass vector.", call. = FALSE)
  fallback / sum(fallback)
}

.ash_hmm_check_data <- function(y, se, sequence_id) {
  if (!is.numeric(y) || !is.numeric(se) || length(y) != length(se)) {
    stop("'y' and 'se' must be numeric vectors of equal length.", call. = FALSE)
  }
  if (!length(y)) stop("At least one observation is required.", call. = FALSE)
  if (anyNA(y) || anyNA(se) || any(!is.finite(y)) || any(!is.finite(se))) {
    stop("'y' and 'se' must be finite and non-missing.", call. = FALSE)
  }
  if (any(se <= 0)) stop("Every standard error must be strictly positive.", call. = FALSE)
  if (length(sequence_id) != length(y) || anyNA(sequence_id)) {
    stop("'sequence_id' must contain one non-missing value per observation.", call. = FALSE)
  }
  invisible(TRUE)
}

# Stephens-style geometric grid, including the point-mass component at zero.
ash_sigma_grid <- function(y, se, multiplier = sqrt(2), sigma_min = NULL,
                           sigma_max = NULL) {
  .ash_hmm_check_data(y, se, rep(1L, length(y)))
  if (!is.numeric(multiplier) || length(multiplier) != 1L || multiplier <= 1) {
    stop("'multiplier' must be a scalar greater than one.", call. = FALSE)
  }
  if (is.null(sigma_min)) sigma_min <- min(se) / 10
  if (is.null(sigma_max)) {
    signal_var <- max(y^2 - se^2)
    sigma_max <- if (signal_var > 0) 2 * sqrt(signal_var) else 8 * sigma_min
  }
  if (sigma_min <= 0 || sigma_max < sigma_min) {
    stop("Require 0 < sigma_min <= sigma_max.", call. = FALSE)
  }
  positive <- sigma_min * multiplier^(0:ceiling(log(sigma_max / sigma_min) /
                                                  log(multiplier)))
  positive[length(positive)] <- sigma_max
  sort(unique(c(0, positive)))
}

# Symmetric state-mean grid with the zero/hub state first.
ash_mu_grid <- function(y, half_grid = 100L, shape = 3, expansion = 1.5,
                        max_abs = NULL, positive_only = FALSE) {
  if (!is.numeric(y) || !length(y) || anyNA(y) || any(!is.finite(y))) {
    stop("'y' must be a finite numeric vector.", call. = FALSE)
  }
  half_grid <- as.integer(half_grid)
  if (half_grid < 1L || shape <= 0 || expansion <= 0) {
    stop("Require half_grid >= 1, shape > 0, and expansion > 0.", call. = FALSE)
  }
  if (is.null(max_abs)) max_abs <- max(abs(y))
  if (!is.finite(max_abs) || max_abs <= 0) max_abs <- 1
  positive <- seq(0, 1, length.out = half_grid)^(1 / shape) * expansion * max_abs
  if (isTRUE(positive_only)) return(positive)
  c(positive, -positive[-1L])
}

# Screen a candidate state-mean grid before any HMM forward-backward pass.
# Each observation is assigned soft iid weights under point emissions centered
# on the candidate means. Their column sums are effective observation counts.
.ash_hmm_prefilter_mu <- function(y, se, mu,
                                  min_state_count = 0.5,
                                  min_state_fraction = 1e-5,
                                  screening_prior_sd = 0,
                                  block_size = 5000L) {
  if (!is.numeric(mu) || !length(mu) || anyNA(mu) || any(!is.finite(mu))) {
    stop("'mu' must be a finite numeric vector.", call. = FALSE)
  }
  if (length(min_state_count) != 1L || !is.finite(min_state_count) ||
      min_state_count < 0 || length(min_state_fraction) != 1L ||
      !is.finite(min_state_fraction) || min_state_fraction < 0 ||
      min_state_fraction > 1) {
    stop("Grid-screening count and fraction thresholds must be valid scalars.",
         call. = FALSE)
  }
  if (length(screening_prior_sd) != 1L || !is.finite(screening_prior_sd) ||
      screening_prior_sd < 0) {
    stop("'screening_prior_sd' must be a finite nonnegative scalar.",
         call. = FALSE)
  }
  block_size <- as.integer(block_size)
  if (block_size < 1L) stop("'block_size' must be positive.", call. = FALSE)

  n <- length(y)
  m <- length(mu)
  effective_count <- numeric(m)
  starts <- seq.int(1L, n, by = block_size)

  for (first in starts) {
    last <- min(n, first + block_size - 1L)
    index <- first:last
    log_score <- matrix(-Inf, length(index), m)
    for (state in seq_len(m)) {
      log_score[, state] <- stats::dnorm(
        y[index], mean = mu[state],
        sd = sqrt(se[index]^2 + screening_prior_sd^2), log = TRUE)
    }
    log_normalizer <- .ash_hmm_row_logsumexp(log_score)
    effective_count <- effective_count +
      colSums(exp(log_score - log_normalizer))
  }

  cutoff <- max(min_state_count, n * min_state_fraction)
  keep <- which(effective_count >= cutoff)
  # State 1 is the zero/hub state for automatically constructed grids and is
  # retained even when the data are entirely nonzero.
  keep <- sort(unique(c(1L, keep)))

  list(full_mu = mu,
       selected_mu = mu[keep],
       selected_index = keep,
       removed_index = setdiff(seq_len(m), keep),
       effective_count = effective_count,
       mean_weight = effective_count / n,
       cutoff = cutoff,
       settings = list(min_state_count = min_state_count,
                       min_state_fraction = min_state_fraction,
                       screening_prior_sd = screening_prior_sd,
                       block_size = block_size))
}

# Construct the original symmetric power grid and remove unsupported states
# before fitting the HMM. This is also useful for inspecting the automatic
# choice without fitting a model.
ash_hmm_select_grid <- function(y, se,
                                half_grid = 100L,
                                shape = 3,
                                expansion = 1.5,
                                max_abs = NULL,
                                min_state_count = 0.5,
                                min_state_fraction = 1e-5,
                                screening_prior_sd = 0,
                                block_size = 5000L,
                                positive_only = FALSE) {
  .ash_hmm_check_data(y, se, rep(1L, length(y)))
  candidate <- ash_mu_grid(y, half_grid = half_grid, shape = shape,
                           expansion = expansion, max_abs = max_abs,
                           positive_only = positive_only)
  ans <- .ash_hmm_prefilter_mu(
    y, se, candidate,
    min_state_count = min_state_count,
    min_state_fraction = min_state_fraction,
    screening_prior_sd = screening_prior_sd,
    block_size = block_size)
  ans$automatic <- TRUE
  ans
}

ash_hmm_transition_mask <- function(n_states, topology = c("hub", "full")) {
  topology <- match.arg(topology)
  n_states <- as.integer(n_states)
  if (n_states < 1L) stop("'n_states' must be positive.", call. = FALSE)
  if (topology == "full") return(matrix(TRUE, n_states, n_states))

  # Hub-and-spoke: the zero/hub state is state 1. A non-hub state may
  # persist or return to the hub, but may not jump directly to another spoke.
  mask <- diag(TRUE, n_states)
  mask[1L, ] <- TRUE
  mask[, 1L] <- TRUE
  mask
}

.ash_hmm_initial_transition <- function(mask, stay_probability = 0.95) {
  m <- nrow(mask)
  stay_probability <- rep(stay_probability, length.out = m)
  if (any(stay_probability < 0 | stay_probability > 1)) {
    stop("'stay_probability' must lie in [0, 1].", call. = FALSE)
  }
  A <- matrix(0, m, m)
  for (i in seq_len(m)) {
    allowed <- which(mask[i, ])
    if (length(allowed) == 1L) {
      A[i, allowed] <- 1
    } else {
      A[i, i] <- stay_probability[i]
      other <- setdiff(allowed, i)
      A[i, other] <- (1 - stay_probability[i]) / length(other)
    }
  }
  A
}

.ash_hmm_stationary <- function(A, tolerance = 1e-13, maxiter = 100000L) {
  p <- rep(1 / nrow(A), nrow(A))
  for (iter in seq_len(maxiter)) {
    p_new <- drop(p %*% A)
    if (max(abs(p_new - p)) < tolerance) return(.ash_hmm_normalize(p_new))
    p <- p_new
  }
  warning("Stationary-distribution iteration did not converge; using its last iterate.")
  .ash_hmm_normalize(p)
}

.ash_hmm_validate_transition <- function(A, mask, tolerance = 1e-10) {
  m <- nrow(mask)
  if (!is.matrix(A) || any(dim(A) != c(m, m)) || anyNA(A) ||
      any(!is.finite(A)) || any(A < 0)) {
    stop("'init_transition' must be a finite, nonnegative M by M matrix.", call. = FALSE)
  }
  if (any(A[!mask] > tolerance)) {
    stop("'init_transition' has positive mass on a structurally forbidden edge.", call. = FALSE)
  }
  rs <- rowSums(A)
  if (any(abs(rs - 1) > tolerance)) {
    stop("Every row of 'init_transition' must sum to one.", call. = FALSE)
  }
  A[!mask] <- 0
  A / rowSums(A)
}

# Log marginal density for one ash component. In the one-sided model the
# continuous prior is N(mu, tau^2) truncated to [0, Inf). The identity
#
# p(y | beta >= 0) = p_untruncated(y) *
#   Pr(beta >= 0 | y) / Pr(beta >= 0)
#
# gives an exact and stable normal-CDF expression for the convolution with
# N(beta, se^2) observation noise.
.ash_hmm_component_log_emission <- function(
    y, se, mu, tau, effect_support = c("real", "nonnegative")) {
  effect_support <- match.arg(effect_support)
  if (tau == 0) {
    if (effect_support == "nonnegative" && mu < 0) {
      return(rep(-Inf, length(y)))
    }
    return(stats::dnorm(y, mean = mu, sd = se, log = TRUE))
  }

  observation_variance <- se^2 + tau^2
  ans <- stats::dnorm(
    y, mean = mu, sd = sqrt(observation_variance), log = TRUE)
  if (effect_support == "real") return(ans)

  posterior_mean <- (se^2 * mu + tau^2 * y) / observation_variance
  posterior_sd <- tau * se / sqrt(observation_variance)
  ans + stats::pnorm(
    posterior_mean / posterior_sd, log.p = TRUE) -
    stats::pnorm(mu / tau, log.p = TRUE)
}

.ash_hmm_component_log_emission_matrix <- function(
    y, se, mu, prior_sd, effect_support = c("real", "nonnegative")) {
  effect_support <- match.arg(effect_support)
  ans <- matrix(-Inf, length(y), length(prior_sd))
  for (component in seq_along(prior_sd)) {
    ans[, component] <- .ash_hmm_component_log_emission(
      y, se, mu, prior_sd[component], effect_support)
  }
  ans
}

.ash_hmm_log_emission <- function(
    y, se, mu, prior_sd, rho,
    effect_support = c("real", "nonnegative")) {
  effect_support <- match.arg(effect_support)
  n <- length(y)
  m <- length(mu)
  l <- length(prior_sd)
  ans <- matrix(-Inf, n, m)
  log_rho <- .ash_hmm_safe_log(rho)
  for (state in seq_len(m)) {
    terms <- .ash_hmm_component_log_emission_matrix(
      y, se, mu[state], prior_sd, effect_support)
    terms <- sweep(terms, 2L, log_rho[state, ], "+")
    ans[, state] <- .ash_hmm_row_logsumexp(terms)
  }
  ans
}

.ash_hmm_segments <- function(sequence_id) {
  starts <- c(1L, which(sequence_id[-1L] != sequence_id[-length(sequence_id)]) + 1L)
  ends <- c(starts[-1L] - 1L, length(sequence_id))
  list(starts = starts, ends = ends)
}

.ash_hmm_forward_backward <- function(log_emission, A, init_prob, mask,
                                      sequence_id) {
  n <- nrow(log_emission)
  m <- ncol(log_emission)
  log_A <- .ash_hmm_safe_log(A)
  log_pi <- .ash_hmm_safe_log(init_prob)
  edge <- which(mask, arr.ind = TRUE)
  edge_from <- edge[, 1L]
  edge_to <- edge[, 2L]
  incoming <- lapply(seq_len(m), function(j) which(edge_to == j))
  outgoing <- lapply(seq_len(m), function(i) which(edge_from == i))
  segments <- .ash_hmm_segments(sequence_id)

  log_alpha <- matrix(-Inf, n, m)
  log_beta <- matrix(-Inf, n, m)
  gamma <- matrix(0, n, m)
  transition_counts <- matrix(0, m, m)
  boundary_probability <- rep(NA_real_, max(0L, n - 1L))
  total_log_likelihood <- 0

  for (segment in seq_along(segments$starts)) {
    first <- segments$starts[segment]
    last <- segments$ends[segment]
    log_alpha[first, ] <- log_pi + log_emission[first, ]

    if (first < last) {
      for (t in (first + 1L):last) {
        for (state in seq_len(m)) {
          e <- incoming[[state]]
          log_alpha[t, state] <- log_emission[t, state] +
            .ash_hmm_logsumexp(log_alpha[t - 1L, edge_from[e]] +
                                 log_A[cbind(edge_from[e], edge_to[e])])
        }
      }
    }

    segment_ll <- .ash_hmm_logsumexp(log_alpha[last, ])
    if (!is.finite(segment_ll)) {
      stop("The HMM assigns zero probability to an observed sequence.", call. = FALSE)
    }
    total_log_likelihood <- total_log_likelihood + segment_ll
    log_beta[last, ] <- 0

    if (first < last) {
      for (t in (last - 1L):first) {
        for (state in seq_len(m)) {
          e <- outgoing[[state]]
          log_beta[t, state] <- .ash_hmm_logsumexp(
            log_A[cbind(edge_from[e], edge_to[e])] +
              log_emission[t + 1L, edge_to[e]] +
              log_beta[t + 1L, edge_to[e]])
        }
      }
    }

    for (t in first:last) {
      lg <- log_alpha[t, ] + log_beta[t, ] - segment_ll
      lg <- lg - .ash_hmm_logsumexp(lg)
      gamma[t, ] <- exp(lg)
    }

    if (first < last) {
      for (t in first:(last - 1L)) {
        log_xi <- log_alpha[t, edge_from] +
          log_A[cbind(edge_from, edge_to)] +
          log_emission[t + 1L, edge_to] +
          log_beta[t + 1L, edge_to] - segment_ll
        log_xi <- log_xi - .ash_hmm_logsumexp(log_xi)
        xi <- exp(log_xi)
        for (e in seq_along(xi)) {
          transition_counts[edge_from[e], edge_to[e]] <-
            transition_counts[edge_from[e], edge_to[e]] + xi[e]
        }
        boundary_probability[t] <- sum(xi[edge_from != edge_to])
      }
    }
  }

  list(log_likelihood = total_log_likelihood,
       state_probability = gamma,
       transition_counts = transition_counts,
       boundary_probability = boundary_probability,
       initial_counts = colSums(gamma[segments$starts, , drop = FALSE]),
       starts = segments$starts)
}

.ash_hmm_component_statistics <- function(
    y, se, mu, prior_sd, rho, log_emission, gamma,
    effect_support = c("real", "nonnegative")) {
  effect_support <- match.arg(effect_support)
  n <- length(y)
  m <- length(mu)
  l <- length(prior_sd)
  counts <- matrix(0, m, l)
  mean_numerator <- mean_precision <- numeric(m)
  log_rho <- .ash_hmm_safe_log(rho)
  for (state in seq_len(m)) {
    terms <- .ash_hmm_component_log_emission_matrix(
      y, se, mu[state], prior_sd, effect_support)
    terms <- sweep(terms, 2L, log_rho[state, ], "+")
    conditional <- exp(terms - log_emission[, state])
    joint <- conditional * gamma[, state]
    counts[state, ] <- colSums(joint)
    if (effect_support == "real") {
      for (component in seq_len(l)) {
        observation_variance <- se^2 + prior_sd[component]^2
        mean_numerator[state] <- mean_numerator[state] +
          sum(joint[, component] * y / observation_variance)
        mean_precision[state] <- mean_precision[state] +
          sum(joint[, component] / observation_variance)
      }
    }
  }
  list(counts = counts,
       occupancy = rowSums(counts),
       mean_numerator = mean_numerator,
       mean_precision = mean_precision)
}

.ash_hmm_component_counts <- function(
    y, se, mu, prior_sd, rho, log_emission, gamma,
    effect_support = c("real", "nonnegative")) {
  .ash_hmm_component_statistics(y, se, mu, prior_sd, rho,
                                log_emission, gamma,
                                effect_support)$counts
}

# Nonoverlapping intervals around the current set of anchor centers. Restricting
# each learned center to its own interval prevents label crossing and stops a
# dense collection of states from collapsing onto the same plateau. When states
# are pruned, recomputing these intervals automatically gives the survivors more
# room to move.
.ash_hmm_voronoi_bounds <- function(anchor) {
  m <- length(anchor)
  lower <- rep(-Inf, m)
  upper <- rep(Inf, m)
  if (m <= 1L) return(cbind(lower = lower, upper = upper))
  if (anyDuplicated(anchor)) {
    stop("State-center anchors must be distinct when using Voronoi bounds.",
         call. = FALSE)
  }
  order_mu <- order(anchor)
  sorted <- anchor[order_mu]
  midpoint <- (sorted[-m] + sorted[-1L]) / 2
  lower_sorted <- c(-Inf, midpoint)
  upper_sorted <- c(midpoint, Inf)
  lower[order_mu] <- lower_sorted
  upper[order_mu] <- upper_sorted
  cbind(lower = lower, upper = upper)
}

.ash_hmm_update_means <- function(mu, anchor, statistics,
                                  learn, fixed, minimum_count = 2,
                                  eligible = rep(TRUE, length(mu)),
                                  damping = 1,
                                  minimum = rep(-Inf, length(mu)),
                                  bounds = c("voronoi", "none")) {
  bounds <- match.arg(bounds)
  m <- length(mu)
  if (length(eligible) != m || anyNA(eligible)) {
    stop("'eligible' must contain one non-missing value per state.",
         call. = FALSE)
  }
  eligible <- as.logical(eligible)
  damping <- rep(damping, length.out = m)
  if (anyNA(damping) || any(!is.finite(damping)) ||
      any(damping < 0 | damping > 1)) {
    stop("Every mean damping value must lie in [0, 1].", call. = FALSE)
  }
  minimum <- rep(minimum, length.out = m)
  if (anyNA(minimum)) {
    stop("Mean lower bounds must not be missing.", call. = FALSE)
  }
  if (!learn || !m) {
    return(list(mu = mu, raw = mu, updated = rep(FALSE, m),
                lower = rep(-Inf, m), upper = rep(Inf, m)))
  }
  interval <- if (bounds == "voronoi") {
    .ash_hmm_voronoi_bounds(anchor)
  } else {
    cbind(lower = rep(-Inf, m), upper = rep(Inf, m))
  }
  interval[, "lower"] <- pmax(interval[, "lower"], minimum)
  raw <- mu
  estimable <- eligible & statistics$occupancy >= minimum_count &
    is.finite(statistics$mean_precision) & statistics$mean_precision > 0
  estimable[fixed] <- FALSE
  raw[estimable] <- statistics$mean_numerator[estimable] /
    statistics$mean_precision[estimable]
  target <- pmax(interval[, "lower"], pmin(interval[, "upper"], raw))
  updated <- estimable & is.finite(target) & damping > 0
  ans <- mu
  ans[updated] <- mu[updated] +
    damping[updated] * (target[updated] - mu[updated])
  ans[fixed] <- mu[fixed]
  list(mu = ans, raw = raw, updated = updated,
       lower = interval[, "lower"], upper = interval[, "upper"])
}

# Conditional maximization for centers in the one-sided ash model. Component
# responsibilities are held at the current parameter snapshot, so each state's
# criterion is a one-dimensional weighted sum of exact truncated-normal
# component log emissions. Candidate moves are accepted only when this
# conditional criterion does not decrease; damping is backtracked if needed.
.ash_hmm_update_truncated_means <- function(
    y, se, mu, anchor, prior_sd, rho, log_emission, gamma, statistics,
    learn, fixed, minimum_count = 2,
    eligible = rep(TRUE, length(mu)), damping = 1,
    minimum = rep(0, length(mu)), bounds = c("voronoi", "none")) {
  bounds <- match.arg(bounds)
  m <- length(mu)
  if (length(eligible) != m || anyNA(eligible)) {
    stop("'eligible' must contain one non-missing value per state.",
         call. = FALSE)
  }
  eligible <- as.logical(eligible)
  damping <- rep(damping, length.out = m)
  if (anyNA(damping) || any(!is.finite(damping)) ||
      any(damping < 0 | damping > 1)) {
    stop("Every mean damping value must lie in [0, 1].", call. = FALSE)
  }
  minimum <- rep(minimum, length.out = m)
  if (anyNA(minimum) || any(!is.finite(minimum))) {
    stop("One-sided mean lower bounds must be finite.", call. = FALSE)
  }

  interval <- if (bounds == "voronoi") {
    .ash_hmm_voronoi_bounds(anchor)
  } else {
    cbind(lower = rep(-Inf, m), upper = rep(Inf, m))
  }
  interval[, "lower"] <- pmax(interval[, "lower"], minimum)

  # stats::optimize requires finite intervals. This upper cap is deliberately
  # conservative relative to both the starting grid and plausible Gaussian
  # observations; it only affects the outermost anchor cell.
  scale_pad <- max(1e-8, stats::median(se))
  search_upper <- max(c(mu, anchor, y + 8 * se), na.rm = TRUE)
  search_upper <- max(search_upper, max(interval[, "lower"]) + scale_pad)
  infinite_upper <- !is.finite(interval[, "upper"])
  interval[infinite_upper, "upper"] <- pmax(
    search_upper, interval[infinite_upper, "lower"] + scale_pad)

  raw <- ans <- mu
  updated <- rep(FALSE, m)
  criterion_before <- criterion_after <- rep(NA_real_, m)
  if (!learn || !m) {
    return(list(mu = ans, raw = raw, updated = updated,
                lower = interval[, "lower"], upper = interval[, "upper"],
                criterion_before = criterion_before,
                criterion_after = criterion_after))
  }

  estimable <- eligible & statistics$occupancy >= minimum_count & damping > 0
  estimable[fixed] <- FALSE
  log_rho <- .ash_hmm_safe_log(rho)

  for (state in which(estimable)) {
    component_log_current <- .ash_hmm_component_log_emission_matrix(
      y, se, mu[state], prior_sd, "nonnegative")
    terms <- sweep(
      component_log_current, 2L, log_rho[state, ], "+")
    component_probability <- exp(terms - log_emission[, state])
    joint <- component_probability * gamma[, state]
    positive_weight <- joint > 0
    if (!any(positive_weight)) next

    conditional_objective <- function(center) {
      component_log <- .ash_hmm_component_log_emission_matrix(
        y, se, center, prior_sd, "nonnegative")
      sum(joint[positive_weight] * component_log[positive_weight])
    }

    lower <- interval[state, "lower"]
    upper <- interval[state, "upper"]
    if (!is.finite(lower) || !is.finite(upper) || upper <= lower) next
    current <- min(upper, max(lower, mu[state]))
    current_value <- conditional_objective(current)
    optimum <- stats::optimize(
      function(center) -conditional_objective(center),
      interval = c(lower, upper),
      tol = sqrt(.Machine$double.eps))
    candidates <- unique(c(current, lower, optimum$minimum, upper))
    values <- vapply(candidates, conditional_objective, numeric(1L))
    target <- candidates[which.max(values)]
    target_value <- max(values)
    raw[state] <- target
    criterion_before[state] <- current_value

    if (!is.finite(target_value) ||
        target_value <= current_value + 1e-10 * (1 + abs(current_value))) {
      criterion_after[state] <- current_value
      next
    }

    fraction <- damping[state]
    proposed <- current + fraction * (target - current)
    proposed_value <- conditional_objective(proposed)
    while (fraction > 2^-20 &&
           proposed_value < current_value - 1e-10 * (1 + abs(current_value))) {
      fraction <- fraction / 2
      proposed <- current + fraction * (target - current)
      proposed_value <- conditional_objective(proposed)
    }
    if (is.finite(proposed_value) && proposed_value >= current_value -
        1e-10 * (1 + abs(current_value))) {
      ans[state] <- proposed
      updated[state] <- abs(proposed - mu[state]) > 1e-12
      criterion_after[state] <- proposed_value
    } else {
      criterion_after[state] <- current_value
    }
  }

  ans[fixed] <- mu[fixed]
  list(mu = ans, raw = raw, updated = updated,
       lower = interval[, "lower"], upper = interval[, "upper"],
       criterion_before = criterion_before,
       criterion_after = criterion_after)
}

# Effective dimension used by the optional BIC/extended-BIC comparison with
# the strict all-null model. Only numerically active simplex coordinates count;
# this is a pragmatic diagnostic because ordinary BIC regularity fails on
# mixture boundaries.
.ash_hmm_effective_parameter_count <- function(
    A, mask, rho, fixed_pointmass_states = integer(),
    learned_mean_states = integer(), estimate_init = FALSE,
    init_prob = NULL, tolerance = 1e-6) {
  active_transition <- mask & A > tolerance
  transition_df <- sum(pmax(rowSums(active_transition) - 1L, 0L))
  mixture_df <- 0L
  free <- setdiff(seq_len(nrow(rho)), fixed_pointmass_states)
  if (length(free)) {
    mixture_df <- sum(pmax(rowSums(rho[free, , drop = FALSE] > tolerance) - 1L, 0L))
  }
  mean_df <- length(unique(learned_mean_states))
  init_df <- 0L
  if (estimate_init && !is.null(init_prob)) {
    init_df <- max(sum(init_prob > tolerance) - 1L, 0L)
  }
  as.integer(transition_df + mixture_df + mean_df + init_df)
}

# Decode a state path with an explicit l0 penalty for every state change. This
# separates step-count selection from the fitted Markov transition matrix.
.ash_hmm_penalized_viterbi <- function(log_emission, mask, sequence_id,
                                       change_penalty) {
  n <- nrow(log_emission)
  m <- ncol(log_emission)
  if (!is.numeric(change_penalty) || length(change_penalty) != 1L ||
      !is.finite(change_penalty) || change_penalty < 0) {
    stop("'change_penalty' must be a finite nonnegative scalar.", call. = FALSE)
  }
  segments <- .ash_hmm_segments(sequence_id)
  path <- integer(n)
  for (segment in seq_along(segments$starts)) {
    first <- segments$starts[segment]
    last <- segments$ends[segment]
    len <- last - first + 1L
    delta <- matrix(-Inf, len, m)
    back <- matrix(NA_integer_, len, m)
    delta[1L, ] <- log_emission[first, ]
    if (len > 1L) {
      for (tt in 2L:len) {
        for (state in seq_len(m)) {
          from <- which(mask[, state])
          value <- delta[tt - 1L, from] -
            change_penalty * as.numeric(from != state)
          winner <- which.max(value)
          delta[tt, state] <- log_emission[first + tt - 1L, state] + value[winner]
          back[tt, state] <- from[winner]
        }
      }
    }
    path[last] <- which.max(delta[len, ])
    if (len > 1L) {
      for (t in last:(first + 1L)) {
        path[t - 1L] <- back[t - first + 1L, path[t]]
      }
    }
  }
  path
}

.ash_hmm_step_summary <- function(path, sequence_id, mu = NULL) {
  n <- length(path)
  segments <- .ash_hmm_segments(sequence_id)
  rows <- vector("list", 0L)
  per_sequence <- data.frame(
    sequence = character(), steps = integer(), changes = integer(),
    stringsAsFactors = FALSE)
  for (s in seq_along(segments$starts)) {
    first <- segments$starts[s]
    last <- segments$ends[s]
    local_change <- if (first < last) {
      which(path[(first + 1L):last] != path[first:(last - 1L)]) + first
    } else integer()
    starts <- c(first, local_change)
    ends <- c(local_change - 1L, last)
    for (j in seq_along(starts)) {
      state <- path[starts[j]]
      rows[[length(rows) + 1L]] <- data.frame(
        sequence = as.character(sequence_id[first]),
        start = starts[j], end = ends[j], length = ends[j] - starts[j] + 1L,
        state = state,
        state_mean = if (is.null(mu)) NA_real_ else mu[state],
        stringsAsFactors = FALSE)
    }
    per_sequence <- rbind(
      per_sequence,
      data.frame(sequence = as.character(sequence_id[first]),
                 steps = length(starts), changes = length(starts) - 1L,
                 stringsAsFactors = FALSE))
  }
  segment_table <- if (length(rows)) do.call(rbind, rows) else data.frame()
  list(segments = segment_table,
       per_sequence = per_sequence,
       step_count = sum(per_sequence$steps),
       change_count = sum(per_sequence$changes),
       occupied_state_count = length(unique(path)))
}

# Mark the lower-occupancy member of each very close pair as a possible pruning
# candidate. Actual deletion is accepted only after the caller checks the HMM
# marginal likelihood.
.ash_hmm_close_state_candidates <- function(mu, occupancy, distance,
                                            protected = integer()) {
  if (length(mu) < 2L || !is.finite(distance) || distance <= 0) {
    return(integer())
  }
  ord <- order(mu)
  close <- which(diff(mu[ord]) <= distance)
  if (!length(close)) return(integer())
  drop <- integer()
  for (j in close) {
    pair <- ord[c(j, j + 1L)]
    eligible <- setdiff(pair, protected)
    if (!length(eligible)) next
    if (length(eligible) == 1L) {
      drop <- c(drop, eligible)
    } else {
      drop <- c(drop, eligible[which.min(occupancy[eligible])])
    }
  }
  unique(drop)
}

.ash_hmm_expand_prior <- function(x, nr, nc, name) {
  if (length(x) == 1L) x <- matrix(x, nr, nc)
  else if (length(x) == nc) x <- matrix(rep(x, each = nr), nr, nc)
  else if (is.matrix(x) && all(dim(x) == c(nr, nc))) x <- x
  else stop(sprintf("'%s' must be scalar, length %d, or a %d by %d matrix.",
                    name, nc, nr, nc), call. = FALSE)
  if (anyNA(x) || any(!is.finite(x)) || any(x < 1)) {
    stop(sprintf("Every '%s' entry must be finite and at least one.", name),
         call. = FALSE)
  }
  x
}

.ash_hmm_dirichlet_penalty <- function(probability, alpha, active = NULL) {
  if (is.null(active)) active <- matrix(TRUE, nrow(probability), ncol(probability))
  use <- active & alpha > 1
  if (!any(use)) return(0)
  if (any(probability[use] <= 0)) return(-Inf)
  sum((alpha[use] - 1) * log(probability[use]))
}

.ash_hmm_viterbi <- function(log_emission, A, init_prob, mask, sequence_id) {
  n <- nrow(log_emission)
  m <- ncol(log_emission)
  log_A <- .ash_hmm_safe_log(A)
  log_pi <- .ash_hmm_safe_log(init_prob)
  edge <- which(mask, arr.ind = TRUE)
  incoming <- lapply(seq_len(m), function(j) which(edge[, 2L] == j))
  segments <- .ash_hmm_segments(sequence_id)
  path <- integer(n)

  for (segment in seq_along(segments$starts)) {
    first <- segments$starts[segment]
    last <- segments$ends[segment]
    delta <- matrix(-Inf, last - first + 1L, m)
    back <- matrix(NA_integer_, last - first + 1L, m)
    delta[1L, ] <- log_pi + log_emission[first, ]
    if (first < last) {
      for (tt in 2L:(last - first + 1L)) {
        for (state in seq_len(m)) {
          e <- incoming[[state]]
          from <- edge[e, 1L]
          values <- delta[tt - 1L, from] + log_A[cbind(from, rep(state, length(from)))]
          winner <- which.max(values)
          delta[tt, state] <- log_emission[first + tt - 1L, state] + values[winner]
          back[tt, state] <- from[winner]
        }
      }
    }
    path[last] <- which.max(delta[nrow(delta), ])
    if (first < last) {
      for (t in last:(first + 1L)) {
        path[t - 1L] <- back[t - first + 1L, path[t]]
      }
    }
  }
  path
}

# Mean and variance of a standard normal truncated below at -z. For very
# negative z, direct inverse-Mills calculations suffer catastrophic
# cancellation; the displayed tail expansions retain accuracy.
.ash_hmm_lower_truncated_standard_moments <- function(z) {
  mean_shift <- variance_factor <- numeric(length(z))
  extreme <- z < -10

  if (any(!extreme)) {
    zz <- z[!extreme]
    inverse_mills <- exp(
      stats::dnorm(zz, log = TRUE) - stats::pnorm(zz, log.p = TRUE))
    mean_shift[!extreme] <- inverse_mills
    variance_factor[!extreme] <- pmax(
      0, 1 - zz * inverse_mills - inverse_mills^2)
  }

  if (any(extreme)) {
    a <- -z[extreme]
    inverse_a <- 1 / a
    # lambda(-a) - a and Var{Z | Z >= a}; terms through a^-7/a^-6.
    residual <- inverse_a - 2 * inverse_a^3 +
      10 * inverse_a^5 - 74 * inverse_a^7
    mean_shift[extreme] <- a + residual
    variance_factor[extreme] <- pmax(
      0, inverse_a^2 - 6 * inverse_a^4 + 50 * inverse_a^6)
  }

  list(mean_shift = mean_shift, variance_factor = variance_factor)
}

.ash_hmm_posterior_summary <- function(
    y, se, mu, prior_sd, rho, log_emission, gamma,
    effect_support = c("real", "nonnegative")) {
  effect_support <- match.arg(effect_support)
  n <- length(y)
  m <- length(mu)
  l <- length(prior_sd)
  posterior_mean <- posterior_second <- prob_ge_zero <- prob_le_zero <-
    prob_zero <- numeric(n)
  log_rho <- .ash_hmm_safe_log(rho)

  for (state in seq_len(m)) {
    terms <- .ash_hmm_component_log_emission_matrix(
      y, se, mu[state], prior_sd, effect_support)
    terms <- sweep(terms, 2L, log_rho[state, ], "+")
    component_probability <- exp(terms - log_emission[, state])
    state_mean <- state_second <- state_ge <- state_le <- state_zero <- numeric(n)

    for (component in seq_len(l)) {
      tau <- prior_sd[component]
      if (tau == 0) {
        component_mean <- rep(mu[state], n)
        component_variance <- rep(0, n)
        ge <- rep(as.numeric(mu[state] >= 0), n)
        le <- rep(as.numeric(mu[state] <= 0), n)
        p0 <- rep(as.numeric(mu[state] == 0), n)
      } else {
        component_variance <- tau^2 * se^2 / (tau^2 + se^2)
        component_mean <- (se^2 * mu[state] + tau^2 * y) / (tau^2 + se^2)
        component_sd <- sqrt(component_variance)
        if (effect_support == "nonnegative") {
          z <- component_mean / component_sd
          truncated <- .ash_hmm_lower_truncated_standard_moments(z)
          component_mean <- component_mean +
            component_sd * truncated$mean_shift
          component_variance <- component_variance *
            truncated$variance_factor
          component_mean <- pmax(0, component_mean)
          ge <- rep(1, n)
          le <- rep(0, n)
        } else {
          ge <- stats::pnorm(component_mean / component_sd)
          le <- stats::pnorm(-component_mean / component_sd)
        }
        p0 <- rep(0, n)
      }
      w <- component_probability[, component]
      state_mean <- state_mean + w * component_mean
      state_second <- state_second + w * (component_variance + component_mean^2)
      state_ge <- state_ge + w * ge
      state_le <- state_le + w * le
      state_zero <- state_zero + w * p0
    }

    g <- gamma[, state]
    posterior_mean <- posterior_mean + g * state_mean
    posterior_second <- posterior_second + g * state_second
    prob_ge_zero <- prob_ge_zero + g * state_ge
    prob_le_zero <- prob_le_zero + g * state_le
    prob_zero <- prob_zero + g * state_zero
  }

  if (effect_support == "nonnegative") {
    posterior_mean <- pmax(0, posterior_mean)
    prob_ge_zero <- rep(1, n)
  }
  posterior_variance <- pmax(0, posterior_second - posterior_mean^2)
  list(mean = posterior_mean,
       sd = sqrt(posterior_variance),
       probability_ge_zero = pmin(1, pmax(0, prob_ge_zero)),
       probability_le_zero = pmin(1, pmax(0, prob_le_zero)),
       probability_zero = pmin(1, pmax(0, prob_zero)),
       lfsr = pmin(prob_ge_zero, prob_le_zero))
}




#' @export
print.ash_hmm_fit <- function(x, ...) {
  cat("Adaptive-shrinkage HMM fit\n")
  cat("  observations:", nrow(x$state_probability), "\n")
  cat("  states:", ncol(x$state_probability))
  if (isTRUE(x$grid$automatic)) {
    cat(" (retained from", length(x$grid$full_mu), "automatic candidates)")
  }
  cat("\n")
  if (!is.null(x$fitted$effect_support)) {
    cat("  effect support:", x$fitted$effect_support, "\n")
  }
  if (isTRUE(x$fitted$learn_state_means)) {
    moved <- sum(abs(x$fitted$mu - x$fitted$mean_anchor) > 1e-10)
    cat("  learned state means:", moved, "of", length(x$fitted$mu), "\n")
  }
  if (!is.null(x$pruning_history) && nrow(x$pruning_history)) {
    cat("  dynamically pruned states:", nrow(x$pruning_history), "\n")
  }
  if (!is.null(x$model_selection) &&
      isTRUE(x$model_selection$collapsed_to_null)) {
    cat("  global model selection: strict null\n")
  }
  if (!is.null(x$step_selection)) {
    cat("  penalized steps:", x$step_selection$step_count,
        "(changes:", x$step_selection$change_count, ")\n")
  }
  cat("  iterations:", x$iterations,
      if (x$converged) "(converged)" else "(maximum reached)", "\n")
  cat("  log likelihood:", format(x$log_likelihood, digits = 8), "\n")
  invisible(x)
}
