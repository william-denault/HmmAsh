library(ashr)
'%!in%' <- function(x,y)!('%in%'(x,y))
fit_hmm <- function (x, sd,
                     halfK=100,
                     mult=3,
                     smooth=FALSE,
                     thresh=0.00001,
                     prefilter=TRUE,
                     thresh_prefilter=1e-30,
                     maxiter=3,
                     max_zscore=20,
                     thresh_sd=1e-30,
                     epsilon=1e-2) {

  # Deal with case where very close to zero sds
  if( length(which(sd < thresh_sd)) > 0){
    sd[ which(sd < thresh_sd)] <- thresh_sd
  }
  if (sum(is.na(sd)) > 0){
    x [ which( is.na(sd))] <- 0
    sd[ which( is.na(x))] <- 1
  }
  if(sum(!is.finite(sd)) > 0){
    x [which(!is.finite(sd))] = 0
    sd[which(!is.finite(sd))] = 1
  }

  if( length(which(abs(x/sd) > max_zscore)) > 0){ # avoid underflow a z-score of 20 => pv < e-90
    sd[which(abs(x/sd) > max_zscore)] <- abs(x[which(abs(x/sd) > max_zscore)]) / max_zscore
  }

  K = 2*halfK-1
  sd = sd
  X <- x

  pos <- seq(0, 1, length.out=halfK)

  # Define the mean states
  mu <- (pos^(1/mult))*1.5*max(abs(X)) # put 0 state at the first place
  mu <- c(mu, -mu[-1] )

  min_delta <- abs(mu[2]-mu[1])
  if( prefilter ){
    tt <- apply(
      do.call( rbind, lapply(1:length(x), function( i){
        tt <- dnorm(x[i], mean = mu, sd=sd[i])
        return( tt / sum(tt))
      } )),
      2,
      mean, na.rm=TRUE)

    temp_idx <- which(tt > thresh_prefilter)
    if( 1 %!in% temp_idx){
      temp_idx <- c(1, temp_idx)
    }
    mu <- mu[temp_idx]
    K <- length(mu)
  }

  # Enforce hub-and-spoke topology strictly in initial P
  P <- matrix(0, ncol=K, nrow=K)
  diag(P) <- 0.5 # Self-transitions

  if (K > 1) {
    P[1, 2:K] <- 0.5 / (K - 1) # Null to non-nulls
    P[2:K, 1] <- 0.5           # Non-nulls to null
  }

  P <- P + matrix(epsilon, ncol=K, nrow=K)
  # Keep structural zeros for non-null to non-null
  if (K > 1) {
    for (i in 2:K) {
      for (j in 2:K) {
        if (i != j) P[i, j] <- 0
      }
    }
  }
  P <- P / rowSums(P) # Exact row normalization

  # Strong prior on starting in the null state to prevent boundary false positives
  pi <- rep(epsilon, K)
  if (K > 1) {
    pi[1] <- 1 - sum(pi[-1])
  } else {
    pi[1] <- 1
  }

  emit = function(k, x, t){
    dnorm(x, mean=mu[k], sd=sd[t])
  }

  alpha_hat = matrix(nrow = length(X), ncol=K)
  alpha_tilde = matrix(nrow = length(X), ncol=K)
  G_t <- rep(NA, length(X))

  for(k in 1:K){
    alpha_hat[1, ] = pi * emit(1:K, x=X[1], t=1)
    alpha_tilde[1, ] = pi * emit(1:K, x=X[1], t=1)
  }

  # Initial Forward algorithm
  for(t in 1:(length(X)-1)){
    m = alpha_hat[t,] %*% P
    alpha_tilde[t+1, ] = m * emit(1:K, x=X[t+1], t= t+1 )
    G_t[t+1] <- sum( alpha_tilde[t+1,])
    alpha_hat[t+1,] <-  alpha_tilde[t+1,] / ( G_t[t+1])
  }

  beta_hat = matrix(nrow = length(X), ncol=K)
  beta_tilde = matrix(nrow = length(X), ncol=K)
  C_t <- rep(NA, length(X))

  # Initialize beta
  for(k in 1:K){
    beta_hat[ length(X), k] = 1
    beta_tilde [ length(X), k] = 1
  }

  # Initial Backwards algorithm
  for(t in (length(X)-1):1){
    emissio_p <- emit(1:K, X[t+1], t=t+1)
    beta_tilde [t, ] = apply( sweep( P, 2, beta_hat[t+1,]*emissio_p, "*" ), 1, sum)
    C_t[t] <- max(beta_tilde[t,])
    beta_hat[t,] <- beta_tilde [t, ] / C_t[t]
  }

  ab = alpha_hat*beta_hat
  prob = ab/rowSums(ab)

  xi <- array(0, dim = c(K, K))
  for (t in 1:(length(X) - 1)) {
    xi_t <- outer(alpha_hat[t, ], beta_hat[t+1, ] * emit(1:K, X[t+1], t+1)) * P
    xi_t <- xi_t / sum(xi_t)
    xi <- xi + xi_t
  }

  # Enforce structural zeros in first pre-ASH P update
  if (K > 1) {
    for (i in 2:K) {
      for (j in 2:K) {
        if (i != j) xi[i, j] <- 0
      }
    }
  }
  row_sums <- rowSums(xi)
  row_sums[row_sums == 0] <- 1  # prevent division by zero
  P <- xi / row_sums
  P[P < epsilon & P > 0] <- epsilon
  P <- P / rowSums(P)

  idx_comp <- which( apply(prob, 2, mean) > thresh )
  if ( !(1 %in% idx_comp) ){ # ensure 0 is in the model
    idx_comp <- c(1, idx_comp)
  }

  ash_obj <- list()
  x_post <- 0*x

  for (i in 2:length(idx_comp)){
    mu_ash <- mu[idx_comp[i]]
    weight <- prob[, idx_comp[i]]

    ash_obj[[i]] <- ash(x, sd,
                        weight=weight,
                        mode=mu_ash,
                        mixcompdist = "normal")
    x_post <- x_post + weight*ash_obj[[i]]$result$PosteriorMean
  }

  prob <- prob[ ,idx_comp, drop=FALSE]
  iter = 1
  K = length(idx_comp)

  # Re-normalize P after subsetting
  P = P[idx_comp, idx_comp, drop=FALSE]
  P <- P / rowSums(P)

  while( iter < maxiter ){

    alpha_hat = matrix(nrow = length(X), ncol=K)
    alpha_tilde = matrix(nrow = length(X), ncol=K)
    G_t <- rep(NA, length(X))

    # Proper E-step initialization for t=1
    pi_iter <- prob[1, ]
    pi_iter[pi_iter < epsilon] <- epsilon # prevent log(0) issues
    pi_iter <- pi_iter / sum(pi_iter)

    data0_t1 <- set_data(X[1], sd[1])
    emissio_p1 <- c(dnorm(X[1], mean=0, sd=sd[1]),
                    sapply(2:K, function(k) exp(ashr::calc_loglik(ash_obj[[k]], data0_t1))))

    alpha_tilde[1, ] <- pi_iter * emissio_p1
    G_t[1] <- sum(alpha_tilde[1, ])
    alpha_hat[1, ] <- alpha_tilde[1, ] / G_t[1]

    # Forward algorithm
    for(t in 1:(length(X)-1)){
      m = alpha_hat[t,] %*% P
      data0 <- set_data(X[t+1], sd[t+1])

      alpha_tilde[t+1, ] = m * c(dnorm(X[t+1], mean=0, sd=sd[t+1]),
                                 sapply(2:K, function(k) exp(ashr::calc_loglik(ash_obj[[k]], data0))))

      G_t[t+1] <- sum( alpha_tilde[t+1,])
      alpha_hat[t+1,] <- alpha_tilde[t+1,] / ( G_t[t+1])
    }

    beta_hat = matrix(nrow = length(X), ncol=K)
    beta_tilde = matrix(nrow = length(X), ncol=K)
    C_t <- rep(NA, length(X))

    # Initialize beta
    for(k in 1:K){
      beta_hat[ length(X), k] = 1
      beta_tilde [ length(X), k] = 1
    }

    # Backwards algorithm
    for(t in (length(X)-1):1){
      data0 <- set_data(X[t+1], sd[t+1])
      emissio_p <- c(dnorm(X[t+1], mean=0, sd=sd[t+1]),
                     sapply(2:K, function(k) exp(ashr::calc_loglik(ash_obj[[k]], data0))))

      beta_tilde [t, ] = apply( sweep( P, 2, beta_hat[t+1,]*emissio_p, "*" ), 1, sum)
      C_t[t] <- max(beta_tilde[t,])
      beta_hat[t,] <- beta_tilde [t, ] / C_t[t]
    }

    ab = alpha_hat*beta_hat
    prob = ab/rowSums(ab)

    ash_obj <- list()
    x_post <- 0*x

    for ( k in 2:K){
      mu_ash <- mu[k]
      weight <- prob[, k]

      ash_obj[[k]] <- ash(x, sd,
                          weight=weight,
                          mode=mu_ash,
                          mixcompdist = "normal")
      x_post <- x_post + weight*ash_obj[[k]]$result$PosteriorMean
    }

    # Baum_Welch updates for transition matrix
    ab = alpha_hat*beta_hat
    prob = ab/rowSums(ab)

    xi <- array(0, dim = c(K, K))
    for (t in 1:(length(X) - 1)) {
      data0 <- set_data(X[t+1], sd[t+1])
      emissio_p_xi <- c(dnorm(X[t+1], mean=0, sd=sd[t+1]),
                        sapply(2:K, function(k) exp(ashr::calc_loglik(ash_obj[[k]], data0))))

      xi_t <- outer(alpha_hat[t, ], beta_hat[t+1, ] * emissio_p_xi) * P
      xi_t <- xi_t / sum(xi_t)
      xi <- xi + xi_t
    }

    # Enforce structural zeros in main Baum-Welch P update
    if (K > 1) {
      for (i in 2:K) {
        for (j in 2:K) {
          if (i != j) xi[i, j] <- 0
        }
      }
    }

    row_sums <- rowSums(xi)
    row_sums[row_sums == 0] <- 1  # prevent division by zero
    P <- xi / row_sums

    # Ensure no probabilities drop exactly to 0 (except structural zeros)
    P[P < epsilon & P > 0] <- epsilon
    P <- P / rowSums(P)

    iter = iter + 1
  }

  ## ---------------------------------------------------------------
  ## Local false sign rate under the HMM mixture posterior.
  ##
  ## At position t, the posterior on beta_t is the mixture
  ##   p(beta_t | x) = prob[t,1] * delta_0
  ##                 + sum_{k>=2} prob[t,k] * g_k_post(beta_t | x_t)
  ## where g_k_post is the per-observation ash posterior under prior g_k.
  ##
  ## The lfsr is
  ##   lfsr(t) = P(beta_t = 0 | x) + min(P(beta_t > 0 | x), P(beta_t < 0 | x))
  ##           = 1 - max(P(beta_t > 0 | x), P(beta_t < 0 | x))
  ## with
  ##   P(beta_t > 0 | x) = sum_{k>=2} prob[t,k] * ash_obj[[k]]$PositiveProb[t]
  ##   P(beta_t < 0 | x) = sum_{k>=2} prob[t,k] * ash_obj[[k]]$NegativeProb[t]
  ## (state 1 is delta_0, which contributes 0 to both.)
  ##
  ## The previous version
  ##     lfsr <- prob[,1] + sum_k prob[,k] * ash_obj[[k]]$lfsr
  ## sums per-state lfsr values.  Each per-state lfsr is computed against
  ## the natural sign of g_k's mode, so when posterior mass at position t
  ## is split across non-null states with opposite signs the per-state
  ## values are individually small (confident in opposite signs) and the
  ## sum underestimates the true lfsr.  The reformulation above is exact
  ## under the HMM mixture and matches the convention used by ashr.
  ## ---------------------------------------------------------------
  P_pos <- rep(0, length(x))
  P_neg <- rep(0, length(x))
  for ( k in 2:K){
    pp <- ash_obj[[k]]$result$PositiveProb
    pn <- ash_obj[[k]]$result$NegativeProb
    if (is.null(pp) || is.null(pn)) {
      ## Fallback for older ashr that doesn't expose PositiveProb/NegativeProb:
      ## derive them from PosteriorMean / PosteriorSD assuming a normal post.
      pm <- ash_obj[[k]]$result$PosteriorMean
      ps <- pmax(ash_obj[[k]]$result$PosteriorSD, .Machine$double.eps)
      pp <- 1 - pnorm(0, mean = pm, sd = ps)
      pn <- pnorm(0, mean = pm, sd = ps)
    }
    P_pos <- P_pos + prob[, k] * pp
    P_neg <- P_neg + prob[, k] * pn
  }
  lfsr_est <- pmin(1, pmax(0, 1 - pmax(P_pos, P_neg)))

  # --- NEW: Log Bayes Factor Calculation for Global Testing ---

  # 1. Log-likelihood of the full HMM model
  # The scaling factors (G_t) from the forward algorithm give the marginal likelihood of the sequence
  ll_hmm <- sum(log(G_t))

  # 2. Log-likelihood of the strict null model (mean = 0 everywhere)
  ll_null <- sum(dnorm(X, mean = 0, sd = sd, log = TRUE))

  # 3. Log Bayes Factor
  log_BF <- ll_hmm - ll_null
  # ------------------------------------------------------------
  #print(log_BF)
  out <- list(prob = prob,
              x_post = x_post,
              lfsr = lfsr_est,
              mu = mu,
              ll_hmm = ll_hmm,      # NEW
              ll_null = ll_null,    # NEW
              log_BF = log_BF)      # NEW

  return(out)
}
