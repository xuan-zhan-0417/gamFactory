########################
# Multi-start backfitting initialization for single-index ("si") nested effects.
#
# For a gam_nl() model, each "special" (non-standard) effect has coefficients split into
# an inner block (info$iec[[idx]][1:na], the single-index direction "alpha") and an outer
# block (the outer spline coefficients "beta"). For a single index (bundle_trans_si, R/bundle_trans_si.R) the
# effective single-index value is ax = Xi %*% (alpha + a0), Xi and a0 being stored in
# info$extra[[idx]]$si$margin[[1]], with a0 a fixed offset. So for any desired
# *effective* direction v, the coefficient values to place are simply v - a0. We exploit
# this to treat "resolved via search", "resolved via the historical default", and "the
# current candidate under test" uniformly (see .fit_reduced_si below).
#
# Only nested effects made of one single index margin (bundle_trans_si) get the sphere/eigenvector
# multi-start search; all other nested effects (exp, mgks, si_nexp, positive single indices, two margins)
# keep the historical behaviour (use the construction-time default alpha as-is, leaving their
# already-constructed outer basis columns untouched).

# .sd_n (population sd) is defined in R/I_nest_bundle.R; the si bundle rescales alpha so that it equals 1
# for the default direction.

# a0 (the fixed additive offset in alpha + a0) only exists for single index margins; it is zero for the
# parameters of the other margins, so that v = alpha + a0 / alpha = v - a0 both correctly degenerate to the
# raw alpha for them. Returns one value per inner parameter of the effect, margin after margin.
.si_a0 <- function(extra) {
  unlist(lapply(extra$si$margin, function(mk) if (is.null(mk$a0)) numeric(length(mk$idx)) else mk$a0))
}

########################
# Build a matrix of n_init candidate directions for a single-index effect with inner
# design matrix Xi (n x r). Every row v is rescaled so that Var(Xi %*% v) == 1.
# Testing code
# Xi <- matrix(rnorm(1000 * 20), nrow = 1000, ncol = 20)
# Xi <- scale(Xi, scale = F)
# alpha <- rep(1, 20)
# 
# # Should give a warning
# v <- .si_search_candidates(Xi, alpha, n_init = 1000, n_eigen = min(ncol(Xi), 10), oversample = 10)
# 
# range(apply(v, 1, function(x) sd(Xi %*% x) * sqrt(999/1000)))
#
.si_search_candidates <- function(Xi, alpha, n_init = 1000, n_eigen = min(ncol(Xi), 10), oversample = 10, seed = 1) {

  d <- ncol(Xi)
  
  if(abs(.sd_n(Xi %*% alpha) - 1) > 1e-6) {
    warning("Initial alpha does not have unit variance; rescaling to unit variance.")
    alpha <- alpha / .sd_n(Xi %*% alpha)
  }
  
  if( any(abs(colMeans(Xi)) > 1e-6) ) {
   stop("Covariates must be centered (mean zero) for single-index search.")
  }

  # Initialisation provided by bundle_trans_si (R/bundle_trans_si.R). It is such that 
  # Var(Xi %*% alpha) == 1, so we can use it as-is as the first candidate.
  cand <- t(as.vector(alpha))
  
  if(n_init > 1){
    XtX <- crossprod(Xi)
    
    # Add eigenvectors of the covariates cross-product. These will be rescaled, but not rotated!
    n_eigen <- min(n_eigen, d)
    ev <- eigen(XtX, symmetric = TRUE)$vectors[ , 1:n_eigen, drop = FALSE]
    ev <- t(ev) / sqrt(colMeans( (Xi %*% ev)^2 ))
    
    cand <- rbind(cand, ev)
    
    # Add vectors on a half-sphere (to avoid the sign ambiguity of the single-index direction)
    # Makes sure that the variance of the single index is 1 for all candidates, i.e. Var(Xi %*% v) == 1.
    # This is done using the Cholesky factor of the covariates cross-product, so that the half-sphere is 
    # stretched to match the covariates' geometry. This rotates the simulated vectors, while we do not want
    # to rotate the eigenvectors. 
    n_left <- n_init - nrow(cand)
    if(n_left > 0) {
      vsim <- .grid_on_half_sphere(N = n_left, r = d, oversample = oversample, seed = seed)
      ch <- chol(crossprod(Xi))
      vsim <- t( backsolve(ch, t(vsim)) ) * sqrt(nrow(Xi))
      cand <- rbind(cand, vsim)
    }

  }

  return( cand )
}

########################
# Given the set of all special-effect indices (nested_idx), a named list `resolved`
# (as.character(idx) -> effective direction v) for the effects currently fixed / under
# test, and `reeval_idx` (the subset of `resolved` whose v may differ from the
# construction-time default and therefore need their outer basis recomputed - i.e. the
# "si" effects currently being searched, never other nested effects), build the
# reduced x/E/lpi handed to family$initialize_bundle():
#  - effects in `reeval_idx` get their outer/basis columns overwritten via
#    basis$evalX(Xi %*% v)$X0, and their inner columns dropped;
#  - effects in `resolved` but not in `reeval_idx` (all other nested effects) keep their
#    already-constructed outer columns untouched (they were built at the
#    construction-time default, which is exactly what `resolved` holds for them), and
#    only their inner columns are dropped;
#  - effects not yet in `resolved` are dropped entirely (both inner and outer columns).
.build_reduced_design_for_si_init <- function(x, E, lpi, info, nested_idx, resolved, reeval_idx) {

  drop_idx <- integer(0)

  for (idx in nested_idx) {

    iec <- info$iec[[idx]]
    extra <- info$extra[[idx]]
    na <- length(extra$si$alpha)
    inner <- iec[1:na]
    outer <- iec[-(1:na)]

    key <- as.character(idx)
    if (is.null(resolved[[key]])) {
      drop_idx <- c(drop_idx, iec)
      next
    }

    drop_idx <- c(drop_idx, inner)

    if (idx %in% reeval_idx) {
      v <- resolved[[key]]
      ax <- drop( extra$si$margin[[1]]$X %*% v )
      x[ , outer] <- extra$basis$evalX(z1 = ax, deriv = 0)$X0
    }
  }

  keep_idx <- setdiff(seq_len(ncol(x)), drop_idx)

  lpi_new <- lapply(lpi, function(.x) match(.x[!(.x %in% drop_idx)], keep_idx))

  list(x = x[ , keep_idx, drop = FALSE],
       E = E[ , keep_idx, drop = FALSE],
       lpi = lpi_new,
       keep_idx = keep_idx)
}

########################
# Build the reduced design for the current `resolved` set, run family$initialize_bundle
# on it, and scatter the result into a full-length coefficient vector (zero everywhere
# by default, so any not-yet-resolved "si" effect contributes exactly zero to the linear
# predictor - see eff_nest.R / linpreds.R).
.fit_reduced_si <- function(y, nobs, x, E, lpi, info, family, offset, weights, unscaled, nested_idx, resolved, reeval_idx) {

  reduced <- .build_reduced_design_for_si_init(x = x, E = E, lpi = lpi, info = info,
                                                nested_idx = nested_idx, resolved = resolved,
                                                reeval_idx = reeval_idx)

  attr(reduced$x, "lpi") <- reduced$lpi

  start_reduced <- family$initialize_bundle(y = y, nobs = nobs, E = reduced$E, x = reduced$x,
                                             family = family, offset = offset, jj = reduced$lpi,
                                             unscaled = unscaled, weights = weights)

  start <- numeric( ncol(x) )
  start[reduced$keep_idx] <- start_reduced

  # Need to subtract a0 from the resolved effective direction to get the actual alpha to place in the start vector.
  for (idx in nested_idx) {
    key <- as.character(idx)
    if (is.null(resolved[[key]])) next
    iec <- info$iec[[idx]]
    na <- length(info$extra[[idx]]$si$alpha)
    start[iec[1:na]] <- resolved[[key]] - .si_a0(info$extra[[idx]])
  }

  start
}

########################
# Top-level orchestrator: sequentially (one pass, no re-visits) resolve (i.e. find good initialitation for) 
# every effect made of one single index margin via a ~n_init-candidate search, scored by
# the exact model log-likelihood after fitting the remaining coefficients via family$initialize_bundle(). 
# Other nested effects (exp, mgks, si_nexp, positive single indices, two margins) use the construction-time initialisation, alpha, 
# provided by the init function of their margin bundles (R/bundle_trans_*.R).
# NOTE: code is looking for initialisation alpha* = alpha + a0, so the returned start vector has alpha = alpha* - a0 for each effect.
#       The helper function .si_a0() returns a0 = 0 for the parameters of non single index margins, so the returned start vector has 
#       alpha = alpha* for those effect types.
.backfit_si_initialize <- function(y, nobs, E, x, family, offset, weights, info,
                                    n_init = 1000, n_eigen = 10, oversample = 10, seed = 1) {

  lpi <- attr(x, "lpi")
  unscaled <- attr(E, "use.unscaled")

  is_special <- !sapply(info$type, function(.z) identical(.z[1], "stand"))
  nested_idx <- which(is_special)

  if (length(nested_idx) == 0) {
    return( family$initialize_bundle(y = y, nobs = nobs, E = E, x = x, family = family,
                                      offset = offset, jj = lpi, unscaled = unscaled, weights = weights) )
  }

  # effects made of one single index margin are searched, the others keep their initial alpha
  is_si <- sapply(nested_idx, function(idx) {
    mg <- info$extra[[idx]]$si$margin
    length(mg) == 1 && mg[[1]]$bundle_nam == "bundle_trans_si"
  })

  si_search <- nested_idx[is_si]
  si_fixed <- nested_idx[!is_si]

  resolved <- list()
  for (idx in si_fixed) {
    resolved[[as.character(idx)]] <- info$extra[[idx]]$si$alpha + .si_a0(info$extra[[idx]])
  }

  for (idx in si_search) {

    Xi <- info$extra[[idx]]$si$margin[[1]]$X
    a0 <- info$extra[[idx]]$si$margin[[1]]$a0
    alpha <- info$extra[[idx]]$si$alpha
    cands <- .si_search_candidates(Xi = Xi, alpha = alpha + a0, n_init = n_init, n_eigen = n_eigen, oversample = oversample, seed = seed)

    best_score <- -Inf
    best_v <- cands[1, ]

    for (ii in seq_len(nrow(cands))) {

      resolved_try <- resolved
      resolved_try[[as.character(idx)]] <- cands[ii, ]

      start_try <- .fit_reduced_si(y = y, nobs = nobs, x = x, E = E, lpi = lpi, info = info,
                                    family = family, offset = offset, weights = weights,
                                    unscaled = unscaled, nested_idx = nested_idx, resolved = resolved_try,
                                    reeval_idx = si_search)

      score <- family$ll(y, x, start_try, weights, family, offset = offset, deriv = 0)$l

      if ( is.finite(score) && score > best_score ) {
        best_score <- score
        best_v <- cands[ii, ]
      }
    }

    resolved[[as.character(idx)]] <- best_v

  }

  .fit_reduced_si(y = y, nobs = nobs, x = x, E = E, lpi = lpi, info = info, family = family,
                   offset = offset, weights = weights, unscaled = unscaled,
                   nested_idx = nested_idx, resolved = resolved, reeval_idx = si_search)

}
