# Helpers shared by the margin bundles (R/bundle_trans_*.R) and by the generic code of nested effects.
#   .nest_bundle        the bundle of a margin, rebuilt from its stored name (as a family from bundle_nam)
#   .nest_eval_centred  z = g - xm and its derivatives, g being the output of expsmooth() / mgks() / ...
#   .nest_eval_scaled   z = exp(scale) (g - xm) and its derivatives
#   .nest_with_par      which margins have inner parameters (all but plain ones)
#   .nest_z             the indices of all margins as an n x K matrix
#   .nest_evalX         the outer basis (1D or 2D) evaluated at those indices
#   .sd_n               standard deviation dividing by n, as the penalty on the variance of the index does
#   .nest_kex_default   range of the outer knots for margins without a range of their own (plain, si_nexp)

.nest_bundle <- function(bundle_nam) do.call(bundle_nam, list())

# xm = NULL (fitting): centre g with its current mean, which moves with the parameters, so the derivatives are
# centred too. Otherwise (prediction) centre with the frozen constant xm.
.nest_eval_centred <- function(inner, deriv, xm){
  g <- drop(inner$d0)
  centre <- is.null(xm)
  if( centre ){ xm <- mean(g) }
  ctr <- function(D) if( centre ) t(t(D) - colMeans(D)) else D
  out <- list(z = g - xm, xm = xm)
  if( deriv >= 1 ){ out$g1 <- ctr(inner$d1) }
  if( deriv >= 2 ){ out$g2 <- ctr(inner$d2) }
  if( deriv >= 3 ){ out$g3 <- ctr(inner$d3) }
  out
}

# The first parameter is the log-scale. Every derivative involving it collapses onto a lower order one
# (dz/dscale = z, d2z/dscale dpar_j = dz/dpar_j, ...), so g1 starts with z, g2 with g1 and g3 with g2.
.nest_eval_scaled <- function(inner, scale, deriv, xm){
  a <- exp(scale)
  o <- .nest_eval_centred(inner, deriv, xm)
  out <- list(z = a * o$z, xm = o$xm)
  if( deriv >= 1 ){ out$g1 <- cbind(out$z, a * o$g1) }
  if( deriv >= 2 ){ out$g2 <- cbind(out$g1, a * o$g2) }
  if( deriv >= 3 ){ out$g3 <- cbind(out$g2, a * o$g3) }
  out
}

.nest_with_par <- function(margin) which(lengths(lapply(margin, "[[", "idx")) > 0)

# m: outputs of the bundles' eval, one per margin
.nest_z <- function(m){
  z <- do.call("cbind", lapply(m, "[[", "z"))
  colnames(z) <- paste0("z", seq_along(m))
  z
}

# basis$evalX(z1, deriv) for one margin, basis$evalX(z1, z2, deriv) for two
.nest_evalX <- function(basis, z, deriv){
  arg <- lapply(seq_len(ncol(z)), function(k) z[ , k])
  names(arg) <- colnames(z)
  do.call(basis$evalX, c(arg, list(deriv = deriv)))
}

.sd_n <- function(x) {
  n <- length(x)
  sd(x) * sqrt(n - 1) / sqrt(n)
}

# the range of z widened by 1, and at least (-6, 6)
.nest_kex_default <- function(z) {
  r <- range(z, na.rm = TRUE) + c(-1, 1)
  if( r[1] < -6 || r[2] > 6 ) r else c(-6, 6)
}
