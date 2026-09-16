######
# Penalty on the variance of both inner indices of an inter_le effect.
#
#   P = (var(z1) - v)^2 + (var(z2) - v)^2
#
# z1 = Xi_1 %*% (alpha_1 + a0_1) is linear in alpha_1, so its slice is handled
# by .pen_var_si() exactly as for a plain "si" margin.
# z2 = exp(alpha_scale) * (g(alpha_2) - mean(g(alpha_2))) is nonlinear in its
# own (alpha_scale, alpha_2) block (via expsmooth), so its slice needs the
# generic g/g1/g2 treatment of .pen_var_gen() -- the same machinery used for
# a standalone nexpsm effect.
#
# The two slices act on disjoint parameter blocks, so the penalty (and its
# Hessian) is block-separable, as in .pen_var_inter_dlinear().
#
#' @noRd
.pen_var_inter_le <- function(o, v, deriv = 0){

  na1 <- o$na1
  na2 <- o$na2
  na  <- o$na
  idx1 <- 1:na1
  idx2 <- (na1 + 1):na   # alpha_scale + alpha_2

  shim1 <- list(na = na1, param = o$param[idx1], a0 = o$a0_1,
                store = list(Xi = o$store$Xi_1))
  p1 <- .pen_var_si(o = shim1, v = v, deriv = min(deriv, 2))

  shim2 <- list(na = length(idx2), deriv = o$deriv,
                param = o$param[idx2],
                store = list(g  = o$store$g[ , "z2"],
                             g1 = o$store$g1$g1_2,
                             g2 = if( !is.null(o$store$g2) ) o$store$g2$g2_2 else NULL))
  p2 <- .pen_var_gen(o = shim2, v = v, deriv = min(deriv, 2))

  l0 <- p1$d0 + p2$d0

  l1 <- l2 <- NULL
  if( deriv ){
    l1 <- numeric(na)
    l1[idx1] <- drop(p1$d1)
    l1[idx2] <- drop(p2$d1)

    if( deriv > 1 ){
      l2 <- matrix(0, na, na)
      l2[idx1, idx1] <- p1$d2
      l2[idx2, idx2] <- p2$d2
    }
  }

  return( list("d0" = l0, "d1" = l1, "d2" = l2, "d3" = NULL) )
}
