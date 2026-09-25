######
# Penalty on the variance of both inner indices of an inter_ee effect.
#
#   P = (var(z1) - v)^2 + (var(z2) - v)^2
#
# Both z_k = exp(scale_k) * (g_k - mean(g_k)) are nonlinear in their own
# (scale_k, rate coefficients) block, so each slice uses the generic
# g / g1 / g2 treatment of .pen_var_gen(). The two slices act on disjoint
# parameter blocks, hence the penalty and its Hessian are block-separable.
#
#' @noRd
.pen_var_inter_ee <- function(o, v, deriv = 0){

  d1 <- o$d1
  na <- o$na
  idx <- list(1:d1, (d1 + 1):na)

  p <- lapply(1:2, function(k){
    shim <- list(na = length(idx[[k]]), deriv = o$deriv,
                 param = o$param[idx[[k]]],
                 store = list(g  = o$store$g[ , k],
                              g1 = o$store$g1[[k]],
                              g2 = if( !is.null(o$store$g2) ) o$store$g2[[k]] else NULL))
    .pen_var_gen(o = shim, v = v, deriv = min(deriv, 2))
  })

  l0 <- p[[1]]$d0 + p[[2]]$d0

  l1 <- l2 <- NULL
  if( deriv ){
    l1 <- numeric(na)
    l1[idx[[1]]] <- drop(p[[1]]$d1)
    l1[idx[[2]]] <- drop(p[[2]]$d1)

    if( deriv > 1 ){
      l2 <- matrix(0, na, na)
      l2[idx[[1]], idx[[1]]] <- p[[1]]$d2
      l2[idx[[2]], idx[[2]]] <- p[[2]]$d2
    }
  }

  return( list("d0" = l0, "d1" = l1, "d2" = l2, "d3" = NULL) )
}
