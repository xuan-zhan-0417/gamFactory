# Minimal per-margin view of an inter_dlinear effect, shaped so that the
# single-index penalty routines can consume it. They read exactly
#   o$na, o$param[1:o$na], o$a0, o$store$Xi
# and nothing else; if that contract changes upstream, this is the only
# place that needs fixing.
#' @noRd
.marg_si_shim <- function(o, ik){
  list(na    = length(ik),
       param = o$param[ik],
       a0    = o$a0[ik],
       store = list(Xi = o$store$Xi[, ik, drop = FALSE]))
}

#' @noRd
.marg_idx_inter_dlinear <- function(o){
  na1 <- o$na1; na2 <- o$na2
  if( is.null(na1) || is.null(na2) ){
    stop("o$na1 / o$na2 not set by eff_inter_dlinear().")
  }
  if( is.null(o$a0) || length(o$a0) != na1 + na2 ){
    stop("o$a0 must be a numeric vector of length na1 + na2 ",
         "(normalise it in eff_inter_dlinear()).")
  }
  list(1:na1, (na1 + 1):(na1 + na2))
}


######
# Penalty on the variance of BOTH inner single indices.
#
#   P = (var(z1) - v)^2 + (var(z2) - v)^2
#
# z_k depends only on alpha_k, so the penalty is separable: it is literally two
# single-index penalties side by side, and its Hessian is block diagonal with an
# exactly zero alpha_1 / alpha_2 cross block. Hence .pen_var_si() once per margin
# rather than a second copy of its algebra.
#
.pen_var_inter_dlinear <- function(o, v, deriv = 0){
  
  na  <- o$na
  idx <- .marg_idx_inter_dlinear(o)
  
  # build_family_nl() never reads pen$d3, and .eval_penalties() can be called
  # with deriv = 3 on the outer-derivative path, so cap it and skip the p^3 loop.
  p <- lapply(idx, function(ik){
    .pen_var_si(o = .marg_si_shim(o, ik), v = v, deriv = min(deriv, 2))
  })
  
  l0 <- p[[1]]$d0 + p[[2]]$d0
  
  l1 <- l2 <- NULL
  if( deriv ){
    l1 <- numeric(na)
    l1[ idx[[1]] ] <- drop(p[[1]]$d1)
    l1[ idx[[2]] ] <- drop(p[[2]]$d1)
    
    if( deriv > 1 ){
      # base R, not Matrix::bdiag: build_family_nl() does
      #   ret$lbb[zz, zz] - lamVar * pen$d2
      # on a dense matrix, and a sparse d2 would silently promote it.
      l2 <- matrix(0, na, na)
      l2[ idx[[1]], idx[[1]] ] <- p[[1]]$d2
      l2[ idx[[2]], idx[[2]] ] <- p[[2]]$d2
    }
  }
  
  return( list("d0" = l0, "d1" = l1, "d2" = l2, "d3" = NULL) )
}


# Derivatives of the penalty's Hessian w.r.t. the smoothing parameters.
# Same separability, so the same reuse.
.pen_var_inter_dlinear_outer <- function(o, v, DaDr){
  
  na  <- o$na
  idx <- .marg_idx_inter_dlinear(o)
  
  h <- lapply(idx, function(ik){
    .pen_var_si_outer(o    = .marg_si_shim(o, ik),
                      v    = v,
                      DaDr = DaDr[ik, , drop = FALSE])
  })
  
  lapply(seq_len(ncol(DaDr)), function(ii){
    H <- matrix(0, na, na)
    H[ idx[[1]], idx[[1]] ] <- h[[1]][[ii]]
    H[ idx[[2]], idx[[2]] ] <- h[[2]][[ii]]
    H
  })
}