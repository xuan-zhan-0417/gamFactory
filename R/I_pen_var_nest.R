######
# Penalty on the variance of the index of each margin (with inner parameters) of a nested effect
#
#   P = sum_k (var(z_k) - v)^2
#
# z_k depends only on the parameters of margin k, so the penalty is block-separable: each block is the
# generic penalty .pen_var_gen() on the g / g1 / g2 (g3) of that margin (for a si margin g2 = g3 = 0).
#
.pen_var_nest <- function(o, v, deriv = 0){

  ia <- lapply(o$margin[o$nest], "[[", "idx")

  # build_family_nl() never reads pen$d3
  p <- lapply(seq_along(ia), function(k) .pen_var_gen(o = .nest_margin_view(o, k), v = v, deriv = min(deriv, 2)))

  l1 <- l2 <- NULL
  if( deriv ){
    l1 <- numeric(o$na)
    for(k in seq_along(ia)){ l1[ ia[[k]] ] <- drop(p[[k]]$d1) }

    if( deriv > 1 ){
      l2 <- matrix(0, o$na, o$na)
      for(k in seq_along(ia)){ l2[ ia[[k]], ia[[k]] ] <- p[[k]]$d2 }
    }
  }

  return( list("d0" = sum(sapply(p, "[[", "d0")), "d1" = l1, "d2" = l2, "d3" = NULL) )
}

# Derivatives of the penalty's Hessian w.r.t. the smoothing parameters, block by block as above
.pen_var_nest_outer <- function(o, v, DaDr){

  ia <- lapply(o$margin[o$nest], "[[", "idx")

  h <- lapply(seq_along(ia), function(k){
    .pen_var_gen_outer(o = .nest_margin_view(o, k), v = v, DaDr = DaDr[ia[[k]], , drop = FALSE])
  })

  lapply(seq_len(ncol(DaDr)), function(ii){
    H <- matrix(0, o$na, o$na)
    for(k in seq_along(ia)){ H[ ia[[k]], ia[[k]] ] <- h[[k]][[ii]] }
    H
  })
}

# Margin k (k-th margin with parameters) of a nested effect seen on its own: its parameters and its index with
# derivatives, in the layout read by .pen_var_gen() / .pen_var_gen_outer()
.nest_margin_view <- function(o, k){
  ik <- o$margin[[ o$nest[k] ]]$idx
  list(na = length(ik), param = o$param[ik], deriv = o$deriv,
       store = list(g = o$store$g[ , o$nest[k]], g1 = o$store$g1[[k]], g2 = o$store$g2[[k]], g3 = o$store$g3[[k]]))
}
