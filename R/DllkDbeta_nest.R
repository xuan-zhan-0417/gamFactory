#'
#' Derivatives of nested effects
#'
#' @rdname DllkDbeta.nest
#' @export DllkDbeta.nest
#' @export
#'
DllkDbeta.nest <- function(o, llk, deriv = 1, param = NULL){

  if( deriv == 0 ){ return( list() ) }

  if( is.null(param) ){ param <- o$param }

  # Need to update the object
  if( is.null(o$param) || !identical(param, o$param) || o$deriv < deriv ){
    o <- o$eval(param = param, deriv = deriv)
  }

  s  <- o$store
  le <- llk$d1
  J  <- .jac_nest(o)                                  # d eta / d(alpha, beta)

  d1 <- drop( crossprod(J, le) )
  d2 <- NULL

  if( deriv > 1 ){

    # Gauss-Newton part, plus the terms where le multiplies the second derivatives of eta
    d2 <- crossprod(J, llk$d2 * J)

    ia <- lapply(o$margin[o$nest], "[[", "idx")        # alpha block of each nested margin
    ib <- (o$na + 1):length(param)                     # beta

    for(k in seq_along(ia)){
      # z_k nonlinear in its own parameters (zero for a linear margin)
      d2[ia[[k]], ia[[k]]] <- d2[ia[[k]], ia[[k]]] +
                              .vec_to_sym_mat(colSums(s$g2[[k]] * (le * s$f1[[k]])), length(ia[[k]]))
      # d2 eta / dbeta dalpha_k = (dX/dz_k) dz_k/dalpha_k
      H_bk <- crossprod(s$X1[[k]], le * s$g1[[k]])
      d2[ib, ia[[k]]] <- d2[ib, ia[[k]]] + H_bk
      d2[ia[[k]], ib] <- d2[ia[[k]], ib] + t(H_bk)
      # f nonlinear in (z_k, z_j)
      for(j in seq_along(ia)){
        d2[ia[[k]], ia[[j]]] <- d2[ia[[k]], ia[[j]]] + crossprod(s$g1[[k]], (le * s$f2[[k]][[j]]) * s$g1[[j]])
      }
    }
  }

  return( list("d1" = d1, "d2" = d2, "d3" = NULL) )
}
