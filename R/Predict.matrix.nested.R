#'
#' Predict using nested effects
#'
#' @param object smooth object of class "nest" (see \link{smooth.construct.nest.smooth.spec}).
#' @param data data frame / list holding the terms.
#' @param get.xa if TRUE, return the inner indices \code{xa} (n x number of margins) and, for each margin with
#'               inner parameters, the derivatives of its index w.r.t. them (\code{xa_da}), instead of the model matrix.
#' @name Predict.matrix.nested
#' @rdname Predict.matrix.nested
#' @export
#'
Predict.matrix.nested <- function(object, data, get.xa = FALSE){

  si <- object$xt$si

  # each margin at the new data, centred with the constant frozen at its fitted value
  m <- lapply(seq_along(si$margin), function(k){
    b  <- .nest_bundle(si$margin[[k]]$bundle_nam)
    mk <- b$newdata(si$margin[[k]], data[[ object$term[k] ]])
    b$eval(mk, si$alpha[mk$idx], deriv = as.numeric(get.xa), xm = mk$xm)
  })
  z <- .nest_z(m)

  if( get.xa ){
    return( list(xa = z, xa_da = lapply(m[ .nest_with_par(si$margin) ], "[[", "g1")) )
  }

  X0 <- .nest_evalX(object$xt$basis, z, deriv = 0)$X0

  # [ 0_alpha , X0 ] : predict.gam multiplies the first block by alpha, which has no effect
  cbind(matrix(0, nrow(X0), length(si$alpha)), X0)
}
