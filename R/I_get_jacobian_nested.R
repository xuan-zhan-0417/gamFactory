#'
#' Jacobian of a nested effect w.r.t. its parameters c(alpha, beta), at the inner indices of data
#'
#' @noRd
get_jacobian.nested <- function(object, data, param){

  beta <- param[ -seq_along(object$xt$si$alpha) ]

  x_nest <- Predict.matrix.nested(object, data = data, get.xa = TRUE)
  store <- .nest_evalX(object$xt$basis, x_nest$xa, deriv = 1)

  # df/dalpha_k = (dX/dz_k %*% beta) * dz_k/dalpha_k for each margin k with parameters ; df/dbeta = X0
  JJ <- cbind(do.call("cbind", Map(function(X1k, g1k) drop(X1k %*% beta) * g1k, store$X1, x_nest$xa_da)),
              store$X0)

  return( list("JJ" = JJ, "xa" = x_nest$xa) )
}
