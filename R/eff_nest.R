#'
#' Build a nested effect
#'
#' @param margin list of the margin descriptions \code{si$margin} (one or two) built by
#'               \link{smooth.construct.nest.smooth.spec}: bundle name, positions \code{idx} of the margin's
#'               parameters in the inner parameter vector, data and centring constant.
#' @param basis outer basis, whose \code{evalX(z1, [z2,] deriv)} returns the model matrix and its
#'              derivatives w.r.t. the indices of the margins with inner parameters.
#' @details The parameter vector is \code{c(alpha, beta)}: \code{alpha} stacks the inner parameters of
#'          the margins, \code{beta} the outer spline coefficients. The index \code{z_k} of each margin and its
#'          derivatives are computed by the \code{eval} function of its bundle (\link{trans_bundles}).
#' @name eff_nest
#' @rdname eff_nest
#' @export eff_nest
#'
eff_nest <- function(margin, basis){

  force(margin); force(basis)

  na   <- length(unlist(lapply(margin, "[[", "idx")))
  nest <- .nest_with_par(margin)                       # margins with inner parameters
  bdl  <- lapply(margin, function(mk) .nest_bundle(mk$bundle_nam))

  eval <- function(param, deriv = 0){

    beta <- param[ -(1:na) ]

    m <- lapply(seq_along(margin), function(k) bdl[[k]]$eval(margin[[k]], param[margin[[k]]$idx], deriv))

    store <- .nest_evalX(basis, .nest_z(m), deriv)
    store$g <- .nest_z(m)

    # one element per nested margin k (same order as store$X1): derivatives of z_k w.r.t. its own
    # parameters (g1, g2, g3) and of f w.r.t. z_k (f1) and (z_k, z_j) (f2)
    if( deriv >= 1 ){
      store$g1 <- lapply(m[nest], "[[", "g1")
      store$f1 <- lapply(store$X1, function(X1) drop(X1 %*% beta))
      if( deriv >= 2 ){
        store$g2 <- lapply(m[nest], "[[", "g2")
        store$f2 <- lapply(nest, function(k) lapply(nest, function(j) drop(store$X2[[paste0("dz", k, "_z", j)]] %*% beta)))
        if( deriv >= 3 ){
          store$g3 <- lapply(m[nest], "[[", "g3")
        }
      }
    }

    o <- eff_nest(margin = margin, basis = basis)
    o$f     <- drop( store$X0 %*% beta )
    o$param <- param
    o$na    <- na
    o$store <- store
    o$deriv <- deriv

    return( o )

  }

  out <- structure(list("eval" = eval, "margin" = margin, "nest" = nest), class = c("nest", "nested"))

  return( out )

}
