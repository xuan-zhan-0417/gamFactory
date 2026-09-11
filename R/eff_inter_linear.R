#'
#' Build interaction_linear effect
#'
#' @param Xi matrix to be projected via single index vector(s). When margin 2
#'           is also nested (\code{s(si(x1), si(x2))}), \code{Xi} is the
#'           column-bound margin-1/margin-2 matrix carrying attributes
#'           \code{"na1"} and \code{"na2"}; otherwise (\code{s(si(x1), t)})
#'           \code{Xi} is the margin-1 matrix only and \code{t} must be
#'           supplied.
#' @param basis function which takes \code{z1} (and \code{z2}) as input and
#'                  returns the model matrix and its derivatives.
#' @param a0 A numeric vector of initialization shifts (length \code{ncol(Xi)}).
#' @param t Vector of the plain (non-nested) margin-2 covariate. Only used
#'          when \code{Xi} carries no \code{"na2"} attribute.
#' @name eff_inter_linear
#' @rdname eff_inter_linear
#' @export eff_inter_linear
#'
eff_inter_linear <- function(Xi, basis, a0 = NULL, t = NULL){

  force(Xi); force(basis); force(a0); force(t);

  na1 <- attr(Xi, "na1")
  na2 <- attr(Xi, "na2")
  nested_2 <- !is.null(na2)

  if( nested_2 ){
    if( is.null(na1) ) {
      stop("Attribute 'na1' must be defined in Xi when 'na2' is set.")
    }
    na <- na1 + na2
    if (ncol(Xi) != na) stop("ncol(Xi) must equal na1 + na2.")
    if (is.null(a0)) { a0 <- numeric(na) }
    if (length(a0) != na) stop("length(a0) must equal na1 + na2.")
    idx1 <- 1:na1
    idx2 <- (na1 + 1):na
  } else {
    na <- ncol(Xi)
    if (is.null(a0)) { a0 <- numeric(na) }
  }

  eval <- function(param, deriv = 0){

    n <- nrow(Xi)
    beta <- param[ -(1:na) ]

    if( nested_2 ){

      alpha_1 <- param[ idx1 ]
      alpha_2 <- param[ idx2 ]

      # inner indices z1, z2
      ax_1 <- drop( Xi[, idx1, drop = FALSE] %*% (alpha_1 + a0[idx1]) )
      ax_2 <- drop( Xi[, idx2, drop = FALSE] %*% (alpha_2 + a0[idx2]) )

      store <- basis$evalX(z1 = ax_1, z2 = ax_2, deriv = deriv)
      store$Xi <- Xi
      store$g  <- cbind(z1 = ax_1, z2 = ax_2)

      if( deriv >= 1 ){

        # dg_k / dalpha_k ; the cross blocks dg_1/dalpha_2, dg_2/dalpha_1 are zero
        store$g1 <- list(
          g1_1 = Xi[, idx1, drop = FALSE],
          g1_2 = Xi[, idx2, drop = FALSE]
        )

        # df / dz_k
        store$f1 <- list(
          f1_1 = drop( store$X1$dz1 %*% beta ),
          f1_2 = drop( store$X1$dz2 %*% beta )
        )

        if( deriv >= 2 ){

          # d2f / dz_j dz_k
          store$f2 <- list(
            f2_11 = drop( store$X2$dz1_z1 %*% beta ),
            f2_22 = drop( store$X2$dz2_z2 %*% beta ),
            f2_12 = drop( store$X2$dz1_z2 %*% beta )
          )

          # z_k is linear in alpha_k, so every second derivative of g vanishes.
          # Kept per component so that a non-linear inner map can drop in later.
          store$g2 <- list(
            g2_1 = matrix(0, nrow = n, ncol = na1 * (na1 + 1) / 2),
            g2_2 = matrix(0, nrow = n, ncol = na2 * (na2 + 1) / 2)
          )
        }
      }

    } else {

      alpha <- param[ 1:na ]

      # Project covariates on single index vector
      ax <- drop( Xi %*% (alpha + a0) )

      store <- basis$evalX(z1 = ax, z2 = t, deriv = deriv)
      store$Xi <- Xi

      if( deriv >= 1 ){
        store$f1 <- drop( store$X1 %*% beta )
        store$g1 <- Xi
        if( deriv >= 2 ){
          store$f2 <- drop( store$X2 %*% beta )
          if( deriv >= 3 ){
            store$f3 <- drop( store$X3 %*% beta )
          }
        }
      }
    }

    o <- eff_inter_linear(Xi = Xi, basis = basis, a0 = a0, t = t)
    o$f     <- drop( store$X0 %*% beta )
    o$param <- param
    o$a0    <- a0
    o$na    <- na
    if( nested_2 ){ o$na1 <- na1; o$na2 <- na2 }
    o$store <- store
    o$deriv <- deriv

    return( o )

  }

  out <- structure(list("eval" = eval),
                    class = if(nested_2) c("inter_dlinear", "nested") else c("si", "nested"))

  return( out )

}
