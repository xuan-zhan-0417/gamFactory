#'
#' Derivatives of the exp(x1) / exp(x2) interaction effect
#'
#' @rdname DllkDbeta.inter_ee
#' @export DllkDbeta.inter_ee
#' @export
#'
DllkDbeta.inter_ee <- function(o, llk, deriv = 1, param = NULL){

  if( deriv == 0 ){ return( list() ) }
  if( deriv >  2 ){
    stop("DllkDbeta.inter_ee does not implement deriv > 2. ",
         "Third order derivatives are only needed by DHessDrho (outDer = TRUE).")
  }

  if( is.null(param) ){
    param <- o$param
    if( is.null(param) ){ stop("param vector not provided!") }
  }

  # Need to update the object
  if( is.null(o$param) || !identical(param, o$param) || o$deriv < deriv ){
    o <- o$eval(param = param, deriv = deriv)
  }

  d1 <- o$d1                          # size of the margin-1 block: alpha_scale_1 + rate coefficients
  d2 <- o$d2

  g1_1 <- o$store$g1$g1_1             # n x d1   dz1 / d(margin-1 params)
  g1_2 <- o$store$g1$g1_2             # n x d2   dz2 / d(margin-2 params)

  X <- o$store$X0                     # outer design matrix

  der1 <- der2 <- der3 <- NULL

  # ---------------------------------------------------------------
  # 1. Gradient
  # ---------------------------------------------------------------
  le   <- llk$d1
  f1_1 <- o$store$f1$f1_1             # df / dz1
  f1_2 <- o$store$f1$f1_2             # df / dz2

  der1 <- c( crossprod(g1_1, le * f1_1),
             crossprod(g1_2, le * f1_2),
             crossprod(X, le) )

  # ---------------------------------------------------------------
  # 2. Hessian
  # ---------------------------------------------------------------
  if( deriv > 1 ){

    lee   <- llk$d2
    f2_11 <- o$store$f2$f2_11
    f2_22 <- o$store$f2$f2_22
    f2_12 <- o$store$f2$f2_12

    X1_dz1 <- o$store$X1$dz1
    X1_dz2 <- o$store$X1$dz2

    # Own block of a nonlinear margin: (g2-weighted curvature term) + t(g1) %*% (lgg * g1),
    # same construction as the single-margin generic DllkDbeta.nested.
    own_block <- function(g1, g2, lg, lgg, p){
      tmp <- colSums(g2 * lg)
      M <- matrix(0, p, p)
      M[lower.tri(M, diag = TRUE)] <- tmp
      M <- t(M)
      M[lower.tri(M, diag = TRUE)] <- tmp
      M + crossprod(g1, lgg * g1)
    }

    ll_aa11 <- own_block(g1_1, o$store$g2$g2_1, le * f1_1, le * f2_11 + lee * f1_1^2, d1)
    ll_aa22 <- own_block(g1_2, o$store$g2$g2_2, le * f1_2, le * f2_22 + lee * f1_2^2, d2)

    # cross block: z1 and z2 depend on disjoint parameters, only f2_12 couples them
    ll_aa12 <- crossprod(g1_1, (le * f2_12 + lee * f1_1 * f1_2) * g1_2)

    ll_aa <- rbind(cbind(ll_aa11,    ll_aa12),
                   cbind(t(ll_aa12), ll_aa22))

    ll_bb <- crossprod(X, lee * X)

    ll_ba <- cbind(
      crossprod(X, (lee * f1_1) * g1_1) + crossprod(X1_dz1, le * g1_1),
      crossprod(X, (lee * f1_2) * g1_2) + crossprod(X1_dz2, le * g1_2)
    )

    der2 <- rbind(cbind(ll_aa, t(ll_ba)),
                  cbind(ll_ba, ll_bb))
  }

  return( list("d1" = der1, "d2" = der2, "d3" = der3) )
}
