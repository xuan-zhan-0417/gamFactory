#'
#' Derivatives of the si(x) / exp(x) interaction effect
#'
#' @rdname DllkDbeta.inter_le
#' @export DllkDbeta.inter_le
#' @export
#'
DllkDbeta.inter_le <- function(o, llk, deriv = 1, param = NULL){

  if( deriv == 0 ){ return( list() ) }
  if( deriv >  2 ){
    stop("DllkDbeta.inter_le does not implement deriv > 2. ",
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

  na  <- o$na
  na1 <- o$na1
  na2 <- o$na2
  p2  <- 1 + na2                      # margin-2 block size: alpha_scale + alpha_2

  g1_1 <- o$store$g1$g1_1             # n x na1  (= dz1/dalpha_1, linear -> Xi_1)
  g1_2 <- o$store$g1$g1_2             # n x p2   (= dz2/d(alpha_scale, alpha_2))

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

    g2_2 <- o$store$g2$g2_2           # n x p2*(p2+1)/2, d2(z2)/d(alpha_scale,alpha_2)^2

    X1_dz1 <- o$store$X1$dz1          # dX / dz1
    X1_dz2 <- o$store$X1$dz2          # dX / dz2

    # --- alpha1 - alpha1 (z1 linear in alpha_1: its own 2nd derivative is 0) ---
    lgg_11  <- le * f2_11 + lee * f1_1^2
    ll_aa11 <- crossprod(g1_1, lgg_11 * g1_1)

    # --- margin-2 own block (z2 nonlinear in alpha_scale/alpha_2) -------------
    # Same construction as the single-margin generic (DllkDbeta.nested): the
    # g2_2-weighted term accounts for z2's own curvature, on top of the
    # "linear-in-g1" term shared with the alpha1-alpha1 block above.
    lgg_22 <- le * f2_22 + lee * f1_2^2
    lg_2   <- le * f1_2

    tmp <- colSums(g2_2 * lg_2)
    ll_aa22 <- matrix(0, p2, p2)
    ll_aa22[lower.tri(ll_aa22, diag = TRUE)] <- tmp
    ll_aa22 <- t(ll_aa22)
    ll_aa22[lower.tri(ll_aa22, diag = TRUE)] <- tmp
    ll_aa22 <- ll_aa22 + crossprod(g1_2, lgg_22 * g1_2)

    # --- alpha1 / margin-2 cross block (only through f2_12, z1 and z2 depend
    #     on disjoint parameters so there is no cross "g2" term) -------------
    lgg_12  <- le * f2_12 + lee * f1_1 * f1_2
    ll_aa12 <- crossprod(g1_1, lgg_12 * g1_2)

    ll_aa <- rbind(cbind(ll_aa11,    ll_aa12),
                   cbind(t(ll_aa12), ll_aa22))

    # --- beta-beta -------------------------------------------------------
    ll_bb <- crossprod(X, lee * X)

    # --- beta-alpha (nb x na) ---------------------------------------------
    ll_ba <- cbind(
      crossprod(X, (lee * f1_1) * g1_1) + crossprod(X1_dz1, le * g1_1),
      crossprod(X, (lee * f1_2) * g1_2) + crossprod(X1_dz2, le * g1_2)
    )

    der2 <- rbind(cbind(ll_aa, t(ll_ba)),
                  cbind(ll_ba, ll_bb))
  }

  out <- list("d1" = der1, "d2" = der2, "d3" = der3)

  return( out )
}
