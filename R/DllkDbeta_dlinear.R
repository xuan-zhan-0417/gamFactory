#'
#' Derivatives of nested interaction (double single-index) effects
#' 
#' @rdname DllkDbeta.inter_dlinear
#' @export DllkDbeta.inter_dlinear
#' @export
#'
DllkDbeta.inter_dlinear <- function(o, llk, deriv = 1, param = NULL){
  
  if( deriv == 0 ){ return( list() ) }
  if( deriv >  2 ){
    stop("DllkDbeta.inter_dlinear does not implement deriv > 2. ",
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
  
  Xi   <- o$store$Xi
  Xi_1 <- Xi[, 1:na1, drop = FALSE]
  Xi_2 <- Xi[, (na1 + 1):na, drop = FALSE]
  
  X <- o$store$X0                      # outer design matrix
  
  der1 <- der2 <- der3 <- NULL
  
  # ---------------------------------------------------------------
  # 1. Gradient
  # ---------------------------------------------------------------
  le   <- llk$d1
  f1_1 <- o$store$f1$f1_1              # df / dz1
  f1_2 <- o$store$f1$f1_2              # df / dz2
  
  # dz_k / dalpha_k = Xi_k ; the cross blocks are zero, hence the block form
  der1 <- c( crossprod(Xi_1, le * f1_1),
             crossprod(Xi_2, le * f1_2),
             crossprod(X, le) )
  
  # ---------------------------------------------------------------
  # 2. Hessian
  # ---------------------------------------------------------------
  if( deriv > 1 ){
    
    lee   <- llk$d2
    f2_11 <- o$store$f2$f2_11
    f2_22 <- o$store$f2$f2_22
    f2_12 <- o$store$f2$f2_12
    
    X1_dz1 <- o$store$X1$dz1           # dX / dz1
    X1_dz2 <- o$store$X1$dz2           # dX / dz2
    
    # --- alpha-alpha -------------------------------------------------
    # NOTE: the  sum_i le_i * f1_i * d2(z)/dalpha2  term is absent because
    # z_k = Xi_k %*% alpha_k is linear, i.e. store$g2 == 0.  Reinstate it
    # here if a non-linear inner map (e.g. positive_si) is ever added.
    lgg_11 <- le * f2_11 + lee * f1_1^2
    lgg_22 <- le * f2_22 + lee * f1_2^2
    lgg_12 <- le * f2_12 + lee * f1_1 * f1_2
    
    ll_aa11 <- crossprod(Xi_1, lgg_11 * Xi_1)
    ll_aa22 <- crossprod(Xi_2, lgg_22 * Xi_2)
    ll_aa12 <- crossprod(Xi_1, lgg_12 * Xi_2)
    
    ll_aa <- rbind(cbind(ll_aa11,      ll_aa12),
                   cbind(t(ll_aa12),   ll_aa22))
    
    # --- beta-beta ---------------------------------------------------
    ll_bb <- crossprod(X, lee * X)
    
    # --- beta-alpha  (nb x na) ---------------------------------------
    # d2 l / dbeta dalpha_k = X' diag(lee * f1_k) Xi_k + (dX/dz_k)' diag(le) Xi_k
    ll_ba <- cbind(
      crossprod(X, (lee * f1_1) * Xi_1) + crossprod(X1_dz1, le * Xi_1),
      crossprod(X, (lee * f1_2) * Xi_2) + crossprod(X1_dz2, le * Xi_2)
    )
    
    der2 <- rbind(cbind(ll_aa, t(ll_ba)),
                  cbind(ll_ba, ll_bb))
  }
  
  out <- list("d1" = der1, "d2" = der2, "d3" = der3)
  
  return( out )
}