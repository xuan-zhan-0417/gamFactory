#'
#' Build interaction effect
#' 
#' @param Xi matrix to be projected via single index vector \code{alpha}.
#' @param basis function which takes \code{si = Xi\%*\%alpha} as input and returns model
#'                  matrix and its derivatives w.r.t. \code{si}.
#' @name eff_inter
#' @rdname eff_inter
#' @export eff_inter
#'
eff_inter <- function(Xi, t, basis, a0 = NULL){
  
  force(Xi); force(basis); force(a0)
  
  eval <- function(param, deriv = 0){

    na <- ncol( Xi )
    nb <- length(param) - na 
    
    # Single index and spline coefficients
    alpha <- param[ 1:na ]
    beta <- param[ -(1:na) ]
    
    if( is.null(a0) ){
      a0 <- alpha * 0
    }
    
    # Project covariates on single index vector 
    ax <- drop( Xi %*% (alpha + a0) )
    
    # Build P-spline basis and its derivatives
    # The error is probably due to the fact that no observations falls within range
    store <- basis$evalX(z = ax, t = t,deriv = deriv)
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
    
    # # =================== check f1,f2,f3 ===========================================
    # 
    # chk <- gamFactory:::check_wrap_z_deriv_numDeriv_mean(
    #   Xi = Xi,
    #   t = t,
    #   basis = basis,
    #   param = param,
    #   a0 = a0,
    #   deriv = 3,
    #   method = "Richardson"
    # )
    # 
    # chk
    # # ==============================================================

    o <- eff_inter(Xi = Xi, t = t, basis = basis, a0 = a0)
    o$f <- drop( store$X0 %*% beta )
    o$param <- param
    o$a0 <- a0
    o$na <- na
    o$store <- store
    o$deriv <- deriv
    
    return( o )
    
  }
  
  out <- structure(list("eval" = eval), class = c("si", "nested"))
  
  return( out )
  
}









