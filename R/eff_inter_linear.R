#'
#' Build interaction_linear effect
#' 
#' @param Xi matrix to be projected via single index vector \code{alpha}.
#' @param basis function which takes \code{si = Xi\%*\%alpha} as input and returns model
#'                  matrix and its derivatives w.r.t. \code{si}.
#' @name eff_inter_linear
#' @rdname eff_inter_linear
#' @export eff_inter_linear
#'
eff_inter_linear <- function(Xi, t, basis, a0 = NULL){
  
  force(Xi); force(basis); force(a0); force(t);
  
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
    store <- basis$evalX(z1 = ax, z2 = t,deriv = deriv)
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
    
    # # # =================== check f1,f2,f3 ===================
    # n_check <- min(20, length(ax))
    # 
    # f1_results <- matrix(NA, nrow = n_check, ncol = 2, dimnames = list(NULL, c("EX", "FD")))
    # f2_results <- matrix(NA, nrow = n_check, ncol = 2, dimnames = list(NULL, c("EX", "FD")))
    # f3_results <- matrix(NA, nrow = n_check, ncol = 2, dimnames = list(NULL, c("EX", "FD")))
    # 
    # for (i in 1:n_check) {
    # 
    #   mock_obj <- list(
    #     d0 = function(zi) { drop(basis$evalX(z = zi, t = t[i], deriv = 0)$X0 %*% beta) },
    #     d1 = function(zi) { drop(basis$evalX(z = zi, t = t[i], deriv = 1)$X1 %*% beta) },
    #     d2 = function(zi) { drop(basis$evalX(z = zi, t = t[i], deriv = 2)$X2 %*% beta) },
    #     d3 = function(zi) { drop(basis$evalX(z = zi, t = t[i], deriv = 3)$X3 %*% beta) }
    #   )
    # 
    #   res <- check_deriv(obj = mock_obj, param = ax[i], ord = 1:deriv)
    #   
    #   if (deriv >= 1) f1_results[i, ] <- res$fd1
    #   if (deriv >= 2) f2_results[i, ] <- res$fd2
    #   if (deriv >= 3) f3_results[i, ] <- res$fd3
    # }
    # 
    # 
    # if (deriv >= 1) {
    #   print(f1_results)
    # }
    # if (deriv >= 2) {
    #   print(f2_results)
    # }
    # if (deriv >= 3) {
    #   print(f3_results)
    # }
    # # ==============================================================

    o <- eff_inter_linear(Xi = Xi, t = t, basis = basis, a0 = a0)
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









