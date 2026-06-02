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
    # browser()
    # check_wrap_z_deriv_numDeriv_mean <- function(Xi, t, basis, param, a0 = NULL,
    #                                              deriv = 3,
    #                                              method = "Richardson",
    #                                              method.args = list()) {
    #   
    #   na <- ncol(Xi)
    #   alpha <- param[1:na]
    #   beta  <- param[-(1:na)]
    #   
    #   if (is.null(a0)) {
    #     a0 <- alpha * 0
    #   }
    #   
    #   z0 <- drop(Xi %*% (alpha + a0))
    #   
    #   # analytical derivatives
    #   st <- basis$evalX(z = z0, t = t, deriv = deriv)
    #   
    #   # scalar function: common shift c applied to all z_i
    #   eta_mean_at <- function(c) {
    #     z_eval <- z0 + c
    #     st_eval <- basis$evalX(z = z_eval, t = t, deriv = 0)
    #     eta_eval <- drop(st_eval$X0 %*% beta)
    #     mean(eta_eval)
    #   }
    #   
    #   out <- list()
    #   
    #   if (deriv >= 1) {
    #     f1_ex <- mean(drop(st$X1 %*% beta))
    #     
    #     f1_fd <- numDeriv::grad(
    #       func = eta_mean_at,
    #       x = 0,
    #       method = method,
    #       method.args = method.args
    #     )
    #     
    #     out$f1 <- c(
    #       EX = f1_ex,
    #       FD = f1_fd,
    #       abs_err = abs(f1_ex - f1_fd),
    #       rel_err = abs(f1_ex - f1_fd) / max(1, abs(f1_fd))
    #     )
    #   }
    #   
    #   if (deriv >= 2) {
    #     f2_ex <- mean(drop(st$X2 %*% beta))
    #     
    #     f2_fd <- as.numeric(numDeriv::hessian(
    #       func = eta_mean_at,
    #       x = 0,
    #       method = method,
    #       method.args = method.args
    #     ))
    #     
    #     out$f2 <- c(
    #       EX = f2_ex,
    #       FD = f2_fd,
    #       abs_err = abs(f2_ex - f2_fd),
    #       rel_err = abs(f2_ex - f2_fd) / max(1, abs(f2_fd))
    #     )
    #   }
    #   
    #   if (deriv >= 3) {
    #     
    #     d2_mean_at <- function(c) {
    #       as.numeric(numDeriv::hessian(
    #         func = eta_mean_at,
    #         x = c,
    #         method = method,
    #         method.args = method.args
    #       ))
    #     }
    #     
    #     f3_ex <- mean(drop(st$X3 %*% beta))
    #     
    #     f3_fd <- numDeriv::grad(
    #       func = d2_mean_at,
    #       x = 0,
    #       method = method,
    #       method.args = method.args
    #     )
    #     
    #     out$f3 <- c(
    #       EX = f3_ex,
    #       FD = f3_fd,
    #       abs_err = abs(f3_ex - f3_fd),
    #       rel_err = abs(f3_ex - f3_fd) / max(1, abs(f3_fd))
    #     )
    #   }
    #   
    #   out
    # }
    # 
    # chk <- check_wrap_z_deriv_numDeriv_mean(
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









