#' Predict using 2D inter linear effects (ANOVA Tensor Product)
#' 
#' @param object 
#' @param data 
#' @param get.xa 
#' @noRd
.predict.matrix.inter_linear <- function(object, data, get.xa = FALSE){
  
  term_x <- object$term[1]  # e.g "X_mat"
  term_t <- object$term[2]  # e.g "t_doy"
  
  si <- object$xt$si
  alpha <- si$alpha
  a0 <- si$a0
  di <- length(alpha)
  
  Xi <- data[[term_x]]
  # if(!is.matrix(Xi)){
  #   Xi <- matrix(Xi, ncol = length(alpha))
  # }

  Xi <- t(t(Xi) - si$xm) %*% si$B
  xa <- Xi %*% (alpha + a0)  # z
  
  if(get.xa){ 
    return(list(xa = xa, xa_da = Xi))
  }
  
  t_var <- data[[term_t]]  - si$tm  
  
  X0 <- object$xt$basis$evalX(z1 = xa, z2 = t_var, deriv = 0)$X0
  
  # [ 0_di, X_2D ]
  Xtot <- cbind(matrix(0, nrow(X0), di), X0) 
  
  attr(Xtot, "inner_linpred_unscaled") <- xa
  
  return(Xtot)
}