#'
#' Compute jacobian of nested effects
#' 
#' @noRd
get_jacobian.nested <- function(object,data, param){
  
  if(class(object)[1] == "si"){
    return( .get.jacobian.si(object, data, param) ) 
  }
  if(class(object)[1] == "nexpsm" || class(object)[1] == "mgks"){
    return( .get.jacobian.nexpsm.mgks(object, data, param) )
  }
  if(class(object)[1] == "si_nexpsm"){ 
    return(.get.jacobian.si_nexpsm(object, data, param))
  }
  if(class(object)[1] == "inter_linear"){
    return(.get.jacobian.inter_linear(object, data, param))
  }
  if(class(object)[1] == "inter_nexp"){ 
    return(.get.jacobian.inter_nexp(object, data, param))
  }
  if(class(object)[1] == "inter_mgks"){ 
    return(.get.jacobian.inter_mgks(object, data, param))
  }
  if(class(object)[1] == "inter_dlinear"){ 
    return(.get.jacobian.inter_dlinear(object, data, param))
  }
  stop("I do not know this effect type")
  
}


# ------------------------------subfunction------------------------------
.get.jacobian.si <- function(object, data, param){
  
  na <- length(object$xt$si$alpha)
  # Single index and spline coefficients
  beta <- param[ -(1:na) ]
  
  x_nest <- Predict.matrix.nested(object, data = data, get.xa = TRUE)
  store <- object$xt$basis$evalX(x = x_nest$xa, deriv = 1)
  
  JJ <- cbind(drop(store$X1 %*% beta) * x_nest$xa_da, # df/da = M1%*%b * ds/da
              store$X0) # df/db = Ma
  return(list("JJ" = JJ, "xa" = NULL) )
}

.get.jacobian.nexpsm.mgks <- function(object, data, param){
  
  na <- length(object$xt$si$alpha)
  # Single index and spline coefficients
  beta <- param[ -(1:na) ]
  
  x_nest <- Predict.matrix.nested(object, data = data, get.xa = TRUE)
  
  store <- object$xt$basis$evalX(x = x_nest$xa, deriv = 1)
  X1beta <- drop(store$X1 %*% beta)
  
  JJ <- cbind(X1beta * x_nest$xa, # df/da = M1%*%b * ds/da0 (where ds/da0 = s because s = exp(a0) xa)
              X1beta * x_nest$xa_da, # df/da = M1%*%b * ds/da
              store$X0) # df/db = Ma 
  
  return(list("JJ" = JJ, "xa" = NULL) )
}

.get.jacobian.si_nexpsm <- function(object, data, param){
  na <- length(object$xt$si$alpha)
  beta <- param[ -(1:na) ]
  
  x_nest <- gamFactory:::Predict.matrix.nested(object, data = data, get.xa = TRUE)
  store <- object$xt$basis$evalX(x = x_nest$xa, deriv = 1)
  X1beta <- drop(store$X1 %*% beta)  # derivative to beta
  
  # ∂eta/∂alpha_si = (dM/dz %*% beta) * ∂z/∂alpha_si
  J_alpha_si <- X1beta * x_nest$xa_dalpha_si  # n × na_si
  # ∂eta/∂alpha_nexp = (dM/dz %*% beta) * ∂z/∂alpha_nexp
  J_alpha_nexp <- X1beta * x_nest$xa_dalpha_nexp  # n × na_nexp
  # ∂eta/∂beta = M
  J_beta <- store$X0
  
  JJ <- cbind(J_alpha_nexp, J_alpha_si, J_beta)
  
  return( list("JJ" = JJ, "xa" = x_nest$xa) )
}

.get.jacobian.inter_linear <- function(object, data, param){
  
  na <- length(object$xt$si$alpha)
  # nested index and spline coefficients
  beta <- param[ -(1:na) ]
  
  x_nest <- Predict.matrix.nested(object, data = data, get.xa = TRUE)
  store <- object$xt$basis$evalX(z1 = x_nest$xa, z2 = object$xt$si$t , deriv = 1)
  
  JJ <- cbind(drop(store$X1 %*% beta) * x_nest$xa_da, # df/da = M1%*%b * ds/da
              store$X0) # df/db = Ma
  return(list("JJ" = JJ, "xa" = x_nest$xa) )
}

.get.jacobian.inter_nexp <- function(object, data, param){
  
  na <- length(object$xt$si$alpha)
  # Single index and spline coefficients
  beta <- param[ -(1:na) ]
  
  x_nest <- Predict.matrix.nested(object, data = data, get.xa = TRUE)
  
  store <- object$xt$basis$evalX(z = x_nest$xa, t = object$xt$si$t, deriv = 1)
  X1beta <- drop(store$X1 %*% beta)
  
  JJ <- cbind(X1beta * x_nest$xa, # df/da = M1%*%b * ds/da0 (where ds/da0 = s because s = exp(a0) xa)
              X1beta * x_nest$xa_da, # df/da = M1%*%b * ds/da
              store$X0) # df/db = Ma 
  
  return(list("JJ" = JJ, "xa" = x_nest$xa) )
}

.get.jacobian.inter_mgks <- function(object, data, param){
  
  na <- length(object$xt$si$alpha)
  # Single index and spline coefficients
  beta <- param[ -(1:na) ]
  
  x_nest <- Predict.matrix.nested(object, data = data, get.xa = TRUE)
  
  store <- object$xt$basis$evalX(z = x_nest$xa, t = object$xt$si$t, deriv = 1)
  X1beta <- drop(store$X1 %*% beta)
  
  JJ <- cbind(X1beta * x_nest$xa, # df/da = M1%*%b * ds/da0 (where ds/da0 = s because s = exp(a0) xa)
              X1beta * x_nest$xa_da, # df/da = M1%*%b * ds/da
              store$X0) # df/db = Ma 
  
  return(list("JJ" = JJ, "xa" = x_nest$xa) )
}

.get.jacobian.inter_dlinear <- function(object, data, param){
  
  na   <- length(object$xt$si$alpha)     # na1 + na2
  beta <- param[ -(1:na) ]
  
  x_nest <- Predict.matrix.nested(object, data = data, get.xa = TRUE)
  
  store <- object$xt$basis$evalX(z1 = x_nest$z1, z2 = x_nest$z2, deriv = 1)
  
  # df/dz_k = (dX/dz_k) %*% beta
  f1_1 <- drop( store$X1$dz1 %*% beta )
  f1_2 <- drop( store$X1$dz2 %*% beta )
  
  # df/dalpha_k = (df/dz_k) * (dz_k/dalpha_k) ; cross blocks are exactly zero
  JJ <- cbind(f1_1 * x_nest$xa_da[[1]],
              f1_2 * x_nest$xa_da[[2]],
              store$X0)
  
  if( ncol(JJ) != length(param) ){
    stop("Jacobian has ", ncol(JJ), " columns but param has length ",
         length(param), ".")
  }
  
  return( list("JJ" = JJ, "xa" = x_nest$xa) )
}
