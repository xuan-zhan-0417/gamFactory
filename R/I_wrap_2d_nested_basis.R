.wrap_2d_nested_basis <- function(basis1, basis2, di_total) {
  force(basis1); force(basis2); force(di_total)
  
  evalX <- function(x1, x2, deriv = 0) {
    withCallingHandlers({
      # marginal value
      o1 <- basis1$evalX(x = x1, deriv = deriv)
      o2 <- basis2$evalX(x = x2, deriv = deriv)
      
      #Row-wise Kronecker
      X_2D_pred <- mgcv::tensor.prod.model.matrix(list(o1$X0, o2$X0))
      X_2D_pred <- cbind(o1$X0, o2$X0, X_2D_pred)

      
      # only return X0, should have X1, X2,... ideally
      return(list(X0 = X_2D_pred))
      
    }, warning = function(w) {
      if (length(grep("there is \\*no\\* information about some basis coefficients", conditionMessage(w)))) {
        invokeRestart("muffleWarning")
      }
    })
  }
  
  return(list("evalX" = evalX))
}