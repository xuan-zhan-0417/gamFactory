.wrap_2d_nested_basis <- function(b1, b2){
  
  force(b1); force(b2)

  evalX <- function(z, t, deriv = 0){

    out1 <- b1$evalX(x = z, deriv = deriv) 
    out2 <- b2$evalX(x = t, deriv = 0)
    X2_0 <- out2$X0 
    
    n <- length(z)
    p2 <- ncol(X2_0)
    
    out_2d <- list()
    
    for(ii in 0:deriv){
      name <- paste0("X", ii)
      X1_d <- out1[[name]]
      
      if(ii == 0){
        # ---------------------------------------------------------
        # deriv = 0: [X1_0, X2_0, X1_0 ⊙ X2_0]
        # ---------------------------------------------------------
        X_inter <- mgcv::tensor.prod.model.matrix(list(X1_d, X2_0))
        out_2d[[name]] <- cbind(X1_d, X2_0, X_inter)
        
      } else {
        # ---------------------------------------------------------
        # [X1_d, 0, X1_d ⊙ X2_0]
        # X2_d = 0
        # ---------------------------------------------------------
        X2_d <- matrix(0, nrow = n, ncol = p2)
        X_inter_d <- mgcv::tensor.prod.model.matrix(list(X1_d, X2_0))
        out_2d[[name]] <- cbind(X1_d, X2_d, X_inter_d)
      }
    }
    
    return(out_2d)
  }
  
  out <- list("evalX" = evalX)
  return(out)
}
