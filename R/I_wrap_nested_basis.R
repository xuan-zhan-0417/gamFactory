
.wrap_nested_basis <- function(b, P, Xth, add_slope){
  
  force(b); force(P); force(Xth)
  
  evalX <- function(x, deriv){
    withCallingHandlers({
      o <- b$evalX(x = x, deriv = deriv) # Get raw basis & derivatives
      if(add_slope){ # Add a slope in the last column: 1st derivative is 1, rest 0.
        o$X0 <- cbind(o$X0, x)
        if(deriv){
          o$X1 <- cbind(o$X1, 1)
          if(deriv >= 2){
            for(ii in 2:deriv){
              o[[paste0("X", ii)]] <- cbind(o[[paste0("X", ii)]], 0)
            }
          }
        }
      }
      o <- lapply(o, function(X) X %*% P)       # Reparametrise
      o <- linextr(x = x, b = o, th = b$krange, # Linearly extrapolate: NOTE method = "smooth" won't work if add_slope == TRUE
                   Xbo = Xth$X0%*%P, Xbo1 = Xth$X1%*%P, method = "simple")
      o
    }, warning = function(w) {
      if (length(grep("there is \\*no\\* information about some basis coefficients", conditionMessage(w)))){
        invokeRestart("muffleWarning")
      }
    })
  }
  
  out <- list("evalX" = evalX)

  return(out)
}

# The 1D outer basis in the layout of the 2D one (.wrap_2d_nested_basis), so that the generic code treats both
# the same way: evalX(z1, deriv) returns X0 and the lists X1 = list(dz1), X2 = list(dz1_z1), X3 = list(dz1_z1_z1)
.wrap_1d_nested_basis <- function(b){

  force(b)

  evalX <- function(z1, deriv = 0){
    o <- b$evalX(x = z1, deriv = deriv)
    out <- list(X0 = o$X0)
    if( deriv >= 1 ){ out$X1 <- list(dz1 = o$X1) }
    if( deriv >= 2 ){ out$X2 <- list(dz1_z1 = o$X2) }
    if( deriv >= 3 ){ out$X3 <- list(dz1_z1_z1 = o$X3) }
    out
  }

  list(evalX = evalX)
}
