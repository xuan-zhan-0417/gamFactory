#' Predict using 2D double-nested effects (both margins are single indices)
#' 
#' @param object smooth object of class "inter_dlinear".
#' @param data data frame / list holding the two inner covariate matrices.
#' @param get.xa if TRUE, return the inner indices and their Jacobians
#'               instead of the model matrix.
#' @noRd
.predict.matrix.inter_dlinear <- function(object, data, get.xa = FALSE){
  
  term_x_1 <- object$term[1]   # e.g. "X1_mat"
  term_x_2 <- object$term[2]   # e.g. "X2_mat"
  
  si <- object$xt$si
  
  na1 <- attr(si$X, "na1")
  na2 <- attr(si$X, "na2")
  if( is.null(na1) || is.null(na2) ){
    stop("Attributes 'na1' / 'na2' are missing from object$xt$si$X.")
  }
  di <- na1 + na2                     # total number of inner coefficients
  
  if( length(si$alpha) != di || length(si$a0) != di || length(si$xm) != di ){
    stop("si$alpha / si$a0 / si$xm must all have length na1 + na2 = ", di, ".")
  }
  if( !is.list(si$B) || length(si$B) != 2 ||
      any(vapply(si$B, is.null, logical(1))) ){
    stop("si$B must be a list of the two marginal reparametrisation matrices ",
         "(use the identity when a margin is unpenalised, never NULL).")
  }
  
  ik <- list(1:na1, (na1 + 1):di)     # alpha slices
  tm <- c(term_x_1, term_x_2)
  
  Xi <- vector("list", 2)             # marginal inner matrices = dz_k / dalpha_k
  za <- vector("list", 2)             # marginal inner indices z_k
  
  for(k in 1:2){
    
    Xk <- as.matrix( data[[ tm[k] ]] )
    if( ncol(Xk) != length(ik[[k]]) ){
      stop("Term '", tm[k], "' has ", ncol(Xk), " columns but na", k,
           " = ", length(ik[[k]]), ".")
    }
    
    # same reparametrisation as .predict.matrix.si: centre with xm, then rotate with B
    Xk <- t(t(Xk) - si$xm[ ik[[k]] ]) %*% si$B[[k]]
    
    Xi[[k]] <- Xk
    za[[k]] <- drop( Xk %*% (si$alpha[ ik[[k]] ] + si$a0[ ik[[k]] ]) )
  }
  
  if( get.xa ){
    return( list(xa    = cbind(z1 = za[[1]], z2 = za[[2]]),  # n x 2
                 xa_da = Xi,                                  # list(n x na1, n x na2)
                 z1 = za[[1]], z2 = za[[2]],
                 Xi_1 = Xi[[1]], Xi_2 = Xi[[2]],
                 na1 = na1, na2 = na2) )
  }
  
  X0 <- object$xt$basis$evalX(z1 = za[[1]], z2 = za[[2]], deriv = 0)$X0
  
  # [ 0_di , X_2D ] : predict.gam multiplies the leading block by alpha -> 0
  Xtot <- cbind(matrix(0, nrow(X0), di), X0)
  
  if( !is.null(object$bs.dim) && ncol(Xtot) != object$bs.dim ){
    stop("Predict matrix has ", ncol(Xtot), " columns but bs.dim = ",
         object$bs.dim, ". Check that bs.dim includes the di alpha columns.")
  }
  
  attr(Xtot, "inner_linpred_unscaled") <- cbind(z1 = za[[1]], z2 = za[[2]])
  
  return( Xtot )
}