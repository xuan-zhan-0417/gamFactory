######
# Penalty on the variance of BOTH inner single indices.
# P = (var(z1) - v)^2 + (var(z2) - v)^2 ; separable, hence a block-diagonal Hessian.
#
.pen_var_inter_dlinear <- function(o, v, deriv = 0){
  
  na  <- o$na
  na1 <- o$na1
  na2 <- o$na2
  if( is.null(na1) || is.null(na2) ){
    stop("o$na1 / o$na2 not set by eff_inter_dlinear().")
  }
  if( is.null(o$a0) || length(o$a0) != na ){
    stop("o$a0 must be a numeric vector of length na1 + na2 (normalise it in eff_inter_dlinear()).")
  }
  
  a <- o$param[1:na] + o$a0
  x <- o$store$Xi
  n <- nrow(x)
  
  idx <- list(1:na1, (na1 + 1):na)
  
  l0 <- 0
  l1 <- if( deriv )     numeric(na)          else NULL
  l2 <- if( deriv > 1 ) matrix(0, na, na)    else NULL
  
  for(k in 1:2){
    
    ik <- idx[[k]]
    xk <- x[, ik, drop = FALSE]
    ak <- a[ik]
    
    axk  <- drop( xk %*% ak )
    vhat <- sum(axk^2)/n - mean(axk)^2      # same expression as .pen_var_si
    
    l0 <- l0 + (vhat - v)^2
    
    if( deriv ){
      Sk <- cov(xk) * (n - 1) / n
      Sa <- Sk %*% ak
      l1[ik] <- 4 * (vhat - v) * drop(Sa)
      if( deriv > 1 ){
        # cross block d2/(da_1 da_2) is exactly zero: the penalty is separable
        l2[ik, ik] <- 8 * tcrossprod(Sa) + 4 * (vhat - v) * Sk
      }
    }
  }
  
  # NOTE: pen[[ii]]$d3 is never consumed by build_family_nl(), so we leave it NULL.
  return( list("d0" = l0, "d1" = l1, "d2" = l2, "d3" = NULL) )
}

# Derivatives of the penalty's Hessian w.r.t. the smoothing parameters
.pen_var_inter_dlinear_outer <- function(o, v, DaDr){
  
  na  <- o$na
  na1 <- o$na1
  na2 <- o$na2
  
  a <- o$param[1:na] + o$a0
  x <- o$store$Xi
  n <- nrow(x)
  m <- ncol(DaDr)
  
  idx <- list(1:na1, (na1 + 1):na)
  d1H <- lapply(1:m, function(nouse) matrix(0, na, na))
  
  for(k in 1:2){
    ik <- idx[[k]]
    xk <- x[, ik, drop = FALSE]
    ak <- a[ik]
    Sk <- cov(xk) * (n - 1) / n
    Sa <- Sk %*% ak
    for(ii in 1:m){
      SD <- Sk %*% DaDr[ik, ii]
      d1H[[ii]][ik, ik] <- 8 * ( tcrossprod(Sa, SD) + tcrossprod(SD, Sa) +
                                   Sk * drop(crossprod(ak, SD)) )
    }
  }
  
  return( d1H )
}