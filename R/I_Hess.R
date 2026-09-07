
.Hess.stand_stand <- function(o1, o2, llk){
  
  out <- crossprod(o2$store$X, llk$d2 * o1$store$X)
  
  return( out )
  
}


.Hess.si_stand <- .Hess.nexpsm_stand <- function(o1, o2, llk){

  out <- cbind(crossprod(o2$store$X, llk$d2 * o1$store$f1 * o1$store$g1), 
               crossprod(o2$store$X, llk$d2 * o1$store$X0))
    
 return( out )
  
}

.Hess.stand_si <- .Hess.stand_nexpsm <- function(o1, o2, llk){

  # Swapped o1 with o2 and lp1 with lp2
  t( .Hess.si_stand(o2, o1, llk) ) 
  
}

.Hess.si_nexpsm <- function(o1, o2, llk){
  
  lgg1 <- llk$d2 * o1$store$f1 * o2$store$f1
  
  # o2 vertical, o1 horizontal
  out <- rbind(cbind(crossprod(o2$store$g1, lgg1 * o1$store$g1), 
                     crossprod(o2$store$g1, llk$d2 * o2$store$f1 * o1$store$X0)),
               cbind(crossprod(o2$store$X0, (llk$d2 * o1$store$f1) * o1$store$g1), 
                     crossprod(o2$store$X0, llk$d2 * o1$store$X0)))
  
  
  return( out )
  
}

.Hess.nexpsm_si <- function(o1, o2, llk){
  
  # Swapped o1 with o2 and lp1 with lp2
  out <- t( .Hess.si_nexpsm(o1 = o2, o2 = o1, llk) )
  
  return( out )
  
}

.Hess.si_si <- .Hess.nexpsm_nexpsm <- function(o1, o2, llk){
  
  lgg1 <- llk$d2 * o1$store$f1 * o2$store$f1
  
  # o2 vertical, o1 horizontal
  out <- rbind(cbind(crossprod(o2$store$g1, lgg1 * o1$store$g1), 
                     crossprod(o2$store$g1, llk$d2 * o2$store$f1 * o1$store$X0)),
               cbind(crossprod(o2$store$X0, (llk$d2 * o1$store$f1) * o1$store$g1), 
                     crossprod(o2$store$X0, llk$d2 * o1$store$X0)))
  
  
  return( out )
  
}

#' @noRd
.jac_inter_dlinear <- function(o){
  Xi  <- o$store$Xi
  na1 <- o$na1; na2 <- o$na2
  cbind(o$store$f1$f1_1 * Xi[, 1:na1, drop = FALSE],
        o$store$f1$f1_2 * Xi[, (na1 + 1):(na1 + na2), drop = FALSE],
        o$store$X0)
}

#' @noRd
.jac_si_like <- function(o) cbind(o$store$f1 * o$store$g1, o$store$X0)

.Hess.inter_dlinear_stand <- function(o1, o2, llk){
  crossprod(o2$store$X, llk$d2 * .jac_inter_dlinear(o1))
}
.Hess.stand_inter_dlinear <- function(o1, o2, llk){
  t( .Hess.inter_dlinear_stand(o1 = o2, o2 = o1, llk) )
}

.Hess.inter_dlinear_inter_dlinear <- function(o1, o2, llk){
  crossprod(.jac_inter_dlinear(o2), llk$d2 * .jac_inter_dlinear(o1))
}

.Hess.inter_dlinear_si <- .Hess.inter_dlinear_nexpsm <- function(o1, o2, llk){
  crossprod(.jac_si_like(o2), llk$d2 * .jac_inter_dlinear(o1))
}
.Hess.si_inter_dlinear <- .Hess.nexpsm_inter_dlinear <- function(o1, o2, llk){
  t( .Hess.inter_dlinear_si(o1 = o2, o2 = o1, llk) )
}