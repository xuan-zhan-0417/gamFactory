#' @export
check_wrap_z_deriv_numDeriv_mean <- function(Xi, t, basis, param, a0 = NULL,
                                             deriv = 3,
                                             method = "Richardson",
                                             method.args = list()) {
  
  na <- ncol(Xi)
  alpha <- param[1:na]
  beta  <- param[-(1:na)]
  
  if (is.null(a0)) {
    a0 <- alpha * 0
  }
  
  z0 <- drop(Xi %*% (alpha + a0))
  
  # analytical derivatives
  st <- basis$evalX(z = z0, t = t, deriv = deriv)
  
  # scalar function: common shift c applied to all z_i
  eta_mean_at <- function(c) {
    z_eval <- z0 + c
    st_eval <- basis$evalX(z = z_eval, t = t, deriv = 0)
    eta_eval <- drop(st_eval$X0 %*% beta)
    mean(eta_eval)
  }
  
  out <- list()
  
  if (deriv >= 1) {
    f1_ex <- mean(drop(st$X1 %*% beta))
    
    f1_fd <- numDeriv::grad(
      func = eta_mean_at,
      x = 0,
      method = method,
      method.args = method.args
    )
    
    out$f1 <- c(
      EX = f1_ex,
      FD = f1_fd,
      abs_err = abs(f1_ex - f1_fd),
      rel_err = abs(f1_ex - f1_fd) / max(1, abs(f1_fd))
    )
  }
  
  if (deriv >= 2) {
    f2_ex <- mean(drop(st$X2 %*% beta))
    
    f2_fd <- as.numeric(numDeriv::hessian(
      func = eta_mean_at,
      x = 0,
      method = method,
      method.args = method.args
    ))
    
    out$f2 <- c(
      EX = f2_ex,
      FD = f2_fd,
      abs_err = abs(f2_ex - f2_fd),
      rel_err = abs(f2_ex - f2_fd) / max(1, abs(f2_fd))
    )
  }
  
  if (deriv >= 3) {
    
    d2_mean_at <- function(c) {
      as.numeric(numDeriv::hessian(
        func = eta_mean_at,
        x = c,
        method = method,
        method.args = method.args
      ))
    }
    
    f3_ex <- mean(drop(st$X3 %*% beta))
    
    f3_fd <- numDeriv::grad(
      func = d2_mean_at,
      x = 0,
      method = method,
      method.args = method.args
    )
    
    out$f3 <- c(
      EX = f3_ex,
      FD = f3_fd,
      abs_err = abs(f3_ex - f3_fd),
      rel_err = abs(f3_ex - f3_fd) / max(1, abs(f3_fd))
    )
  }
  
  out
}