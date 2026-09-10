#' Predict using 2D inter linear effects (ANOVA Tensor Product)
#' 
#' @param object 
#' @param data 
#' @param get.xa 
#' Predict using nested interactive exponential smoothing effects
#' 
#' @noRd
.predict.matrix.inter_nexp <- function(object, data, get.xa = FALSE){
  
  # =========================================================================
  # 1. Extract smooth info and alpha parameters
  # =========================================================================
  si <- object$xt$si
  
  # Ensure xm exists for postproc_gam_nl compatibility
  if( is.null(si$xm) ){
    si$xm <- 0
  }
  if( is.null(si$tm) ){
    si$tm <- 0
  }
  
  alpha <- si$alpha
  a0 <- alpha[1]      # alpha_scale
  a1 <- alpha[-1]     # [alpha_intercept, alpha_w]
  
  term_x <- object$term[1]
  term_t <- object$term[2]
  
  # =========================================================================
  # 2. Extract and split data
  # =========================================================================
  X_mat <- data[[term_x]]
  n <- nrow(X_mat)
  nms <- colnames(X_mat)
  x_raw <- as.vector( t(X_mat[ , which(nms == "y")]) )
  times <- NULL
  tmp <- which(nms == "times")
  if( length(tmp) ){
    times <- X_mat[ , tmp]
  }
  W_mat <- X_mat[ , which(nms == "x"), drop = FALSE]
  
  nrep <- ceiling( length(x_raw)/n )
  dXi <- ncol(W_mat)/nrep
  if(nrep > 1){
    tmp <- rep(1:dXi, nrep)
    W_mat <- apply(W_mat, 1, function(x) do.call("cbind", tapply(x, tmp, I)), simplify = FALSE)
    W_mat <- do.call("rbind", W_mat)
  }
  
  # =========================================================================
  # 3. Inner model: Exponential smoothing and derivatives
  # =========================================================================
  # Reparameterize W_mat using the basis transformation matrix B from constructor
  W_mat_rot <- W_mat %*% si$B
  
  # Evaluate the C++ backend (assuming expsmooth handles the sigmoid internally or 
  # beta weights linear predictor. If expsmooth is standard, ensure it matches your omega definition).
  xsm_list <- expsmooth(y = x_raw, Xi = W_mat_rot, beta = a1, times = times, deriv = get.xa)
  
  # Center and scale the smoothed variable
  xsm_unscaled <- xsm_list$d0 - si$xm
  xsm <- exp(a0) * xsm_unscaled
  
  # =========================================================================
  # 4. Return inner parameters and Jacobian if requested by optimizer
  # =========================================================================
  if(get.xa){ 
    # The Jacobian xa_da needs to account for ALL elements in alpha.
    # d(xsm)/d(a0) = exp(a0) * xsm_unscaled = xsm
    # d(xsm)/d(a1) = exp(a0) * xsm_list$d1
    
    return(list(xa = xsm,
                xa_da = exp(a0) * xsm_list$d1))
  }
  
  # =========================================================================
  # 5. Build outer 2D model matrix X_2D
  # =========================================================================
  t_var <- data[[term_t]]  - si$tm  
  
  # Compute outer model matrix using the basis evaluator created in constructor
  # For the custom 2D tensor, evalX typically expects the updated data environment
  X0 <- object$xt$basis$evalX(z1 = xsm, z2 = t_var, deriv = 0)$X0
  
  # Total model matrix is X0 preceded by a matrix of zeros. 
  # predict.gam will multiply the latter by alpha, which will have no effect (the standard gamFactory trick).
  Xtot <- cbind(matrix(0, nrow(X0), length(alpha)), X0) 
  
  attr(Xtot, "inner_linpred_unscaled") <- xsm_unscaled
  
  return(Xtot)
  
}