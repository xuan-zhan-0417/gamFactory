#' Predict using 2D inter mgks effects (ANOVA Tensor Product)
#' 
#' @param object Smooth spec object containing stored smooth info
#' @param data Data frame containing prediction variables
#' @param get.xa Logical; if TRUE, returns inner linear predictor and its Jacobian matrix wrt alpha
#' 
#' @noRd
.predict.matrix.inter_mgks <- function(object, data, get.xa = FALSE){
  
  # =========================================================================
  # 1. Extract smooth info and alpha parameters
  # =========================================================================
  si <- object$xt$si
  
  # Ensure xm and tm exist for postproc_gam_nl compatibility
  if( is.null(si$xm) ){
    si$xm <- 0
  }
  if( is.null(si$tm) ){
    si$tm <- 0
  }
  
  alpha <- si$alpha
  a0 <- alpha[1]       # Variance scale parameter (log scale)
  a1 <- alpha[-1]      # MGKS kernel bandwidth parameters (beta)
  
  term_x <- object$term[1]
  term_t <- object$term[2]
  
  # =========================================================================
  # 2. Extract and split data (X_mat and t_doy)
  # =========================================================================
  Xi <- data[[term_x]]
  nms <- colnames(Xi)
  
  # Extract y0 (reference observations at historical locations)
  y0 <- Xi[ , which(nms == "y"), drop = FALSE]
  if( !ncol(y0) ){
    y0 <- si$y0
  }
  if( is.null(y0) ){
    y0 <- si$x
  }
  
  # Extract pairwise distance matrices (d1, d2, ...) sequentially
  Dist <- list()
  kk <- 1
  while( TRUE ){
    idx <- which(startsWith(nms, "d") & endsWith(nms, as.character(kk)) & sapply(nms, function(.x) nchar(.x) == 2))
    if( !length(idx) ){
      break
    }
    Dist[[kk]] <- Xi[ , idx, drop = FALSE]
    kk <- kk + 1
  }
  
  t_vec <- data[[term_t]]
  
  # =========================================================================
  # 3. Inner model: MGKS kernel smoothing and derivatives
  # =========================================================================
  # Evaluate MGKS C++ backend / core function
  xsm_list <- mgks(y = y0, dist = Dist, beta = a1, deriv = get.xa)
  
  # Center using training mean (si$xm) and scale using exp(a0)
  g_unscaled <- xsm_list$d0 - si$xm
  xsm <- exp(a0) * g_unscaled
  
  # =========================================================================
  # 4. Return inner parameters and Jacobian if requested by optimizer
  # =========================================================================
  if(get.xa){ 
    # Calculate Jacobian xa_da for full alpha = [a0, a1]:
    # d(xsm) / d(a0) = exp(a0) * g_unscaled = xsm
    # d(xsm) / d(a1) = exp(a0) * (dg / d_beta)
    d_a0 <- xsm
    d_a1 <- exp(a0) * xsm_list$d1
    
    return(list(xa = xsm,
                xa_da = cbind(d_a0, d_a1)))
  }
  
  # =========================================================================
  # 5. Build outer 2D model matrix X_2D
  # =========================================================================
  # Center time variable using stored training mean (si$tm)
  t_var <- t_vec - si$tm  
  
  # Compute outer 2D model matrix using the ANOVA tensor product basis evaluator
  X0 <- object$xt$basis$evalX(z1 = xsm, z2 = t_var, deriv = 0)$X0
  
  # Total design matrix: pad with zero columns corresponding to non-linear parameters alpha
  # (Standard gamFactory trick to allow predict.gam to align dimensions)
  Xtot <- cbind(matrix(0, nrow(X0), length(alpha)), X0) 
  
  attr(Xtot, "inner_linpred_unscaled") <- g_unscaled
  
  return(Xtot)
}