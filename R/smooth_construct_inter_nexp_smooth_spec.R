#' Nested interactive adaptive exponential smoothing effect
#' 
#' @name smooth.construct.inter_nexp.smooth.spec
#' @rdname smooth.construct.inter_nexp.smooth.spec
#' @importFrom MASS Null
#' @importFrom Matrix rankMatrix
#' @export
#'
smooth.construct.inter_nexp.smooth.spec <- function(object, data, knots){
  
  # =========================================================================
  # 1. Validate and extract variables
  # =========================================================================
  if(length(object$term) != 2) {
    stop("The smooth effect must contain exactly two terms: the multivariate matrix (X_mat) and the time variable (t_doy).")
  }
  term_x <- object$term[1] # X_mat (x, w1, w2)
  term_t <- object$term[2] # t_doy
  
  si <- object$xt$si
  if( is.null(si) ){ si <- object$xt$si <- list() }
  
  # =========================================================================
  # 2. Center the time marginal variable
  # =========================================================================
  t_vec <- data[[term_t]]
  t_mean <- mean(t_vec)
  data[[term_t]] <- t_vec - t_mean
  si$t <- data[[term_t]]
  si$tm <- t_mean 
  
  # =========================================================================
  # 3. Split X_mat into raw x and weight matrix W (w1, w2)
  # =========================================================================
  X_mat <- data[[term_x]]
  if (ncol(X_mat) < 2) {
    stop("X_mat must have at least 2 columns: x (data to smooth) and w variables.")
  }
  
  x_raw <- X_mat[, 1]
  W_mat <- X_mat[, -1, drop = FALSE] # Contains w1 (intercept) and w2
  n <- nrow(X_mat)
  
  # di is the total number of inner parameters: alpha_scale + alpha_intercept + alpha_w
  di <- ncol(W_mat) + 1 
  
  # =========================================================================
  # 4. Handle inner penalty (for alpha_w) and reparameterization
  # =========================================================================
  Si <- si$S
  no_pen <- is.null(Si) && is.null(si$pord)
  
  if( no_pen ){ 
    si$X <- W_mat
    si$B <- diag(nrow = ncol(W_mat))
    si$rank <- 0 
  } else {
    # not 100% sure this part
    rankSi <- ifelse(is.null(Si), ncol(W_mat) - si$pord, rankMatrix(Si))
    if(is.null(Si)) Si <- .psp(d = ncol(W_mat), ord = si$pord)
    si <- append(si, gamFactory:::.diagPen(X = W_mat, S = Si, r = rankSi))
  }
  
  # =========================================================================
  # 5. Initialize alpha,perform initial smoothing, center and scale
  # =========================================================================
  # Need to initialize inner coefficients?
  alpha <- si$alpha
  if( is.null(alpha) ){ 
    # alpha[1] s.t. sd(inner_lin_pred) = 1 (target variance)
    g <- expsmooth(y = x_raw, Xi = si$X, beta = rep(0, di-1))$d0
    alpha <- si$alpha <- c(log(1/sd(g)), rep(0, di-1))
  } else {
    alpha <- solve(si$B) %*% alpha
    g <- expsmooth(y = x_raw, Xi = si$X, beta = alpha, times = times)$d0
    alpha <- si$alpha <- c(log(1/sd(g)), alpha)
  }
  
  # Center and scale the initialized inner linear predictor
  data[[term_x]] <- exp(alpha[1]) * (g - mean(g))
  
  # si$xm <- mean(g) # Save mean for prediction
  si$x_raw <- x_raw
  si$W_mat <- W_mat
  
  # =========================================================================
  # 7. Build custom 2D B-spline basis (X_2D)
  # =========================================================================
  out <- .build_n_inter_bspline_basis(object = object, data = data, knots = knots, si = si)
  
  # =========================================================================
  # 8. Deal with overlap in the external penalty matrices
  # =========================================================================
  n_mats <- length(out$S)
  matrix_dim <- nrow(out$S[[1]])
  
  diag_matrix <- sapply(out$S, diag)
  ol_ind <- which(rowSums(diag_matrix == 1) > 1)
  out$S[[n_mats + 1]] <- matrix(0, nrow = matrix_dim, ncol = matrix_dim)
  
  if (length(ol_ind) > 0) {
    diag(out$S[[n_mats + 1]])[ol_ind] <- 1
    for (i in 1:n_mats) {
      diag(out$S[[i]])[ol_ind] <- 0
    }
  }
  out$rank <- sapply(out$S, function(Sm) as.numeric(Matrix::rankMatrix(Sm)))
  
  # =========================================================================
  # 9. Assemble final block-diagonal penalty matrix
  # =========================================================================
  dsmo <- out$bs.dim - di 
  si <- out$xt$si
  
  if ( !no_pen ) {
    inner_pen <- rbind(0, cbind(0, si$S))
    
    alpha_penalty_padded <- rbind(
      cbind(inner_pen, matrix(0, di, dsmo)),
      cbind(matrix(0, dsmo, di), matrix(0, dsmo, dsmo))
    )
    
    out$S[[length(out$S) + 1]] <- alpha_penalty_padded
    
    out$null.space.dim <- out$null.space.dim + (out$bs.dim - si$rank)
    out$rank <- c(out$rank, si$rank)
  }
  
  # # debug test  
  # out$S <- NULL
  # out$S[[1]] <-  as.matrix(Matrix::bdiag(diag(0,nrow = 3, ncol = 3), diag(1,nrow = 15, ncol = 15)))
  # out$rank <- 15
  # # test end   
    
  class(out) <- c("inter_nexp", "nested")
  
  return( out )
}