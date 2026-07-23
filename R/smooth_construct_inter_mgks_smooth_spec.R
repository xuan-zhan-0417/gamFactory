#' Nested interactive mgks smoothing effect
#' 
#' @name smooth.construct.inter_mgks.smooth.spec
#' @rdname smooth.construct.inter_mgks.smooth.spec
#' @importFrom MASS Null
#' @importFrom Matrix rankMatrix
#' @importFrom stats sd
#' @export
#'
smooth.construct.inter_mgks.smooth.spec <- function(object, data, knots){
  
  # =========================================================================
  # 1. Check and initialize terms
  # =========================================================================
  if(length(object$term) != 2) {
    stop("The inter_mgks smooth effect must contain exactly two terms: the multivariate matrix (X_mat) and the time variable (t_doy).")
  }
  term_x <- object$term[1]
  term_t <- object$term[2]
  
  si <- object$xt$si
  if( is.null(si) ){ si <- object$xt$si <- list() }
  
  # =========================================================================
  # 2. Centering time variable t_doy and record mean for prediction
  # =========================================================================
  t_vec <- data[[term_t]]
  t_mean <- mean(t_vec)
  data[[term_t]] <- t_vec - t_mean
  si$t <- data[[term_t]]
  si$tm <- t_mean 
  
  # =========================================================================
  # 3. Parse composite matrix Xi (Extract y0 and pairwise distance matrices)
  # =========================================================================
  Xi <- data[[term_x]]
  nms <- colnames(Xi)
  
  # Extract response variable y0 at historical/reference locations
  y0 <- Xi[ , which(nms == "y"), drop = FALSE]
  if( !ncol(y0) ){
    y0 <- si$y0
  }
  
  # Extract distance matrices (d1, d2, ...) sequentially
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
  
  si$x <- y0
  si$dist <- Dist
  
  # =========================================================================
  # 4. Initialize and compute inner MGKS smooth predictor & alpha
  # =========================================================================
  alpha <- si$alpha
  if( is.null(alpha) ){ 
    # Cold start: bandwidths set to -log(sd(Dist)/10)
    init_beta <- -log(sapply(si$dist, sd) / 10)
    g <- mgks(y = si$x, dist = Dist, beta = init_beta)$d0
    alpha <- si$alpha <- c(log(1 / sd(g)), init_beta) 
  } else {
    # Warm start / Iterative step inside efs optimizer
    if( length(alpha) == length(Dist) + 1 ){
      g <- mgks(y = si$x, dist = Dist, beta = alpha[-1])$d0
      alpha[1] <- log(1 / sd(g))
      si$alpha <- alpha
    } else {
      g <- mgks(y = si$x, dist = Dist, beta = alpha)$d0
      alpha <- si$alpha <- c(log(1 / sd(g)), alpha)
    }
  }
  
  # Center and scale inner linear predictor: sd constraint applied
  ax <- exp(si$alpha[1]) * (g - mean(g))
  data[[term_x]] <- ax
  
  # =========================================================================
  # 5. Build 2D interactive B-spline basis matrix X_2D
  # =========================================================================
  out <- .build_n_inter_bspline_basis(object = object, data = data, knots = knots, si = si)
  
  # =========================================================================
  # 6. Resolve overlapping constraints in outer penalty matrices
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
  # 7. Assign classes and return
  # =========================================================================
  class(out) <- c("inter_mgks", "nested")
  
  return( out )
}