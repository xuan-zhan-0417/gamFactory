#' Nested single-index in tensor product
#' 
#' @name smooth.construct.inter_linear.smooth.spec
#' @rdname smooth.construct.inter_linear.smooth.spec
#' @importFrom MASS Null
#' @importFrom Matrix rankMatrix
#' @export
#'
smooth.construct.inter_linear.smooth.spec <- function(object, data, knots){
  
  # =========================================================================
  # initialise si
  # =========================================================================
  if(length(object$term) != 2) {
    stop("The smooth effect must contain exactly two terms: the multivariate matrix and the time variable.")
  }
  term_x <- object$term[1]
  term_t <- object$term[2]
  
  si <- object$xt$si
  if( is.null(si) ){ si <- object$xt$si <- list() }
  
  # =========================================================================
  # centering to two marginal variable, and save the mean for prediction
  # =========================================================================
  Xi <- data[[term_x]]
  t_vec <- data[[term_t]]
  
  t_mean <- mean(t_vec)
  data[[term_t]] <- t_vec - t_mean
  si$t <- data[[term_t]]
  si$tm <- t_mean 
  
  Xi <- scale(Xi, scale = FALSE)
  si$xm <- attr(Xi, "scaled:center") #overwrite data[[term_x]] later
  
  di <- ncol(Xi)
  n <- nrow(Xi)
  
  # =========================================================================
  # # diag penalty matrix
  # =========================================================================
  Si <- si$S
  no_pen <- is.null(Si) && is.null(si$pord)
  
  if( no_pen ){ 
    # case [a]: no penalty
    si$X <- Xi
    si$B <- diag(nrow = ncol(Xi))
    si$rank <- 0 
  } else {
    if( is.null(Si) ){ 
      # case [b]: P-splines penalty
      Si <- .psp(d = di, ord = si$pord)
      rankSi <- ncol(Xi) - si$pord
    } else { 
      # case [c]: designed penalty
      rankSi <- rankMatrix(Si)
    }
    si <- append(si, gamFactory:::.diagPen(X = Xi, S = Si, r = rankSi))
  }
  
  # =========================================================================
  # initialize alpha and a0 with B matrix
  # =========================================================================
  # full_alpha = alpha + a0
  if( is.null(si$a0) ){
    if( no_pen ){
      si$a0 <- rep(0, di)
    } else {
      si$a0 <- rep(1, di)
    }
  }
  
  if( is.null(si$alpha) ){ 
    if( is.null(si$a0) || all(si$a0 == 0) ){
      si$alpha <- rep(1, di) 
    } else {
      si$alpha <- rep(0, di) 
    }
  }
  
  si$alpha <- solve(si$B, si$alpha)
  si$a0 <- solve(si$B, si$a0)
  
  # impose variance constraint: var(X * (alpha + a0)) = 1
  tmp <- sd(si$X %*% (si$alpha + si$a0))
  si$alpha <- si$alpha / tmp
  si$a0 <- si$a0 / tmp
  
  # =========================================================================
  # z1 = x * t(alpha) and build design matrix X_2D
  # =========================================================================
  ax <- drop( si$X %*% (si$alpha + si$a0) ) 
  data[[term_x]] <- ax
  out <- .build_n_inter_bspline_basis(object = object, data = data, knots = knots, si = si, nested = c(TRUE, FALSE))
  
  # # =======================================================================
  # assemble penalty matrix
  # =========================================================================
  if( !no_pen ){
    # add alpha penalty
    dsmo <- out$bs.dim - di 
    si <- out$xt$si
    
    alpha_penalty_padded <- rbind(
      cbind(si$S, matrix(0, di, dsmo)),
      cbind(matrix(0, dsmo, di), matrix(0, dsmo, dsmo))
    )
    
    # out$S[[3]] <- alpha_penalty_padded
    out$S[[length(out$S) + 1]] <- alpha_penalty_padded
  
    out$null.space.dim <- out$null.space.dim + (out$bs.dim - si$rank)
    out$rank <- c(out$rank, si$rank)
  }
  
  class(out) <- c("inter_linear", "nested")
  
  # out$repara = FALSE
  # out$nl.reg <- TRUE

  return( out )
}