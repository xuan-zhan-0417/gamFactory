#' Nested single-index in tensor product
#' 
#' @name smooth.construct.inter_linear.smooth.spec
#' @rdname smooth.construct.inter_linear.smooth.spec
#' @importFrom MASS Null
#' @importFrom Matrix rankMatrix
#' @export
#'
smooth.construct.inter_linear.smooth.spec <- function(object, data, knots){
  
  if(length(object$term) != 2) {
    stop("The smooth effect must contain exactly two terms.")
  }
  term_1 <- object$term[1]
  term_2 <- object$term[2]
  
  si <- object$xt$si
  if( is.null(si) ){ si <- object$xt$si <- list() }

  # margin 2 is a plain covariate (t) when it is a vector, or a matrix with a
  # single column; it is treated as a nested single index (si) otherwise.
  d2_raw <- if( is.null(dim(data[[term_2]])) ) 1L else ncol(data[[term_2]])
  nested_2 <- d2_raw > 1L

  d1 <- ncol(as.matrix(data[[term_1]]))
  
  # sanity checks on user-supplied dims (mirrors old inter_dlinear checks)
  if (!is.null(si$alpha_1) && length(si$alpha_1) != d1) stop(...)
  if (!is.null(si$a0_1)    && length(si$a0_1)    != d1) stop(...)
  if (!is.null(si$S_1) && (nrow(si$S_1) != d1 || ncol(si$S_1) != d1)) stop(...)
  
  # ---- margin 1: always a nested single index ------------------------------
  res1 <- .init_marginal_si(Xi = as.matrix(data[[term_1]]),
                            S = si$S_1, pord = si$pord_1,
                            a0 = si$a0_1, alpha = si$alpha_1)
  data[[term_1]] <- res1$ax
  
  if( nested_2 ){
    d2 <- ncol(as.matrix(data[[term_2]]))
    if (!is.null(si$alpha_2) && length(si$alpha_2) != d2) stop(...)
    if (!is.null(si$a0_2)    && length(si$a0_2)    != d2) stop(...)
    if (!is.null(si$S_2) && (nrow(si$S_2) != d2 || ncol(si$S_2) != d2)) stop(...)
    
    res2 <- .init_marginal_si(Xi = as.matrix(data[[term_2]]),
                              S = si$S_2, pord = si$pord_2,
                              a0 = si$a0_2, alpha = si$alpha_2)
    data[[term_2]] <- res2$ax
    
    si$X <- list(res1$X, res2$X)
    si$B <- list(res1$B, res2$B)
    si$S <- list(res1$S, res2$S)
    si$xm    <- c(res1$xm, res2$xm)
    si$alpha <- c(res1$alpha, res2$alpha)
    si$a0    <- c(res1$a0, res2$a0)
    si$rank  <- c(res1$rank, res2$rank)
    si$na1 <- d1; si$na2 <- d2; si$na <- d1 + d2
    
  } else {
    t_vec <- data[[term_2]]
    if( !is.null(dim(t_vec)) && ncol(as.matrix(t_vec)) > 1 ){
      stop("term[2] must be a single numeric vector when nested_2 = FALSE.")
    }
    t_mean <- mean(t_vec)
    data[[term_2]] <- t_vec - t_mean
    si$t  <- data[[term_2]]
    si$tm <- t_mean
    
    si$X <- list(res1$X, NULL)
    si$B <- list(res1$B, NULL)
    si$S <- list(res1$S, NULL)
    si$xm    <- res1$xm
    si$alpha <- res1$alpha
    si$a0    <- res1$a0
    si$rank  <- c(res1$rank, 0)
    si$na1 <- d1; si$na2 <- 0; si$na <- d1
  }
  
  out <- .build_n_inter_bspline_basis(object = object, data = data, knots = knots,
                                      si = si, nested = c(TRUE, nested_2))
  
  # ---- penalty: pad each margin's alpha penalty into bs.dim space ----------
  .pad_alpha <- function(Sk, off, dk){
    P <- matrix(0, out$bs.dim, out$bs.dim)
    ii <- off + seq_len(dk)
    P[ii, ii] <- Sk
    P
  }
  added_rank <- 0
  if( !is.null(si$S[[1]]) && isTRUE(si$rank[1] > 0) ){
    out$S[[length(out$S)+1]] <- .pad_alpha(si$S[[1]], 0, d1)
    out$rank <- c(out$rank, si$rank[1]); added_rank <- added_rank + si$rank[1]
  }
  if( nested_2 && !is.null(si$S[[2]]) && isTRUE(si$rank[2] > 0) ){
    out$S[[length(out$S)+1]] <- .pad_alpha(si$S[[2]], d1, si$na2)
    out$rank <- c(out$rank, si$rank[2]); added_rank <- added_rank + si$rank[2]
  }
  if( added_rank > 0 ){
    out$null.space.dim <- out$null.space.dim + (out$bs.dim - added_rank)
  }
  
  # ---- stack marginal si design matrices ------------------------------------
  if( nested_2 ){
    out$xt$si$X <- cbind(out$xt$si$X[[1]], out$xt$si$X[[2]])
    attr(out$xt$si$X, "na1") <- d1
    attr(out$xt$si$X, "na2") <- si$na2
  } else {
    out$xt$si$X <- out$xt$si$X[[1]]
  }
  
  class(out) <- c("inter_linear", "nested")
  return( out )
}