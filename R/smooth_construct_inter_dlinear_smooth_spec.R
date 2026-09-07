#' Nested single-index in tensor product
#' 
#' @name smooth.construct.inter_dlinear.smooth.spec
#' @rdname smooth.construct.inter_dlinear.smooth.spec
#' @importFrom MASS Null
#' @importFrom Matrix rankMatrix bdiag
#' @export
smooth.construct.inter_dlinear.smooth.spec <- function(object, data, knots){  
  # =========================================================================
  # 1. check and initialization
  # =========================================================================
  if(length(object$term) != 2) {
    stop("The smooth effect must contain exactly two terms.")
  }
  term_x_1 <- object$term[1]
  term_x_2 <- object$term[2]
  
  si <- object$xt$si
  if( is.null(si) ){ si <- object$xt$si <- list() }
  
  d1 <- ncol(as.matrix(data[[term_x_1]]))
  d2 <- ncol(as.matrix(data[[term_x_2]]))
  
  # ---- two-sided inputs: kept separate, gathered into per-side lists -------
  alpha_in <- list(si$alpha_1, si$alpha_2)
  a0_in    <- list(si$a0_1,    si$a0_2)
  S_in     <- list(si$S_1,     si$S_2)
  pord_in  <- list(si$pord_1,  si$pord_2)
  d_in     <- c(d1, d2)
  
  # sanity check on the supplied dimensions
  for (k in 1:2) {
    if (!is.null(alpha_in[[k]]) && length(alpha_in[[k]]) != d_in[k]) {
      stop(sprintf("length(alpha_%d) must equal ncol(%s) = %d.",
                   k, object$term[k], d_in[k]))
    }
    if (!is.null(a0_in[[k]]) && length(a0_in[[k]]) != d_in[k]) {
      stop(sprintf("length(a0_%d) must equal ncol(%s) = %d.",
                   k, object$term[k], d_in[k]))
    }
    if (!is.null(S_in[[k]]) &&
        (nrow(S_in[[k]]) != d_in[k] || ncol(S_in[[k]]) != d_in[k])) {
      stop(sprintf("S_%d must be a %d x %d matrix.", k, d_in[k], d_in[k]))
    }
  }
  
  si$X  <- vector("list", 2)
  si$B  <- vector("list", 2)
  si$S  <- vector("list", 2)
  
  xm_out    <- numeric(0)
  alpha_out <- numeric(0)
  a0_out    <- numeric(0)
  rank_out  <- numeric(0)
  
  # =========================================================================
  # 2. loop marginal Single-Index
  # =========================================================================
  for (k in 1:2) {
    term_name <- object$term[k]
    
    # initialize marginal single-index for the k-th term
    res <- .init_marginal_si(
      Xi    = as.matrix(data[[term_name]]),
      S     = S_in[[k]],
      pord  = pord_in[[k]],
      a0    = a0_in[[k]],
      alpha = alpha_in[[k]]
    )
    
    # results ( use x[k] <- list(.) so that NULL keeps its slot )
    si$X[k]  <- list(res$X)
    si$B[k]  <- list(res$B)
    si$S[k]  <- list(res$S)
    
    xm_out    <- c(xm_out, res$xm)
    alpha_out <- c(alpha_out, res$alpha)
    a0_out    <- c(a0_out, res$a0)
    rank_out  <- c(rank_out, res$rank)
    
    # update data
    data[[term_name]] <- res$ax
  }
  
  # A NULL rank would be dropped by c(), leaving si$rank too short; downstream
  # si$rank[2] would then be NA and `if (NA > 0)` errors. Fail here instead.
  if (length(rank_out) != 2L) {
    stop(".init_marginal_si must return a numeric rank for each margin ",
         "(0 when the margin is unpenalised).")
  }
  if (length(xm_out) != d1 + d2 || length(alpha_out) != d1 + d2 ||
      length(a0_out) != d1 + d2) {
    stop("xm / alpha / a0 returned by .init_marginal_si must have length ",
         "ncol(Xi) for each margin.")
  }
  
  # update si : alpha = [alpha_1, alpha_2] as a single stacked vector
  si$xm    <- xm_out
  si$alpha <- alpha_out
  si$a0    <- a0_out
  si$rank  <- rank_out
  
  # keep the block sizes so that alpha can be split downstream
  si$na1 <- d1
  si$na2 <- d2
  si$na  <- d1 + d2
  
  # =========================================================================
  # 3. 2D B-spline design matrix
  # =========================================================================
  out <- .build_n_inter_bspline_basis(object = object, data = data,
                                      knots = knots, si = si,
                                      nested = c(TRUE, TRUE))
  
  # =========================================================================
  # 4. Penalty Matrix
  # =========================================================================
  no_pen <- is.null(si$S[[1]]) && is.null(si$S[[2]]) &&
    is.null(pord_in[[1]]) && is.null(pord_in[[2]])
  
  if( !no_pen ){
    
    # pad a marginal alpha penalty into the full bs.dim x bs.dim parameter
    # space. alpha occupies the leading di columns, matching
    # Parameter = [alpha, beta] and the zero padding in Predict.matrix.
    .pad_alpha <- function(Sk, off, dk){
      P  <- matrix(0, out$bs.dim, out$bs.dim)
      ii <- off + seq_len(dk)
      P[ii, ii] <- Sk
      P
    }
    
    # ONE smoothing parameter per margin: alpha_1 and alpha_2 get smoothed
    # independently. A margin with no penalty contributes nothing at all --
    # never a zero matrix, which would give mgcv a lambda that does nothing
    # and a rank-0 block in ldetS.
    added_rank <- 0
    
    if( !is.null(si$S[[1]]) && isTRUE(si$rank[1] > 0) ){
      out$S[[ length(out$S) + 1 ]] <- .pad_alpha(si$S[[1]], 0,  d1)
      out$rank   <- c(out$rank, si$rank[1])
      added_rank <- added_rank + si$rank[1]
    }
    if( !is.null(si$S[[2]]) && isTRUE(si$rank[2] > 0) ){
      out$S[[ length(out$S) + 1 ]] <- .pad_alpha(si$S[[2]], d1, d2)
      out$rank   <- c(out$rank, si$rank[2])
      added_rank <- added_rank + si$rank[2]
    }
    
    # Null space of the COMBINED alpha penalty: accumulated once, and only for
    # what was actually appended. Keying this off sum(si$rank) instead would
    # bump Mp even when no penalty made it into out$S, shifting REML by a
    # constant with no error raised.
    if( added_rank > 0 ){
      out$null.space.dim <- out$null.space.dim + (out$bs.dim - added_rank)
    }
  }
  
  # =========================================================================
  # 5. stack the two marginal single-index design matrices
  # =========================================================================
  if (is.list(out$xt$si$X)){
    out$xt$si$X <- cbind(out$xt$si$X[[1]], out$xt$si$X[[2]])
    # eff_inter_dlinear only receives this matrix (via .build_effects), not si,
    # so the block sizes have to travel with it as attributes. Elsewhere prefer
    # si$na1 / si$na2 -- attributes do not survive subsetting.
    attr(out$xt$si$X, "na1") <- d1
    attr(out$xt$si$X, "na2") <- d2
  }
  
  class(out) <- c("inter_dlinear", "nested")
  
  return( out )
}