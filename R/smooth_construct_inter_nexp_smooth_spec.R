#' Nested interactive adaptive exponential smoothing effect
#' 
#' @name smooth.construct.inter_nexp.smooth.spec
#' @rdname smooth.construct.inter_nexp.smooth.spec
#' @description Builds either
#' \itemize{
#'   \item{\code{s(exp(x), t)}}{ (class \code{"inter_nexp"}): \code{object$term[1]} is a matrix with
#'         columns named \code{"y"}, \code{"x"} (and optionally \code{"times"}), and \code{object$term[2]}
#'         is a plain (non-nested) vector \code{t}.}
#'   \item{\code{s(exp(x1), exp(x2))}}{ (class \code{"inter_ee"}): both \code{object$term[1]} and
#'         \code{object$term[2]} are matrices with columns named \code{"y"}, \code{"x"} (and optionally
#'         \code{"times"}). The two-sided case is detected automatically from the column names of
#'         \code{object$term[2]}.}
#' }
#' In the two-sided case, margin-specific initial values / penalties are passed via
#' \code{pord_1, S_1, alpha_1} (margin 1) and \code{pord_2, S_2, alpha_2} (margin 2), see \link{trans_inter_nexp}.
#' As in \link{trans_nexpsm}, each margin has a free scaling parameter (no design column) which is
#' followed by the smoothing-rate coefficients; an intercept for the rate is just one of the \code{"x"} columns.
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
  
  # Two-sided case s(exp(x1), exp(x2)): the second term is itself a matrix with "y"/"x" columns
  nms_t <- colnames(data[[term_t]])
  if( !is.null(nms_t) && all(c("y", "x") %in% nms_t) ){
    return( .smooth.construct.inter_ee(object, data, knots) )
  }
  
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
  nms <- colnames(X_mat)
  
  # Extract times and remove from X_mat to not interfere with x and W extraction
  times <- NULL
  tmp <- which(nms == "times")
  if (length(tmp)) {
    times <- X_mat[, tmp]
  }
  
  if (ncol(X_mat) < 2) {
    stop("X_mat must have at least 2 columns: x (data to smooth) and w variables.")
  }
  n <- nrow(X_mat)
  
  x_raw <- as.vector( t(X_mat[ , which(nms == "y")]) ) #data to be smoothed
  W_mat <- X_mat[ , which(nms == "x"), drop = FALSE] # model matrix used to model the exp smoothing rate/weight
  
  nrep <- ceiling( length(x_raw)/n )
  dXi <- ncol(W_mat)/nrep
  if(nrep > 1){
    tmp <- rep(1:dXi, nrep)
    W_mat <- apply(W_mat, 1, function(x) do.call("cbind", tapply(x, tmp, I)), simplify = FALSE)
    W_mat <- do.call("rbind", W_mat)
  }
  
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
    g <- expsmooth(y = x_raw, Xi = si$X, beta = rep(0, di-1), times = times)$d0
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
  si$times <- times
  
  # =========================================================================
  # 6. Build custom 2D B-spline basis (X_2D)
  # =========================================================================
  out <- .build_n_inter_bspline_basis(object = object, data = data, knots = knots, si = si, nested = c(TRUE, FALSE))
  
  # =========================================================================
  # 7. Assemble final block-diagonal penalty matrix
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


# ---------------------------------------------------------------------------
# Marginal initialisation of a single-level adaptive exponential smooth
# (same "y" / "x" / "times" conventions and nrep replication as
# smooth.construct.nexpsm.smooth.spec). Used by the two-sided branch below.
# ---------------------------------------------------------------------------
.init_marginal_nexp <- function(Xmat, S = NULL, pord = NULL, alpha = NULL, label = "term"){
  
  nms <- colnames(Xmat)
  if( is.null(nms) || sum(nms == "y") < 1 || sum(nms == "x") < 1 ){
    stop(label, " must be a matrix with columns named 'y' (data to smooth), 'x' (smoothing-rate covariates) ",
         "and, optionally, 'times'.")
  }
  
  n <- nrow(Xmat)
  y_raw <- as.vector( t(Xmat[ , which(nms == "y"), drop = FALSE]) )
  
  times <- NULL
  tmp <- which(nms == "times")
  if( length(tmp) ){ times <- Xmat[ , tmp] }
  
  W <- Xmat[ , which(nms == "x"), drop = FALSE]
  nrep <- ceiling( length(y_raw) / n )
  na <- ncol(W) / nrep
  if( na != round(na) ){
    stop(label, ": the number of 'x' columns (", ncol(W), ") must be a multiple of the number of 'y' columns (", nrep, ").")
  }
  if( nrep > 1 ){
    tmp <- rep(1:na, nrep)
    W <- apply(W, 1, function(z) do.call("cbind", tapply(z, tmp, I)), simplify = FALSE)
    W <- do.call("rbind", W)
  }
  
  if( !is.null(alpha) && length(alpha) != na ){
    stop("length of the initial smoothing-rate coefficients of ", label, " must be ", na, ".")
  }
  if( !is.null(S) && (nrow(S) != na || ncol(S) != na) ){
    stop("The penalty matrix of ", label, " must be ", na, "x", na, ".")
  }
  
  # Penalty on the rate coefficients (the scaling parameter is never penalised) + reparametrisation
  no_pen <- is.null(S) && is.null(pord)
  if( no_pen ){
    W_rot <- W; B <- diag(nrow = na); rank_S <- 0; S_out <- NULL
  } else {
    if( is.null(S) ){
      S <- .psp(d = na, ord = pord)
      rankS <- na - pord
    } else {
      rankS <- Matrix::rankMatrix(S)
    }
    dp <- gamFactory:::.diagPen(X = W, S = S, r = rankS)
    W_rot <- dp$X; B <- dp$B; rank_S <- dp$rank; S_out <- dp$S
  }
  
  alpha_w <- if( is.null(alpha) ) rep(0, na) else drop(solve(B, alpha))
  
  g <- drop( expsmooth(y = y_raw, Xi = W_rot, beta = alpha_w, times = times)$d0 )
  alpha_scale <- log(1 / sd(g))       # sd(inner index) = 1 at the start
  gm <- mean(g)
  
  list(z = exp(alpha_scale) * (g - gm), y_raw = y_raw, W = W_rot, B = B, S = S_out, rank = rank_S,
       times = times, na = na, alpha_scale = alpha_scale, alpha_w = alpha_w, xm = gm)
}


# ---------------------------------------------------------------------------
# Two-sided case: s(exp(x1), exp(x2))  (class "inter_ee")
# Inner parameter vector: c(alpha_scale_1, alpha_w_1, alpha_scale_2, alpha_w_2)
# ---------------------------------------------------------------------------
.smooth.construct.inter_ee <- function(object, data, knots){
  
  term_1 <- object$term[1]
  term_2 <- object$term[2]
  
  si <- object$xt$si
  if( is.null(si) ){ si <- object$xt$si <- list() }
  
  # Margin-specific inputs; pord/S/alpha (as in the one-sided case) are accepted for margin 1
  pick <- function(a, b){ if( !is.null(a) ) a else b }
  m1 <- .init_marginal_nexp(data[[term_1]], S = pick(si$S_1, si$S), pord = pick(si$pord_1, si$pord),
                            alpha = pick(si$alpha_1, si$alpha), label = paste0("term[1] (", term_1, ")"))
  m2 <- .init_marginal_nexp(data[[term_2]], S = si$S_2, pord = si$pord_2, alpha = si$alpha_2,
                            label = paste0("term[2] (", term_2, ")"))
  
  data[[term_1]] <- m1$z
  data[[term_2]] <- m2$z
  
  na1 <- m1$na; na2 <- m2$na
  d1  <- 1 + na1                       # size of the margin-1 block: scaling + rate coefficients
  d2  <- 1 + na2
  
  # Drop the user-supplied inputs so that only the internal (reparametrised) versions remain
  si$pord <- si$S <- si$alpha <- si$pord_1 <- si$pord_2 <- si$alpha_1 <- si$alpha_2 <- NULL
  
  si$na1 <- na1; si$na2 <- na2; si$d1 <- d1; si$d2 <- d2; si$na <- d1 + d2
  
  si$y_raw_1 <- m1$y_raw; si$W_1 <- m1$W; si$B_1 <- m1$B; si$S_1 <- m1$S; si$rank_1 <- m1$rank; si$times_1 <- m1$times
  si$y_raw_2 <- m2$y_raw; si$W_2 <- m2$W; si$B_2 <- m2$B; si$S_2 <- m2$S; si$rank_2 <- m2$rank; si$times_2 <- m2$times
  si$xm1 <- m1$xm; si$xm2 <- m2$xm
  si$alpha_scale_1 <- m1$alpha_scale; si$alpha_w_1 <- m1$alpha_w
  si$alpha_scale_2 <- m2$alpha_scale; si$alpha_w_2 <- m2$alpha_w
  
  si$alpha <- c(m1$alpha_scale, m1$alpha_w, m2$alpha_scale, m2$alpha_w)
  
  out <- .build_n_inter_bspline_basis(object = object, data = data, knots = knots,
                                       si = si, nested = c(TRUE, TRUE))
  
  # ---- penalties on the rate coefficients, padded into bs.dim space -------
  .pad_alpha <- function(Sk, off, dk){
    P <- matrix(0, out$bs.dim, out$bs.dim)
    ii <- off + seq_len(dk)
    P[ii, ii] <- Sk
    P
  }
  added_rank <- 0
  if( !is.null(si$S_1) && isTRUE(si$rank_1 > 0) ){
    out$S[[length(out$S) + 1]] <- .pad_alpha(si$S_1, 1, na1)              # skip alpha_scale_1
    out$rank <- c(out$rank, si$rank_1); added_rank <- added_rank + si$rank_1
  }
  if( !is.null(si$S_2) && isTRUE(si$rank_2 > 0) ){
    out$S[[length(out$S) + 1]] <- .pad_alpha(si$S_2, d1 + 1, na2)         # skip alpha_scale_2
    out$rank <- c(out$rank, si$rank_2); added_rank <- added_rank + si$rank_2
  }
  if( added_rank > 0 ){
    # each added penalty lives on alpha columns that were fully in the null space of S1/S2/S_inter
    out$null.space.dim <- out$null.space.dim - added_rank
  }
  
  class(out) <- c("inter_ee", "nested")
  return( out )
}
