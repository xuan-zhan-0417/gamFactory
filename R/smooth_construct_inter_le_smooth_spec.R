#' Nested single-index / adaptive exponential-smoothing interaction
#'
#' @name smooth.construct.inter_le.smooth.spec
#' @rdname smooth.construct.inter_le.smooth.spec
#' @description Builds the effect \code{s(si(x), exp(x))}: a tensor-product interaction
#'              between a nested single index (margin 1) and a nested, single-level
#'              adaptive exponential smooth (margin 2). Margin 2 is "single-level"
#'              in the sense that, unlike \link{trans_linear_nexpsm}, the covariate being
#'              smoothed is not itself the output of a nested single index: it is a plain
#'              data column, and only the smoothing rate is modelled via covariates.
#' @details \code{object$term[1]} must be a matrix with \code{na_1} columns (the data used
#'          to build the margin-1 single index, as in \link{trans_linear}). \code{object$term[2]}
#'          must be a matrix with columns named
#'          \itemize{
#'            \item{\code{"y"}}{ exactly one column: the data to be exponentially smoothed.}
#'            \item{\code{"intercept"}}{ exactly one column, entirely equal to 1: it multiplies
#'                  the scaling parameter \code{alpha_scale} (see \link{trans_inter_le}).}
#'            \item{\code{"x"}}{ \code{na_2} columns used to model the exponential smoothing
#'                  rate, via coefficients \code{alpha_2}.}
#'          }
#'          An optional \code{"times"} column can be included, as in \link{trans_nexpsm}, when
#'          the response is observed at a different frequency than the covariates driving the
#'          smoothing rate.
#' @importFrom MASS Null
#' @importFrom Matrix rankMatrix
#' @export
#'
smooth.construct.inter_le.smooth.spec <- function(object, data, knots){

  if(length(object$term) != 2) {
    stop("The smooth effect must contain exactly two terms: si(x) and exp(x).")
  }
  term_1 <- object$term[1]
  term_2 <- object$term[2]

  si <- object$xt$si
  if( is.null(si) ){ si <- object$xt$si <- list() }

  # =========================================================================
  # Margin 1: si(x) -- standard nested single index (as in inter_linear)
  # =========================================================================
  X1  <- as.matrix(data[[term_1]])
  na1 <- ncol(X1)

  if( !is.null(si$alpha_1) && length(si$alpha_1) != na1 ){
    stop("length(si$alpha_1) must equal ncol(", term_1, ") = ", na1, ".")
  }
  if( !is.null(si$a0_1) && length(si$a0_1) != na1 ){
    stop("length(si$a0_1) must equal ncol(", term_1, ") = ", na1, ".")
  }
  if( !is.null(si$S_1) && (nrow(si$S_1) != na1 || ncol(si$S_1) != na1) ){
    stop("si$S_1 must be a ", na1, "x", na1, " matrix.")
  }

  res1 <- .init_marginal_si(Xi = X1, S = si$S_1, pord = si$pord_1,
                             a0 = si$a0_1, alpha = si$alpha_1)
  data[[term_1]] <- res1$ax

  # =========================================================================
  # Margin 2: exp(x) -- single-level adaptive exponential smoothing
  # =========================================================================
  X2 <- data[[term_2]]
  if( is.null(X2) || is.null(colnames(X2)) ){
    stop("term[2] must be a matrix with columns named 'y', 'intercept' and 'x'.")
  }
  nms2 <- colnames(X2)

  if( sum(nms2 == "y") != 1 ){
    stop("term[2] must contain exactly one column named 'y' (the data to be smoothed).")
  }
  if( sum(nms2 == "intercept") != 1 ){
    stop("term[2] must contain exactly one column named 'intercept' (a column of 1s).")
  }
  if( sum(nms2 == "x") < 1 ){
    stop("term[2] must contain at least one column named 'x' (covariates for the smoothing rate).")
  }

  y_raw  <- as.vector(X2[ , which(nms2 == "y")])
  interc <- X2[ , which(nms2 == "intercept")]
  if( !all(interc == 1) ){
    stop("term[2] column 'intercept' must be a column of 1s.")
  }
  W2  <- X2[ , which(nms2 == "x"), drop = FALSE]
  na2 <- ncol(W2)

  times <- NULL
  tmp <- which(nms2 == "times")
  if( length(tmp) ){ times <- X2[ , tmp] }

  if( !is.null(si$alpha_2) && length(si$alpha_2) != na2 ){
    stop("length(si$alpha_2) must equal the number of 'x' columns in term[2] (", na2, ").")
  }
  if( !is.null(si$S_2) && (nrow(si$S_2) != na2 || ncol(si$S_2) != na2) ){
    stop("si$S_2 must be a ", na2, "x", na2, " matrix.")
  }

  # ---- penalty & diagonalisation on alpha_2 (the scaling parameter is never penalised) ----
  S2 <- si$S_2
  no_pen_2 <- is.null(S2) && is.null(si$pord_2)

  if( no_pen_2 ){
    W2_rot <- W2
    B2     <- diag(nrow = na2)
    rank_2 <- 0
    S2_out <- NULL
  } else {
    if( is.null(S2) ){
      S2     <- .psp(d = na2, ord = si$pord_2)
      rankS2 <- na2 - si$pord_2
    } else {
      rankS2 <- Matrix::rankMatrix(S2)
    }
    diag_pen_2 <- gamFactory:::.diagPen(X = W2, S = S2, r = rankS2)
    rank_2 <- diag_pen_2$rank
    W2_rot <- diag_pen_2$X
    B2     <- diag_pen_2$B
    S2_out <- diag_pen_2$S
  }

  # ---- initialise alpha_2 (smoothing-rate coefficients) and alpha_scale ----
  alpha_2 <- si$alpha_2
  if( is.null(alpha_2) ){
    alpha_2 <- rep(0, na2)
  } else {
    alpha_2 <- drop(solve(B2, alpha_2))
  }
  g <- expsmooth(y = y_raw, Xi = W2_rot, beta = alpha_2, times = times)$d0

  alpha_scale <- si$alpha_scale
  if( is.null(alpha_scale) ){
    # alpha_scale s.t. sd(inner_lin_pred) = 1 (target variance), as in trans_nexpsm
    alpha_scale <- log(1 / sd(g))
  }

  gm <- mean(g)
  data[[term_2]] <- exp(alpha_scale) * (g - gm)

  # =========================================================================
  # Store everything needed downstream (prediction / jacobians / penalties)
  # =========================================================================
  si$na1 <- na1; si$na2 <- na2; si$na <- na1 + 1 + na2

  si$X_1 <- res1$X; si$B_1 <- res1$B; si$S_1 <- res1$S; si$rank_1 <- res1$rank
  si$xm1 <- res1$xm; si$alpha_1 <- res1$alpha; si$a0_1 <- res1$a0

  si$y_raw  <- y_raw; si$W_2 <- W2_rot; si$B_2 <- B2; si$S_2 <- S2_out
  si$rank_2 <- rank_2; si$times <- times
  si$xm2 <- gm; si$alpha_scale <- alpha_scale; si$alpha_2 <- alpha_2

  # Full inner parameter vector: c(alpha_1, alpha_scale, alpha_2)
  si$alpha <- c(res1$alpha, alpha_scale, alpha_2)

  out <- .build_n_inter_bspline_basis(object = object, data = data, knots = knots,
                                       si = si, nested = c(TRUE, TRUE))

  # ---- penalty: pad each margin's alpha penalty into bs.dim space ----------
  .pad_alpha <- function(Sk, off, dk){
    P <- matrix(0, out$bs.dim, out$bs.dim)
    ii <- off + seq_len(dk)
    P[ii, ii] <- Sk
    P
  }
  added_rank <- 0
  if( !is.null(si$S_1) && isTRUE(si$rank_1 > 0) ){
    out$S[[length(out$S) + 1]] <- .pad_alpha(si$S_1, 0, na1)
    out$rank <- c(out$rank, si$rank_1); added_rank <- added_rank + si$rank_1
  }
  if( !is.null(si$S_2) && isTRUE(si$rank_2 > 0) ){
    # offset by na1 + 1 to skip margin-1 alphas and the (unpenalised) scaling parameter
    out$S[[length(out$S) + 1]] <- .pad_alpha(si$S_2, na1 + 1, na2)
    out$rank <- c(out$rank, si$rank_2); added_rank <- added_rank + si$rank_2
  }
  if( added_rank > 0 ){
    # Each added penalty acts on a block of alpha coefficients that was
    # previously entirely in the null space of S1/S2/S_inter (those are zero
    # on the leading "di" alpha columns), so it removes exactly "added_rank"
    # dimensions from the null space -- not "bs.dim - added_rank".
    out$null.space.dim <- out$null.space.dim - added_rank
  }

  class(out) <- c("inter_le", "nested")
  return( out )
}
