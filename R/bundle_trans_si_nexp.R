#' @rdname trans_bundles
#' @export bundle_trans_si_nexp
#'
bundle_trans_si_nexp <- function(){
  list(
    bundle_nam = "bundle_trans_si_nexp",
    scale = FALSE,
    linear = FALSE,

    # parameters: c(alpha_nexp, alpha_si). z = g - mean(g), g = exponential smooth of the single index X_si alpha_si
    # with smoothing rate given by X_nexp alpha_nexp (see deriv_si_nexp)
    init = function(X, trans, label){
      positive_si <- isTRUE(trans$positive_si)

      # number of single index and of exp-smooth covariates
      n_si <- trans$n_si
      n_nexp <- trans$n_nexp
      if( any(colnames(X) == "times") ){
        if( is.null(n_si) ){ stop("`n_si` is NULL: the number of single-index covariates must be specified.", call. = FALSE) }
        if( is.null(n_nexp) ){ stop("`n_nexp` is NULL: the number of exponential-smoothing covariates must be specified.", call. = FALSE) }
      } else {
        if( !is.null(n_nexp) ){ n_si <- ncol(X) - n_nexp } else { n_nexp <- ncol(X) - n_si }
        if( n_si + n_nexp != ncol(X) ){
          stop(sprintf("Inconsistent number of variables: n_si + n_nexp = %d + %d != %d (ncol of Xall)", n_si, n_nexp, ncol(X)))
        }
      }
      ex <- .extract_si_nexp_margin(X, n_si, n_nexp)

      # optionally centre the single index covariates
      X_si <- ex$X_si
      xm_si <- rep(0, n_si)
      if( isTRUE(trans$center) ){
        X_si <- scale(X_si, scale = FALSE)
        xm_si <- attr(X_si, "scaled:center")
      }

      # penalties: reparametrise each block so that its penalty is diagonal
      no_pen_si <- is.null(trans$S_si) && is.null(trans$pord)
      if( no_pen_si ){
        B_si <- diag(n_si); S_si <- NULL; rank_si <- 0
      } else {
        S_si <- trans$S_si
        if( is.null(S_si) ){ S_si <- .psp(d = n_si, ord = trans$pord); rank_si <- n_si - trans$pord }
        else { rank_si <- Matrix::rankMatrix(S_si) }
        dp <- .diagPen(X = X_si, S = S_si, r = rank_si)
        X_si <- dp$X; B_si <- dp$B; S_si <- dp$S; rank_si <- dp$rank
      }
      X_nexp <- ex$X_nexp
      if( is.null(trans$S_nexp) ){
        B_nexp <- diag(n_nexp); S_nexp <- NULL; rank_nexp <- 0
      } else {
        dp <- .diagPen(X = X_nexp, S = trans$S_nexp, r = Matrix::rankMatrix(trans$S_nexp))
        X_nexp <- dp$X; B_nexp <- dp$B; S_nexp <- dp$S; rank_nexp <- dp$rank
      }

      # initial parameters (alpha_center is a fixed shift of alpha_si, as a0 of a si margin)
      alpha_center <- trans$alpha_center
      if( is.null(alpha_center) ){ alpha_center <- if( no_pen_si ) rep(0, n_si) else rep(1, n_si) }
      alpha_si <- trans$alpha_si
      if( is.null(alpha_si) ){
        alpha_si <- if( all(alpha_center == 0) ) rep(1, n_si) else if( positive_si ) rep(0.01, n_si) else rep(0, n_si)
      }
      alpha_si <- solve(B_si, alpha_si)
      alpha_center <- solve(B_si, alpha_center)
      if( positive_si ){ alpha_si <- log(pmax(alpha_si, 1e-8)) }       # positive weights: exp(alpha_si)
      alpha_nexp <- if( is.null(trans$alpha_nexp) ) rep(0, n_nexp) else solve(B_nexp) %*% trans$alpha_nexp

      g <- deriv_si_nexp(X_si = X_si, X_nexp = X_nexp, param = c(alpha_nexp, alpha_si), times = ex$times,
                         alpha_center = alpha_center, Z0 = trans$Z0, positive_si = positive_si)$d0

      # sd(z) = 1 at the start: g is linear in the single index, so rescale alpha_si and alpha_center
      alpha_si <- if( positive_si ) alpha_si - log(sd(g)) else alpha_si / sd(g)
      alpha_center <- alpha_center / sd(g)

      pen <- list()
      if( !no_pen_si ){ pen[[length(pen) + 1]] <- list(S = S_si, rank = rank_si, offset = n_nexp) }
      if( !is.null(S_nexp) ){ pen[[length(pen) + 1]] <- list(S = S_nexp, rank = rank_nexp, offset = 0) }

      list(z = (g - mean(g)) / sd(g), alpha = c(alpha_nexp, alpha_si), pen = pen,
           margin = list(X_si = X_si, X_nexp = X_nexp, times = ex$times, B_si = B_si, B_nexp = B_nexp,
                         n_si = n_si, n_nexp = n_nexp, center = isTRUE(trans$center), xm_si = xm_si,
                         alpha_center = alpha_center, Z0 = trans$Z0, positive_si = positive_si,
                         xm = mean(g) / sd(g)))
    },

    eval = function(mk, par, deriv = 0, xm = NULL){
      .nest_eval_centred(deriv_si_nexp(X_si = mk$X_si, X_nexp = mk$X_nexp, param = par, times = mk$times,
                                       deriv = deriv, alpha_center = mk$alpha_center, Z0 = mk$Z0,
                                       positive_si = mk$positive_si),
                         deriv = deriv, xm = xm)
    },

    newdata = function(mk, X){
      ex <- .extract_si_nexp_margin(X, mk$n_si, mk$n_nexp)
      X_si <- ex$X_si
      if( mk$center ){ X_si <- scale(X_si, center = mk$xm_si, scale = FALSE) }
      mk$X_si <- X_si %*% mk$B_si
      mk$X_nexp <- ex$X_nexp %*% mk$B_nexp
      mk$times <- ex$times
      mk
    }
  )
}

# Split the term of a si_nexp margin into the single index covariates X_si, the exp-smooth covariates X_nexp
# and the times. With a "times" column each row holds nrep sub-observations, [X_si x nrep | X_nexp x nrep | times]:
# they are put one per row, and the rows after the last fully finite one are dropped.
.extract_si_nexp_margin <- function(X, n_si, n_nexp){

  i_t <- which(colnames(X) == "times")
  if( !length(i_t) ){
    return( list(X_si = X[ , 1:n_si, drop = FALSE], X_nexp = X[ , n_si + 1:n_nexp, drop = FALSE], times = NULL) )
  }

  times <- X[ , i_t]
  X <- X[ , -i_t, drop = FALSE]
  nrep <- ceiling(ncol(X) / (n_si + n_nexp))
  stopifnot(ncol(X) == nrep * (n_si + n_nexp))

  # row d of a block -> nrep rows of p columns
  long <- function(B, p) matrix(as.vector(t(B)), ncol = p, byrow = TRUE)
  X_si   <- long(X[ , 1:(nrep * n_si), drop = FALSE], n_si)
  X_nexp <- long(X[ , nrep * n_si + 1:(nrep * n_nexp), drop = FALSE], n_nexp)

  ok <- rowSums(is.finite(X_si)) == n_si & rowSums(is.finite(X_nexp)) == n_nexp
  last_good <- max(which(ok))

  list(X_si = X_si[1:last_good, , drop = FALSE], X_nexp = X_nexp[1:last_good, , drop = FALSE], times = times)
}
