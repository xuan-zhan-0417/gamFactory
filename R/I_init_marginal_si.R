
# used by bundle_trans_si and bundle_trans_si_pos (R/bundle_trans_si.R, R/bundle_trans_si_pos.R)
# build the marginal si for each variable
# including centering, diagonalization, and initialization of alpha and a0
# positive = TRUE: the weights are exp(alpha + a0), so alpha is on the log scale

.init_marginal_si <- function(Xi, S = NULL, pord = NULL, a0 = NULL, alpha = NULL, positive = FALSE) {
  # (a) Centering
  Xi_centered <- scale(Xi, scale = FALSE)
  xm <- attr(Xi_centered, "scaled:center")
  di <- ncol(Xi_centered)
  
  # (b) Penalty & Diagonalization (.diagPen)
  no_pen <- is.null(S) && is.null(pord)
  
  if (no_pen) {
    X <- Xi_centered
    B <- diag(nrow = di)
    rank_S <- 0
    S_out <- NULL
  } else {
    if (is.null(S)) {
      S <- .psp(d = di, ord = pord)
      rankS <- di - pord
    } else {
      rankS <- Matrix::rankMatrix(S)
    }
    diag_pen <- gamFactory:::.diagPen(X = Xi_centered, S = S, r = rankS)
    rank_S   <- diag_pen$rank
    X        <- diag_pen$X
    B        <- diag_pen$B
    S_out    <- diag_pen$S
  }
  
  # (c) Initialize alpha & a0
  if (is.null(a0)) {
    a0 <- if (no_pen) rep(0, di) else rep(1, di)
  }
  
  if (is.null(alpha)) {
    alpha <- if (all(a0 == 0)) rep(1, di) else rep(0, di)
  }
  
  # Basis transformation using B matrix
  alpha <- solve(B, alpha)
  a0    <- solve(B, a0)
  
  # Impose variance constraint: var(z) = 1, dividing by n as the penalty on the variance does
  if (positive) {
    alpha <- alpha - log(.sd_n(X %*% exp(alpha + a0)))
    ax    <- drop(X %*% exp(alpha + a0))
  } else {
    scale_factor <- .sd_n(X %*% (alpha + a0))
    alpha <- alpha / scale_factor
    a0    <- a0 / scale_factor
    ax    <- drop(X %*% (alpha + a0))
  }
  
  list(
    xm = xm, X = X, B = B, S = S_out, rank = rank_S,
    a0 = a0, alpha = alpha, ax = ax
  )
}