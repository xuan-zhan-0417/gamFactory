# Sibling of R/I_init_marginal_si.R: initialisation of one exp margin (used by .inter_init_exp).

# ---------------------------------------------------------------------------
# One adaptive exponential smooth margin: same "y" / "x" / "times" conventions and nrep
# replication as smooth.construct.nexpsm.smooth.spec. alpha_scale (a free scalar, no
# design column) is followed by the smoothing-rate coefficients alpha_w.
# ---------------------------------------------------------------------------
.init_marginal_nexp <- function(Xmat, S = NULL, pord = NULL, alpha = NULL, alpha_scale = NULL, label = "term"){

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
  if( is.null(alpha_scale) ){ alpha_scale <- log(1 / sd(g)) }   # sd(inner index) = 1 at the start
  gm <- mean(g)

  list(z = exp(alpha_scale) * (g - gm), y_raw = y_raw, W = W_rot, W_raw = W, B = B, S = S_out, rank = rank_S,
       times = times, na = na, alpha_scale = alpha_scale, alpha_w = alpha_w, xm = gm, no_pen = no_pen)
}
