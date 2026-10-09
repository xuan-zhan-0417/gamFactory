#' @rdname trans_bundles
#' @export bundle_trans_exp
#'
bundle_trans_exp <- function(){
  list(
    bundle_nam = "bundle_trans_exp",
    scale = TRUE,
    linear = FALSE,

    init = function(X, trans, label){
      res <- .init_marginal_nexp(X, S = trans$S, pord = trans$pord, alpha = trans$alpha,
                                 alpha_scale = trans$alpha_scale, label = label)
      pen <- if( !is.null(res$S) ) list(list(S = res$S, rank = res$rank, offset = 1))   # the scale is not penalised
      list(z = res$z, alpha = c(res$alpha_scale, res$alpha_w), pen = pen, kex = .kex_exp(res$y_raw, nrow(res$W), res$times),
           margin = list(y = res$y_raw, W = res$W, times = res$times, B = res$B, xm = res$xm))
    },

    eval = function(mk, par, deriv = 0, xm = NULL){
      .nest_eval_scaled(expsmooth(y = mk$y, Xi = mk$W, beta = par[-1], times = mk$times, deriv = deriv),
                         scale = par[1], deriv = deriv, xm = xm)
    },

    newdata = function(mk, X){
      ex <- .extract_exp_margin(X, na = ncol(mk$B), B = mk$B)
      mk$y <- ex$y; mk$W <- ex$W; mk$times <- ex$times
      mk
    }
  )
}

# Extract the "y" / "x" / "times" pieces of one exp(x) margin from its term matrix
# (same conventions and nrep replication as .predict.matrix.nexpsm) and rotate the
# rate design with the margin's reparametrisation matrix B.
.extract_exp_margin <- function(X, na, B){

  nms <- colnames(X)
  n   <- nrow(X)

  y <- as.vector( t(X[ , which(nms == "y"), drop = FALSE]) )
  times <- NULL
  tmp <- which(nms == "times")
  if( length(tmp) ){ times <- X[ , tmp] }

  W <- X[ , which(nms == "x"), drop = FALSE]
  nrep <- ceiling( length(y) / n )
  if( nrep > 1 ){
    tmp <- rep(1:na, nrep)
    W <- apply(W, 1, function(z) do.call("cbind", tapply(z, tmp, I)), simplify = FALSE)
    W <- do.call("rbind", W)
  }

  list(y = y, W = W %*% B, times = times)
}

# Range of the outer knots of an exponential smooth: the smallest and largest value of the standardised smooth
# over constant smoothing rates (nr = number of rows of the rate design)
.kex_exp <- function(y, nr, times){
  one <- matrix(1, nrow = nr)
  z_std <- function(a){
    g <- expsmooth(y = y, Xi = one, beta = a, times = times)$d0
    (g - mean(g)) / .sd_n(g)
  }
  rates <- qlogis(c(1e-4, 1 - 1e-4))
  lo <- optimize(function(a) -min(z_std(a))^2, interval = rates)$objective
  hi <- optimize(function(a) -max(z_std(a))^2, interval = rates)$objective
  1.1 * c(-sqrt(-lo), sqrt(-hi))
}
