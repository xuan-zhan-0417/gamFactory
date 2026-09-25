# Extract the "y" / "x" / "times" pieces of one exp(x) margin from its term matrix
# (same conventions and nrep replication as .predict.matrix.nexpsm) and rotate the
# rate design with the margin's reparametrisation matrix B.
#' @noRd
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

#' Predict using the exp(x1) / exp(x2) interaction effect (ANOVA tensor product)
#'
#' @param object smooth object of class "inter_ee".
#' @param data data frame / list holding the two inner term matrices.
#' @param get.xa if TRUE, return the inner indices and their Jacobians
#'               instead of the model matrix.
#' @noRd
.predict.matrix.inter_ee <- function(object, data, get.xa = FALSE){

  term <- object$term
  si   <- object$xt$si
  di   <- si$na

  z <- zu <- xa_da <- vector("list", 2)
  Xi <- vector("list", 2)

  for(k in 1:2){
    ex <- .extract_exp_margin(data[[ term[k] ]], na = si[[paste0("na", k)]], B = si[[paste0("B_", k)]])
    a_scale <- si[[paste0("alpha_scale_", k)]]
    xm <- si[[paste0("xm", k)]]; if( is.null(xm) ){ xm <- 0 }

    sm <- expsmooth(y = ex$y, Xi = ex$W, beta = si[[paste0("alpha_w_", k)]], times = ex$times, deriv = get.xa)
    zu[[k]] <- drop(sm$d0) - xm
    z[[k]]  <- exp(a_scale) * zu[[k]]
    Xi[[k]] <- ex$W
    if( get.xa ){ xa_da[[k]] <- cbind(z[[k]], exp(a_scale) * sm$d1) }
  }

  if( get.xa ){
    return( list(xa    = cbind(z1 = z[[1]], z2 = z[[2]]),
                 xa_da = xa_da,
                 z1 = z[[1]], z2 = z[[2]],
                 Xi_1 = Xi[[1]], Xi_2 = Xi[[2]]) )
  }

  X0 <- object$xt$basis$evalX(z1 = z[[1]], z2 = z[[2]], deriv = 0)$X0

  # [ 0_di , X_2D ] : predict.gam multiplies the leading block by alpha -> 0
  Xtot <- cbind(matrix(0, nrow(X0), di), X0)

  # margin values stored *unscaled* (pre exp(alpha_scale)) so that postproc_gam_nl can
  # refresh xm1/xm2 to mean(g) under the final fitted parameters
  attr(Xtot, "inner_linpred_unscaled") <- cbind(z1_unscaled = zu[[1]], z2_unscaled = zu[[2]])

  return( Xtot )
}
