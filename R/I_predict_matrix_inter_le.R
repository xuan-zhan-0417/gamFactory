#' Predict using the si(x) / exp(x) interaction effect (ANOVA tensor product)
#'
#' @param object smooth object of class "inter_le".
#' @param data data frame / list holding the two inner term matrices.
#' @param get.xa if TRUE, return the inner indices and their Jacobians
#'               instead of the model matrix.
#' @noRd
.predict.matrix.inter_le <- function(object, data, get.xa = FALSE){

  term_1 <- object$term[1]
  term_2 <- object$term[2]

  si  <- object$xt$si
  na1 <- si$na1
  na2 <- si$na2
  di  <- si$na

  if( is.null(si$xm2) ){ si$xm2 <- 0 }

  # ---- margin 1: si(x1) -- same reparametrisation as .predict.matrix.si ----
  X1 <- as.matrix(data[[term_1]])
  X1 <- t(t(X1) - si$xm1) %*% si$B_1
  z1 <- drop( X1 %*% (si$alpha_1 + si$a0_1) )

  # ---- margin 2: exp(x) -- same extraction as .predict.matrix.nexpsm -------
  X2   <- data[[term_2]]
  n2   <- nrow(X2)
  nms2 <- colnames(X2)

  y_raw <- as.vector( t(X2[ , which(nms2 == "y"), drop = FALSE]) )
  times <- NULL
  tmp <- which(nms2 == "times")
  if( length(tmp) ){ times <- X2[ , tmp] }

  W2   <- X2[ , which(nms2 == "x"), drop = FALSE]
  nrep <- ceiling( length(y_raw) / n2 )
  if( nrep > 1 ){
    tmp <- rep(1:na2, nrep)
    W2 <- apply(W2, 1, function(z) do.call("cbind", tapply(z, tmp, I)), simplify = FALSE)
    W2 <- do.call("rbind", W2)
  }
  W2 <- W2 %*% si$B_2

  xsm_list <- expsmooth(y = y_raw, Xi = W2, beta = si$alpha_2, times = times, deriv = get.xa)
  g2_unscaled <- drop(xsm_list$d0) - si$xm2
  z2 <- exp(si$alpha_scale) * g2_unscaled

  if( get.xa ){
    return( list(xa    = cbind(z1 = z1, z2 = z2),
                 xa_da = list(X1, exp(si$alpha_scale) * cbind(g2_unscaled, xsm_list$d1)),
                 z1 = z1, z2 = z2,
                 Xi_1 = X1, Xi_2 = W2,
                 na1 = na1, na2 = na2) )
  }

  X0 <- object$xt$basis$evalX(z1 = z1, z2 = z2, deriv = 0)$X0

  # [ 0_di , X_2D ] : predict.gam multiplies the leading block by alpha -> 0
  Xtot <- cbind(matrix(0, nrow(X0), di), X0)

  # margin-1 value is already "final" (si-style); margin-2 value is stored
  # *unscaled* (nexp-style, pre exp(alpha_scale)) so that postproc_gam_nl can
  # refresh xm2 to mean(g) under the final fitted alpha_2, exactly as for
  # standalone nexpsm/inter_nexp.
  attr(Xtot, "inner_linpred_unscaled") <- cbind(z1 = z1, z2_unscaled = g2_unscaled)

  return( Xtot )
}
