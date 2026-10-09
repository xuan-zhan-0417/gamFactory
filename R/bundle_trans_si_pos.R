#' @rdname trans_bundles
#' @export bundle_trans_si_pos
#'
bundle_trans_si_pos <- function(){

  # same data and prediction as a single index, but with positive weights exp(alpha + a0)
  b <- bundle_trans_si()
  b$bundle_nam <- "bundle_trans_si_pos"
  b$linear <- FALSE

  b$init <- function(X, trans, label){
    res <- .init_marginal_si(Xi = as.matrix(X), S = trans$S, pord = trans$pord, a0 = trans$a0, alpha = trans$alpha,
                             positive = TRUE)
    pen <- if( res$rank > 0 ) list(list(S = res$S, rank = res$rank, offset = 0))
    list(z = res$ax, alpha = res$alpha, pen = pen, kex = .kex_si(res$X),
         margin = list(X = res$X, B = res$B, xm = res$xm, a0 = res$a0))
  }

  # z = X exp(alpha + a0). dz/dalpha_j = X_j exp(alpha_j + a0_j), and so are all the derivatives w.r.t. alpha_j
  # only (j, j) and (j, j, j); the mixed ones are zero
  b$eval <- function(mk, par, deriv = 0, xm = NULL){
    p  <- length(par)
    g1 <- t(t(mk$X) * exp(par + mk$a0))
    out <- list(z = rowSums(g1), xm = mk$xm, g1 = g1)
    ind <- trind.generator(p)
    if( deriv >= 2 ){
      out$g2 <- matrix(0, nrow(g1), p * (p + 1) / 2)
      out$g2[ , ind$i2[cbind(1:p, 1:p)]] <- g1
    }
    if( deriv >= 3 ){
      out$g3 <- matrix(0, nrow(g1), p * (p + 1) * (p + 2) / 6)
      out$g3[ , ind$i3[cbind(1:p, 1:p, 1:p)]] <- g1
    }
    out
  }

  b
}
