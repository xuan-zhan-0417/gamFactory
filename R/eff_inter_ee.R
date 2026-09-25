#'
#' Build exp(x1) / exp(x2) interaction effect
#'
#' @param y1,y2 data to be exponentially smoothed for margin 1 / margin 2.
#' @param Xi_1,Xi_2 smoothing-rate design matrices of margin 1 / margin 2 (already rotated).
#' @param basis function which takes \code{z1, z2} as input and returns the
#'              model matrix and its derivatives.
#' @param times_1,times_2 optional vectors of times, as in \link{trans_nexpsm}.
#' @details The inner parameter vector is \code{c(alpha_scale_1, alpha_w_1, alpha_scale_2, alpha_w_2)},
#'          and \code{z_k = exp(alpha_scale_k) * (g_k - mean(g_k))} with \code{g_k} the exponential smooth of margin \code{k}.
#' @name eff_inter_ee
#' @rdname eff_inter_ee
#' @export eff_inter_ee
#'
eff_inter_ee <- function(y1, Xi_1, y2, Xi_2, basis, times_1 = NULL, times_2 = NULL){

  force(y1); force(Xi_1); force(y2); force(Xi_2); force(basis); force(times_1); force(times_2)

  na1 <- ncol(Xi_1)
  na2 <- ncol(Xi_2)
  d1  <- 1 + na1
  d2  <- 1 + na2
  na  <- d1 + d2

  idx1 <- 1:d1
  idx2 <- (d1 + 1):na

  # One exp(x) margin: z = exp(par[1]) * (g - mean(g)), and its derivatives w.r.t. par = (scale, rate coefs).
  # Because of the exp(scale) factor, every derivative involving the scale collapses onto a lower-order
  # one (d z/d scale = z, d2 z/d scale d par_j = d z / d par_j), so the first "d" columns of g2 are g1.
  margin <- function(y, Xi, times, par, deriv){
    a0 <- exp(par[1])
    inner <- expsmooth(y = y, Xi = Xi, beta = par[-1], times = times, deriv = min(deriv, 2))
    g <- drop(inner$d0)
    z <- a0 * (g - mean(g))
    out <- list(z = z)
    if( deriv >= 1 ){
      out$g1 <- cbind(z, a0 * t(t(inner$d1) - colMeans(inner$d1)))
      if( deriv >= 2 ){
        out$g2 <- cbind(out$g1, a0 * t(t(inner$d2) - colMeans(inner$d2)))
      }
    }
    out
  }

  eval <- function(param, deriv = 0){

    beta <- param[ -(1:na) ]

    m1 <- margin(y1, Xi_1, times_1, param[idx1], deriv)
    m2 <- margin(y2, Xi_2, times_2, param[idx2], deriv)

    store <- basis$evalX(z1 = m1$z, z2 = m2$z, deriv = deriv)
    store$g <- cbind(z1 = m1$z, z2 = m2$z)

    if( deriv >= 1 ){

      store$f1 <- list(
        f1_1 = drop( store$X1$dz1 %*% beta ),
        f1_2 = drop( store$X1$dz2 %*% beta )
      )
      store$g1 <- list(g1_1 = m1$g1, g1_2 = m2$g1)

      if( deriv >= 2 ){

        store$f2 <- list(
          f2_11 = drop( store$X2$dz1_z1 %*% beta ),
          f2_22 = drop( store$X2$dz2_z2 %*% beta ),
          f2_12 = drop( store$X2$dz1_z2 %*% beta )
        )
        store$g2 <- list(g2_1 = m1$g2, g2_2 = m2$g2)
      }
    }

    o <- eff_inter_ee(y1 = y1, Xi_1 = Xi_1, y2 = y2, Xi_2 = Xi_2, basis = basis,
                      times_1 = times_1, times_2 = times_2)
    o$f     <- drop( store$X0 %*% beta )
    o$param <- param
    o$na    <- na; o$na1 <- na1; o$na2 <- na2; o$d1 <- d1; o$d2 <- d2
    o$store <- store
    o$deriv <- deriv

    return( o )

  }

  out <- structure(list("eval" = eval), class = c("inter_ee", "nested"))

  return( out )

}
