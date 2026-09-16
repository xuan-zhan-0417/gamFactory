#'
#' Build si(x) / exp(x) interaction effect
#'
#' @param Xi_1 margin-1 (single index) design matrix, already centred and
#'             rotated (as produced by \link{smooth.construct.inter_le.smooth.spec}).
#' @param a0_1 fixed initialisation shift for margin 1 (length \code{ncol(Xi_1)}).
#' @param y2 data to be exponentially smoothed for margin 2.
#' @param Xi_2 margin-2 design matrix (rate covariates), already rotated.
#' @param basis function which takes \code{z1, z2} as input and returns the
#'              model matrix and its derivatives.
#' @param times optional vector of times, as in \link{trans_nexpsm}.
#' @name eff_inter_le
#' @rdname eff_inter_le
#' @export eff_inter_le
#'
eff_inter_le <- function(Xi_1, a0_1 = NULL, y2, Xi_2, basis, times = NULL){

  force(Xi_1); force(a0_1); force(y2); force(Xi_2); force(basis); force(times);

  na1 <- ncol(Xi_1)
  na2 <- ncol(Xi_2)
  na  <- na1 + 1 + na2

  if( is.null(a0_1) ){ a0_1 <- numeric(na1) }

  idx1  <- 1:na1
  idx_s <- na1 + 1
  idx2  <- (na1 + 2):na

  eval <- function(param, deriv = 0){

    n <- nrow(Xi_1)
    beta <- param[ -(1:na) ]

    alpha_1     <- param[idx1]
    alpha_scale <- param[idx_s]
    alpha_2     <- param[idx2]
    a0s <- exp(alpha_scale)

    # margin 1: si(x1) -- linear in alpha_1
    z1 <- drop( Xi_1 %*% (alpha_1 + a0_1) )

    # margin 2: exp(x) -- single-level adaptive exponential smoothing
    inner <- expsmooth(y = y2, Xi = Xi_2, beta = alpha_2, times = times, deriv = deriv)
    g2raw <- drop(inner$d0)
    z2    <- a0s * (g2raw - mean(g2raw))

    store <- basis$evalX(z1 = z1, z2 = z2, deriv = deriv)
    store$Xi_1 <- Xi_1
    store$Xi_2 <- Xi_2
    store$g    <- cbind(z1 = z1, z2 = z2)

    if( deriv >= 1 ){

      store$f1 <- list(
        f1_1 = drop( store$X1$dz1 %*% beta ),
        f1_2 = drop( store$X1$dz2 %*% beta )
      )

      # dz2 / d(alpha_scale, alpha_2): the exp(alpha_scale) trick makes the
      # cross/self partial w.r.t. alpha_scale collapse to z2 itself (see
      # DllkDbeta.inter_le / .get_eff_eval_general for the same trick).
      g1_2 <- cbind(z2, a0s * (t(t(inner$d1) - colMeans(inner$d1))))

      store$g1 <- list(
        g1_1 = Xi_1,  # dz1 / dalpha_1
        g1_2 = g1_2   # dz2 / d(alpha_scale, alpha_2)
      )

      if( deriv >= 2 ){

        store$f2 <- list(
          f2_11 = drop( store$X2$dz1_z1 %*% beta ),
          f2_22 = drop( store$X2$dz2_z2 %*% beta ),
          f2_12 = drop( store$X2$dz1_z2 %*% beta )
        )

        g2_2 <- cbind(g1_2, a0s * (t(t(inner$d2) - colMeans(inner$d2))))

        store$g2 <- list(
          g2_1 = matrix(0, n, na1 * (na1 + 1) / 2),  # z1 linear in alpha_1
          g2_2 = g2_2
        )
      }
    }

    o <- eff_inter_le(Xi_1 = Xi_1, a0_1 = a0_1, y2 = y2, Xi_2 = Xi_2, basis = basis, times = times)
    o$f     <- drop( store$X0 %*% beta )
    o$param <- param
    o$a0_1  <- a0_1
    o$na    <- na; o$na1 <- na1; o$na2 <- na2
    o$store <- store
    o$deriv <- deriv

    return( o )

  }

  out <- structure(list("eval" = eval), class = c("inter_le", "nested"))

  return( out )

}
