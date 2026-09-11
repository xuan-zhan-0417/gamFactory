#' Predict using 2D inter linear effects (ANOVA Tensor Product)
#'
#' Handles both sub-cases produced by \code{smooth.construct.inter_linear.smooth.spec}:
#' margin 2 nested (\code{s(si(x1), si(x2))}, \code{si$na2 > 0}) and margin 2
#' plain (\code{s(si(x1), t)}, \code{si$na2 == 0}).
#'
#' @param object smooth object of class "inter_linear".
#' @param data data frame / list holding the inner covariate matrix (matrices).
#' @param get.xa if TRUE, return the inner indices and their Jacobians
#'               instead of the model matrix.
#' @noRd
.predict.matrix.inter_linear <- function(object, data, get.xa = FALSE){

  term_1 <- object$term[1]
  term_2 <- object$term[2]

  si <- object$xt$si
  na1 <- si$na1
  na2 <- si$na2
  di  <- si$na
  nested_2 <- isTRUE(na2 > 0)

  if( nested_2 ){

    if( length(si$alpha) != di || length(si$a0) != di || length(si$xm) != di ){
      stop("si$alpha / si$a0 / si$xm must all have length na1 + na2 = ", di, ".")
    }
    if( !is.list(si$B) || length(si$B) != 2 ||
        any(vapply(si$B, is.null, logical(1))) ){
      stop("si$B must be a list of the two marginal reparametrisation matrices ",
           "(use the identity when a margin is unpenalised, never NULL).")
    }

    ik <- list(1:na1, (na1 + 1):di)     # alpha slices
    tm <- c(term_1, term_2)

    Xi <- vector("list", 2)             # marginal inner matrices = dz_k / dalpha_k
    za <- vector("list", 2)             # marginal inner indices z_k

    for(k in 1:2){

      Xk <- as.matrix( data[[ tm[k] ]] )
      if( ncol(Xk) != length(ik[[k]]) ){
        stop("Term '", tm[k], "' has ", ncol(Xk), " columns but na", k,
             " = ", length(ik[[k]]), ".")
      }

      # same reparametrisation as .predict.matrix.si: centre with xm, then rotate with B
      Xk <- t(t(Xk) - si$xm[ ik[[k]] ]) %*% si$B[[k]]

      Xi[[k]] <- Xk
      za[[k]] <- drop( Xk %*% (si$alpha[ ik[[k]] ] + si$a0[ ik[[k]] ]) )
    }

    if( get.xa ){
      return( list(xa    = cbind(z1 = za[[1]], z2 = za[[2]]),  # n x 2
                   xa_da = Xi,                                  # list(n x na1, n x na2)
                   z1 = za[[1]], z2 = za[[2]],
                   Xi_1 = Xi[[1]], Xi_2 = Xi[[2]],
                   na1 = na1, na2 = na2) )
    }

    X0 <- object$xt$basis$evalX(z1 = za[[1]], z2 = za[[2]], deriv = 0)$X0

    # [ 0_di , X_2D ] : predict.gam multiplies the leading block by alpha -> 0
    Xtot <- cbind(matrix(0, nrow(X0), di), X0)

    if( !is.null(object$bs.dim) && ncol(Xtot) != object$bs.dim ){
      stop("Predict matrix has ", ncol(Xtot), " columns but bs.dim = ",
           object$bs.dim, ". Check that bs.dim includes the di alpha columns.")
    }

    attr(Xtot, "inner_linpred_unscaled") <- cbind(z1 = za[[1]], z2 = za[[2]])

  } else {

    alpha <- si$alpha
    a0 <- si$a0

    Xi <- data[[term_1]]
    Xi <- t(t(Xi) - si$xm) %*% si$B[[1]]
    xa <- Xi %*% (alpha + a0)  # z

    if(get.xa){
      return(list(xa = xa, xa_da = Xi))
    }

    t_var <- data[[term_2]] - si$tm

    X0 <- object$xt$basis$evalX(z1 = xa, z2 = t_var, deriv = 0)$X0

    # [ 0_di, X_2D ]
    Xtot <- cbind(matrix(0, nrow(X0), di), X0)

    attr(Xtot, "inner_linpred_unscaled") <- xa
  }

  return( Xtot )
}
