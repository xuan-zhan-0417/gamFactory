#' @rdname trans_bundles
#' @export bundle_trans_mgks
#'
bundle_trans_mgks <- function(){
  list(
    bundle_nam = "bundle_trans_mgks",
    scale = TRUE,
    linear = FALSE,

    init = function(X, trans, label){
      ex <- .extract_mgks_margin(X)
      y0 <- ex$y                                     # responses at the reference locations
      if( !ncol(y0) ){ y0 <- trans$y0 }
      Dist <- ex$dist                                # distance matrices d1, d2, ...

      alpha <- trans$alpha
      if( is.null(alpha) ){
        # cold start: bandwidths set to -log(sd(Dist)/10)
        init_beta <- -log(sapply(Dist, sd) / 10)
        g <- mgks(y = y0, dist = Dist, beta = init_beta)$d0
        alpha <- c(log(1 / sd(g)), init_beta)
      } else if( length(alpha) == length(Dist) + 1 ){
        # warm start with the scaling parameter included
        g <- mgks(y = y0, dist = Dist, beta = alpha[-1])$d0
        alpha[1] <- log(1 / sd(g))
      } else {
        g <- mgks(y = y0, dist = Dist, beta = alpha)$d0
        alpha <- c(log(1 / sd(g)), alpha)
      }

      list(z = exp(alpha[1]) * (g - mean(g)), alpha = alpha,
           margin = list(y = y0, dist = Dist, xm = mean(g)))
    },

    eval = function(mk, par, deriv = 0, xm = NULL){
      .nest_eval_scaled(mgks(y = mk$y, dist = mk$dist, beta = par[-1], deriv = deriv),
                         scale = par[1], deriv = deriv, xm = xm)
    },

    newdata = function(mk, X){
      ex <- .extract_mgks_margin(X)
      if( ncol(ex$y) ){ mk$y <- ex$y }               # otherwise keep the training responses
      mk$dist <- ex$dist
      mk
    }
  )
}

# Extract the responses ("y" columns) and the distance matrices ("d1", "d2", ..., in order) of one mgks margin
.extract_mgks_margin <- function(X){

  nms <- colnames(X)

  dist <- list()
  kk <- 1
  while( any(nms == paste0("d", kk)) ){
    dist[[kk]] <- X[ , nms == paste0("d", kk), drop = FALSE]
    kk <- kk + 1
  }

  list(y = X[ , nms == "y", drop = FALSE], dist = dist)
}
