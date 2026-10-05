#' @rdname nest_bundles
#' @export bundle_plain
#'
bundle_plain <- function(){
  list(
    bundle_nam = "bundle_plain",
    scale = FALSE,
    linear = TRUE,

    init = function(X, trans, label){
      if( !is.null(dim(X)) && ncol(as.matrix(X)) > 1 ){
        stop(label, " must be a single numeric vector when it is a plain (non-nested) margin.")
      }
      tm <- mean(X)
      list(z = X - tm, alpha = numeric(0), margin = list(x = X - tm, xm = tm))
    },

    eval = function(mk, par, deriv = 0, xm = NULL){
      list(z = mk$x, xm = mk$xm)
    },

    newdata = function(mk, X){
      mk$x <- X - mk$xm
      mk
    }
  )
}
