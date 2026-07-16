#'
#' Build interaction_linear effect
#' 
#' @param Xi matrix to be projected via single index vector \code{alpha}.
#' @param basis function which takes \code{si = Xi\%*\%alpha} as input and returns model
#'                  matrix and its derivatives w.r.t. \code{si}.
#' @name eff_inter_nexp
#' @rdname eff_inter_nexp
#' @export eff_inter_nexp
#'
eff_inter_nexp <- function(y, Xi, t, basis, x0 = NULL, times = NULL){
  
  force(y); force(Xi); force(t); force(basis); force(x0); force(times);
  
  na <- ncol(Xi) + 1
  
  incall <- as.expression(quote(do.call("expsmooth", list("y" = y, "Xi" = Xi, "beta" = alpha, "times" = times, "deriv" = deriv), quote = TRUE)))
  efcall <- as.expression(quote(do.call("eff_inter_nexp", list("y" = y, "Xi" = Xi, "t" = t, "basis" = basis, "x0" = x0, "times" = times), quote = TRUE)))
  
  .eval <- .get_eff_eval_inter()
  
  environment(.eval) <- as.environment(environment())
  
  out <- structure(list("eval" = .eval), class = c("nexpsm", "nested"))
  
  return( out )
  
}









