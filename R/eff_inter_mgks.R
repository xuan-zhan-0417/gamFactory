#'
#' Build nested inter mgks smooth effect
#' 
#' @name eff_inter_mgks
#' @rdname eff_inter_mgks
#' @export eff_inter_mgks
#'
eff_inter_mgks <- function(y, t, dist, basis){
  
  force(y); force(t); force(dist); force(basis);
  
  na <- length(dist) + 1
  Xi <- NULL
  
  incall <- as.expression(quote(do.call("mgks", list("y" = y, "dist" = dist, "beta" = alpha, "deriv" = deriv), quote = TRUE)))
  efcall <- as.expression(quote(do.call("eff_inter_mgks", list("y" = y, "t" = t, "dist" = dist, "basis" = basis), quote = TRUE)))
  
  .eval <- .get_eff_eval_inter()
  
  environment(.eval) <- as.environment(environment())
  
  out <- structure(list("eval" = .eval), class = c("mgks", "nested"))
  
  return( out )
  
}

