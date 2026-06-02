#'❗️Maybe we can remove this function if we return "si" class in the end of "eff_inter"
#'
#' Derivative of inter_linear Hessian w.r.t smoothing parameters
#' 
#' @name DHessDrho.inter_linear
#' @rdname DHessDrho.inter_linear
#' @export 
#'
DHessDrho.inter_linear <- function(o, llk, DbDr){
  
  DHessDrho.si(o = o, llk = llk, DbDr = DbDr)
  
}
