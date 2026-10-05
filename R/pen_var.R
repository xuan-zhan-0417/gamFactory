#'
#' Penalty on variance of the inner linear predictor
#' 
#' @description Penalty on variance of a nested effect. 
#' @param o the effect object.
#' @param v target variance of the inner linear predictor.
#' @param deriv number of derivatives to be calculated.
#' @param DaDr derivative for parameters of the inner linear predictor 
#'             w.r.t. the smoothing parameter.
#' @name pen_var
#' @rdname pen_var
#' @export pen_var
#'
pen_var <- function(o, v, deriv = 0){
  .pen_var_nest(o = o, v = v, deriv = deriv)
}

#' @name pen_var
#' @rdname pen_var
#' @export pen_var_outer
pen_var_outer <- function(o, v, DaDr){
  .pen_var_nest_outer(o = o, v = v, DaDr = DaDr)
}
