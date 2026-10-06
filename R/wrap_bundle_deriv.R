#'
#' Check the derivatives of a margin bundle
#'
#' @description Wraps the \code{eval} function of a margin bundle (\link{trans_bundles}) so that its derivatives can
#'              be compared with finite differences by \link{check_deriv}, as done for a family with the derivatives of
#'              \link{llk_gaussian}. The function differentiated is \eqn{\sum_i w_i z_i}, where \eqn{z} is the index of
#'              the margin computed by the bundle at the term \code{X}.
#' @param bundle a margin bundle, e.g. \code{bundle_trans_exp()}.
#' @param X the term of the margin, as passed to \link{s_nest}.
#' @param trans optional arguments of the margin, e.g. \code{trans_linear(pord = 1)}.
#' @param w weights of the observations. Random by default: \eqn{\sum_i z_i} alone is always 0 because \eqn{z} is centred.
#' @return A list with the functions \code{d0}, \code{d1}, \code{d2}, \code{d3} of the inner parameters, to be passed
#'         to \link{check_deriv}, and \code{param}, the initial inner parameters set by the bundle.
#' @examples
#' # see ?trans_bundles
#' @name wrap_bundle_deriv
#' @rdname wrap_bundle_deriv
#' @export wrap_bundle_deriv
#'
wrap_bundle_deriv <- function(bundle, X, trans = NULL, w = NULL){

  init <- bundle$init(X, trans, "X")
  mk <- init$margin
  if( is.null(w) ){ w <- rnorm(length(init$z)) }

  # derivative of order k of sum(w * z): z for k = 0, then the columns of g1, g2, g3
  d <- function(k){
    function(param){
      o <- bundle$eval(mk, param, deriv = k)
      if( k == 0 ) sum(w * o$z) else colSums(w * o[[paste0("g", k)]])
    }
  }

  list(d0 = d(0), d1 = d(1), d2 = d(2), d3 = d(3), param = init$alpha)
}
