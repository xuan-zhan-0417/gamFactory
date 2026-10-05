#'
#' Specifying inner transformations
#' 
#' @description Functions to specify nested covariate transformations, to be used
#'              within \link{s_nest}.

#' @rdname trans_xxx
#' 
#' @param pord The order of the differences of the transformation parameters
#'             that will be penalised. E.g., \code{pord = 1} corresponds to penalising the squared 
#'             first order differences \eqn{(alpha_k - alpha_{k-1})^2}. Ignored if argument \code{S}
#'             is supplied. 
#' @param S Positive semi-definite matrix used to specify a generalised ridge
#'          penalty on the transformation parameters. In \code{trans_linear}, if \code{S} is not supplied, the \code{pord}-order differences 
#'          between parameters are penalised, see the \code{pord} argument.
#' @param S_si Penalty matrix used in trans_linear_nexpsm. For the linear transformation parameters.
#' @param S_nexp Penalty matrix used in expsmooth. For the linear transformation parameters.
#' @param alpha Vector containing the initial values for the parameters of the transformation. If provided they can be used
#'              by the optimiser.
#' @param alpha_si Initial values used in trans_linear_nexpsm. For linear transformation parameters.
#' @param alpha_nexp Initial values used in trans_linear_nexpsm. For expsmooth parameters.
#' @param y0 Vector of \eqn{n} observations corresponding to the rows of the variables in \code{s_nest}.
#' @param n_si Length of \code{alpha_si} if \code{alpha_si} is not provided.
#' @param n_nexp Length of \code{alpha_nexp} if \code{alpha_nexp} is not provided.
#' 
#' @details The types of transformations currently provided are:
#' \itemize{
#'  \item{\code{trans_linear}{ a linear transformation \eqn{X^\top \alpha}, which can be used to specify, e.g, a single index vector (i.e. a projection) or
#'                       a linear effect.}}
#'  \item{\code{trans_mgks}{ a multivariate kernel smooth transformation based on the response variable observation 
#'                          vector \code{y0} and corresponding distance matrix.}}
#'  \item{\code{trans_nexpsm}{ an exponential smoothing transformation (\code{trans_exp} is the same function).}}
#'  \item{\code{trans_linear_nexpsm}{ a combination of a linear transformation and an exponential smoothing transformation.}}
#' }
#' @export trans_linear
#'
trans_linear <- function(pord, S, alpha, a0){
  
  out <- lapply(as.list(match.call())[-1], eval, envir = parent.frame())
  out$type <- "si"
  
  return(out)
  
}

#
# Specifying a MGKS transformation
# 
#' @rdname trans_xxx
#' @export trans_mgks
#'
trans_mgks <- function(y0, alpha){
  
  out <- lapply(as.list(match.call())[-1], eval, envir = parent.frame())
  out$type <- "mgks"
  
  return(out)
  
} 
# trans_mgks <- function(X0, y0, alpha){
#  
#   if( missing(X0) ){ stop("Argument \"X0\" is missing but a value is required") }
#    
#   out <- lapply(as.list(match.call())[-1], eval, envir = parent.frame())
#   out$type <- "mgks"
#  
#   return(out)
#   
# } 

#
# Specifying exponential smooth transformation 
# 
#' @rdname trans_xxx
#' @export trans_nexpsm
#'
trans_nexpsm <- function(pord, S, alpha, alpha_scale){

  out <- lapply(as.list(match.call())[-1], eval, envir = parent.frame())
  out$type <- "exp"

  return(out)

}

#' @rdname trans_xxx
#' @export trans_exp
#'
trans_exp <- trans_nexpsm

#
# Specifying linear transform + exponential smooth transformation 
# 
#' @rdname trans_xxx
#' @export trans_linear_nexpsm
#'
trans_linear_nexpsm <- function(
    pord_1       = NULL,
    pord_2       = NULL,
    S_si         = NULL,
    S_nexp       = NULL,
    alpha_nexp   = NULL,
    alpha_si     = NULL,
    alpha_scale  = NULL,
    center       = FALSE,
    alpha_center = NULL,
    n_si         = NULL,
    n_nexp       = NULL,
    Z0           = NULL,   # initialization for expsmooth
    positive_si = FALSE
){
  # if missing n_si/n_nexp，get it from alpha_si/alpha_nexp
  if (is.null(n_si)   && !is.null(alpha_si))   n_si   <- length(alpha_si)
  if (is.null(n_nexp) && !is.null(alpha_nexp)) n_nexp <- length(alpha_nexp)
  if (is.null(alpha_nexp) && is.null(alpha_si) && is.null(n_si) && is.null(n_nexp)) {
    stop("You must provide at least one of: alpha_si, alpha_nexp, n_si, n_nexp.")
  }
  
  out <- as.list(environment())
  out$type <- "si_nexp"
  return(out)
}


#
# Specifying a two-margin nested interaction, s(m_1, m_2)
# (see smooth.construct.nest.smooth.spec)
#
#' @rdname trans_xxx
#'
#' @param margin1,margin2 Transformation of margin 1 / margin 2 of an interactive transformation: the output of
#'                        \code{trans_linear} (single index), \code{trans_exp} (adaptive exponential smooth),
#'                        \code{trans_mgks} (kernel smooth) or \code{trans_plain} (ordinary covariate), whose
#'                        \code{type} selects the margin's bundle (\link{nest_bundles}). \code{NULL} (default):
#'                        inferred from the shape / column names of the term. Supported combinations are listed in
#'                        \link{smooth.construct.nest.smooth.spec}.
#' @param alpha_scale Initial value for the scaling parameter of an exponential smooth (which multiplies the
#'                    centred output of the exponential smooth). In \code{trans_nexpsm}, \code{pord} and \code{S}
#'                    penalise the smoothing-rate coefficients only.
#'
#' @details \code{trans_inter} specifies any of the two-margin nested interactions
#'          \eqn{s(si(x), t)}, \eqn{s(si(x_1), si(x_2))}, \eqn{s(si(x), exp(x))}, \eqn{s(exp(x), t)},
#'          \eqn{s(exp(x_1), exp(x_2))} and \eqn{s(mgks(x), t)}, e.g.
#'          \code{s_nest(X, E, trans = trans_inter(trans_linear(pord = 1), trans_exp()))}. A margin left to
#'          \code{NULL} is inferred from its term: a matrix with columns named \code{"y"} and \code{"x"} is an
#'          exponential smooth, one with \code{"y"} and \code{"d1"}, \code{"d2"}, ... a kernel smooth, any other
#'          matrix a single index and a vector a plain covariate. The arguments of \code{trans_linear},
#'          \code{trans_exp} and \code{trans_mgks} apply to that margin only; for \code{trans_exp} the penalty
#'          acts on the smoothing-rate coefficients and \code{alpha} has one element per column \code{"x"}
#'          within a replicate.
#' @export trans_inter
#'
trans_inter <- function(margin1 = NULL, margin2 = NULL){
  list(type = "inter", trans = list(margin1, margin2))
}

#' @rdname trans_xxx
#' @export trans_plain
#'
trans_plain <- function(){
  list(type = "plain")
}
