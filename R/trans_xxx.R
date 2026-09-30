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
#'  \item{\code{trans_nexpsm}{ an exponential smoothing transformation.}}
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
trans_nexpsm <- function(S, alpha){
  
  out <- lapply(as.list(match.call())[-1], eval, envir = parent.frame())
  out$type <- "nexpsm"
  
  return(out)
  
} 

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
  out$type <- "si_nexpsm"
  return(out)
}


#
# Specifying a two-margin nested interaction, s(m_1, m_2)
# (see smooth.construct.inter.smooth.spec)
#
#' @rdname trans_xxx
#'
#' @param margins Character vector of length 2 giving the type of each margin of an interactive
#'                transformation: one of \code{"si"} (single index), \code{"exp"} (adaptive exponential smooth),
#'                \code{"mgks"} (kernel smooth), \code{"plain"} (ordinary covariate) or \code{"auto"}
#'                (inferred from the shape / column names of the term, the default). Supported combinations are
#'                listed in \link{smooth.construct.inter.smooth.spec}.
#' @param pord_1,pord_2 As \code{pord}, but for margin 1 / margin 2 of an interactive transformation. For an
#'                      \code{"exp"} margin the penalty acts on the smoothing-rate coefficients only.
#' @param S_1,S_2 As \code{S}, but for margin 1 / margin 2.
#' @param alpha_1,alpha_2 Initial values for the inner coefficients of margin 1 / margin 2 (for an \code{"exp"}
#'                        margin: the smoothing-rate coefficients, one per column of \code{"x"} within a replicate).
#' @param a0_1,a0_2 Fixed initialisation shifts of the coefficients of a \code{"si"} margin.
#' @param alpha_scale Initial value for the scaling parameter of the \code{"exp"} margin 2 (which multiplies the
#'                    centred output of the exponential smooth).
#'
#' @details \code{trans_inter} specifies any of the two-margin nested interactions
#'          \eqn{s(si(x), t)}, \eqn{s(si(x_1), si(x_2))}, \eqn{s(si(x), exp(x))}, \eqn{s(exp(x), t)},
#'          \eqn{s(exp(x_1), exp(x_2))} and \eqn{s(mgks(x), t)}. With the default \code{margins = c("auto", "auto")}
#'          the structure is inferred from the two terms passed to \link{s_nest}: a matrix with columns named
#'          \code{"y"} and \code{"x"} is an \code{"exp"} margin, one with \code{"y"} and \code{"d1"}, \code{"d2"}, ...
#'          is an \code{"mgks"} margin, any other matrix is a \code{"si"} margin, and a vector is a \code{"plain"}
#'          margin. For backward compatibility \code{pord}, \code{S} and \code{alpha} are accepted for margin 1
#'          of an \code{"exp"} or \code{"mgks"} margin.
#' @export trans_inter
#'
trans_inter <- function(
    margins = c("auto", "auto"),
    pord = NULL, S = NULL, alpha = NULL,
    pord_1 = NULL, pord_2 = NULL,
    S_1 = NULL, S_2 = NULL,
    alpha_1 = NULL, alpha_2 = NULL,
    a0_1 = NULL, a0_2 = NULL,
    alpha_scale = NULL
){
  out <- lapply(as.list(match.call())[-1], eval, envir = parent.frame())
  out$type <- "inter"
  out$margins <- margins
  return(out)
}
