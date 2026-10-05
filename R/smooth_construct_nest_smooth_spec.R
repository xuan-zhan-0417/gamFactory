#' Nested smooth effects
#'
#' @name smooth.construct.nest.smooth.spec
#' @rdname smooth.construct.nest.smooth.spec
#' @description One constructor for all nested effects built by \link{s_nest}: \eqn{s(m_1)} with one term
#'              and \eqn{s(m_1, m_2)} (ANOVA tensor product) with two terms. Each margin \eqn{m_k} is an inner
#'              transformation of its term, described by a bundle (\link{nest_bundles}):
#' \itemize{
#'   \item{\code{"si"}}{ single index \eqn{X\alpha}, set by \link{trans_linear}. The term is a numeric matrix.}
#'   \item{\code{"exp"}}{ adaptive exponential smooth, set by \link{trans_nexpsm} (or \code{trans_exp}). The term is a
#'         matrix with columns named \code{"y"} (data to smooth), \code{"x"} (smoothing-rate covariates) and,
#'         optionally, \code{"times"}.}
#'   \item{\code{"mgks"}}{ multivariate kernel smooth, set by \link{trans_mgks}. The term is a matrix with columns
#'         named \code{"y"} and \code{"d1"}, \code{"d2"}, ... (distances).}
#'   \item{\code{"si_nexp"}}{ exponential smooth of a single index, set by \link{trans_linear_nexpsm}.}
#'   \item{\code{"plain"}}{ ordinary covariate (second margin only), set by \link{trans_plain}. The term is a vector.}
#' }
#' Supported structures: one margin of type si, exp, mgks or si_nexp; two margins si|plain, si|si, si|exp,
#' exp|plain, exp|exp, mgks|plain (\code{trans = trans_inter(margin1, margin2)}). A margin of
#' \code{trans_inter} without a trans object is inferred from its term.
#' @return A smooth of class \code{c("nest", "nested")}. \code{xt$si} holds \code{alpha} (the inner parameters,
#'         margin after margin) and \code{margin}, one description per margin (\code{bundle_nam}, positions
#'         \code{idx} of its parameters in \code{alpha}, position \code{iscale} of its scaling parameter, data
#'         and centring constant \code{xm}), which is all the downstream code reads.
#' @importFrom MASS Null
#' @importFrom Matrix rankMatrix
#' @importFrom stats sd cov
#' @importFrom mgcv smooth.construct Predict.matrix.Bspline.smooth sdiag bandchol psum.chisq
#' @export
#'
smooth.construct.nest.smooth.spec <- function(object, data, knots){

  term  <- object$term
  trans <- object$xt$si$trans

  # one bundle per margin, chosen by the type of the margin
  mt <- .nest_margin_types(trans, data, term)
  struct <- paste(mt, collapse = "|")
  if( !(struct %in% .nest_structures) ){
    stop("Unsupported structure: ", struct, ". Supported: ", paste(.nest_structures, collapse = ", "), ".")
  }
  b <- lapply(mt, function(type) .nest_bundle(paste0("bundle_", type)))

  # initialise each margin: its inner parameters, and its index z on which the outer basis is built
  m <- lapply(seq_along(term), function(k) b[[k]]$init(data[[ term[k] ]], trans[[k]], term[k]))
  for(k in seq_along(term)){ data[[ term[k] ]] <- m[[k]]$z }
  si <- .nest_assemble_si(m, b)

  # outer basis: B-splines of z1, or ANOVA tensor product of B-splines of z1 and z2
  if( length(term) == 1 ){
    out <- .build_nested_bspline_basis(object = object, data = data, knots = knots, si = si)
    out$xt$basis <- .wrap_1d_nested_basis(out$xt$basis)
  } else {
    out <- .build_n_inter_bspline_basis(object = object, data = data, knots = knots, si = si,
                                         nested = sapply(m, function(mk) length(mk$alpha) > 0))
  }

  # penalties on the inner parameters, padded into the full coefficient space
  off <- 0
  for(k in seq_along(m)){
    for(pk in m[[k]]$pen){
      out$S[[length(out$S) + 1]] <- .nest_pad_penalty(pk$S, off + pk$offset, out$bs.dim)
      out$rank <- c(out$rank, pk$rank)
      out$null.space.dim <- out$null.space.dim - pk$rank   # these parameters were unpenalised so far
    }
    off <- off + length(m[[k]]$alpha)
  }

  class(out) <- c("nest", "nested")
  return( out )
}
