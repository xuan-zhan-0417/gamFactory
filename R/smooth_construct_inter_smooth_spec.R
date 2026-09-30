
#' Nested interactive (tensor-product) smooth effects
#'
#' @name smooth.construct.inter.smooth.spec
#' @rdname smooth.construct.inter.smooth.spec
#' @description One constructor for all two-margin nested interaction effects
#' \eqn{s(m_1, m_2)}, where each margin \eqn{m_k} is one of
#' \itemize{
#'   \item{\code{"si"}}{ a nested single index \eqn{X\alpha}: a numeric matrix.}
#'   \item{\code{"exp"}}{ a nested adaptive exponential smooth: a matrix with columns named
#'         \code{"y"} (data to smooth), \code{"x"} (smoothing-rate covariates) and, optionally, \code{"times"}.}
#'   \item{\code{"mgks"}}{ a nested multivariate kernel smooth: a matrix with columns named \code{"y"}
#'         and \code{"d1"}, \code{"d2"}, ... (distances). Margin 1 only.}
#'   \item{\code{"plain"}}{ an ordinary (non-nested) covariate: a numeric vector. Margin 2 only.}
#' }
#' Supported combinations (which decide the class of the resulting smooth, and hence the downstream code used):
#' \tabular{lll}{
#'   margins \tab class \tab meaning \cr
#'   si | plain \tab \code{inter_linear} \tab \eqn{s(si(x), t)} \cr
#'   si | si    \tab \code{inter_linear} \tab \eqn{s(si(x_1), si(x_2))} \cr
#'   si | exp   \tab \code{inter_le}     \tab \eqn{s(si(x), exp(x))} \cr
#'   exp | plain \tab \code{inter_nexp}  \tab \eqn{s(exp(x), t)} \cr
#'   exp | exp  \tab \code{inter_ee}     \tab \eqn{s(exp(x_1), exp(x_2))} \cr
#'   mgks | plain \tab \code{inter_mgks} \tab \eqn{s(mgks(x), t)}
#' }
#' The margin types are read from \code{object$xt$si$margins} (a length-2 character vector, set by the
#' \code{trans_inter} function); an element equal to \code{"auto"} (or a missing \code{margins}) is
#' inferred from the column names / shape of the corresponding term.
#' Margin-specific arguments are \code{pord_k}, \code{S_k}, \code{alpha_k}, \code{a0_k} (\code{k = 1, 2}); for
#' backward compatibility \code{pord}, \code{S}, \code{alpha} are accepted for margin 1 of \code{exp}/\code{mgks}
#' margins. \code{alpha_scale} (or \code{alpha_scale_k}) initialises the scaling parameter of an \code{exp} margin.
#' @importFrom MASS Null
#' @importFrom Matrix rankMatrix
#' @importFrom stats sd
#' @export
#'
smooth.construct.inter.smooth.spec <- function(object, data, knots){

  if( length(object$term) != 2 ){
    stop("An 'inter' smooth effect must contain exactly two terms.")
  }
  term <- object$term

  si <- object$xt$si
  if( is.null(si) ){ si <- object$xt$si <- list() }

  # name of structure
  mt <- .inter_margin_types(si$margins, data, term)
  struct <- paste(mt, collapse = "|")
  if( !(struct %in% names(.inter_classes)) ){
    stop("Unsupported combination of margins: ", struct, ". Supported: ",
         paste(names(.inter_classes), collapse = ", "), ".")
  }

  # initialise each margin (its inner index z and inner parameters)
  m <- vector("list", 2)
  for(k in 1:2){ m[[k]] <- .inter_init_margin(mt[k], k, data[[ term[k] ]], si, term[k]) }
  for(k in 1:2){ data[[ term[k] ]] <- m[[k]]$z }

  # assemble si
  si <- .inter_assemble_si(struct, si, m)

  # outer 2D basis (a margin is "nested" if it has inner parameters) 
  out <- .build_n_inter_bspline_basis(object = object, data = data, knots = knots, si = si,
                                       nested = c(m[[1]]$nested, m[[2]]$nested))

  # penalties on the inner parameters, padded into bs.dim space 
  off <- 0
  added_rank <- 0
  for(k in 1:2){
    mk <- m[[k]]
    if( !is.null(mk$S) && isTRUE(mk$rank > 0) ){
      out$S[[length(out$S) + 1]] <- .inter_pad_penalty(mk$S, off + mk$pen_offset, mk$pen_dim, out$bs.dim)
      out$rank <- c(out$rank, mk$rank); added_rank <- added_rank + mk$rank
    }
    off <- off + length(mk$alpha)
  }
  if( added_rank > 0 ){
    # each added penalty lives on alpha columns that were fully in the null space of S1/S2/S_inter
    out$null.space.dim <- out$null.space.dim - added_rank
  }

  # structure-specific finishing touches 
  if( struct == "si|si" ){
    out$xt$si$X <- cbind(out$xt$si$X[[1]], out$xt$si$X[[2]])
    attr(out$xt$si$X, "na1") <- si$na1
    attr(out$xt$si$X, "na2") <- si$na2
  } else if( struct == "si|plain" ){
    out$xt$si$X <- out$xt$si$X[[1]]
  }

  class(out) <- c(.inter_classes[[struct]], "nested")
  return( out )
}
