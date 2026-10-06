#'
#' Bundles for the margins of nested effects
#'
#' @description Every nested effect (\link{smooth.construct.nest.smooth.spec}) is made of one or two margins, and
#'              each type of margin is described by a bundle, in the same way as a family is described by e.g.
#'              \link{bundle_fam_gaussian}. The generic code (constructor, \code{eff_nest}, predict, postproc) only
#'              calls the bundle; the smooth object stores its name (\code{si$margin[[k]]$bundle_nam}), never its
#'              functions. Bundles of inner transformations are named \code{bundle_trans_<type>}, bundles of families
#'              \code{bundle_fam_<family>}.
#' @return A list with elements
#' \itemize{
#'   \item{\code{bundle_nam}}{ name of the bundle function.}
#'   \item{\code{scale}}{ TRUE if the first inner parameter is a log-scale, \eqn{z = exp(scale) (g - xm)}.}
#'   \item{\code{linear}}{ TRUE if the index is linear in the inner parameters.}
#'   \item{\code{init(X, trans, label)}}{ initialisation from the term \code{X} and the user arguments \code{trans}
#'         (e.g. the output of \link{trans_linear}). Returns the initial index \code{z}, the initial inner parameters
#'         \code{alpha}, the penalties \code{pen} on them (each a list with \code{S}, \code{rank} and the position
#'         \code{offset} of its block) and \code{margin}, the data and centring constant \code{xm} stored in the
#'         smooth object.}
#'   \item{\code{eval(mk, par, deriv, xm)}}{ index \code{z} of margin \code{mk} at parameters \code{par} and its
#'         derivatives \code{g1}, \code{g2}, \code{g3} w.r.t. \code{par} (columns in \code{trind.generator} order).
#'         \code{xm = NULL} (fitting) centres with the current mean, otherwise (prediction) with the frozen \code{xm}.
#'         Also returns the centring constant \code{xm}.}
#'   \item{\code{newdata(mk, X)}}{ margin \code{mk} with its data replaced by the new term \code{X}.}
#' }
#' @details \code{bundle_trans_si}: single index \eqn{z = X (\alpha + a_0)}, X centred and rotated.
#'          \code{bundle_trans_exp}: adaptive exponential smooth. \code{bundle_trans_mgks}: kernel smooth.
#'          \code{bundle_trans_si_nexp}: exponential smooth of a single index. \code{bundle_trans_plain}: ordinary covariate,
#'          no inner parameters. The derivatives returned by \code{eval} can be checked with \link{check_deriv},
#'          see \link{wrap_bundle_deriv}.
#' @examples
#' library(gamFactory)
#' set.seed(1)
#' n <- 200
#'
#' # exponential smooth margin: compare exact and finite-difference derivatives up to order 3
#' E <- cbind(rnorm(n), 1, runif(n))
#' colnames(E) <- c("y", "x", "x")
#' obj <- wrap_bundle_deriv(bundle_trans_exp(), E)
#' der <- check_deriv(obj = obj, param = obj$param + 0.1, ord = 1:3)
#' sapply(der, function(d) max(abs(d[ , 1] - d[ , 2])))   # column 1: exact, column 2: finite differences
#'
#' # single index margin with a penalty
#' X <- matrix(rnorm(n * 3), n, 3)
#' obj <- wrap_bundle_deriv(bundle_trans_si(), X, trans = trans_linear(pord = 1))
#' der <- check_deriv(obj = obj, param = obj$param + 0.1, ord = 1:3)
#' sapply(der, function(d) max(abs(d[ , 1] - d[ , 2])))   # column 1: exact, column 2: finite differences
#' @name trans_bundles
#' @rdname trans_bundles
#' @export bundle_trans_si
#'
bundle_trans_si <- function(){
  list(
    bundle_nam = "bundle_trans_si",
    scale = FALSE,
    linear = TRUE,

    init = function(X, trans, label){
      X <- as.matrix(X)
      d <- ncol(X)
      if( !is.null(trans$alpha) && length(trans$alpha) != d ){
        stop("length(alpha) of the si margin ", label, " must equal ncol(", label, ") = ", d, ".")
      }
      if( !is.null(trans$a0) && length(trans$a0) != d ){
        stop("length(a0) of the si margin ", label, " must equal ncol(", label, ") = ", d, ".")
      }
      if( !is.null(trans$S) && (nrow(trans$S) != d || ncol(trans$S) != d) ){
        stop("S of the si margin ", label, " must be a ", d, "x", d, " matrix.")
      }
      res <- .init_marginal_si(Xi = X, S = trans$S, pord = trans$pord, a0 = trans$a0, alpha = trans$alpha)
      pen <- if( res$rank > 0 ) list(list(S = res$S, rank = res$rank, offset = 0))
      list(z = res$ax, alpha = res$alpha, pen = pen,
           margin = list(X = res$X, B = res$B, xm = res$xm, a0 = res$a0))
    },

    # z is linear in par, so all its second and third derivatives are zero
    eval = function(mk, par, deriv = 0, xm = NULL){
      n <- nrow(mk$X)
      p <- length(par)
      out <- list(z = drop(mk$X %*% (par + mk$a0)), xm = mk$xm, g1 = mk$X)
      if( deriv >= 2 ){ out$g2 <- matrix(0, n, p * (p + 1) / 2) }
      if( deriv >= 3 ){ out$g3 <- matrix(0, n, p * (p + 1) * (p + 2) / 6) }
      out
    },

    newdata = function(mk, X){
      mk$X <- t(t(as.matrix(X)) - mk$xm) %*% mk$B
      mk
    }
  )
}
