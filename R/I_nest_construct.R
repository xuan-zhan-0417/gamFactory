# Helpers of the constructor of nested effects (R/smooth_construct_nest_smooth_spec.R).
#   .nest_structures    supported structures: one margin, or "type of margin 1|type of margin 2"
#   .nest_margin_types  type of each margin, from its trans object or inferred from its term
#   .nest_assemble_si   the si list read by all the downstream code
#   .nest_pad_penalty   penalty on a block of inner parameters, padded into the full coefficient space

# The downstream code works margin by margin, so it does not depend on the structure; these are the tested ones
.nest_structures <- c("si", "si_pos", "exp", "mgks", "si_nexp",
                      "si|plain", "si|si", "si|exp", "exp|plain", "exp|exp", "mgks|plain")

.nest_margin_types <- function(trans, data, term){

  # no trans object: a matrix with columns "y" and "x" is exp, one with "y" and "d1", ... is mgks,
  # another matrix is si, a vector is plain (margin 2 only)
  detect <- function(X, k){
    nms <- colnames(X)
    if( !is.null(nms) && all(c("y", "x") %in% nms) ){ return( "exp" ) }
    if( !is.null(nms) && "y" %in% nms && any(grepl("^d[0-9]$", nms)) ){ return( "mgks" ) }
    d <- if( is.null(dim(X)) ) 1L else ncol(X)
    if( k == 1L || d > 1L ) "si" else "plain"
  }

  mt <- vapply(seq_along(term), function(k){
    if( is.null(trans[[k]]) ) detect(data[[ term[k] ]], k) else trans[[k]]$type
  }, character(1))

  bad <- setdiff(mt, c("si", "si_pos", "exp", "mgks", "plain", "si_nexp"))
  if( length(bad) ){ stop("Unknown margin type(s): ", paste(bad, collapse = ", "), ".") }
  mt
}

# m: outputs of the bundles' init, b: the bundles. si holds
#   alpha    inner parameters: those of margin 1, then those of margin 2
#   margin   per margin: bundle_nam, positions idx of its parameters in alpha, position iscale of its log-scale
#            parameter (bundles with scale = TRUE), centring constant xm and data (whatever the bundle's init stored)
.nest_assemble_si <- function(m, b){

  margin <- lapply(m, "[[", "margin")
  off <- 0
  for(k in seq_along(m)){
    margin[[k]]$bundle_nam <- b[[k]]$bundle_nam
    margin[[k]]$idx <- off + seq_along(m[[k]]$alpha)
    if( b[[k]]$scale ){ margin[[k]]$iscale <- off + 1 }
    off <- off + length(m[[k]]$alpha)
  }

  list(margin = margin, alpha = unlist(lapply(m, "[[", "alpha")))
}

.nest_pad_penalty <- function(Sk, off, bs.dim){
  P <- matrix(0, bs.dim, bs.dim)
  ii <- off + seq_len(nrow(Sk))
  P[ii, ii] <- Sk
  P
}
