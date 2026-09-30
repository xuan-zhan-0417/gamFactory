# Penalty of an 'inter' smooth on the inner (alpha) parameters, padded into the full coefficient space.
# Used by smooth.construct.inter.smooth.spec (R/smooth_construct_inter_smooth_spec.R).
#   Sk       penalty matrix (dk x dk) on a block of dk inner parameters
#   off      number of coefficients preceding that block
#   dk       size of the block
#   bs.dim   total number of coefficients of the smooth (inner parameters + outer basis)
.inter_pad_penalty <- function(Sk, off, dk, bs.dim){
  P <- matrix(0, bs.dim, bs.dim)
  ii <- off + seq_len(dk)
  P[ii, ii] <- Sk
  P
}
