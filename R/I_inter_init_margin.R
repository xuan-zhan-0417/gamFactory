# Margin initialisers of an 'inter' smooth: one function per margin type, all returning the same kind of list.
# Used by smooth.construct.inter.smooth.spec (R/smooth_construct_inter_smooth_spec.R).
#   .inter_init_margin   dispatcher on the margin type
#   .inter_arg           user-supplied argument <name>_k of margin k
#   .inter_init_si / .inter_init_exp / .inter_init_mgks / .inter_init_plain
# The computation behind an exp margin lives in R/I_init_marginal_nexp.R, the one behind a si margin in
# R/I_init_marginal_si.R.

# ---------------------------------------------------------------------------
# Margin initialisers. Each returns
#   type, nested, z (inner index, centred), alpha (inner parameter block),
#   S / rank / pen_offset / pen_dim (penalty on the inner parameters, if any),
#   plus type-specific pieces used by .inter_assemble_si().
# ---------------------------------------------------------------------------
.inter_init_margin <- function(type, k, X, si, label){
  switch(type,
         "si"    = .inter_init_si(X, k, si, label),
         "exp"   = .inter_init_exp(X, k, si, label),
         "mgks"  = .inter_init_mgks(X, k, si, label),
         "plain" = .inter_init_plain(X, label))
}

# user-supplied argument of margin k: <nm>_k (margin 1 can also use the legacy un-suffixed name)
.inter_arg <- function(si, nm, k, legacy = FALSE){
  a <- si[[paste0(nm, "_", k)]]
  if( is.null(a) && legacy && k == 1L ){ a <- si[[nm]] }
  a
}

.inter_init_si <- function(X, k, si, label){

  X <- as.matrix(X)
  d <- ncol(X)
  a <- function(nm) .inter_arg(si, nm, k)

  if( !is.null(a("alpha")) && length(a("alpha")) != d ){
    stop("length(si$alpha_", k, ") must equal ncol(", label, ") = ", d, ".")
  }
  if( !is.null(a("a0")) && length(a("a0")) != d ){
    stop("length(si$a0_", k, ") must equal ncol(", label, ") = ", d, ".")
  }
  if( !is.null(a("S")) && (nrow(a("S")) != d || ncol(a("S")) != d) ){
    stop("si$S_", k, " must be a ", d, "x", d, " matrix.")
  }

  res <- .init_marginal_si(Xi = X, S = a("S"), pord = a("pord"), a0 = a("a0"), alpha = a("alpha"))

  list(type = "si", nested = TRUE, z = res$ax, res = res, alpha = res$alpha, na = d,
       S = res$S, rank = res$rank, pen_offset = 0, pen_dim = d)
}

.inter_init_exp <- function(X, k, si, label){

  a <- function(nm) .inter_arg(si, nm, k, legacy = TRUE)
  scale0 <- si[[paste0("alpha_scale_", k)]]
  if( is.null(scale0) && k == 2L ){ scale0 <- si$alpha_scale }

  res <- .init_marginal_nexp(X, S = a("S"), pord = a("pord"), alpha = a("alpha"),
                             alpha_scale = scale0, label = label)

  list(type = "exp", nested = TRUE, z = res$z, res = res, alpha = c(res$alpha_scale, res$alpha_w), na = res$na,
       S = res$S, rank = res$rank, pen_offset = 1, pen_dim = res$na)      # offset 1: the scaling parameter is never penalised
}

.inter_init_mgks <- function(X, k, si, label){

  nms <- colnames(X)

  # response at the reference locations
  y0 <- X[ , which(nms == "y"), drop = FALSE]
  if( !ncol(y0) ){ y0 <- si$y0 }

  # distance matrices d1, d2, ... (in order)
  Dist <- list()
  kk <- 1
  while( TRUE ){
    idx <- which(startsWith(nms, "d") & endsWith(nms, as.character(kk)) & sapply(nms, function(.x) nchar(.x) == 2))
    if( !length(idx) ){ break }
    Dist[[kk]] <- X[ , idx, drop = FALSE]
    kk <- kk + 1
  }

  alpha <- .inter_arg(si, "alpha", k, legacy = TRUE)
  if( is.null(alpha) ){
    # cold start: bandwidths set to -log(sd(Dist)/10)
    init_beta <- -log(sapply(Dist, sd) / 10)
    g <- mgks(y = y0, dist = Dist, beta = init_beta)$d0
    alpha <- c(log(1 / sd(g)), init_beta)
  } else if( length(alpha) == length(Dist) + 1 ){
    # warm start with the scaling parameter included
    g <- mgks(y = y0, dist = Dist, beta = alpha[-1])$d0
    alpha[1] <- log(1 / sd(g))
  } else {
    g <- mgks(y = y0, dist = Dist, beta = alpha)$d0
    alpha <- c(log(1 / sd(g)), alpha)
  }

  list(type = "mgks", nested = TRUE, z = exp(alpha[1]) * (g - mean(g)), alpha = alpha, y0 = y0, dist = Dist,
       S = NULL, rank = 0, pen_offset = 0, pen_dim = 0)
}

.inter_init_plain <- function(X, label){

  if( !is.null(dim(X)) && ncol(as.matrix(X)) > 1 ){
    stop(label, " must be a single numeric vector when it is a plain (non-nested) margin.")
  }
  tm <- mean(X)

  list(type = "plain", nested = FALSE, z = X - tm, tm = tm, alpha = numeric(0),
       S = NULL, rank = 0, pen_offset = 0, pen_dim = 0)
}
