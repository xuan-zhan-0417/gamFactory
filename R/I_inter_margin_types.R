# Structure of an 'inter' smooth: which combinations of margins are supported (and the class each one gets),
# and how the type of each margin is declared or inferred.
# Used by smooth.construct.inter.smooth.spec (R/smooth_construct_inter_smooth_spec.R).
#   .inter_classes       structure -> class table (the class is the key of all downstream dispatch)
#   .inter_margin_types  types of the two margins, from si$margins or inferred from the terms

# class of the smooth (= key for all downstream dispatch) for each supported structure
.inter_classes <- list("si|plain"   = "inter_linear",
                       "si|si"      = "inter_linear",
                       "si|exp"     = "inter_le",
                       "exp|plain"  = "inter_nexp",
                       "exp|exp"    = "inter_ee",
                       "mgks|plain" = "inter_mgks")


# ---------------------------------------------------------------------------
# Margin types: declared in si$margins, or inferred from the term itself
# ---------------------------------------------------------------------------
.inter_margin_types <- function(margins, data, term){

  if( is.null(margins) ){ margins <- c("auto", "auto") }
  if( length(margins) != 2 ){ stop("si$margins must have length 2.") }

  detect <- function(X, k){
    nms <- colnames(X)
    if( !is.null(nms) && all(c("y", "x") %in% nms) ){ return( "exp" ) }
    if( !is.null(nms) && "y" %in% nms && any(grepl("^d[0-9]$", nms)) ){ return( "mgks" ) }
    d <- if( is.null(dim(X)) ) 1L else ncol(X)
    if( k == 1L || d > 1L ) "si" else "plain"
  }

  mt <- vapply(1:2, function(k){
    if( identical(margins[k], "auto") ) detect(data[[ term[k] ]], k) else as.character(margins[k])
  }, character(1))

  bad <- setdiff(mt, c("si", "exp", "mgks", "plain"))
  if( length(bad) ){ stop("Unknown margin type(s): ", paste(bad, collapse = ", "), ".") }
  mt
}
