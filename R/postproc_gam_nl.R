#' Post-process a GAM containing non-standard effects
#' @name postproc_gam_nl
#' @rdname postproc_gam_nl
#' @export postproc_gam_nl
#'
postproc_gam_nl <- function(o, info) {
  
  ne <- length(info$type)
  
  for (ii in seq_len(ne)) {
    types <- info$type[[ii]]
    
    # Coefficients of nested smooth must be saved in smooth object
    if ("nested" %in% types) {
      
      ism <- info$extra[[ii]]$ism
      sii <- o$smooth[[ism]]
      si  <- sii$xt$si
      alpha_len <- length(si$alpha)
      si$alpha <- coef(o)[info$iec[[ii]][seq_len(alpha_len)]]
      
      # refresh the centring constant of each margin under the final parameters
      # (exp / mgks / si_nexp margins are centred with mean(g); si / plain ones keep the constant of construction)
      si$margin <- lapply(si$margin, function(mk){ mk$xm <- .nest_bundle(mk$bundle_nam)$eval(mk, si$alpha[mk$idx])$xm; mk })
      sii$xt$si <- si

      jacobian <- get_jacobian.nested(sii, data = o$model, param = coef(o)[info$iec[[ii]]])
      sii$xt$jacobian <- jacobian$JJ
      sii$xt$xa       <- jacobian$xa  # inner indices (one column per margin), useful for plot
      o$smooth[[ism]] <- sii
    }
  }
  
  o$linear.predictors <- predict.gamnl(o, type = "link", newdata = o$model)
  o$fitted.values     <- predict.gamnl(o, type = "response", newdata = o$model)
  
  return(o)
}