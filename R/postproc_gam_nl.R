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
      
      # special case si_nexpsm, update alpha with positive constrain
      has_si_nexpsm <- "si_nexpsm" %in% types
      if (has_si_nexpsm) {
        n_nexp <- si$n_nexp
        n_si   <- si$n_si
        si$alpha_nexp <- as.vector(si$alpha[seq_len(n_nexp)])
        si$alpha_si   <- as.vector(si$alpha[(n_nexp + 1):(n_si + n_nexp)])

        if (isTRUE(si$positive_si)) {
          si$alpha_si_inner <- si$alpha_si             # optimized alpha_si
          si$alpha_si_true  <- exp(si$alpha_si)        # positive alpha_si
        }
      }

      # special case inter_le, split the combined alpha back into its
      # margin-1 (si) and margin-2 (nexp) pieces, which is what
      # .predict.matrix.inter_le / eff_inter_le actually read.
      has_inter_le <- "inter_le" %in% types
      if (has_inter_le) {
        na1 <- si$na1; na2 <- si$na2
        si$alpha_1     <- as.vector(si$alpha[1:na1])
        si$alpha_scale <- si$alpha[na1 + 1]
        si$alpha_2     <- as.vector(si$alpha[(na1 + 2):(na1 + 1 + na2)])
      }

      # special case inter_ee: split the combined alpha back into the two margin blocks
      # c(alpha_scale_1, alpha_w_1, alpha_scale_2, alpha_w_2), which is what
      # .predict.matrix.inter_ee actually reads.
      has_inter_ee <- "inter_ee" %in% types
      if (has_inter_ee) {
        d1 <- si$d1; d2 <- si$d2
        si$alpha_scale_1 <- si$alpha[1]
        si$alpha_w_1     <- as.vector(si$alpha[2:d1])
        si$alpha_scale_2 <- si$alpha[d1 + 1]
        si$alpha_w_2     <- as.vector(si$alpha[(d1 + 2):(d1 + d2)])
      }

      sii$xt$si <- si

      # Inner smooth must be centered using original data
      needs_base_xm <- is.null(si$xm)
      if (has_inter_le) {
        # margin-1 (si) xm1 is fixed raw column means -- never needs refreshing.
        # margin-2 (nexp) xm2 must be refreshed to mean(g) under the final
        # fitted alpha_2, exactly as for standalone nexpsm/inter_nexp.
        ip <- attr(Predict.matrix.nested(sii, data = o$model), "inner_linpred_unscaled")
        si$xm2 <- si$xm2 + mean(ip[ , "z2_unscaled"])
        sii$xt$si$xm2 <- si$xm2
      } else if (has_inter_ee) {
        # refresh xm1/xm2 to mean(g_k) under the final fitted parameters, as for nexpsm
        ip <- attr(Predict.matrix.nested(sii, data = o$model), "inner_linpred_unscaled")
        si$xm1 <- si$xm1 + mean(ip[ , "z1_unscaled"])
        si$xm2 <- si$xm2 + mean(ip[ , "z2_unscaled"])
        sii$xt$si$xm1 <- si$xm1
        sii$xt$si$xm2 <- si$xm2
      } else if (isTRUE(si$na2 > 0)) {
        na_tot <- length(si$alpha)
        if (needs_base_xm || length(si$xm) != na_tot) {
          stop("inter_linear (nested_2 = TRUE): si$xm must be a length-", na_tot,
               " vector of marginal column means (set in smooth.construct).")
        }
      } else if (needs_base_xm || has_si_nexpsm) {
        inner_mean <- mean(attr(Predict.matrix.nested(sii, data = o$model),
                                "inner_linpred_unscaled"))
        if (needs_base_xm) si$xm <- inner_mean
        if (has_si_nexpsm) si$xm <- c(si$xm, nexp = inner_mean)
        sii$xt$si$xm <- si$xm
      }
      
      jacobian <- get_jacobian.nested(sii, data = o$model, param = coef(o)[info$iec[[ii]]])
      sii$xt$jacobian <- jacobian$JJ
      sii$xt$xa       <- jacobian$xa  # smoothed s_t, useful for plot
      o$smooth[[ism]] <- sii
    }
  }
  
  o$linear.predictors <- predict.gamnl(o, type = "link", newdata = o$model)
  o$fitted.values     <- predict.gamnl(o, type = "response", newdata = o$model)
  
  return(o)
}