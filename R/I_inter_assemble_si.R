# Last step before building the outer basis of an 'inter' smooth: put into si exactly the fields that the
# downstream code of each class (Predict.matrix.nested, get_jacobian.nested, pen_var, eff_*, postproc_gam_nl, ...) reads.
# Used by smooth.construct.inter.smooth.spec (R/smooth_construct_inter_smooth_spec.R).

# ---------------------------------------------------------------------------
# Put into si exactly what the downstream code of each class reads
# ---------------------------------------------------------------------------
.inter_assemble_si <- function(struct, si, m){

  m1 <- m[[1]]; m2 <- m[[2]]
  r1 <- m1$res; r2 <- m2$res

  if( struct == "si|plain" ){
    d1 <- m1$na
    si$t <- m2$z; si$tm <- m2$tm
    si$X <- list(r1$X, NULL); si$B <- list(r1$B, NULL); si$S <- list(r1$S, NULL)
    si$xm <- r1$xm; si$alpha <- r1$alpha; si$a0 <- r1$a0; si$rank <- c(r1$rank, 0)
    si$na1 <- d1; si$na2 <- 0; si$na <- d1

  } else if( struct == "si|si" ){
    d1 <- m1$na; d2 <- m2$na
    si$X <- list(r1$X, r2$X); si$B <- list(r1$B, r2$B); si$S <- list(r1$S, r2$S)
    si$xm <- c(r1$xm, r2$xm); si$alpha <- c(r1$alpha, r2$alpha); si$a0 <- c(r1$a0, r2$a0)
    si$rank <- c(r1$rank, r2$rank)
    si$na1 <- d1; si$na2 <- d2; si$na <- d1 + d2

  } else if( struct == "si|exp" ){                                    # inter_le
    na1 <- m1$na; na2 <- m2$na
    si$na1 <- na1; si$na2 <- na2; si$na <- na1 + 1 + na2
    si$X_1 <- r1$X; si$B_1 <- r1$B; si$S_1 <- r1$S; si$rank_1 <- r1$rank
    si$xm1 <- r1$xm; si$alpha_1 <- r1$alpha; si$a0_1 <- r1$a0
    si$y_raw <- r2$y_raw; si$W_2 <- r2$W; si$B_2 <- r2$B; si$S_2 <- r2$S
    si$rank_2 <- r2$rank; si$times <- r2$times
    si$xm2 <- r2$xm; si$alpha_scale <- r2$alpha_scale; si$alpha_2 <- r2$alpha_w
    si$alpha <- c(r1$alpha, r2$alpha_scale, r2$alpha_w)

  } else if( struct == "exp|plain" ){                                 # inter_nexp
    si$t <- m2$z; si$tm <- m2$tm
    si$X <- r1$W; si$B <- r1$B; si$rank <- r1$rank
    if( !r1$no_pen ){ si$S <- r1$S }
    si$alpha <- m1$alpha
    si$x_raw <- r1$y_raw; si$W_mat <- r1$W_raw; si$times <- r1$times

  } else if( struct == "exp|exp" ){                                   # inter_ee
    na1 <- m1$na; na2 <- m2$na
    # drop the user-supplied inputs so that only the internal (reparametrised) versions remain
    si$pord <- si$S <- si$alpha <- si$pord_1 <- si$pord_2 <- si$alpha_1 <- si$alpha_2 <- NULL
    si$alpha_scale <- NULL
    si$na1 <- na1; si$na2 <- na2; si$d1 <- 1 + na1; si$d2 <- 1 + na2; si$na <- si$d1 + si$d2
    for(k in 1:2){
      r <- m[[k]]$res
      si[[paste0("y_raw_", k)]] <- r$y_raw; si[[paste0("W_", k)]] <- r$W; si[[paste0("B_", k)]] <- r$B
      si[[paste0("rank_", k)]] <- r$rank; si[[paste0("times_", k)]] <- r$times
      si[[paste0("xm", k)]] <- r$xm
      si[[paste0("alpha_scale_", k)]] <- r$alpha_scale; si[[paste0("alpha_w_", k)]] <- r$alpha_w
    }
    si$S_1 <- r1$S; si$S_2 <- r2$S
    si$alpha <- c(m1$alpha, m2$alpha)

  } else if( struct == "mgks|plain" ){                                # inter_mgks
    si$t <- m2$z; si$tm <- m2$tm
    si$x <- m1$y0; si$dist <- m1$dist
    si$alpha <- m1$alpha
  }

  si
}
