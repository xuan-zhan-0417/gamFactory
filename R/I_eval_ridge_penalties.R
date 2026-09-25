########
# Evaluate additional ridge penalties and their derivatives
#
.eval_ridge_penalties <- function(eff, info, deriv){
  effType <- sapply(info$type, paste0, collapse = '')
  nstand <- which(effType != "stand")
  pen <- vector(mode = "list", length = length(nstand))
  kk <- 1
  for(ii in nstand){
    iec <- info$iec[[ii]]
    extra <- info$extra[[ii]]
    aii <- iec[1:eff[[ii]]$na]
    cl <- class(eff[[ii]])
    if("nested" %in% cl){
      if("si" %in% cl || "inter_dlinear" %in% cl){
        ipc <- 1:eff[[ii]]$na
      } else if("inter_le" %in% cl){
        # inner parameters are c(alpha_1, alpha_scale, alpha_2): exclude only the scaling parameter
        ipc <- (1:eff[[ii]]$na)[ -(eff[[ii]]$na1 + 1) ]
      } else if("inter_ee" %in% cl){
        # exclude the two scaling parameters (first element of each margin block)
        ipc <- (1:eff[[ii]]$na)[ -c(1, eff[[ii]]$d1 + 1) ]
      } else {
        ipc <- 2:eff[[ii]]$na
      }
        pen[[kk]] <- pen_ridge_var(o = eff[[ii]], extra = extra, ipc = ipc, deriv = deriv)
    }
    
    pen[[kk]]$iec <- aii 
    kk <- kk + 1
  }
  return( pen )
}

pen_ridge_var <- function(o, extra, ipc, deriv){
  
  a <- o$param[ipc]
  a_init <- extra$si$alpha[ipc]

  l0 <- .5 * sum((a-a_init)^2)
  
  l1 <- l2 <- NULL
  if(deriv){
    l1 <- rep(0, o$na)
    l1[ipc] <- a - a_init
    
    if(deriv > 1){
      l2 <- matrix(0, o$na, o$na)
      diag(l2)[ipc] <- 1
    }
  }
  return(list("d0" = l0, "d1" = l1, "d2" = l2))
}
