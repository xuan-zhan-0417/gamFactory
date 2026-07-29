.build_n_inter_bspline_basis <- function(object, data, knots, si) {
  
  # k <- object$bs.dim
  # if (is.null(k) || all(k < 0) || all(is.na(k))) {
  #   k1 <- 10
  #   k2 <- 10
  # } else {
  #   k1 <- if (length(k) >= 1 && !is.na(k[1]) && k[1] > 0) k[1] else 10
  #   k2 <- if (length(k) >= 2 && !is.na(k[2]) && k[2] > 0) k[2] else k1
  # }
  k1 <- 10
  k2 <- 10

  m <- object$p.order
  if (is.null(m) || any(is.na(m))) {
    # default：m[1]=2(order of penalty), m[2]=2(slope constrain)
    m <- c(2, 2) 
  }
  
  # fake maginal obj
  obj1 <- object
  obj1$term <- object$term[1]
  obj1$bs.dim <- k1
  obj1$p.order <- m
  
  obj2 <- object
  obj2$term <- object$term[2]
  obj2$bs.dim <- k2
  obj2$p.order <- m
  
  # si = null to avoid pad 0 in the left of X
  out1 <- .build_nested_bspline_basis(obj1, data, knots, si = NULL)
  out2 <- .build_nested_bspline_basis(obj2, data, knots, si = NULL) 
  
  X1 <- out1$X; X2 <- out2$X
  n <- nrow(X1)
  p1 <- ncol(X1); p2 <- ncol(X2)
  
  X_2D <- mgcv::tensor.prod.model.matrix(list(X1, X2))
  X_spline <- cbind(X1, X2, X_2D)
  
  # penalty matrix
  S1 <- out1$S[[1]]
  S2 <- out2$S[[1]]
  S_1_inter <- kronecker(S1, diag(p1))
  S_2_inter <- kronecker(S2, diag(p2))

  overlap <- (S_1_inter != 0) & (S_2_inter != 0)
  
  S3 <- matrix(0, nrow = nrow(S_1_inter), ncol = ncol(S_1_inter))
  S3[overlap] <- 1
  S_1_inter[overlap] <- 0
  S_2_inter[overlap] <- 0
  
  S1_spline <- Matrix::bdiag(S1, matrix(0, p2, p2), S_1_inter)
  S2_spline <- Matrix::bdiag(matrix(0, p1, p1), S2, S_2_inter)
  S_inter_spline <- Matrix::bdiag(matrix(0, p1, p1), matrix(0, p2, p2), S3)
  
  # pad 0 to X_2D and S_inter
  di <- length(si$alpha)
  dsmo_total <- ncol(X_spline)
  X_final <- cbind(matrix(0, n, di), X_spline)
  pad_mat   <- matrix(0, di, di)
  pad_cross <- matrix(0, di, dsmo_total)
  S1_final <- rbind(cbind(pad_mat, pad_cross), cbind(t(pad_cross), as.matrix(S1_spline)))
  S2_final <- rbind(cbind(pad_mat, pad_cross), cbind(t(pad_cross), as.matrix(S2_spline)))
  S_inter_final <-  rbind(cbind(pad_mat, pad_cross), cbind(t(pad_cross), as.matrix(S_inter_spline)))
    
  # out <- out1
  out <- object
  out$X <- as.matrix(X_final)
  out$S <- list(as.matrix(S1_final), as.matrix(S2_final), as.matrix(S_inter_final))
  out$bs.dim <- ncol(X_final)
  out$rank <- c(as.numeric(Matrix::rankMatrix(S1_final)), 
                as.numeric(Matrix::rankMatrix(S2_final)),
                as.numeric(Matrix::rankMatrix(S_inter_final)))
  out$null.space.dim <- out$bs.dim - sum(out$rank)
  out$df <- out$bs.dim      
  out$C <- matrix(0, 0, out$bs.dim)
  out$side.constrain <- FALSE
  out$no.rescale <- TRUE
  out$plot.me <- FALSE
  out$repara <- TRUE
  out$xt <- list()
  out$xt$si <- si # only one si for linear effect
  out$xt$basis <- .wrap_2d_nested_basis(b1 = out1$xt$basis, b2 = out2$xt$basis)
  out$xt$sumConv <- FALSE
  
  return(out)
}