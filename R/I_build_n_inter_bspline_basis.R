.build_n_inter_bspline_basis <- function(object, data, knots, si) {
  
  # 1. 拆分 object
  k1 <- if(!is.null(object$bs.dim)) object$bs.dim[1] else 10
  k2 <- if(!is.null(object$bs.dim)) object$bs.dim[2] else 10
  obj1 <- object; obj1$term <- object$term[1]; obj1$bs.dim <- k1
  obj2 <- object; obj2$term <- object$term[2]; obj2$bs.dim <- k2
  
  # =====================================================================
  # 2. 核心适配：传入 si = NULL 提取纯净的样条设计矩阵
  # =====================================================================
  # 传入 NULL 使得 di = 0，从而避免原函数在 X 的左侧垫零
  out1 <- .build_nested_bspline_basis(obj1, data, knots[[1]], si = NULL)
  out2 <- .build_nested_bspline_basis(obj2, data, knots[[2]], si = NULL) 
  
  X1 <- out1$X; X2 <- out2$X
  n <- nrow(X1)
  p1 <- ncol(X1); p2 <- ncol(X2)
  
  # 3. 构造 2D 设计矩阵 (此时 X1 和 X2 是纯净的)
  X_2D <- mgcv::tensor.prod.model.matrix(list(X1, X2))
  X_spline <- cbind(X1, X2, X_2D)
  
  # 4. 严谨地构造块对角惩罚矩阵
  S1 <- out1$S[[1]]
  S2 <- out2$S[[1]]
  S1_inter <- kronecker(S1, diag(p2))
  S2_inter <- kronecker(diag(p1), S2)
  
  S1_spline <- Matrix::bdiag(S1, matrix(0, p2, p2), S1_inter)
  S2_spline <- Matrix::bdiag(matrix(0, p1, p1), S2, S2_inter)
  
  # =====================================================================
  # 5. 全局补零：在这里统一为 alpha 系数腾出位置
  # =====================================================================
  di <- length(si$alpha)
  dsmo_total <- ncol(X_spline)
  
  # 对 X 进行全局补零
  X_final <- cbind(matrix(0, n, di), X_spline)
  
  # 对 S1 和 S2 进行全局补零
  pad_mat   <- matrix(0, di, di)
  pad_cross <- matrix(0, di, dsmo_total)
  
  S1_final <- rbind(cbind(pad_mat, pad_cross), cbind(t(pad_cross), as.matrix(S1_spline)))
  S2_final <- rbind(cbind(pad_mat, pad_cross), cbind(t(pad_cross), as.matrix(S2_spline)))
  
  # 6. 组装返回列表
  out <- list()
  out$X <- as.matrix(X_final)
  out$S <- list(as.matrix(S1_final), as.matrix(S2_final))
  
  out$bs.dim <- ncol(X_final)
  out$rank <- c(Matrix::rankMatrix(S1_final), Matrix::rankMatrix(S2_final))
  out$null.space.dim <- out$bs.dim - sum(out$rank)
  out$df <- out$bs.dim      
  
  out$C <- matrix(0, 0, out$bs.dim)
  out$side.constrain <- FALSE
  out$no.rescale <- TRUE
  out$plot.me <- FALSE
  out$repara <- TRUE
  
  out$xt <- list()
  out$xt$si <- si # 将真实的 si 放回 xt 供 Predict.matrix 使用
  out$xt$basis <- .wrap_2d_nested_basis(out1$xt$basis, out2$xt$basis, di)
  
  return(out)
}