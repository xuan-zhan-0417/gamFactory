#' 交互型 Single Index 平滑效应的 construct 函数
#' 
#' @name smooth.construct.inter_linear.smooth.spec
#' @rdname smooth.construct.inter_linear.smooth.spec
#' @importFrom MASS Null
#' @importFrom Matrix rankMatrix
#' @export
#'
smooth.construct.inter_linear.smooth.spec <- function(object, data, knots){
  
  browser()
  # =========================================================================
  # 块 1：提取变量名与初始化 si 列表
  # =========================================================================
  # 识别交互项的两个变量：term_x 是多维矩阵，term_t 是时间变量(如 doy)
  if(length(object$term) != 2) {
    stop("The smooth effect must contain exactly two terms: the multivariate matrix and the time variable.")
  }
  term_x <- object$term[1]
  term_t <- object$term[2]
  
  # 提取 xt 中关于 single index (si) 的设定
  si <- object$xt$si
  if( is.null(si) ){ si <- object$xt$si <- list() }
  
  # =========================================================================
  # 块 2：处理内部多维数据 X 与 时间变量 t (预处理)
  # =========================================================================
  Xi <- data[[term_x]]
  t_vec <- data[[term_t]]
  
  # 1. 对时间变量 t 进行中心化 (保持其作为 Vector 的属性)
  t_mean <- mean(t_vec)
  data[[term_t]] <- t_vec - t_mean
  si$tm <- t_mean  # 保存 t 的均值，供预测时使用
  
  # 2. 对多维输入 Xi 进行中心化 (Xi 是 Matrix，用 scale 很合适)
  Xi <- scale(Xi, scale = FALSE)
  si$xm <- attr(Xi, "scaled:center") # 保存 X 各列的均值，供预测时使用
  
  di <- ncol(Xi)
  n <- nrow(Xi)
  
  # =========================================================================
  # 块 3：处理内部系数 alpha 的惩罚矩阵 (Inner Penalty)
  # =========================================================================
  Si <- si$S
  no_pen <- is.null(Si) && is.null(si$pord)
  
  if( no_pen ){ 
    # 情况 [a]: 内部系数 alpha 不受罚
    si$X <- Xi
    si$B <- diag(nrow = ncol(Xi))
    si$rank <- 0 
  } else {
    if( is.null(Si) ){ 
      # 情况 [b]: 默认使用 P-splines 惩罚（需调用内部函数 .psp）
      Si <- .psp(d = di, ord = si$pord)
      rankSi <- ncol(Xi) - si$pord
    } else { 
      # 情况 [c]: 用户自定义惩罚矩阵
      rankSi <- rankMatrix(Si)
    }
    # 对 Xi 进行重参数化，使得关于 single index 向量的惩罚矩阵是对角化的
    si <- append(si, gamFactory:::.diagPen(X = Xi, S = Si, r = rankSi))
  }
  
  # =========================================================================
  # 块 4：初始化及标准化 alpha 系数
  # =========================================================================
  # alpha 是内部系数向量，a0 是一个 offset (full_alpha = alpha + a0)
  if( is.null(si$a0) ){
    if( no_pen ){
      si$a0 <- rep(0, di)
    } else {
      si$a0 <- rep(1, di)
    }
  }
  
  if( is.null(si$alpha) ){ 
    if( is.null(si$a0) || all(si$a0 == 0) ){
      si$alpha <- rep(1, di) 
    } else {
      si$alpha <- rep(0, di) 
    }
  }
  
  # 应用块 3 中的 B 矩阵进行重参数化转换
  si$alpha <- solve(si$B, si$alpha)
  si$a0 <- solve(si$B, si$a0)
  
  # 施加方差为 1 的约束（为了使得模型可识别，因为外部平滑项能吸收尺度 variation）
  tmp <- sd(si$X %*% (si$alpha + si$a0))
  si$alpha <- si$alpha / tmp
  si$a0 <- si$a0 / tmp
  
  # =========================================================================
  # 块 5：计算一维投影，交接给基函数构建模块
  # =========================================================================
  # 计算当前的降维一维变量 z1 = x * t(alpha)
  ax <- drop( si$X %*% (si$alpha + si$a0) ) 
  
  # 核心设计：将 data 中的多维矩阵替换为计算出的一维向量 ax
  # 这样 .build_n_inter_bspline_basis 只需要处理两个一维变量（ax 和 t）的边际矩阵和张量积
  data[[term_x]] <- ax
  
  # 调用你提到的专门处理基函数、零空间约束、外推以及按行张量积构建 X_2D 的核心函数
  # 此时 data 里的变量都是一维的了
  out <- .build_n_inter_bspline_basis(object = object, data = data, knots = knots, si = si)
  
  # =========================================================================
  # 块 6：组装全局惩罚矩阵列表
  # =========================================================================
  # 在 mgcv 架构中，如果有内部 alpha 惩罚，我们需要将其拼接到外层参数的惩罚矩阵中
  if( !no_pen ){
    # 假设 out$bs.dim 是包含了内层和外层所有参数的总维度
    # 如果外层参数（样条系数）个数 = 总参数量 - 内部alpha的个数 (di)
    dsmo <- out$bs.dim - di 
    si <- out$xt$si
    
    # 构建包含 alpha 惩罚的扩展块对角矩阵。
    # 这里的顺序假设参数向量是 c(alpha_coeffs, spline_coeffs)
    alpha_penalty_padded <- rbind(
      cbind(si$S, matrix(0, di, dsmo)),
      cbind(matrix(0, dsmo, di), matrix(0, dsmo, dsmo))
    )
    
    # 将其加入到 out$S 的列表中（具体放置位置视 mgcv 优化器需要，通常加在最后）
    out$S[[length(out$S) + 1]] <- alpha_penalty_padded
    
    # 更新平滑项整体的零空间维度和每个惩罚矩阵对应的 rank
    out$null.space.dim <- out$null.space.dim + (out$bs.dim - si$rank)
    out$rank <- c(out$rank, si$rank)
  }
  
  # =========================================================================
  # 块 7：定义返回对象的 class
  # =========================================================================
  # 定义 class 使得 mgcv 的 Predict.matrix 和其他辅助函数能够识别这种新效应
  class(out) <- c("inter_linear", "nested")
  
  return( out )
}