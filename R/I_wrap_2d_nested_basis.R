#' 2D nested basis wrapper (one or both margins nested)
#'
#' @param b1,b2 marginal bases, each exposing evalX(x, deriv).
#' @param nested logical, length 2: which margins are single indices and
#'               therefore need derivatives. Recycled if length 1.
#' @param simplify if TRUE and exactly one margin is nested, X1/X2/X3 are
#'                 returned as bare matrices instead of one-element lists,
#'                 preserving the single-sided API (store$X1 %*% beta).
#' @noRd
.wrap_2d_nested_basis <- function(b1, b2, nested = c(TRUE, TRUE), simplify = TRUE) {
  
  force(b1); force(b2);
  nested <- as.logical(nested)
  if (length(nested) == 1L) nested <- rep(nested, 2L)
  stopifnot(length(nested) == 2L, any(nested))
  
  idx <- which(nested)   # 需要求导的边际
  
  evalX <- function(z1, z2, deriv = 0) {
    
    ## 只对嵌套的一侧要高阶导：非嵌套的边际只算 X0，
    ## 那一侧的基也就不需要支持高阶
    out1 <- b1$evalX(x = z1, deriv = if (nested[1]) deriv else 0)
    out2 <- b2$evalX(x = z2, deriv = if (nested[2]) deriv else 0)
    
    n  <- nrow(out1$X0)
    Z1 <- matrix(0, n, ncol(out1$X0))
    Z2 <- matrix(0, n, ncol(out2$X0))
    
    ## ---------------------------------------------------------------------
    ## blk(a, b) = d^(a+b) X_2D / dz1^a dz2^b ，其中
    ##   X_2D = [ X1(z1) , X2(z2) , X1(z1) (x)_row X2(z2) ]
    ##
    ##   主效应 1 只依赖 z1  ->  b > 0 时恒零，否则 X1^(a)
    ##   主效应 2 只依赖 z2  ->  a > 0 时恒零，否则 X2^(b)
    ##   交互块是行向 Kronecker 积，对两个因子双线性，而两个因子各自只依赖
    ##   一个变量，所以没有乘积法则的交叉项：d^(a+b)(A (x) B) = A^(a) (x) B^(b)
    ##
    ## 混合偏导只取决于每个方向被求导的「次数」，与顺序无关（Schwarz），
    ## 所以按 (a, b) 记忆化就能让所有因对称而重合的高阶块只算一次。
    ## cache 是 evalX 的局部变量，每次调用重建 —— 绝不能跨不同的 z 复用。
    ## ---------------------------------------------------------------------
    cache <- list()
    blk <- function(a, b) {
      key <- paste0(a, "_", b)
      if (!is.null(cache[[key]])) return(cache[[key]])
      A <- out1[[paste0("X", a)]]
      B <- out2[[paste0("X", b)]]
      res <- cbind(if (b == 0) A else Z1,
                   if (a == 0) B else Z2,
                   mgcv::tensor.prod.model.matrix(list(A, B)))
      cache[[key]] <<- res
      res
    }
    
    ## 把一串求导方向（每个元素是 1 或 2）折成多重指标 (a, b)
    blk_dirs <- function(dirs) blk(sum(dirs == 1L), sum(dirs == 2L))
    
    out_2d <- list(X0 = blk(0, 0))
    
    ## 一阶：每个嵌套方向一项
    if (deriv >= 1) {
      g <- lapply(idx, blk_dirs)
      names(g) <- paste0("dz", idx)
      out_2d$X1 <- if (simplify && length(idx) == 1L) g[[1]] else g
    }
    
    ## 二阶：嵌套方向的所有有序配对；对称的名字指向同一个缓存块
    if (deriv >= 2) {
      h <- list()
      for (i in idx) for (j in idx) {
        h[[paste0("dz", i, "_z", j)]] <- blk_dirs(c(i, j))
      }
      out_2d$X2 <- if (simplify && length(idx) == 1L) h[[1]] else h
    }
    
    ## 三阶：同理，两侧嵌套时 8 个命名条目最多只对应 4 个不同的块
    if (deriv >= 3) {
      t3 <- list()
      for (i in idx) for (j in idx) for (k in idx) {
        t3[[paste0("dz", i, "_z", j, "_z", k)]] <- blk_dirs(c(i, j, k))
      }
      out_2d$X3 <- if (simplify && length(idx) == 1L) t3[[1]] else t3
    }
    
    out_2d
  }
  
  list(evalX = evalX)
}