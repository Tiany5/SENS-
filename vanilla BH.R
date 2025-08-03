vanilla_bh.func <- function(X,q,B){
  # 将输入转换为统一列表格式
  if (!is.list(X)) {
    X_list <- split(X, col(X))  # 按列分割矩阵为列表
  } else {
    X_list <- X
  }
  m <- length(X_list)
  
  # 计算单样本t统计量 (稳健化处理)
  compute_z_stat <- function(x, mu0 = 0) {
    x <- na.omit(x)  # 处理缺失值
    n <- length(x)
    if (n < 2) return(NA)  # 样本量不足时返回NA
    
    # 计算t统计量
    sample_mean <- mean(x)
    sample_sd <- sd(x)  # 样本标准差 (分母为n-1)
    t_stat <- (sample_mean - mu0) / (sample_sd / sqrt(n))
    
    # 计算自由度并应用转换公式
    df <- n - 1
    t_cdf <- pt(t_stat, df)        # t分布的CDF值
    z_value <- qnorm(t_cdf)        # 转换为标准正态分位数
    
    return(z_value)
  }
  
  # 并行初始化
  library(parallel)
  cl <- makeCluster(detectCores())
  on.exit(stopCluster(cl))  # 确保退出时关闭集群
  clusterExport(cl, c("compute_z_stat", "X_list"), envir = environment())
  
  # 并行生成置换统计量 (确保返回长度一致)
  z_perm <- parLapply(cl, 1:B, function(b) {
    tryCatch({
      sapply(X_list, function(x) {
        if (length(x) < 2) return(NA)
        signs <- sample(c(-1, 1), length(x), replace = TRUE)
        flipped_x <- x * signs
        compute_z_stat(flipped_x)
      })
    }, error = function(e) {
      warning("Permutation ", b, " failed: ", e$message)
      rep(NA, m)  # 确保始终返回长度m的向量
    })
  })
  
  # 转换为矩阵并转置
  z_perm <- matrix(unlist(z_perm), nrow = B, ncol = m, byrow = TRUE)
  
  # 计算原始统计量
  z_original <- sapply(X_list, compute_z_stat)
  
  # 计算双侧p值 (带平滑修正)
  p_values <- sapply(1:m, function(i) {
    extreme_count <- sum(abs(z_perm[, i]) >= abs(z_original[i]), na.rm = TRUE)
    (extreme_count + 1)/(B + 1)
  })
  y <- bh.func(p_values,q)
  return (y)
}
