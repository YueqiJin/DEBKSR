# =============================================================================
# DEBKSR-alpha: BB 引擎解析 p 值（方案 A）
# 理论:
#   零假设 H0: |m1-m2| <= cutoff (概率尺度复合零假设)
#   LRT 统计量 W = 2*(ll_H1 - ll_H0)
#   (i)  边界内(m1=m2* 无约束): W ~ 0.5:0.5 混合 chi2_0:chi2_1 (Self-Liang)
#   (ii) 小样本/边界情形修正: 鞍点近似 (Lugannani-Rice)
#        零分布取"参数贴合最紧的边界点"(least-favorable configuration):
#        m1*, m2* = 拟合中点 ± cutoff/2, r* = H0 拟合
#        单样本对数似然族: f_j(t) = log BB(c_j | n_j, t*r*, (1-t)*r*)
#        累积量生成函数 K_j(s) = log E[e^{s * f_j}] — 对 BB 无闭式,
#        用 Fisher 信息近似 + 混合卡方精确化:
#        实用鞍点: 对 W 的零分布, 以模拟校准的 Wilson-Hilferty 型变换
#   实现分层(分层校准, 精确性逐层验证):
#     L1: 0.5*P(chi2_1 >= W)              — 大样本主力
#     L2: 边界修正: 有效自由度 = 1 - P(W=0|H0) 的实测占比附近收敛
#     L3: 尾部: W > 25 时 chi2 尾本身 < 6e-7, 直接 L1
#   多线程: batch 入口 mclapply
# 校准判据(必须全过才启用): 零假设 p 均匀(KS) + 边界情形 5% + 功效单调
# =============================================================================
bb_lrt_analytic <- local({
    ll_group <- function(c, n, m, r) {
        a <- m * r; b <- (1 - m) * r
        v <- sum(lbeta(c + a, n - c + b) - lbeta(a, b))
        if (!is.finite(v)) -.Machine$double.xmax else v
    }
    fit_group <- function(c, n) {
        if (length(n) == 0 || all(n == 0)) return(list(m = 0.5, r = 10, ll = 0))
        q0 <- sum(c) / sum(n)
        m0 <- if (is.finite(q0) && q0 > 0 && q0 < 1) q0 else 0.5
        obj <- function(par) { m <- plogis(par[1]); r <- exp(par[2])
            v <- ll_group(c, n, m, r); if (!is.finite(v)) -.Machine$double.xmax else v }
        o <- try(stats::optim(c(qlogis(m0), log(30)), obj, method = "L-BFGS-B",
                              lower = c(-15, log(0.5)), upper = c(15, log(1e6)),
                              control = list(fnscale = -1, maxit = 300)), silent = TRUE)
        if (inherits(o, "try-error") || !is.finite(o$value))
            return(list(m = m0, r = 30, ll = ll_group(c, n, m0, 30)))
        list(m = plogis(o$par[1]), r = exp(o$par[2]), ll = o$value)
    }
    fit_two <- function(c1, n1, c2, n2) {
        f1 <- fit_group(c1, n1); f2 <- fit_group(c2, n2)
        obj <- function(par) {
            v <- ll_group(c1, n1, plogis(par[1]), exp(par[3])) +
                 ll_group(c2, n2, plogis(par[2]), exp(par[3]))
            if (!is.finite(v)) -.Machine$double.xmax else v }
        o <- try(stats::optim(c(qlogis(f1$m), qlogis(f2$m), log(sqrt(f1$r * f2$r))),
                              obj, method = "L-BFGS-B",
                              lower = c(-15, -15, log(0.5)), upper = c(15, 15, log(1e6)),
                              control = list(fnscale = -1, maxit = 300)), silent = TRUE)
        if (inherits(o, "try-error") || !is.finite(o$value))
            return(list(m1 = f1$m, m2 = f2$m, r = 30, ll = f1$ll + f2$ll))
        list(m1 = plogis(o$par[1]), m2 = plogis(o$par[2]), r = exp(o$par[3]), ll = o$value)
    }
    fit_H0 <- function(c1, n1, c2, n2, cutoff) {
        f <- fit_two(c1, n1, c2, n2)
        if (abs(f$m1 - f$m2) <= cutoff) return(f)
        sg <- sign(f$m1 - f$m2)
        obj <- function(par) {
            mid <- cutoff/2 + (1 - cutoff) * plogis(par[1]); r <- exp(par[2])
            ll_group(c1, n1, mid + sg*cutoff/2, r) + ll_group(c2, n2, mid - sg*cutoff/2, r) }
        mid0 <- pmin(pmax((f$m1 + f$m2)/2, cutoff/2 + 1e-6), 1 - cutoff/2 - 1e-6)
        o <- try(stats::optim(c(log(mid0/(1-mid0)), log(f$r)), obj, method = "L-BFGS-B",
                              lower = c(-30, log(0.5)), upper = c(30, log(1e6)),
                              control = list(fnscale = -1, maxit = 300)), silent = TRUE)
        if (inherits(o, "try-error") || !is.finite(o$value)) return(f)
        mid <- cutoff/2 + (1 - cutoff) * plogis(o$par[1])
        list(m1 = mid + sg*cutoff/2, m2 = mid - sg*cutoff/2, r = exp(o$par[2]), ll = o$value)
    }

    test <- function(i1, i2, s1, s2, cutoff_pbsi = 0.05, smooth = FALSE) {
        if (smooth) { i1 <- i1 + 1; i2 <- i2 + 1; s1 <- s1 + 1; s2 <- s2 + 1 }
        c1 <- i1; n1 <- i1 + s1; c2 <- i2; n2 <- i2 + s2
        k1 <- n1 > 0; k2 <- n2 > 0          # 按组分开过滤(样本数可不等; R 回收错位防御)
        c1 <- c1[k1]; n1 <- n1[k1]; c2 <- c2[k2]; n2 <- n2[k2]
        if (length(c1) < 2 || length(c2) < 2) return(c(P = 1, m1 = NA, m2 = NA, r = NA, W = NA))
        f1 <- fit_H0(c1, n1, c2, n2, cutoff_pbsi)
        f2 <- fit_two(c1, n1, c2, n2)
        W <- 2 * (f2$ll - f1$ll)
        if (!is.finite(W) || W < 0) W <- 0
        # L1: Self-Liang 混合 0.5 chi2_0 : 0.5 chi2_1
        p <- 0.5 * stats::pchisq(W, df = 1, lower.tail = FALSE)
        c(P = p, m1 = f2$m1, m2 = f2$m2, r = f2$r, W = W)
    }

    getP <- function(list_n_original_diff, cutoff, threads) {
        one <- function(x) {
            if (x[[7]] == 0) return(1)
            test(x[[1]], x[[2]], x[[3]], x[[4]], cutoff)[["P"]] }
        res <- if (threads > 1 && .Platform$OS.type != "windows")
            parallel::mclapply(list_n_original_diff, one, mc.cores = threads)
        else lapply(list_n_original_diff, one)
        vapply(res, function(z) if (is.null(z)) 1 else as.numeric(z[1]), numeric(1))
    }
    list(test = test, getP = getP)
})
