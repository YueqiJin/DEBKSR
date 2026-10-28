# =============================================================================
# DEBKSR-alpha: 贝塔-二项引擎（beta-binomial engine, engine = "bb"）
# =============================================================================
# 背景: study-3 发现 rMATS 继承的方差尺度错配(概率尺度方差 ÷ logit 尺度偏差);
#       study-7 pilot 模型选择(rmats-core 的模型族 = logit-normal+binomial)在
#       1,976 个拆半零假设位点上: 边界层参数不可辨识、模拟方差高估 1.4-2.1x、
#       覆盖率 12-50%; 而贝塔-二项全分层方差重现 0.95-1.00、覆盖 94-100%。
#
# 模型: 对每个位点, 每样本的环状接頭计数 c_s ~ Binomial(n_s, q_s),
#       q_s ~ Beta(m*r, (1-m)*r)   (均值 m × 精度 r 参数化)
#       组 = 两个 m (m1, m2); 检验 H0: |logit(m1)-logit(m2)| <= cutoff_logit
#       vs H1: 不受限。LRT 统计量 2*(ll_H1 - ll_H0), p = 1-pchisq(stat, df=1)。
#
# 似然: 给定 (m, r), 每样本边际似然有闭式(beta-binomial pmf):
#       log BB(c | n, a, b) = lchoose(n,c) + lbeta(c+a, n-c+b) - lbeta(a, b)
#       组对数似然 = sum_s log BB(c_s | n_s, a, b)   (各样本独立)
#       — 无 Laplace、无数值积分、无迭代 psi 优化。
#
# 精度参数 r: 逐位点 MLE(组间共享一个 r —— 与原引擎的 var 语义对应:
#       重复间过离散是位点属性, 不是组属性)。r 的估计在 H0 与 H1 下各自优化。
#
# 有效长度: 组间检验用原始计数分割 — 有效长度对同一位点是跨样本/跨组常数,
#       在二项分割似然中数学上约掉(改变的是计数的"单次试验成功概率"的尺度,
#       即报告的 PBSI, 不改变两组间的相对比较)。报告的 pbsi1/pbsi2 仍按
#       DEBKSR 原有效长度口径计算(adjustPsi), 保证输出表语义不变。
# =============================================================================

debks_bb_core <- local({
    env <- new.env(parent = baseenv())

    # ---- 每组对数似然(闭式): 向量 c(环状), n(总数), 共享 (m, r) ----
    ll_group <- function(c, n, m, r) {
        a <- m * r; b <- (1 - m) * r
        sum(lbeta(c + a, n - c + b) - lbeta(a, b))   # lchoose 项与参数无关, 检验中约掉
    }

    # ---- 单组 (m, r) MLE ----
    fit_group <- function(c, n) {
        if (length(n) == 0 || all(n == 0)) return(list(m = 0.5, r = 10, ll = 0))
        q0 <- sum(c) / sum(n)
        m0 <- if (is.finite(q0) && q0 > 0 && q0 < 1) q0 else 0.5
        obj <- function(par) {
            m <- plogis(par[1]); r <- exp(par[2])
            v <- ll_group(c, n, m, r)
            if (!is.finite(v)) -.Machine$double.xmax else v
        }
        o <- try(stats::optim(c(qlogis(m0), log(30)), obj, method = "L-BFGS-B",
                              lower = c(-15, log(0.5)), upper = c(15, log(1e6)),
                              control = list(fnscale = -1, maxit = 300)), silent = TRUE)
        if (inherits(o, "try-error") || !is.finite(o$value))
            return(list(m = m0, r = 30, ll = ll_group(c, n, m0, 30)))
        list(m = plogis(o$par[1]), r = exp(o$par[2]), ll = o$value)
    }

    # ---- 两组共享 r 的对数似然(profile: 给定 r, m 各自闭式一维优化不易, 直接二维) ----
    ll_two <- function(c1, n1, c2, n2, m1, m2, r) {
        ll_group(c1, n1, m1, r) + ll_group(c2, n2, m2, r)
    }
    fit_two_shared_r <- function(c1, n1, c2, n2) {
        f1 <- fit_group(c1, n1); f2 <- fit_group(c2, n2)
        obj <- function(par) {
            m1 <- plogis(par[1]); m2 <- plogis(par[2]); r <- exp(par[3])
            v <- ll_two(c1, n1, c2, n2, m1, m2, r)
            if (!is.finite(v)) -.Machine$double.xmax else v
        }
        o <- try(stats::optim(c(qlogis(f1$m), qlogis(f2$m), log(sqrt(f1$r * f2$r))),
                              obj, method = "L-BFGS-B",
                              lower = c(-15, -15, log(0.5)), upper = c(15, 15, log(1e6)),
                              control = list(fnscale = -1, maxit = 300)), silent = TRUE)
        if (inherits(o, "try-error") || !is.finite(o$value))
            return(list(m1 = f1$m, m2 = f2$m, r = 30, ll = f1$ll + f2$ll))
        list(m1 = plogis(o$par[1]), m2 = plogis(o$par[2]), r = exp(o$par[3]), ll = o$value)
    }

    # ---- H0: |m1 - m2| <= cutoff(概率尺度 — 与原引擎 cutoff 语义一致:
    #        原引擎 theta2 + cutoff 的约束就在概率尺度) ----
    fit_H0 <- function(c1, n1, c2, n2, cutoff) {
        f <- fit_two_shared_r(c1, n1, c2, n2)
        if (abs(f$m1 - f$m2) <= cutoff) return(f)
        sg <- sign(f$m1 - f$m2)
        obj <- function(par) {   # par: 均值中点(概率尺度, logit 参数化保证界内) + log r
            mid <- plogis(par[1]); r <- exp(par[2])
            m1 <- mid + sg * cutoff / 2; m2 <- mid - sg * cutoff / 2
            if (m1 <= 0 || m1 >= 1 || m2 <= 0 || m2 >= 1) return(-.Machine$double.xmax)
            v <- ll_two(c1, n1, c2, n2, m1, m2, r)
            if (!is.finite(v)) -.Machine$double.xmax else v
        }
        o <- try(stats::optim(c(qlogis((f$m1 + f$m2) / 2), log(f$r)),
                              obj, method = "L-BFGS-B",
                              lower = c(-15, log(0.5)), upper = c(15, log(1e6)),
                              control = list(fnscale = -1, maxit = 300)), silent = TRUE)
        if (inherits(o, "try-error") || !is.finite(o$value)) return(f)
        mid <- plogis(o$par[1]); r <- exp(o$par[2])
        list(m1 = mid + sg * cutoff / 2, m2 = mid - sg * cutoff / 2, r = r, ll = o$value)
    }

    # ---- 主检验入口 ----
    # 输入遵循 stat_r 约定: i/s 已经过 vecAddOne(+1 平滑)与无效样本过滤
    assign("bb_test_one", function(i1, i2, s1, s2, cutoff_pbsi = 0.05) {
        c1 <- i1; n1 <- i1 + s1
        c2 <- i2; n2 <- i2 + s2
        keep <- (n1 > 0) & (n2 > 0)
        c1 <- c1[keep]; n1 <- n1[keep]; c2 <- c2[keep]; n2 <- n2[keep]
        if (length(c1) < 2 || length(c2) < 2) return(c(P = 1, m1 = NA, m2 = NA, r = NA))
        f1 <- fit_H0(c1, n1, c2, n2, cutoff_pbsi)
        f2 <- fit_two_shared_r(c1, n1, c2, n2)
        stat <- 2 * (f2$ll - f1$ll)
        if (!is.finite(stat) || stat < 0) stat <- 0
        p <- stats::pchisq(stat, df = 1, lower.tail = FALSE)
        c(P = p, m1 = f2$m1, m2 = f2$m2, r = f2$r)
    }, envir = env)

    # ---- 批量入口(与 getrMATS_P 同签名, 可直接注入 stat_r) ----
    assign("bb_getP", function(list_n_original_diff, cutoff, threads) {
        one <- function(x) {
            if (x[[7]] == 0) return(c(P = 1, m1 = NA, m2 = NA, r = NA))
            bb_test_one(x[[1]], x[[2]], x[[3]], x[[4]], cutoff)
        }
        res <- if (threads > 1 && .Platform$OS.type != "windows") {
            parallel::mclapply(list_n_original_diff, one, mc.cores = threads)
        } else lapply(list_n_original_diff, one)
        vapply(res, function(z) z[["P"]], numeric(1))
    }, envir = env)

    env
})
