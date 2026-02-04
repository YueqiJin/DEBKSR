.debks_rmats_env <- new.env(parent = baseenv())

local({
    env <- environment()
    # Try to load C++ helpers for hot paths (optional, skips if Rcpp missing)
    if (requireNamespace("Rcpp", quietly = TRUE)) {
        cpp_path <- "DEBKSR/src/rmats_helpers.cpp"
        if (file.exists(cpp_path)) {
            try(Rcpp::sourceCpp(cpp_path, env), silent = TRUE)
        }
    }
    assign("rho", 0.9, envir = env)
if (exists("logit_clamped_cpp", envir = env, inherits = FALSE)) {
    assign("logit", function(x) logit_clamped_cpp(x), envir = env)
} else {
    assign("logit", function(x) {
            x <- pmin(pmax(x, 0.01), 0.99)
            log(x / (1 - x))
        }, envir = env)
}

if (exists("pop_var_cpp", envir = env, inherits = FALSE)) {
    assign("pop_var", function(x) pop_var_cpp(x), envir = env)
} else {
    assign("pop_var", function(x) {
            m <- mean(x)
            mean((x - m)^2)
        }, envir = env)
}

assign("optim_ctrl", list(factr = 1e7, pgtol = 1e-5, maxit = 15000), envir = env)

assign("myfunc_multivar", function(x, psi1, psi2, var1, var2) {
        sum1 <- sum((logit(psi1) - logit(x[1]))^2) / var1 / 2
        sum2 <- sum((logit(psi2) - logit(x[2]))^2) / var2 / 2
        sum1 + sum2 + 0.1 * 0.5 * (rho^2) / (1 - rho^2) *
            (stats::qnorm(x[1])^2 + stats::qnorm(x[2])^2 - 2 * rho * stats::qnorm(x[1]) * stats::qnorm(x[2]))
    }, envir = env)

assign("myfunc_multivar_der", function(x, psi1, psi2, var1, var2) {
        sum1 <- sum(-2 * (logit(psi1) - logit(x[1])) / x[1] / (1 - x[1])) / var1 / 2
        sum2 <- sum(-2 * (logit(psi2) - logit(x[2])) / x[2] / (1 - x[2])) / var2 / 2
        res1 <- sum1 + 0.1 * 0.5 * (rho^2) / (1 - rho^2) *
            (2 * stats::qnorm(x[1]) - 2 * rho * stats::qnorm(x[2])) / stats::dnorm(stats::qnorm(x[1]))
        res2 <- sum2 + 0.1 * 0.5 * (rho^2) / (1 - rho^2) *
            (2 * stats::qnorm(x[2]) - 2 * rho * stats::qnorm(x[1])) / stats::dnorm(stats::qnorm(x[2]))
        c(res1, res2)
    }, envir = env)

assign("myfunc_1", function(x, psi1, psi2, var1, var2, cutoff) {
        sum1 <- sum((logit(psi1) - logit(x + cutoff))^2) / var1 / 2
        sum2 <- sum((logit(psi2) - logit(x))^2) / var2 / 2
        sum1 + sum2 + 0.1 * 0.5 * (rho^2) / (1 - rho^2) *
            (stats::qnorm(x + cutoff)^2 + stats::qnorm(x)^2 - 2 * rho * stats::qnorm(x + cutoff) * stats::qnorm(x))
    }, envir = env)

assign("myfunc_der_1", function(x, psi1, psi2, var1, var2, cutoff) {
        sum1 <- sum(-2 * (logit(psi1) - logit(x + cutoff)) / (x + cutoff) / (1 - x - cutoff)) / var1 / 2
        res1 <- sum1 + 0.1 * 0.5 * (rho^2) / (1 - rho^2) *
            (2 * stats::qnorm(x + cutoff) - 2 * rho * stats::qnorm(x)) / stats::dnorm(stats::qnorm(x + cutoff))
        sum2 <- sum(-2 * (logit(psi2) - logit(x)) / x / (1 - x)) / var2 / 2
        res2 <- sum2 + 0.1 * 0.5 * (rho^2) / (1 - rho^2) *
            (2 * stats::qnorm(x) - 2 * rho * stats::qnorm(x + cutoff)) / stats::dnorm(stats::qnorm(x))
        res1 + res2
    }, envir = env)

assign("myfunc_2", function(x, psi1, psi2, var1, var2, cutoff) {
        sum1 <- sum((logit(psi1) - logit(x))^2) / var1 / 2
        sum2 <- sum((logit(psi2) - logit(x + cutoff))^2) / var2 / 2
        sum1 + sum2 + 0.1 * 0.5 * (rho^2) / (1 - rho^2) *
            (stats::qnorm(x + cutoff)^2 + stats::qnorm(x)^2 - 2 * rho * stats::qnorm(x + cutoff) * stats::qnorm(x))
    }, envir = env)

assign("myfunc_der_2", function(x, psi1, psi2, var1, var2, cutoff) {
        sum1 <- sum(-2 * (logit(psi1) - logit(x)) / x / (1 - x)) / var1 / 2
        res1 <- sum1 + 0.1 * 0.5 * (rho^2) / (1 - rho^2) *
            (2 * stats::qnorm(x + cutoff) - 2 * rho * stats::qnorm(x)) / stats::dnorm(stats::qnorm(x + cutoff))
        sum2 <- sum(-2 * (logit(psi2) - logit(x + cutoff)) / (x + cutoff) / (1 - x - cutoff)) / var2 / 2
        res2 <- sum2 + 0.1 * 0.5 * (rho^2) / (1 - rho^2) *
            (2 * stats::qnorm(x) - 2 * rho * stats::qnorm(x + cutoff)) / stats::dnorm(stats::qnorm(x))
        res1 + res2
    }, envir = env)

assign("myfunc_marginal_2der", function(x, I, S, beta, var, inclusion_length, skipping_length) {
        temp1 <- 1 / x^2 / (1 - x)^2 * (((2 * x - 1) * (logit(beta) - logit(x)) - 1) / var - 1)
        temp2 <- I * skipping_length * ((2 * inclusion_length + skipping_length) * x + skipping_length * (1 - x)) /
            x^2 / (inclusion_length * x + skipping_length * (1 - x))^2
        temp3 <- S * inclusion_length * ((inclusion_length + 2 * skipping_length) * (1 - x) + inclusion_length * x) /
            (1 - x)^2 / (inclusion_length * x + skipping_length * (1 - x))^2
        temp1 - temp2 - temp3
    }, envir = env)

assign("myfunc_marginal", function(x, I, S, psi, var, inclusion_length, skipping_length) {
        beta <- x
        sumv <- 0
        for (i in seq_along(psi)) {
            new_psi <- inclusion_length * psi[i] / (inclusion_length * psi[i] + skipping_length * (1 - psi[i]))
            f1 <- I[i] * log(new_psi) + S[i] * log(1 - new_psi) -
                (logit(psi[i]) - logit(beta))^2 / 2 / var - log(psi[i]) - log(1 - psi[i])
            f1_2der <- abs(myfunc_marginal_2der(psi[i], I[i], S[i], beta, var, inclusion_length, skipping_length))
            sumv <- sumv + (0.5 * log(abs(f1_2der) + 1e-05) - f1)
        }
        sumv
    }, envir = env)

assign("myfunc_marginal_der", function(x, I, S, psi, var, inclusion_length, skipping_length) {
        beta <- x
        sumv <- 0
        for (i in seq_along(psi)) {
            new_psi <- inclusion_length * psi[i] / (inclusion_length * psi[i] + skipping_length * (1 - psi[i]))
            f1_3der <- (2 * psi[i] - 1) / psi[i]^2 / (1 - psi[i])^2 / beta / (1 - beta) / var
            f1_2der <- myfunc_marginal_2der(psi[i], I[i], S[i], beta, var, inclusion_length, skipping_length)
            f1_der <- (logit(psi[i]) - logit(beta)) / beta / (1 - beta) / var
            f1 <- I[i] * log(new_psi) + S[i] * log(1 - new_psi) -
                (logit(psi[i]) - logit(beta))^2 / 2 / var + log(psi[i]) + log(1 - psi[i])
            sumv <- sumv + (0.5 * f1_3der / f1_2der - f1_der)
        }
        sumv
    }, envir = env)

assign("myfunc_marginal_1", function(x, I1, S1, psi1, var1, I2, S2, psi2, var2, inclusion_length, skipping_length, cutoff) {
        beta2 <- x
        beta1 <- x + cutoff
        myfunc_marginal(beta1, I1, S1, psi1, var1, inclusion_length, skipping_length) +
            myfunc_marginal(beta2, I2, S2, psi2, var2, inclusion_length, skipping_length)
    }, envir = env)

assign("myfunc_marginal_1_der", function(x, I1, S1, psi1, var1, I2, S2, psi2, var2, inclusion_length, skipping_length, cutoff) {
        beta2 <- x
        beta1 <- x + cutoff
        myfunc_marginal_der(beta1, I1, S1, psi1, var1, inclusion_length, skipping_length) +
            myfunc_marginal_der(beta2, I2, S2, psi2, var2, inclusion_length, skipping_length)
    }, envir = env)

assign("myfunc_marginal_2", function(x, I1, S1, psi1, var1, I2, S2, psi2, var2, inclusion_length, skipping_length, cutoff) {
        beta2 <- x + cutoff
        beta1 <- x
        myfunc_marginal(beta1, I1, S1, psi1, var1, inclusion_length, skipping_length) +
            myfunc_marginal(beta2, I2, S2, psi2, var2, inclusion_length, skipping_length)
    }, envir = env)

assign("myfunc_marginal_2_der", function(x, I1, S1, psi1, var1, I2, S2, psi2, var2, inclusion_length, skipping_length, cutoff) {
        beta2 <- x + cutoff
        beta1 <- x
        myfunc_marginal_der(beta1, I1, S1, psi1, var1, inclusion_length, skipping_length) +
            myfunc_marginal_der(beta2, I2, S2, psi2, var2, inclusion_length, skipping_length)
    }, envir = env)

if (exists("myfunc_individual_cpp", envir = env, inherits = FALSE)) {
    assign("myfunc_individual", function(x, I, S, beta, var, inclusion_length, skipping_length) {
        myfunc_individual_cpp(x, I, S, beta, var, inclusion_length, skipping_length)
    }, envir = env)
    assign("myfunc_individual_der", function(x, I, S, beta, var, inclusion_length, skipping_length) {
        myfunc_individual_der_cpp(x, I, S, beta, var, inclusion_length, skipping_length)
    }, envir = env)
} else {
    assign("myfunc_individual", function(x, I, S, beta, var, inclusion_length, skipping_length) {
        new_psi <- inclusion_length * x / (inclusion_length * x + skipping_length * (1 - x))
        -1 * (I * log(new_psi) + S * log(1 - new_psi) -
              (logit(x) - logit(beta)) * (logit(x) - logit(beta)) / 2 / var - log(x) - log(1 - x))
    }, envir = env)

    assign("myfunc_individual_der", function(x, I, S, beta, var, inclusion_length, skipping_length) {
        new_psi <- inclusion_length * x / (inclusion_length * x + skipping_length * (1 - x))
        new_psi_der <- inclusion_length * skipping_length / (inclusion_length * x + skipping_length * (1 - x))^2
        -1 * (I / new_psi * new_psi_der - S / (1 - new_psi) * new_psi_der -
              (logit(x) - logit(beta)) / var / x / (1 - x) - 1 / x + 1 / (1 - x))
    }, envir = env)
}

assign("myfunc_likelihood", function(x, I, S, beta, var) {
        N <- I + S
        if (N == 0) return(0)
        -0.5 * ((I - N * x) * (I - N * x) / (N * x) + (S - N * (1 - x)) * (S - N * (1 - x)) / (N * (1 - x))) -
            log(sqrt(var)) - (x - beta) * (x - beta) / 2 / var
    }, envir = env)

    if (exists("psi_optimize_golden_cpp", envir = env, inherits = FALSE)) {
        assign("optimize_psi_batch", function(psi_init, I, S, beta, var, inclusion_length, skipping_length) {
            res <- psi_optimize_golden_cpp(I, S, beta, var, inclusion_length, skipping_length,
                                           lower = 0.01, upper = 0.99)
            list(psi = as.numeric(res$psi), value = as.numeric(res$value))
        }, envir = env)
    } else {
        assign("optimize_psi_batch", function(psi_init, I, S, beta, var, inclusion_length, skipping_length) {
            psi <- numeric(length(psi_init))
            value <- numeric(length(psi_init))
            for (i in seq_along(psi_init)) {
                xopt <- stats::optim(par = psi_init[i], fn = myfunc_individual, gr = myfunc_individual_der,
                                     I = I[i], S = S[i], beta = beta, var = var,
                                     inclusion_length = inclusion_length, skipping_length = skipping_length,
                                     method = "L-BFGS-B", lower = 0.01, upper = 0.99,
                                     control = optim_ctrl)
                psi[i] <- as.numeric(xopt$par)
                value[i] <- as.numeric(xopt$value)
            }
            list(psi = psi, value = value)
        }, envir = env)
    }

    MLE_iteration_constrain <- function(i1, i2, s1, s2, effective_inclusion_length, effective_skipping_length, cutoff) {
        psi1 <- vec2psi(i1, s1, effective_inclusion_length, effective_skipping_length)
        psi2 <- vec2psi(i2, s2, effective_inclusion_length, effective_skipping_length)
        iter_cutoff <- 1
        iter_maxrun <- 100
        count <- 0
        previous_sum <- 0
        beta_0 <- mean(psi1)
        beta_1 <- mean(psi2)
        var1 <- 10 * pop_var(psi1 - beta_0)
        var2 <- 10 * pop_var(psi2 - beta_1)
        if (is.na(var1) || var1 <= 0.01) var1 <- 0.01
        if (is.na(var2) || var2 <= 0.01) var2 <- 0.01

        while ((iter_cutoff > 0.01) && (count <= iter_maxrun)) {
            count <- count + 1
            beta_0 <- mean(psi1)
            beta_1 <- mean(psi2)
            if (mean(psi1) > mean(psi2)) {
                xopt <- stats::optim(par = mean(psi2), fn = myfunc_1, gr = myfunc_der_1,
                                     psi1 = psi1, psi2 = psi2, var1 = var1, var2 = var2, cutoff = cutoff,
                                     method = "L-BFGS-B", lower = 0.001, upper = 0.999 - cutoff,
                                     control = optim_ctrl)
                theta2 <- max(min(as.numeric(xopt$par), 1 - cutoff), 0)
                theta1 <- theta2 + cutoff
            } else {
                xopt <- stats::optim(par = mean(psi1), fn = myfunc_2, gr = myfunc_der_2,
                                     psi1 = psi1, psi2 = psi2, var1 = var1, var2 = var2, cutoff = cutoff,
                                     method = "L-BFGS-B", lower = 0.001, upper = 0.999 - cutoff,
                                     control = optim_ctrl)
                theta1 <- max(min(as.numeric(xopt$par), 1 - cutoff), 0)
                theta2 <- theta1 + cutoff
            }
            beta_0 <- theta1
            beta_1 <- theta2
            res1 <- optimize_psi_batch(psi1, i1, s1, beta_0, var1, effective_inclusion_length, effective_skipping_length)
            res2 <- optimize_psi_batch(psi2, i2, s2, beta_1, var2, effective_inclusion_length, effective_skipping_length)
            new_psi1 <- res1$psi
            new_psi2 <- res2$psi
            current_sum <- sum(res1$value) + sum(res2$value)
            psi1 <- new_psi1
            psi2 <- new_psi2
            if (count > 1) iter_cutoff <- abs(previous_sum - current_sum)
            previous_sum <- current_sum
        }
        list(current_sum, list(psi1, psi2, beta_0, beta_1, var1, var2))
    }

    MLE_iteration <- function(i1, i2, s1, s2, effective_inclusion_length, effective_skipping_length) {
        psi1 <- vec2psi(i1, s1, effective_inclusion_length, effective_skipping_length)
        psi2 <- vec2psi(i2, s2, effective_inclusion_length, effective_skipping_length)
        iter_cutoff <- 1
        iter_maxrun <- 100
        count <- 0
        previous_sum <- 0
        beta_0 <- mean(psi1)
        beta_1 <- mean(psi2)
        var1 <- 10 * pop_var(psi1 - beta_0)
        var2 <- 10 * pop_var(psi2 - beta_1)
        if (is.na(var1) || var1 <= 0.01) var1 <- 0.01
        if (is.na(var2) || var2 <= 0.01) var2 <- 0.01

        while ((iter_cutoff > 0.01) && (count <= iter_maxrun)) {
            count <- count + 1
            beta_0 <- mean(psi1)
            beta_1 <- mean(psi2)
            xopt <- stats::optim(par = c(beta_0, beta_1), fn = myfunc_multivar, gr = myfunc_multivar_der,
                                 psi1 = psi1, psi2 = psi2, var1 = var1, var2 = var2,
                                 method = "L-BFGS-B", lower = c(0.01, 0.01), upper = c(0.99, 0.99),
                                 control = optim_ctrl)
            beta_0 <- as.numeric(xopt$par[1])
            beta_1 <- as.numeric(xopt$par[2])
            res1 <- optimize_psi_batch(psi1, i1, s1, beta_0, var1, effective_inclusion_length, effective_skipping_length)
            res2 <- optimize_psi_batch(psi2, i2, s2, beta_1, var2, effective_inclusion_length, effective_skipping_length)
            psi1 <- res1$psi
            psi2 <- res2$psi
            current_sum <- sum(res1$value) + sum(res2$value)
            if (count > 1) iter_cutoff <- abs(previous_sum - current_sum)
            previous_sum <- current_sum
        }
        if (count > iter_maxrun) {
            return(list(current_sum, list(psi1, psi2, 0, 0, var1, var2)))
        }
        list(current_sum, list(psi1, psi2, beta_0, beta_1, var1, var2))
    }

    MLE_marginal_iteration <- function(i1, i2, s1, s2, effective_inclusion_length, effective_skipping_length) {
        psi1 <- vec2psi(i1, s1, effective_inclusion_length, effective_skipping_length)
        psi2 <- vec2psi(i2, s2, effective_inclusion_length, effective_skipping_length)
        iter_cutoff <- 1
        iter_maxrun <- 100
        count <- 0
        previous_sum <- 0
        beta_0 <- mean(psi1)
        beta_1 <- mean(psi2)
        var1 <- 10 * pop_var(psi1 - beta_0)
        var2 <- 10 * pop_var(psi2 - beta_1)
        if (is.na(var1) || var1 <= 0.01) var1 <- 0.01
        if (is.na(var2) || var2 <= 0.01) var2 <- 0.01

        while ((iter_cutoff > 0.01) && (count <= iter_maxrun)) {
            count <- count + 1
            beta_0 <- mean(psi1)
            beta_1 <- mean(psi2)
            xopt <- stats::optim(par = c(beta_0, beta_1), fn = myfunc_multivar, gr = myfunc_multivar_der,
                                 psi1 = psi1, psi2 = psi2, var1 = var1, var2 = var2,
                                 method = "L-BFGS-B", lower = c(0.01, 0.01), upper = c(0.99, 0.99),
                                 control = optim_ctrl)
            beta_0 <- as.numeric(xopt$par[1])
            beta_1 <- as.numeric(xopt$par[2])
            res1 <- optimize_psi_batch(psi1, i1, s1, beta_0, var1, effective_inclusion_length, effective_skipping_length)
            res2 <- optimize_psi_batch(psi2, i2, s2, beta_1, var2, effective_inclusion_length, effective_skipping_length)
            psi1 <- res1$psi
            psi2 <- res2$psi
            current_sum <- sum(res1$value) + sum(res2$value)
            if (count > 1) iter_cutoff <- abs(previous_sum - current_sum)
            previous_sum <- current_sum
        }
        if (count > iter_maxrun) {
            return(list(current_sum, list(psi1, psi2, 0, 0, var1, var2)))
        }

        iter_cutoff <- 1
        iter_maxrun <- 100
        count <- 0
        previous_sum <- 0
        while ((iter_cutoff > 0.01) && (count <= iter_maxrun)) {
            count <- count + 1
            xopt1 <- stats::optim(par = beta_0, fn = myfunc_marginal, gr = myfunc_marginal_der,
                                  I = i1, S = s1, psi = psi1, var = var1,
                                  inclusion_length = effective_inclusion_length, skipping_length = effective_skipping_length,
                                  method = "L-BFGS-B", lower = 0.01, upper = 0.99,
                                  control = optim_ctrl)
            xopt2 <- stats::optim(par = beta_1, fn = myfunc_marginal, gr = myfunc_marginal_der,
                                  I = i2, S = s2, psi = psi2, var = var2,
                                  inclusion_length = effective_inclusion_length, skipping_length = effective_skipping_length,
                                  method = "L-BFGS-B", lower = 0.01, upper = 0.99,
                                  control = optim_ctrl)
            beta_0 <- as.numeric(xopt1$par)
            beta_1 <- as.numeric(xopt2$par)
            res1 <- optimize_psi_batch(psi1, i1, s1, beta_0, var1, effective_inclusion_length, effective_skipping_length)
            res2 <- optimize_psi_batch(psi2, i2, s2, beta_1, var2, effective_inclusion_length, effective_skipping_length)
            psi1 <- res1$psi
            psi2 <- res2$psi
            current_sum <- sum(res1$value) + sum(res2$value)
            if (count > 1) iter_cutoff <- abs(previous_sum - current_sum)
            previous_sum <- current_sum
        }
        if (count > iter_maxrun) {
            return(list(current_sum, list(psi1, psi2, 0, 0, var1, var2)))
        }
        list(current_sum, list(psi1, psi2, beta_0, beta_1, var1, var2))
    }

    MLE_marginal_iteration_constrain <- function(i1, i2, s1, s2, effective_inclusion_length, effective_skipping_length, cutoff) {
        psi1 <- vec2psi(i1, s1, effective_inclusion_length, effective_skipping_length)
        psi2 <- vec2psi(i2, s2, effective_inclusion_length, effective_skipping_length)
        iter_cutoff <- 1
        iter_maxrun <- 100
        count <- 0
        previous_sum <- 0
        beta_0 <- mean(psi1)
        beta_1 <- mean(psi2)
        var1 <- 10 * pop_var(psi1 - beta_0)
        var2 <- 10 * pop_var(psi2 - beta_1)
        if (is.na(var1) || var1 <= 0.01) var1 <- 0.01
        if (is.na(var2) || var2 <= 0.01) var2 <- 0.01

        while ((iter_cutoff > 0.01) && (count <= iter_maxrun)) {
            count <- count + 1
            beta_0 <- mean(psi1)
            beta_1 <- mean(psi2)
            if (mean(psi1) > mean(psi2)) {
                xopt <- stats::optim(par = mean(psi2), fn = myfunc_1, gr = myfunc_der_1,
                                     psi1 = psi1, psi2 = psi2, var1 = var1, var2 = var2, cutoff = cutoff,
                                     method = "L-BFGS-B", lower = 0.001, upper = 0.999 - cutoff,
                                     control = optim_ctrl)
                theta2 <- max(min(as.numeric(xopt$par), 1 - cutoff), 0)
                theta1 <- theta2 + cutoff
            } else {
                xopt <- stats::optim(par = mean(psi1), fn = myfunc_2, gr = myfunc_der_2,
                                     psi1 = psi1, psi2 = psi2, var1 = var1, var2 = var2, cutoff = cutoff,
                                     method = "L-BFGS-B", lower = 0.001, upper = 0.999 - cutoff,
                                     control = optim_ctrl)
                theta1 <- max(min(as.numeric(xopt$par), 1 - cutoff), 0)
                theta2 <- theta1 + cutoff
            }
            beta_0 <- theta1
            beta_1 <- theta2
            res1 <- optimize_psi_batch(psi1, i1, s1, beta_0, var1, effective_inclusion_length, effective_skipping_length)
            res2 <- optimize_psi_batch(psi2, i2, s2, beta_1, var2, effective_inclusion_length, effective_skipping_length)
            psi1 <- res1$psi
            psi2 <- res2$psi
            current_sum <- sum(res1$value) + sum(res2$value)
            if (count > 1) iter_cutoff <- abs(previous_sum - current_sum)
            previous_sum <- current_sum
        }
        if (count > iter_maxrun) {
            return(list(current_sum, list(psi1, psi2, 0, 0, var1, var2)))
        }

        iter_cutoff <- 1
        iter_maxrun <- 100
        count <- 0
        previous_sum <- 0
        while ((iter_cutoff > 0.01) && (count <= iter_maxrun)) {
            count <- count + 1
            if (mean(psi1) > mean(psi2)) {
                xopt <- stats::optim(par = beta_1, fn = myfunc_marginal_1, gr = myfunc_marginal_1_der,
                                     I1 = i1, S1 = s1, psi1 = psi1, var1 = var1,
                                     I2 = i2, S2 = s2, psi2 = psi2, var2 = var2,
                                     inclusion_length = effective_inclusion_length, skipping_length = effective_skipping_length,
                                     cutoff = cutoff, method = "L-BFGS-B", lower = 0.001, upper = 0.999 - cutoff,
                                     control = optim_ctrl)
                beta_1 <- max(min(as.numeric(xopt$par), 1 - cutoff), 0)
                beta_0 <- beta_1 + cutoff
            } else {
                xopt <- stats::optim(par = beta_0, fn = myfunc_marginal_2, gr = myfunc_marginal_2_der,
                                     I1 = i1, S1 = s1, psi1 = psi1, var1 = var1,
                                     I2 = i2, S2 = s2, psi2 = psi2, var2 = var2,
                                     inclusion_length = effective_inclusion_length, skipping_length = effective_skipping_length,
                                     cutoff = cutoff, method = "L-BFGS-B", lower = 0.001, upper = 0.999 - cutoff,
                                     control = optim_ctrl)
                beta_0 <- max(min(as.numeric(xopt$par), 1 - cutoff), 0)
                beta_1 <- beta_0 + cutoff
            }
            res1 <- optimize_psi_batch(psi1, i1, s1, beta_0, var1, effective_inclusion_length, effective_skipping_length)
            res2 <- optimize_psi_batch(psi2, i2, s2, beta_1, var2, effective_inclusion_length, effective_skipping_length)
            psi1 <- res1$psi
            psi2 <- res2$psi
            current_sum <- sum(res1$value) + sum(res2$value)
            if (count > 1) iter_cutoff <- abs(previous_sum - current_sum)
            previous_sum <- current_sum
        }
        if (count > iter_maxrun) {
            return(list(current_sum, list(psi1, psi2, 0, 0, var1, var2)))
        }

        list(current_sum, list(psi1, psi2, beta_0, beta_1, var1, var2))
    }

    likelihood_test <- function(i1, i2, s1, s2, effective_inclusion_length, effective_skipping_length, flag, id, cutoff) {
        if (flag == 0) {
            return(c(1, 1))
        }
        res <- MLE_marginal_iteration(i1, i2, s1, s2, effective_inclusion_length, effective_skipping_length)
        if (abs(res[[2]][[3]] - res[[2]][[4]]) <= cutoff) {
            return(c(1, res[[1]]))
        }
        res_constrain <- MLE_marginal_iteration_constrain(i1, i2, s1, s2, effective_inclusion_length, effective_skipping_length, cutoff)
        stat <- 2 * abs(res_constrain[[1]] - res[[1]])
        pval <- 1 - stats::pchisq(stat, df = 1)
        c(pval)
    }

    MultiProcessorPool <- function(n_original_diff, cutoff) {
        i1 <- n_original_diff[[1]]; i2 <- n_original_diff[[2]]
        s1 <- n_original_diff[[3]]; s2 <- n_original_diff[[4]]
        effective_inclusion_length <- n_original_diff[[5]]
        effective_skipping_length <- n_original_diff[[6]]
        flag <- n_original_diff[[7]]
        id <- n_original_diff[[8]]
        likelihood_test(i1, i2, s1, s2, effective_inclusion_length, effective_skipping_length, flag, id, cutoff)
    }

    vec2float <- function(vec) {
        as.numeric(vec)
    }

    vecAddOne <- function(vec) {
        vec + 1
    }

    vecprod <- function(vec) {
        prod(vec)
    }

    vecadd <- function(vec1, vec2) {
        vec1 + vec2
    }

    vec2remove0psi <- function(inc, skp) {
        keep <- (inc != 0) | (skp != 0)
        list(inc[keep], skp[keep])
    }

    if (exists("vec2psi_cpp", envir = env, inherits = FALSE)) {
        vec2psi <- function(inc, skp, effective_inclusion_length, effective_skipping_length) {
            vec2psi_cpp(inc, skp, effective_inclusion_length, effective_skipping_length)
        }
    } else {
        vec2psi <- function(inc, skp, effective_inclusion_length, effective_skipping_length) {
            psi <- numeric(length(inc))
            for (i in seq_along(inc)) {
                if ((inc[i] + skp[i]) == 0) {
                    psi[i] <- 0.5
                } else {
                    psi[i] <- (inc[i] / effective_inclusion_length) /
                        ((inc[i] / effective_inclusion_length) + (skp[i] / effective_skipping_length))
                }
            }
            psi
        }
    }

    vec210 <- function(vec) {
        ifelse(vec > 0, 1, -1)
    }

    myttest <- function(vec1, vec2) {
        if (length(vec1) == 1 && length(vec2) == 1) {
            stats::t.test(c(vec1[1], vec1[1]), c(vec2[1], vec2[1]))
        } else {
            stats::t.test(vec1, vec2)
        }
    }

    myorder <- function(p, reverse) {
        p_dict <- list()
        for (index in seq_along(p)) {
            key <- as.character(p[index])
            if (is.null(p_dict[[key]])) p_dict[[key]] <- integer(0)
            p_dict[[key]] <- c(p_dict[[key]], index)
        }
        o_index <- sort(as.numeric(names(p_dict)), decreasing = reverse)
        o <- integer(0)
        for (idx in o_index) {
            o <- c(o, p_dict[[as.character(idx)]])
        }
        o
    }

    mycummin <- function(p) {
        cummin(p)
    }

    myFDR <- function(p) {
        lp <- length(p)
        i <- seq(lp, 1)
        o <- myorder(p, TRUE)
        ro <- myorder(o, FALSE)
        p_new <- numeric(length(o))
        for (index in seq_along(o)) {
            p_new[index] <- p[o[index]] * lp / i[index]
        }
        p_mycummin <- mycummin(p_new)
        p_mycummin_new <- pmin(p_mycummin, 1)
        res <- numeric(length(p_mycummin_new))
        for (index in seq_along(p_mycummin_new)) {
            res[index] <- p_mycummin_new[ro[index]]
        }
        res
    }

    adjustPsi <- function(inc, skp, effective_inclusion_length, effective_skipping_length) {
        inc_adj <- inc / effective_inclusion_length
        skp_adj <- skp / effective_skipping_length
        paste(inc_adj / (inc_adj + skp_adj), collapse = ",")
    }

    myand <- function(l1, l2) {
        which((l1 & l2) == TRUE)
    }

    getrMATS_P <- function(list_n_original_diff, cutoff, threads) {
        if (threads > 1) {
            if (.Platform$OS.type == "windows") {
                cl <- parallel::makeCluster(threads)
                on.exit(parallel::stopCluster(cl))
                probability <- parallel::parLapply(cl, list_n_original_diff, MultiProcessorPool, cutoff = cutoff)
            } else {
                probability <- parallel::mclapply(list_n_original_diff, MultiProcessorPool, cutoff = cutoff, mc.cores = threads)
            }
        } else {
            probability <- lapply(list_n_original_diff, MultiProcessorPool, cutoff = cutoff)
        }
        vapply(probability, function(x) x[1], numeric(1))
    }
}, envir = .debks_rmats_env)
