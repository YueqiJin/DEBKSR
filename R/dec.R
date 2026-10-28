DEBKS <- function(dbs,
                  long_read = FALSE,
                  paired = FALSE,
                  group_var = NULL,
                  reference_group = NULL,
                  observe_group = NULL,
                  paired_id = NULL,
                  circRNA_filter = 0,
                  sample_filter = 0,
                  threads = 4,
                  cutoff = 0.05,
                  read_len = 150,
                  hangover_len = 6,
                  engine = c("rmats", "bb"),
                  output = NULL) {
    if (!inherits(dbs, "DEBKS")) stop("dbs must be a DEBKS object.")
    engine <- match.arg(engine)

    threads <- as.integer(threads)
    cutoff <- min(as.numeric(cutoff), 1.00)
    read_len <- as.integer(read_len)
    hangover_len <- as.integer(hangover_len)
    circRNA_filter <- max(as.integer(circRNA_filter), 0)
    sample_filter <- max(as.integer(sample_filter), 0)

    cjc <- dbs$cjc
    ljc <- dbs$ljc
    circ_len <- dbs$circ_len

    if (is.null(group_var)) {
        stop("group_var (coldata) must be provided.")
    }

    if (!is.null(group_var)) {
        if (!is.data.frame(dbs$coldata)) stop("dbs$coldata must be a data.frame.")
        if (!group_var %in% colnames(dbs$coldata)) {
            stop("group_var must be a column name in coldata.")
        }
        group_vals <- dbs$coldata[[group_var]]
        group_levels <- unique(group_vals)
        if (is.null(reference_group) || is.null(observe_group)) {
            if (length(group_levels) != 2) {
                stop("When group_var has more than two groups, reference_group and observe_group are required.")
            }
            reference_group <- group_levels[1]
            observe_group <- group_levels[2]
        }
        if (!reference_group %in% group_levels) stop("reference_group not found in group_var.")
        if (!observe_group %in% group_levels) stop("observe_group not found in group_var.")
        if (reference_group == observe_group) stop("reference_group and observe_group must differ.")

        idx_ref <- which(group_vals == reference_group)
        idx_obs <- which(group_vals == observe_group)
        if (length(idx_ref) == 0 || length(idx_obs) == 0) {
            stop("reference_group/observe_group selection returned zero samples.")
        }

        if (paired) {
            if (is.null(paired_id)) stop("paired_id must be provided for paired analysis.")
            if (!paired_id %in% colnames(dbs$coldata)) {
                stop("paired_id must be a column name in coldata.")
            }
            pair_vals <- dbs$coldata[[paired_id]]
            ref_pairs <- pair_vals[idx_ref]
            obs_pairs <- pair_vals[idx_obs]
            if (any(is.na(ref_pairs)) || any(is.na(obs_pairs))) {
                stop("paired_id contains NA values for selected samples.")
            }
            if (any(duplicated(ref_pairs)) || any(duplicated(obs_pairs))) {
                stop("paired_id must be unique within each group.")
            }
            common_pairs <- ref_pairs[ref_pairs %in% obs_pairs]
            if (length(common_pairs) == 0) {
                stop("No overlapping paired_id values between reference_group and observe_group.")
            }
            if (length(common_pairs) != length(ref_pairs) || length(common_pairs) != length(obs_pairs)) {
                stop("Paired analysis requires one sample per group for each paired_id.")
            }
            ref_ordered <- idx_ref[match(common_pairs, ref_pairs)]
            obs_ordered <- idx_obs[match(common_pairs, obs_pairs)]
            new_order <- c(ref_ordered, obs_ordered)
        } else {
            new_order <- c(idx_ref, idx_obs)
        }

        cjc <- cjc[, new_order, drop = FALSE]
        ljc <- ljc[, new_order, drop = FALSE]
        if (!is.null(dbs$cjc2)) {
            dbs$cjc2 <- dbs$cjc2[, new_order, drop = FALSE]
        }
        dbs$coldata <- dbs$coldata[new_order, , drop = FALSE]
    }

    if (paired) {
        rmats_env <- .debks_rmats_paired_env
        getrMATS_P <- rmats_env$getrMATS_P
        vec2float <- rmats_env$vec2float
        vecAddOne <- rmats_env$vecAddOne
        adjustPsi <- rmats_env$adjustPsi
        myFDR <- rmats_env$myFDR
        groupPair <- TRUE
    } else {
        rmats_env <- .debks_rmats_env
        getrMATS_P <- rmats_env$getrMATS_P
        vec2float <- rmats_env$vec2float
        vecAddOne <- rmats_env$vecAddOne
        adjustPsi <- rmats_env$adjustPsi
        myFDR <- rmats_env$myFDR
        groupPair <- FALSE
    }

    sampleNum <- ncol(cjc)
    group1_num <- length(idx_ref)
    group2_num <- length(idx_obs)
    if (group1_num < 2 || group2_num < 2) stop("Require replicates in both groups.")

    adjustReadLength <- read_len - hangover_len
    if (adjustReadLength <= 0) stop("read_len must be larger than hangover_len.")

    adjust_inc_len <- rep(4 * adjustReadLength, nrow(cjc))
    adjust_skp_len <- rep(2 * adjustReadLength, nrow(cjc))

    if (!is.null(circ_len)) {
        adj <- getCircLen_from_df(circ_len, adjustReadLength)
        adjust_inc_len <- adj$adjust_inc_len
        adjust_skp_len <- adj$adjust_skp_len
    }

    if (!is.null(dbs$cjc2)) {
        if (!is.null(dbs$circ_len2)) {
            adj2 <- getCircLen_from_df(dbs$circ_len2, adjustReadLength)
            adjust_inc_len <- adj2$adjust_skp_len
        } else {
            adjust_inc_len <- rep(2 * adjustReadLength, nrow(cjc))
        }
    }

    if (long_read) {
        adjust_inc_len <- adjust_skp_len
    }

    if (circRNA_filter > 0) {
        passID <- filterCirc_matrix(cjc, circRNA_filter)
        if (length(passID) > 0) {
            cjc <- cjc[passID, , drop = FALSE]
            ljc <- ljc[passID, , drop = FALSE]
            circ_len <- circ_len[passID, , drop = FALSE]
            adjust_inc_len <- adjust_inc_len[passID]
            adjust_skp_len <- adjust_skp_len[passID]
            if (!is.null(dbs$cjc2)) {
                dbs$cjc2 <- dbs$cjc2[passID, , drop = FALSE]
            }
        } else {
            stop("Filter kept zero circRNAs for dec.")
        }
    }

    if (engine == "bb" && paired) {
        stop("engine='bb' currently supports unpaired designs only.")
    }

    res <- stat_r(cjc, ljc, group1_num, sampleNum, adjust_inc_len, adjust_skp_len,
                 cutoff, groupPair, threads, output, !is.null(dbs$cjc2), sample_filter,
                 read_len, circ_len, getrMATS_P, vec2float, vecAddOne, adjustPsi, myFDR,
                 engine = engine)
    dbs$results <- res$matStat
    dbs
}

filterCirc_matrix <- function(mat, filterNum) {
    rowSum <- rowSums(mat)
    which(rowSum >= filterNum)
}

getCircLen_from_df <- function(circ_len, adjustReadLength) {
    if (any(circ_len$length <= 0)) stop("circ_len contains zero length.")
    circLen <- as.integer(circ_len$length)
    adjust_inc_len <- integer(length(circLen))
    adjust_skp_len <- integer(length(circLen))
    for (i in seq_along(circLen)) {
        if (circLen[i] < 2 * adjustReadLength) {
            adjust_inc_len[i] <- 2 * adjustReadLength + circLen[i]
            adjust_skp_len[i] <- circLen[i]
        } else {
            adjust_inc_len[i] <- 4 * adjustReadLength
            adjust_skp_len[i] <- 2 * adjustReadLength
        }
    }
    list(adjust_inc_len = adjust_inc_len, adjust_skp_len = adjust_skp_len)
}

psi2PBSI <- function(psi) {
    pbsi <- character(length(psi))
    for (i in seq_along(psi)) {
        eachPsi <- strsplit(psi[i], ",", fixed = TRUE)[[1]]
        eachPbsi <- vapply(eachPsi, function(j) {
            if (j != "") {
                as.character(1 - as.numeric(j))
            } else {
                ""
            }
        }, character(1))
        pbsi[i] <- paste(eachPbsi, collapse = ",")
    }
    pbsi
}

getMean <- function(x) {
    if (length(x) == 0) return(0)
    vals <- suppressWarnings(as.numeric(x))
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0) return(0)
    mean(vals)
}

get_delta <- function(PBSI1, PBSI2) {
    delta_pbsi <- numeric(length(PBSI1))
    for (i in seq_along(PBSI1)) {
        eachPbsi1_mean <- getMean(strsplit(PBSI1[i], ",", fixed = TRUE)[[1]])
        eachPbsi2_mean <- getMean(strsplit(PBSI2[i], ",", fixed = TRUE)[[1]])
        delta_pbsi[i] <- eachPbsi1_mean - eachPbsi2_mean
    }
    delta_pbsi
}

vec_remove_na_pair <- function(inc1, skp1, inc2, skp2) {
    keep <- (inc1 + skp1) > 0 & (inc2 + skp2) > 0
    list(inc1[keep], skp1[keep], inc2[keep], skp2[keep])
}

vec_remove_na <- function(inc, skp) {
    keep <- (inc + skp) > 0
    list(inc[keep], skp[keep])
}

vec_remove_low_pair <- function(inc1, skp1, inc2, skp2, num) {
    keep <- (skp1 >= num) & (skp2 >= num)
    list(inc1[keep], skp1[keep], inc2[keep], skp2[keep])
}

vec_remove_low <- function(inc, skp, num) {
    keep <- skp >= num
    list(inc[keep], skp[keep])
}

checkInf <- function(x) {
    x[is.infinite(x)] <- NaN
    x
}

stat_r <- function(circMat, linearMat, group1_num, sampleNum,
                   adjust_inc_len, adjust_skp_len, cutoff, groupPair,
                   threads, output, isc2, filterSample, readLength,
                   circ_len, getrMATS_P, vec2float, vecAddOne, adjustPsi, myFDR,
                   engine = "rmats") {
    list_n_original_diff <- list()
    psi_list_1 <- character(nrow(circMat))
    psi_list_2 <- character(nrow(circMat))

    for (i in seq_len(nrow(circMat))) {
        inc1 <- linearMat[i, 1:group1_num]
        inc2 <- linearMat[i, (group1_num + 1):sampleNum]
        skp1 <- circMat[i, 1:group1_num]
        skp2 <- circMat[i, (group1_num + 1):sampleNum]
        effective_inclusion_length <- adjust_inc_len[i]
        effective_skipping_length <- adjust_skp_len[i]
        inc1 <- vec2float(inc1); skp1 <- vec2float(skp1)
        inc2 <- vec2float(inc2); skp2 <- vec2float(skp2)
        inc1 <- checkInf(inc1); skp1 <- checkInf(skp1)
        inc2 <- checkInf(inc2); skp2 <- checkInf(skp2)
        if (effective_inclusion_length <= 0) {
            effective_inclusion_length <- readLength
        }
        if (groupPair) {
            temp <- vec_remove_na_pair(inc1, skp1, inc2, skp2)
            inc1_nona <- temp[[1]]; skp1_nona <- temp[[2]]
            inc2_nona <- temp[[3]]; skp2_nona <- temp[[4]]
            temp <- vec_remove_low_pair(inc1_nona, skp1_nona, inc2_nona, skp2_nona, filterSample)
            inc1_nona <- temp[[1]]; skp1_nona <- temp[[2]]
            inc2_nona <- temp[[3]]; skp2_nona <- temp[[4]]
            if (length(inc1_nona) == 0 || length(inc2_nona) == 0) {
                psi_list_1[i] <- ""
                psi_list_2[i] <- ""
                list_n_original_diff[[length(list_n_original_diff) + 1]] <- list(inc1, inc2, skp1, skp2,
                                                                                effective_inclusion_length, effective_skipping_length, 0, i)
            } else {
                psi_list_1[i] <- adjustPsi(inc1, skp1, effective_inclusion_length, effective_skipping_length)
                psi_list_2[i] <- adjustPsi(inc2, skp2, effective_inclusion_length, effective_skipping_length)
                inc1 <- inc1_nona; skp1 <- skp1_nona
                inc2 <- inc2_nona; skp2 <- skp2_nona
                inc1 <- vecAddOne(inc1); skp1 <- vecAddOne(skp1)
                inc2 <- vecAddOne(inc2); skp2 <- vecAddOne(skp2)
                if (length(inc1) == 1 || length(inc2) == 1) {
                    list_n_original_diff[[length(list_n_original_diff) + 1]] <- list(inc1, inc2, skp1, skp2,
                                                                                    effective_inclusion_length, effective_skipping_length, 0, i)
                } else {
                    list_n_original_diff[[length(list_n_original_diff) + 1]] <- list(inc1, inc2, skp1, skp2,
                                                                                    effective_inclusion_length, effective_skipping_length, 1, i)
                }
            }
        } else {
            temp1 <- vec_remove_na(inc1, skp1)
            temp2 <- vec_remove_na(inc2, skp2)
            inc1_nona <- temp1[[1]]; skp1_nona <- temp1[[2]]
            inc2_nona <- temp2[[1]]; skp2_nona <- temp2[[2]]
            temp1 <- vec_remove_low(inc1_nona, skp1_nona, filterSample)
            temp2 <- vec_remove_low(inc2_nona, skp2_nona, filterSample)
            inc1_nona <- temp1[[1]]; skp1_nona <- temp1[[2]]
            inc2_nona <- temp2[[1]]; skp2_nona <- temp2[[2]]
            if (length(inc1_nona) == 0 || length(inc2_nona) == 0 || effective_inclusion_length == 0) {
                psi_list_1[i] <- ""
                psi_list_2[i] <- ""
                list_n_original_diff[[length(list_n_original_diff) + 1]] <- list(inc1, inc2, skp1, skp2,
                                                                                effective_inclusion_length, effective_skipping_length, 0, i)
            } else {
                psi_list_1[i] <- adjustPsi(inc1, skp1, effective_inclusion_length, effective_skipping_length)
                psi_list_2[i] <- adjustPsi(inc2, skp2, effective_inclusion_length, effective_skipping_length)
                inc1 <- inc1_nona; skp1 <- skp1_nona
                inc2 <- inc2_nona; skp2 <- skp2_nona
                inc1 <- vecAddOne(inc1); skp1 <- vecAddOne(skp1)
                inc2 <- vecAddOne(inc2); skp2 <- vecAddOne(skp2)
                list_n_original_diff[[length(list_n_original_diff) + 1]] <- list(inc1, inc2, skp1, skp2,
                                                                                effective_inclusion_length, effective_skipping_length, 1, i)
            }
        }
    }

    if (engine == "bb") {
        P <- debks_bb_core$bb_getP(list_n_original_diff, cutoff, threads)
    } else {
        P <- getrMATS_P(list_n_original_diff, cutoff, threads)
    }
    FDR <- myFDR(P)
    pbsi1 <- psi2PBSI(psi_list_1)
    pbsi2 <- psi2PBSI(psi_list_2)
    deltPBSI <- get_delta(pbsi1, pbsi2)

    cjc1 <- apply(circMat[, 1:group1_num, drop = FALSE], 1, function(x) paste(x, collapse = ","))
    cjc2 <- apply(circMat[, (group1_num + 1):sampleNum, drop = FALSE], 1, function(x) paste(x, collapse = ","))
    ljc1 <- apply(linearMat[, 1:group1_num, drop = FALSE], 1, function(x) paste(x, collapse = ","))
    ljc2 <- apply(linearMat[, (group1_num + 1):sampleNum, drop = FALSE], 1, function(x) paste(x, collapse = ","))

    if (isc2) {
        matStat <- data.frame(
            chr = circ_len$chr,
            start = circ_len$start,
            end = circ_len$end,
            cjc_1 = cjc1,
            cjc_2 = cjc2,
            cjc2_1 = ljc1,
            cjc2_2 = ljc2,
            adj_cjc_len = adjust_skp_len,
            adj_cjc2_len = adjust_inc_len,
            pbsi1 = pbsi1,
            pbsi2 = pbsi2,
            delta_pbsi = deltPBSI,
            P = P,
            FDR = FDR
        )
    } else {
        matStat <- data.frame(
            chr = circ_len$chr,
            start = circ_len$start,
            end = circ_len$end,
            cjc_1 = cjc1,
            cjc_2 = cjc2,
            ljc_1 = ljc1,
            ljc_2 = ljc2,
            adj_cjc_len = adjust_skp_len,
            adj_ljc_len = adjust_inc_len,
            pbsi1 = pbsi1,
            pbsi2 = pbsi2,
            delta_pbsi = deltPBSI,
            P = P,
            FDR = FDR
        )
    }

    if (!is.null(output)) {
        utils::write.table(matStat, file = output, sep = "\t", quote = FALSE, row.names = FALSE)
    }

    list(matStat = matStat, decNum = sum(FDR < 0.05, na.rm = TRUE))
}
