DEBKSDataSet <- function(cjc, ljc, circ_len, coldata, cjc2 = NULL, circ_len2 = NULL) {
    if (is.data.frame(cjc)) cjc <- as.matrix(cjc)
    if (is.data.frame(ljc)) ljc <- as.matrix(ljc)
    if (!is.matrix(cjc) || !is.matrix(ljc)) {
        stop("cjc and ljc must be matrices or data frames.")
    }
    if (nrow(cjc) != nrow(ljc)) {
        stop("cjc and ljc must have the same number of rows.")
    }
    if (ncol(cjc) != ncol(ljc)) {
        stop("cjc and ljc must have the same number of columns.")
    }
    if (!is.data.frame(circ_len)) {
        stop("circ_len must be a data.frame with columns chr, start, end, length.")
    }
    required_cols <- c("chr", "start", "end", "length")
    if (!all(required_cols %in% colnames(circ_len))) {
        stop("circ_len must have columns: chr, start, end, length.")
    }
    if (nrow(circ_len) != nrow(cjc)) {
        stop("circ_len row count must match cjc/ljc row count.")
    }
    if (!is.data.frame(coldata)) {
        stop("coldata must be a data.frame.")
    }
    if (nrow(coldata) != ncol(cjc)) {
        stop("coldata rows must match sample columns in cjc/ljc.")
    }

    if (!is.null(cjc2)) {
        if (is.data.frame(cjc2)) cjc2 <- as.matrix(cjc2)
        if (!is.matrix(cjc2)) stop("cjc2 must be a matrix or data.frame.")
        if (nrow(cjc2) != nrow(cjc)) stop("cjc2 row count must match cjc.")
        if (ncol(cjc2) != ncol(cjc)) stop("cjc2 column count must match cjc.")
        if (is.null(circ_len2)) stop("circ_len2 is required when cjc2 is provided.")
        if (!is.data.frame(circ_len2)) stop("circ_len2 must be a data.frame.")
        if (!all(required_cols %in% colnames(circ_len2))) {
            stop("circ_len2 must have columns: chr, start, end, length.")
        }
        if (nrow(circ_len2) != nrow(cjc)) {
            stop("circ_len2 row count must match cjc.")
        }
    }

    obj <- list(
        cjc = cjc,
        ljc = ljc,
        circ_len = circ_len,
        coldata = coldata,
        cjc2 = cjc2,
        circ_len2 = circ_len2,
        results = NULL
    )
    class(obj) <- "DEBKS"
    obj
}

sizeFactors <- function(x) {
    if (!inherits(x, "DEBKS")) stop("x must be a DEBKS object.")
    rep(1, ncol(x$cjc))
}

psiMatrix <- function(x) {
    if (!inherits(x, "DEBKS")) stop("x must be a DEBKS object.")
    inc <- x$ljc
    skp <- x$cjc
    inc_adj <- sweep(inc, 1, x$circ_len$length, "/")
    skp_adj <- sweep(skp, 1, x$circ_len$length, "/")
    inc_adj / (inc_adj + skp_adj)
}

results <- function(x) {
    if (!inherits(x, "DEBKS")) stop("x must be a DEBKS object.")
    x$results
}
