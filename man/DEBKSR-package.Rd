\name{DEBKSR-package}
\alias{DEBKSR-package}
\alias{DEBKSR}
\docType{package}
	\title{DEBKSR: Differential Expression for circRNA (DEBKS) in R}
\description{
An R implementation of the DEBKS "dec" analysis for differential circular RNA
expression between two groups. The package implements constructors and helpers to
create a DEBKS dataset and to run differential analysis (paired or unpaired designs).
}
\details{
Use the `DEBKSDataSet` constructor to assemble junction-count matrices and sample
metadata, then call `DEBKS()` to perform the statistical test. Results (PBSI,
delta PBSI, p-values, and FDR) are stored in the returned object and accessible
via the `results()` accessor. Use `psiMatrix()` or `pbsiMatrix()` to compute
per-sample PSI or PBSI matrices directly from the object.
For paired designs, include a pairing identifier column in `coldata` and pass it
via `paired_id` when calling `DEBKS()`.
}
\author{Yueqi Jin, Ence Yang}
\keyword{package}
\seealso{\code{DEBKSDataSet}, \code{DEBKS}}
\examples{
# Minimal runnable example (toy data)
cjc <- matrix(c(10,5,12,6, 3,2,4,1), nrow = 2, byrow = TRUE) # skipped counts
ljc <- matrix(c(5,6,4,5, 2,1,3,2), nrow = 2, byrow = TRUE)  # inclusion counts
circ_len <- data.frame(chr = c("chr1","chr1"), start = c(100,200), end = c(150,250), length = c(500,400))
coldata <- data.frame(group = c("A","A","B","B"))
dbs <- DEBKSDataSet(cjc = cjc, ljc = ljc, circ_len = circ_len, coldata = coldata)
# Run a quick DEC analysis (single-thread, return all results)
res_obj <- DEBKS(dbs, paired = FALSE, threads = 1, cutoff = 1,
				 group_var = "group", reference_group = "A", observe_group = "B")
head(results(res_obj))
}
