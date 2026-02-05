# DEBKSR

DEBKSR is an R package implementing the DEBKS "dec" analysis for detecting
differentially expressed circular RNAs between two groups. It supports paired
and unpaired designs and exposes a small API for building datasets and running
the analysis.

## Installation

First ensure prerequisites are installed (R >= 4.x and Rtools on Windows). The
package uses `Rcpp`; install it from CRAN if you don't have it:

```r
install.packages("Rcpp", repos = "https://cloud.r-project.org")
```

To build and install the package from github:

```r
devtools::install_github("YueqiJin/DEBKSR")
```


## Quick usage (package-style)

After installation, load the package and use the constructor and functions:

```r
library(DEBKSR)

# toy example
cjc <- matrix(c(10,5,12,6, 3,2,4,1), nrow = 2, byrow = TRUE)
ljc <- matrix(c(5,6,4,5, 2,1,3,2), nrow = 2, byrow = TRUE)
circ_len <- data.frame(chr = c("chr1","chr1"), start = c(100,200), end = c(150,250), length = c(500,400))
coldata <- data.frame(group = c("A","A","B","B"))

dbs <- DEBKSDataSet(cjc = cjc, ljc = ljc, circ_len = circ_len, coldata = coldata)
res_obj <- DEBKS(dbs, paired = FALSE, threads = 1, cutoff = 1,
				 group_var = "group", reference_group = "A", observe_group = "B")
head(results(res_obj))

# For paired designs, add a column like `pair_id` in coldata and pass
# paired = TRUE with paired_id = "pair_id".
```

## Output

The main result table contains coordinates, junction counts, PBSI estimates,
delta PBSI, p-values and FDR-adjusted p-values. Use `results()` to access the
table from the returned `DEBKS` object.

## Citation

If you use DEBKSR in your research, please cite the following publication:
Liu Z, Ding H, She J, Chen C, Zhang W, Yang E. DEBKS: A Tool to Detect Differentially Expressed Circular RNAs. Genomics Proteomics Bioinformatics. 2022;20(3):549-556. doi:10.1016/j.gpb.2021.01.003