# CalcCPM.R <count_matrix.tsv> <out_cpm.tsv>
# TMM-normalized CPM via edgeR. Reads the matrix written by ReadCount.sh.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)

suppressPackageStartupMessages(library(edgeR))

x <- read.delim(args[1], row.names = 1, check.names = FALSE)
stopifnot(ncol(x) > 0, all(sapply(x, is.numeric)))

y <- calcNormFactors(DGEList(counts = x))
write.table(cpm(y, normalized.lib.sizes = TRUE), args[2],
            sep = "\t", quote = FALSE, col.names = NA)

cat(sprintf("[CPM ] %d genes x %d samples -> %s\n", nrow(x), ncol(x), args[2]))
