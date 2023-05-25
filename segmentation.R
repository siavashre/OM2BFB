library('DNAcopy')
args = commandArgs(trailingOnly=TRUE)
data <- read.csv(args[1])
genomedata <- data[[4]]
chrom <-rep(1,length(genomedata))
w <- data[[5]]
maploc <- c(1:length(genomedata))
ans <- segment(CNA(genomedata, chrom, maploc), weights=w)
write.table(ans[2], args[2], append = FALSE, sep = "\t", dec = ".",
            row.names = TRUE, col.names = TRUE)
