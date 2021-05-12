#! /usr/bin/env Rscript-4

# This script takes any number of bigwig files as input and quantile normalizes the count distributions to the one of the first bigwig in the argument list.
# The distributions are scaled to the range [0,1].
# Values will be capped at a particular percentile, such that extreme outliers (often the case for ATAC-seq) are be capped to that value.
# The numeric target_column argument specifies which samplet to use as a target distribution to which the others are quantile normalized.

# parse arguments
args <- commandArgs(trailingOnly=T)
if (length(args) < 6) stop('Usage: Rscript quantile_normalize_bws.R <genome.sizes> <outdir> <percentile-cap [0-1]> <target_column> <file_1.bw> ... <file_n.bw>')

sizes_file <- args[1]
outdir <- normalizePath(args[2])
percentile_cap <- as.numeric(args[3])
target_col <- as.numeric(args[4])
bwfiles <- file.path(outdir, args[-(1:4)])
cat('Loading packages\n')
packages <- c('preprocessCore', 'rtracklayer', 'tictoc')
for (pkg in packages) suppressMessages(library(pkg, character.only=T))

tic()
cat('Tiling genome\n')
sizes <- read.table(sizes_file, col.names=c('chrom','size'))
sizes$size <- as.integer(sizes$size / 10) * 10 # crop chromosome lengths to a multiple of the binsize (10)
sizes <- sizes[!grepl('_|M', sizes$chrom),]
seq_lengths <- sizes$size
names(seq_lengths) <- sizes$chr
genome_gr <- unlist(tile(GRanges(paste0(sizes$chrom, ':1-', sizes$size)), width=10))
seqlengths(genome_gr) <- seq_lengths
toc()

tic()
cat('Reading bigwigs\n')
print(bwfiles)
gr_list <- mclapply(bwfiles, function(bw) import.bw(bw, which=genome_gr), mc.cores=min(10,length(bwfiles)))
print(gr_list)
df <- sapply(gr_list, function(gr) rep(gr$score, ceiling(width(gr)/10))) # GRanges contain multi-bin regions with the same score. unfold them by repeating their scores with the number of bins their regions span
toc()

tic()
cat('Quantile normalizing counts\n')
df_qnorm <- normalize.quantiles.use.target(x=df, target=df[,target_col])
df_maxnorm <- df_qnorm / quantile(df_qnorm, percentile_cap)
df_maxnorm[df_maxnorm > 1] <- 1
toc()

tic()
cat('Writing to files\n')
for (i in seq_len(ncol(df_maxnorm))) {
    scores <- c(df_maxnorm[,i], rep(0, length(genome_gr)-length(df_maxnorm[,i]))) # in some cases the df has one element less than genome_gr, can't figure out why. this extends the scores vector to the same length as the genome_gr
    genome_gr$score <- scores # use genome_gr here again, not genome_gr_narrow, as we only used the narrowed regions for importing the bigwigs
    outfile <- file.path(outdir, gsub('\\.bw', '\\.qnorm\\.bw', basename(bwfiles[[i]])))
    cat(sprintf('%s\n', outfile))
    export.bw(genome_gr, outfile)
}
toc()
cat('Done\n')
