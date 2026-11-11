#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(GenomicRanges)
})

opt = parse_args(OptionParser(option_list = list(
  make_option("--sample",              type = "character"),
  make_option("--consensus_segments",  type = "character"),
  make_option("--regions",             type = "character"),
  make_option("--min_callers",         type = "integer", default = 2L),
  make_option("--outdir",              type = "character", default = ".")
)))

stopifnot(!is.null(opt$sample), !is.null(opt$consensus_segments), !is.null(opt$regions))

reg = readr::read_tsv(opt$regions, show_col_types = FALSE)
seg = readr::read_tsv(opt$consensus_segments, show_col_types = FALSE)
if ("role" %in% names(seg)) seg = seg[seg$role != "cycle", ]

keep = reg$ecdna_region_id[reg$n_callers >= opt$min_callers]
message("Regions with >= ", opt$min_callers, " callers: ", length(keep), " of ", nrow(reg))

for (id in keep) {
  s = dplyr::filter(seg, ecdna_region_id == id)
  gr = GenomicRanges::reduce(GenomicRanges::GRanges(s$chr, IRanges::IRanges(s$start, s$end)))
  bed = data.frame(chr = as.character(GenomicRanges::seqnames(gr)), start = GenomicRanges::start(gr) - 1L, end = GenomicRanges::end(gr))
  out_file = file.path(opt$outdir, paste0(opt$sample, "_region", id, ".bed"))
  write.table(bed, out_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  message("Region ", id, ": ", nrow(bed), " intervals, ", sum(bed$end - bed$start), " bp -> ", out_file)
}
