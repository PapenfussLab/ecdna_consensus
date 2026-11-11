#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
})

opt = parse_args(OptionParser(option_list = list(
  make_option("--sample",             type = "character"),
  make_option("--consensus_segments", type = "character"),
  make_option("--svclone",            type = "character"),
  make_option("--tol",                type = "integer", default = 50L),
  make_option("--outdir",             type = "character", default = ".")
)))

stopifnot(!is.null(opt$sample), !is.null(opt$consensus_segments), !is.null(opt$svclone))

seg = readr::read_tsv(opt$consensus_segments, show_col_types = FALSE)
if ("role" %in% names(seg)) seg = seg[seg$role != "footprint", ]
sv = read.delim(opt$svclone)
seg = seg[order(seg$call_id, seg$segment_order), ]

seg$exit_pos = ifelse(seg$strand == "+", seg$end, seg$start)
seg$exit_dir = seg$strand
seg$entry_pos = ifelse(seg$strand == "+", seg$start, seg$end)
seg$entry_dir = ifelse(seg$strand == "+", "-", "+")

n = nrow(seg)
nxt = seq_len(n) + 1
nxt[nxt > n] = 1
first_row = match(seg$call_id, seg$call_id)
is_last = seg$call_id[nxt] != seg$call_id
nxt[is_last] = first_row[is_last]

jn = data.frame(
  sample_id = rep(opt$sample, n),
  ecdna_region_id = seg$ecdna_region_id,
  caller = seg$caller,
  call_id = seg$call_id,
  chrA = seg$chr,
  posA = seg$exit_pos,
  dirA = seg$exit_dir,
  chrB = seg$chr[nxt],
  posB = seg$entry_pos[nxt],
  dirB = seg$entry_dir[nxt]
)
jn$contiguous = jn$chrA == jn$chrB & abs(jn$posA - jn$posB) <= 1 & jn$dirA != jn$dirB

jn$sv_row = rep(NA_integer_, nrow(jn))
jn$dist = rep(Inf, nrow(jn))

for (i in seq_len(nrow(jn))) {
  j = jn[i, ]

  fw = sv$chr1 == j$chrA & sv$chr2 == j$chrB & sv$dir1 == j$dirA & sv$dir2 == j$dirB
  d_fw = pmax(abs(sv$pos1 - j$posA), abs(sv$pos2 - j$posB))
  d_fw[!fw] = Inf

  rv = sv$chr1 == j$chrB & sv$chr2 == j$chrA & sv$dir1 == j$dirB & sv$dir2 == j$dirA
  d_rv = pmax(abs(sv$pos1 - j$posB), abs(sv$pos2 - j$posA))
  d_rv[!rv] = Inf

  d = pmin(d_fw, d_rv)
  k = which.min(d)
  if (length(k) == 1) {
    jn$sv_row[i] = k
    jn$dist[i] = d[k]
  }
}

jn$matched = !jn$contiguous & is.finite(jn$dist) & jn$dist <= opt$tol

jn$svclone_event = as.character(sv$gridss_event[jn$sv_row])
jn$ccf = sv$ccf[jn$sv_row]
jn$cluster = sv$cluster[jn$sv_row]
jn$clonal = sv$clonal[jn$sv_row]

jn$svclone_event[!jn$matched] = NA
jn$ccf[!jn$matched] = NA
jn$cluster[!jn$matched] = NA
jn$clonal[!jn$matched] = NA
jn$dist[is.infinite(jn$dist)] = NA

jn = jn[, c("sample_id", "ecdna_region_id", "caller", "call_id", "chrA", "posA", "dirA", "chrB", "posB", "dirB",
            "contiguous", "dist", "matched", "svclone_event", "ccf", "cluster", "clonal")]

summ = jn[!duplicated(jn$call_id), c("sample_id", "ecdna_region_id", "caller", "call_id")]
summ$n_junctions = rep(0L, nrow(summ))
summ$n_matched = rep(0L, nrow(summ))
summ$median_ccf = rep(NA_real_, nrow(summ))

for (i in seq_len(nrow(summ))) {
  rows = jn$call_id == summ$call_id[i]
  summ$n_junctions[i] = sum(!jn$contiguous[rows])
  summ$n_matched[i] = sum(jn$matched[rows])
  summ$median_ccf[i] = median(jn$ccf[rows], na.rm = TRUE)
}

readr::write_tsv(jn, file.path(opt$outdir, paste0(opt$sample, "_junction_ccf.tsv")), na = "")
readr::write_tsv(summ, file.path(opt$outdir, paste0(opt$sample, "_call_ccf.tsv")), na = "")

message("Junctions: ", nrow(jn), " (non-contiguous: ", sum(!jn$contiguous), ", matched to SVclone: ", sum(jn$matched), ")")
