#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(GenomicRanges)
  library(readr)
})

opt = parse_args(OptionParser(option_list = list(
  make_option("--sample",      type = "character"),
  make_option("--segments",    type = "character"),
  make_option("--min_overlap", type = "double", default = 0.5),
  make_option("--nested_frac", type = "double", default = 0.9),
  make_option("--outdir",      type = "character", default = ".")
)))

stopifnot(!is.null(opt$sample), !is.null(opt$segments))

shared_bp = function(gr1, gr2) {
  sum(IRanges::width(suppressWarnings(GenomicRanges::intersect(gr1, gr2, ignore.strand = TRUE))))
}

cluster_calls_complete = function(segments_df, min_overlap, nested_frac) {
  stopifnot(min_overlap > 0, min_overlap <= 1)
  calls = segments_df %>%
    dplyr::group_by(sample_id, caller, call_id) %>%
    dplyr::summarise(
      chr_set = paste(sort(unique(chr)), collapse = ","),
      role = dplyr::first(role),
      gr = list(GenomicRanges::GRanges(chr, IRanges::IRanges(start, end))),
      .groups = "drop"
    ) %>%
    dplyr::mutate(is_interchr = grepl(",", chr_set), ecdna_region_id = NA_integer_)

  n = nrow(calls)
  grs = lapply(calls$gr, GenomicRanges::reduce, ignore.strand = TRUE)
  size = sapply(grs, function(g) sum(IRanges::width(g)))

  inter = matrix(0, n, n)
  if (n > 1) {
    for (a in seq_len(n - 1)) {
      for (b in (a + 1):n) {
        inter[a, b] = inter[b, a] = shared_bp(grs[[a]], grs[[b]])
      }
    }
  }
  larger = outer(size, size, pmin)
  score = ifelse(larger > 0, inter / larger, 0)
  diag(score) = 1

  if (n == 1) {
    calls$ecdna_region_id = 1L
  } else {
    hc = hclust(as.dist(1 - score), method = "average")
    cl = cutree(hc, h = 1 - min_overlap)
    calls$ecdna_region_id = match(cl, unique(cl))
  }

  calls$nested_in_region = NA_integer_
  for (i in seq_len(n)) {
    cand = which(size > size[i] & inter[i, ] / size[i] >= nested_frac &
                   calls$ecdna_region_id != calls$ecdna_region_id[i])
    if (length(cand) > 0) calls$nested_in_region[i] = calls$ecdna_region_id[cand[which.max(size[cand])]]
  }

  idx = which(upper.tri(inter) & inter > 0, arr.ind = TRUE)
  pairs = data.frame(
    sample_id = calls$sample_id[idx[, 1]],
    call_id_a = calls$call_id[idx[, 1]],
    caller_a = calls$caller[idx[, 1]],
    call_id_b = calls$call_id[idx[, 2]],
    caller_b = calls$caller[idx[, 2]],
    region_a = calls$ecdna_region_id[idx[, 1]],
    region_b = calls$ecdna_region_id[idx[, 2]],
    size_a_bp = size[idx[, 1]],
    size_b_bp = size[idx[, 2]],
    overlap_bp = inter[idx],
    frac_a_covered = inter[idx] / size[idx[, 1]],
    frac_b_covered = inter[idx] / size[idx[, 2]],
    overlap_score = score[idx],
    jaccard = inter[idx] / (size[idx[, 1]] + size[idx[, 2]] - inter[idx]),
    stringsAsFactors = FALSE
  )
  pairs$same_region = pairs$region_a == pairs$region_b
  pairs = pairs[, c("sample_id", "call_id_a", "caller_a", "call_id_b", "caller_b",
                    "region_a", "region_b", "same_region", "size_a_bp", "size_b_bp",
                    "overlap_bp", "frac_a_covered", "frac_b_covered", "overlap_score","jaccard")]

  list(calls = dplyr::select(calls, -gr), pairs = pairs)
}

seg = readr::read_tsv(opt$segments, show_col_types = FALSE)
if (!"role" %in% names(seg)) seg$role = rep("call", nrow(seg))

calls_file = file.path(opt$outdir, paste0(opt$sample, "_consensus_calls.tsv"))
segs_file = file.path(opt$outdir, paste0(opt$sample, "_consensus_segments.tsv"))
pairs_file = file.path(opt$outdir, paste0(opt$sample, "_call_overlaps.tsv"))
regions_file = file.path(opt$outdir, paste0(opt$sample, "_regions.tsv"))

if (nrow(seg) == 0) {
  message("No segments for ", opt$sample, ", writing empty outputs")
  write.table(
    data.frame(
      ecdna_region_id = integer(), sample_id = character(), caller = character(),
      role = character(), call_id = character(), chr_set = character(), is_interchr = logical(),
      nested_in_region = integer()
    ),
    calls_file, sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.table(
    data.frame(
      sample_id = character(), chr = character(), start = integer(), end = integer(),
      strand = character(), segment_order = integer(), caller = character(),
      call_id = character(), role = character(), ecdna_region_id = integer()
    ),
    segs_file, sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.table(
    data.frame(
      sample_id = character(), call_id_a = character(), caller_a = character(),
      call_id_b = character(), caller_b = character(), region_a = integer(),
      region_b = integer(), same_region = logical(), size_a_bp = numeric(),
      size_b_bp = numeric(), overlap_bp = numeric(), frac_a_covered = numeric(),
      frac_b_covered = numeric(), overlap_score = numeric(), jaccard = numeric()
    ),
    pairs_file, sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.table(
    data.frame(
      sample_id = character(), ecdna_region_id = integer(), n_calls = integer(),
      n_callers = integer(), callers = character(), chr_set = character(),
      nested_regions = character()
    ),
    regions_file, sep = "\t", quote = FALSE, row.names = FALSE
  )
  quit(status = 0)
}

core = seg[seg$role != "cycle", ]
clustered = cluster_calls_complete(core, opt$min_overlap, opt$nested_frac)
consensus_calls = clustered$calls

cycle_calls = unique(seg[seg$role == "cycle", c("sample_id", "caller", "call_id", "role")])
if (nrow(cycle_calls) > 0) {
  region_ids = sort(unique(consensus_calls$ecdna_region_id))
  region_grs = list()
  for (r in region_ids) {
    s = core[core$call_id %in% consensus_calls$call_id[consensus_calls$ecdna_region_id == r], ]
    region_grs[[length(region_grs) + 1]] = GenomicRanges::reduce(
      GenomicRanges::GRanges(s$chr, IRanges::IRanges(s$start, s$end)), ignore.strand = TRUE
    )
  }
  cycle_calls$ecdna_region_id = rep(NA_integer_, nrow(cycle_calls))
  cycle_calls$chr_set = rep("", nrow(cycle_calls))
  for (i in seq_len(nrow(cycle_calls))) {
    s = seg[seg$call_id == cycle_calls$call_id[i], ]
    gr_cyc = GenomicRanges::reduce(GenomicRanges::GRanges(s$chr, IRanges::IRanges(s$start, s$end)), ignore.strand = TRUE)
    shared = sapply(region_grs, function(g) shared_bp(gr_cyc, g))
    if (max(shared) > 0) cycle_calls$ecdna_region_id[i] = region_ids[which.max(shared)]
    cycle_calls$chr_set[i] = paste(sort(unique(s$chr)), collapse = ",")
  }
  cycle_calls$is_interchr = grepl(",", cycle_calls$chr_set)
  cycle_calls$nested_in_region = rep(NA_integer_, nrow(cycle_calls))
  message("Cycles assigned to a region: ", sum(!is.na(cycle_calls$ecdna_region_id)), " of ", nrow(cycle_calls))
  consensus_calls = dplyr::bind_rows(consensus_calls, cycle_calls)
}

core_calls = consensus_calls[consensus_calls$role != "cycle", ]

message("Total calls: ", nrow(core_calls))
message("Total regions: ", dplyr::n_distinct(core_calls$ecdna_region_id))
print(
  core_calls %>%
    dplyr::group_by(ecdna_region_id) %>%
    dplyr::summarise(
      n_calls = dplyr::n(),
      n_callers = dplyr::n_distinct(caller),
      callers = paste(sort(unique(caller)), collapse = ","),
      .groups = "drop"
    ),
  n = Inf
)

consensus_calls_out = consensus_calls %>%
  dplyr::select(ecdna_region_id, sample_id, caller, role, call_id, chr_set, is_interchr, nested_in_region) %>%
  dplyr::arrange(ecdna_region_id, caller, call_id)
write_tsv(consensus_calls_out, calls_file)

consensus_segments = seg %>%
  dplyr::left_join(
    dplyr::select(consensus_calls, call_id, ecdna_region_id),
    by = "call_id"
  ) %>%
  dplyr::select(sample_id, chr, start, end, strand, segment_order, caller, call_id, role, ecdna_region_id) %>%
  dplyr::arrange(ecdna_region_id, call_id, segment_order)
write_tsv(consensus_segments, segs_file)

pairs_out = dplyr::arrange(clustered$pairs, dplyr::desc(same_region), dplyr::desc(overlap_score))
write_tsv(pairs_out, pairs_file)

nested_map = consensus_calls %>%
  dplyr::filter(!is.na(nested_in_region)) %>%
  dplyr::group_by(parent = nested_in_region) %>%
  dplyr::summarise(nested_regions = paste(sort(unique(ecdna_region_id)), collapse = ","), .groups = "drop")

regions_out = core_calls %>%
  dplyr::group_by(sample_id, ecdna_region_id) %>%
  dplyr::summarise(
    n_calls = dplyr::n(),
    n_callers = dplyr::n_distinct(caller),
    callers = paste(sort(unique(caller)), collapse = ","),
    chr_set = paste(sort(unique(unlist(strsplit(chr_set, ",")))), collapse = ","),
    .groups = "drop"
  ) %>%
  dplyr::left_join(nested_map, by = c("ecdna_region_id" = "parent")) %>%
  dplyr::arrange(ecdna_region_id)
write_tsv(regions_out, regions_file, na = "")

message("Written: ", calls_file)
message("Written: ", segs_file)
message("Written: ", pairs_file)
message("Written: ", regions_file)
