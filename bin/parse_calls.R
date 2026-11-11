#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(GenomicRanges)
  library(stringr)
})

opt = parse_args(OptionParser(option_list = list(
  make_option("--sample",                type = "character"),
  make_option("--aa_profiles",           type = "character", default = "NO_FILE"),
  make_option("--aa_beds",               type = "character", default = "NO_FILE"),
  make_option("--aa_cycles",             type = "character", default = "NO_FILE"),
  make_option("--coral_profiles",        type = "character", default = "NO_FILE"),
  make_option("--coral_beds",            type = "character", default = "NO_FILE"),
  make_option("--coral_cycles",          type = "character", default = "NO_FILE"),
  make_option("--decoil_bed",            type = "character", default = "NO_FILE"),
  make_option("--outdir",                type = "character", default = ".")
)))

stopifnot(!is.null(opt$sample))

ecdna_like_cycles = function(cycles_file) {
  lines = readLines(cycles_file)
  seg_lines = grep("^Segment\\s+", lines, value = TRUE)
  segments = dplyr::bind_rows(lapply(seg_lines, function(l) {
    p = strsplit(trimws(l), "\\s+")[[1]]
    data.frame(segment_id = p[2], chr = p[3], start = as.integer(p[4]), end = as.integer(p[5]))
  }))
  cyc_lines = grep("^Cycle=", lines, value = TRUE)
  cyc_lines = grep("CycleClass=ecDNA-like(;|$)", cyc_lines, value = TRUE)
  cyc_lines = grep("IsCyclicPath=True(;|$)", cyc_lines, value = TRUE)
  if (length(cyc_lines) == 0) return(NULL)
  cycles = dplyr::bind_rows(lapply(cyc_lines, function(l) {
    path = strsplit(sub(".*Segments=([^;]+).*", "\\1", l), ",")[[1]]
    if (any(grepl("^0[+-]?$", path))) return(NULL)
    data.frame(
      cycle_id = sub(".*Cycle=([^;]+).*", "\\1", l),
      segment_order = seq_along(path),
      segment_id = sub("[+-]$", "", path),
      strand = sub("^.*([+-])$", "\\1", path)
    )
  }))
  if (nrow(cycles) == 0) return(NULL)
  cycles = dplyr::left_join(cycles, segments, by = "segment_id")
  if (anyNA(cycles$chr) || anyNA(cycles$start) || anyNA(cycles$end)) {
    stop("Some cycle segments have missing genomic coordinates in ", cycles_file)
  }
  cycles
}

read_ecdna_amplicons = function(prof_file, classes = "ecDNA") {
  prof = read.delim(prof_file, stringsAsFactors = FALSE)
  cols = paste0(classes, ".")
  cols = cols[cols %in% names(prof)]
  if (!"amplicon_number" %in% names(prof) || length(cols) == 0) {
    stop("Unexpected columns in ", prof_file, ": ", paste(names(prof), collapse = ", "))
  }
  is_pos = rep(FALSE, nrow(prof))
  for (cl in cols) is_pos = is_pos | prof[[cl]] %in% "Positive"
  pos = prof$amplicon_number[is_pos]
  if (length(pos) == 0) return(character(0))
  paste0("amplicon", as.integer(gsub("\\D", "", pos)))
}

load_classified_calls = function(sample, caller, tag, prof_file, bed_dir, cycles_dir, classes = "ecDNA") {
  ecdna_amps = read_ecdna_amplicons(prof_file, classes)
  class_re = paste0("(", paste(classes, collapse = "|"), ")")
  message(caller, " ecDNA-positive amplicons: ", paste(ecdna_amps, collapse = ", "))
  if (length(ecdna_amps) == 0) return(NULL)
  beds = list.files(
    bed_dir,
    pattern = paste0("_amplicon[0-9]+_", class_re, "_[0-9]+_intervals\\.bed$"),
    full.names = TRUE
  )
  beds = beds[file.info(beds)$size > 0]
  beds = beds[stringr::str_extract(basename(beds), "amplicon[0-9]+") %in% ecdna_amps]
  if (length(beds) == 0) {
    message("No ecDNA interval BEDs found in ", bed_dir, " for ", sample)
    return(NULL)
  }

  footprints = dplyr::bind_rows(lapply(beds, function(f) {
    id = stringr::str_extract(basename(f), paste0("amplicon[0-9]+_", class_re, "_[0-9]+"))
    bed = read.table(f, sep = "\t", header = FALSE)[, 1:3]
    colnames(bed) = c("chr", "start", "end")
    data.frame(
      sample_id = sample, chr = bed$chr, start = bed$start, end = bed$end,
      strand = "*", segment_order = seq_len(nrow(bed)),
      caller = caller, call_id = paste(sample, tag, id, "footprint", sep = "|"), role = "footprint"
    )
  }))

  cycle_files = list.files(cycles_dir, pattern = "_annotated_cycles\\.txt$", full.names = TRUE)
  cycle_amps = stringr::str_extract(basename(cycle_files), "amplicon[0-9]+")

  cycles = dplyr::bind_rows(lapply(beds, function(f) {
    a = stringr::str_extract(basename(f), "amplicon[0-9]+")
    id = stringr::str_extract(basename(f), paste0("amplicon[0-9]+_", class_re, "_[0-9]+"))
    cf = cycle_files[cycle_amps %in% a][1]
    if (is.na(cf)) {
      message(id, ": cycles file not found, no cycle added")
      return(NULL)
    }
    cyc = ecdna_like_cycles(cf)
    if (is.null(cyc)) return(NULL)
    bed = read.table(f, sep = "\t", header = FALSE)[, 1:3]
    bed_gr = GenomicRanges::reduce(GenomicRanges::GRanges(bed[[1]], IRanges::IRanges(bed[[2]] + 1L, bed[[3]])), ignore.strand = TRUE)
    cycle_ids = unique(cyc$cycle_id)
    jac = sapply(cycle_ids, function(cid) {
      s = cyc[cyc$cycle_id == cid, ]
      g = GenomicRanges::reduce(GenomicRanges::GRanges(s$chr, IRanges::IRanges(s$start, s$end)), ignore.strand = TRUE)
      inter = sum(IRanges::width(suppressWarnings(GenomicRanges::intersect(g, bed_gr, ignore.strand = TRUE))))
      inter / (sum(IRanges::width(g)) + sum(IRanges::width(bed_gr)) - inter)
    })
    best = cycle_ids[which.max(jac)]
    message(id, ": ", length(cycle_ids), " cycles, best Jaccard with BED = ", round(max(jac), 3), " (cycle ", best, ")")
    cyc[cyc$cycle_id == best, ] %>%
      dplyr::transmute(
        sample_id = sample, chr, start, end, strand, segment_order,
        caller = caller,
        call_id = paste(sample, tag, a, cycle_id, sep = "|"),
        role = "cycle"
      )
  }))
  if (nrow(cycles) > 0) cycles = cycles[!duplicated(cycles[, c("call_id", "segment_order")]), ]

  dplyr::bind_rows(footprints, cycles)
}

load_aa_calls = function(sample, prof_file, bed_dir, cycles_dir) {
  load_classified_calls(sample, "AmpliconArchitect", "AA", prof_file, bed_dir, cycles_dir, classes = c("ecDNA", "FAN"))
}

load_coral_calls = function(sample, prof_file, bed_dir, cycles_dir) {
  load_classified_calls(sample, "CoRAL", "CoRAL", prof_file, bed_dir, cycles_dir, classes = c("ecDNA", "FAN"))
}

load_decoil_calls = function(sample, f) {
  if (!file.exists(f)) return(NULL)
  bed = read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "")
  colnames(bed)[1] = "chr"
  if (nrow(bed) == 0) return(NULL)
  bed %>%
    dplyr::group_by(circ_id) %>%
    dplyr::mutate(segment_order = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      sample_id = sample, chr, start, end, strand, segment_order,
      caller = "Decoil",
      call_id = paste(sample, "Decoil", circ_id, sep = "|"),
      role = "call"
    )
}

parse_calls = function(sample, aa_profiles, aa_beds, aa_cycles, coral_profiles, coral_beds, coral_cycles, decoil_bed) {
  present_file = function(f) !is.na(f) && basename(f) != "NO_FILE" && file.exists(f) && !dir.exists(f) && file.info(f)$size > 0
  present_dir = function(d) !is.na(d) && dir.exists(d) && any(list.files(d) != "NO_FILE")
  parts = list()
  if (present_file(coral_profiles) && present_dir(coral_beds) && present_dir(coral_cycles)) {
    parts$coral = load_coral_calls(sample, coral_profiles, coral_beds, coral_cycles)
  }
  if (present_file(decoil_bed)) parts$decoil = load_decoil_calls(sample, decoil_bed)
  if (present_file(aa_profiles) && present_dir(aa_beds) && present_dir(aa_cycles)) {
    parts$aa = load_aa_calls(sample, aa_profiles, aa_beds, aa_cycles)
  }
  out = dplyr::bind_rows(parts)
  if (nrow(out) == 0) return(NULL)
  out %>%
    dplyr::mutate(chr = ifelse(startsWith(chr, "chr"), chr, paste0("chr", chr))) %>%
    dplyr::select(sample_id, chr, start, end, strand, segment_order, caller, call_id, role)
}

res = parse_calls(
  opt$sample, opt$aa_profiles, opt$aa_beds, opt$aa_cycles,
  opt$coral_profiles, opt$coral_beds, opt$coral_cycles, opt$decoil_bed
)

if (is.null(res)) {
  message("No calls found for ", opt$sample)
  res = data.frame(
    sample_id = character(), chr = character(), start = integer(), end = integer(),
    strand = character(), segment_order = integer(), caller = character(), call_id = character(),
    role = character()
  )
}

message("Calls per caller and role:")
print(table(res$caller[!duplicated(res$call_id)], res$role[!duplicated(res$call_id)]))

out_file = file.path(opt$outdir, paste0(opt$sample, "_segments.tsv"))
write.table(res, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
message("Written: ", out_file)
