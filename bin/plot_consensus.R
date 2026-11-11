#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(GenomicRanges)
  library(gTrack)
  library(gGnome)
  library(dplyr)
  library(readr)
  library(scales)
  library(rtracklayer)
})

opt = parse_args(OptionParser(option_list = list(
  make_option("--sample",             type = "character"),
  make_option("--consensus_segments", type = "character"),
  make_option("--consensus_calls",    type = "character"),
  make_option("--decoil_bed",         type = "character", default = "NO_FILE"),
  make_option("--decoil_links",       type = "character", default = "NO_FILE"),
  make_option("--gtf",                type = "character"),
  make_option("--oncogenes",          type = "character", default = ""),
  make_option("--max_genes",          type = "integer", default = 10L),
  make_option("--flye_pafs",          type = "character", default = ""),
  make_option("--flye_min_overlap",   type = "double", default = 0.5),
  make_option("--svclone_junctions",  type = "character", default = ""),
  make_option("--outdir",             type = "character", default = ".")
)))

stopifnot(!is.null(opt$sample), !is.null(opt$consensus_segments), !is.null(opt$consensus_calls), !is.null(opt$gtf))


tools = c("AmpliconArchitect", "Decoil", "CoRAL", "Flye")

build_graph_track = function(segments, track_name, node_col = "gray40",
                             edge_col_start = "orange", edge_col_end = "red",
                             links = NULL, edge_labels = NULL) {

  nodes_gr = makeGRangesFromDataFrame(
    segments,
    seqnames.field     = "chr",
    start.field        = "start",
    end.field          = "end",
    strand.field       = "strand",
    keep.extra.columns = TRUE
  )
  nodes_gr$label        = ""
  nodes_gr$col          = alpha(node_col, 0.7)
  nodes_gr$border.col   = alpha(node_col, 0.7)
  nodes_gr$terminal.col = alpha(node_col, 0.7)
  nodes_gr$left.col     = alpha(node_col, 0.7)
  nodes_gr$right.col    = alpha(node_col, 0.7)

  n_segs = nrow(segments)

  if (!is.null(links) && nrow(links) > 0) {
    node_ids = sort(unique(c(links$n1, links$n2)))
    if (length(node_ids) != n_segs) {
      stop("Decoil circle ", track_name, ": links have ", length(node_ids),
           " nodes but BED has ", n_segs, " segments")
    }
    edge_colors = colorRampPalette(c(edge_col_start, edge_col_end))(nrow(links))
    edges_df = data.frame(
      circle_id = track_name,
      n1      = match(links$n1, node_ids),
      n2      = match(links$n2, node_ids),
      n1.side = as.integer(links$n1.side),
      n2.side = as.integer(links$n2.side),
      col     = alpha(edge_colors, 1),
      lwd     = 2
    )
    if (anyNA(edges_df$n1) || anyNA(edges_df$n2)) {
      stop("Decoil link node could not be mapped to a segment (", track_name, ")")
    }

  } else {
    if (n_segs == 1) {
      edges_df = data.frame(
        circle_id = track_name, n1 = 1L, n2 = 1L,
        n1.side = 0L, n2.side = 1L,
        col = alpha(edge_col_start, 1), lwd = 2
      )
    } else {
      edge_colors = colorRampPalette(c(edge_col_start, edge_col_end))(n_segs)
      edges_df = data.frame(
        circle_id = track_name,
        n1      = seq_len(n_segs),
        n2      = c(seq(2, n_segs), 1L),
        n1.side = ifelse(segments$strand == "-", 0L, 1L),
        n2.side = ifelse(c(segments$strand[-1], segments$strand[1]) == "-", 1L, 0L),
        col     = alpha(edge_colors, 1),
        lwd     = 2
      )
    }

    non_self   = edges_df[edges_df$n1 != edges_df$n2, ]
    self_loops = unique(edges_df[edges_df$n1 == edges_df$n2, ])
    self_loops$n1.side = ifelse(self_loops$n1.side == 1L & self_loops$n2.side == 1L, 0L, self_loops$n1.side)
    self_loops$n2.side = ifelse(self_loops$n1.side == 0L, 1L, self_loops$n2.side)
    nodes_with_any_edge = unique(c(non_self$n1, non_self$n2))
    self_loops = self_loops[!self_loops$n1 %in% nodes_with_any_edge, ]
    edges_df = rbind(non_self, self_loops)
  }

  if (!is.null(edge_labels) && length(edge_labels) == n_segs) {
    nxt_idx = c(seq_len(n_segs)[-1], 1L)
    fwd = edges_df$n2 == nxt_idx[edges_df$n1]
    rev = edges_df$n1 == nxt_idx[edges_df$n2]
    edges_df$label = ifelse(fwd, edge_labels[edges_df$n1], ifelse(rev, edge_labels[edges_df$n2], ""))
    edges_df$label = sub("^0\\.", ".", sub("^1\\.00$", "1", edges_df$label))
    edges_df$lty = ifelse(grepl("\\*$", edges_df$label), 2L, 1L)
  }

  edges_df = as.data.frame(edges_df)

  graph = gG(nodes = nodes_gr, edges = edges_df,
             meta  = list(sep.bg.col    = alpha("white", 0),
                          bundle        = FALSE,
                          gr.labelfield = "label",
                          gr.cex.label  = 1,
                          draw.nodes    = FALSE))
  gt            = graph$gt
  gt$name       = track_name
  gt$sep.bg.col = NA
  gt$xaxis.cex  = 0.2
  formatting(gt)$node.col       = alpha(node_col, 0.7)
  formatting(gt)$node.border    = alpha(node_col, 0.7)
  formatting(gt)$nodeborder.col = alpha(node_col, 0.7)
  return(gt)
}

crop_pdf = function(f, m = 20) {
  gs = Sys.which("gs")
  if (gs == "") return(invisible(NULL))
  bb = suppressWarnings(system2(gs, c("-q", "-dBATCH", "-dNOPAUSE", "-sDEVICE=bbox", shQuote(f)), stdout = TRUE, stderr = TRUE))
  b = grep("^%%BoundingBox:", bb, value = TRUE)
  if (length(b) == 0) return(invisible(NULL))
  v = as.numeric(strsplit(trimws(sub("^%%BoundingBox:", "", b[1])), " +")[[1]])
  tmp = paste0(f, ".crop.pdf")
  ps = sprintf("<</PageOffset [%f %f]>> setpagedevice", m - v[1], m - v[2])
  st = system2(gs, c("-q", "-dBATCH", "-dNOPAUSE", "-sDEVICE=pdfwrite",
                     paste0("-dDEVICEWIDTHPOINTS=", v[3] - v[1] + 2 * m),
                     paste0("-dDEVICEHEIGHTPOINTS=", v[4] - v[2] + 2 * m),
                     "-dFIXEDMEDIA", "-dAutoRotatePages=/None", paste0("-sOutputFile=", shQuote(tmp)),
                     "-c", shQuote(ps), "-f", shQuote(f)))
  if (st == 0 && file.exists(tmp)) file.rename(tmp, f)
  invisible(NULL)
}

draw_ecdna_cycle = function(calls_df, genes_gr, out_dir = ".", decoil_links = NULL, decoil_prop = NULL, only_regions = NULL, suffix = "ecdna",
                            oncogenes = character(0), max_genes = 10, svclone_junctions = NULL) {
  caller_styles = list(
    AmpliconArchitect = list(node_col = "#1a6600", edge_col_start = "#c2ffab", edge_col_end = "#39FF14"),
    Decoil            = list(node_col = "#7a007a", edge_col_start = "#ffabff", edge_col_end = "#FF00FF"),
    CoRAL             = list(node_col = "#00607a", edge_col_start = "#abefff", edge_col_end = "#00EFFF"),
    Flye              = list(node_col = "#7a4400", edge_col_start = "#ffdc99", edge_col_end = "#FF8C00")
  )

  apply_fmt = function(gt, n_segs = 2) {
    formatting(gt)$ywid     = 0.3
    formatting(gt)$height   = max(20, 3 * n_segs)
    formatting(gt)$ygap     = 20
    formatting(gt)$sep.draw = FALSE
    gt
  }

  decoil_ids = sort(unique(as.character(calls_df$call_id[calls_df$caller == "Decoil"])))

  decoil_colors = setNames(
    grDevices::hcl.colors(max(length(decoil_ids), 1), palette = "Dark 3")[seq_along(decoil_ids)],
    decoil_ids
  )

  region_ids = sort(unique(calls_df$ecdna_region_id))
  if (!is.null(only_regions)) region_ids = intersect(region_ids, only_regions)
  message("Total regions to plot: ", length(region_ids))

  for (region_id in region_ids) {
    message("Starting region ", region_id)
    withCallingHandlers(
      tryCatch({
        region_rows = dplyr::filter(calls_df, ecdna_region_id == region_id)
        region_rows = region_rows[order(match(region_rows$caller, tools), !(region_rows$role %in% "footprint")), ]
        sample_id   = region_rows$sample_id[1]
        message("Plotting ecdna_region_id ", region_id, " (sample: ", sample_id,
                ", ", nrow(region_rows), " calls)")
        amp_tracks    = list()
        all_segs_list = list()

        for (k in seq_len(nrow(region_rows))) {
          row     = region_rows[k, ]
          segs_gr = row$gr[[1]]
          strands = row$strands[[1]]
          message("caller=", row$caller, " n_segs=", length(segs_gr))

          if (length(segs_gr) == 0) {
            message("  Row ", k, ": no segments, skipping")
            next
          }
          raw_strands = as.character(strands)
          raw_strands[is.na(raw_strands) | !raw_strands %in% c("+", "-")] = "*"
          cycle_segs = data.frame(
            chr        = as.character(seqnames(segs_gr)),
            start      = start(segs_gr),
            end        = end(segs_gr),
            strand     = raw_strands,
            segment_id = if (!is.null(segs_gr$segment_id)) segs_gr$segment_id else seq_len(length(segs_gr)),
            stringsAsFactors = FALSE
          )
          message("  Row ", k, " caller=", row$caller, " segs=", nrow(cycle_segs))
          all_segs_list[[k]] = segs_gr

          style = caller_styles[[row$caller]]
          short_caller = ifelse(row$caller == "AmpliconArchitect", "AA", row$caller)
          track_name = short_caller
          if (row$caller == "Flye") track_name = paste0("Flye ", row$call_id)
          if (row$caller %in% c("AmpliconArchitect", "CoRAL")) {
            id_parts = strsplit(as.character(row$call_id), "\\|")[[1]]
            if (length(id_parts) >= 4) track_name = paste0(short_caller, " ", sub("amplicon", "amp", id_parts[3]), " cyc", id_parts[4])
          }

          if (!is.na(row$role) && row$role == "footprint") {
            fp_gr = segs_gr
            fp_gr$col = alpha(style$node_col, 0.5)
            fp_kind = ifelse(grepl("_FAN_", as.character(row$call_id)), "NC", "footprint")
            fp_feature = sub("^amplicon([0-9]+)_", "amp\\1 ", strsplit(as.character(row$call_id), "\\|")[[1]][3])
            fp_track = gTrack(fp_gr, name = paste0(short_caller, " ", fp_kind, "\n", fp_feature))
            formatting(fp_track)$height = 24
            formatting(fp_track)$ygap   = 20
            amp_tracks[[k]] = fp_track
            next
          }

          links_k = NULL

          if (row$caller == "Decoil") {
            circle_id  = sub("^.*\\|", "", as.character(row$call_id))
            circle_col = unname(decoil_colors[as.character(row$call_id)])
            style = list(node_col = circle_col, edge_col_start = circle_col, edge_col_end = circle_col)
            track_name = paste0("Decoil ", circle_id)
            if (!is.null(decoil_prop) && !is.na(decoil_prop[circle_id])) track_name = paste0(track_name, "\nprop = ", decoil_prop[circle_id])

            if (!is.null(decoil_links)) {
              links_k = dplyr::filter(decoil_links, circ_id == circle_id)
              cycle_segs$strand = "*"
            }
          }

          if (row$ecdna_region_id != region_id) {
            style = list(node_col = "#4d4d4d", edge_col_start = "#bdbdbd", edge_col_end = "#1a1a1a")
          }

          edge_labels = NULL
          if (!is.null(svclone_junctions)) {
            jr = svclone_junctions[svclone_junctions$call_id == as.character(row$call_id), ]
            if (nrow(jr) == nrow(cycle_segs)) {
              edge_labels = ifelse(jr$matched, paste0(sprintf("%.2f", jr$ccf), ifelse(jr$clonal %in% FALSE, "*", "")), "")
            }
          }

          gt = build_graph_track(
            cycle_segs,
            track_name     = track_name,
            node_col       = style$node_col,
            edge_col_start = style$edge_col_start,
            edge_col_end   = style$edge_col_end,
            links          = links_k,
            edge_labels    = edge_labels
          )
          amp_tracks[[k]] = apply_fmt(gt, n_segs = nrow(cycle_segs))
        }

        amp_tracks = Filter(Negate(is.null), amp_tracks)
        if (length(amp_tracks) == 0) stop("no tracks for region")

        all_segs_gr = do.call(c, Filter(Negate(is.null), all_segs_list))
        windows     = GenomicRanges::reduce(all_segs_gr + 200000)

        genes_in = subsetByOverlaps(genes_gr, windows, ignore.strand = TRUE)
        if (length(genes_in) > max_genes) {
          message(length(genes_in), " genes in window, showing oncogenes only")
          genes_in = genes_in[toupper(genes_in$gene_name) %in% oncogenes]
        }

        n_amp_tracks    = length(amp_tracks)
        track_list      = amp_tracks
        if (length(genes_in) > 0) {
          local_gene_track = make_gene_track(genes_in, oncogenes)
          formatting(local_gene_track)$ygap     = 2
          formatting(local_gene_track)$sep.draw = FALSE
          formatting(local_gene_track)$height   = 18
          track_list = c(list(local_gene_track), amp_tracks)
        }
        combined_tracks = Reduce(c, track_list)

        safe_region = gsub("[|/]", "_", region_id)
        out_file   = file.path(out_dir, paste0(sample_id, "_region", safe_region, "_", suffix, ".pdf"))
        pdf_height = (1200 * n_amp_tracks + 1600) / 300

        pdf(out_file, width = 20, height = pdf_height)
        par(mar = c(8, 2, 0, 10))
        gTrack::plot.gTrack(combined_tracks, windows,
                            sep.draw  = FALSE, cex.ylabel = 1.1)
        dev.off()
        tryCatch(crop_pdf(out_file), error = function(e) message("crop skipped: ", conditionMessage(e)))
        message("Saved: ", out_file)

      }, error = function(e) {
        if (dev.cur() > 1) dev.off()
        message("ERROR in region ", region_id, ": ", conditionMessage(e))
      }),
      warning = function(w) {
        message("WARNING in region ", region_id, ": ", conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
  }
}

load_genes = function(gtf_file) {
  message("Loading gene annotation from: ", gtf_file)
  gtf      = rtracklayer::import(gtf_file, feature.type = "gene")
  genes_gr = gtf[gtf$type == "gene" & gtf$gene_type == "protein_coding"]
  if (length(genes_gr) == 0) genes_gr = gtf[gtf$type == "gene"]
  genes_gr
}

make_gene_track = function(genes_gr, oncogenes = character(0)) {
  is_onc = toupper(genes_gr$gene_name) %in% oncogenes
  genes_gr$col = alpha("#7F7F7F", 0.35)
  genes_gr$col[is_onc] = alpha("#0B6623", 0.3)
  gt = gTrack(genes_gr,
              name          = "Genes",
              gr.labelfield = "gene_name",
              gr.cex.label  = 1.5)
  gt$height = 40
  gt
}

seg = readr::read_tsv(opt$consensus_segments, show_col_types = FALSE)
if (!"role" %in% names(seg)) seg$role = rep("call", nrow(seg))

if (nrow(seg) == 0) {
  message("No consensus segments for ", opt$sample, ", nothing to plot")
  quit(status = 0)
}

calls_df = seg %>%
  dplyr::arrange(ecdna_region_id, caller, call_id, segment_order) %>%
  dplyr::group_by(sample_id, ecdna_region_id, caller, call_id) %>%
  dplyr::summarise(
    chr_set = paste(sort(unique(chr)), collapse = ","),
    role    = dplyr::first(role),
    strands = list(strand),
    gr      = list(GenomicRanges::GRanges(chr, IRanges::IRanges(start, end))),
    .groups = "drop"
  ) %>%
  dplyr::filter(caller %in% tools) %>%
  dplyr::left_join(
    dplyr::select(readr::read_tsv(opt$consensus_calls, show_col_types = FALSE), call_id, nested_in_region),
    by = "call_id"
  )

decoil_links_df = NULL
if (basename(opt$decoil_links) != "NO_FILE" && file.exists(opt$decoil_links)) {
  decoil_links_df = read.table(opt$decoil_links, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
    dplyr::mutate(circ_id = as.character(circ_id))
}

decoil_prop = NULL
if (basename(opt$decoil_bed) != "NO_FILE" && file.exists(opt$decoil_bed)) {
  dbed = read.table(opt$decoil_bed, header = TRUE, sep = "\t", comment.char = "", stringsAsFactors = FALSE)
  if ("estimated_proportions" %in% names(dbed)) {
    decoil_prop = tapply(dbed$estimated_proportions, as.character(dbed$circ_id), function(x) x[1])
  }
}

flye_regions = NULL
if (nzchar(opt$flye_pafs)) {
  flye_rows = list()
  for (f in strsplit(opt$flye_pafs, ",")[[1]]) {
    if (file.info(f)$size == 0) next
    rid = as.integer(sub("^.*_region([0-9]+)_circular_to_ref\\.paf$", "\\1", basename(f)))
    fields = strsplit(readLines(f), "\t", fixed = TRUE)
    paf = data.frame(
      contig = vapply(fields, `[`, "", 1),
      query_start = as.integer(vapply(fields, `[`, "", 3)),
      strand = vapply(fields, `[`, "", 5),
      chr = vapply(fields, `[`, "", 6),
      start = as.integer(vapply(fields, `[`, "", 8)) + 1L,
      end = as.integer(vapply(fields, `[`, "", 9)),
      stringsAsFactors = FALSE
    )
    paf = paf[paf$end >= paf$start, ]
    foot_df = dplyr::filter(seg, ecdna_region_id == rid)
    foot = GenomicRanges::reduce(GenomicRanges::GRanges(foot_df$chr, IRanges::IRanges(foot_df$start, foot_df$end)))
    for (ctg in unique(paf$contig)) {
      p = paf[paf$contig == ctg, ]
      p = p[order(p$query_start), ]
      ctg_gr = GenomicRanges::reduce(GenomicRanges::GRanges(p$chr, IRanges::IRanges(p$start, p$end)))
      ctg_bp = sum(IRanges::width(ctg_gr))
      ov_bp = sum(IRanges::width(suppressWarnings(GenomicRanges::intersect(ctg_gr, foot, ignore.strand = TRUE))))
      ov_frac = ov_bp / ctg_bp
      message("Region ", rid, " Flye contig ", ctg, ": ", round(100 * ov_frac, 1), "% of aligned bases overlap the region")
      if (ov_frac < opt$flye_min_overlap) {
        message("  dropped (below --flye_min_overlap ", opt$flye_min_overlap, ")")
        next
      }
      flye_rows[[length(flye_rows) + 1]] = dplyr::tibble(
        sample_id = opt$sample, ecdna_region_id = rid, caller = "Flye", call_id = ctg,
        chr_set = paste(sort(unique(p$chr)), collapse = ","),
        strands = list(p$strand),
        gr = list(GenomicRanges::GRanges(p$chr, IRanges::IRanges(p$start, p$end))),
        role = "call",
        nested_in_region = NA_integer_
      )
    }
  }
  if (length(flye_rows) == 0) {
    message("No Flye alignments for ", opt$sample, ", nothing to plot")
    quit(status = 0)
  }
  flye_df = dplyr::bind_rows(flye_rows)
  flye_regions = unique(flye_df$ecdna_region_id)
  calls_df = dplyr::bind_rows(calls_df, flye_df)
}

oncogenes = NULL
if (nzchar(opt$oncogenes) && file.exists(opt$oncogenes) && file.info(opt$oncogenes)$size > 0) {
  oncogenes = toupper(trimws(readLines(opt$oncogenes)))
  oncogenes = oncogenes[nzchar(oncogenes)]
  message("Oncogene list: ", length(oncogenes), " genes")
}

genes_gr = load_genes(opt$gtf)
if (is.null(oncogenes)) oncogenes = character(0)

svclone_df = NULL
if (nzchar(opt$svclone_junctions) && file.exists(opt$svclone_junctions) && file.info(opt$svclone_junctions)$size > 0) {
  svclone_df = readr::read_tsv(opt$svclone_junctions, show_col_types = FALSE)
  message("SVclone junctions: ", nrow(svclone_df), " rows, ", sum(svclone_df$matched), " matched")
}

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

plot_suffix = "ecdna"
if (!is.null(flye_regions)) plot_suffix = "flye"
if (is.null(flye_regions) && !is.null(svclone_df)) plot_suffix = "svclone"

draw_ecdna_cycle(calls_df, genes_gr, opt$outdir, decoil_links = decoil_links_df, decoil_prop = decoil_prop,
                 only_regions = flye_regions, suffix = plot_suffix,
                 oncogenes = oncogenes, max_genes = opt$max_genes, svclone_junctions = svclone_df)
