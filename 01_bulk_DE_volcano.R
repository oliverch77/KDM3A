#!/usr/bin/env Rscript
# KDM3A: original edgeR differential analysis and volcano plots.
# Run from the repository root; see README.md for input formats and caveats.
# No packages are installed or updated by this script.
# Usage: Rscript 01_bulk_DE_volcano.R hsMDM|mmBMM|RAW counts.tsv output_directory

options(stringsAsFactors = FALSE)

# Required libraries (use the original analysis package versions).
suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(statmod)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
})

# 1. Fixed analysis settings -------------------------------------------------
# Keep these values unchanged when reproducing the original analyses.
settings <- list(
  cpm = 0.5, raw_sample_fraction = 0.33,
  p = 0.05, fdr = 0.05, absolute_logFC = 0,
  raw_color_logFC = 1,             # Approved change: PLOT COLOR ONLY.
  robust_QL = TRUE
)

mouse_comparisons <- data.frame(
  comparison = c("KO_(without_LPS)", "KO_(with_LPS)"),
  base = c("WT_Ctrl", "WT_LPS"),
  target = c("KO_Ctrl", "KO_LPS"),
  row.names = c("KO_noLPS", "KO_withLPS")
)
# Human comparisons and their order were supplied by the author (2026-09-26).
human_comparisons <- data.frame(
  comparison = c("KDM3Akd_(without_LPS)", "KDM3Akd_(with_LPS)"),
  base = c("scrambled_Ctrl", "scrambled_LPS"),
  target = c("KDM3Akd_Ctrl", "KDM3Akd_LPS"),
  row.names = c("KDM3Akd_noLPS", "KDM3Akd_withLPS")
)

human_labels <- c("IFIT1", "MX1", "CCL5", "IL8", "IL10", "TNF", "IL6",
                  "MARCO", "KDM3A", "CXCL10", "SOCS3", "CXCL1")
mouse_labels <- c("Ifit1", "Ifit2", "Ifit3", "Mx1", "Oas1g", "Oas2",
                  "Ifnb1", "Irf7", "Irf3", "Irf5", "Irf8", "Ccl5", "Ccl8",
                  "Cxcl1", "Cxcl2", "Cxcl9", "Cxcl10", "Marco", "Cd163",
                  "Mrc1", "Il10", "Il1rn", "Lgals3", "Aldoa", "Trem2", "Cd9",
                  "Spp1", "Ctsb", "Tnf", "Il6", "Cxcl14", "Mmp9", "Il1b", "Il12b")

dataset_config <- function(dataset) {
  if (!dataset %in% c("hsMDM", "mmBMM", "RAW"))
    stop("dataset must be hsMDM, mmBMM, or RAW.")
  list(dataset = dataset, species = if (dataset == "hsMDM") "Hs" else "Mm",
       comparisons = if (dataset == "hsMDM") human_comparisons else mouse_comparisons,
       labels = if (dataset == "hsMDM") human_labels else
         if (dataset == "RAW") append(mouse_labels, "Ifih1", after = 6) else mouse_labels,
       top_labels = if (dataset == "hsMDM") 3 else 0,
       # These two legacy values were defined but UNUSED by the plot code.
       # They remain inactive here; activating them would change labels.
       label_fdr_unused = 0.05,
       label_logFC_unused = if (dataset == "hsMDM") 2.5 else 0)
}

# 2. Count import and original sample naming ---------------------------------
# From the author's 00_f_base.R. Keep replacements sequential and literal:
# a later pattern may match text introduced by an earlier replacement.
renameinbulk <- function(vector, pattern, replacement) {
  for (i in seq_along(pattern)) {
    vector <- gsub(pattern[i], replacement[i], vector, fixed = TRUE)
  }
  vector
}

parse_sample_names <- function(original, dataset) {
  canonical_pattern <- switch(dataset,
    hsMDM = "^.{3}_(scrambled|KDM3Akd)_(Ctrl|LPS)$",
    mmBMM = "^(WT|KO)_(Ctrl|LPS)_mm[1-4]$",
    RAW = "^(WT|KO)_(Ctrl|LPS)_rep[1-3]$")
  if (all(grepl(canonical_pattern, original))) return(original)

  if (dataset == "hsMDM") {
    end <- regexpr("_R1_", original) - 5
    if (any(end < 14)) stop("Human headers do not match the original _R1_ naming.")
    n <- substr(original, 14, end)
    n <- renameinbulk(n, c("K", "S", "C", "P"),
                       c("_KDM3Akd_", "_scrambled_", "Ctrl", "LPS"))
  } else if (dataset == "mmBMM") {
    n <- renameinbulk(substr(original, 1, 3), c("ko", "wt", "1", "2", "4", "5"),
                       c("KO_", "WT_", "Ctrl_mm1", "Ctrl_mm2", "LPS_mm1", "LPS_mm2"))
    wt <- substr(n, 1, 1) == "W"
    n[wt] <- renameinbulk(n[wt], c("1", "2"), c("3", "4"))
  } else {
    if (length(original) != 12L)
      stop("Original RAW naming assumes 12 columns in their original order.")
    n <- substr(original, 5, 10)
    lps <- substr(n, 1, 1) == "L"
    n[!lps] <- paste0(substr(n[!lps], 1, 2), "_Ctrl")
    n[lps] <- paste0(substr(n[lps], 5, 6), "_LPS")
    n <- paste0(n, "_rep", rep(c(1, 2, 3), times = 4))
  }
  if (anyNA(n) || !all(grepl(canonical_pattern, n)) || anyDuplicated(n))
    stop("Unrecognized/duplicated sample labels. Use original headers or the ",
         "already-renamed headers specified in README.md; do not guess sample identities.")
  n
}

sample_groups <- function(sample_names, dataset) {
  switch(dataset,
    hsMDM = substring(sample_names, 5),
    mmBMM = substr(sample_names, 1, regexpr("mm", sample_names) - 2),
    RAW = substr(sample_names, 1, regexpr("rep", sample_names) - 2))
}

read_counts <- function(path, cfg) {
  if (!file.exists(path)) stop("Count file not found: ", path)
  x <- read.delim(path, header = TRUE, row.names = 1)
  annotation <- c("chr", "start", "end", "Length", "strand", "Copies",
                  "Annotation.Divergence")
  x <- x[, !colnames(x) %in% annotation, drop = FALSE]
  if (!nrow(x) || !ncol(x) || !all(vapply(x, is.numeric, logical(1))))
    stop("Expect a RefSeq-by-sample numeric count table, optionally with HOMER annotations.")
  if (anyNA(x) || any(!is.finite(as.matrix(x))) || any(as.matrix(x) < 0))
    stop("Counts must be finite, nonmissing, and nonnegative.")
  original_names <- colnames(x)
  colnames(x) <- parse_sample_names(original_names, cfg$dataset)
  if (anyDuplicated(colnames(x))) stop("Duplicated sample names.")
  groups <- sample_groups(colnames(x), cfg$dataset)
  expected <- unique(c(cfg$comparisons$base, cfg$comparisons$target))
  if (!setequal(unique(groups), expected))
    stop("Count-table groups differ from the embedded comparisons: ",
         paste(unique(groups), collapse = ", "))
  mapping <- data.frame(original_name = original_names, sample = colnames(x), group = groups)
  mapping$legacy_zero_filter_block <- if (cfg$dataset == "RAW") NA_character_ else
    substr(colnames(x), 2, 3)
  mapping$model_block <- switch(cfg$dataset,
    hsMDM = substr(colnames(x), 2, 3),
    mmBMM = substring(colnames(x), nchar(colnames(x)) - 2),
    RAW = rep(NA_character_, ncol(x)))
  list(counts = x, mapping = mapping)
}

# 3. Original annotation, filters, and model ---------------------------------
annotate_and_filter <- function(counts, group, cfg) {
  y <- edgeR::DGEList(counts = counts, genes = rownames(counts), group = group)
  pkg <- paste0("org.", cfg$species, ".eg.db")
  prefix <- paste0("org.", cfg$species, ".eg")
  refseq <- getExportedValue(pkg, paste0(prefix, "REFSEQ"))
  refseq_table <- AnnotationDbi::toTable(refseq)
  y$genes$Entrez <- refseq_table$gene_id[match(y$genes$genes, refseq_table$accession)]
  symbol_table <- AnnotationDbi::toTable(getExportedValue(pkg, paste0(prefix, "SYMBOL")))
  y$genes$Sym <- symbol_table$symbol[match(y$genes$Entrez, symbol_table$gene_id)]
  chr_table <- AnnotationDbi::toTable(getExportedValue(pkg, paste0(prefix, "CHR")))
  y$genes$Chr <- chr_table$chromosome[match(y$genes$Entrez, chr_table$gene_id)]
  protein_table <- AnnotationDbi::toTable(getExportedValue(pkg, paste0(prefix, "UNIPROT")))
  y$genes$Uniport <- protein_table$uniprot_id[match(y$genes$Entrez, protein_table$gene_id)]
  # Keep the first mapping and the highest-total-count transcript PER CONTRAST.
  # Do not aggregate transcripts, update symbols, or prefilter across contrasts.
  y <- y[y$genes$genes %in% AnnotationDbi::mappedRkeys(refseq), ]
  y <- y[order(rowSums(y$counts), decreasing = TRUE), ]
  y <- y[!duplicated(y$genes$Entrez), ]
  rownames(y$counts) <- rownames(y$genes) <- y$genes$Entrez

  min_samples <- if (cfg$dataset == "RAW")
    settings$raw_sample_fraction * ncol(edgeR::cpm(y)) else length(unique(group))
  keep <- rowSums(edgeR::cpm(y) > settings$cpm) >= min_samples
  y <- y[keep, , keep.lib.sizes = FALSE]

  if (cfg$dataset != "RAW") {
    # Preserve the literal legacy substring grouping, INCLUDING for mmBMM.
    # For mmBMM this yields O_/T_, not mouse IDs; changing it changes results.
    block <- substr(colnames(y), 2, 3)
    notkeep <- rep(FALSE, nrow(y))
    for (j in unique(block)) {
      notkeep <- as.logical(notkeep + (rowSums(y[, block == j]$counts) == 0))
    }
    y <- y[!notkeep, , keep.lib.sizes = FALSE]
  }
  if (!nrow(y)) stop("No genes remain after the original filters.")
  y
}

make_design <- function(sample_names, group, dataset) {
  if (dataset == "hsMDM") {
    donor <- substr(sample_names, 2, 3)
    design <- model.matrix(~ donor + group)
  } else if (dataset == "mmBMM") {
    mouse_n <- substring(sample_names, nchar(sample_names) - 2)
    design <- model.matrix(~ mouse_n + group)
    design <- design[, !(colSums(design) == 1), drop = FALSE]
  } else {
    design <- model.matrix(~ group)
  }
  rownames(design) <- sample_names
  if (qr(design)$rank != ncol(design) || nrow(design) <= ncol(design))
    stop("Original design is not estimable for these samples; no automatic redesign was applied.")
  if (!startsWith(tail(colnames(design), 1), "group"))
    stop("Last model coefficient is not the intended group comparison.")
  design
}

run_contrast <- function(counts, cfg, i, output_dir) {
  comp <- cfg$comparisons[i, ]
  use <- sample_groups(colnames(counts), cfg$dataset) %in% c(comp$base, comp$target)
  counts <- counts[, use, drop = FALSE]
  group <- relevel(factor(sample_groups(colnames(counts), cfg$dataset)), ref = comp$base)
  y <- annotate_and_filter(counts, group, cfg)
  design <- make_design(colnames(y), group, cfg$dataset)
  y$samples$lib.size <- colSums(y$counts)
  y <- edgeR::calcNormFactors(y)
  y <- edgeR::estimateDisp(y, design = design)
  fit <- edgeR::glmQLFit(y, design, robust = settings$robust_QL)
  coef <- tail(colnames(design), 1)
  qlf <- edgeR::glmQLFTest(fit, coef = coef)
  tab <- edgeR::topTags(qlf, n = nrow(qlf$table))$table

  d <- file.path(output_dir, "differential", i)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  stem <- paste0(comp$target, "_vs_", comp$base)
  result_file <- file.path(d, paste0(stem, "_DEG.csv"))
  write.csv(tab, result_file, row.names = FALSE)
  write.csv(design, file.path(d, "design_matrix.csv"))
  write.csv(y$samples, file.path(d, "normalization_factors.csv"))
  pdf(file.path(d, "BCVplot.pdf"), width = 5, height = 4)
  edgeR::plotBCV(y)
  dev.off()
  pdf(file.path(d, "QLDisp.pdf"), width = 5, height = 4)
  edgeR::plotQLDisp(fit)
  dev.off()
  # The original plotting stage read the exported CSV back in.
  plot_volcano(read.csv(result_file), cfg, comp, file.path(d, paste0(stem, "_volcano.pdf")))
  result_file
}

# 4. Volcano plots ----------------------------------------------------------
volcano_data <- function(tab, cfg) {
  v <- tab[, c("Sym", "logFC", "PValue", "FDR")]
  v$legacy_category <- dplyr::case_when(
    v$FDR < settings$fdr & v$logFC > settings$absolute_logFC ~ "upFDR",
    v$FDR < settings$fdr & v$logFC < -settings$absolute_logFC ~ "dnFDR",
    v$PValue < settings$p & v$FDR > settings$fdr & v$logFC > settings$absolute_logFC ~ "upPV",
    v$PValue < settings$p & v$FDR > settings$fdr & v$logFC < -settings$absolute_logFC ~ "dnPV",
    TRUE ~ "NS")
  v$cols <- v$legacy_category
  if (cfg$dataset == "RAW") {
    # Only the displayed categories change. Tests and label selection do not.
    v$cols <- dplyr::case_when(
      v$FDR < settings$fdr & v$logFC > settings$raw_color_logFC ~ "upStrong",
      v$FDR < settings$fdr & v$logFC < -settings$raw_color_logFC ~ "dnStrong",
      v$FDR < settings$fdr & v$logFC > 0 ~ "upOther",
      v$FDR < settings$fdr & v$logFC < 0 ~ "dnOther",
      TRUE ~ "NS")
  }
  # Preserve the source's label behavior, including RAW's legacy 1:0 indexing.
  # mmBMM alone had a top_labels > 0 guard. Human top_labels remains 3.
  fdr_for_labels <- if (cfg$dataset == "hsMDM") v$FDR else -log10(v$FDR)
  fdr_for_labels[is.infinite(fdr_for_labels)] <- 325
  v$anlab <- "NO"
  v$anlab[match(v$Sym, cfg$labels) & v$legacy_category != "NS"] <- "Yes"
  n <- cfg$top_labels
  if (cfg$dataset != "mmBMM" || n > 0) {
    v$anlab[v$logFC > 0][1:n] <- "Yes"
    v$anlab[v$logFC < 0][1:n] <- "Yes"
    v$anlab[v$logFC >= sort(v$logFC, decreasing = TRUE)[n] & fdr_for_labels < settings$fdr] <- "Yes"
    v$anlab[v$logFC <= sort(v$logFC)[n] & fdr_for_labels < settings$fdr] <- "Yes"
  }
  v$plot_y <- -log10(v$PValue)
  if (cfg$dataset != "hsMDM") v$plot_y[is.infinite(v$plot_y)] <- 325
  v
}

plot_volcano <- function(tab, cfg, comp, path) {
  v <- volcano_data(tab, cfg)
  raw <- cfg$dataset == "RAW"
  mm <- cfg$dataset == "mmBMM"
  palette <- if (raw) c(dnStrong = "#1c8038", dnOther = "#68c74a", NS = "#C0C0C0",
                        upOther = "#f2a85e", upStrong = "#e07204") else
    c(dnFDR = "#1c8038", dnPV = "#68c74a", NS = "#C0C0C0", upPV = "#f2a85e", upFDR = "#e07204")
  legend <- if (raw) c(
    dnStrong = "decreased (FDR < 0.05; log2FC < -1)",
    dnOther = "decreased (FDR < 0.05; -1 <= log2FC < 0)", NS = "not significant",
    upOther = "increased (FDR < 0.05; 0 < log2FC <= 1)",
    upStrong = "increased (FDR < 0.05; log2FC > 1)") else
    c(dnFDR = "decreased (FDR < 0.05)    ", dnPV = "decreased (P < 0.05)    ",
      NS = "not significant    ", upPV = "increased (P < 0.05)    ", upFDR = "increased (FDR < 0.05)    ")
  v$cols <- factor(v$cols, levels = names(palette))
  xlimvalue <- if (mm) max(abs(v$logFC), na.rm = TRUE) * 1.05 else max(abs(v$logFC)) + 0.1
  ylimuplimit <- if (mm) max(v$plot_y, na.rm = TRUE) * 1.05 else
    if (raw) round(max(v$plot_y, na.rm = TRUE) * 1.05) + 0.5 else round(max(v$plot_y)) + 0.5
  # Preserve original FDR / nominal-P annotation counts, including RAW.
  # The authorized RAW change affects point colors and their legend only.
  counts <- table(factor(v$legacy_category,
                         levels = c("dnFDR", "dnPV", "NS", "upPV", "upFDR")))
  up_counts <- c(counts["upFDR"], counts["upFDR"] + counts["upPV"])
  dn_counts <- c(counts["dnFDR"], counts["dnFDR"] + counts["dnPV"])
  stat_labels <- c(paste0("decreased: ", dn_counts[1], " / ", dn_counts[2]),
                   paste0("increased: ", up_counts[1], " / ", up_counts[2]))
  p <- ggplot2::ggplot(v, ggplot2::aes(x = logFC, y = plot_y)) +
    ggplot2::geom_point(alpha = 0.5, pch = 16, size = 2.5, ggplot2::aes(color = cols)) +
    ggplot2::xlim(-xlimvalue, xlimvalue) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(-0.16, ylimuplimit)) +
    ggplot2::theme_minimal(base_size = 24) +
    ggplot2::scale_color_manual(values = palette, labels = legend, name = " ") +
    ggplot2::ggtitle(paste0(comp$target, " vs ", comp$base)) +
    ggplot2::labs(x = expression(paste(log[2], "(fold change)")),
                  y = expression(paste("-", log[10], "(", italic(P), " value)")))
  label_args <- list(data = v[v$anlab == "Yes", ], mapping = ggplot2::aes(label = Sym),
                     fontface = "italic", size = 3, segment.alpha = 0.5, force = 1,
                     point.padding = 0.16, min.segment.length = 0.25)
  if (mm) label_args$max.overlaps <- 12
  p <- p + do.call(ggrepel::geom_text_repel, label_args) +
    ggplot2::theme(legend.position = "top", legend.key = ggplot2::element_blank(),
      legend.box.background = ggplot2::element_blank(), legend.text = ggplot2::element_text(size = 6),
      legend.title = ggplot2::element_text(size = 2, hjust = 0.5),
      legend.background = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 12, hjust = 0.5),
      axis.text = ggplot2::element_text(size = 8), axis.title = ggplot2::element_text(size = 10)) +
    ggplot2::annotate("label", x = c(-1, 1) * xlimvalue * (if (mm) 0.75 else 0.85),
      y = ylimuplimit * (if (mm) 0.05 else 0.075), fill = NA, size = 2.5,
      label = stat_labels, col = scales::alpha(c("#186300", "#b55b00"), 0.8))
  if (mm) p <- p + ggplot2::theme(plot.background = ggplot2::element_blank(),
                                  panel.background = ggplot2::element_blank(), aspect.ratio = 0.75)
  if (cfg$dataset != "hsMDM" && any(is.infinite(-log10(v$PValue)))) {
    p <- p + ggplot2::geom_hline(yintercept = 325, color = "#C0C0C0", linetype = "dashed") +
      ggplot2::annotate("text", y = 333, x = 0,
        label = expression(paste(italic(P), " value below R minimal limit")), colour = "#C0C0C0")
  }
  if (mm) {
    ggplot2::ggsave(plot = p, filename = path, width = 6, height = 5)
  } else {
    pdf(path, width = 6, height = 5)
    print(p)
    dev.off()
  }
}

# 5. Run one dataset --------------------------------------------------------
main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 3L) stop(
    "Usage: Rscript 01_bulk_DE_volcano.R hsMDM|mmBMM|RAW counts.tsv output_directory")
  cfg <- dataset_config(args[1])
  input <- read_counts(args[2], cfg)
  output_dir <- args[3]
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt")), add = TRUE)
  write.csv(input$mapping, file.path(output_dir, "sample_mapping.csv"), row.names = FALSE)
  write.csv(cfg$comparisons, file.path(output_dir, "comparisons_used.csv"))
  dput(list(settings = settings, config = cfg), file.path(output_dir, "settings_used.R"))
  result_files <- setNames(vector("list", nrow(cfg$comparisons)), rownames(cfg$comparisons))
  for (i in rownames(cfg$comparisons)) {
    message("Running ", cfg$dataset, ": ", i)
    result_files[[i]] <- run_contrast(input$counts, cfg, i, output_dir)
  }
  message("Completed: ", normalizePath(output_dir))
  invisible(result_files)
}

if (sys.nframe() == 0L) main()
