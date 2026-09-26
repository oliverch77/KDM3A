#!/usr/bin/env Rscript
# KDM3A: top-100 fold-change signatures in an ALREADY PREPARED Seurat object.
# No normalization, integration, clustering, or embedding is performed here.
# Usage: Rscript 02_single_cell_scores.R DEG_overview.csv prepared.rds output_dir [cluster_column]
# Omit cluster_column to use the object's current Idents, as in the source.

options(stringsAsFactors = FALSE)

# Required libraries; preserve lme4 -> lmerTest loading order.
suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(lme4)
  library(lmerTest)
  library(ggplot2)
})

# 1. Fixed settings ---------------------------------------------------------
settings <- list(
  top_n = 100, signature_p = 0.05, mitochondrial_percent = 7.5,
  ctrl = 100, k = FALSE,
  # These are the Seurat defaults used implicitly by the original calls.
  nbin = 24, seed = 1, search = FALSE, slot = "data",
  patient_column = "Patient", adjustment = "BH", significance = 0.05
)
contrast <- "KDM3Akd_LPS_vs_scrambled_LPS"

# 2. Original top-100 signature selection -----------------------------------
select_signatures <- function(DEG_KDM) {
  fc_column <- paste0(contrast, "_logFC")
  pv_column <- paste0(contrast, "_PV")
  required <- c("Sym", fc_column, pv_column)
  if (!all(required %in% names(DEG_KDM)))
    stop("Use the intersection-based human DEG_overview.csv; required columns: ",
         paste(required, collapse = ", "))
  # Exact source sequence: order rows by FC, restrict to nominal P < 0.05,
  # then take positions 1:100. No new FDR gate, sign gate, deduplication,
  # symbol remapping, or replacement of missing features is introduced.
  dn <- DEG_KDM[order(DEG_KDM[[fc_column]], decreasing = FALSE), ]
  up <- DEG_KDM[order(DEG_KDM[[fc_column]], decreasing = TRUE), ]
  dn_genes <- dn[dn[[pv_column]] < settings$signature_p, "Sym"][1:settings$top_n]
  up_genes <- up[up[[pv_column]] < settings$signature_p, "Sym"][1:settings$top_n]
  if (anyNA(dn_genes) || anyNA(up_genes) || any(!nzchar(c(dn_genes, up_genes))))
    stop("The original top-100 selection contains missing symbols or fewer than 100 valid entries. ",
         "Check the input; no genes were substituted.")
  list(dnFC100 = dn_genes, upFC100 = up_genes)
}

# 3. Original linear mixed model and cell-level Cohen's d --------------------
score_cluster <- function(cluster_name, score, cell_clusters, sample_id) {
  group <- ifelse(cell_clusters == cluster_name, cluster_name, "rest")
  group <- factor(group, levels = c(cluster_name, "rest"))
  in_cluster <- cell_clusters == cluster_name
  out_cluster <- !in_cluster
  mean_in <- mean(score[in_cluster], na.rm = TRUE)
  mean_out <- mean(score[out_cluster], na.rm = TRUE)
  n_in <- sum(in_cluster)
  n_out <- sum(out_cluster)
  df <- data.frame(score = score, group = group, sample_id = factor(sample_id))
  # lmerTest was loaded after lme4 in the source; this is the same fitter.
  # Keep its default REML fit, optimizer, and ANOVA defaults (type III,
  # Satterthwaite denominator degrees of freedom). Do not switch to glmer.
  fit <- lmerTest::lmer(score ~ group + (1 | sample_id), data = df)
  aov_tab <- stats::anova(fit)
  p_val <- aov_tab["group", "Pr(>F)"]
  sd_pooled <- sqrt(
    ((n_in - 1) * var(score[in_cluster], na.rm = TRUE) +
       (n_out - 1) * var(score[out_cluster], na.rm = TRUE)) /
      (n_in + n_out - 2))
  cohen_d <- (mean_in - mean_out) / sd_pooled
  data.frame(cluster = cluster_name, n_cluster = n_in, mean_cluster = mean_in,
              mean_rest = mean_out, cohen_d = cohen_d, p_val = p_val,
              stringsAsFactors = FALSE)
}

test_signature <- function(object, signature_name) {
  cell_subset <- Seurat::Idents(object)
  # Iteration over clusters is necessary for the original cluster-vs-rest
  # tests. There is no loop over signature definitions or parameter variants.
  res <- lapply(levels(cell_subset), score_cluster,
                score = object[[]][[paste0(signature_name, "1")]],
                cell_clusters = cell_subset,
                sample_id = object[[]][[settings$patient_column]])
  res <- do.call(rbind, res)
  # BH correction across clusters separately for EACH signature, as before.
  res$p_val_adj <- p.adjust(res$p_val, method = settings$adjustment)
  res$score <- signature_name
  res[, c("score", "cluster", "n_cluster", "mean_cluster", "mean_rest",
           "cohen_d", "p_val", "p_val_adj")]
}

# 4. Original effect-size heatmap styling -----------------------------------
effect_heatmap <- function(results) {
  results$border_size <- ifelse(results$p_val_adj < settings$significance, "sig", "ns")
  results$cluster <- factor(results$cluster, levels = unique(results$cluster))
  ggplot2::ggplot(results, ggplot2::aes(x = cluster, y = score, fill = cohen_d)) +
    ggplot2::geom_tile(ggplot2::aes(linewidth = border_size), color = "black") +
    ggplot2::scale_linewidth_manual(values = c(sig = 1.2, ns = 0)) +
    ggplot2::scale_fill_gradient2(low = "#060486", mid = "#FFFFFF", high = "#860604",
                                  midpoint = 0, name = "Cohen's d") +
    ggplot2::ggtitle("Effect Size (Cohen's d) Across Macrophage Subsets") +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(plot.background = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank(), panel.background = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 20),
      panel.grid.major = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank(), legend.position = "right")
}

# 5. Run explicitly: down signature, up signature, then mt filtering ---------
main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (!(length(args) %in% c(3L, 4L))) stop(
    "Usage: Rscript 02_single_cell_scores.R DEG_overview.csv prepared.rds output_dir [cluster_column]")
  if (!file.exists(args[1]) || !file.exists(args[2])) stop("Input file not found.")
  output_dir <- args[3]
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt")), add = TRUE)

  DEG_KDM <- read.csv(args[1])
  signatures <- select_signatures(DEG_KDM)
  cntdata <- readRDS(args[2])
  if (!inherits(cntdata, "Seurat")) stop("The supplied RDS must contain a prepared Seurat object.")
  metadata <- cntdata[[]]
  if (!all(c(settings$patient_column, "percent.mt") %in% names(metadata)))
    stop("The prepared Seurat object must contain Patient and percent.mt metadata.")
  if (!is.numeric(metadata$percent.mt) || anyNA(metadata$percent.mt) ||
      any(!is.finite(metadata$percent.mt)))
    stop("percent.mt must be finite and nonmissing, on the 0-100 percent scale.")
  if (anyNA(metadata[[settings$patient_column]]) ||
      any(!nzchar(as.character(metadata[[settings$patient_column]]))))
    stop("Patient identifiers must be nonmissing.")
  if (length(args) == 4L && !args[4] %in% names(metadata))
    stop("Cluster metadata column not found: ", args[4])

  # Preserve the input's default assay. Forcing RNA would change an analysis
  # that originally ran on another default assay. No assay/data is recomputed.
  assay_used <- Seurat::DefaultAssay(cntdata)
  options(lifecycle_verbosity = "warning")
  # Do not require an S3 method here: AddModuleScore is an ordinary function
  # in Seurat v4. The original scoring settings are passed explicitly below.

  # Separate calls are intentional: combining the two lists in ONE call can
  # change control-gene sampling. Each source call reset the default seed to 1.
  cntdata <- Seurat::AddModuleScore(
    cntdata, features = list(signatures$dnFC100), name = "dnFC100",
    ctrl = settings$ctrl, k = settings$k, nbin = settings$nbin,
    seed = settings$seed, search = settings$search, slot = settings$slot,
    assay = assay_used, pool = NULL)
  cntdata <- Seurat::AddModuleScore(
    cntdata, features = list(signatures$upFC100), name = "upFC100",
    ctrl = settings$ctrl, k = settings$k, nbin = settings$nbin,
    seed = settings$seed, search = settings$search, slot = settings$slot,
    assay = assay_used, pool = NULL)

  # Optional selection of a supplied annotation; no annotations are inferred.
  if (length(args) == 4L) Seurat::Idents(cntdata) <- cntdata[[]][[args[4]]]
  if (anyNA(Seurat::Idents(cntdata))) stop("Cluster identities must be nonmissing.")
  write.csv(cntdata[[]], file.path(output_dir, "scores.csv"))
  available <- rownames(cntdata[[assay_used]])
  gene_membership <- rbind(
    data.frame(signature = "dnFC100", rank = 1:100, Sym = signatures$dnFC100,
                in_assay = signatures$dnFC100 %in% available),
    data.frame(signature = "upFC100", rank = 1:100, Sym = signatures$upFC100,
                in_assay = signatures$upFC100 %in% available))
  write.csv(gene_membership, file.path(output_dir, "signature_genes.csv"), row.names = FALSE)
  dput(list(settings = settings, assay = assay_used,
             cluster_source = if (length(args) == 4L) args[4] else "existing Idents",
             contrast = contrast), file.path(output_dir, "settings_used.R"))

  # IMPORTANT: filter AFTER scoring, matching the original code exactly.
  fcnt <- subset(x = cntdata, subset = percent.mt < 7.5)
  cluster_sizes <- table(Seurat::Idents(fcnt))
  if (length(cluster_sizes) < 2L || any(cluster_sizes < 2L) ||
      length(unique(fcnt[[]][[settings$patient_column]])) < 2L)
    stop("The filtered object needs at least two clusters, two cells per cluster, and two patients.")
  if ("rest" %in% levels(Seurat::Idents(fcnt)))
    stop("A cluster is named 'rest', which conflicts with the original model's reference label.")
  write.csv(fcnt[[]], file.path(output_dir, "scores_filtered.csv"))

  down_results <- test_signature(fcnt, "dnFC100")
  up_results <- test_signature(fcnt, "upFC100")
  write.csv(rbind(down_results, up_results), file.path(output_dir, "GLMM_filtered.csv"),
             row.names = FALSE)

  p <- Seurat::VlnPlot(fcnt, features = c("dnFC1001", "upFC1001"),
                       pt.size = 0.5, log = FALSE, group.by = NULL, split.by = NULL, idents = NULL)
  ggplot2::ggsave(filename = file.path(output_dir, "KDMscore_filtered.pdf"),
                  plot = p, width = 16, height = 16)
  # Keep the source's CSV roundtrip before plotting effect sizes.
  plot_results <- read.csv(file.path(output_dir, "GLMM_filtered.csv"))
  p_down <- effect_heatmap(plot_results[plot_results$score == "dnFC100", ])
  p_up <- effect_heatmap(plot_results[plot_results$score == "upFC100", ])
  ggplot2::ggsave(filename = file.path(output_dir, "CD_filtered_dn.pdf"),
                  plot = p_down, width = 16, height = 7)
  ggplot2::ggsave(filename = file.path(output_dir, "CD_filtered_up.pdf"),
                  plot = p_up, width = 16, height = 7)
  if ("umap" %in% Seurat::Reductions(fcnt)) {
    p_umap <- Seurat::DimPlot(fcnt, reduction = "umap") + ggplot2::theme(aspect.ratio = 1)
    ggplot2::ggsave(filename = file.path(output_dir, "Umap_filtered.pdf"),
                    plot = p_umap, width = 20, height = 4)
  }
  message("Completed: ", normalizePath(output_dir))
  invisible(list(down = down_results, up = up_results))
}

if (sys.nframe() == 0L) main()
