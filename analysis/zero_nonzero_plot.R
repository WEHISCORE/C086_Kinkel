# Luke G’s zero-nonzero plot.
# Peter Hickey
# 2021-05-06

# Setup ------------------------------------------------------------------------

library(here)
library(SingleCellExperiment)
library(pheatmap)
library(edgeR)

sce <- readRDS(here("data", "SCEs", "C086_Kinkel.preprocessed.SCE.rds"))
sce <- sce[, order(sce$smchd1_genotype_updated, sce$sex, sce$mouse_number)]

sce_aggr <- scuttle::aggregateAcrossCells(
  sce,
  ids = colData(sce)[, c("genotype.mouse")],
  use.altexps = FALSE,
  use.dimred = FALSE,
  use.assay.type = c("UMI_counts", "read_counts"))

# Some useful colours
sex_colours <- setNames(
  unique(sce$sex_colours),
  unique(names(sce$sex_colours)))
smchd1_genotype_updated_colours <- setNames(
  unique(sce$smchd1_genotype_updated_colours),
  unique(names(sce$smchd1_genotype_updated_colours)))
mouse_number_colours <- setNames(
  unique(sce$mouse_number_colours),
  unique(names(sce$mouse_number_colours)))

zeroNonZeroPlot <- function(sce, exprs_values) {
  mat <- assay(sce, exprs_values)
  mat <- mat[rowSums(mat) > 0, ]
  ave_lcpm <- cpmByGroup(mat, log = TRUE)
  I <- as(mat > 0, "dgCMatrix")

  # Re-order data by aveLogCpm
  o <- order(ave_lcpm, decreasing = TRUE)
  ave_lcpm <- ave_lcpm[o]
  I <- I[o, ]
  mat <- mat[o, ]

  keep_exprs <- filterByExpr(
    mat,
    group = sce$smchd1_genotype_updated,
    min.count = 5)

  pheatmap(
    mat = as.matrix(I),
    color = c("yellow", "black"),
    scale = "none",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = FALSE,
    legend = FALSE,
    annotation_col = as.data.frame(
      colData(sce)[, c("mouse_number", "sex", "smchd1_genotype_updated")]),
    annotation_row = data.frame(
      aveLogCPM = ave_lcpm,
      filterByExpr = as.factor(keep_exprs),
      row.names = rownames(I)),
    annotation_colors = list(
      "mouse_number" = mouse_number_colours,
      "sex" = sex_colours,
      "smchd1_genotype_updated" = smchd1_genotype_updated_colours,
      "filterByExpr" = c("FALSE" = "white", "TRUE" = "forestgreen")),
    main = exprs_values,
    fontsize = 8)
}

# Plot -------------------------------------------------------------------------

dir.create(here("output/figures/for_Terry"), recursive = TRUE)
pdf(
  here("output/figures/for_Terry/zero_nonzero_plot.pdf"),
  width = 6,
  height = 9)
zeroNonZeroPlot(sce, "read_counts")
zeroNonZeroPlot(sce, "UMI_counts")
zeroNonZeroPlot(sce_aggr, "read_counts")
zeroNonZeroPlot(sce_aggr, "UMI_counts")
dev.off()
