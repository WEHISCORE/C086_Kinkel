# Deterministic imputation of count matrix with diagnostic plots.
# Peter Hickey
# 2021-05-06

# Setup ------------------------------------------------------------------------

library(here)
library(SingleCellExperiment)
library(pheatmap)
library(edgeR)
library(RUVSeq)

source(here("code/helper_functions.R"))

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

deterministicImputation <- function(sce, exprs_values) {
  mat <- assay(sce, exprs_values)
  pi <- rowSums(mat) / sum(mat)
  nonzero_idx <- which(mat > 0)
  imputed <- as.matrix(pi) %*%
    matrix(colSums(mat), ncol = ncol(mat), dimnames = list(NULL, colnames(mat)))
  imputed[nonzero_idx] <- mat[nonzero_idx]
  assay(sce, paste0("imputed_", exprs_values)) <- imputed
  sce
}

diagnosticPlots <- function(sce, exprs_values) {

  # Prepare data
  dgel <- DGEList(
    counts = as.matrix(assay(sce, exprs_values)),
    samples = colData(sce),
    group = factor(sce$smchd1_genotype_updated),
    genes = flattenDF(rowData(sce)))
  dgel <- dgel[rowSums(dgel$counts) > 0, ]
  keep_exprs <- filterByExpr(
    dgel,
    group = dgel$samples$smchd1_genotype_updated,
    min.count = 5)
  dgel <- dgel[keep_exprs, , keep.lib.sizes = FALSE]
  dgel <- calcNormFactors(dgel)
  # TODO: With or without mouse_number?
  design <- model.matrix(~0 + smchd1_genotype_updated + sex, dgel$samples)
  dgel <- estimateDisp(dgel, design)
  dgeglm <- glmQLFit(dgel, design)
  contrasts <- makeContrasts(
    smchd1_genotype_updatedDel - smchd1_genotype_updatedHet,
    levels = design)
  dgelrt <- glmQLFTest(dgeglm, contrast = contrasts)
  tt <- topTags(dgelrt, n = Inf)

  # Plot
  par(mfrow = c(2, 2))
  plotPCA(
    cpm(dgel$counts),
    isLog = FALSE,
    col = dgel$samples$smchd1_genotype_updated_colours,
    main = exprs_values)
  plotRLE(
    dgel$counts,
    outline = FALSE,
    col = dgel$samples$smchd1_genotype_updated_colours,
    las = 2)
  legend(
    "bottomright",
    fill = unique(dgel$samples$smchd1_genotype_updated_colours),
    legend = unique(dgel$samples$smchd1_genotype_updated),
    cex = 0.4)
  plotBCV(dgel, main = exprs_values)
  hist(tt$table$FDR, breaks = 0:20 / 20, xlab = "FDR", main = exprs_values)
}

# Imputation -------------------------------------------------------------------

sce <- deterministicImputation(sce, "read_counts")
sce <- deterministicImputation(sce, "UMI_counts")
sce_aggr <- deterministicImputation(sce_aggr, "read_counts")
sce_aggr <- deterministicImputation(sce_aggr, "UMI_counts")

# Plot -------------------------------------------------------------------------

dir.create(here("output/figures/for_Terry"), recursive = TRUE)
pdf(
  here("output/figures/for_Terry/imputation_diagnostics.pdf"),
  width = 9,
  height = 9)
diagnosticPlots(sce, "read_counts")
diagnosticPlots(sce, "imputed_read_counts")
diagnosticPlots(sce, "UMI_counts")
diagnosticPlots(sce, "imputed_UMI_counts")
diagnosticPlots(sce_aggr, "read_counts")
diagnosticPlots(sce_aggr, "imputed_read_counts")
diagnosticPlots(sce_aggr, "UMI_counts")
diagnosticPlots(sce_aggr, "imputed_UMI_counts")
dev.off()
