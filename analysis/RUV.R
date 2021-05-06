# RUVs
# Peter Hickey
# 2021-05-06

# Setup ------------------------------------------------------------------------

library(here)
library(SingleCellExperiment)
library(RUVSeq)
library(edgeR)

source(here("code/helper_functions.R"))

sce <- readRDS(here("data", "SCEs", "C086_Kinkel.preprocessed.SCE.rds"))
sce <- sce[, order(sce$smchd1_genotype_updated, sce$sex, sce$mouse_number)]
sce$genotype_sex <- paste0(sce$smchd1_genotype_updated, "_", sce$sex)

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

# TODO: In calls to plotPCA(), should counts be (log)cpm-ed?
runRUVs <- function(sce, exprs_values, g, ks) {
  message("Constructing DGEList")
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

  message("Running RUVs with various k")
  differences <- makeGroups(dgel$samples[[g]])
  ruvs <- lapply(ks, function(k) {
    RUVs(dgel$counts, rownames(dgel), k = k, scIdx = differences)
  })

  message("RLE coloured by genotype")
  nr <- nc <- ceiling(sqrt(1 + length(ks)))
  par(mfrow = c(nr, nc))
  plotRLE(
    dgel$counts,
    outline = FALSE,
    col = dgel$samples$smchd1_genotype_updated_colours,
    las = 2,
    main = "Raw",
    cex = 0.5)
  lapply(ks, function(k) {
    set <- ruvs[[as.character(k)]]
    plotRLE(
      set$normalizedCounts,
      outline = FALSE,
      col = dgel$samples$smchd1_genotype_updated_colours,
      las = 2,
      main = paste0("RUVs: k = ", k),
      cex = 0.5)
  })

  message("PCA coloured by genotype")
  par(mfrow = c(nr, nc))

  plotPCA(
    cpm(dgel$counts),
    isLog = FALSE,
    col = dgel$samples$smchd1_genotype_updated_colours,
    main = "Raw")
  lapply(ks, function(k) {
    set <- ruvs[[as.character(k)]]
    plotPCA(
      set$normalizedCounts,
      isLog = FALSE,
      col = dgel$samples$smchd1_genotype_updated_colours,
      main = paste0("RUVs: k = ", k))
  })

  message("PCA coloured by sex")
  par(mfrow = c(nr, nc))
  plotPCA(
    cpm(dgel$counts),
    isLog = FALSE,
    col = dgel$samples$sex_colours,
    main = "Raw")
  lapply(ks, function(k) {
    set <- ruvs[[as.character(k)]]
    plotPCA(
      set$normalizedCounts,
      isLog = FALSE,
      col = dgel$samples$sex_colours,
      main = paste0("RUVs: k = ", k))
  })

  message("PCA coloured by mouse")
  par(mfrow = c(nr, nc))
  plotPCA(
    cpm(dgel$counts),
    isLog = FALSE,
    col = dgel$samples$mouse_number_colours,
    main = "Raw")
  lapply(ks, function(k) {
    set <- ruvs[[as.character(k)]]
    plotPCA(
      set$normalizedCounts,
      isLog = FALSE,
      col = dgel$samples$mouse_number_colours,
      main = paste0("RUVs: k = ", k))
  })

  message("Diagnostic plots")
  lapply(ks, function(k) {
    message("\tk = ", k)
    set <- ruvs[[as.character(k)]]
    # TODO: With or without mouse_number?
    design <- model.matrix(
      ~0 + smchd1_genotype_updated + sex + set$W,
      dgel$samples)
    colnames(design) <- sub(
      "set\\$|smchd1_genotype_updated|sex",
      "",
      colnames(design))
    dgel <- estimateDisp(dgel, design)
    dgeglm <- glmQLFit(dgel, design)
    contrasts <- makeContrasts(Del - Het, levels = design)
    dgelrt <- glmQLFTest(dgeglm, contrast = contrasts)
    tt <- topTags(dgelrt, n = Inf)

    # Plots
    par(mfrow = c(2, 3))
    plotPCA(
      set$normalizedCounts,
      isLog = FALSE,
      col = dgel$samples$smchd1_genotype_updated_colours,
      main = paste0("RUVs: k = ", k))
    plotRLE(
      set$normalizedCounts,
      outline = FALSE,
      col = dgel$samples$smchd1_genotype_updated_colours,
      las = 2,
      main = paste0("RUVs: k = ", k))
    legend(
      "bottomright",
      fill = unique(dgel$samples$smchd1_genotype_updated_colours),
      legend = unique(dgel$samples$smchd1_genotype_updated),
      cex = 0.4)
    plotBCV(dgel, main = paste0("RUVs: k = ", k))
    hist(
      tt$table$FDR,
      breaks = 0:20 / 20,
      xlab = "FDR",
      main = paste0("RUVs: k = ", k))
    plotMD(dgelrt)
    abline(h = 0, lty = 2, col = "dodgerBlue")
    plotMD(
      dgelrt,
      status = ifelse(rownames(dgelrt) == "Smchd1", "Smchd1", "Other"),
      hl.col = "orange")
    abline(h = 0, lty = 2, col = "dodgerBlue")
  })
}

dir.create(here("output/figures/for_Terry"), recursive = TRUE)

# Individual tech reps ---------------------------------------------------------

ks <- seq(1, 15, 2)
names(ks) <- ks
# Read counts
pdf(
  here("output/figures/for_Terry/RUVs.read_counts.mouse_number.pdf"),
  width = 9,
  height = 9)
runRUVs(
  sce = sce,
  exprs_values = "read_counts",
  g = "mouse_number",
  ks = ks)
dev.off()
pdf(
  here("output/figures/for_Terry/RUVs.read_counts.genotype_sex.pdf"),
  width = 9,
  height = 9)
runRUVs(
  sce = sce,
  exprs_values = "read_counts",
  g = "genotype_sex",
  ks = ks)
dev.off()
pdf(
  here("output/figures/for_Terry/RUVs.read_counts.smchd1_genotype_updated.pdf"),
  width = 9,
  height = 9)
runRUVs(
  sce = sce,
  exprs_values = "read_counts",
  g = "smchd1_genotype_updated",
  ks = ks)
dev.off()

# UMI counts
pdf(
  here("output/figures/for_Terry/RUVs.UMI_counts.mouse_number.pdf"),
  width = 9,
  height = 9)
runRUVs(
  sce = sce,
  exprs_values = "UMI_counts",
  g = "mouse_number",
  ks = ks)
dev.off()
pdf(
  here("output/figures/for_Terry/RUVs.UMI_counts.genotype_sex.pdf"),
  width = 9,
  height = 9)
runRUVs(
  sce = sce,
  exprs_values = "UMI_counts",
  g = "genotype_sex",
  ks = ks)
dev.off()
pdf(
  here("output/figures/for_Terry/RUVs.UMI_counts.smchd1_genotype_updated.pdf"),
  width = 9,
  height = 9)
runRUVs(
  sce = sce,
  exprs_values = "UMI_counts",
  g = "smchd1_genotype_updated",
  ks = ks)
dev.off()

# TODO: Imputed counts; do imputed counts need to be rounded?

# Summed tech reps -------------------------------------------------------------

ks <- seq(1, 7, 1)
names(ks) <- ks
# Read counts
pdf(
  here("output/figures/for_Terry/RUVs.read_counts.genotype_sex.summedTechReps.pdf"),
  width = 9,
  height = 9)
runRUVs(
  sce = sce_aggr,
  exprs_values = "read_counts",
  g = "genotype_sex",
  ks = ks)
dev.off()
pdf(
  here("output/figures/for_Terry/RUVs.read_counts.smchd1_genotype_updated.summedTechReps.pdf"),
  width = 9,
  height = 9)
runRUVs(
  sce = sce_aggr,
  exprs_values = "read_counts",
  g = "smchd1_genotype_updated",
  ks = ks)
dev.off()

# UMI_counts
pdf(
  here("output/figures/for_Terry/RUVs.UMI_counts.genotype_sex.summedTechReps.pdf"),
  width = 9,
  height = 9)
runRUVs(
  sce = sce_aggr,
  exprs_values = "UMI_counts",
  g = "genotype_sex",
  ks = ks)
dev.off()
pdf(
  here("output/figures/for_Terry/RUVs.UMI_counts.smchd1_genotype_updated.summedTechReps.pdf"),
  width = 9,
  height = 9)
runRUVs(
  sce = sce_aggr,
  exprs_values = "UMI_counts",
  g = "smchd1_genotype_updated",
  # NOTE: Can't estimate k > 8.
  ks = ks)
dev.off()

# TODO: Imputed counts.

# TODO: Some of the boxplots are basically squashed, e.g.:
#       - RUVs.read_counts.smchd1_genotype_updated.pdf
#       - RUVs.UMI_counts.mouse_number.pdf (extreme)
#       - RUVs.UMI_counts.genotype_sex.pdf (extreme)
#       - RUVs.UMI_counts.smchd1_genotype_updated.pdf (extreme)
#       - RUVs.UMI_counts.genotype_sex.summedTechReps.pdf
#       - RUVs.UMI_counts.smchd1_genotype_updated.summedTechReps.pdf (extreme)
