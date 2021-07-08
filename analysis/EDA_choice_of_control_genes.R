# Comparison of running RUVs with all genes vs. house-keeping genes
# Peter Hickey
# 2021-07-08

# One the basis of results from RUV.R, I am using:
# - Read counts rather than UMI counts
# - Summed counts
# - Raw counts rather than imputed counts
# - Standard filtering
# - Group by genotype:sex
# - K = 3

# Setup ------------------------------------------------------------------------

library(here)
library(SingleCellExperiment)
library(RUVSeq)
library(edgeR)

source(here("code/helper_functions.R"))

sce <- readRDS(here("data", "SCEs", "C086_Kinkel.preprocessed.SCE.rds"))
sce <- sce[, order(sce$smchd1_genotype_updated, sce$sex, sce$mouse_number)]

# NOTE: Filter to only retain protein coding genes.
sce <- sce[any(grepl("protein_coding", rowData(sce)$ENSEMBL.GENEBIOTYPE)), ]

# Drop WT sample
sce <- sce[, sce$smchd1_genotype_updated != "WT"]
colData(sce) <- droplevels(colData(sce))

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

# Output directory
dir.create(here("output/figures/for_Terry"), recursive = TRUE)

# Mouse housekeeping genes from Ramyar (via Luke).
hk <- readRDS(here("data/scHK_mouse.rds"))

# Filtering --------------------------------------------------------------------

sce_aggr_reads_standard <- filterCountMatrix(
  sce_aggr,
  "read_counts",
  min.count = 5,
  method = "standard")

# RUV --------------------------------------------------------------------------

sce <- sce_aggr_reads_standard
exprs_values <- "read_counts"
g <- "genotype_sex"
dgel <- DGEList(
  counts = as.matrix(assay(sce, exprs_values)),
  samples = colData(sce),
  group = factor(sce$smchd1_genotype_updated),
  genes = flattenDF(rowData(sce)))
dgel <- calcNormFactors(dgel, method = "upperquartile")
k <- 3
differences <- makeGroups(dgel$samples[[g]])
set_hk <- RUVs(
  x = betweenLaneNormalization(dgel$counts, "upper"),
  cIdx = toupper(rownames(dgel)) %in% hk,
  k = k,
  scIdx = differences)
set_all <- RUVs(
  x = betweenLaneNormalization(dgel$counts, "upper"),
  cIdx = rownames(dgel),
  k = k,
  scIdx = differences)

par(mfrow = c(1, 2))
plotRLE(
  set_hk$normalizedCounts,
  outline = FALSE,
  col = dgel$samples$smchd1_genotype_updated_colours,
  las = 2,
  main = "HK",
  cex = 0.5)
plotRLE(
  set_all$normalizedCounts,
  outline = FALSE,
  col = dgel$samples$smchd1_genotype_updated_colours,
  las = 2,
  main = "All",
  cex = 0.5)

par(mfrow = c(1, 2))
plotPCA(
  set_hk$normalizedCounts,
  col = dgel$samples$smchd1_genotype_updated_colours,
  main = "HK")
plotPCA(
  set_all$normalizedCounts,
  col = dgel$samples$smchd1_genotype_updated_colours,
  main = "All")

par(mfrow = c(2, 2))
plot(
  set_hk$W[, 1],
  set_all$W[, 1],
  main = "W1",
  xlab = "HK",
  ylab = "All",
  sub = paste0("cor = ", round(cor(set_hk$W[, 1], set_all$W[, 1]), 2)),
  xlim = range(c(set_hk$W[, 1], set_all$W[, 1])),
  ylim = range(c(set_hk$W[, 1], set_all$W[, 1])))
abline(a = 0, b = 1, lty = 2)
plot(
  set_hk$W[, 2],
  set_all$W[, 2],
  main = "W1",
  xlab = "HK",
  ylab = "All",
  sub = paste0("cor = ", round(cor(set_hk$W[, 2], set_all$W[, 2]), 2)),
  xlim = range(c(set_hk$W[, 2], set_all$W[, 2])),
  ylim = range(c(set_hk$W[, 2], set_all$W[, 2])))
abline(a = 0, b = 1, lty = 2)
plot(
  set_hk$W[, 3],
  set_all$W[, 3],
  main = "W1",
  xlab = "HK",
  ylab = "All",
  sub = paste0("cor = ", round(cor(set_hk$W[, 3], set_all$W[, 3]), 2)),
  xlim = range(c(set_hk$W[, 3], set_all$W[, 3])),
  ylim = range(c(set_hk$W[, 3], set_all$W[, 3])))
abline(a = 0, b = 1, lty = 2)

sets <- list(HK = set_hk, All = set_all)
lapply(names(sets), function(n) {
  set <- sets[[n]]
  design <- model.matrix(
    ~0 + smchd1_genotype_updated + sex + set$W,
    dgel$samples)
  colnames(design) <- sub(
    "set\\$|smchd1_genotype_updated|sex",
    "",
    colnames(design))
  dgel <- estimateDisp(dgel, design)
  dgeglm <- glmFit(dgel, design)
  contrasts <- makeContrasts(Het - Del, levels = design)
  dgelrt <- glmLRT(dgeglm, contrast = contrasts)
  tt <- topTags(dgelrt, n = Inf)
  print(summary(decideTests(dgelrt)))

  # Plots
  par(mfrow = c(2, 3))
  plotPCA(
    set$normalizedCounts,
    col = dgel$samples$smchd1_genotype_updated_colours,
    main = n)
  plotRLE(
    set$normalizedCounts,
    outline = FALSE,
    col = dgel$samples$smchd1_genotype_updated_colours,
    las = 2,
    main = n)
  legend(
    "bottomright",
    fill = unique(dgel$samples$smchd1_genotype_updated_colours),
    legend = unique(dgel$samples$smchd1_genotype_updated),
    cex = 0.4)
  plotBCV(dgel, main = n)
  hist(
    tt$table$PValue,
    breaks = 0:20 / 20,
    xlab = "PValue",
    main = n)
  plotMD(dgelrt)
  abline(h = 0, lty = 2, col = "dodgerBlue")
  plotMD(
    dgelrt,
    status = ifelse(rownames(dgelrt) == "Smchd1", "Smchd1", "Other"),
    hl.col = "orange")
  abline(h = 0, lty = 2, col = "dodgerBlue")
})
