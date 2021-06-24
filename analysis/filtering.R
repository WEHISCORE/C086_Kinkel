# 'Standard' vs. 'careful' filtering of count matrix
# Peter Hickey
# 2021-06-24

# Setup ------------------------------------------------------------------------

library(here)
library(SingleCellExperiment)
library(edgeR)

sce <- readRDS(here("data", "SCEs", "C086_Kinkel.preprocessed.SCE.rds"))
sce <- sce[, order(sce$smchd1_genotype_updated, sce$sex, sce$mouse_number)]

# NOTE: Filter to only retain protein coding genes.
sce <- sce[any(grepl("protein_coding", rowData(sce)$ENSEMBL.GENEBIOTYPE)), ]

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

# Read counts ------------------------------------------------------------------

y <- assay(sce, "read_counts")
min.count <- 10

keep_standard <- filterByExpr(
  y,
  group = sce$smchd1_genotype_updated,
  min.count = min.count)

keep_del <- rowSums(y[, sce$smchd1_genotype_updated == "Del"] >= min.count) >=
  sum(sce$smchd1_genotype_updated == "Del")
keep_het <- rowSums(y[, sce$smchd1_genotype_updated == "Het"] >= min.count) >=
  sum(sce$smchd1_genotype_updated == "Het")
keep_wt <- rowSums(y[, sce$smchd1_genotype_updated == "WT"] >= min.count) >=
  sum(sce$smchd1_genotype_updated == "WT")
keep_careful <- keep_del | keep_het | keep_wt

table(standard = keep_standard, careful = keep_careful)

# UMI counts -------------------------------------------------------------------

y <- assay(sce, "UMI_counts")
min.count <- 5

keep_standard <- filterByExpr(
  y,
  group = sce$smchd1_genotype_updated,
  min.count = min.count)

keep_del <- rowSums(y[, sce$smchd1_genotype_updated == "Del"] >= min.count) >=
  sum(sce$smchd1_genotype_updated == "Del")
keep_het <- rowSums(y[, sce$smchd1_genotype_updated == "Het"] >= min.count) >=
  sum(sce$smchd1_genotype_updated == "Het")
keep_wt <- rowSums(y[, sce$smchd1_genotype_updated == "WT"] >= min.count) >=
  sum(sce$smchd1_genotype_updated == "WT")
keep_careful <- keep_del | keep_het | keep_wt

table(standard = keep_standard, careful = keep_careful)

# Summary ----------------------------------------------------------------------

# Roughly, keep_careful is a subset of keep_standard
