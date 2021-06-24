# Deterministic imputation of count matrix with diagnostic plots.
# Peter Hickey
# 2021-06-24

# Setup ------------------------------------------------------------------------

library(here)
library(SingleCellExperiment)
library(pheatmap)
library(edgeR)
library(RUVSeq)

source(here("code/helper_functions.R"))

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

# Output directory
dir.create(here("output/figures/for_Terry"), recursive = TRUE)

# NOTE: This function assumes that the genes have already been filtered.
diagnosticPlots <- function(sce, exprs_values) {

  # Prepare data
  dgel <- DGEList(
    counts = assay(sce, exprs_values),
    samples = colData(sce),
    group = factor(sce$smchd1_genotype_updated))

  # DE analysis
  dgel <- calcNormFactors(dgel, method = "upperquartile")
  design <- model.matrix(~0 + smchd1_genotype_updated + sex, dgel$samples)
  dgel <- estimateDisp(dgel, design)
  dgeglm <- glmFit(dgel, design)
  contrasts <- makeContrasts(
    smchd1_genotype_updatedHet - smchd1_genotype_updatedDel,
    levels = design)
  dgelrt <- glmLRT(dgeglm, contrast = contrasts)
  tt <- topTags(dgelrt, n = Inf)

  # Plot
  par(mfrow = c(2, 2))
  plotPCA(
    cpm(dgel),
    col = dgel$samples$smchd1_genotype_updated_colours,
    main = exprs_values)
  plotRLE(
    cpm(dgel),
    outline = FALSE,
    col = dgel$samples$smchd1_genotype_updated_colours,
    las = 2,
    main = exprs_values)
  legend(
    "bottomright",
    fill = unique(dgel$samples$smchd1_genotype_updated_colours),
    legend = unique(dgel$samples$smchd1_genotype_updated),
    cex = 0.4)
  plotBCV(dgel, main = exprs_values)
  hist(tt$table$PValue, breaks = 0:20 / 20, xlab = "P", main = exprs_values)
}

# Filtering --------------------------------------------------------------------

sce_reads_standard <- filterCountMatrix(
  sce,
  "read_counts",
  min.count = 5,
  method = "standard")
sce_reads_careful <- filterCountMatrix(
  sce,
  "read_counts",
  min.count = 5,
  method = "careful")
sce_umis_standard <- filterCountMatrix(
  sce,
  "UMI_counts",
  min.count = 5,
  method = "standard")
sce_umis_careful <- filterCountMatrix(
  sce,
  "UMI_counts",
  min.count = 5,
  method = "careful")

# Imputation -------------------------------------------------------------------

sce_reads_standard_all <- deterministicImputation(
  sce = sce_reads_standard,
  exprs_values = "read_counts",
  method = "all")
sce_reads_standard_group <- deterministicImputation(
  sce = sce_reads_standard,
  exprs_values = "read_counts",
  method = "group-specific")
sce_reads_careful_all <- deterministicImputation(
  sce = sce_reads_careful,
  exprs_values = "read_counts",
  method = "all")
sce_reads_careful_group <- deterministicImputation(
  sce = sce_reads_careful,
  exprs_values = "read_counts",
  method = "group-specific")
sce_umis_standard_all <- deterministicImputation(
  sce = sce_umis_standard,
  exprs_values = "UMI_counts",
  method = "all")
sce_umis_standard_group <- deterministicImputation(
  sce = sce_umis_standard,
  exprs_values = "UMI_counts",
  method = "group-specific")
sce_umis_careful_all <- deterministicImputation(
  sce = sce_umis_careful,
  exprs_values = "UMI_counts",
  method = "all")
sce_umis_careful_group <- deterministicImputation(
  sce = sce_umis_careful,
  exprs_values = "UMI_counts",
  method = "group-specific")

# Plots: Individual replicates, read counts ------------------------------------

pdf(
  here("output/figures/for_Terry/imputation_diagnostics.read_counts.standard_filtering.pdf"),
  width = 9,
  height = 9)
diagnosticPlots(sce_reads_standard, "read_counts")
dev.off()

pdf(
  here("output/figures/for_Terry/imputation_diagnostics.read_counts.careful_filtering.pdf"),
  width = 9,
  height = 9)
diagnosticPlots(sce_reads_careful, "read_counts")
dev.off()

# Plots: Individual replicates, all imputed read counts ------------------------

pdf(
  here("output/figures/for_Terry/imputation_diagnostics.all_imputed_read_counts.standard_filtering.pdf"),
  width = 9,
  height = 9)
diagnosticPlots(sce_reads_standard_all, "imputed_read_counts")
dev.off()

pdf(
  here("output/figures/for_Terry/imputation_diagnostics.all_imputed_read_counts.careful_filtering.pdf"),
  width = 9,
  height = 9)
diagnosticPlots(sce_reads_careful_all, "imputed_read_counts")
dev.off()

# Plots: Individual replicates, group imputed read counts ----------------------

pdf(
  here("output/figures/for_Terry/imputation_diagnostics.group_imputed_read_counts.standard_filtering.pdf"),
  width = 9,
  height = 9)
diagnosticPlots(sce_reads_standard_group, "imputed_read_counts")
dev.off()

pdf(
  here("output/figures/for_Terry/imputation_diagnostics.group_imputed_read_counts.careful_filtering.pdf"),
  width = 9,
  height = 9)
diagnosticPlots(sce_reads_careful_group, "imputed_read_counts")
dev.off()

# Plots: Individual replicates, UMI counts -------------------------------------

pdf(
  here("output/figures/for_Terry/imputation_diagnostics.UMI_counts.standard_filtering.pdf"),
  width = 9,
  height = 9)
diagnosticPlots(sce_umis_standard, "UMI_counts")
dev.off()

pdf(
  here("output/figures/for_Terry/imputation_diagnostics.UMI_counts.careful_filtering.pdf"),
  width = 9,
  height = 9)
diagnosticPlots(sce_umis_careful, "UMI_counts")
dev.off()

# Plots: Individual replicates, all imputed UMI counts -------------------------

pdf(
  here("output/figures/for_Terry/imputation_diagnostics.all_imputed_UMI_counts.standard_filtering.pdf"),
  width = 9,
  height = 9)
diagnosticPlots(sce_umis_standard_all, "imputed_UMI_counts")
dev.off()

pdf(
  here("output/figures/for_Terry/imputation_diagnostics.all_imputed_UMI_counts.careful_filtering.pdf"),
  width = 9,
  height = 9)
diagnosticPlots(sce_umis_careful_all, "imputed_UMI_counts")
dev.off()

# Plots: Individual replicates, group imputed UMI counts -----------------------

pdf(
  here("output/figures/for_Terry/imputation_diagnostics.group_imputed_UMI_counts.standard_filtering.pdf"),
  width = 9,
  height = 9)
diagnosticPlots(sce_umis_standard_group, "imputed_UMI_counts")
dev.off()

pdf(
  here("output/figures/for_Terry/imputation_diagnostics.group_imputed_UMI_counts.careful_filtering.pdf"),
  width = 9,
  height = 9)
diagnosticPlots(sce_umis_careful_group, "imputed_UMI_counts")
dev.off()
