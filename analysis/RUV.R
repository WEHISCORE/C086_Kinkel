# RUVs
# Peter Hickey
# 2021-06-24

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

runRUVs <- function(sce, exprs_values, g, ks) {
  message("Constructing DGEList")
  dgel <- DGEList(
    counts = as.matrix(assay(sce, exprs_values)),
    samples = colData(sce),
    group = factor(sce$smchd1_genotype_updated),
    genes = flattenDF(rowData(sce)))
  dgel <- calcNormFactors(dgel, method = "upperquartile")

  message("Running RUVs with various k")
  differences <- makeGroups(dgel$samples[[g]])
  ruvs <- lapply(ks, function(k) {
    # TODO: Seems you need to pass normalized counts to RUVs; Check with Luke.
    RUVs(
      # TODO: Not sure why, but can't simply do `x = cpm(dgel)`
      x = betweenLaneNormalization(dgel$counts, "upper"),
      # TODO: Investigate using mouse house keeping control genes (see Peter.Rmd
      #       from Luke).
      cIdx = rownames(dgel),
      k = k,
      scIdx = differences)
  })

  message("RLE coloured by genotype")
  nr <- nc <- ceiling(sqrt(1 + length(ks)))
  par(mfrow = c(nr, nc))
  plotRLE(
    cpm(dgel),
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
    cpm(dgel),
    col = dgel$samples$smchd1_genotype_updated_colours,
    main = "Raw")
  lapply(ks, function(k) {
    set <- ruvs[[as.character(k)]]
    plotPCA(
      set$normalizedCounts,
      col = dgel$samples$smchd1_genotype_updated_colours,
      main = paste0("RUVs: k = ", k))
  })

  message("PCA coloured by sex")
  par(mfrow = c(nr, nc))
  plotPCA(
    cpm(dgel),
    col = dgel$samples$sex_colours,
    main = "Raw")
  lapply(ks, function(k) {
    set <- ruvs[[as.character(k)]]
    plotPCA(
      set$normalizedCounts,
      col = dgel$samples$sex_colours,
      main = paste0("RUVs: k = ", k))
  })

  message("PCA coloured by mouse")
  par(mfrow = c(nr, nc))
  plotPCA(
    cpm(dgel),
    col = dgel$samples$mouse_number_colours,
    main = "Raw")
  lapply(ks, function(k) {
    set <- ruvs[[as.character(k)]]
    plotPCA(
      set$normalizedCounts,
      col = dgel$samples$mouse_number_colours,
      main = paste0("RUVs: k = ", k))
  })

  message("Diagnostic plots")
  lapply(ks, function(k) {
    message("\tk = ", k)
    set <- ruvs[[as.character(k)]]
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

    # Plots
    par(mfrow = c(2, 3))
    plotPCA(
      set$normalizedCounts,
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
      tt$table$PValue,
      breaks = 0:20 / 20,
      xlab = "PValue",
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

ks <- seq(1, 15, 2)
names(ks) <- ks

pdf(
  here("output/figures/for_Terry/RUV.read_counts.standard_filtering.mouse_number.pdf"),
  width = 9,
  height = 9)
runRUVs(
  sce = sce_reads_standard,
  exprs_values = "read_counts",
  g = "mouse_number",
  ks = ks)
dev.off()

pdf(
  here("output/figures/for_Terry/RUV.read_counts.standard_filtering.genotype_sex.pdf"),
  width = 9,
  height = 9)
runRUVs(
  sce = sce_reads_standard,
  exprs_values = "read_counts",
  g = "genotype_sex",
  ks = ks)
dev.off()

pdf(
  here("output/figures/for_Terry/RUV.read_counts.standard_filtering.smchd1_genotype_updated.pdf"),
  width = 9,
  height = 9)
runRUVs(
  sce = sce_reads_standard,
  exprs_values = "read_counts",
  g = "smchd1_genotype_updated",
  ks = ks)
dev.off()

pdf(
  here("output/figures/for_Terry/RUV.read_counts.careful_filtering.mouse_number.pdf"),
  width = 9,
  height = 9)
runRUVs(
  sce = sce_reads_careful,
  exprs_values = "read_counts",
  g = "mouse_number",
  ks = ks)
dev.off()

pdf(
  here("output/figures/for_Terry/RUV.read_counts.careful_filtering.genotype_sex.pdf"),
  width = 9,
  height = 9)
runRUVs(
  sce = sce_reads_careful,
  exprs_values = "read_counts",
  g = "genotype_sex",
  ks = ks)
dev.off()

pdf(
  here("output/figures/for_Terry/RUV.read_counts.careful_filtering.genotype.pdf"),
  width = 9,
  height = 9)
runRUVs(
  sce = sce_reads_careful,
  exprs_values = "read_counts",
  g = "genotype",
  ks = ks)
dev.off()

# Plots: Individual replicates, all imputed read counts ------------------------

# TODO
ks <- seq(1, 15, 2)
names(ks) <- ks


# Plots: Individual replicates, group imputed read counts ----------------------

# TODO
ks <- seq(1, 15, 2)
names(ks) <- ks


# Plots: Individual replicates, UMI counts -------------------------------------

ks <- seq(1, 15, 2)
names(ks) <- ks


# Plots: Individual replicates, all imputed UMI counts -------------------------

# TODO
ks <- seq(1, 15, 2)
names(ks) <- ks

# Plots: Individual replicates, group imputed UMI counts -----------------------

# TODO
ks <- seq(1, 15, 2)
names(ks) <- ks

# Plots: Summed replicates, read counts ----------------------------------------

# TODO
ks <- seq(1, 7, 1)
names(ks) <- ks

# Plots: Summed replicates, all imputed read counts ----------------------------

# TODO
ks <- seq(1, 7, 1)
names(ks) <- ks

# Plots: Summed replicates, group imputed read counts --------------------------

# TODO
ks <- seq(1, 7, 1)
names(ks) <- ks

# Plots: Summed replicates, UMI counts -----------------------------------------

# TODO
ks <- seq(1, 7, 1)
names(ks) <- ks

# Plots: Summed replicates, all imputed UMI counts -----------------------------

# TODO
ks <- seq(1, 7, 1)
names(ks) <- ks

# Plots: Summed replicates, group imputed UMI counts ---------------------------

# TODO
ks <- seq(1, 7, 1)
names(ks) <- ks
