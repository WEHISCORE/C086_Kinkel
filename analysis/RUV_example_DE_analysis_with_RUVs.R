# DE analysis following lessons learned from discussions with Terry and Luke
# (explored in EDA_*.R scripts and EDA_summary.md).
# Peter Hickey
# 2021-07-08

# Setup ------------------------------------------------------------------------

library(here)
library(SingleCellExperiment)
library(scater)
library(RUVSeq)
library(edgeR)
library(pheatmap)

source(here("code/helper_functions.R"))

sce <- readRDS(here("data", "SCEs", "C086_Kinkel.preprocessed.SCE.rds"))
# NOTE: Reorder samples for plots (e.g., RLE plots).
sce <- sce[, order(sce$smchd1_genotype_updated, sce$sex, sce$mouse_number)]

# Restrict analysis to read counts ---------------------------------------------

assays(sce) <- list(counts = assay(sce, "read_counts"))

# Filter genes and samples -----------------------------------------------------

# Filter to only retain protein coding genes.
sce <- sce[any(grepl("protein_coding", rowData(sce)$ENSEMBL.GENEBIOTYPE)), ]

# Filter out WT samples as these are not useful for DE analysis.
sce <- sce[, sce$smchd1_genotype_updated != "WT"]

# Update metadata
colData(sce) <- droplevels(colData(sce))

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

# Create DGEList ---------------------------------------------------------------

x <- DGEList(
  counts = as.matrix(counts(sce)),
  samples = colData(sce),
  group = sce$smchd1_genotype_updated,
  genes = flattenDF(rowData(sce)))

# Variation of technical replicates --------------------------------------------

par(mfrow = c(2, 2))
plotRLE(
  x$counts,
  outline = FALSE,
  col = x$samples$smchd1_genotype_updated_colours,
  las = 2,
  main = "Raw",
  cex = 0.5)
plotPCA(
  x$counts,
  col = x$samples$smchd1_genotype_updated_colours,
  main = "Raw")
plotPCA(
  x$counts,
  col = x$samples$sex_colours,
  main = "Raw")
plotPCA(
  x$counts,
  col = x$samples$mouse_number_colours,
  main = "Raw")

x <- calcNormFactors(x, method = "upperquartile")
par(mfrow = c(2, 2))
plotRLE(
  cpm(x),
  outline = FALSE,
  col = x$samples$smchd1_genotype_updated_colours,
  las = 2,
  main = "UQ",
  cex = 0.5)
plotPCA(
  cpm(x),
  col = x$samples$smchd1_genotype_updated_colours,
  main = "UQ")
plotPCA(
  cpm(x),
  col = x$samples$sex_colours,
  main = "UQ")
plotPCA(
  cpm(x),
  col = x$samples$mouse_number_colours,
  main = "UQ")

# Summary: Even with UQ normalization, there is considerable intra-sample
# variation between the technical replicates (at least for some samples).

# Aggregate (sum) technical replicates -----------------------------------------

y <- sumTechReps(x, x$samples$genotype.mouse)

# Gene filtering ---------------------------------------------------------------

keep <- filterByExpr(
  y,
  group = y$samples$smchd1_genotype_updated,
  min.count = 5)
y <- y[keep, , keep.lib.sizes = FALSE]

# Variation of biological replicates -------------------------------------------

par(mfrow = c(2, 2))
plotRLE(
  y$counts,
  outline = FALSE,
  col = y$samples$smchd1_genotype_updated_colours,
  las = 2,
  main = "Raw",
  cex = 0.5,
  )
plotPCA(
  y$counts,
  col = y$samples$smchd1_genotype_updated_colours,
  main = "Raw")
plotPCA(
  y$counts,
  col = y$samples$sex_colours,
  main = "Raw")
plotPCA(
  y$counts,
  col = y$samples$mouse_number_colours,
  main = "Raw")

y <- calcNormFactors(y, method = "upperquartile")
par(mfrow = c(2, 2))
plotRLE(
  cpm(y),
  outline = FALSE,
  col = y$samples$smchd1_genotype_updated_colours,
  las = 2,
  main = "UQ",
  cex = 0.5)
plotPCA(
  cpm(y),
  col = y$samples$smchd1_genotype_updated_colours,
  main = "UQ")
plotPCA(
  cpm(y),
  col = y$samples$sex_colours,
  main = "UQ")
plotPCA(
  cpm(y),
  col = y$samples$mouse_number_colours,
  main = "UQ")

# Vanilla DE analysis ----------------------------------------------------------

design <- model.matrix(
  ~0 + smchd1_genotype_updated + sex,
  y$samples)
colnames(design) <- sub(
  "smchd1_genotype_updated|sex",
  "",
  colnames(design))

y <- estimateDisp(y, design)
par(mfrow = c(1, 1))
plotBCV(y)

dgeglm <- glmFit(y, design)
contrasts <- makeContrasts(Het - Del, levels = design)
dgelrt <- glmLRT(dgeglm, contrast = contrasts)
tt <- topTags(dgelrt, n = Inf)
summary(decideTests(dgelrt))

par(mfrow = c(2, 2))
hist(
  tt$table$PValue,
  breaks = 0:20 / 20,
  xlab = "PValue")
plotMD(dgelrt)
abline(h = 0, lty = 2, col = "dodgerBlue")
plotMD(
  dgelrt,
  status = ifelse(rownames(dgelrt) == "Smchd1", "Smchd1", "Other"),
  hl.col = "orange")
abline(h = 0, lty = 2, col = "dodgerBlue")

# RUVs DE analysis -------------------------------------------------------------

differences <- makeGroups(
  interaction(y$samples$smchd1_genotype_updated, y$samples$sex))
ruv <- RUVs(
    x = betweenLaneNormalization(y$counts, "upper"),
    cIdx = rownames(y),
    k = 3,
    scIdx = differences)

design <- model.matrix(
  ~0 + smchd1_genotype_updated + sex + ruv$W,
  y$samples)
colnames(design) <- sub(
  "ruv\\$|smchd1_genotype_updated|sex",
  "",
  colnames(design))

y <- estimateDisp(y, design)
par(mfrow = c(1, 1))
plotBCV(y)

par(mfrow = c(2, 2))
plotRLE(
  ruv$normalizedCounts,
  outline = FALSE,
  col = y$samples$smchd1_genotype_updated_colours,
  las = 2)
plotPCA(
  ruv$normalizedCounts,
  col = y$samples$smchd1_genotype_updated_colours)
plotPCA(
  ruv$normalizedCounts,
  col = y$samples$sex_colours)

dgeglm <- glmFit(y, design)
contrasts <- makeContrasts(Het - Del, levels = design)
dgelrt <- glmLRT(dgeglm, contrast = contrasts)
tt <- topTags(dgelrt, n = Inf)
summary(decideTests(dgelrt))

par(mfrow = c(2, 2))
hist(
  tt$table$PValue,
  breaks = 0:20 / 20,
  xlab = "PValue")
plotMD(dgelrt)
abline(h = 0, lty = 2, col = "dodgerBlue")
plotMD(
  dgelrt,
  status = ifelse(rownames(dgelrt) == "Smchd1", "Smchd1", "Other"),
  hl.col = "orange")
abline(h = 0, lty = 2, col = "dodgerBlue")

# Visualizations of DE results -------------------------------------------------

pheatmap(
  cpm(y, log = TRUE)[head(rownames(tt[tt$table$FDR < 0.05, ]), 100), ],
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = hcl.colors(101, "Blue-Red 3"),
  breaks = seq(-3, 3, length.out = 101),
  scale = "row",
  annotation_col = data.frame(
    genotype = y$samples$smchd1_genotype_updated,
    sex = y$samples$sex,
    row.names = colnames(y)),
  annotation_colors = list(
    genotype = smchd1_genotype_updated_colours,
    sex = sex_colours),
  main = "UQ")

pheatmap(
  ruv$normalizedCounts[head(rownames(tt[tt$table$FDR < 0.05, ]), 100), ],
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = hcl.colors(101, "Blue-Red 3"),
  breaks = seq(-3, 3, length.out = 101),
  scale = "row",
  annotation_col = data.frame(
    genotype = y$samples$smchd1_genotype_updated,
    sex = y$samples$sex,
    row.names = colnames(y)),
  annotation_colors = list(
    genotype = smchd1_genotype_updated_colours,
    sex = sex_colours),
  main = "RUV")

# Gene set analysis ------------------------------------------------------------

# TODO: These need to incorporate directionality.

male_het_vs_del <- readxl::read_excel(
  here("data", "sarahs_gene_lists", "MaleHetvDel_EdgeR.xlsx"))
female_het_vs_del <- readxl::read_excel(
  here("data", "sarahs_gene_lists", "Bcell_female_HetvDel_EdgeR_mRNA.xlsx"))

y_X <- y$genes$ENSEMBL.SEQNAME == "X"
y_protocadherin <- grepl("Pcdh", rownames(y))
y_pw <- rownames(y) %in% c("Ndn", "Mkrn3", "Peg12")
# TODO: Something more sophisticated than taking genes with FDR < 0.05?
y_male <- rownames(y) %in% male_het_vs_del$Probe[male_het_vs_del$FDR < 0.05]
y_female <- rownames(y) %in%
  female_het_vs_del$Probe[female_het_vs_del$FDR < 0.05]

y_index <- list(
  `X-linked` = y_X,
  Protocadherin = y_protocadherin,
  `Prader-Willi` = y_pw,
  `Male: Het vs. Del` = y_male,
  `Female: Het vs. Del` = y_female)

fry(
  y,
  index = y_index,
  design = design,
  contrast = contrasts)

par(mfrow = c(length(y_index), 2))
for (n in names(y_index)) {
  y_status <- rep("Other", nrow(dgelrt))
  y_status[y_index[[n]]] <- n
  plotMD(
    dgelrt,
    status = y_status,
    values = n,
    hl.col = "red",
    legend = "bottomright")
  abline(h = 0, col = "darkgrey")

  barcodeplot(
    dgelrt$table$logFC,
    y_index[[n]],
    labels = c("Down in Del", "Up in Del"))
}

# Vanilla DE analysis (females) ------------------------------------------------

y <- sumTechReps(x, x$samples$genotype.mouse)
y <- y[, y$samples$sex == "F"]
keep <- filterByExpr(
  y,
  group = y$samples$smchd1_genotype_updated,
  min.count = 5)
y <- y[keep, , keep.lib.sizes = FALSE]

design <- model.matrix(
  ~0 + smchd1_genotype_updated,
  y$samples)
colnames(design) <- sub(
  "smchd1_genotype_updated",
  "",
  colnames(design))

y <- estimateDisp(y, design)
par(mfrow = c(1, 1))
plotBCV(y)

dgeglm <- glmFit(y, design)
contrasts <- makeContrasts(Het - Del, levels = design)
dgelrt <- glmLRT(dgeglm, contrast = contrasts)
tt <- topTags(dgelrt, n = Inf)
summary(decideTests(dgelrt))

par(mfrow = c(2, 2))
hist(
  tt$table$PValue,
  breaks = 0:20 / 20,
  xlab = "PValue")
plotMD(dgelrt)
abline(h = 0, lty = 2, col = "dodgerBlue")
plotMD(
  dgelrt,
  status = ifelse(rownames(dgelrt) == "Smchd1", "Smchd1", "Other"),
  hl.col = "orange")
abline(h = 0, lty = 2, col = "dodgerBlue")

# RUVs DE analysis (females) ---------------------------------------------------

differences <- makeGroups(y$samples$smchd1_genotype_updated)
ruv <- RUVs(
  x = betweenLaneNormalization(y$counts, "upper"),
  cIdx = rownames(y),
  k = 1,
  scIdx = differences)

design <- model.matrix(
  ~0 + smchd1_genotype_updated + ruv$W,
  y$samples)
colnames(design) <- sub(
  "ruv\\$|smchd1_genotype_updated",
  "",
  colnames(design))

y <- estimateDisp(y, design)
par(mfrow = c(1, 1))
plotBCV(y)

par(mfrow = c(2, 2))
plotRLE(
  ruv$normalizedCounts,
  outline = FALSE,
  col = y$samples$smchd1_genotype_updated_colours,
  las = 2)
plotPCA(
  ruv$normalizedCounts,
  col = y$samples$smchd1_genotype_updated_colours)
plotPCA(
  ruv$normalizedCounts,
  col = y$samples$sex_colours)

dgeglm <- glmFit(y, design)
contrasts <- makeContrasts(Het - Del, levels = design)
dgelrt <- glmLRT(dgeglm, contrast = contrasts)
tt <- topTags(dgelrt, n = Inf)
summary(decideTests(dgelrt))

par(mfrow = c(2, 2))
hist(
  tt$table$PValue,
  breaks = 0:20 / 20,
  xlab = "PValue")
plotMD(dgelrt)
abline(h = 0, lty = 2, col = "dodgerBlue")
plotMD(
  dgelrt,
  status = ifelse(rownames(dgelrt) == "Smchd1", "Smchd1", "Other"),
  hl.col = "orange")
abline(h = 0, lty = 2, col = "dodgerBlue")

# Vanilla DE analysis (males) ------------------------------------------------

y <- sumTechReps(x, x$samples$genotype.mouse)
y <- y[, y$samples$sex == "M"]
keep <- filterByExpr(
  y,
  group = y$samples$smchd1_genotype_updated,
  min.count = 5)
y <- y[keep, , keep.lib.sizes = FALSE]

design <- model.matrix(
  ~0 + smchd1_genotype_updated,
  y$samples)
colnames(design) <- sub(
  "smchd1_genotype_updated",
  "",
  colnames(design))

y <- estimateDisp(y, design)
par(mfrow = c(1, 1))
plotBCV(y)

dgeglm <- glmFit(y, design)
contrasts <- makeContrasts(Het - Del, levels = design)
dgelrt <- glmLRT(dgeglm, contrast = contrasts)
tt <- topTags(dgelrt, n = Inf)
summary(decideTests(dgelrt))

par(mfrow = c(2, 2))
hist(
  tt$table$PValue,
  breaks = 0:20 / 20,
  xlab = "PValue")
plotMD(dgelrt)
abline(h = 0, lty = 2, col = "dodgerBlue")
plotMD(
  dgelrt,
  status = ifelse(rownames(dgelrt) == "Smchd1", "Smchd1", "Other"),
  hl.col = "orange")
abline(h = 0, lty = 2, col = "dodgerBlue")

# RUVs DE analysis (males) -----------------------------------------------------

differences <- makeGroups(y$samples$smchd1_genotype_updated)
ruv <- RUVs(
  x = betweenLaneNormalization(y$counts, "upper"),
  cIdx = rownames(y),
  k = 1,
  scIdx = differences)

design <- model.matrix(
  ~0 + smchd1_genotype_updated + ruv$W,
  y$samples)
colnames(design) <- sub(
  "ruv\\$|smchd1_genotype_updated",
  "",
  colnames(design))

y <- estimateDisp(y, design)
par(mfrow = c(1, 1))
plotBCV(y)

par(mfrow = c(2, 2))
plotRLE(
  ruv$normalizedCounts,
  outline = FALSE,
  col = y$samples$smchd1_genotype_updated_colours,
  las = 2)
plotPCA(
  ruv$normalizedCounts,
  col = y$samples$smchd1_genotype_updated_colours)
plotPCA(
  ruv$normalizedCounts,
  col = y$samples$sex_colours)

dgeglm <- glmFit(y, design)
contrasts <- makeContrasts(Het - Del, levels = design)
dgelrt <- glmLRT(dgeglm, contrast = contrasts)
tt <- topTags(dgelrt, n = Inf)
summary(decideTests(dgelrt))

par(mfrow = c(2, 2))
hist(
  tt$table$PValue,
  breaks = 0:20 / 20,
  xlab = "PValue")
plotMD(dgelrt)
abline(h = 0, lty = 2, col = "dodgerBlue")
plotMD(
  dgelrt,
  status = ifelse(rownames(dgelrt) == "Smchd1", "Smchd1", "Other"),
  hl.col = "orange")
abline(h = 0, lty = 2, col = "dodgerBlue")

