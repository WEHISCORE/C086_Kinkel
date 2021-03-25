# Setup ------------------------------------------------------------------------

library(SingleCellExperiment)
library(here)
library(edgeR)
library(ggplot2)
library(cowplot)
library(patchwork)
library(ggrepel)

source(here("code/helper_functions.R"))

sce <- readRDS(
  here("data", "SCEs", "C086_Kinkel.preprocessed.SCE.rds"))

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

x <- DGEList(
  counts = as.matrix(assay(sce, "read_counts")),
  samples = colData(sce),
  group = factor(sce$smchd1_genotype_updated),
  genes = flattenDF(rowData(sce)))
x <- x[!matrixStats::rowAlls(x$counts == 0), ]
x0 <- x
y <- sumTechReps(x0, x0$samples$genotype.mouse)
ym <- y[, y$samples$sex == "M"]
yf <- y[, y$samples$sex == "F"]

# TODO: Experiment with minimum count
x_keep_exprs <- filterByExpr(x, group = x$samples$group)
table(x_keep_exprs)
x <- x[x_keep_exprs, , keep.lib.sizes = FALSE]
y_keep_exprs <- filterByExpr(y, group = y$samples$group)
table(y_keep_exprs)
y <- y[y_keep_exprs, , keep.lib.sizes = FALSE]
ym_keep_exprs <- filterByExpr(ym, group = ym$samples$group)
table(ym_keep_exprs)
ym <- ym[ym_keep_exprs, , keep.lib.sizes = FALSE]
yf_keep_exprs <- filterByExpr(yf, group = yf$samples$group)
table(yf_keep_exprs)
yf <- yf[yf_keep_exprs, , keep.lib.sizes = FALSE]

x <- calcNormFactors(x, method = "TMM")
y <- calcNormFactors(y, method = "TMM")
ym <- calcNormFactors(ym, method = "TMM")
yf <- calcNormFactors(yf, method = "TMM")

# MDS of tech reps -------------------------------------------------------------

x_mds <- plotMDS(x, plot = FALSE)

x_df <- cbind(data.frame(x = x_mds$x, y = x_mds$y), x$samples)
rownames(x_df) <- colnames(x)

p1 <- ggplot(
  aes(x = x, y = y, colour = smchd1_genotype_updated, label = rownames(x_df)),
  data = x_df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2") +
  scale_colour_manual(values = smchd1_genotype_updated_colours)
p2 <- ggplot(
  aes(x = x, y = y, colour = sex, label = rownames(x_df)),
  data = x_df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2") +
  scale_colour_manual(values = sex_colours)
p3 <- ggplot(
  aes(x = x, y = y, colour = mouse_number, label = rownames(x_df)),
  data = x_df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2") +
  scale_colour_manual(values = mouse_number_colours)
p4 <- ggplot(
  aes(x = x, y = y, colour = log10(lib.size), label = rownames(x_df)),
  data = x_df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  scale_colour_viridis_c() +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2")
p5 <- ggplot(
  aes(x = log10(lib.size), y = x, colour = log10(lib.size), label = rownames(x_df)),
  data = x_df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  scale_colour_viridis_c() +
  ylab("Leading logFC dim 1")
p6 <- ggplot(
  aes(x = log10(lib.size), y = y, colour = log10(lib.size), label = rownames(x_df)),
  data = x_df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  scale_colour_viridis_c() +
  ylab("Leading logFC dim 2")

p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 2)

# MDS of tech reps (excluding outliers) ----------------------------------------

outliers <- c("Del.1282.2", "Het.1248.2")
x2 <- x[, setdiff(colnames(x), outliers)]

x2_mds <- plotMDS(x2, plot = FALSE)

x2_df <- cbind(data.frame(x = x2_mds$x, y = x2_mds$y), x2$samples)
rownames(x2_df) <- colnames(x2)

p1 <- ggplot(
  aes(x = x, y = y, colour = smchd1_genotype_updated, label = rownames(x2_df)),
  data = x2_df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2") +
  scale_colour_manual(values = smchd1_genotype_updated_colours)
p2 <- ggplot(
  aes(x = x, y = y, colour = sex, label = rownames(x2_df)),
  data = x2_df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2") +
  scale_colour_manual(values = sex_colours)
p3 <- ggplot(
  aes(x = x, y = y, colour = mouse_number, label = rownames(x2_df)),
  data = x2_df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2") +
  scale_colour_manual(values = mouse_number_colours)
p4 <- ggplot(
  aes(x = x, y = y, colour = log10(lib.size), label = rownames(x2_df)),
  data = x2_df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  scale_colour_viridis_c() +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2")
p5 <- ggplot(
  aes(x = log10(lib.size), y = x, colour = log10(lib.size), label = rownames(x2_df)),
  data = x2_df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  scale_colour_viridis_c() +
  ylab("Leading logFC dim 1")
p6 <- ggplot(
  aes(x = log10(lib.size), y = y, colour = log10(lib.size), label = rownames(x2_df)),
  data = x2_df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  scale_colour_viridis_c() +
  ylab("Leading logFC dim 2")

p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 2)

# MDS of aggregated tech reps --------------------------------------------------

y_mds <- plotMDS(y, plot = FALSE)

y_df <- cbind(data.frame(x = y_mds$x, y = y_mds$y), y$samples)
rownames(y_df) <- colnames(y)

p1 <- ggplot(
  aes(x = x, y = y, colour = smchd1_genotype_updated, label = rownames(y_df)),
  data = y_df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2") +
  scale_colour_manual(values = smchd1_genotype_updated_colours)
p2 <- ggplot(
  aes(x = x, y = y, colour = sex, label = rownames(y_df)),
  data = y_df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2") +
  scale_colour_manual(values = sex_colours)
p3 <- ggplot(
  aes(x = x, y = y, colour = mouse_number, label = rownames(y_df)),
  data = y_df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2") +
  scale_colour_manual(values = mouse_number_colours)
p4 <- ggplot(
  aes(x = x, y = y, colour = log10(lib.size), label = rownames(y_df)),
  data = y_df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  scale_colour_viridis_c() +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2")
p5 <- ggplot(
  aes(x = log10(lib.size), y = x, colour = log10(lib.size), label = rownames(y_df)),
  data = y_df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  scale_colour_viridis_c() +
  ylab("Leading logFC dim 1")
p6 <- ggplot(
  aes(x = log10(lib.size), y = y, colour = log10(lib.size), label = rownames(y_df)),
  data = y_df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  scale_colour_viridis_c() +
  ylab("Leading logFC dim 2")

p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 2)

# arrayWeights of tech reps ----------------------------------------------------

x_o <- order(colnames(x))
x_design <- model.matrix(~0 + smchd1_genotype_updated, x$samples)
colnames(x_design) <- sub("smchd1_genotype_updated", "", colnames(x_design))
# TODO: How to interpret small intra-block correlation (0.02)
# TODO: Interpret voom plot
x_fit <- voomLmFit(
  x,
  design = x_design,
  block = x$samples$mouse_number,
  sample.weights = TRUE,
  plot = TRUE)
x_aw <- x_fit$targets$sample.weight
names(x_aw) <- colnames(x)
x_aw <- x_aw[x_o]
barplot(x_aw, col = x$samples$smchd1_genotype_updated_colours[x_o], las = 2)

# arrayWeights of aggregated tech reps -----------------------------------------

y_o <- order(colnames(y))
y_design <- model.matrix(~0 + smchd1_genotype_updated, y$samples)
colnames(y_design) <- sub("smchd1_genotype_updated", "", colnames(y_design))
y_fit <- voomLmFit(
  y,
  design = y_design,
  sample.weights = TRUE,
  plot = TRUE)
y_aw <- y_fit$targets$sample.weight
names(y_aw) <- colnames(y)
y_aw <- y_aw[y_o]
barplot(y_aw, col = y$samples$smchd1_genotype_updated_colours[y_o], las = 2)

# Comparing arrayWeights of tech reps vs. aggregated tech reps -----------------

par(mfrow = c(2, 1))
barplot(x_aw, col = x$samples$smchd1_genotype_updated_colours[x_o], las = 2)
barplot(y_aw, col = y$samples$smchd1_genotype_updated_colours[y_o], las = 2)

# DE of tech reps --------------------------------------------------------------

x_contrasts <- makeContrasts(Del - Het, levels = x_design)
x_fit <- contrasts.fit(x_fit, x_contrasts)
x_fit <- eBayes(x_fit)

par(mfrow = c(1, 1))
plotSA(x_fit, main="Final model: Mean-variance trend")

x_dt <- decideTests(x_fit)
summary(x_dt)
par(mfrow = c(1, 1))
plotMD(x_fit, status = x_dt[, 1])

topTable(x_fit)
topTable(x_fit, n = Inf)["Smchd1", ]

# DE of aggregated tech reps ---------------------------------------------------

y_contrasts <- makeContrasts(Del - Het, levels = y_design)
y_fit <- contrasts.fit(y_fit, y_contrasts)
y_fit <- eBayes(y_fit)

par(mfrow = c(1, 1))
plotSA(y_fit, main="Final model: Mean-variance trend")

y_dt <- decideTests(y_fit)
summary(y_dt)
par(mfrow = c(1, 1))
plotMD(y_fit, status = y_dt[, 1])

topTable(y_fit)
topTable(y_fit, n = Inf)["Smchd1", ]

par(mfrow = c(1, 2))
hist(topTable(y_fit, n = Inf)$P.Value, xlim = c(0, 1), breaks = 0:20 / 20)
hist(topTable(y_fit, n = Inf)$adj.P.Val, xlim = c(0, 1), breaks = 0:20 / 20)

y_Y <- y$genes$ENSEMBL.SEQNAME == "Y"
y_X <- y$genes$ENSEMBL.SEQNAME == "X"

y_index <- list(Y = y_Y, X = y_X)
y <- estimateDisp(y, design = y_design)
fry(y, index = y_index, design = y_design, contrast = y_contrasts)

v_y <- voomWithQualityWeights(y, y_design)
fry(v_y, index = y_index, design = y_design, contrast = y_contrasts)

y_status <- rep("Other", nrow(y_fit))
y_status[y_Y] <- "Y"
y_status[y_X] <- "X"

par(mfrow = c(1, 2))
plotMD(
  y_fit,
  status = y_status,
  values = c("X","Y"),
  hl.col = c("red", "blue"),
  legend = "bottomright")
abline(h = 0, col = "darkgrey")

barcodeplot(
  y_fit$t[, 1],
  y_X,
  y_Y,
  labels = c("Down in Del", "Up in Del"))
legend("top", legend = c("Y", "X"), lty = 1, col = c("blue", "red"))

# MDS, arrayWeights, and DE analysis of male aggregated tech reps --------------

ym_mds <- plotMDS(ym, plot = FALSE)

ym_df <- cbind(data.frame(x = ym_mds$x, y = ym_mds$y), ym$samples)
rownames(ym_df) <- colnames(ym)

p1 <- ggplot(
  aes(x = x, y = y, colour = smchd1_genotype_updated, label = rownames(ym_df)),
  data = ym_df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2") +
  scale_colour_manual(values = smchd1_genotype_updated_colours)
p2 <- ggplot(
  aes(x = x, y = y, colour = sex, label = rownames(ym_df)),
  data = ym_df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2") +
  scale_colour_manual(values = sex_colours)
p3 <- ggplot(
  aes(x = x, y = y, colour = mouse_number, label = rownames(ym_df)),
  data = ym_df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2") +
  scale_colour_manual(values = mouse_number_colours)
p4 <- ggplot(
  aes(x = x, y = y, colour = log10(lib.size), label = rownames(ym_df)),
  data = ym_df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  scale_colour_viridis_c() +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2")
p5 <- ggplot(
  aes(x = log10(lib.size), y = x, colour = log10(lib.size), label = rownames(ym_df)),
  data = ym_df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  scale_colour_viridis_c() +
  ylab("Leading logFC dim 1")
p6 <- ggplot(
  aes(x = log10(lib.size), y = y, colour = log10(lib.size), label = rownames(ym_df)),
  data = ym_df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  scale_colour_viridis_c() +
  ylab("Leading logFC dim 2")

p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 2)

ym_o <- order(colnames(ym))
ym_design <- model.matrix(~0 + smchd1_genotype_updated, ym$samples)
colnames(ym_design) <- sub("smchd1_genotype_updated", "", colnames(ym_design))
ym_fit <- voomLmFit(
  ym,
  design = ym_design,
  sample.weights = TRUE,
  plot = TRUE)
ym_aw <- ym_fit$targets$sample.weight
names(ym_aw) <- colnames(ym)
ym_aw <- ym_aw[ym_o]
barplot(ym_aw, col = ym$samples$smchd1_genotype_updated_colours[ym_o], las = 2)

ym_contrasts <- makeContrasts(Del - Het, levels = ym_design)
ym_fit <- contrasts.fit(ym_fit, ym_contrasts)
ym_fit <- eBayes(ym_fit)

par(mfrow = c(1, 1))
plotSA(ym_fit, main="Final model: Mean-variance trend")

ym_dt <- decideTests(ym_fit)
summary(ym_dt)
par(mfrow = c(1, 1))
plotMD(ym_fit, status = ym_dt[, 1])

topTable(ym_fit)
topTable(ym_fit, n = Inf)["Smchd1", ]

par(mfrow = c(1, 2))
hist(topTable(ym_fit, n = Inf)$P.Value, xlim = c(0, 1), breaks = 0:20 / 20)
hist(topTable(ym_fit, n = Inf)$adj.P.Val, xlim = c(0, 1), breaks = 0:20 / 20)

ym_Y <- ym$genes$ENSEMBL.SEQNAME == "Y"
ym_X <- ym$genes$ENSEMBL.SEQNAME == "X"

ym_index <- list(Y = ym_Y, X = ym_X)
ym <- estimateDisp(ym, design = ym_design)
fry(ym, index = ym_index, design = ym_design, contrast = ym_contrasts)

v_ym <- voomWithQualityWeights(ym, ym_design)
fry(v_ym, index = ym_index, design = ym_design, contrast = ym_contrasts)

ym_status <- rep("Other", nrow(ym_fit))
ym_status[ym_Y] <- "Y"
ym_status[ym_X] <- "X"

par(mfrow = c(1, 2))
plotMD(
  ym_fit,
  status = ym_status,
  values = c("X","Y"),
  hl.col = c("red", "blue"),
  legend = "bottomright")
abline(h = 0, col = "darkgrey")

barcodeplot(
  ym_fit$t[, 1],
  ym_X,
  ym_Y,
  labels = c("Down in Del", "Up in Del"))
legend("top", legend = c("Y", "X"), lty = 1, col = c("blue", "red"))

# MDS, arrayWeights, and DE analysis of female aggregated tech reps ------------

yf_mds <- plotMDS(yf, plot = FALSE)

yf_df <- cbind(data.frame(x = yf_mds$x, y = yf_mds$y), yf$samples)
rownames(yf_df) <- colnames(yf)

p1 <- ggplot(
  aes(x = x, y = y, colour = smchd1_genotype_updated, label = rownames(yf_df)),
  data = yf_df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2") +
  scale_colour_manual(values = smchd1_genotype_updated_colours)
p2 <- ggplot(
  aes(x = x, y = y, colour = sex, label = rownames(yf_df)),
  data = yf_df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2") +
  scale_colour_manual(values = sex_colours)
p3 <- ggplot(
  aes(x = x, y = y, colour = mouse_number, label = rownames(yf_df)),
  data = yf_df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2") +
  scale_colour_manual(values = mouse_number_colours)
p4 <- ggplot(
  aes(x = x, y = y, colour = log10(lib.size), label = rownames(yf_df)),
  data = yf_df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  scale_colour_viridis_c() +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2")
p5 <- ggplot(
  aes(x = log10(lib.size), y = x, colour = log10(lib.size), label = rownames(yf_df)),
  data = yf_df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  scale_colour_viridis_c() +
  ylab("Leading logFC dim 1")
p6 <- ggplot(
  aes(x = log10(lib.size), y = y, colour = log10(lib.size), label = rownames(yf_df)),
  data = yf_df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  scale_colour_viridis_c() +
  ylab("Leading logFC dim 2")

p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 2)

yf_o <- order(colnames(yf))
yf_design <- model.matrix(~0 + smchd1_genotype_updated, yf$samples)
# NOTE: Drop WT column as there are no female WT samples.
yf_design <- yf_design[, colSums(yf_design) > 0]
colnames(yf_design) <- sub("smchd1_genotype_updated", "", colnames(yf_design))
yf_fit <- voomLmFit(
  yf,
  design = yf_design,
  sample.weights = TRUE,
  plot = TRUE)
yf_aw <- yf_fit$targets$sample.weight
names(yf_aw) <- colnames(yf)
yf_aw <- yf_aw[yf_o]
barplot(yf_aw, col = yf$samples$smchd1_genotype_updated_colours[yf_o], las = 2)

yf_contrasts <- makeContrasts(Del - Het, levels = yf_design)
yf_fit <- contrasts.fit(yf_fit, yf_contrasts)
yf_fit <- eBayes(yf_fit)

par(mfrow = c(1, 1))
plotSA(yf_fit, main="Final model: Mean-variance trend")

yf_dt <- decideTests(yf_fit)
summary(yf_dt)
par(mfrow = c(1, 1))
plotMD(yf_fit, status = yf_dt[, 1])

topTable(yf_fit)
topTable(yf_fit, n = Inf)["Smchd1", ]

par(mfrow = c(1, 2))
hist(topTable(yf_fit, n = Inf)$P.Value, xlim = c(0, 1), breaks = 0:20 / 20)
hist(topTable(yf_fit, n = Inf)$adj.P.Val, xlim = c(0, 1), breaks = 0:20 / 20)

yf_Y <- yf$genes$ENSEMBL.SEQNAME == "Y"
yf_X <- yf$genes$ENSEMBL.SEQNAME == "X"

yf_index <- list(Y = yf_Y, X = yf_X)
yf <- estimateDisp(yf, design = yf_design)
fry(yf, index = yf_index, design = yf_design, contrast = yf_contrasts)

v_yf <- voomWithQualityWeights(yf, yf_design)
fry(v_yf, index = yf_index, design = yf_design, contrast = yf_contrasts)

yf_status <- rep("Other", nrow(yf_fit))
yf_status[yf_Y] <- "Y"
yf_status[yf_X] <- "X"

par(mfrow = c(1, 2))
plotMD(
  yf_fit,
  status = yf_status,
  values = c("X","Y"),
  hl.col = c("red", "blue"),
  legend = "bottomright")
abline(h = 0, col = "darkgrey")

barcodeplot(
  yf_fit$t[, 1],
  yf_X,
  yf_Y,
  labels = c("Down in Del", "Up in Del"))
legend("top", legend = c("Y", "X"), lty = 1, col = c("blue", "red"))

# TODOs ------------------------------------------------------------------------

# UP TO HERE

# - [x] DE analysis using x_v (voomLmFit blocking on mouse_number) and y_v
#   - [ ] Compare gene lists
# - [x] Remember, if using `y` then should re-filter genes.
# - [ ] Exclude pseudogenes and all that crap?
