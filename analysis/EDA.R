# Run after multi-sample comparisons.

library(ggrepel)

# MDS of tech reps -------------------------------------------------------------

mds <- plotMDS(x, plot = FALSE)

df <- cbind(data.frame(x = mds$x, y = mds$y), x$samples)
rownames(df) <- colnames(x)

p1 <- ggplot(
  aes(x = x, y = y, colour = smchd1_genotype_updated, label = rownames(df)),
  data = df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2") +
  scale_colour_manual(values = smchd1_genotype_updated_colours)
p2 <- ggplot(
  aes(x = x, y = y, colour = sex, label = rownames(df)),
  data = df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2") +
  scale_colour_manual(values = sex_colours)
p3 <- ggplot(
  aes(x = x, y = y, colour = mouse_number, label = rownames(df)),
  data = df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2") +
  scale_colour_manual(values = mouse_number_colours)
p4 <- ggplot(
  aes(x = x, y = y, colour = log10(lib.size), label = rownames(df)),
  data = df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  scale_colour_viridis_c() +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2")
p5 <- ggplot(
  aes(x = log10(lib.size), y = x, colour = log10(lib.size), label = rownames(df)),
  data = df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  scale_colour_viridis_c() +
  ylab("Leading logFC dim 1")
p6 <- ggplot(
  aes(x = log10(lib.size), y = y, colour = log10(lib.size), label = rownames(df)),
  data = df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  scale_colour_viridis_c() +
  ylab("Leading logFC dim 2")

p1 + p2 + p3 + p4 + p5 + p6

# MDS of tech reps (excluding outliers) ----------------------------------------

outliers <- c("Del.1282.2", "Het.1248.2")
x2 <- x[, setdiff(colnames(x), outliers)]

mds <- plotMDS(x2, plot = FALSE)

df <- cbind(data.frame(x = mds$x, y = mds$y), x2$samples)
rownames(df) <- colnames(x2)

p1 <- ggplot(
  aes(x = x, y = y, colour = smchd1_genotype_updated, label = rownames(df)),
  data = df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2") +
  scale_colour_manual(values = smchd1_genotype_updated_colours)
p2 <- ggplot(
  aes(x = x, y = y, colour = sex, label = rownames(df)),
  data = df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2") +
  scale_colour_manual(values = sex_colours)
p3 <- ggplot(
  aes(x = x, y = y, colour = mouse_number, label = rownames(df)),
  data = df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2") +
  scale_colour_manual(values = mouse_number_colours)
p4 <- ggplot(
  aes(x = x, y = y, colour = log10(lib.size), label = rownames(df)),
  data = df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  scale_colour_viridis_c() +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2")
p5 <- ggplot(
  aes(x = log10(lib.size), y = x, colour = log10(lib.size), label = rownames(df)),
  data = df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  scale_colour_viridis_c() +
  ylab("Leading logFC dim 1")
p6 <- ggplot(
  aes(x = log10(lib.size), y = y, colour = log10(lib.size), label = rownames(df)),
  data = df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  scale_colour_viridis_c() +
  ylab("Leading logFC dim 2")

p1 + p2 + p3 + p4 + p5 + p6


# MDS of aggregated tech reps --------------------------------------------------

y <- sumTechReps(x, x$samples$genotype.mouse)

mds <- plotMDS(y, plot = FALSE)

df <- cbind(data.frame(x = mds$x, y = mds$y), y$samples)
rownames(df) <- colnames(y)

p1 <- ggplot(
  aes(x = x, y = y, colour = smchd1_genotype_updated, label = rownames(df)),
  data = df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2") +
  scale_colour_manual(values = smchd1_genotype_updated_colours)
p2 <- ggplot(
  aes(x = x, y = y, colour = sex, label = rownames(df)),
  data = df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2") +
  scale_colour_manual(values = sex_colours)
p3 <- ggplot(
  aes(x = x, y = y, colour = mouse_number, label = rownames(df)),
  data = df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2") +
  scale_colour_manual(values = mouse_number_colours)
p4 <- ggplot(
  aes(x = x, y = y, colour = log10(lib.size), label = rownames(df)),
  data = df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  scale_colour_viridis_c() +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2")
p5 <- ggplot(
  aes(x = log10(lib.size), y = x, colour = log10(lib.size), label = rownames(df)),
  data = df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  scale_colour_viridis_c() +
  ylab("Leading logFC dim 1")
p6 <- ggplot(
  aes(x = log10(lib.size), y = y, colour = log10(lib.size), label = rownames(df)),
  data = df) +
  geom_point() +
  geom_text_repel(size = 4) +
  theme_cowplot() +
  scale_colour_viridis_c() +
  ylab("Leading logFC dim 2")

p1 + p2 + p3 + p4 + p5 + p6

# arrayWeights of tech reps ----------------------------------------------------

x_o <- order(colnames(x))
# TODO: Should I be using voomLmFit() and supplying mouse_number as block?
# x_design <- model.matrix(~smchd1_genotype_updated + mouse_number, x$samples)
# x_v <- voom(x, x_design)
# x_aw <- arrayWeights(x_v)
# names(x_aw) <- colnames(x)
x_design <- model.matrix(~0 + smchd1_genotype_updated, x$samples)
colnames(x_design) <- sub("smchd1_genotype_updated", "", colnames(x_design))
# TODO: How to interpret small intra-block correlation (0.02)
# TODO: Interpret voom plot
x_fit <- voomLmFit(x, x_design, block = x$samples$mouse_number, sample.weights = TRUE, plot = TRUE)
x_aw <- x_fit$targets$sample.weight
names(x_aw) <- colnames(x)
x_aw <- x_aw[x_o]
barplot(x_aw, col = x$samples$smchd1_genotype_updated_colours[x_o], las = 2)

# arrayWeights of aggregated tech reps -----------------------------------------

y_o <- order(colnames(y))
y_design <- model.matrix(~0 + smchd1_genotype_updated, y$samples)
colnames(y_design) <- sub("smchd1_genotype_updated", "", colnames(y_design))
y_v <- voom(y, y_design)
y_aw <- arrayWeights(y_v)
names(y_aw) <- colnames(y)
y_aw <- y_aw[y_o]
barplot(y_aw, col = y$samples$smchd1_genotype_updated_colours[y_o], las = 2)

# Comparing arrayWeights of tech reps vs. aggregated tech reps -----------------

par(mfrow = c(2, 1))
barplot(x_aw, col = x$samples$smchd1_genotype_updated_colours[x_o], las = 2)
barplot(y_aw, col = y$samples$smchd1_genotype_updated_colours[y_o], las = 2)

# DE of tech reps --------------------------------------------------------------

x_contrasts <- makeContrasts(Del - Het, levels = x_design)
# TODO: Rename x_v as x_fit (above)
x_fit <- contrasts.fit(x_v, x_contrasts)
x_fit <- eBayes(x_fit)
plotSA(x_fit, main="Final model: Mean-variance trend")
x_dt <- decideTests(x_fit)
summary(x_dt)
par(mfrow = c(1, 1))
plotMD(x_fit, status = x_dt[, 1])
topTable(x_fit)

# DE of aggregated tech reps ---------------------------------------------------

y_contrasts <- makeContrasts(Del - Het, levels = y_design)
y_fit <- lmFit(y_v, y_design)
y_fit <- contrasts.fit(y_fit, y_contrasts)
y_fit <- eBayes(y_fit)
plotSA(y_fit, main="Final model: Mean-variance trend")
y_dt <- decideTests(y_fit)
summary(y_dt)
par(mfrow = c(1, 1))
plotMD(y_fit, status = y_dt[, 1])
topTable(y_fit)

# TODOs ------------------------------------------------------------------------

# UP TO HERE

# - [x] DE analysis using x_v (voomLmFit blocking on mouse_number) and y_v
#   - [ ] Compare gene lists
# - [ ] Remember, if using `y` then should re-filter genes.

# OLD STUFF --------------------------------------------------------------------

par(mfrow = c(1, 3))
plotMDS(y, col = y$samples$sex_colours, main = "Sex", cex = 0.7)
legend("topright", legend = names(sex_colours), col = sex_colours, pch = 16, bty = "n", cex = 0.7)
plotMDS(y, col = y$samples$mouse_number_colours, main = "Mouse", cex = 0.7)
plotMDS(y, col = y$samples$smchd1_genotype_updated_colours, main = "smchd1_genotype_updated", cex = 0.7)
legend("topright", legend = names(smchd1_genotype_updated_colours), col = smchd1_genotype_updated_colours, pch = 16, bty = "n", cex = 0.7)

y2 <- y[setdiff(rownames(y), sex_set), ]
par(mfrow = c(1, 3))
plotMDS(y2, col = y2$samples$sex_colours, main = "Sex", cex = 0.7)
legend("topright", legend = names(sex_colours), col = sex_colours, pch = 16, bty = "n", cex = 0.7)
plotMDS(y2, col = y2$samples$mouse_number_colours, main = "Mouse", cex = 0.7)
plotMDS(y2, col = y2$samples$smchd1_genotype_updated_colours, main = "smchd1_genotype_updated", cex = 0.7)
legend("topright", legend = names(smchd1_genotype_updated_colours), col = smchd1_genotype_updated_colours, pch = 16, bty = "n", cex = 0.7)

z <- y[, !y$samples$mouse_number %in% c("1282", "1248")]
par(mfrow = c(1, 3))
plotMDS(z, col = z$samples$sex_colours, main = "Sex", cex = 0.7)
legend("topright", legend = names(sex_colours), col = sex_colours, pch = 16, bty = "n", cex = 0.7)
plotMDS(z, col = z$samples$mouse_number_colours, main = "Mouse", cex = 0.7)
plotMDS(z, col = z$samples$smchd1_genotype_updated_colours, main = "smchd1_genotype_updated", cex = 0.7)
legend("topright", legend = names(smchd1_genotype_updated_colours), col = smchd1_genotype_updated_colours, pch = 16, bty = "n", cex = 0.7)

z2 <- y2[, !y2$samples$mouse_number %in% c("1282", "1248")]
par(mfrow = c(1, 3))
plotMDS(z2, col = z2$samples$sex_colours, main = "Sex", cex = 0.7)
legend("bottomright", legend = names(sex_colours), col = sex_colours, pch = 16, bty = "n", cex = 0.7)
plotMDS(z2, col = z2$samples$mouse_number_colours, main = "Mouse", cex = 0.7)
plotMDS(z2, col = z2$samples$smchd1_genotype_updated_colours, main = "smchd1_genotype_updated", cex = 0.7)
legend("bottomright", legend = names(smchd1_genotype_updated_colours), col = smchd1_genotype_updated_colours, pch = 16, bty = "n", cex = 0.7)

# DE of summed tech reps -------------------------------------------------------

y <- y[, order(colnames(y))]
y$samples$tmp <- relevel(y$samples$smchd1_genotype_updated, "WT")
design <- model.matrix(~tmp + sex, y$samples)
colnames(design) <- sub("tmp|sex", "", colnames(design))

v <- voomWithQualityWeights(y, design = design, plot = TRUE, col = y$samples$smchd1_genotype_updated_colours)
colnames(v)

con <- limma::makeContrasts(Del - Het, levels = design)
fit <- lmFit(v, design)
fit <- contrasts.fit(fit, con)
fit <- eBayes(fit)
plotSA(fit, main="Final model: Mean-variance trend")
dt <- decideTests(fit)
summary(dt)
par(mfrow = c(1, 1))
plotMD(fit, status = dt[, 1])
topTable(fit)

# DE of summed tech reps (excl. WT) --------------------------------------------

y3 <- y[, y$samples$smchd1_genotype_updated != "WT"]
y3$samples <- droplevels(y3$samples)
design <- model.matrix(~0 + tmp, y3$samples)
colnames(design) <- sub("tmp|sex", "", colnames(design))
v <- voomWithQualityWeights(y3, design = design, plot = TRUE, col = y$samples$smchd1_genotype_updated_colours)
colnames(v)

con <- limma::makeContrasts(Del - Het, levels = design)
fit <- lmFit(v, design)
fit <- contrasts.fit(fit, con)
fit <- eBayes(fit)
plotSA(fit, main="Final model: Mean-variance trend")
dt <- decideTests(fit)
summary(dt)
par(mfrow = c(1, 1))
plotMD(fit, status = dt[, 1])
topTable(fit)

# DE of tech reps --------------------------------------------------------------

x$samples$tmp <- relevel(x$samples$smchd1_genotype_updated, "WT")
design <- model.matrix(~tmp + sex, x$samples)
colnames(design) <- sub("tmp|sex", "", colnames(design))
fit <- voomLmFit(x, design, block = x$samples$mouse_number, plot = TRUE)

con <- limma::makeContrasts(Del - Het, levels = design)
fit <- contrasts.fit(fit, con)
fit <- eBayes(fit)
plotSA(fit, main="Final model: Mean-variance trend")
dt <- decideTests(fit)
summary(dt)
par(mfrow = c(1, 1))
plotMD(fit, status = dt[, 1])
topTable(fit)
