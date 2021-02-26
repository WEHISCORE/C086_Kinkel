# Run after multi-sample comparisons.

par(mfrow = c(1, 3))
plotMDS(x, col = x$samples$sex_colours, main = "Sex", cex = 0.7)
legend("topright", legend = names(sex_colours), col = sex_colours, pch = 16, bty = "n", cex = 0.7)
plotMDS(x, col = x$samples$mouse_number_colours, main = "Mouse", cex = 0.7)
plotMDS(x, col = x$samples$smchd1_genotype_updated_colours, main = "smchd1_genotype_updated", cex = 0.7)
legend("topright", legend = names(smchd1_genotype_updated_colours), col = smchd1_genotype_updated_colours, pch = 16, bty = "n", cex = 0.7)

y <- sumTechReps(x, x$samples$genotype.mouse)

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
