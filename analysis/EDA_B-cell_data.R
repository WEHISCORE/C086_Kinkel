# DE analysis of Sarah's B cell RNA-seq data
# Peter Hickey
# 2021-10-06

# Setup ------------------------------------------------------------------------

library(here)
library(edgeR)

x <- data.table::fread(
  here("data/sarahs_gene_lists/Smchd1_femBcells_rawcounts.txt"),
  data.table = FALSE)
# TODO: What are the unique gene IDs? There are genes with duplicated `Probe`
#       or `Feature` or `ID`.
# This could explain how we end up with probes with the same name because:
#   Almost all rows in the count matrix have a unique ID​ value and Ensembl gene IDs are unambiguous and highly stable
# The exception are 56 features where the ID​ is blank but there is a I'm not sure what happened there.
# There are probes with the same name because multiple Ensembl genes can map to the same gene symbol.
# An example is Nav1​ which has 2 entries in the count matrix.
# We need to filter/aggregate these 2 entries prior to analysis
# > x[x$Feature == "Nav1", c("Probe", "ID", grep("bam", colnames(x), value = TRUE))]
#     Probe                 ID Bcell_B167het_trim.bam Bcell_B189het_trim.bam
# 847  Nav1 ENSMUSG00000009418                     36                     30
# 849  Nav1 ENSMUSG00000090399                      2                      0
#
# ENSMUSG00000009418 is 'active' whereas ENSMUSG00000090399 has been 'archived'.

y <- DGEList(
  counts = as.matrix(x[, grep("bam", colnames(x))]),
  genes = x[
    ,
    c("Probe", "Chromosome", "Start", "End", "Probe Strand", "Feature", "ID",
      "Description", "Feature Strand", "Type", "Feature Orientation")])
y$samples$genotype <- ifelse(grepl("het", colnames(y)), "Het", "Del")
y$samples$mouse <- substr(colnames(y), 7, 10)
y$samples$group <- y$samples$genotype
library(Mus.musculus)
y$genes$NCBI.ENTREZID <- unstrsplit(
  mapIds(
    x = Mus.musculus,
    # NOTE: Need to remove gene version number prior to lookup.
    keys = y$genes$ID,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "CharacterList"),
  sep = "; ")

# Gene filtering ---------------------------------------------------------------

y_ <- y
y_keep_exprs <- filterByExpr(y, group = y$samples$group, min.count = 5)
y <- y[y_keep_exprs, , keep.lib.sizes = FALSE]
table(y_keep_exprs)

# Normalization ----------------------------------------------------------------

y <- calcNormFactors(y, method = "TMM")

# EDA --------------------------------------------------------------------------

plotMDS(y, col = as.integer(factor(y$samples$genotype)))

# DE analysis ------------------------------------------------------------------

y_design <- model.matrix(
  ~0 + genotype,
  y$samples)
colnames(y_design) <- sub("genotype", "", colnames(y_design))
y <- estimateDisp(y, y_design)
plotBCV(y)

y_dgeglm <- glmFit(y, y_design)
y_contrast <- makeContrasts(Del - Het, levels = y_design)
y_dgelrt <- glmLRT(y_dgeglm, contrast = y_contrast)
y_dt <- decideTests(y_dgelrt)

summary(y_dt)

plotMD(y_dgelrt)

# Gene set tests ---------------------------------------------------------------

# goana
list_of_go <- lapply(colnames(y_contrast), function(j) {
  dgelrt <- glmLRT(y_dgeglm, contrast = y_contrast[, j])

  go <- goana(
    de = dgelrt,
    # TODO: Not sure why, but need to use `[` rather than `[[`
    geneid = sapply(strsplit(dgelrt$genes$NCBI.ENTREZID, "; "), "[", 1),
    species = "Mm")
  # gzout <- gzfile(
  #   description = file.path(outdir, "aggregated_tech_reps.goana.csv.gz"),
  #   open = "wb")
  # write.csv(
  #   topGO(go, number = Inf),
  #   gzout,
  #   row.names = TRUE)
  # close(gzout)
  go
})
names(list_of_go) <- colnames(y_contrast)
head(topGO(list_of_go[[1]]))

# kegga
list_of_kegg <- lapply(colnames(y_contrast), function(j) {
  dgelrt <- glmLRT(y_dgeglm, contrast = y_contrast[, j])

  kegg <- kegga(
    de = dgelrt,
    # TODO: Not sure why, but need to use `[` rather than `[[`
    geneid = sapply(strsplit(dgelrt$genes$NCBI.ENTREZID, "; "), "[", 1),
    species = "Mm")
  # gzout <- gzfile(
  #   description = file.path(outdir, "aggregated_tech_reps.kegga.csv.gz"),
  #   open = "wb")
  # write.csv(
  #   topKEGG(kegg, number = Inf),
  #   gzout,
  #   row.names = TRUE)
  # close(gzout)
  kegg
})
names(list_of_kegg) <- colnames(y_contrast)
head(topKEGG(list_of_kegg[[1]]))

# camera
# NOTE: Using BiocFileCache to avoid re-downloading these gene sets everytime
#       the report is rendered.
library(BiocFileCache)
bfc <- BiocFileCache()
# NOTE: Creating list of gene sets in this slightly convoluted way so that each
#       gene set name is prepended by its origin (e.g. H, C2, or C7).
msigdb <- do.call(
  c,
  list(
    H = readRDS(
      bfcrpath(
        bfc,
        "http://bioinf.wehi.edu.au/MSigDB/v7.1/Mm.h.all.v7.1.entrez.rds")),
    C2 = readRDS(
      bfcrpath(
        bfc,
        "http://bioinf.wehi.edu.au/MSigDB/v7.1/Mm.c2.all.v7.1.entrez.rds"))))

y_idx <- ids2indices(
  msigdb,
  # TODO: Not sure why, but need to use `[` rather than `[[`
  id = sapply(strsplit(y$genes$NCBI.ENTREZID, "; "), "[", 1))
list_of_camera <- lapply(colnames(y_contrast), function(j) {
  cam <- camera(
    y = y,
    index = y_idx,
    design = y_design,
    contrast = y_contrast[, j])
  # gzout <- gzfile(
  #   description = file.path(outdir, "aggregated_tech_reps.camera.csv.gz"),
  #   open = "wb")
  # write.csv(
  #   cam,
  #   gzout,
  #   row.names = TRUE)
  # close(gzout)
  cam
})
names(list_of_camera) <- colnames(y_contrast)
head(list_of_camera[[1]])
head(list_of_camera[[1]][grep("^H", rownames(list_of_camera[[1]])), ])

# chrX gene set test
y_X <- y$genes$Chromosome == "X"
y_index <- list(`X-linked` = y_X)

fry(
  y,
  index = y_index,
  design = y_design,
  contrast = y_contrast,
  sort = "none")

par(mfrow = c(1, 2))
n <- names(y_index)[1]
y_status <- rep("Other", nrow(y_dgelrt))
y_status[y_index[[n]]] <- n
plotMD(
  y_dgelrt,
  status = y_status,
  values = n,
  hl.col = "red",
  legend = "topright",
  main = n)
abline(h = 0, col = "darkgrey")
barcodeplot(
  y_dgelrt$table$logFC,  y_index[[n]])

# Comparison to results from Sarah's analysis ----------------------------------

# TODO: The DE results are slightly different to what Sarah got.

female_het_vs_del <- readxl::read_excel(
  here("data", "sarahs_gene_lists", "Bcell_female_HetvDel_EdgeR_mRNA.xlsx"))
female_het_vs_del <- as.data.frame(female_het_vs_del)

y_tt <- topTags(y_dgelrt, n = Inf)
# DEGs only found by Sarah
setdiff(female_het_vs_del$Probe, y_tt$table$Probe[y_tt$table$FDR < 0.05])

i <- intersect(
  female_het_vs_del$Probe,
  y_tt$table$Probe[y_tt$table$FDR < 0.05])
tmp <- data.frame(
  probe = i,
  sarah = female_het_vs_del[
    match(i, female_het_vs_del$Probe), "Differential expression"],
  me = y_tt[match(i, y_tt$table$Probe), ]$table$logFC)
par(mfrow = c(1, 1))
plot(tmp$sarah, tmp$me, xlab = "Sarah logFC", ylab = "Pete logFC")
abline(a = 0, b = 1, col = "red")

# Results for Sarah's DEGs that aren't DE (but are tested) in my analysis
y_tt[
  y_tt$table$Probe %in%
    intersect(
      female_het_vs_del$Probe,
      y_tt$table$Probe[y_tt$table$FDR >= 0.05]), ]

# DEGs only found by me.
setdiff(y_tt$table$Probe[y_tt$table$FDR < 0.05], female_het_vs_del$Probe)

a <- data.table::fread(
  "data/sarahs_gene_lists/210923_EdgeR stats p1.0 log2.txt",
  data.table = FALSE)
# TODO: There are ~8000 genes with logFC == 0 in Sarah's results. This is rather
#       suspicious and suggests that genes have not been filtered (e.g.,
#       genes with all zero counts would result in logFC = 0).
aa <- a[a$`Log2 Fold Change (EdgeR stats p<1.0 after correction)` != 0, ]
table(
  rowSums(y_[y_$genes$Probe %in% a[a$`Log2 Fold Change (EdgeR stats p<1.0 after correction)` == 0, "Probe"], ]$counts))

par(mfrow = c(1, 2))
barcodeplot(
  aa$`Log2 Fold Change (EdgeR stats p<1.0 after correction)`,
  aa$Chromosome == "X",
  sub = "Sarah's analysis",
  main = n)
barcodeplot(
  y_dgelrt$table$logFC,
  y_index[[n]],
  sub = "My analysis",
  main = n)

# TODO: SeqMonk uses exactTest() and no gene filtering; https://github.com/s-andrews/SeqMonk/blob/master/uk/ac/babraham/SeqMonk/Filters/EdgeRFilter/edger_template.r
