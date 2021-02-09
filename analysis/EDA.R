# EDA of C086_Kinkel
# 2021-02-09
# Peter Hickey

library(SingleCellExperiment)
library(here)

# Using the non-deduped data
sce <- readRDS(here("data/SCEs/C086_Kinkel.not_UMI_deduped.scPipe.SCE.rds"))

# Add gene metadata
# Extract rownames (Ensembl IDs) to use as key in database lookups.
ensembl <- rownames(sce)
# Pull out useful gene-based annotations from the Ensembl-based database.
library(EnsDb.Mmusculus.v79)
library(ensembldb)
# NOTE: These columns were customised for this project.
ensdb_columns <- c(
  "GENEBIOTYPE", "GENENAME", "GENESEQSTART", "GENESEQEND", "SEQNAME", "SYMBOL")
names(ensdb_columns) <- paste0("ENSEMBL.", ensdb_columns)
stopifnot(all(ensdb_columns %in% columns(EnsDb.Mmusculus.v79)))
ensdb_df <- DataFrame(
  lapply(ensdb_columns, function(column) {
    mapIds(
      x = EnsDb.Mmusculus.v79,
      # NOTE: Need to remove gene version number prior to lookup.
      keys = gsub("\\.[0-9]+$", "", ensembl),
      keytype = "GENEID",
      column = column,
      multiVals = "CharacterList")
  }),
  row.names = ensembl)
# NOTE: Can't look up GENEID column with GENEID key, so have to add manually.
ensdb_df$ENSEMBL.GENEID <- ensembl
# NOTE: Mus.musculus combines org.Mm.eg.db and
#       TxDb.Mmusculus.UCSC.mm10.knownGene (as well as others) and therefore
#       uses entrez gene and RefSeq based data.
library(Mus.musculus)
# NOTE: These columns were customised for this project.
ncbi_columns <- c(
  # From TxDB: None required
  # From OrgDB
  "ALIAS", "ENTREZID", "GENENAME", "REFSEQ", "SYMBOL")
names(ncbi_columns) <- paste0("NCBI.", ncbi_columns)
stopifnot(all(ncbi_columns %in% columns(Mus.musculus)))
ncbi_df <- DataFrame(
  lapply(ncbi_columns, function(column) {
    mapIds(
      x = Mus.musculus,
      # NOTE: Need to remove gene version number prior to lookup.
      keys = gsub("\\.[0-9]+$", "", ensembl),
      keytype = "ENSEMBL",
      column = column,
      multiVals = "CharacterList")
  }),
  row.names = ensembl)
rowData(sce) <- cbind(ensdb_df, ncbi_df)
# Replace the row names of the SCE by the gene symbols (where available).
rownames(sce) <- uniquifyFeatureNames(
  ID = rownames(sce),
  # NOTE: An Ensembl ID may map to 0, 1, 2, 3, ... gene symbols.
  #       When there are multiple matches only the 1st match is used.
  names = sapply(rowData(sce)$ENSEMBL.SYMBOL, function(x) {
    if (length(x)) {
      x[[1]]
    } else {
      NA_character_
    }
  }))

# Add sample metadata
library(readxl)
sample_sheet <- read_excel(
  here("data/sample_sheets/C086_Kinkel_slamLSK_MB_NNXXX_Seqprimer11Dec2020.xlsx"),
  sheet = "Pool info",
  range = "A21:G33")
library(janitor)
sample_sheet <- clean_names(sample_sheet)
sample_sheet$mouse_number <- as.character(sample_sheet$mouse_number)
library(dplyr)
colData(sce) <- left_join(
  as.data.frame(colData(sce)),
  sample_sheet,
  by = c("mouse_number" = "mouse_number")) %>%
  DataFrame()
sce$smchd1_genotype_3[is.na(sce$smchd1_genotype_3)] <- "NA"
sce$genotype.mouse <- paste0(sce$smchd1_genotype_3, ".", sce$mouse_number)

count(as.data.frame(colData(sce)), mouse_number, sex, smchd1_genotype_3)

# Quick QC
library(scater)
sce <- addPerCellQC(sce)

plotColData(
  sce,
  "sum",
  x = "genotype.mouse",
  colour_by = "smchd1_genotype_3") +
  geom_hline(yintercept = 5e6, lty = 2, col = "red")
plotColData(
  sce,
  "detected",
  x = "genotype.mouse",
  colour_by = "smchd1_genotype_3")
plotExpression(sce, "Xist", x = "mouse_number", colour_by = "sex", exprs_values = "counts")
plotExpression(sce, "Smchd1", x = "genotype.mouse", colour_by = "smchd1_genotype_3", exprs_values = "counts")

sce <- logNormCounts(sce)

keep <- !isOutlier(sce$sum, log = TRUE)
sce <- sce[, keep]

library(scran)
stats <- modelGeneVar(sce)
hvgs <- getTopHVGs(stats, n = 1000)
str(hvgs)

sce <- runMDS(sce, subset_row = hvgs)
plotMDS(sce, colour_by = "smchd1_genotype_3")
plotMDS(sce, colour_by = "sex")
plotMDS(sce, colour_by = "mouse_number")

sex_set <- rownames(sce)[any(rowData(sce)$ENSEMBL.SEQNAME %in% c("X", "Y"))]
sce <- runMDS(sce, subset_row = setdiff(hvgs, sex_set))
plotMDS(sce, colour_by = "smchd1_genotype_3")
plotMDS(sce, colour_by = "sex")
plotMDS(sce, colour_by = "mouse_number")
