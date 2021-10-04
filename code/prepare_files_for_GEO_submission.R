# Prepare C086_Kinkel data for GEO submission
# Peter Hickey
# 2021-10-04

library(here)

outdir <- here("GEO")
dir.create(outdir, recursive = TRUE)

# FASTQs -----------------------------------------------------------------------

dir.create(file.path(outdir, "FASTQ"))
# NOTE: These plate-level FASTQ files are created by code/scPipe.R
file.copy(
  from = here("extdata/NN206/scPipe/LCE504/LCE504.R1.fastq.gz"),
  to = file.path(outdir, "FASTQ"),
  recursive = FALSE,
  overwrite = FALSE)
file.copy(
  from = here("extdata/NN206/scPipe/LCE504/LCE504.R2.fastq.gz"),
  to = file.path(outdir, "FASTQ"),
  recursive = FALSE,
  overwrite = FALSE)

# SCEs -------------------------------------------------------------------------

dir.create(file.path(outdir, "SCE"))
# NOTE: Starting from the SCE used for the DE analysis.
sce <- readRDS(here("data", "SCEs", "C086_Kinkel.preprocessed.SCE.rds"))
# NOTE: Revert some of the changes to the rowData.
rownames(sce) <- rowData(sce)$ENSEMBL.GENEID
rowData(sce) <- S4Vectors::make_zero_col_DFrame(nrow(sce))

# Gene counts
# NOTE: Exporting read counts.
write.csv(
  x = as.data.frame(as.matrix(assay(sce, "read_counts"))),
  file = gzfile(file.path(outdir, "SCE", "gene_counts.csv.gz")),
  row.names = TRUE)

# ERCC counts
# NOTE: Exporting read counts.
write.csv(
  x = as.data.frame(as.matrix(assay(altExp(sce, "ERCC"), "read_counts"))),
  file = gzfile(file.path(outdir, "SCE", "ERCC_counts.csv.gz")),
  row.names = TRUE)

# colData
write.csv(
  x = as.data.frame(colData(sce)),
  file = gzfile(file.path(outdir, "SCE", "sample_sheet.csv.gz")),
  row.names = TRUE)
