# Prepare C086_Kinkel data for GEO submission
# Peter Hickey
# 2021-09-21

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
sce <- readRDS(here("data", "SCEs", "C086_Kinkel.preprocessed.SCE.rds"))
# Restrict analysis to read counts
assays(sce) <- list(counts = assay(sce, "read_counts"))

# TODO: Any sample-level filtering? If so, will need to modify FASTQ files.

# Gene counts
write.csv(
  x = as.data.frame(as.matrix(counts(sce))),
  file = gzfile(file.path(outdir, "SCE", "gene_counts.csv.gz")),
  row.names = TRUE)

# ERCC counts
write.csv(
  x = as.data.frame(as.matrix(counts(altExp(sce)))),
  file = gzfile(file.path(outdir, "SCE", "ERCC_counts.csv.gz")),
  row.names = TRUE)

# colData
write.csv(
  x = as.data.frame(colData(sce)),
  file = gzfile(file.path(outdir, "SCE", "sample_sheet.csv.gz")),
  row.names = TRUE)
