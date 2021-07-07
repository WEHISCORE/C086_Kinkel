# Helper function to Combine data from 2 SCEs using gene names.
# NOTE: This assumes more than I'd like about the rowData and doesn't do much
#       checking of these assumptions.
.combine <- function(x, y, rowData_by = c("ENSEMBL", "SYMBOL", "CHR")) {
  if (is.null(rowData_by)) {
    rowData <- dplyr::full_join(
      as.data.frame(rowData(x)) %>%
        tibble::rownames_to_column(var = "gene"),
      as.data.frame(rowData(y)) %>%
        tibble::rownames_to_column(var = "gene")) %>%
      tibble::column_to_rownames("gene") %>%
      DataFrame(., row.names = rownames(.))
  } else {
    rowData <- dplyr::full_join(
      as.data.frame(rowData(x)[, rowData_by, drop = FALSE]),
      as.data.frame(rowData(y)[, rowData_by, drop = FALSE]),
      by = rowData_by) %>%
      DataFrame(row.names = scater::uniquifyFeatureNames(
        .$ENSEMBL,
        .$SYMBOL))
    rownames(x) <- rownames(rowData)[match(rowData(x)$ENSEMBL, rowData$ENSEMBL)]
    rownames(y) <- rownames(rowData)[match(rowData(y)$ENSEMBL, rowData$ENSEMBL)]
  }

  colData <- rbind(colData(x), colData(y))

  counts <- matrix(
    data = 0L,
    nrow = nrow(rowData), ncol = nrow(colData),
    dimnames = list(rownames(rowData), rownames(colData)))
  counts[rownames(x), colnames(x)] <- counts(
    x,
    withDimnames = FALSE)
  counts[rownames(y), colnames(y)] <- counts(
    y,
    withDimnames = FALSE)

  stopifnot(
    identical(
      metadata(x)$scPipe$version,
      metadata(y)$scPipe$version))
  stopifnot(
    identical(
      metadata(x)$scPipe$QC_cols,
      metadata(y)$scPipe$QC_cols))
  stopifnot(
    identical(
      metadata(x)$scPipe$demultiplex_info$status,
      metadata(y)$scPipe$demultiplex_info$status))
  stopifnot(
    identical(
      metadata(x)$scPipe$UMI_dup_info$duplication.number,
      metadata(y)$scPipe$UMI_dup_info$duplication.number))
  stopifnot(identical(metadata(x)$Biomart, metadata(y)$Biomart))
  metadata <- list(
    scPipe = list(
      version = metadata(x)$scPipe$version,
      QC_cols = metadata(x)$scPipe$QC_cols,
      demultiplex_info = data.frame(
        status = metadata(x)$scPipe$demultiplex_info$status,
        count = metadata(x)$scPipe$demultiplex_info$count +
          metadata(y)$scPipe$demultiplex_info$count),
      UMI_dup_info = data.frame(
        duplication.number = metadata(
          x)$scPipe$UMI_dup_info$duplication.number,
        count = metadata(x)$scPipe$UMI_dup_info$count +
          metadata(y)$scPipe$UMI_dup_info$count)),
    Biomart = metadata(x)$Biomart)

  sce <- SingleCellExperiment(
    rowData = rowData,
    colData = colData,
    assays = list(counts = counts),
    metadata = metadata)

  stopifnot(identical(int_metadata(x), int_metadata(y)))
  int_metadata(sce) <- int_metadata(x)

  # NOTE: Not trying to combine int_elementMetadata of objects. Each is a
  #       DataFrame with a zero-column DataFrame as a `rowPairs` column. This
  #       is effectively no data and the SCE constructor makes one, anyway.

  stopifnot(validObject(sce))
  sce
}

# NOTE: Need to use my own gene counting function because not using UMI
#       deduplication.
geneCountingNoUMIDedup <- function(outdir, bc_anno) {
  files <- list.files(file.path(outdir, "count"), full.names = TRUE)
  names(files) <- sub("\\.csv", "", basename(files))
  counts <- lapply(files, function(file) {
    message(basename(file))
    data.table::fread(file, select = 1)[, table(gene_id)]
  })
  genes <- Reduce(union, lapply(counts, names))
  x <- matrix(
    0L,
    nrow = length(genes),
    ncol = length(files),
    dimnames = list(genes, names(counts)))
  for (j in names(counts)) {
    xx <- counts[[j]]
    x[names(xx), j] <- xx
  }
  z <- cbind(
    data.frame(gene_id = rownames(x)),
    as.data.frame(x))
  data.table::fwrite(
    x = z,
    file = file.path(paste0(outdir, "_no_dedup"), "gene_count.csv"),
    row.names = FALSE,
    nThread = 1)
}

# Take a DataFrame with AtomicList columns and return a DataFrame where these
# columns have been flattened by paste-ing together the elements separated by
# `sep`.
flattenDF <- function(x, sep = "; ") {
  DataFrame(
    endoapply(x, function(xx) {
      if (!is(xx, "AtomicList")) {
        return(xx)
      }
      unstrsplit(as(xx, "CharacterList"), sep = sep)
    }),
    row.names = rownames(x))
}

filterCountMatrix <- function(sce, exprs_values, min.count, method) {
  y <- assay(sce, exprs_values)

  if (method == "standard") {
    keep <- filterByExpr(
      y,
      group = sce$smchd1_genotype_updated,
      min.count = min.count)
  } else if (method == "careful") {
    keeps <- lapply(levels(sce$smchd1_genotype_updated), function(j) {
      rowSums(
        y[, sce$smchd1_genotype_updated == j, drop = FALSE] >= min.count) >=
        sum(sce$smchd1_genotype_updated == j)
    })
    keep <- Reduce("|", keeps)
  } else {
    stop("Unknown method")
  }
  sce[keep, ]
}

deterministicImputation <- function(sce, exprs_values, method) {
  y <- assay(sce, exprs_values)

  # NOTE: Feel like there's a more efficient way do this (perhaps matrix algebra
  #       or statmod::vecmat()/statmod::matvec()). In any case, still much more
  #       efficient than Luke's implementation that explicitly loops over rows
  #       of the matrix.
  if (method == "all") {
    # NOTE: This treats all samples as coming from a single biological group.
    num <- rowSums(y)
    N <- matrix(
      colSums(y),
      nrow = nrow(y),
      ncol = ncol(y),
      byrow = TRUE,
      dimnames = dimnames(y))
    I <- y > 0
    denom <- rowSums(N * I)
    pi_nonzero <- matrix(
      num / denom,
      ncol = nrow(y),
      dimnames = list(NULL, rownames(y)))
    # NOTE: Don't know why I need `[,]` but I do.
    y_imp <- pi_nonzero[, ] * N[, ]
    y2 <- y_imp
    y2[which(I > 0)] <- y[which(I > 0)]
  } else if (method == "group-specific") {
    # NOTE: Adapted from Luke's code.
    # NOTE: This imputes samples in each biological group separately.
    # NOTE: This could be further sped up by pre-allocating a matrix, `y2`,
    #       and filling it directly.
    y2s <- lapply(levels(sce$smchd1_genotype_updated), function(j) {
      y <- y[, sce$smchd1_genotype_updated == j, drop = FALSE]
      num <- rowSums(y)
      N <- matrix(
        colSums(y),
        nrow = nrow(y),
        ncol = ncol(y),
        byrow = TRUE,
        dimnames = dimnames(y))
      I <- y > 0
      denom <- rowSums(N * I)
      # NOTE: If all samples have zero counts for a particular gene then that
      #       value of denom will be 0. To avoid 0 / 0 (i.e. NaN) values in
      #       pi_nonzero we replace these by an arbitrary value.
      denom[denom == 0] <- 1L
      pi_nonzero <- matrix(
        num / denom,
        ncol = nrow(y),
        dimnames = list(NULL, rownames(y)))
      # NOTE: Don't know why I need `[,]` but I do.
      y_imp <- pi_nonzero[,] * N[,]
      if (is.null(dim(y_imp))) {
        # NOTE: Needed for WT when using summed tech reps (n=1).
        dim(y_imp) <- dim(y)
        dimnames(y_imp) <- dimnames(y)
      }
      y2 <- y_imp
      y2[which(I > 0)] <- y[which(I > 0)]
      y2
    })
    y2 <- do.call(cbind, y2s)
    y2 <- y2[, colnames(y)]
  } else {
    stop("Unknown 'method'")
  }
  assay(sce, paste0("imputed_", exprs_values)) <- y2
  sce
}
