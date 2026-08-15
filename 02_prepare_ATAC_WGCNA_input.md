# Part 2 — Preparation of ATAC-seq matrices for WGCNA

This document describes the preparation of the chromatin-accessibility matrices used for the *H. erato* and *H. melpomene* co-accessibility analyses. It begins with peak-level count tables generated after read trimming, alignment, alignment filtering, duplicate removal, and peak calling.

Before the R analysis, peaks were annotated to the closest RNA-expressed gene with HOMER and their coordinates were converted into the common pan-genome coordinate system. The matrices used for the final analysis were restricted to peaks associated with the forewing-development gene set identified by RNA-seq WGCNA.

## Required R packages

DESeq2 was used for variance-stabilizing transformation, RUVSeq for estimation of unwanted variation, edgeR for the residual model required by `RUVr()`, limma for regression of the estimated unwanted factors, and WGCNA for the final sample and peak quality check.

```r
library(DESeq2)
library(RUVSeq)
library(edgeR)
library(limma)
library(WGCNA)
library(Biobase)
library(dplyr)
library(matrixStats)
```

## Input files and retained parameters

The two peak-count tables contain the pan-genome coordinates, original genome coordinates, closest expressed gene, and sample-level ATAC-seq counts. The metadata files contain the sample identifier, developmental stage, and sampled wing region.

```r
erato_atac_file       <- "erato.full.ATAC.genes.counts.header.txt"
melpomene_atac_file   <- "melp.full.ATAC.genes.counts.header.txt"
erato_metadata_file   <- "header_samples_dem_hyd.txt"
melpomene_metadata_file <- "header_samples_melp_ros.txt"

minimum_mean_accessibility <- 30
ruv_k_atac <- 7
within_stage_variance_quantile <- 0.75

forewing_regions <- c("forewing", "FD", "FM", "FP")

```


## Load the annotated peak-count tables

The first three columns contain the pan-genome interval, the following three columns contain the corresponding original-genome interval, and `gene` contains the closest expressed gene. These seven annotation columns were removed from the numeric count matrix after constructing a unique peak identifier.

```r
read_atac_counts <- function(file) {
  atac_table <- read.table(
    file,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  required_columns <- c(
    "chrom", "start", "end",
    "chrom_ori", "start_ori", "end_ori",
    "gene"
  )

  missing_columns <- setdiff(required_columns, colnames(atac_table))
  if (length(missing_columns) > 0) {
    stop(
      "Missing ATAC annotation columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  peak_id <- paste(
    atac_table$chrom,
    atac_table$start,
    atac_table$end,
    atac_table$gene,
    sep = "_"
  )

  if (anyDuplicated(peak_id)) {
    stop("Duplicated pan-genome peak identifiers were detected")
  }

  peak_annotation <- atac_table[, required_columns, drop = FALSE]
  rownames(peak_annotation) <- peak_id

  sample_columns <- setdiff(colnames(atac_table), required_columns)
  count_matrix <- atac_table[, sample_columns, drop = FALSE]
  rownames(count_matrix) <- peak_id

  list(
    counts = count_matrix,
    annotation = peak_annotation
  )
}

erato_atac <- read_atac_counts(erato_atac_file)
melpomene_atac <- read_atac_counts(melpomene_atac_file)
```

## Load and standardize the sample metadata

The original metadata files were read without column headers and referred to the first three columns as `V1`, `V2`, and `V3`. They correspond to sample identifier, developmental stage, and wing region, respectively.

```r
read_atac_metadata <- function(file) {
  metadata <- read.table(
    file,
    header = FALSE,
    stringsAsFactors = FALSE
  )

  colnames(metadata)[1:3] <- c(
    "sample_id",
    "stage",
    "region"
  )

  metadata$stage_region <- paste(
    metadata$stage,
    metadata$region,
    sep = "_"
  )

  metadata
}

erato_metadata <- read_atac_metadata(erato_metadata_file)
melpomene_metadata <- read_atac_metadata(melpomene_metadata_file)

erato_metadata <- erato_metadata %>%
  filter(region %in% forewing_regions)

melpomene_metadata <- melpomene_metadata %>%
  filter(
    region %in% forewing_regions,
    !sample_id %in% excluded_melpomene_samples
  )
```

## Select samples and filter low-count peaks

Sample selection was performed from the metadata rather than by column position. Peaks with a mean count below 30 across the selected samples were removed

```r
select_atac_samples <- function(count_matrix, metadata) {
  missing_samples <- setdiff(
    metadata$sample_id,
    colnames(count_matrix)
  )

  if (length(missing_samples) > 0) {
    stop(
      "Samples missing from the ATAC count matrix: ",
      paste(missing_samples, collapse = ", ")
    )
  }

  selected_counts <- count_matrix[
    , metadata$sample_id,
    drop = FALSE
  ]

  metadata <- metadata[
    match(colnames(selected_counts), metadata$sample_id),
    ,
    drop = FALSE
  ]
  rownames(metadata) <- metadata$sample_id

  selected_counts <- round(selected_counts)
  keep_peak <- rowMeans(selected_counts, na.rm = TRUE) >=
    minimum_mean_accessibility

  list(
    counts = selected_counts[keep_peak, , drop = FALSE],
    metadata = metadata
  )
}

erato_input <- select_atac_samples(
  erato_atac$counts,
  erato_metadata
)

melpomene_input <- select_atac_samples(
  melpomene_atac$counts,
  melpomene_metadata
)
```

## Variance-stabilizing transformation and RUVSeq correction

Each species was processed separately. Counts were transformed with DESeq2 using `blind = FALSE`. As in the original script, the transformed values were shifted to be non-negative before constructing the RUVSeq object. `RUVr()` was run with seven unwanted factors in both species. The estimated factors were then regressed from the transformed matrix with `removeBatchEffect()`, while developmental-stage effects were retained in the design matrix.

```r
prepare_atac_wgcna_matrix <- function(count_matrix, metadata, k = 7) {
  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = metadata,
    design = ~ stage_region
  )

  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  vst_matrix <- assay(vsd)

  # RUVSeq input used in the original analysis.
  shifted_matrix <- vst_matrix + abs(min(vst_matrix))
  ruv_counts <- round(shifted_matrix)
  ruv_keep <- rowSums(ruv_counts > 5) > 2
  ruv_counts <- ruv_counts[ruv_keep, , drop = FALSE]

  ruv_set <- newSeqExpressionSet(
    as.matrix(ruv_counts),
    phenoData = data.frame(
      metadata,
      row.names = colnames(ruv_counts)
    )
  )

  ruv_set <- betweenLaneNormalization(
    ruv_set,
    which = "upper"
  )

  residual_design <- model.matrix(
    ~ stage_region,
    data = pData(ruv_set)
  )

  y <- DGEList(
    counts = counts(ruv_set),
    group = metadata$stage_region
  )
  y <- calcNormFactors(y, method = "upperquartile")
  y <- estimateGLMCommonDisp(y, residual_design)
  y <- estimateGLMTagwiseDisp(y, residual_design)

  fit <- glmFit(y, residual_design)
  residuals_deviance <- residuals(fit, type = "deviance")

  ruv_fit <- RUVr(
    ruv_set,
    rownames(ruv_counts),
    k = k,
    residuals_deviance
  )

  ruv_metadata <- pData(ruv_fit)
  unwanted_columns <- setdiff(
    colnames(ruv_metadata),
    colnames(metadata)
  )
  unwanted_factors <- as.matrix(
    ruv_metadata[, unwanted_columns, drop = FALSE]
  )

  biological_design <- model.matrix(
    ~ 0 + metadata$stage
  )

  corrected_matrix <- removeBatchEffect(
    shifted_matrix,
    covariates = unwanted_factors,
    design = biological_design
  )

  # WGCNA requires samples in rows and peaks in columns.
  corrected_matrix <- t(corrected_matrix)

  quality <- goodSamplesGenes(corrected_matrix, verbose = 3)
  corrected_matrix <- corrected_matrix[
    quality$goodSamples,
    quality$goodGenes,
    drop = FALSE
  ]

  corrected_matrix
}

erato_wgcna_atac <- prepare_atac_wgcna_matrix(
  erato_input$counts,
  erato_input$metadata,
  k = ruv_k_atac
)

melpomene_wgcna_atac <- prepare_atac_wgcna_matrix(
  melpomene_input$counts,
  melpomene_input$metadata,
  k = ruv_k_atac
)
```

## Remove peaks dominated by within-stage variation

The original script calculated peak-wise variance separately within fifth-instar larvae, pupal day 1, and pupal day 2. A peak was retained when its within-stage variance was below the 75th percentile in at least one developmental stage. This removes peaks that show extremely high replicate-to-replicate variation at every sampled stage.

```r
retain_by_within_stage_variance <- function(
  accessibility_matrix,
  metadata,
  probability = 0.75
) {
  metadata <- metadata[
    match(rownames(accessibility_matrix), metadata$sample_id),
    ,
    drop = FALSE
  ]

  retained_peaks <- lapply(unique(metadata$stage), function(stage_name) {
    stage_matrix <- accessibility_matrix[
      metadata$stage == stage_name,
      ,
      drop = FALSE
    ]

    peak_variance <- colVars(as.matrix(stage_matrix))
    variance_limit <- quantile(
      peak_variance,
      probabilities = probability,
      na.rm = TRUE
    )

    colnames(stage_matrix)[peak_variance < variance_limit]
  })

  retained_peaks <- unique(unlist(retained_peaks))
  accessibility_matrix[, retained_peaks, drop = FALSE]
}

erato_wgcna_atac <- retain_by_within_stage_variance(
  erato_wgcna_atac,
  erato_input$metadata,
  probability = within_stage_variance_quantile
)

melpomene_wgcna_atac <- retain_by_within_stage_variance(
  melpomene_wgcna_atac,
  melpomene_input$metadata,
  probability = within_stage_variance_quantile
)
```
