# Part 1 — Preparation of RNA-seq matrices for WGCNA

## Required R packages

DESeq2 was used for variance-stabilizing transformation, RUVSeq for estimation of unwanted variation, edgeR for the residual model required by `RUVr()`, limma for regression of the estimated unwanted factors, and WGCNA for the final sample and gene quality check.

```r
library(DESeq2)
library(RUVSeq)
library(edgeR)
library(limma)
library(WGCNA)
library(dplyr)
library(Biobase)
```

## Input files and retained parameters

The count matrices contain genes in rows and samples in columns. The metadata table contains the sample identifiers, species or race, wing, compartment, and developmental stage. The orthology table contains the reciprocal *H. melpomene*–*H. erato* gene assignments.

```r
metadata_file          <- "sample_info_Joe.txt"
erato_count_file       <- "counts_table_erato_nonStranded.txt"
melpomene_count_file   <- "counts_table_melp_nonStranded.txt"
orthology_file         <- "../homogy.final.txt"

minimum_mean_count <- 10
ruv_k_erato        <- 5
ruv_k_melpomene    <- 4
```

The values of `k` above are the species-specific choices labelled in the original RUVSeq block. The exploratory loop testing values from 1 to 7 is not included because it was used only to inspect the correction by PCA.


## Replace *H. melpomene* identifiers with *H. erato* orthologs

To place the two species in a common gene-name space, *H. melpomene* identifiers with reciprocal orthologs were replaced by the corresponding *H. erato* identifiers. Genes without an identified reciprocal ortholog retained their original names. 

```r
melpomene_ids <- rownames(melpomene_counts)

for (i in seq_len(nrow(orthologs))) {
  matched <- grepl(
    as.character(orthologs$melp.single.prot[i]),
    melpomene_ids,
    fixed = TRUE
  )

  melpomene_ids[matched] <- as.character(
    orthologs$erato.single.prot[i]
  )
}

melpomene_counts$gene_id <- melpomene_ids
melpomene_counts <- melpomene_counts %>%
  distinct(gene_id, .keep_all = TRUE)

rownames(melpomene_counts) <- melpomene_counts$gene_id
melpomene_counts$gene_id <- NULL

# Name used for WntA in the downstream comparative analysis.
rownames(melpomene_counts)[
  rownames(melpomene_counts) == "Wnt5"
] <- "WntA"
```

## Select forewing samples and filter low-count genes

Samples were selected using the metadata rather than fixed column positions. Genes with a mean count below 10 in the selected samples were removed before transformation

```r
select_forewing_data <- function(count_matrix, metadata, race_name) {
  selected_metadata <- metadata %>%
    filter(
      race == race_name,
      wing == "forewing",
      day %in% c("larval", "1", "2")
    )

  missing_samples <- setdiff(
    selected_metadata$sample_id,
    colnames(count_matrix)
  )

  if (length(missing_samples) > 0) {
    stop(
      "Samples missing from the count matrix: ",
      paste(missing_samples, collapse = ", ")
    )
  }

  selected_counts <- count_matrix[
    , selected_metadata$sample_id,
    drop = FALSE
  ]

  selected_metadata <- selected_metadata[
    match(colnames(selected_counts), selected_metadata$sample_id),
    ,
    drop = FALSE
  ]
  rownames(selected_metadata) <- selected_metadata$sample_id

  selected_counts <- round(selected_counts)
  keep_gene <- rowMeans(selected_counts, na.rm = TRUE) >= minimum_mean_count
  selected_counts <- selected_counts[keep_gene, , drop = FALSE]

  list(
    counts = selected_counts,
    metadata = selected_metadata
  )
}

erato_input <- select_forewing_data(
  erato_counts,
  metadata,
  race_name = "demophoon"
)

melpomene_input <- select_forewing_data(
  melpomene_counts,
  metadata,
  race_name = "rosina"
)
```

## Variance-stabilizing transformation and RUVSeq correction

Each species was processed separately. Counts were transformed with the DESeq2 variance-stabilizing transformation using `blind = FALSE`. `RUVr()` then estimated unwanted factors from deviance residuals. The estimated factors were regressed from the transformed matrix with `removeBatchEffect()`, while the tissue design matrix was retained so that biological variation associated with the sampled developmental tissues was not removed.

```r
prepare_wgcna_matrix <- function(count_matrix, metadata, k) {
  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = metadata,
    design = ~ tissue
  )

  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  vst_matrix <- assay(vsd)

  # This filtering was used only to estimate the unwanted factors.
  ruv_counts <- round(vst_matrix)
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
    ~ metadata$tissue,
    data = pData(ruv_set)
  )

  y <- DGEList(
    counts = counts(ruv_set),
    group = metadata$tissue
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
    ~ 0 + metadata$tissue
  )

  corrected_matrix <- removeBatchEffect(
    vst_matrix,
    covariates = unwanted_factors,
    design = biological_design
  )

  # WGCNA requires samples in rows and genes in columns.
  corrected_matrix <- t(corrected_matrix)

  quality <- goodSamplesGenes(corrected_matrix, verbose = 3)
  corrected_matrix <- corrected_matrix[
    quality$goodSamples,
    quality$goodGenes,
    drop = FALSE
  ]

  corrected_matrix
}

erato_wgcna <- prepare_wgcna_matrix(
  erato_input$counts,
  erato_input$metadata,
  k = ruv_k_erato
)

melpomene_wgcna <- prepare_wgcna_matrix(
  melpomene_input$counts,
  melpomene_input$metadata,
  k = ruv_k_melpomene
)
```

## Construct the developmental-stage matrices

Developmental stages were encoded as binary variables for the subsequent module–stage correlation analysis. Rows were placed in exactly the same order as the samples in the corrected expression matrices.

```r
make_stage_matrix <- function(metadata, sample_order) {
  metadata <- metadata[
    match(sample_order, metadata$sample_id),
    ,
    drop = FALSE
  ]

  data.frame(
    larval = as.integer(metadata$day == "larval"),
    day1 = as.integer(metadata$day == "1"),
    day2 = as.integer(metadata$day == "2"),
    row.names = sample_order,
    check.names = FALSE
  )
}

erato_stage <- make_stage_matrix(
  erato_input$metadata,
  rownames(erato_wgcna)
)

melpomene_stage <- make_stage_matrix(
  melpomene_input$metadata,
  rownames(melpomene_wgcna)
)
```


These four files are the inputs to the RNA network-construction and module-preservation analysis described in Part 3.
