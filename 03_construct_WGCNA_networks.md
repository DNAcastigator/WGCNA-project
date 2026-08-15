# Part 3 — WGCNA module construction and developmental-stage associations

This document describes the construction of four independent weighted networks: *H. erato* RNA-seq, *H. melpomene* RNA-seq, *H. erato* ATAC-seq, and *H. melpomene* ATAC-seq. The input matrices are the corrected, forewing-specific matrices produced in Parts 1 and 2, with samples in rows and genes or ATAC-seq peaks in columns.

The code below retains the final network settings recoverable from the supplied scripts. Exploratory parameter searches, overwritten objects, individual candidate-gene queries, and plotting experiments were excluded. Cross-species module preservation is intentionally handled in Part 4.

## Required R package

```r
library(WGCNA)

options(stringsAsFactors = FALSE)
```

## Input files

The RNA-seq matrices use the common gene-name space established through reciprocal orthology assignment in Part 1. The ATAC-seq matrices use the pan-genome peak identifiers and RNA-associated peak restriction described in Part 2.

```r
input_files <- list(
  erato_RNA = list(
    matrix = "erato_FW_filtered_all.txt",
    traits = "erato_FW_day_trait.txt"
  ),
  melpomene_RNA = list(
    matrix = "melpomene_FW_filtered_all.txt",
    traits = "melpomene_FW_day_trait.txt"
  ),
  erato_ATAC = list(
    matrix = "erato_FW_filtered_ATAC_updated_conv-WN.txt",
    traits = "erato_FW_day_trait_ATAC.txt"
  ),
  melpomene_ATAC = list(
    matrix = "melpomene_FW_filtered_ATAC_updated_conv-WN.txt",
    traits = "melpomene_FW_day_trait_ATAC.txt"
  )
)
```

## Dataset-specific WGCNA parameters

Soft-thresholding powers were evaluated over powers 1–10 and the even powers 12–30 with `pickSoftThreshold()`, using a signed-hybrid network and biweight midcorrelation. The selected power was the first value reaching a signed scale-free topology fit of at least 0.80. The powers used in the final downstream scripts were 18 and 20 for the *H. erato* and *H. melpomene* RNA networks, respectively, and 4 for both ATAC networks.

The remaining module-detection settings reproduce the final RNA and ATAC network blocks. These values should also be reported in Supplementary Table S6.

| Dataset | Power | Maximum block size | Minimum module size | Merge cut height | Deep split | Correlation | Network type |
| --- | ---: | ---: | ---: | ---: | ---: | --- | --- |
| *H. erato* RNA | 18 | 6,000 | 30 | 0.25 | 3 | bicor | signed hybrid |
| *H. melpomene* RNA | 20 | 6,000 | 30 | 0.25 | 3 | bicor | signed hybrid |
| *H. erato* ATAC | 4 | 6,000 | 30 | 0.30 | 3 | bicor | signed hybrid |
| *H. melpomene* ATAC | 4 | 6,000 | 30 | 0.30 | 3 | bicor | signed hybrid |

```r
network_parameters <- data.frame(
  dataset = c(
    "erato_RNA",
    "melpomene_RNA",
    "erato_ATAC",
    "melpomene_ATAC"
  ),
  power = c(18, 20, 4, 4),
  max_block_size = c(6000, 6000, 6000, 6000),
  min_module_size = c(30, 30, 30, 30),
  merge_cut_height = c(0.25, 0.25, 0.30, 0.30),
  deep_split = c(3, 3, 3, 3),
  stringsAsFactors = FALSE
)

rownames(network_parameters) <- network_parameters$dataset
```

## Load and align the expression/accessibility and trait matrices

Both the molecular matrix and trait matrix must contain sample identifiers as row names. The checks below prevent traits from being silently correlated with samples in the wrong order.

```r
read_numeric_table <- function(file) {
  data <- read.table(
    file,
    header = TRUE,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  data[] <- lapply(data, as.numeric)
  data
}

load_wgcna_dataset <- function(matrix_file, trait_file) {
  dat_expr <- read_numeric_table(matrix_file)
  dat_traits <- read_numeric_table(trait_file)

  if (anyDuplicated(rownames(dat_expr))) {
    stop("Duplicated sample identifiers in ", matrix_file)
  }

  if (anyDuplicated(rownames(dat_traits))) {
    stop("Duplicated sample identifiers in ", trait_file)
  }

  missing_traits <- setdiff(
    rownames(dat_expr),
    rownames(dat_traits)
  )

  if (length(missing_traits) > 0) {
    stop(
      "Samples missing from the trait matrix: ",
      paste(missing_traits, collapse = ", ")
    )
  }

  dat_traits <- dat_traits[
    match(rownames(dat_expr), rownames(dat_traits)),
    ,
    drop = FALSE
  ]

  quality <- goodSamplesGenes(dat_expr, verbose = 3)

  if (!quality$allOK) {
    dat_expr <- dat_expr[
      quality$goodSamples,
      quality$goodGenes,
      drop = FALSE
    ]

    dat_traits <- dat_traits[
      rownames(dat_expr),
      ,
      drop = FALSE
    ]
  }

  list(
    expression = dat_expr,
    traits = dat_traits
  )
}

wgcna_data <- lapply(input_files, function(files) {
  load_wgcna_dataset(
    matrix_file = files$matrix,
    trait_file = files$traits
  )
})
```

## Evaluate the soft-thresholding power

The signed scale-free topology fit is calculated as in the supplied plotting scripts. The table produced for each dataset documents the full set of tested powers and the first power satisfying the 0.80 criterion.

```r
candidate_powers <- c(1:10, seq(12, 30, by = 2))
scale_free_fit_threshold <- 0.80

evaluate_soft_threshold <- function(dat_expr) {
  fit <- pickSoftThreshold(
    dat_expr,
    powerVector = candidate_powers,
    networkType = "signed hybrid",
    corFnc = "bicor",
    verbose = 5
  )

  fit_table <- as.data.frame(fit$fitIndices)
  fit_table$signed_R2 <- -sign(fit_table$slope) * fit_table$SFT.R.sq

  passing <- fit_table$Power[
    fit_table$signed_R2 >= scale_free_fit_threshold
  ]

  selected_power <- if (length(passing) == 0) {
    NA_real_
  } else {
    passing[1]
  }

  list(
    fit = fit,
    table = fit_table,
    first_passing_power = selected_power
  )
}

soft_threshold_results <- lapply(
  wgcna_data,
  function(dataset) {
    evaluate_soft_threshold(dataset$expression)
  }
)

soft_threshold_summary <- data.frame(
  dataset = names(soft_threshold_results),
  first_passing_power = vapply(
    soft_threshold_results,
    function(result) result$first_passing_power,
    numeric(1)
  ),
  power_used = network_parameters[
    names(soft_threshold_results),
    "power"
  ],
  row.names = NULL
)

soft_threshold_summary
```

The `first_passing_power` and `power_used` columns should agree. Retaining both columns makes the selection criterion auditable instead of reporting only the final power.

## Construct the four networks

Networks were generated independently for each species and assay. `blockwiseModules()` used biweight midcorrelation, a signed-hybrid network, a signed topological overlap matrix, and color module labels. The grey class contains features not assigned to a module and is retained in the network object but excluded from module–trait reporting.

```r
construct_network <- function(dat_expr, parameters) {
  blockwiseModules(
    dat_expr,
    power = parameters$power,
    maxBlockSize = parameters$max_block_size,
    networkType = "signed hybrid",
    TOMType = "signed",
    minModuleSize = parameters$min_module_size,
    reassignThreshold = 0,
    mergeCutHeight = parameters$merge_cut_height,
    numericLabels = FALSE,
    corType = "bicor",
    deepSplit = parameters$deep_split,
    verbose = 3
  )
}

networks <- lapply(names(wgcna_data), function(dataset_name) {
  construct_network(
    dat_expr = wgcna_data[[dataset_name]]$expression,
    parameters = network_parameters[dataset_name, ]
  )
})

names(networks) <- names(wgcna_data)
```

## Summarize module membership

The complete feature-to-module assignments and module sizes are retained for reporting and downstream analyses.

```r
summarize_modules <- function(network) {
  assignments <- data.frame(
    feature = names(network$colors),
    module = unname(network$colors),
    stringsAsFactors = FALSE
  )

  sizes <- as.data.frame(table(assignments$module))
  colnames(sizes) <- c("module", "size")
  sizes <- sizes[order(sizes$size, decreasing = TRUE), ]

  list(
    assignments = assignments,
    sizes = sizes
  )
}

module_summaries <- lapply(networks, summarize_modules)
```

## Calculate module–developmental-stage relationships

Module eigengenes were correlated with the three binary developmental-stage variables using Pearson correlation, as in the supplied scripts. Student correlation p-values were adjusted separately within each developmental-stage column using the Benjamini–Hochberg false-discovery-rate procedure. Associations with an adjusted p-value below 0.01 were treated as statistically supported.

```r
module_trait_alpha <- 0.01

calculate_module_trait_relationships <- function(
  dat_expr,
  dat_traits,
  network,
  alpha = 0.01
) {
  module_colors <- network$colors

  eigengenes <- moduleEigengenes(
    dat_expr,
    colors = module_colors
  )$eigengenes

  eigengenes <- orderMEs(eigengenes)

  module_trait_cor <- cor(
    eigengenes,
    dat_traits,
    use = "pairwise.complete.obs",
    method = "pearson"
  )

  module_trait_p <- corPvalueStudent(
    module_trait_cor,
    nSamples = nrow(dat_expr)
  )

  module_trait_fdr <- apply(
    module_trait_p,
    2,
    p.adjust,
    method = "BH"
  )

  module_trait_fdr <- matrix(
    module_trait_fdr,
    nrow = nrow(module_trait_cor),
    ncol = ncol(module_trait_cor),
    dimnames = dimnames(module_trait_cor)
  )

  # Exclude the unassigned grey class from reported modules.
  keep_module <- rownames(module_trait_cor) != "MEgrey"
  module_trait_cor <- module_trait_cor[keep_module, , drop = FALSE]
  module_trait_p <- module_trait_p[keep_module, , drop = FALSE]
  module_trait_fdr <- module_trait_fdr[keep_module, , drop = FALSE]

  significant <- module_trait_fdr < alpha

  list(
    eigengenes = eigengenes,
    correlation = module_trait_cor,
    p_value = module_trait_p,
    fdr = module_trait_fdr,
    significant = significant
  )
}

module_trait_results <- lapply(names(wgcna_data), function(dataset_name) {
  calculate_module_trait_relationships(
    dat_expr = wgcna_data[[dataset_name]]$expression,
    dat_traits = wgcna_data[[dataset_name]]$traits,
    network = networks[[dataset_name]],
    alpha = module_trait_alpha
  )
})

names(module_trait_results) <- names(wgcna_data)
```

## Create complete module–trait result tables

Each row of the output contains a module–stage correlation, its raw p-value, its FDR-adjusted p-value, and the resulting significance classification. All tested module–stage combinations are saved; non-significant results are not discarded.

```r
make_module_trait_table <- function(result) {
  correlation_table <- as.data.frame(as.table(result$correlation))
  p_value_table <- as.data.frame(as.table(result$p_value))
  fdr_table <- as.data.frame(as.table(result$fdr))

  data.frame(
    module = sub("^ME", "", correlation_table$Var1),
    stage = correlation_table$Var2,
    correlation = correlation_table$Freq,
    p_value = p_value_table$Freq,
    fdr = fdr_table$Freq,
    significant = fdr_table$Freq < module_trait_alpha,
    stringsAsFactors = FALSE
  )
}

module_trait_tables <- lapply(
  module_trait_results,
  make_module_trait_table
)
```

## Plot module–stage relationships

The heatmap displays correlations for every assigned module. Text is shown only for FDR-supported associations, avoiding the visual removal of non-significant modules from the underlying result table.

```r
plot_module_trait_heatmap <- function(
  result,
  module_sizes,
  title
) {
  correlation_matrix <- result$correlation
  fdr_matrix <- result$fdr

  module_names <- sub("^ME", "", rownames(correlation_matrix))
  size_lookup <- setNames(module_sizes$size, module_sizes$module)

  row_labels <- paste0(
    module_names,
    " (",
    size_lookup[module_names],
    ")"
  )

  text_matrix <- matrix(
    "",
    nrow = nrow(correlation_matrix),
    ncol = ncol(correlation_matrix),
    dimnames = dimnames(correlation_matrix)
  )

  supported <- fdr_matrix < module_trait_alpha
  text_matrix[supported] <- paste0(
    signif(correlation_matrix[supported], 2),
    "\n(q=",
    signif(fdr_matrix[supported], 2),
    ")"
  )

  labeledHeatmap(
    Matrix = correlation_matrix,
    xLabels = colnames(correlation_matrix),
    yLabels = row_labels,
    ySymbols = row_labels,
    colorLabels = FALSE,
    colors = blueWhiteRed(50),
    textMatrix = text_matrix,
    setStdMargins = FALSE,
    cex.text = 0.7,
    zlim = c(-1, 1),
    main = title
  )
}
```

## Save network objects and tabular results

The saved R object contains the input matrices, selected parameters, network assignments, module eigengenes, and complete module–trait statistics needed for preservation and WntA-centered analyses.

```r
output_directory <- "wgcna_results"
dir.create(output_directory, showWarnings = FALSE)

for (dataset_name in names(wgcna_data)) {
  write.table(
    soft_threshold_results[[dataset_name]]$table,
    file = file.path(
      output_directory,
      paste0(dataset_name, "_soft_threshold.tsv")
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    module_summaries[[dataset_name]]$assignments,
    file = file.path(
      output_directory,
      paste0(dataset_name, "_module_assignments.tsv")
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    module_summaries[[dataset_name]]$sizes,
    file = file.path(
      output_directory,
      paste0(dataset_name, "_module_sizes.tsv")
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    module_trait_tables[[dataset_name]],
    file = file.path(
      output_directory,
      paste0(dataset_name, "_module_trait_relationships.tsv")
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  pdf(
    file = file.path(
      output_directory,
      paste0(dataset_name, "_module_trait_heatmap.pdf")
    ),
    width = 8,
    height = 10
  )

  plot_module_trait_heatmap(
    result = module_trait_results[[dataset_name]],
    module_sizes = module_summaries[[dataset_name]]$sizes,
    title = paste(
      gsub("_", " ", dataset_name),
      "developmental stage"
    )
  )

  dev.off()
}

saveRDS(
  list(
    data = wgcna_data,
    parameters = network_parameters,
    soft_thresholds = soft_threshold_results,
    networks = networks,
    module_summaries = module_summaries,
    module_trait_results = module_trait_results
  ),
  file = file.path(
    output_directory,
    "four_WGCNA_networks.rds"
  )
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(output_directory, "sessionInfo.txt")
)
```

The `four_WGCNA_networks.rds` object is the direct input for the module-preservation analysis in Part 4. The complete module–trait tables can also be used to construct Supplementary Table S3 without rerunning the networks.
