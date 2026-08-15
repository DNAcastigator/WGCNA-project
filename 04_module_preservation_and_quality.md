# Part 4 — Module preservation and within-reference module quality

This document describes the cross-species module-preservation analyses for the RNA-seq and ATAC-seq networks constructed in Part 3. It also extracts the within-reference quality statistics calculated by `modulePreservation()`, allowing module robustness to be reported directly rather than inferred only from the preservation plots.

Two distinct questions are addressed:

1. **Module quality:** How well defined is each module within the reference dataset in which it was constructed?
2. **Module preservation:** How well is the organization of that reference module retained in the other species?

The quality and preservation statistics are related outputs of the same WGCNA function, but they should not be interpreted as the same quantity. A high-quality module in the reference species is not necessarily preserved in the comparison species.

## Required R packages

```r
library(WGCNA)
library(ggplot2)

options(stringsAsFactors = FALSE)
```

## Load the networks generated in Part 3

```r
wgcna_object_file <- file.path(
  "wgcna_results",
  "four_WGCNA_networks.rds"
)

wgcna_results <- readRDS(wgcna_object_file)

required_datasets <- c(
  "erato_RNA",
  "melpomene_RNA",
  "erato_ATAC",
  "melpomene_ATAC"
)

missing_datasets <- setdiff(
  required_datasets,
  names(wgcna_results$data)
)

if (length(missing_datasets) > 0) {
  stop(
    "Datasets missing from the Part 3 object: ",
    paste(missing_datasets, collapse = ", ")
  )
}
```

## Module-preservation settings

The analyses used signed-hybrid networks, biweight midcorrelation, 100 permutations, a maximum module size of 10,000, and a maximum gold-module size of 500. Both species were used in turn as the reference network.

| Assay | Correlation used for preservation | Permutations | Network type | Maximum module size | Maximum gold-module size |
| --- | --- | ---: | --- | ---: | ---: |
| RNA-seq | bicor | 100 | signed hybrid | 10,000 | 500 |
| ATAC-seq | bicor | 100 | signed hybrid | 10,000 | 500 |

The same correlation and permutation settings were used for the within-reference module-quality checks.

```r
preservation_settings <- list(
  RNA = list(
    cor_function = "bicor",
    permutations = 100
  ),
  ATAC = list(
    cor_function = "bicor",
    permutations = 100
  )
)

quality_settings <- list(
  cor_function = "bicor",
  permutations = 100
)

preservation_random_seed <- 12345
```

## Helper functions for expression matrices and module colors

`modulePreservation()` matches features by column name. Cross-species comparisons therefore use only genes or pan-genome peaks present in both species, while internal quality checks use the complete feature set of each individual network.

```r
get_expression <- function(dataset_name) {
  wgcna_results$data[[dataset_name]]$expression
}

get_module_colors <- function(dataset_name) {
  colors <- wgcna_results$networks[[dataset_name]]$colors

  if (is.null(names(colors))) {
    stop("Unnamed module-color vector for ", dataset_name)
  }

  colors
}

align_colors_to_expression <- function(dat_expr, colors) {
  missing_colors <- setdiff(colnames(dat_expr), names(colors))

  if (length(missing_colors) > 0) {
    stop(
      "Features missing module assignments: ",
      paste(head(missing_colors, 10), collapse = ", ")
    )
  }

  colors[colnames(dat_expr)]
}
```

## Identify the features shared between species

RNA features are matched in the common ortholog-name space established in Part 1. ATAC features are matched using the pan-genome peak identifiers established in Part 2.

```r
prepare_cross_species_data <- function(assay) {
  erato_name <- paste0("erato_", assay)
  melpomene_name <- paste0("melpomene_", assay)

  erato_expression <- get_expression(erato_name)
  melpomene_expression <- get_expression(melpomene_name)

  shared_features <- intersect(
    colnames(erato_expression),
    colnames(melpomene_expression)
  )

  if (length(shared_features) == 0) {
    stop("No shared features were found for ", assay)
  }

  erato_expression <- erato_expression[
    , shared_features,
    drop = FALSE
  ]

  melpomene_expression <- melpomene_expression[
    , shared_features,
    drop = FALSE
  ]

  erato_colors <- align_colors_to_expression(
    erato_expression,
    get_module_colors(erato_name)
  )

  melpomene_colors <- align_colors_to_expression(
    melpomene_expression,
    get_module_colors(melpomene_name)
  )

  list(
    multi_expression = list(
      erato = list(data = erato_expression),
      melpomene = list(data = melpomene_expression)
    ),
    multi_color = list(
      erato = erato_colors,
      melpomene = melpomene_colors
    ),
    shared_features = shared_features
  )
}

cross_species_data <- list(
  RNA = prepare_cross_species_data("RNA"),
  ATAC = prepare_cross_species_data("ATAC")
)
```

## Run the cross-species preservation analyses

Each species is treated as the reference in turn. Consequently, preservation of *H. erato* modules in *H. melpomene* and preservation of *H. melpomene* modules in *H. erato* are both calculated.

```r
run_cross_species_preservation <- function(
  prepared_data,
  cor_function,
  permutations
) {
  modulePreservation(
    multiData = prepared_data$multi_expression,
    multiColor = prepared_data$multi_color,
    referenceNetworks = c(1, 2),
    networkType = "signed hybrid",
    corFnc = cor_function,
    maxModuleSize = 10000,
    maxGoldModuleSize = 500,
    nPermutations = permutations,
    randomSeed = preservation_random_seed,
    verbose = 4
  )
}

cross_species_preservation <- lapply(
  names(cross_species_data),
  function(assay) {
    settings <- preservation_settings[[assay]]

    run_cross_species_preservation(
      prepared_data = cross_species_data[[assay]],
      cor_function = settings$cor_function,
      permutations = settings$permutations
    )
  }
)

names(cross_species_preservation) <- names(cross_species_data)
```

## Run within-reference module-quality checks

For each network, the same molecular matrix is supplied as both the reference and comparison set. This reproduces the dedicated quality-check design in the supplied scripts while avoiding any dependence on the number of features shared between species.

```r
run_internal_quality_check <- function(dataset_name) {
  dat_expr <- get_expression(dataset_name)
  colors <- align_colors_to_expression(
    dat_expr,
    get_module_colors(dataset_name)
  )

  duplicated_expression <- list(
    reference = list(data = dat_expr),
    duplicate = list(data = dat_expr)
  )

  reference_colors <- list(
    reference = colors
  )

  modulePreservation(
    multiData = duplicated_expression,
    multiColor = reference_colors,
    referenceNetworks = 1,
    networkType = "signed hybrid",
    corFnc = quality_settings$cor_function,
    maxModuleSize = 10000,
    maxGoldModuleSize = 500,
    nPermutations = quality_settings$permutations,
    randomSeed = preservation_random_seed,
    verbose = 4
  )
}

internal_quality <- lapply(
  required_datasets,
  run_internal_quality_check
)

names(internal_quality) <- required_datasets
```

## Extract summary quality and preservation statistics

For module quality, `Zsummary.qual` summarizes how strongly the module is defined in its reference network, while `medianRank.qual` ranks module quality relative to the other reference modules. For cross-species preservation, `Zsummary.pres` summarizes permutation-based preservation evidence, whereas `medianRank.pres` provides a relative preservation ranking that is less dependent on module size. Lower median rank indicates stronger relative quality or preservation.

```r
extract_summary_statistics <- function(
  observed_statistics,
  z_statistics,
  statistic_type = c("quality", "preservation")
) {
  statistic_type <- match.arg(statistic_type)

  common_modules <- intersect(
    rownames(observed_statistics),
    rownames(z_statistics)
  )

  observed_statistics <- observed_statistics[
    common_modules,
    ,
    drop = FALSE
  ]

  z_statistics <- z_statistics[
    common_modules,
    ,
    drop = FALSE
  ]

  if (statistic_type == "quality") {
    output <- data.frame(
      module = common_modules,
      analyzed_module_size = observed_statistics$moduleSize,
      median_rank_quality = observed_statistics$medianRank.qual,
      z_summary_quality = z_statistics$Zsummary.qual,
      stringsAsFactors = FALSE
    )
  } else {
    output <- data.frame(
      module = common_modules,
      analyzed_module_size = observed_statistics$moduleSize,
      median_rank_preservation = observed_statistics$medianRank.pres,
      z_summary_preservation = z_statistics$Zsummary.pres,
      z_density_preservation = z_statistics$Zdensity.pres,
      z_connectivity_preservation = z_statistics$Zconnectivity.pres,
      stringsAsFactors = FALSE
    )
  }

  output[!output$module %in% c("grey", "gold"), , drop = FALSE]
}
```

## Build the within-reference quality table


```r
extract_internal_quality <- function(result, dataset_name) {
  observed <- result$quality$observed[[1]][[2]]
  z_scores <- result$quality$Z[[1]][[2]]

  quality_table <- extract_summary_statistics(
    observed_statistics = observed,
    z_statistics = z_scores,
    statistic_type = "quality"
  )

  quality_table$dataset <- dataset_name
  quality_table$quality_above_2 <-
    quality_table$z_summary_quality > 2

  quality_table
}

module_quality_tables <- lapply(
  names(internal_quality),
  function(dataset_name) {
    extract_internal_quality(
      result = internal_quality[[dataset_name]],
      dataset_name = dataset_name
    )
  }
)

module_quality_table <- do.call(
  rbind,
  module_quality_tables
)

rownames(module_quality_table) <- NULL

quality_overview <- do.call(
  rbind,
  lapply(
    split(module_quality_table, module_quality_table$dataset),
    function(dataset_table) {
      data.frame(
        dataset = unique(dataset_table$dataset),
        number_of_modules = nrow(dataset_table),
        minimum_z_summary_quality = min(
          dataset_table$z_summary_quality,
          na.rm = TRUE
        ),
        all_modules_above_2 = all(
          dataset_table$quality_above_2,
          na.rm = TRUE
        )
      )
    }
  )
)
```

## Annotate modules by developmental-stage association

The preservation figures in the supplied scripts focused on modules showing a statistically supported positive association with at least one developmental stage. The annotation below uses the complete module–trait results from Part 3 and assigns each supported module to the stage with its strongest supported positive correlation.

```r
get_stage_annotations <- function(dataset_name) {
  module_trait <- wgcna_results$module_trait_results[[dataset_name]]

  correlations <- module_trait$correlation
  fdr <- module_trait$fdr

  modules <- sub("^ME", "", rownames(correlations))

  annotations <- lapply(seq_along(modules), function(i) {
    supported <-
      fdr[i, ] < 0.01 &
      correlations[i, ] > 0

    supported[is.na(supported)] <- FALSE

    if (any(supported, na.rm = TRUE)) {
      supported_correlations <- correlations[i, ]
      supported_correlations[!supported] <- -Inf

      stage <- names(which.max(supported_correlations))
      stage_correlation <- max(
        supported_correlations,
        na.rm = TRUE
      )
    } else {
      stage <- NA_character_
      stage_correlation <- NA_real_
    }

    data.frame(
      module = modules[i],
      stage_associated = any(supported, na.rm = TRUE),
      strongest_supported_stage = stage,
      stage_correlation = stage_correlation,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, annotations)
}

stage_annotations <- lapply(
  required_datasets,
  get_stage_annotations
)

names(stage_annotations) <- required_datasets
```

## Build complete bidirectional preservation tables

The preservation output is extracted in both directions. `total_reference_module_size` refers to the complete reference module, whereas `analyzed_module_size` is the number of its features present in both species and therefore included in the preservation calculation.

```r
extract_directional_preservation <- function(
  result,
  assay,
  reference_index,
  test_index,
  reference_species,
  test_species,
  shared_features
) {
  observed <- result$preservation$observed[[reference_index]][[test_index]]

  z_scores <- result$preservation$Z[[reference_index]][[test_index]]

  preservation_table <- extract_summary_statistics(
    observed_statistics = observed,
    z_statistics = z_scores,
    statistic_type = "preservation"
  )

  dataset_name <- paste0(reference_species, "_", assay)
  colors <- get_module_colors(dataset_name)

  total_sizes <- table(colors)
  shared_colors <- colors[
    intersect(names(colors), shared_features)
  ]
  shared_sizes <- table(shared_colors)

  preservation_table$total_reference_module_size <- as.integer(
    total_sizes[preservation_table$module]
  )

  preservation_table$shared_features_in_module <- as.integer(
    shared_sizes[preservation_table$module]
  )

  preservation_table$shared_features_in_module[
    is.na(preservation_table$shared_features_in_module)
  ] <- 0

  preservation_table$percent_reference_module_shared <-
    100 * preservation_table$shared_features_in_module /
    preservation_table$total_reference_module_size

  preservation_table$assay <- assay
  preservation_table$reference_species <- reference_species
  preservation_table$test_species <- test_species

  preservation_table <- merge(
    preservation_table,
    stage_annotations[[dataset_name]],
    by = "module",
    all.x = TRUE,
    sort = FALSE
  )

  preservation_table
}

build_assay_preservation_table <- function(assay) {
  result <- cross_species_preservation[[assay]]
  shared_features <- cross_species_data[[assay]]$shared_features

  erato_reference <- extract_directional_preservation(
    result = result,
    assay = assay,
    reference_index = 1,
    test_index = 2,
    reference_species = "erato",
    test_species = "melpomene",
    shared_features = shared_features
  )

  melpomene_reference <- extract_directional_preservation(
    result = result,
    assay = assay,
    reference_index = 2,
    test_index = 1,
    reference_species = "melpomene",
    test_species = "erato",
    shared_features = shared_features
  )

  rbind(erato_reference, melpomene_reference)
}

preservation_tables <- lapply(
  c("RNA", "ATAC"),
  build_assay_preservation_table
)

names(preservation_tables) <- c("RNA", "ATAC")

complete_preservation_table <- do.call(
  rbind,
  preservation_tables
)

rownames(complete_preservation_table) <- NULL
```

## Plot median-rank preservation for stage-associated modules

Median rank was the primary preservation measure used in the manuscript because it is a relative ranking and is less affected by differences in module size than permutation-derived Z statistics. Lower values indicate stronger relative preservation.

```r
plot_median_rank_preservation <- function(
  preservation_table,
  assay_name
) {
  plot_data <- preservation_table[
    preservation_table$stage_associated %in% TRUE,
    ,
    drop = FALSE
  ]

  plot_data$direction <- paste0(
    plot_data$reference_species,
    " modules in ",
    plot_data$test_species
  )

  ggplot(
    plot_data,
    aes(
      x = analyzed_module_size,
      y = median_rank_preservation,
      color = module,
      shape = strongest_supported_stage
    )
  ) +
    geom_point(size = 3) +
    facet_wrap(~ direction, scales = "free_x") +
    scale_color_identity() +
    theme_classic() +
    labs(
      title = paste(assay_name, "module preservation"),
      x = "Module features included in comparison",
      y = "Median-rank preservation",
      shape = "Developmental stage"
    )
}

rna_preservation_plot <- plot_median_rank_preservation(
  preservation_tables$RNA,
  "RNA-seq"
)

atac_preservation_plot <- plot_median_rank_preservation(
  preservation_tables$ATAC,
  "ATAC-seq"
)
```

## Save preservation objects, tables, and figures

```r
output_directory <- file.path(
  "wgcna_results",
  "module_preservation"
)

dir.create(
  output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

write.table(
  module_quality_table,
  file = file.path(
    output_directory,
    "within_reference_module_quality.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  quality_overview,
  file = file.path(
    output_directory,
    "within_reference_quality_overview.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  complete_preservation_table,
  file = file.path(
    output_directory,
    "cross_species_module_preservation.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

ggsave(
  filename = file.path(
    output_directory,
    "RNA_median_rank_preservation.pdf"
  ),
  plot = rna_preservation_plot,
  width = 9,
  height = 5
)

ggsave(
  filename = file.path(
    output_directory,
    "ATAC_median_rank_preservation.pdf"
  ),
  plot = atac_preservation_plot,
  width = 9,
  height = 5
)

saveRDS(
  list(
    settings = list(
      preservation = preservation_settings,
      quality = quality_settings,
      random_seed = preservation_random_seed
    ),
    shared_features = lapply(
      cross_species_data,
      `[[`,
      "shared_features"
    ),
    internal_quality = internal_quality,
    cross_species_preservation = cross_species_preservation,
    module_quality_table = module_quality_table,
    quality_overview = quality_overview,
    preservation_table = complete_preservation_table
  ),
  file = file.path(
    output_directory,
    "module_preservation_results.rds"
  )
)
```

The file `within_reference_module_quality.tsv` contains the per-module quality statistics that can be added to the supplementary module table. Before stating that all modules had `Zsummary.qual > 2`, confirm that `all_modules_above_2` is `TRUE` for all four datasets in `within_reference_quality_overview.tsv`.
