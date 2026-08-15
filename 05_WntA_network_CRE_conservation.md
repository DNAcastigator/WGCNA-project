# Part 5 — WntA network comparison and CRE conservation

This document describes extraction of the *WntA*-centered RNA networks, cross-species comparison of *WntA* correlations, and integration with the structural and temporal conservation of associated cis-regulatory elements (CREs).

For each gene, CRE homology and synchronicity are calculated as follows:

`CRE homology = number of homologous CREs / total number of associated CREs`

`CRE synchronicity = number of synchronous homologous CREs / number of homologous CREs`

Synchronicity is assigned `NA` when a gene has no homologous CREs because no homologous elements are available for temporal comparison.

## Analysis specifications

- The 99th-percentile threshold is calculated from the vector of TOM similarities between *WntA* and all other genes, not from the complete TOM.
- The *WntA* self-comparison is excluded before thresholding.
- Spearman p-values are adjusted once across the tested *WntA*–gene correlations within each species.
- Missing or unsupported correlations remain `NA`; they are not replaced with zero before classifying correlations as identical.
- Correlations are rounded only after statistical support has been determined.
- Peak homology is obtained from `WGCNA_peaks_homology.txt`.
- Each ATAC peak receives one module and therefore at most one developmental-stage assignment.
- Peaks with multiple homology-map rows are counted once in gene-level proportions.
- Pairwise Kolmogorov–Smirnov p-values are reported both raw and Holm-adjusted within each metric and species.

## Required R packages

```r
library(WGCNA)
library(Hmisc)
library(dplyr)
library(tidyr)
library(ggplot2)

options(stringsAsFactors = FALSE)
```

## Input files and analysis parameters

The Part 3 object contains the four molecular matrices, fitted networks, selected soft-thresholding powers, and complete module–trait statistics. The peak-homology table contains the *H. erato* peak identifier, *H. melpomene* peak identifier, and common pan-genome identifier used by the supplied scripts.

```r
wgcna_object_file <- file.path(
  "wgcna_results",
  "four_WGCNA_networks.rds"
)

peak_homology_file <- "WGCNA_peaks_homology.txt"

wnta_gene_candidates <- c("WntA", "WntA.t1", "Wnt5")
wnta_tom_quantile <- 0.99
wnta_correlation_alpha <- 0.01
wnta_p_adjust_method <- "bonferroni"
correlation_rounding_digits <- 2
number_of_overlap_permutations <- 10000
overlap_random_seed <- 12345
module_trait_alpha <- 0.01

wgcna_results <- readRDS(wgcna_object_file)
```



## Calculate RNA topological-overlap matrices

The species-specific powers selected in Part 3 are reused. Network type, TOM type, and correlation function are kept consistent with the fitted RNA networks.

```r
calculate_rna_tom <- function(dataset_name) {
  dat_expr <- wgcna_results$data[[dataset_name]]$expression
  power <- wgcna_results$parameters[dataset_name, "power"]

  tom <- TOMsimilarityFromExpr(
    datExpr = dat_expr,
    power = power,
    networkType = "signed hybrid",
    TOMType = "signed",
    corType = "bicor",
    verbose = 3
  )

  dimnames(tom) <- list(
    colnames(dat_expr),
    colnames(dat_expr)
  )

  tom
}

rna_tom <- list(
  erato = calculate_rna_tom("erato_RNA"),
  melpomene = calculate_rna_tom("melpomene_RNA")
)
```

## Extract the top 1% of WntA-associated genes

The threshold is calculated independently in each species from the distribution of TOM similarities between *WntA* and every other gene. Ties at the 99th percentile are retained, so a network can contain slightly more than exactly 1% of genes.

```r
extract_wnta_tom_network <- function(tom_matrix) {
  wnta_name <- resolve_wnta_name(colnames(tom_matrix))

  wnta_similarity <- tom_matrix[wnta_name, ]
  wnta_similarity <- wnta_similarity[
    names(wnta_similarity) != wnta_name
  ]

  cutoff <- quantile(
    wnta_similarity,
    probabilities = wnta_tom_quantile,
    na.rm = TRUE,
    names = FALSE
  )

  selected <- names(wnta_similarity)[
    !is.na(wnta_similarity) &
      wnta_similarity >= cutoff
  ]

  list(
    wnta_name = wnta_name,
    similarity = wnta_similarity,
    cutoff = cutoff,
    genes = unique(selected)
  )
}

wnta_tom_networks <- lapply(
  rna_tom,
  extract_wnta_tom_network
)

shared_wnta_genes <- intersect(
  wnta_tom_networks$erato$genes,
  wnta_tom_networks$melpomene$genes
)
```

## Test whether network overlap exceeds random expectation

Random networks are sampled without replacement from each species-specific RNA feature universe, using the observed network sizes. The empirical p-value is the proportion of random overlaps at least as large as the observed overlap, with a plus-one correction.

```r
set.seed(overlap_random_seed)

erato_universe <- names(
  wnta_tom_networks$erato$similarity
)

melpomene_universe <- names(
  wnta_tom_networks$melpomene$similarity
)

observed_overlap <- length(shared_wnta_genes)

random_overlap <- replicate(
  number_of_overlap_permutations,
  {
    random_erato <- sample(
      erato_universe,
      size = length(wnta_tom_networks$erato$genes),
      replace = FALSE
    )

    random_melpomene <- sample(
      melpomene_universe,
      size = length(wnta_tom_networks$melpomene$genes),
      replace = FALSE
    )

    length(intersect(random_erato, random_melpomene))
  }
)

overlap_z_score <- if (sd(random_overlap) > 0) {
  (observed_overlap - mean(random_overlap)) /
    sd(random_overlap)
} else {
  NA_real_
}

overlap_empirical_p <-
  (1 + sum(random_overlap >= observed_overlap)) /
  (number_of_overlap_permutations + 1)

overlap_summary <- data.frame(
  erato_network_size = length(wnta_tom_networks$erato$genes),
  melpomene_network_size = length(
    wnta_tom_networks$melpomene$genes
  ),
  observed_overlap = observed_overlap,
  random_mean = mean(random_overlap),
  random_standard_deviation = sd(random_overlap),
  z_score = overlap_z_score,
  empirical_p_value = overlap_empirical_p
)
```

## Recalculate WntA–gene correlations with Spearman correlation

Correlations are recalculated for genes found in the top 1% TOM network in both species. Spearman correlation is used because the inference depends on the monotonic correspondence between expression profiles and does not require normally distributed expression values or a strictly linear relationship. Bonferroni correction is applied across the tested *WntA*–gene correlations separately in each species.

```r
calculate_wnta_spearman <- function(
  dat_expr,
  genes,
  species_name
) {
  wnta_name <- resolve_wnta_name(colnames(dat_expr))

  genes <- intersect(
    genes,
    colnames(dat_expr)
  )

  correlation_input <- dat_expr[
    , c(wnta_name, genes),
    drop = FALSE
  ]

  correlation_result <- Hmisc::rcorr(
    as.matrix(correlation_input),
    type = "spearman"
  )

  rho <- correlation_result$r[wnta_name, genes]
  p_value <- correlation_result$P[wnta_name, genes]
  pairwise_n <- correlation_result$n[wnta_name, genes]

  adjusted_p_value <- p.adjust(
    p_value,
    method = wnta_p_adjust_method
  )

  data.frame(
    gene = genes,
    species = species_name,
    rho = as.numeric(rho),
    p_value = as.numeric(p_value),
    adjusted_p_value = as.numeric(adjusted_p_value),
    pairwise_n = as.numeric(pairwise_n),
    statistically_supported =
      adjusted_p_value < wnta_correlation_alpha,
    stringsAsFactors = FALSE
  )
}

erato_wnta_correlations <- calculate_wnta_spearman(
  dat_expr = wgcna_results$data$erato_RNA$expression,
  genes = shared_wnta_genes,
  species_name = "erato"
)

melpomene_wnta_correlations <- calculate_wnta_spearman(
  dat_expr = wgcna_results$data$melpomene_RNA$expression,
  genes = shared_wnta_genes,
  species_name = "melpomene"
)
```

## Classify conservation of WntA correlations

Only statistically supported correlations enter the Ic, Cc, and Ncc categories. Rounding is used solely to distinguish Ic from Cc after support has been established in both species.

- **Ic:** supported in both species and equal after rounding to two decimal places.
- **Cc:** supported in both species but different after rounding.
- **Ncc:** supported in only one species.
- **Unconfirmed:** unsupported in both species and excluded from the CRE-group tests.

```r
wnta_correlation_comparison <- erato_wnta_correlations %>%
  select(
    gene,
    erato_rho = rho,
    erato_p_value = p_value,
    erato_adjusted_p = adjusted_p_value,
    erato_supported = statistically_supported,
    erato_n = pairwise_n
  ) %>%
  full_join(
    melpomene_wnta_correlations %>%
      select(
        gene,
        melpomene_rho = rho,
        melpomene_p_value = p_value,
        melpomene_adjusted_p = adjusted_p_value,
        melpomene_supported = statistically_supported,
        melpomene_n = pairwise_n
      ),
    by = "gene"
  ) %>%
  mutate(
    erato_supported =
      !is.na(erato_supported) & erato_supported,
    melpomene_supported =
      !is.na(melpomene_supported) & melpomene_supported,
    erato_rho_rounded = round(
      erato_rho,
      correlation_rounding_digits
    ),
    melpomene_rho_rounded = round(
      melpomene_rho,
      correlation_rounding_digits
    ),
    correlation_class = case_when(
      erato_supported &
        melpomene_supported &
        erato_rho_rounded == melpomene_rho_rounded ~ "Ic",
      erato_supported & melpomene_supported ~ "Cc",
      xor(erato_supported, melpomene_supported) ~ "Ncc",
      TRUE ~ "Unconfirmed"
    )
  )

confirmed_wnta_correlations <- wnta_correlation_comparison %>%
  filter(correlation_class %in% c("Ic", "Cc", "Ncc"))
```

## Assign each ATAC peak to one developmental stage

Each peak inherits its module assignment from the fitted ATAC network. A module receives a stage only when it has a positive module–stage correlation with FDR below 0.01; if more than one stage is supported, the stage with the strongest positive correlation is retained. This prevents the same peak from being assigned to multiple developmental stages.

```r
get_module_stage_table <- function(dataset_name) {
  result <- wgcna_results$module_trait_results[[dataset_name]]
  correlations <- result$correlation
  fdr <- result$fdr

  modules <- sub("^ME", "", rownames(correlations))

  stage_rows <- lapply(seq_along(modules), function(i) {
    supported <-
      fdr[i, ] < module_trait_alpha &
      correlations[i, ] > 0

    supported[is.na(supported)] <- FALSE

    if (any(supported)) {
      supported_correlations <- correlations[i, ]
      supported_correlations[!supported] <- -Inf

      selected_stage <- names(
        which.max(supported_correlations)
      )

      selected_correlation <- max(
        supported_correlations
      )
    } else {
      selected_stage <- NA_character_
      selected_correlation <- NA_real_
    }

    data.frame(
      module = modules[i],
      stage = selected_stage,
      module_stage_correlation = selected_correlation,
      stage_associated = any(supported),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, stage_rows)
}

extract_gene_from_peak_id <- function(peak_id) {
  sub(
    "^[^_]+_[^_]+_[^_]+_",
    "",
    peak_id
  )
}

build_peak_stage_table <- function(dataset_name, species_name) {
  peak_colors <- wgcna_results$networks[[dataset_name]]$colors

  if (anyDuplicated(names(peak_colors))) {
    stop("Duplicated peak identifiers in ", dataset_name)
  }

  peak_table <- data.frame(
    peak = names(peak_colors),
    module = unname(peak_colors),
    stringsAsFactors = FALSE
  )

  peak_table <- peak_table %>%
    left_join(
      get_module_stage_table(dataset_name),
      by = "module"
    ) %>%
    mutate(
      gene = extract_gene_from_peak_id(peak),
      species = species_name
    )

  peak_table
}

peak_stage_tables <- list(
  erato = build_peak_stage_table(
    "erato_ATAC",
    "erato"
  ),
  melpomene = build_peak_stage_table(
    "melpomene_ATAC",
    "melpomene"
  )
)
```

## Read and resolve the peak-homology map

The mapping file is read using the three-column organization used in the supplied scripts: *H. erato* peak, *H. melpomene* peak, and pan-genome peak identifier. The resolver supports matrices named either with species-specific peak identifiers or with the shared pan-genome identifier.

```r
read_peak_homology <- function(file) {
  homology <- read.table(
    file,
    header = FALSE,
    sep = "\t",
    quote = "",
    stringsAsFactors = FALSE,
    fill = TRUE
  )

  if (ncol(homology) < 3) {
    stop("The peak-homology table must contain at least three columns")
  }

  homology <- homology[, 1:3, drop = FALSE]
  colnames(homology) <- c(
    "erato_original_peak",
    "melpomene_original_peak",
    "pan_peak"
  )

  homology <- homology %>%
    filter(
      !is.na(erato_original_peak),
      !is.na(melpomene_original_peak),
      !is.na(pan_peak)
    ) %>%
    distinct()

  homology
}

resolve_peak_homology <- function(
  homology,
  erato_peak_ids,
  melpomene_peak_ids
) {
  homology %>%
    mutate(
      erato_peak = case_when(
        pan_peak %in% erato_peak_ids ~ pan_peak,
        erato_original_peak %in% erato_peak_ids ~
          erato_original_peak,
        TRUE ~ NA_character_
      ),
      melpomene_peak = case_when(
        pan_peak %in% melpomene_peak_ids ~ pan_peak,
        melpomene_original_peak %in% melpomene_peak_ids ~
          melpomene_original_peak,
        TRUE ~ NA_character_
      ),
      homology_id = pan_peak
    ) %>%
    filter(!is.na(erato_peak) | !is.na(melpomene_peak)) %>%
    select(
      homology_id,
      erato_peak,
      melpomene_peak
    ) %>%
    distinct()
}

peak_homology_raw <- read_peak_homology(
  peak_homology_file
)

peak_homology <- resolve_peak_homology(
  homology = peak_homology_raw,
  erato_peak_ids = peak_stage_tables$erato$peak,
  melpomene_peak_ids = peak_stage_tables$melpomene$peak
)
```

A peak absent from this explicit homology map is reported below as non-homologous. 

## Classify CRE homology and temporal behavior

The analysis is performed once with each species as the reference. Only reference peaks assigned to a developmental stage are included in the total associated-CRE denominator, matching the stage-associated ATAC modules used in the downstream analysis.

For a homologous reference peak with multiple mapped counterparts, the peak is classified as synchronous when any counterpart has the same stage. It is otherwise asynchronous when at least one counterpart has a different assigned stage, or species-specific in timing when no counterpart has an assigned stage.

```r
classify_reference_cre <- function(
  reference_species,
  test_species,
  target_genes
) {
  reference_peaks <- peak_stage_tables[[reference_species]] %>%
    filter(
      stage_associated %in% TRUE,
      gene %in% target_genes
    )

  test_peak_stages <- peak_stage_tables[[test_species]] %>%
    select(
      counterpart_peak = peak,
      counterpart_stage = stage
    )

  if (reference_species == "erato") {
    directional_map <- peak_homology %>%
      transmute(
        reference_peak = erato_peak,
        counterpart_peak = melpomene_peak,
        homology_id = homology_id
      )
  } else {
    directional_map <- peak_homology %>%
      transmute(
        reference_peak = melpomene_peak,
        counterpart_peak = erato_peak,
        homology_id = homology_id
      )
  }

  mapped_peaks <- reference_peaks %>%
    rename(reference_peak = peak) %>%
    left_join(
      directional_map,
      by = "reference_peak"
    ) %>%
    left_join(
      test_peak_stages,
      by = "counterpart_peak"
    )

  peak_classification <- mapped_peaks %>%
    group_by(
      reference_peak,
      gene,
      module,
      reference_stage = stage
    ) %>%
    summarise(
      homologous = any(!is.na(homology_id)),
      synchronous =
        homologous &
        any(
          !is.na(counterpart_stage) &
            counterpart_stage == reference_stage
        ),
      counterpart_has_stage = any(
        !is.na(counterpart_stage)
      ),
      counterpart_stages = paste(
        sort(unique(na.omit(counterpart_stage))),
        collapse = ";"
      ),
      .groups = "drop"
    ) %>%
    mutate(
      temporal_class = case_when(
        !homologous ~ "non_homologous",
        synchronous ~ "synchronous",
        counterpart_has_stage ~ "asynchronous",
        TRUE ~ "species_specific_timing"
      ),
      reference_species = reference_species,
      test_species = test_species
    )

  gene_summary <- peak_classification %>%
    group_by(gene) %>%
    summarise(
      total_associated_cre = n_distinct(reference_peak),
      homologous_cre = n_distinct(
        reference_peak[homologous]
      ),
      synchronous_homologous_cre = n_distinct(
        reference_peak[temporal_class == "synchronous"]
      ),
      asynchronous_homologous_cre = n_distinct(
        reference_peak[temporal_class == "asynchronous"]
      ),
      species_specific_timing_cre = n_distinct(
        reference_peak[
          temporal_class == "species_specific_timing"
        ]
      ),
      .groups = "drop"
    ) %>%
    mutate(
      homologous_cre_proportion =
        homologous_cre / total_associated_cre,
      synchronous_cre_proportion = ifelse(
        homologous_cre > 0,
        synchronous_homologous_cre / homologous_cre,
        NA_real_
      ),
      reference_species = reference_species,
      test_species = test_species
    )

  all_target_genes <- data.frame(
    gene = unique(target_genes),
    stringsAsFactors = FALSE
  )

  gene_summary <- all_target_genes %>%
    left_join(gene_summary, by = "gene") %>%
    mutate(
      total_associated_cre = replace_na(
        total_associated_cre,
        0L
      ),
      homologous_cre = replace_na(homologous_cre, 0L),
      synchronous_homologous_cre = replace_na(
        synchronous_homologous_cre,
        0L
      ),
      asynchronous_homologous_cre = replace_na(
        asynchronous_homologous_cre,
        0L
      ),
      species_specific_timing_cre = replace_na(
        species_specific_timing_cre,
        0L
      ),
      reference_species = reference_species,
      test_species = test_species
    )

  list(
    peak_classification = peak_classification,
    gene_summary = gene_summary
  )
}

target_wnta_genes <- confirmed_wnta_correlations$gene

erato_cre_results <- classify_reference_cre(
  reference_species = "erato",
  test_species = "melpomene",
  target_genes = target_wnta_genes
)

melpomene_cre_results <- classify_reference_cre(
  reference_species = "melpomene",
  test_species = "erato",
  target_genes = target_wnta_genes
)
```

## Integrate correlation and CRE-conservation classes

Genes without stage-associated CREs are retained in the output for transparency but excluded from distributional tests because their CRE proportions are undefined.

```r
cre_peak_classification <- bind_rows(
  erato_cre_results$peak_classification,
  melpomene_cre_results$peak_classification
)

gene_cre_conservation <- bind_rows(
  erato_cre_results$gene_summary,
  melpomene_cre_results$gene_summary
) %>%
  left_join(
    confirmed_wnta_correlations,
    by = "gene"
  ) %>%
  mutate(
    homologous_cre_proportion = ifelse(
      total_associated_cre > 0,
      homologous_cre / total_associated_cre,
      NA_real_
    ),
    synchronous_cre_proportion = ifelse(
      homologous_cre > 0,
      synchronous_homologous_cre / homologous_cre,
      NA_real_
    )
  )
```

## Compare Ic, Cc, and Ncc distributions

Pairwise two-sample Kolmogorov–Smirnov tests are performed separately for homologous-CRE and synchronous-CRE proportions in each reference species. Holm-adjusted p-values are included to account for the three pairwise comparisons made for each species and metric.

```r
run_pairwise_ks <- function(
  data,
  metric,
  species_name
) {
  groups <- c("Ic", "Cc", "Ncc")
  group_pairs <- combn(groups, 2, simplify = FALSE)

  tests <- lapply(group_pairs, function(group_pair) {
    x <- data[[metric]][
      data$correlation_class == group_pair[1]
    ]

    y <- data[[metric]][
      data$correlation_class == group_pair[2]
    ]

    x <- x[is.finite(x)]
    y <- y[is.finite(y)]

    if (length(x) == 0 || length(y) == 0) {
      statistic <- NA_real_
      p_value <- NA_real_
    } else {
      test <- suppressWarnings(
        ks.test(x, y, exact = FALSE)
      )

      statistic <- unname(test$statistic)
      p_value <- test$p.value
    }

    data.frame(
      reference_species = species_name,
      metric = metric,
      group_1 = group_pair[1],
      group_2 = group_pair[2],
      n_group_1 = length(x),
      n_group_2 = length(y),
      ks_statistic = statistic,
      p_value = p_value,
      stringsAsFactors = FALSE
    )
  })

  test_table <- bind_rows(tests)
  test_table$holm_adjusted_p <- p.adjust(
    test_table$p_value,
    method = "holm"
  )

  test_table
}

ks_test_results <- bind_rows(
  lapply(c("erato", "melpomene"), function(species_name) {
    species_data <- gene_cre_conservation %>%
      filter(reference_species == species_name)

    bind_rows(
      run_pairwise_ks(
        species_data,
        "homologous_cre_proportion",
        species_name
      ),
      run_pairwise_ks(
        species_data,
        "synchronous_cre_proportion",
        species_name
      )
    )
  })
)
```

## Plot CRE conservation by correlation class

```r
cre_plot_data <- gene_cre_conservation %>%
  filter(correlation_class %in% c("Ic", "Cc", "Ncc")) %>%
  pivot_longer(
    cols = c(
      homologous_cre_proportion,
      synchronous_cre_proportion
    ),
    names_to = "metric",
    values_to = "proportion"
  ) %>%
  mutate(
    correlation_class = factor(
      correlation_class,
      levels = c("Ic", "Cc", "Ncc")
    ),
    metric = recode(
      metric,
      homologous_cre_proportion = "Homologous CREs",
      synchronous_cre_proportion =
        "Synchronous homologous CREs"
    )
  )

cre_conservation_plot <- ggplot(
  cre_plot_data,
  aes(
    x = correlation_class,
    y = proportion,
    fill = correlation_class
  )
) +
  geom_boxplot(
    width = 0.65,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    alpha = 0.55,
    size = 1.5
  ) +
  facet_grid(
    metric ~ reference_species
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(
    x = "WntA correlation-conservation class",
    y = "CRE proportion"
  )
```

## Save the complete downstream analysis

```r
output_directory <- file.path(
  "wgcna_results",
  "WntA_CRE_analysis"
)

dir.create(
  output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

wnta_network_membership <- data.frame(
  gene = union(
    wnta_tom_networks$erato$genes,
    wnta_tom_networks$melpomene$genes
  ),
  stringsAsFactors = FALSE
) %>%
  mutate(
    erato_top_1_percent =
      gene %in% wnta_tom_networks$erato$genes,
    melpomene_top_1_percent =
      gene %in% wnta_tom_networks$melpomene$genes,
    shared_top_1_percent =
      erato_top_1_percent & melpomene_top_1_percent,
    erato_tom_similarity =
      wnta_tom_networks$erato$similarity[gene],
    melpomene_tom_similarity =
      wnta_tom_networks$melpomene$similarity[gene]
  )

write.table(
  wnta_network_membership,
  file = file.path(
    output_directory,
    "WntA_top1percent_TOM_networks.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  overlap_summary,
  file = file.path(
    output_directory,
    "WntA_network_overlap_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  data.frame(random_overlap = random_overlap),
  file = file.path(
    output_directory,
    "WntA_network_random_overlap_distribution.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  wnta_correlation_comparison,
  file = file.path(
    output_directory,
    "WntA_spearman_correlation_comparison.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  cre_peak_classification,
  file = file.path(
    output_directory,
    "CRE_homology_and_timing_by_peak.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  gene_cre_conservation,
  file = file.path(
    output_directory,
    "CRE_conservation_by_WntA_network_gene.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  ks_test_results,
  file = file.path(
    output_directory,
    "CRE_conservation_pairwise_KS_tests.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

ggsave(
  filename = file.path(
    output_directory,
    "CRE_conservation_by_correlation_class.pdf"
  ),
  plot = cre_conservation_plot,
  width = 9,
  height = 7
)

saveRDS(
  list(
    parameters = list(
      tom_quantile = wnta_tom_quantile,
      correlation_alpha = wnta_correlation_alpha,
      p_adjust_method = wnta_p_adjust_method,
      rounding_digits = correlation_rounding_digits,
      overlap_permutations =
        number_of_overlap_permutations,
      overlap_random_seed = overlap_random_seed
    ),
    wnta_tom_networks = wnta_tom_networks,
    overlap_summary = overlap_summary,
    random_overlap = random_overlap,
    correlation_comparison =
      wnta_correlation_comparison,
    peak_homology = peak_homology,
    peak_stage_tables = peak_stage_tables,
    cre_peak_classification =
      cre_peak_classification,
    gene_cre_conservation = gene_cre_conservation,
    ks_tests = ks_test_results
  ),
  file = file.path(
    output_directory,
    "WntA_CRE_analysis_results.rds"
  )
)
```
