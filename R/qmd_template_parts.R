#' Create opening YAML for a report
#'
#' @param title Title for the report
#' @param author Author for the report
#' @param email Email of the author for the report
#' @param count_fn Path to the counts file
#' @param sample_fn Path to the sample-sheet file.
#' @param contrast_fn Path to the contrasts file.
#' @param metric_fn Path to the metric file.  Defaul value is NULL.
#' @param outputdir Path to the output directory where all the plots and outputs will be stored.  Will be created if it is not already present.
#' @param thresh Adjusted p-value threshold to call significant genes.
#'
#' @return opening_yaml string to be used to construct the report qmd
#' @export
#'
#' @examples
#' \dontrun{
#' opnening_yaml <- template_yaml()
#' }
template_yaml <- function(title = "RNA-Seq Report",
                          author = "Sameet",
                          email = "sameet.mehta@yale.edu",
                          count_fn = "gene_count_matrix.csv",
                          sample_fn = "sample-sheet.txt",
                          contrast_fn = "contrasts.txt",
                          metric_fn = NULL,
                          outputdir = "rnaseq_report_outputs",
                          thresh = 0.05
                          ) {

  if(is.null(metric_fn)) {
    metrics_fn_part <- stringr::str_glue("metrics: null")
  } else {
    metrics_fn_part <- stringr::str_glue("metrics: \"{metric_fn}\"")
  }

  stringr::str_glue("
---
title: \"{title}\"
author: \"{author}\"
email: \"{email}\"
date: \"`r Sys.Date()`\"
format:
  typst:
    keep-typ: true
  pdf:
    coloredlinks: true
  html:
    code-fold: true
    code-summary: \"Show Code\"
  ipynb:
    markdown-headings: atx
    prefer-html: true
execute:
  warning: false
  message: false
  echo: true
params:
  counts: \"{count_fn}\"
  meta_fn: \"{sample_fn}\"
  contrasts_fn: \"{contrast_fn}\"
  {metrics_fn_part}
  outputs: \"{outputdir}\"
  use_threshold: {thresh}
---

                    ")
}

#' Make setup part of the qmd file
#'
#' @return setup_qmd string that holds some text to setup the analysis.
#' @export
#'
#' @examples
#' \dontrun{
#' setup_qmd <- make_setup_qmd()
#' }
make_setup_qmd <- function() {
  setup_qmd <- stringr::str_glue("

## Introduction

This is analysis report for differential gene expression analysis.
The data is described in the later sections.
In this section we are setting up the requirements for the analysis to happen correcty.
We will import packages and functions.

```{{r}}
#| label: \"setup\"

# library(remotes)
# remotes::install_github(\"sameet/quartoReport\")

suppressPackageStartupMessages(library(dplyr)) # to make sure %>% is always found.
library(quartoReport)
# use_thresh <- 0.05 # global threshold to use as alpha in the analysis
if(!dir.exists(params$outputs)) dir.create(params$outputs)
```

Please note that this package is currently under active development.
This status may change in the future.
Currently this package is **not** ready for production use.
USE AT YOUR OWN RISK.

                                 ")
  # message(class(setup_qmd))
  setup_qmd
}

#' Create the part to plot the QC file
#'
#' @param fn name of the rnaseq metrics file if available.
#'
#' @return str A string that will correctly plot the QC of the rnaseq metrics.
#' @export
#'
#' @examples
#' \dontrun{
#' qc_bit <- make_qc_bit()
#' }
make_qc_bit <- function(fn = NULL) {
  if(is.null(fn) | fn == "null") {
    qc_bit <- stringr::str_glue("
## QC Report
We do not have a metrics file available so we will do not have a QC section.
                                ")
  } else {
    qc_bit <- stringr::str_glue("
## QC Report
Usually (for most standard model organisms) a metrics file is generated.
This file helpfully named as `rnaseq-metrics.txt` is part of the output folder.
The metrics file reports on a number of metrics about the alignment for each sample.
Here we are showing a handful of those parameters for a quick assessment on each of the samples.
These metrics are generally created using the `PICARD` suite of programs.

```{{r}}
#| label: tab-showQC
#| tab-cap: QC parameters and their values/interpretation.

metrics_df <- readr::read_delim(params$metrics, \"\\t\", show_col_types = FALSE)
metrics_df %>%
  dplyr::mutate(`% Unique` = gsub(\"%\", \"\", `% Unique`)) %>%
  dplyr::mutate(`% Unique` = as.numeric(`% Unique`)) -> metrics_df

metrics_df %>%
  dplyr::select(1, 3, 4, 5, 6, 7, 9, 10, 11, 12) %>%
  knitr::kable(booktabs = TRUE, linesep = \"\") %>%
  kableExtra::kable_paper(\"hover\", full_width = TRUE) %>%
  kableExtra::column_spec(1, color = ifelse(metrics_df$`% Unique` > 85, \"green\", \"red\")) %>%
  kableExtra::add_footnote(\"Samples passing filter are colored green.\")
```
                                ")
  }
  qc_bit
}

#' Create the analysis part of the text
#'
#' @return analysis_str analysis string that will put in place all the correct bits and pieces to actually perform the analysis.
#' @export
#'
#' @examples
#' \dontrun{
#' analysis_str <- make_analysis()
#' }
make_analysis <- function(){

  analysis1 <- stringr::str_glue("

## Analysis

The following sections detail how we ingest the data, prepare the data and process the data for analysis.

### Read in meta-data
Meta-data in this context the the description of the files.
A method to inform which samples bear which labels.
Essentially it is a dicionary that has a one -> many relation ship between conditions and samples.
But, it has one -> one relationshp between sample and condition.

This step expects a tab-separated values file with 2 columns.
We may extend this to more conditions per sample in the future, but currently we expect only 2 columns in this file.
In the code we will rename these columns to `sample`, and `condition`.

```{{r}}
#| label: tbl-read_meta_data
#| tbl-cap: The sample-sheet provided and being used in this analysis.

meta_df <- read_meta_data(params$meta_fn)

meta_df |>
  kableExtra::kbl(booktabs = TRUE) |>
  kableExtra::kable_styling(bootstrap_options = c(\"condensed\", \"striped\"),
                latex_options = c(\"striped\"))
```

### Read in Contrasts

We will read in the contrasts here.
These are the combinations that we are going to use.
If a `contrasts.txt` file is provided, we will do calculations for only those comparisons.
Otherwise, we will do all pair-wise comparisons from the `condition` column in the meta-data.

It is usually best practice to think through the analysis.
It is possible for a given analysis all conditions need to be compared against each other.
However, typically only a small sub-set of the comparisons will be relevant.
E.g. in a time series kind of experimental design the goal is to study differences between two consecutive time points serially.
If there are many time points number of comparisons increase exponentially.

Most experimental designs **do not have** hundreds of comparisons.

```{{r}}
#| label: tbl-contrasts
#| tbl-cap: The contrasts being used in this analysis.

contrasts_df <- get_contrasts(fn = params$contrasts, meta_df = meta_df)

contrasts_df |>
  kableExtra::kbl(booktabs = TRUE) |>
  kableExtra::kable_styling(bootstrap_options = c(\"condensed\", \"striped\"),
                latex_options = c(\"striped\"))
```

### Read in the count data

The RNA-Seq pipeline generates two count matrices, the transcript-level count matrix, and gene-level count matrix.
In this analysis we are going to use the gene-level count matrix *only*.
Consequently, for the smRNA-Seq analysis this will be either the \"mature\" smRNA counts, or the \"hairpin\" smRNA counts.

By design the genes are listed according to their Ensembl gene IDs, or other IDs.
For most model organisms most of the genes have a real human readable name, but for many non-standard model organisms not all genes have a human readable name.
We will be initially doing all the analysis using Ensembl gene IDs, and then annotate only the significantly changed genes for most common model organisms e.g. human, mouse, rat, drosophila etc.

```{{r}}
#| label: read_counts

count_mat <- read_count_mat(fn = params$counts, meta_df)
```

## Gene Expression Analysis

With the `DESeq2` method it is easier to generate the analysis for all possible combinations as defined by the `condition` column of the `meta_df`.
We can then extract the results for the comparions of interest rather easily.

```{{r}}
#| label: DEGA

dds <- make_deseq_object(count_mat = count_mat, meta_df = meta_df)
vs_data <- stabilize_count_data(dds)
```

It is advised that you save the `dds` object and the `vs_data` object.
There are relatively expensive objects and require sizable computation to generate.
We want to preserve them in case of further downstream analysis.

```{{r}}
#| label: saveObjects

save_objects(dds, op_dir = params$outputs)
save_objects(vs_data, op_dir = params$outputs)
```

The DESeqDataSet-... file represents the DESeq data-set that will be used in this analysis.
The DESEqTransform-... file represents a `vst`-transformed count data file.
The files will be saved in output directory defined in the `params` section of the `YAML`.

                                 ")
  analysis1
}

#' Make the bits to show graphs for comparisons.
#'
#' @param res_df Result data frame for single result.
#' @param thresh Adjusted p-value cutoff to call significant genes.
#' @param label_n Number of genes to label default is 30
#' @param comp_n Number of comparisons to annotated.
#'
#' @return A bit that will display graphs for one comparison at a time.
#' @export
#'
#' @examples
#' \dontrun{
#' make_single_comparison_bit(res_df)
#' }
make_single_comparison_bit <- function(res_df, thresh = 0.05, label_n = 30, comp_n) {
  res_df |>
    extract_comparisons() %>%
    gsub(" -- ", "--", .) -> comparison_name

  sig_df <- res_df |>
    dplyr::filter(padj <= thresh)

  sig_genes <- nrow(sig_df)

  if(nrow(sig_df) == 0) {
    single_comparison_bit <- stringr::str_glue("
## Results for comparison {comparison_name}

This comparison as 0 (zero) differentially expressted genes at adjusted $p$-value $\\lte$ {thresh}.
The numbers may change if a less stringent threshold is used for cutoff.
                                               ")
  } else {
    single_comparison_bit <- stringr::str_glue("

## Results for comparison {comparison_name}

This comparison has a total of {sig_genes} significant genes.
Of those in this report we will show informatin about top {label_n} genes.

### Changed genes for {comparison_name}

In @tbl-comparison_{comp_n} we see top genes changed between {comparison_name}.

```{{r}}
#| label: tbl-comparison_{comp_n}
#| tbl-cap: Top changed differentially expressed genes in {comparison_name}.

res_df_{comp_n} <- get_single_result(dds, comp_df = contrasts_df[{comp_n}, ])
res_df_{comp_n} |>
  dplyr::filter(padj <= {thresh}) |>
  dplyr::arrange(desc(abs(log2FoldChange))) |>
  dplyr::slice_max(order_by = abs(log2FoldChange), n = {label_n}) |>
  dplyr::arrange(desc(log2FoldChange)) |>
  dplyr::select(gene_id, log2FoldChange, padj) -> use_sig_df_{comp_n}

upreg_rows <- which(use_sig_df_{comp_n}$log2FoldChange > 0)
downreg_rows <- which(use_sig_df_{comp_n}$log2FoldChange < 0)

use_sig_df_{comp_n} |>
  kableExtra::kbl(booktabs = TRUE) |>
  kableExtra::kable_styling(bootstrap_options = c(\"condensed\"),
                           latex_options = c(\"striped\"),
                           font_size = 8
                           ) |>
  kableExtra::row_spec(upreg_rows, color = \"red\") |>
  kableExtra::row_spec(downreg_rows, color = \"blue\")
```

### Data Visualization for {comparison_name}

```{{r}}
#| label: make-plots_{comp_n}
#| echo: false

all_plot_{comp_n} <- make_all_plots(res_df = res_df_{comp_n},
                                    vs_data = vs_data,
                                    thresh = params$use_threshold,
                                    op_dir = params$outputs)

```

```{{r}}
#| label: fig-volcano_{comp_n}
#| fig-cap: Volcano Plot for {comparison_name}.  Genes in blue are down-regulated, and genes in red are up-regulated.

all_plot_{comp_n}$volcano
```

```{{r}}
#| label: fig-boxplot_single_{comp_n}
#| fig-cap: Boxplot for {comparison_name}.

all_plot_{comp_n}$bp_s
```
```{{r}}
#| label: fig-hm_{comp_n}
#| fig-cap: Heatmp for all significant genes for {comparison_name}.  There are total {sig_genes} significant genes.

ggplotify::as.ggplot(all_plot_{comp_n}$hm)
```


                                               ")
  }
  single_comparison_bit
}

#' Make the bit for overall analysis
#'
#' @param contrasts_df Data frame of contrasts.
#' @return str_glue object that will be printed in the report.
#'
#' @examples
#' \dontrun{
#'   scope_bit <- make_overall_scope_bit(contrasts_df)
#' }
make_overall_scope_bit <- function(contrasts_df) {
  # contrasts_df <- get_contrasts(fn = contrasts_fn)
  if(nrow(contrasts_df) == 1) {
    upset_bit <- stringr::str_glue("
### UpSet Plot

With only a single contrast there there are only two possible conditions the genes are either upregulated, or they are downregulated.
This information can be calculated from the results table already saved.
                                     ")
  } else {
    upset_bit <- stringr::str_glue("
### Upset Plot

Additionally when we have mutiple comparisons, real number of comparisons of interest are typically pattern that follow up/down regulation in multiple comparisons.
Effectively the classes of genes considering all the combinations of comparisons + up/down increase exponentially.
Usually, Venn Diagram is a way to represent this information.
However, most Venn Diagram plotting methods cannot go past 5 sets, and even at 5 sets the figure becomes extremely busy to the point of being un-interpretable.
To get past this problem there is new method of visualizaton called the UpSet plot.
This is a deconstructed Venn Diagram that is much more interpretable.
Upset plot for the data in this analysis is given in @fig-upset

```{{r}}
#| label: fig-upset
#| fig-cap: UpSet plot with with significant genes seen across all the comparisons.

upset_df <- make_upset_df(dds = dds, meta_df = meta_df,
                          contrasts_df = contrasts_df, thresh = params$use_threshold)
upset_plot <- make_upset_plot(upset_df)
ofn_upset <- file.path(params$outputs, \"upset-plot.pdf\")
ggplot2::ggsave(filename = ofn_upset, device = \"pdf\", plot = upset_plot,
                width = 9, height = 11, units = \"in\")
upset_plot
```
In @fig-upset the number in the bars denotes number of genes satisfying that condition.

                                   ")
  }
  overall_bit <- stringr::str_glue("
## Over All Analysis

In this analysis there were a total of `r nrow(contrasts_df)` comparisons between a total of `r nrow(meta_df)` samples.
In this section we will visualize how the data looks as a whole considering all the gene expression, and specifically significantly differentially expressed genes.

### PCA Plot

PCA plot is a potent dimensionality reduction method that can show the discriminative power of the detected significantly differentially expressed genes.

```{{r}}
#| label: fig-pca
#| fig-cap: PCA plots from the data. A) PCA plot using all the data, B) PCA plot using only significant genes, C) Variance captured by (A), and D) The amount of variation captured by each PC in (B).

pca_p_l <- make_pca_plot_l(dds = dds, vs_data = vs_data, contrasts_df = contrasts_df, thresh = params$use_threshold)
patchwork::wrap_plots(pca_p_l, nrow = 2, byrow = TRUE) +
  patchwork::plot_layout(guides = \"collect\") +
  patchwork::plot_annotation(tag_levels = \"A\") -> pca_plot

ofn_pca <- file.path(params$outputs, \"all-pca-plot.pdf\")
ggplot2::ggsave(filename = ofn_pca, device = \"pdf\", plot = pca_plot)
pca_plot
```
PCA plot for the data in this report is given in @fig-pca.

{upset_bit}
                                   ")
}
