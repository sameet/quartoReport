
# quartoReport

<!-- badges: start -->
<!-- badges: end -->

The goal of `quartoReport` is to automatically generate a `Quarto` document with differential gene expression analysis using `DESeq2`.  Currently it is limited to straightup comparisons out of the box.
We are going to make this into a `Docker`-ized version so it can be used even when `R` itself is not available or not installed.  We are also going to make available prescription to create an `apptainer` image of this package to use on system where `Docker` is not allowed.

### Important Consideration

This is mostly going to be an internal facing package that is used for Yale (YCGA) - internal workflows.
I do hope that they are useful outside of that context as well.  

## Installation

You can install the development version of quartoReport from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("sameet/quartoReport")
```

## Expected workflow

Although the final aim of this package is to get an automated report, you can use this package without generating a report.
This package will make most of the commonly anticipated visualizations e.g. Volcano Plots, Box Plots, and Heatmaps.
A typical workflow is described below.
This assumes availability of the following files:

1. A sample-sheet -- This is expected to be a 2 column tab-separated-values file.
It is expected to have first column named `sample`, and second column named `condition`.
In principal it could have more columns but in the current version only first two columns will be used.

2. A counts file -- This is generally a comma-separated-values file. 
Each row contains genes, and each column contains count data for each sample.
This package currently has no plans to add gene annotation functionality.  
The gene identifiers in the data will be used for plotting.

3. A contrasts file -- This os optional.
If it is not available a 2 condition combination will be created from the sample-sheet's condition column.

Primary aim of this package is to be able to generate a report quickly.
To that end it works most efficiently if package is being using inside of a `quarto` document.
It can still do the analysis and generate graphics on its own, albeit currently it can generate only 4 most commonly used graphs in differential gene expression analysis.
In a `quarto` document `YAML` it can put in place-holder paths for the inputs and outputs.  
These can be changed appropriately to generate new report from the same data.
This `qmd` file can also be generated automatically from the given data.
Currently there are plans to emplement a command line interface (CLI) to expose this ability.

### Typical Workflow

This workflow is expected to be "interactive" in `R` console.

```r
library(quartoReport)

## Preliminary processing

### Read in the data
meta_df <- read_meta_data("sample-sheet.txt")
count_mat <- read_count_mat("gene_count_matrix.csv")
contrasts_df <- get_contrasts(fn = "contrasts.txt")

# if you do not have a contrasts.txt file you can generate all-against-all contrsts from the meta data as follows
# contrasts_df <- get_contrasts(meta_df = meta_df)

### Make DESeq2 object for analysis
dds <- make_deseq_object(count_mat = count_mat, meta_df = meta_df) # will save the object
vs_data <- stabilize_count_data(dds = dds) # will save normalized expression.

## Get results for a single comparison

### This gets results for a single comparison from the contrasts_df
res_f <- get_single_result(dds = dds, comp_df = contrasts_df[1, ])

### Save the result tables.
save_results(res_df, op_dir = "<path to output directory>", thresh = 0.05)

### Create plots for the result
all_plots_l <- make_all_plots(res_df = res_df, 
                              vs_data = vs_data,
                              label_n = 30,
                              thresh = 0.05
                              op_dir = "<path to output directory>")
### This function call will also save the graphs to the output directory.
```

At this point you should see 2 text files for `all-results`, and `sig-results` in the output directory.
There should also be a few `pdf` files with graphs in the output directory.

In a `quarto` based workflow there will be following lines after the front matter

```
  params:
    counts: "gene_count_matrix.csv"
    meta_fn: "sample-sheet.txt"
    contrasts_fn: null
    metrics: null
    outputs: "report_output_dir"
    use_threshold: 0.05
```

In side the `quarto` documents you we replace the file names by `params$counts`, `params$meta_fn` etc..

This allows us to use same template over and over without having to type a lot of code.

## How to generate the Quarto document

```r
dds <- readRDS("path/to/existing/dds/object.RDS")
contrasts_df <- get_contrasts(fn = "path/to/contrasts/file.txt")
meta_df <- read_meta_data("path/to/sample-sheet.txt")
compose_qmd(ofn = "path/to/output/directory/where/the/qmd_is_saved.qmd",
            contrasts_df = contrasts_df)
```
