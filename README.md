
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
This package will make most of the commonly anticipated visualizations e.g. Volcano Plots, Box Plots, Heatmaps, and Venn Diagrams.
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

```r
library(quartoReport)

## Preliminary processing

### Read in the data
meta_df <- read_meta_data("sample-sheet.txt")
count_mat <- read_count_mat("gene_count_matrix.csv")

### Make DESeq2 object for analysis
dds <- make_deseq_object(count_mat = count_mat, meta_df = meta_df) # will save the object
vs_data <- stabilize_count_data(dds = dds) # will save normalized expression.
```
