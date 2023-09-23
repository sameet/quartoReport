
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

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(quartoReport)
## basic example code
```

