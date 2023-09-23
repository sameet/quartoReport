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

meta_df <- read_meta_data(params$samplesheet)

meta_df |>
  kbl(booktabs = TRUE) |>
  kable_styling(bootstrap_options = c(\"condensed\", \"striped\"),
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

Most experimental designs do not have hundreds of comparisons.

```{{r}}
#| label: tbl-contrasts
#| tbl-cap: The contrasts being used in this analysis.

contrasts_df <- get_contrasts(fn = params$contrasts, meta_df = meta_df)

contrasts_df |>
  kbl(booktabs = TRUE) |>
  kable_styling(bootstrap_options = c(\"condensed\", \"striped\"),
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

save_objects(dds)
save_objects(vs_data)
```
The `DESEqTransform-...` file represents a `vst`-transformed count data file.
The `DESeqDataSet-...` file represents the DESeq data-set that will be used in this analysis.

                                 ")
  analysis1
}
