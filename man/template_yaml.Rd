% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qmd_template_parts.R
\name{template_yaml}
\alias{template_yaml}
\title{Create opening YAML for a report}
\usage{
template_yaml(
  title = "RNA-Seq Report",
  author = "Sameet",
  email = "sameet.mehta@yale.edu",
  count_fn = "gene_count_matrix.csv",
  sample_fn = "sample-sheet.txt",
  contrast_fn = "contrasts.txt",
  metric_fn = NULL,
  outputdir = "rnaseq_report_outputs",
  thresh = 0.05
)
}
\arguments{
\item{title}{Title for the report}

\item{author}{Author for the report}

\item{email}{Email of the author for the report}

\item{count_fn}{Path to the counts file}

\item{sample_fn}{Path to the sample-sheet file.}

\item{contrast_fn}{Path to the contrasts file.}

\item{metric_fn}{Path to the metric file.  Defaul value is NULL.}

\item{outputdir}{Path to the output directory where all the plots and outputs will be stored.  Will be created if it is not already present.}

\item{thresh}{Adjusted p-value threshold to call significant genes.}
}
\value{
opening_yaml string to be used to construct the report qmd
}
\description{
Create opening YAML for a report
}
\examples{
\dontrun{
opnening_yaml <- template_yaml()
}
}
