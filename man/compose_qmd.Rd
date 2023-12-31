% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compose_qmd.R
\name{compose_qmd}
\alias{compose_qmd}
\title{Compose the Quarto markdown file.}
\usage{
compose_qmd(
  ofn = NULL,
  contrasts_df,
  thresh = 0.05,
  label_n = 30,
  metric_fn = NULL,
  ...
)
}
\arguments{
\item{ofn}{Name of the output file.  A temporary file will be generated}

\item{contrasts_df}{The contrasts_df to determine number of contrasts in the data.}

\item{thresh}{The threshold for calculation of default 0.05}

\item{label_n}{Number of genes to label if possible default 30}

\item{metric_fn}{Name of the file with RNA-Seq alignment metrics.}

\item{...}{Other parameters to pass to opening yaml function}
}
\value{
ofn Path to the output file.
}
\description{
Compose the Quarto markdown file.
}
\examples{
\dontrun{
ofn <- compose_qmd()
}
}
