% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deseq-related-functions.R
\name{get_single_result}
\alias{get_single_result}
\title{Get single result from combination df, and dds object}
\usage{
get_single_result(dds, comp_df)
}
\arguments{
\item{dds}{object of class DESeqDataSet}

\item{comp_df}{one row of a combination df}
}
\value{
res_df A data frame with result comparing the two conditions.
}
\description{
Get single result from combination df, and dds object
}
\examples{
\dontrun{
res_df <- get_single_result(dds, comp_df)
}
}
