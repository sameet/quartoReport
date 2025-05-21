library(testthat)
library(ggplot2)
library(plotly) # Explicitly loading for class checks, though ggplotly loads it.
library(SummarizedExperiment)
library(quartoReport) # To make package functions available

# Tests for to_interactive_plot()
test_that("to_interactive_plot converts ggplot to plotly", {
  p_static <- ggplot(data.frame(x = 1:3, y = 1:3), aes(x, y)) + geom_point()
  p_interactive <- to_interactive_plot(p_static)
  expect_s3_class(p_interactive, "plotly")
  expect_s3_class(p_interactive, "htmlwidget")
})

# Tests for make_all_plots()
test_that("make_all_plots handles interactive argument correctly", {
  # Mock Data Creation
  meta_df <- data.frame(
    sample = paste0("sample", 1:4),
    condition = rep(c("A", "B"), each = 2),
    stringsAsFactors = FALSE
  )

  res_df <- data.frame(
    gene_id = paste0("gene", 1:10),
    baseMean = runif(10, 100, 1000),
    log2FoldChange = rnorm(10, 0, 1.5),
    lfcSE = runif(10, 0.1, 0.5),
    stat = rnorm(10, 0, 1),
    pvalue = runif(10, 0, 1),
    padj = runif(10, 0, 1),
    comparison = "A -- B", # Consistent with meta_df
    stringsAsFactors = FALSE
  )
  res_df$padj[1:3] <- 0.01 # Ensure some significant genes for plots

  counts_matrix <- matrix(rpois(40, lambda = 100), ncol = 4, nrow = 10)
  rownames(counts_matrix) <- paste0("gene", 1:10)
  colnames(counts_matrix) <- meta_df$sample
  col_data_for_se <- S4Vectors::DataFrame(meta_df) # Ensure column names match meta_df
  row_data_for_se <- S4Vectors::DataFrame(gene_id = rownames(counts_matrix))

  # Using unnamed list for assays as per the problem description notes
  vs_data_unnamed_assay <- SummarizedExperiment(
    assays = list(counts_matrix),
    colData = col_data_for_se,
    rowData = row_data_for_se
  )

  # Test interactive = TRUE
  # Suppress messages and warnings during plot generation for cleaner test output
  suppressMessages(suppressWarnings({
    plot_l_interactive <- make_all_plots(
      res_df = res_df,
      vs_data = vs_data_unnamed_assay,
      meta_df = meta_df,
      interactive = TRUE,
      label_n = 2,
      thresh = 0.05 # Explicitly pass thresh
    )
  }))
  expect_s3_class(plot_l_interactive$volcano, "plotly")
  expect_s3_class(plot_l_interactive$bp_s, "plotly")
  # pheatmap object itself or its gtable, or even a list containing the gtable
  expect_true(inherits(plot_l_interactive$hm, "pheatmap") || inherits(plot_l_interactive$hm, "gtable") || is.list(plot_l_interactive$hm))


  # Test interactive = FALSE
  suppressMessages(suppressWarnings({
    plot_l_static <- make_all_plots(
      res_df = res_df,
      vs_data = vs_data_unnamed_assay,
      meta_df = meta_df,
      interactive = FALSE,
      label_n = 2,
      thresh = 0.05 # Explicitly pass thresh
    )
  }))
  expect_s3_class(plot_l_static$volcano, "ggplot")
  expect_s3_class(plot_l_static$bp_s, "ggplot")
  expect_true(inherits(plot_l_static$hm, "pheatmap") || inherits(plot_l_static$hm, "gtable") || is.list(plot_l_static$hm))
})
