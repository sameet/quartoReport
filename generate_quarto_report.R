#!/usr/bin/env Rscript

# Script to generate a Quarto markdown (.qmd) file for RNA-Seq analysis report
# using the quartoReport R package.
#
# Usage:
# Rscript generate_quarto_report.R --count_file <path_to_counts.csv> \
#                                --meta_file <path_to_metadata.tsv> \
#                                --output_qmd <output_report.qmd> \
#                                --report_output_dir <output_directory_name> \
#                                [--contrast_file <path_to_contrasts.tsv>] \
#                                [--metric_file <path_to_metrics.txt>] \
#                                [--report_title "My RNA-Seq Analysis"] \
#                                [--report_author "Your Name"] \
#                                [--report_email "your.email@example.com"] \
#                                [--p_value_threshold 0.05] \
#                                [--label_n_genes 30]
#
# Example:
# Rscript generate_quarto_report.R --count_file data/gene_count_matrix.csv.gz \
#                                --meta_file data/sample-sheet.txt \
#                                --output_qmd rnaseq_analysis.qmd \
#                                --report_output_dir rnaseq_analysis_outputs \
#                                --report_title "Differential Gene Expression Analysis"
#
# The script will generate a .qmd file. To produce the final report,
# you will need to render it using the Quarto CLI:
# quarto render <output_report.qmd>
#
# Ensure the 'quartoReport' package and its dependencies (DESeq2, optparse) are installed.
# You might need to install optparse separately if not already present:
# install.packages("optparse")
#
# If quartoReport is not installed globally, you might need to specify its library path
# or ensure your R environment is set up correctly (e.g., using renv with this project).

# Load necessary libraries
if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse")
}
library(optparse)

# Define command-line options
option_list <- list(
  make_option(c("-c", "--count_file"), type = "character", default = NULL,
              help = "Path to the gene count matrix CSV file (required).", metavar = "character"),
  make_option(c("-m", "--meta_file"), type = "character", default = NULL,
              help = "Path to the sample metadata TSV file (required).", metavar = "character"),
  make_option(c("-t", "--contrast_file"), type = "character", default = NULL,
              help = "Path to the contrasts definition file (optional). If not provided, all pairwise comparisons will be generated.", metavar = "character"),
  make_option(c("-e", "--metric_file"), type = "character", default = NULL,
              help = "Path to an RNA-Seq metrics file (optional).", metavar = "character"),
  make_option(c("-q", "--output_qmd"), type = "character", default = NULL,
              help = "Filepath for the output Quarto markdown (.qmd) file (required).", metavar = "character"),
  make_option(c("-d", "--report_output_dir"), type = "character", default = NULL,
              help = "Directory name for outputs from rendering the Quarto file (plots, tables, etc.) (required).", metavar = "character"),
  make_option(c("--report_title"), type = "character", default = "RNA-Seq Analysis Report",
              help = "Title for the generated report [default %default].", metavar = "character"),
  make_option(c("--report_author"), type = "character", default = "quartoReport User",
              help = "Author of the report [default %default].", metavar = "character"),
  make_option(c("--report_email"), type = "character", default = "",
              help = "Email of the author [default %default].", metavar = "character"),
  make_option(c("-p", "--p_value_threshold"), type = "numeric", default = 0.05,
              help = "Adjusted p-value threshold for significance [default %default].", metavar = "numeric"),
  make_option(c("-l", "--label_n_genes"), type = "numeric", default = 30,
              help = "Number of top genes to label in plots [default %default].", metavar = "numeric")
)

# Parse the options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$count_file) || is.null(opt$meta_file) || is.null(opt$output_qmd) || is.null(opt$report_output_dir)) {
  print_help(opt_parser)
  stop("Missing required arguments: --count_file, --meta_file, --output_qmd, and --report_output_dir must be specified.", call. = FALSE)
}

# Load necessary libraries for core functionality
# quartoReport should be installed, e.g. via devtools::install_github("sameet/quartoReport")
# or ensure it's available in the renv environment if used.
if (!requireNamespace("quartoReport", quietly = TRUE)) {
  stop("The 'quartoReport' package is not installed. Please install it to proceed.", call. = FALSE)
}
library(quartoReport)

# DESeq2 is a dependency of quartoReport and should be available.
# If not, it might indicate an issue with quartoReport installation.
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  stop("The 'DESeq2' package is not installed. It is a dependency of 'quartoReport'.", call. = FALSE)
}

message("All required libraries loaded.")

# Main script logic
message("Starting report generation process...")

# 1. Load metadata
message(paste("Loading metadata from:", opt$meta_file))
meta_df <- read_meta_data(opt$meta_file)
message("Metadata loaded successfully.")

# 2. Load count matrix
message(paste("Loading count matrix from:", opt$count_file))
# Assuming read_count_mat can handle .gz files if the path ends with .gz
# The function itself uses read.csv or read_delim which should handle paths directly.
count_mat <- read_count_mat(fn = opt$count_file, meta_df = meta_df)
message("Count matrix loaded successfully.")

# 3. Get contrasts
if (!is.null(opt$contrast_file)) {
  message(paste("Loading contrasts from:", opt$contrast_file))
  contrasts_df <- get_contrasts(fn = opt$contrast_file, meta_df = meta_df)
  message("Contrasts loaded from file.")
} else {
  message("Generating pairwise contrasts from metadata.")
  contrasts_df <- get_contrasts(meta_df = meta_df)
  message("Pairwise contrasts generated.")
}

# 4. Create DESeqDataSet object
message("Creating DESeqDataSet object (dds)...")
# This dds object will be in the global environment of this script,
# which is then accessible by compose_qmd as per its internal implementation.
dds <- make_deseq_object(count_mat = count_mat, meta_df = meta_df)
message("DESeqDataSet object (dds) created successfully.")

# 5. Call quartoReport::compose_qmd()
message(paste("Generating Quarto markdown file:", opt$output_qmd))

# The contrast_fn parameter for template_yaml (called by compose_qmd)
# should be the path to the contrast file if provided, otherwise NULL.
# compose_qmd itself takes contrasts_df.
# template_yaml needs the file path for the params section of the qmd.
actual_contrast_file_param <- if (!is.null(opt$contrast_file)) basename(opt$contrast_file) else NULL

# compose_qmd uses `...` to pass arguments to `template_yaml`.
# We need to ensure the file paths in `params` within the QMD are relative or just filenames
# if they are expected to be in the same directory as the QMD or a known path.
# For simplicity, we'll use basenames for files assuming they might be copied or accessible
# relative to where the qmd is run. The `report_output_dir` is a directory name.

# Ensure paths passed to params in YAML are what's intended for the QMD context.
# Typically, these would be relative paths from the QMD file.
# For this script, we assume the user provides paths accessible from the script's CWD.
# The QMD, when rendered, will need these files. It's often best if these are relative paths
# from the QMD's location, or the QMD is rendered with a CWD where these paths are valid.
# We will use basenames for file params and assume they are in the same dir as the QMD or in `report_output_dir`.
# The `outputdir` param for `template_yaml` is `opt$report_output_dir`.

quarto_doc_path <- compose_qmd(
  ofn = opt$output_qmd,
  contrasts_df = contrasts_df,
  thresh = opt$p_value_threshold,
  label_n = opt$label_n_genes,
  metric_fn = if (!is.null(opt$metric_file)) basename(opt$metric_file) else NULL, # Pass base filename for qmd param
  title = opt$report_title,
  author = opt$report_author,
  email = opt$report_email,
  count_fn = basename(opt$count_file),        # Pass base filename for qmd param
  sample_fn = basename(opt$meta_file),       # Pass base filename for qmd param
  contrast_fn = actual_contrast_file_param, # Pass base filename for qmd param (or NULL)
  outputdir = opt$report_output_dir      # This is a directory name
)

message(paste("Successfully generated Quarto markdown file:", quarto_doc_path))
message("To render the report, navigate to the directory containing the .qmd file and run:")
message(paste("  quarto render", basename(quarto_doc_path)))
message(paste("Ensure that the input files (", basename(opt$count_file), ",", basename(opt$meta_file),
              if(!is.null(actual_contrast_file_param)) paste(",", actual_contrast_file_param),
              if(!is.null(opt$metric_file)) paste(",", basename(opt$metric_file)),
              ") are accessible from the location of", basename(quarto_doc_path),
              "usually by placing them in the same directory or in a subdirectory referenced accordingly."))

message("Report generation script finished.")
