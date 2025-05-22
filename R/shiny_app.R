# This is a Shiny application.
#
# To run this app:
# 1. Open this file (R/shiny_app.R) in RStudio and click the "Run App" button.
# OR
# 2. In the R console, navigate to the project root directory and run:
#    shiny::runApp("R/shiny_app.R")
#
# This app uses placeholder data for demonstration purposes.

library(shiny)
library(ggplot2)
library(ggrepel)
library(plotly)

# Source necessary functions - assuming not in a package structure for now
source("R/plot-deseq.R")
source("R/interactive_plots.R")

ui <- fluidPage(
  titlePanel("Volcano Plot"),
  sidebarLayout(
    sidebarPanel(
      numericInput("padj_cutoff", "Adjusted P-value Cutoff", value = 0.05, min = 0, max = 1, step = 0.01),
      numericInput("label_n", "Number of Genes to Label", value = 10, min = 0),
      downloadButton("download_pdf", "Download PDF"),
      downloadButton("download_png", "Download PNG")
    ),
    mainPanel(
      plotlyOutput("volcano_plot_interactive")
    )
  )
)

server <- function(input, output, session) {

  # Placeholder res_df data frame
  res_df_placeholder <- data.frame(
    padj = runif(100, 0, 0.1),
    log2FoldChange = rnorm(100, 0, 2),
    gene_id = paste0("gene", 1:100),
    comparison = "ConditionA_vs_ConditionB"
  )

  # Reactive expression for volcano plot data
  reactive_volcano_plot_data <- reactive({
    req(input$padj_cutoff, input$label_n) # Ensure inputs are available
    make_volcano(
      res_df = res_df_placeholder,
      thresh = input$padj_cutoff,
      label_n = input$label_n
    )
  })

  # Render interactive plot
  output$volcano_plot_interactive <- renderPlotly({
    p <- reactive_volcano_plot_data()
    to_interactive_plot(p)
  })

  # Download PDF
  output$download_pdf <- downloadHandler(
    filename = function() {
      "volcano_plot.pdf"
    },
    content = function(file) {
      ggsave(file, plot = reactive_volcano_plot_data(), device = "pdf", width = 9, height = 9)
    }
  )

  # Download PNG
  output$download_png <- downloadHandler(
    filename = function() {
      "volcano_plot.png"
    },
    content = function(file) {
      ggsave(file, plot = reactive_volcano_plot_data(), device = "png", width = 9, height = 9)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
