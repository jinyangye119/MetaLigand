library(shiny)
library(shinythemes)
library(DT)
library(Seurat)

ui = fluidPage(
  theme = shinytheme(theme = "flatly"),
    navbarPage(
    "MetaLigand",
      tabPanel("Home",
              h2("MetaLigand - An R toolkit for analyzing non-peptide ligands"),  # Main title for the home page
              br(),
              p("MetaLigand is a powerful tool for analyzing ligand data.
                    Use the tabs above to upload your data, run analyses, and download results."),
              p("This tool supports multiple species, including human, mouse, and zebrafish.
                    Get started by navigating to the 'Analysis' tab."),
              p("The overall design of MetaLigand is shown below:"),
              img(src="Fig1.png", height="100%", width="100%", align = "center")
            ),
      tabPanel("Analysis",
             sidebarPanel(
                width=3,
                fileInput("file1", "Upload Your Data File",
                          multiple = FALSE,
                          accept = c(".csv", ".txt", ".xls", ".xlsx")),
                selectInput("species", "Select Species",
                            choices = c("human", "mouse", "zebrafish"),
                            selected = "mouse"),
                actionButton("run", "Run Analysis"),
                downloadButton("downloadResults", "Download Results")),
             mainPanel(
               h4("Uploaded File Preview"),
               div(style = "height: 30vh; overflow-y: auto;",
                   DT::dataTableOutput("fileContents")
                   )
            ),
            mainPanel(
              width = 12,
              h4("Analysis Results"),
              div(style = "height: 120vh; overflow-y: auto;",
                  DT::dataTableOutput("analysisResults")
              )
            ),
            mainPanel(
              verbatimTextOutput("status")
              )
            
    )
      )
)
