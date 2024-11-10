options(digits = 3)

server <- function(input, output) {
  analysis_data <- reactiveVal(NULL)  # Store the analysis results

  observeEvent(input$run, {
    req(input$file1)  # Ensure a file is uploaded before proceeding

    # Read the uploaded file
    tryCatch({
      file_data <- read.csv(input$file1$datapath, row.names = 1)
      Meta_out <- Meta_matrix(file_data, species = input$species)%>%
        as.data.frame()  # Pass the selected species

      analysis_data(Meta_out)  # Save the analysis results in the reactive variable

      output$fileContents <- DT::renderDataTable(
       {DT::datatable(file_data%>%
           mutate(across(is.numeric, round, digits = 2)),
           rownames = TRUE,
           extensions = c('Buttons', 'ColReorder', 'FixedHeader', 'Scroller'))})

      output$analysisResults <- DT::renderDataTable(
        {DT::datatable(Meta_out%>%
                         mutate(across(is.numeric, round, digits = 2)),
                       rownames = TRUE,
                       extensions = c('Buttons', 'ColReorder', 'FixedHeader', 'Scroller'))})

      output$status <- renderText({
        paste("Analysis completed successfully for", input$species, "data.")
      })

    }, error = function(e) {
      output$status <- renderText({
        paste("An error occurred: ", e$message, ". Please check your file format and try again.")
      })
    })
  })

  # Download handler for analysis results
  output$downloadResults <- downloadHandler(
    filename = function() {
      paste("analysis_results_", input$species, "_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(analysis_data(), file, row.names = TRUE)
    }
  )

  output$mark <- renderUI({
    HTML(markdown::markdownToHTML(knit("./test.Rmd", quiet = TRUE)))
  })
}
