output$uiOutput_data <- renderUI({ 
elements <- list(
  tabPanel("Raw Data", dataTableOutput("table_rawData")),
  tabPanel("Processed Data", dataTableOutput("table_processedData"))
)

do.call(tabsetPanel, elements)
})