library(shiny)

source("prepareData.R")    
source("manhattan_gene.R")    

# Manhattan plot of SNVs with their corresonding (-log) p-values 
# associated with expression of the gene `r topGene`.

# R shiny app 

ui <- fluidPage(
            titlePanel("Manhattan plot"),
            fluidRow(column(6, sliderInput(inputId = 'gene_width',
                                           label = 'Choose width around gene',
                                           value = 25, 
                                           min = 1, 
                                           max = 10)),
                     column(6, selectizeInput(inputId = "gene_name", 
                                              choices = res$gene_name,
                                              selected = NULL,
                                              label = 'Choose gene'))),
            plotOutput('manhattan')
)

server <- function(input, output){
    
    output$manhattan <- renderPlot({plot_manhattan_gene(res, 
                                                        snps, 
                                                        input$gene_name, 
                                                        gtf, 
                                                        w = input$gene_width)})
    }

shinyApp(ui = ui, server = server)
