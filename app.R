options(repos = BiocManager::repositories())
getOption("repos")

library(shiny)
library(ggplot2)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(DT)
library(tidyverse)
library(enrichplot)
library(ReactomePA)
library(svglite)


ui <- fluidPage(
  # titlePanel("some text"),
  
  navbarPage(
    "KEGG pathway analysis and visualization",
    
    tabPanel(
      "Analysis",
      sidebarPanel(
        fileInput(inputId = "genelist", 
                  label = h3("Gene list input (csv format):"),
                  multiple = FALSE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")),
        
        # radioButtons(inputId = "species", 
        #              label = "Choose",
        #              choices = c("Human"= "hsa", "Mouse"="mmu")),
        
        selectInput(inputId = "species",
                    label = "Choose a species:",
                    choices = c("Human"= "hsa", "Mouse"="mmu")),
        br(),
        sliderInput(inputId = "pval", 
                    label = "KEGG p value cutoff",
                    min = 0, max = 1, 
                    value = 0.05),
        sliderInput(inputId = "qval", 
                    label = "KEGG q value cutoff",
                    min = 0, max = 1, 
                    value = 0.05),
        selectInput(inputId = "padjMethod",
                    label = "Choose a p adjustment method:",
                    choices = c("None"= "none", "Benjamini-Hochberg"="BH")),
        br(),
        selectInput(inputId = "x",
                    label = "Choose X axis:",
                    choices = c("GeneRatio"="GeneRatio", "Count"="Count")),
        selectInput(inputId = "color",
                    label = "Color by:",
                    choices = c("Adjusted P Value"="p.adjust",
                                "P Value"="pvalue", 
                                "q Value"="qvalue")),
        sliderInput(inputId = "show", 
                    label = "Number of pathways to show:",
                    min = 0, max = 100, 
                    value = 10),
        br(),
        selectInput(inputId = "format",
                    label = "Choose File format to save:",
                    choices = c("pdf"="pdf", "png"="png", "svg"="svg"))
      ),
      
      mainPanel(
        tabsetPanel(
          tabPanel("DotPlot", 
                   tags$h1("KEGG"),
                   plotOutput("keggdotplot"), 
                   downloadButton("keggdot_down","Download the plot"),
                   br(),
                   tags$h1("Reactome"),
                   plotOutput("reactomedotplot"),
                   downloadButton("reactomedot_down","Download the plot"),
                   br()), 
          tabPanel("Table", 
                   tags$h1("KEGG"),
                   dataTableOutput("kegg_dt"), 
                   downloadButton("keggdf", "Download the table"),
                   br(),
                   tags$h1("Reactome"),
                   dataTableOutput("reactome_dt"),
                   downloadButton("reactomedf", "Download the table"),
                   br(),
                   verbatimTextOutput("message_1"),
                   br()),
          tabPanel("BarPlot", 
                   tags$h1("KEGG"),
                   plotOutput("keggbarplot"), 
                   downloadButton("keggbar_down","Download the plot"),
                   br(),
                   tags$h1("Reactome"),
                   plotOutput("reactomebarplot"),
                   downloadButton("reactomebar_down","Download the plot"),
                   br()),
          tabPanel("Network", 
                   tags$h1("KEGG"),
                   plotOutput("keggcnetplot"), 
                   downloadButton("keggcnet_down","Download the plot"),
                   br(),
                   tags$h1("Reactome"),
                   plotOutput("reactomecnetplot"),
                   downloadButton("reactomecnet_down","Download the plot"),
                   br())
        )
      )
    ),
    
    tabPanel(
      "FAQ",
      includeMarkdown("manual.md")
    ),
    
    tabPanel(
      "Logs",
      textOutput("num"),
      includeMarkdown("updates.md")
    )
  )
  
  
)



server <- function(input, output, session) {
  
  df <- reactive({
    shiny::validate(
      need(input$genelist, "Upload a file!")
    )
    df <- read_csv(input$genelist$datapath, col_names = F)
    
    return(df)
  })
  
  org <- reactive({
    org <- ifelse(input$species=="hsa","org.Hs.eg.db","org.Mm.eg.db")
  })
  
  id <- reactive({
    id <- bitr(df()$X1, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org())
    return(id)
  })
  
  fc <- reactive({
    if(ncol(df())==2){
      fc <- df()$X2[match(id()$SYMBOL, df()$X1)]
      names(fc) <- id()$ENTREZID
    } else 
      fc <- NULL
    return(fc)
    
  })
  
  kegg <- reactive({
    kegg <- enrichKEGG(id()$ENTREZID, organism = input$species, 
                       pvalueCutoff = input$pval, 
                       pAdjustMethod = input$padjMethod,
                       minGSSize = 10, 
                       maxGSSize = 500, 
                       qvalueCutoff = input$qval) %>%
      setReadable(., OrgDb = org(), keyType = "ENTREZID")
    return(kegg)
  })
  
  reactome_sp <- reactive({
    reactome_sp <- ifelse(input$species=="hsa","human","mouse")
  })
  
  reactome <- reactive({
    reactome <- enrichPathway(id()$ENTREZID, organism = reactome_sp(), 
                              pvalueCutoff = input$pval, 
                              pAdjustMethod = input$padjMethod,
                              minGSSize = 10, 
                              maxGSSize = 500, 
                              qvalueCutoff = input$qval, 
                              readable=T)
    return(reactome)
  })
  
  # output$kegg_dt <- renderDataTable(fc() %>% data.frame())
  output$kegg_dt <- renderDataTable({
    kegg() %>% data.frame() %>% datatable() %>%
      formatRound(columns=5:7, digits=3)
  })
  
  output$reactome_dt <- renderDataTable({
    reactome() %>% data.frame() %>% datatable() %>%
      formatRound(columns=5:7, digits=3)
  })
  
  output$message_1 <- renderText({
    paste0(nrow(id()), " out of ", nrow(df()), " genes are mapped!")
  })
  
  output$keggdotplot <- renderPlot({
    dotplot(kegg(), x=input$x, color=input$color,
            showCategory=input$show)
  })
  
  output$reactomedotplot <- renderPlot({
    dotplot(reactome(), x=input$x, color=input$color,
            showCategory=input$show)
  })
  
  output$keggbarplot <- renderPlot({
    barplot(kegg(), x=input$x, color=input$color,
            showCategory=input$show) + ylab(input$x)
  })
  
  output$reactomebarplot <- renderPlot({
    barplot(reactome(), x=input$x, color=input$color,
            showCategory=input$show) + ylab(input$x)
  })
  
  output$keggcnetplot <- renderPlot({
    cnetplot(kegg(), # x=input$x, color=input$color,
             foldChange = fc(),
             showCategory=input$show)
  })
  
  output$reactomecnetplot <- renderPlot({
    cnetplot(reactome(), # x=input$x, color=input$color,
             foldChange = fc(),
             showCategory=input$show)
  })
  
  
  ## from here, create download files
  output$keggdot_down <- downloadHandler(
    ## specify file name
    filename = function(){
      paste0("keggdotplot.", input$format)
    },
    content = function(file){
      ggsave(filename = file, device = input$format,
             plot = dotplot(kegg(), x=input$x, 
                            color=input$color,
                            showCategory=input$show),
             width = 10, height = 5)
    })
  
  output$reactomedot_down <- downloadHandler(
    ## specify file name
    filename = function(){
      paste0("reactomedotplot.", input$format)
    },
    content = function(file){
      ggsave(filename = file, device = input$format,
             plot = dotplot(reactome(), x=input$x, 
                            color=input$color,
                            showCategory=input$show), 
             width = 10, height = 5)
    })
  
  output$keggbar_down <- downloadHandler(
    ## specify file name
    filename = function(){
      paste0("keggbarplot.", input$format)
    },
    content = function(file){
      ggsave(filename = file, device = input$format,
             plot = barplot(kegg(), x=input$x, 
                            color=input$color,
                            showCategory=input$show) + ylab(input$x),
             width = 10, height = 5)
    })
  
  output$reactomebar_down <- downloadHandler(
    ## specify file name
    filename = function(){
      paste0("reactomedbarplot.", input$format)
    },
    content = function(file){
      ggsave(filename = file, device = input$format,
             plot = barplot(reactome(), x=input$x, 
                            color=input$color,
                            showCategory=input$show) + ylab(input$x),
             width = 10, height = 5)
    })
  
  output$keggcnet_down <- downloadHandler(
    ## specify file name
    filename = function(){
      paste0("keggcnetplot.", input$format)
    },
    content = function(file){
      ggsave(filename = file, device = input$format,
             plot = cnetplot(kegg(), # x=input$x, color=input$color,
                             foldChange = fc(),
                             showCategory=input$show),
             width = 10, height = 5)
    })
  
  output$reactomecnet_down <- downloadHandler(
    ## specify file name
    filename = function(){
      paste0("reactomecnetplot.", input$format)
    },
    content = function(file){
      ggsave(filename = file, device = input$format,
             plot = cnetplot(reactome(), # x=input$x, color=input$color,
                             foldChange = fc(),
                             showCategory=input$show),
             width = 10, height = 5)
    })
  
  ## download pathway table
  output$keggdf <- downloadHandler(
    ## specify file name
    filename = "keggtable.csv",
    content = function(file){
      write.csv(kegg() %>% data.frame(), file)
    })
  
  output$reactomedf <- downloadHandler(
    ## specify file name
    filename = "reactometable.csv",
    content = function(file){
      write.csv(reactome() %>% data.frame(), file)
    })
  
  ## maybe create a zip file contain all?
  # output$downloadall <- downloadHandler(
  #   filename = "all.zip",
  #   content = function(file){
  #     
  #   }
  # )
  
  ## counters
  ## permission issue
  output$num <- renderText({
    if (!file.exists("data/counter.RData")) {
      counter <- 0
    } else
      load(file = "data/counter.RData")
    counter <- counter + 1
    save(counter, file = "data/counter.RData")
    paste0("You are the ", counter, " visitors.")
  })

}

shinyApp(ui, server)





