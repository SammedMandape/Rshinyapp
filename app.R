## Author: @Sammed Mandape, 2019

library(shinydashboard)
library(ggplot2)
library(DT)
library(STRINGdb)
library(extrafont)

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Summary", tabName = "summary"),
    menuItem("Analysis", tabName = "analysis")
  )
)

body <- dashboardBody(
  tabItems(
    tabItem(tabName = "summary",
            fluidRow(
              box(
                title = "Studies", width = 12, status = "danger", solidHeader = TRUE,
                div(style = 'overflow-y: scroll; max-height: 600px', DT::dataTableOutput('table'))
              )
              
            ),
            hr(),
            
            fluidRow(
              box(
                title = "Type of cancer", status = "warning", width = 6, plotOutput("hist1")
              ),
              box(
                title = "Data type", status = "warning", color="yellow", width = 6, plotOutput("DataType"))
            )
    ),
  
    # Text Area Input for gene list
    tabItem(tabName = "analysis", 
            tabsetPanel(type = "tabs",
                        tabPanel("Plots",
                          fluidRow(
                            box(
                              title = "Input gene list here", status = "warning", width = 6, textAreaInput("title1","",value = "", width = '100%',rows=5),
                              actionButton("submit", "Submit")
                            )
                          ),
                          fluidRow(
                            box(
                              title = "Stirng DB network", status = "warning", width = 6, plotOutput("stringplot")
                            ),
                            box(
                              title = "Survival plot", status = "warning", width = 6, plotOutput("survivalplot")
                            )
                          )
                        ),
                        tabPanel("Heat map")
            )
    )
  )
)

# User interface --
ui <- dashboardPage(
  dashboardHeader(title = "Cancer research updates", titleWidth = 300),
  sidebar,
  body
)

# Server logic
server <- function(input, output) {
  #######################
  ## Summary sidebarpanel
  #######################
  pcdata <- read.csv("./data/ProstateCancerInput_Final.txt", sep = "\t", encoding = 'UTF-8')
  colnames(pcdata)[1] <- ''
  colnames(pcdata)[6] <- 'Data type'
  colnames(pcdata)[9] <- 'Type of cancer'
  pcdata$Link <- paste0("<a href ='",pcdata$Link,"'>","Link to data","</a>")
  output$table <- DT::renderDataTable({
    datatable(pcdata, rownames = FALSE, escape = FALSE)
  })

  CancerTypehist <- ggplot(data = pcdata, aes(`Type of cancer`)) +
                  geom_bar(stat = 'count', fill = 'peru') +
                  xlab("Cancer type") + ylab("# Cancer studies") +
                  theme(axis.text = element_text(family="sans", size = 12), axis.title = element_text(family="sans", face = "bold", size = 14))
  CancerTypehistFlip <- CancerTypehist + coord_flip()
  output$hist1 <- renderPlot({
    CancerTypehistFlip
  })
  
  
  DataTypehist <- ggplot(data = pcdata, aes(`Data type`, fill = Species)) + 
    geom_bar(stat = "count") +
    theme(legend.position = 'top') +
    scale_fill_manual(values=c('steelblue1','steelblue')) +
    xlab("Data type") + ylab("# Cancer studies") +
    theme(axis.text = element_text(family="sans", size = 12), axis.title = element_text(family="sans", face = "bold", size = 14),
          legend.text = element_text(size = 12, colour = "darkorange4"), legend.title = element_text(size = 12, face="bold"))
  DataTypehistFlip <- DataTypehist + coord_flip()
  output$DataType <- renderPlot({
    DataTypehistFlip
  })
  
  
  #######################
  ## Analysis sidebarPanel
  #######################
  
  #######################
  ## Plot tab: Stringdb plot for input gene set
  #######################
  subset_dataset <- eventReactive(input$submit,
                                  {data.frame(genes = unlist(strsplit(input$title1,split='[\r\n]')))
                                  })
  string_db <- STRINGdb$new(version="10", species=9606, input_directory="./")
  output$stringplot <- renderPlot({
    req(subset_dataset())
    #string_db <- STRINGdb$new()
    mapped_genes <- string_db$map(subset_dataset(),"genes",removeUnmappedRows = TRUE)
    hits_2<-mapped_genes$STRING_id
    string_db$plot_network(hits_2)})
  
  #######################
  ## Plot tab: Stringdb plot for input gene set
  #######################
}

shinyApp(ui = ui, server = server)