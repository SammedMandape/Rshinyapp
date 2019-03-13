## Author: @Sammed Mandape, 2019
## PFI survival analysis code: @Shu Cheng, 2019

library(biomaRt)

#libraries for rshiny
library(shinydashboard)
library(ggplot2)
library(DT)
library(STRINGdb)
library(extrafont)

#libraries for survival analysis
library(dplyr)
library(tidyr)
library(survival)
library(survminer)

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
  ## Plot tab: Survival plot for input gene set
  #######################
  output$survivalplot <- renderPlot({
    
    #######################
    ##read clinical data
    #######################
    cli <- read.table("nationwidechildrens.org_clinical_follow_up_v1.0_prad.txt", as.is = TRUE,sep="\t",header=TRUE)
    cli <- cli[-c(1,2),]
    cli <- cli %>% dplyr::select("bcr_patient_barcode","new_tumor_event_dx_indicator","new_tumor_event_dx_days_to","new_tumor_event_type",
                                 "tumor_status","vital_status",
                                 "cause_of_death","last_contact_days_to","death_days_to")
    cli$stat <- ifelse(cli$new_tumor_event_dx_indicator == "YES" & 
                         !is.na(cli$new_tumor_event_dx_days_to) & cli$new_tumor_event_dx_days_to != "[Not Available]" |
                         cli$cause_of_death == "Prostate Cancer" & cli$new_tumor_event_dx_indicator != "YES", 1, 0)
    bb <- cli %>% filter(stat==1)
    cli$last_contact_days_to <- as.numeric(gsub('\\[Not Applicable\\]', 'NA',cli$last_contact_days_to))
    cli$death_days_to <- as.numeric(gsub('\\[Not Available\\]', 'NA', cli$death_days_to))
    cli$new_tumor_event_dx_days_to <- as.numeric(gsub('\\[Not Available\\]', 'NA', cli$new_tumor_event_dx_days_to))
    cli$Min <- with(cli, pmin(new_tumor_event_dx_days_to, death_days_to,na.rm=TRUE)) ##pick a smaller value
    cli$survtime <- ifelse(cli$stat == 1, cli[,"Min"], 
                           rowMeans(cli[, c("last_contact_days_to","death_days_to")], na.rm=TRUE))
    cli <-subset(cli,cli$last_contact_days_to !="NA" | cli$vital_status == "Dead") 
    clis1 <- cli %>% group_by(bcr_patient_barcode)
    clisort <- clis1 %>% arrange(new_tumor_event_dx_days_to, .by_group = TRUE) %>%
      arrange(desc(tumor_status), .by_group = TRUE) %>%
      arrange(!is.na(last_contact_days_to),desc(last_contact_days_to), .by_group = TRUE) %>%
      arrange(desc(stat), .by_group = TRUE)
    
    clisort <- clisort[!duplicated(clisort$bcr_patient_barcode),]
    clisort <- clisort %>% drop_na(survtime)
    clisort1 <- clisort[with(clisort, order(bcr_patient_barcode)),] %>% dplyr::select(bcr_patient_barcode,stat,survtime)
    
    #######################
    ##read a large table for gene expression FPKM:
    #######################
    
    tab3rows <- read.table("Prostate_FPKM", header = TRUE, nrows = 3)
    classes <- sapply(tab3rows, class)
    fpkm <- read.table("Prostate_FPKM", header = TRUE, row.names=1, colClasses = classes,comment.char = "")
    
    #######################
    ## Biomart id conversion Gene symbol (hgnc) -> Ensembl gene ids
    #######################
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    hgncgenes <- subset_dataset()$genes
    hgnc2ensembl <- getBM(filters = "hgnc_symbol",attributes = c("hgnc_symbol","ensembl_gene_id"),values = hgncgenes,mart = mart)
    colnames(hgnc2ensembl)<-c('SYMBOL','ENSEMBL')
    
    #browser()
    
    #######################
    ## read input gene list
    #######################
    fpkm1 <- fpkm %>% mutate(ENSEMBL=rownames(fpkm)) %>%
      #filter(rownames(fpkm) %in% as.vector(map.ensembl()$ENSEMBL))
      filter(rownames(fpkm) %in% as.vector(hgnc2ensembl$ENSEMBL))
    fpkm1 <- merge(hgnc2ensembl,fpkm1, by = "ENSEMBL") 
    
    ## Keep gene symbols only
    fpkm1 <- fpkm1[,-1]
    
    colnames(fpkm1) <- substr(colnames(fpkm1), 1, 16)
    fpkmT <- fpkm1[,grepl("01A$", colnames(fpkm1))] ##Tumor samples only
    rownames(fpkmT) <- fpkm1$SYMBOL
    colnames(fpkmT) <- gsub('\\.', '-', colnames(fpkmT)) ##same bar code as cli2 
    
    colnames(fpkmT) <- substr(colnames(fpkmT), 1, 12 ) ##1077
    fpkmT <- fpkmT[,!duplicated(colnames(fpkmT))]  ##1072
    
    fpkmTt <- t(fpkmT)
    fpkmTt <- data.frame(log2(fpkmTt+0.01)) %>% mutate(submitter_id = rownames(fpkmTt))
    
    #######################
    ## merge clinical and FPKM
    #######################
    clifpkm <- merge(clisort1, fpkmTt, by.x ="bcr_patient_barcode",by.y ="submitter_id", all=FALSE)  ##481 
    clifpkm_clu <- clifpkm[,-c(2,3)]  ##no clinical information
    
    #######################
    ## clustering here
    #######################
    clifpkm_clu[,-1] <- apply(clifpkm_clu[,-1], 2, as.numeric)  #convert to a numeric matrix
    summary(clifpkm_clu[,-1])
    #K-Means Cluster Analysis
    #library(fpc)
    fit1 <- kmeans(scale(clifpkm_clu[,-1]), 2) # 2 cluster solution
    clifpkm_clu1 <- clifpkm_clu %>% mutate(cluster1 = fit1$cluster)
    mydata <- left_join(clifpkm[,c("bcr_patient_barcode","survtime","stat")], clifpkm_clu1) 
    
    #######################
    ## output
    #######################
    #pdf("PRAD_PFI_KICSTOR.pdf",  width = 5, height = 4, onefile=FALSE)
    #browser()
    #surv_object <- Surv(time = mydata$survtime, event = mydata$stat)
    #browser()
    fitX <- survfit(Surv(time = survtime, event = stat) ~ cluster1, data = mydata)
    #browser()
    
    ggsurvplot(fitX, data = mydata, 
               pval = TRUE, pval.coord=c(3000,0.25), pval.size = 6, 
               legend.labs=c("Cluster1", "Cluster2"), legend.title="KICSTOR complex",legend = c(0.25,0.25),
               font.legend = c(14, "plain", "black"), palette=c("orchid2","dodgerblue2"),
               xlab="Days",ylab="Probability of Progression-free\n Survival")
    #return(survivalplot)
    #req(subset_dataset())
    #string_db <- STRINGdb$new()
    #mapped_genes <- string_db$map(subset_dataset(),"genes",removeUnmappedRows = TRUE)
    #hits_2<-mapped_genes$STRING_id
    #string_db$plot_network(hits_2)
    
  })
  
}

shinyApp(ui = ui, server = server)