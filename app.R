## Rshiny app to perform data analysis and visualization
## Author: @Sammed Mandape, 2019
## smandape@email.arizona.edu
## PFI survival analysis code: @Shu Cheng, 2019
## Jan 2019

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

library(pheatmap)


library('reutils')
library('hash')
library(pracma)
library(stringr)


sidebar <- dashboardSidebar(
  sidebarMenu(id = "sidebarid", #id important for updateTabItems
    menuItem("Home", tabName = "home"
             #, icon = icon("home")
             ),
    menuItem("Data", tabName = "data"),
    
    menuItem("TCGA - Analysis", tabName = "analysis")
  )
)

body <- dashboardBody(
  tabItems(
    tabItem(tabName = "home",
            fluidRow(
              box(
                width = 12, tags$h2(tags$b("CDAT: Cancer Data Analysis Toolkit")), 
                tags$br(),
                tags$br(),
                tags$p("This tool can be used to perform sruvival analysis, plot biological network, and
                       do pathway enrichment analysis.")
              )
            )        
    ),
  
  
    # Text Area Input for gene list
    tabItem(tabName = "analysis", 
            tabsetPanel(type = "tabs",
                        tabPanel("Analysis",
                          fluidRow(
                            box(
                              title = "Input gene list here", status = "warning", width = 6, textAreaInput("title1","",value = "", width = '100%',rows=5),
                              actionButton("submit", "Submit")
                            )
                          ),
                          fluidRow(
                            # tabBox(
                            #   #title = "Plots",
                            #   width = 6,
                            #   tabPanel(title = "Heat map", status = "warning", width = 6, plotOutput("heatmap")),
                            #   tabPanel(title = "Progression Free Survival", status = "warning", width = 6, plotOutput("survivalplot")),
                            #   tabPanel(title = "String DB network", status = "warning", width = 6, plotOutput("stringplot"))
                            #   
                            # )
                            box(
                              title = "Stirng DB network", status = "warning", width = 4, plotOutput("stringplot")
                            ),
                            box(
                              title = "Survival plot", status = "warning", width = 4, plotOutput("survivalplot")
                            ),
                            box(
                              title = "Heat map", status = "warning", width = 4, plotOutput("heatmap")
                            )
                          )
                        )
                        #tabPanel("Heat map")
            )
    ),
    
    tabItem(tabName = "data",
            fluidRow(
              #column(width = 4,
              box(
                dateRangeInput('dateRange',
                               label = 'Date range input: yyyy/mm/dd',
                               start = Sys.Date() - 30, end = Sys.Date(),
                               format = "yyyy/mm/dd"
                               
                )),
              
              box(
                selectInput("speciesSelect","Species:", c("Homo sapiens" = "Homo sapiens",
                                                          "Mus musculus" = "Mus musculus",
                                                          "Quercus alba" = "Quercus alba",
                                                          "Gossypium hirsutum" = "Gossypium hirsutum",
                                                          "Rhodotorula" = "Rhodotorula",
                                                          "Pandalus borealis" = "Pandalus borealis",
                                                          "Gallus gallus" = "Gallus gallus",
                                                          "Bos taurus" = "Bos taurus",
                                                          "Ovis aries" = "Ovis aries",
                                                          "Megathura crenulata" = "Megathura crenulata",
                                                          "Bothrops moojeni" = "Bothrops moojeni",
                                                          "Frangula alnus" = "Frangula alnus",
                                                          "Senegalia senegal" = "Senegalia senegal",
                                                          "Pachymenia carnosa" = "Pachymenia carnosa",
                                                          "Dermestes lardarius" = "Dermestes lardarius",
                                                          "synthetic construct" = "synthetic construct",
                                                          "Rattus norvegicus" = "Rattus norvegicus",
                                                          "Human alphaherpesvirus 1" = "Human alphaherpesvirus 1",
                                                          "Human betaherpesvirus 5" = "Human betaherpesvirus 5",
                                                          "Human gammaherpesvirus 4" ="Human gammaherpesvirus 4",
                                                          "JC polyomavirus" = "JC polyomavirus",
                                                          "Human immunodeficiency virus 1" = "Human immunodeficiency virus 1",
                                                          "Human gammaherpesvirus 8" = "Human gammaherpesvirus 8",
                                                          "Human polyomavirus 1" = "Human polyomavirus 1",
                                                          "Macaca mulatta polyomavirus 1" = "Macaca mulatta polyomavirus 1",
                                                          "Mus spretus" = "Mus spretus"
                                                          
                ),multiple = TRUE)
                
              ),
              
              hr(),
              box(
                selectInput("DataTypeSelect","DataType:", c("Expression profiling" = "Expression profiling",
                                                            "Genome binding/occupancy profiling" = "Genome binding/occupancy profiling",
                                                            "Genome variation profiling" = "Genome variation profiling",
                                                            "Methylation profiling" = "Methylation profiling",
                                                            "Non coding rna profiling"  = "Non coding rna profiling",
                                                            "Protein profiling" = "Protein profiling",
                                                            "Snp genotyping by snp array" = "Snp genotyping by snp array",
                                                            "Third party reanalysis" = "Third party reanalysis",
                                                            "Other" = "Other"
                                                            
                                                            
                                                            
                ),multiple = TRUE)),
              
              
              
              
              box(textInput("summaryBox", "KeyWord For Study Details:", "")),
              
              
              box(
                title = "Studies", width = 12, status = "danger", solidHeader = TRUE,
                #DT::dataTableOutput("Table1"),
                div(style = 'overflow-y: scroll; max-height: 600px', dataTableOutput('table'))
                #div(style = 'overflow-y: scroll; max-height: 300px', tableOutput('table'))
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
            
            
            
    )
  )
)

title <- tags$a(href='https://cancercenter.arizona.edu/researchers/shared-resources/bioinformatics',tags$img(src="blocka-logo.png", height = '40', width = '40'), tags$b('CDAT', style="color:black"))
#title <- tags$img(src="blocka-logo.png", height = '40', width = '40')

# User interface --
ui <- dashboardPage(
  dashboardHeader(title = title,
                  tags$li(class = "dropdown", actionButton("home1","", icon=icon("home"))
                    )
                  # dropdownMenu(type = "messages",
                  #              messageItem(
                  #                from = "Sales Dept",
                  #                message = "Sales are steady this month."
                  #              ))
                  #,titleWidth = 300
                  ),
  sidebar,
  body
)

# Server logic
server <- function(input, output, session) {
  #######################
  ## Observe event to go back to home
  ######################
  observeEvent(input$home1, {
    updateTabItems(session, "sidebarid", "home")
  })
  
  
  #######################
  ## Data menuItem
  #######################
  dataFrame <-reactive({
    
    startDate <-format(input$dateRange[1],"%Y/%m/%d")
    
    
    endDate <- format(input$dateRange[2],"%Y/%m/%d")
    
    
    query <- paste0("prostate cancer and \"gse\"[Filter] and ",  startDate) 
    query <-paste0(query, "[PDAT]:")
    query <-paste0(query,endDate)
    query <-paste0(query,"[PDAT]")
    
    
    pmids <- esearch(query,"gds",usehistory = TRUE)
    
    pmids
    
    articles <-esummary(pmids)
    articles
    
    #get just the gse number
    Accession <- articles$xmlValue("//Accession")
    Accession <-grep('GSE',Accession,value =TRUE)
    
    
    Study <- articles$xmlValue("//title")
    
    date <-articles$xmlValue("//PDAT")
    date
    
    year <- substring(date,1,4)
    year
    
    
    
    #link?
    Link <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",Accession)
    Link
    
    
    #dataType?
    dataType <-articles$xmlValue("//gdsType")
    dataType
    
    samples <- articles$xmlValue("//n_samples")
    samples
    
    
    #more than one value 
    Species <-articles$xmlValue("//taxon")
    Species
    
    
    #TypeOfCancer
    CancerType <-articles$xmlValue("//CancerType")
    CancerType
    
    
    a <-data.frame(Accession,Study,year,Link,dataType,samples,Species,CancerType)
    a$samples = as.integer(as.character(a$samples))
    
    colnames(a)[1] <- ''
    colnames(a)[5] <- 'Data type'
    colnames(a)[8] <- 'Type of cancer'
    a$Link <- paste0("<a href ='",a$Link,"'>","Link to data","</a>")
    
    a
    
    
  })
  
  DatatableTemp <-reactive({
    
    h <-hash()
    h[["expression profiling by array"]] <- "Expression profiling"
    h[["expression profiling by genome tiling array"]] <- "Expression profiling"
    h[["expression profiling by high throughput sequencing"]] <-"Expression profiling"
    h[["expression profiling by mpss"]] <-"Expression profiling"
    h[["expression profiling by rt-pcr"]] <-"Expression profiling"
    h[["expression profiling by sage"]] <-"Expression profiling"
    h[["expression profiling by snp array"]] <-"Expression profiling"
    
    
    
    h[["genome binding/occupancy profiling by array"]] <-"Genome binding/occupancy profiling"
    h[["genome binding/occupancy profiling by genome tiling array"]] <-"Genome binding/occupancy profiling"
    h[["genome binding/occupancy profiling by high throughput sequencing"]] <-"Genome binding/occupancy profiling"
    h[["genome binding/occupancy profiling by snp array"]] <-"Genome binding/occupancy profiling"
    
    h[["genome variation profiling by array"]] <-"Genome variation profiling"
    h[["genome variation profiling by genome tiling array"]] <-"Genome variation profiling"
    h[["genome variation profiling by high throughput sequencing"]] <-"Genome variation profiling"
    h[["genome variation profiling by snp array"]] <-"Genome variation profiling"
    
    h[["methylation profiling by array"]] <-"Methylation profiling"
    h[["methylation profiling by genome tiling array"]] <-"Methylation profiling"
    h[["methylation profiling by high throughput sequencing"]] <-"Methylation profiling"
    h[["methylation profiling by snp array"]] <-"Methylation profiling"
    
    h[["non-coding rna profiling by array"]] <-"Non coding rna profiling"
    h[["non-coding rna profiling by genome tiling array"]] <-"Non coding rna profiling"
    h[["non-coding rna profiling by high throughput sequencing"]] <-"Non coding rna profiling"
    
    h[["other"]] <-"Other"
    h[["protein profiling by mass spec"]] <-"Protein profiling"
    h[["protein profiling by protein array"]] <-"Protein profiling"
    
    h[["snp genotyping by snp array"]] <-"Snp genotyping by snp array"
    h[["third party reanalysis"]] <-"Third party reanalysis"
    
    
    #Test whether the hash has this key
    #View(has.key('genome variation profiling by array',h))
    
    tableBefore <- dataFrame()
    names(tableBefore)[names(tableBefore) == "Data type"] <- "DataType"
    
    datatypeFrame <- tableBefore[,5, drop=FALSE]
    
    
    
    names(datatypeFrame)[names(datatypeFrame) == "Data type"] <- "DataType"
    
    datatypeFrame$DataType <- as.character(datatypeFrame$DataType)
    
    
    for (row in 1:nrow(datatypeFrame)) {
      
      listTemp <- strsplit(datatypeFrame[row,"DataType"], ';')
      tempValue <-""
      
      
      
      listTemp <-sapply(listTemp, tolower)
      
      a <- as.integer(0)
      for( i in listTemp){
        i <- str_trim(i, side = c("both", "left", "right"))
        
        if(isTRUE(has.key(i,h))){
          
          if(a >=1){
            if(grepl(h[[i]], tempValue)){
              
            }
            else{
              tempValue <- paste(tempValue, " AND ")
            }
          }
          
          if(grepl(h[[i]], tempValue)){
            
          }
          else{
            tempValue <- paste(tempValue, h[[i]])
            a <- sum(a,1)
          }
          
          
          
        }
        
      }
      
      datatypeFrame[row, "DataType"] <- tempValue
      
    }  
    
    
    
    tableBefore$DataType <- datatypeFrame$DataType
    
    
    temp <- tableBefore[1,]
    
    temp <- temp[-1,]
    
    
    ##adding Species filtes to the table below:
    
    for (row in 1:nrow(tableBefore)) {
      
      boolean <- 'TRUE'
      
      for(i in input$speciesSelect){
        i <- str_trim(i, side = c("both", "left", "right"))
        if(grepl(i,tableBefore[row, "Species"])){
          
        }
        else{
          boolean <- 'FALSE'
        }
        
      }
      
      if(boolean == 'TRUE'){
        temp[nrow(temp)+1, ] <- tableBefore[row,]
      }
    }
    
    
    
    
    #temp
    
    ##adding DataType filtes to the table below:
    temp2 <- temp[1,]
    
    temp2 <- temp2[-1,]
    
    
    for (row in 1:nrow(temp)) {
      
      boolean <- 'TRUE'
      
      for(i in input$DataTypeSelect){
        i <- str_trim(i, side = c("both", "left", "right"))
        if(grepl(i,temp[row, "DataType"])){
          
        }
        else{
          boolean <- 'FALSE'
        }
        
      }
      
      if(boolean == 'TRUE'){
        temp2[nrow(temp2)+1, ] <- temp[row,]
      }
    }
    
    
    temp2
    
    
    #adding Summary filtes to the table below
    
    temp3 <- temp2[1,]
    
    temp3 <- temp3[-1,]
    
    
    for (row in 1:nrow(temp2)) {
      
      i <- input$summaryBox
      i <- str_trim(i, side = c("both", "left", "right"))
      i <- tolower(i)
      if(grepl(i,tolower(temp2[row, "Study"]))){
        temp3[nrow(temp3)+1, ] <- temp2[row,]
      }
    }
    
    
    temp3;
    
    
  })
  
  
  output$table<-DT::renderDataTable({
    #test.table
    
    
    
    datatable(DatatableTemp(), rownames = FALSE, escape = FALSE,
              
              options = list(
                
                autoWidth = TRUE,
                columnDefs = list(list(width = '120px', targets = c(0,6))),
                scrollX = TRUE
              )
    )}
  )
  
  DataTypehist <- reactive({
    data = dataFrame()
    speciesFrame <- data[,7, drop=FALSE]
    #View(speciesFrame)
    
    speciesFrame$Species <- as.character(speciesFrame$Species)
    for (row in 1:nrow(speciesFrame)) {
      if(grepl(";",speciesFrame[row, "Species"] , fixed=TRUE)) {
        speciesFrame[row, "Species"] <- 'Combined'
      }
    }
    
    
    #View(speciesFrame)
    
    DataTypehist_Temp <-ggplot(data = speciesFrame, aes(x = `Species`,fill = Species))  + 
      geom_bar(stat = "count",fill = 'peru') +
      theme(legend.position = 'top') +
      
      geom_text(stat='count', aes(label=..count..), position = position_stack(vjust = 0.5),size=6)+
      
      xlab("Species") + ylab("# of Studies") +
      theme(axis.text = element_text(family="sans", size = 12), axis.title = element_text(family="sans", face = "bold", size = 14),
            legend.text = element_text(size = 12, colour = "darkorange4"), legend.title = element_text(size = 12, face="bold"))
    
    DataTypehistFlip <- DataTypehist_Temp + coord_flip()
    DataTypehistFlip
    
  })
  
  
  output$DataType <- renderPlot({
    DataTypehist()
  })
  
  
  
  #######################
  ## Analysis menuItem
  #######################
  
  #######################
  ## Plot tab: Stringdb plot for input gene set
  #######################
  # reactive event to plot survival and network when user hits submit button
  subset_dataset <- eventReactive(input$submit,
                                  {data.frame(genes = unlist(strsplit(input$title1,split='[\r\n]')))
                                  })
  
  string_db <- STRINGdb$new(version="10", species=9606, input_directory="./")
  output$stringplot <- renderPlot({
    req(subset_dataset())
    #took it from below and wrote it outside loop so that the stringdb local copy is downloaded only once. 
    #string_db <- STRINGdb$new()
    mapped_genes <- string_db$map(subset_dataset(),"genes",removeUnmappedRows = TRUE)
    hits_2<-mapped_genes$STRING_id
    string_db$plot_network(hits_2)
    })
  
  
  
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
  
  tab3rows <- read.table("Prostate_FPKM_Ensembl2GeneSymbol_unique", header = TRUE, nrows = 2)
  #browser()
  classes <- sapply(tab3rows, class)
  
  ## UNCOMMENT THIS - NEW
  fpkm <- read.table("Prostate_FPKM_Ensembl2GeneSymbol_unique", header = TRUE, row.names= NULL, colClasses = classes,comment.char = "")
  
  ## COMMENT THIS - OLD
  #fpkm <- read.table("Prostate_FPKM", header = TRUE, row.names=1, colClasses = classes,comment.char = "")
  
  #######################
  ## read input gene list
  #######################
  ## COMMENT THIS -OLD
  # fpkm1 <- fpkm %>% mutate(ENSEMBL=rownames(fpkm)) %>%
  #   filter(rownames(fpkm) %in% as.vector(hgnc2ensembl$ENSEMBL))
  # fpkm1 <- merge(hgnc2ensembl,fpkm1, by = "ENSEMBL")
  
  ## UNCOMMENT THIS - NEW
  mydata1 <- reactive({
    fpkm1 <- fpkm %>% filter(fpkm[,1] %in% as.vector(subset_dataset()$genes)) #after keeping just the tumor samples, the code should start at this line.
  #fpkm2 <- reactive({
    fpkmM <- merge(subset_dataset(),fpkm1, by.x = "genes", by.y = 'hgnc_symbol') # this will do same as above line when read the modifide input file.
    colnames(fpkmM) <- substr(colnames(fpkmM), 1, 16)
    
    fpkmT <- fpkmM[,grepl("01A$", colnames(fpkmM))] ##Tumor samples only
    rownames(fpkmT) <- fpkmM$genes #### UNCOMMENT THIS - NEW
    #rownames(fpkmT) <- fpkm1$SYMBOL ## ## COMMENT THIS -OLD
    colnames(fpkmT) <- gsub('\\.', '-', colnames(fpkmT)) ##same bar code as cli2 
    
    colnames(fpkmT) <- substr(colnames(fpkmT), 1, 12 ) ##1077
    fpkmT <- fpkmT[,!duplicated(colnames(fpkmT))]  ##1072
    
    fpkmTt <- t(fpkmT)
    fpkmTt <- data.frame(log2(fpkmTt+0.01)) %>% mutate(submitter_id = rownames(fpkmTt))
    
    #######################
    ## merge clinical and FPKM
    #######################
    ## with fpkmTt
    clifpkm <- merge(clisort1, fpkmTt, by.x ="bcr_patient_barcode",by.y ="submitter_id", all=FALSE)  ##481 
    #browser()
    clifpkm_clu <- clifpkm[,-c(2,3)]  ##no clinical information
    
    #######################
    ## clustering here
    #######################
    clifpkm_clu[,-1] <- apply(clifpkm_clu[,-1], 2, as.numeric)  #convert to a numeric matrix
    summary(clifpkm_clu[,-1])
    fit1 <- kmeans(scale(clifpkm_clu[,-1]), 2) # 2 cluster solution
    #browser()
    clifpkm_clu1 <- clifpkm_clu %>% mutate(cluster1 = fit1$cluster)
    mydata <- left_join(clifpkm[,c("bcr_patient_barcode","survtime","stat")], clifpkm_clu1) 
    mydata
  })
  
  
  #######################
  ## Plot tab: Survival plot for input gene set
  #######################
  output$survivalplot <- renderPlot({
    fitX <- survfit(Surv(time = survtime, event = stat) ~ cluster1, data = mydata1())
    ggsurvplot(fitX, data = mydata1(), 
               pval = TRUE, pval.coord=c(3000,0.25), pval.size = 6, 
               legend.labs=c("Cluster1", "Cluster2"), 
               #legend.title="KICSTOR complex",legend = c(0.25,0.25),
               font.legend = c(14, "plain", "black"), palette=c("orchid2","dodgerblue2"),
               xlab="Days",ylab="Probability of Progression-free\n Survival")
    
    # sometimes there is need to print ggsurvplot explicitly and was trying out if return would work.
    #return(survivalplot)
    #req(subset_dataset())
  })
  
  #######################
  ## Plot tab: Heatmap
  #######################
  output$heatmap <- renderPlot({
    data_input <- mydata1()[,-c(1:3)]
    rownames(data_input) <- mydata1()[,1]
    data_inputSort <- data_input[order(data_input$cluster1),]
    data2matrix <- as.matrix(data_inputSort[,-ncol(data_inputSort)])
    anno_row <- data.frame("Clusters" = data_inputSort[,ncol(data_inputSort)])
    rownames(anno_row) <- rownames(data_inputSort) # check this again
    pheatmap(data2matrix, scale = "column", annotation_row = anno_row, cluster_row = F)
  })
  
}

shinyApp(ui = ui, server = server)