####----Install and load packages----####

packages <- c("shiny","shinythemes","shinyjqui","stringi","rhoR","ggrepel",
              "dplyr","DT","ggplot2","ggpubr","plotly","ggVennDiagram",
              "readr","tidyr","stringr","shinycssloaders","reshape2")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
#bioconductor packages
bioCpacks <- c("clusterProfiler","GSVA","limma","enrichplot")
installed_packages_BIOC <- bioCpacks %in% rownames(installed.packages())
if (any(installed_packages_BIOC == FALSE)) {
  BiocManager::install(bioCpacks[!installed_packages_BIOC], ask = F)
}

invisible(lapply(packages, library, character.only = TRUE))
invisible(lapply(bioCpacks, library, character.only = TRUE))


####----Project Name----####

ProjName <- 'CCLE'


####----File Names----####

##--Database Files--##
#Meta
db_meta_file <- '~/R/DRPPM-EASY-Database-Integration-main/CCLE_data/CCLE_meta.zip'
#Meta Selector File
db_meta_selec_file <- '~/R/DRPPM-EASY-Database-Integration-main/CCLE_data/CCLE_meta_selector.tsv'
#Expression Data
db_expr_file <- '~/R/DRPPM-EASY-Database-Integration-main/CCLE_data/CCLE_expr.zip'
#Name Map File
db_namemap_file <- '~/R/DRPPM-EASY-Database-Integration-main/CCLE_data/CCLE_NameMap.tsv'




####
##
## Everything below should stay unchanged
##
####

##--Example Files--##
examp_file3.R <- "~/R/DRPPM-EASY-Database-Integration-main/Example_data/USP7_RNAseq_expr.txt"
examp_file3.RM <- "~/R/DRPPM-EASY-Database-Integration-main/Example_data/USP7_RNAseq_meta.tsv"

##--Gene Set Files--##
#MSigDB GS and Category table
MSigDBcat_file <- '~/R/DRPPM-EASY-Database-Integration-main/GeneSet_data/msigdb_gsNcat_HS.tsv'
#gene set list for ssGSEA
GSlist_file <- '~/R/DRPPM-EASY-Database-Integration-main/GeneSet_data/msigdb_gs_HS.RData'



####----Read and Reformat Files----####

##--Gene Set Files--##
#MSigDB gene sets FOR UI
MSigDBcat <- read.delim(MSigDBcat_file, header = T, sep = '\t')
#gene set list for ssGSEA
load(GSlist_file)
#subset category names
MSigDBcat2 <- MSigDBcat[,-3]
MSigDBcat3 <- unique(MSigDBcat2)
rownames(MSigDBcat3) <- 1:nrow(MSigDBcat3)
colnames(MSigDBcat3) <- c("MSigDB Categories","MSigDB Sub-Categories")


##--Database Files--##
#Read Meta File
db_meta <- as.data.frame(read_delim(db_meta_file, delim = '\t'))

#read meta selector file if it exists
if (file.exists(db_meta_selec_file) == T) {
  db_meta_selec <- as.data.frame(read_tsv(db_meta_selec_file))
}
if (file.exists(db_meta_selec_file) == F) {
  db_meta_selec <- NULL
}

#read expression
db_expr <- as.data.frame(read_delim(db_expr_file, delim = '\t'))
rownames(db_expr) <- db_expr[,1]
db_expr <- db_expr[,-1]
colnames(db_expr) <- gsub("[_.-]", ".", colnames(db_expr))

#read name map file if exists
if (file.exists(db_namemap_file) == TRUE) {
  
  db_namemap <- as.data.frame(read_tsv(db_namemap_file))
  
}
if (file.exists(db_namemap_file) == FALSE) {
  
  db_namemap <- NULL
  
}

#lineage/disease selection choices
subsetter_list <- list()
#if the selector file is given
if (is.null(db_meta_selec) == F) {
  
  #variables to subset data
  subsetters <- db_meta_selec[which(db_meta_selec[,2] == "Selector"),1]
  #split meta into selectors and phenotypes
  meta_subsetters <- db_meta[which(db_meta[,2] %in% subsetters),]
  meta_phenos <- db_meta[which(!db_meta[,2] %in% subsetters),]
  
  for (i in levels(factor(meta_subsetters[,2]))) {
    
    subsetter_list[[paste(i,"_Choices",sep = "")]] <- unique(meta_subsetters[which(meta_subsetters[,2] == i),3])
    
  }
  list2env(subsetter_list,globalenv())
  
  #phenotype selection choices
  db_pheno_choice <- unique(meta_phenos[,2])
  
}
#if the selector file is not given
if (is.null(db_meta_selec) == T) {
  
  #there are no choices to subset the expression
  subsetter_list[["No_Choices_Avaliable"]] <- NA
  
  #all choices in the meta are phenotypes
  meta_phenos <- db_meta
  
  #phenotype selection choices
  db_pheno_choice <- unique(meta_phenos[,2])
  
}


shinytheme("sandstone")

ui <-
  
  navbarPage(paste("{ ",ProjName," Integrative Expression Analysis }", sep = ""),
             
             ####----Intro Tab----####
             
             tabPanel("Intro and Methods",
                      fluidPage(
                        mainPanel(
                          tabsetPanel(
                            id = 'intropage',
                            tabPanel("Intoduction",
                                     h3("Introduction"),
                                     p("This is an extention of the DRPPM Expression Analysis ShinY (EASY) App which allows a user to integrate data from the Cancer Cell Line Encyclopedia (CCLE) with a dataset of their choice. With the gathered data they may perform differential gene expression comparison and reciprocal gene set enrichment analyses. This R Shiny app is very similar in features to the original EASY Integration app, expect for the addition of the Data Input and Sample Selection tab which allows for the user to upload their data and perform CCLE sample selection as well as preparation for comparison analysis downstream."),
                                     h3("Differential Gene Expression Methods"),
                                     p("Differential gene expression analysis is performed on the expression data between two groups defined by the provided metadata file. The samples chosen are log-transformed (log2 + 1). A model matrix is then designed for LIMMA linear modeling followed by empirical Bayes statistics to moderate gene variance and modeling of the global characteristic of the data. Of note, the current pipeline expects the matrix input quantification are pre-normalized, such as in CPM, TMM, or FPKM."),
                                     h3("Gene Set Enrichment Analysis (GSEA) Methods"),
                                     p("Gene Set Enrichment Analysis (GSEA) is performed through signal-to-noise calculations on the expression matrix. This calculation requires at least two types of sample groups and at least three samples for each grouping. The signal-to-noise calculation scales the difference of means by standard deviation and generates a ranked list of genes, effectively capturing the differences between phenotypes. The GSEA identifies the enrichment score (ES) for a gene set, which indicates the extent that the gene set is overrepresented at the beginning or end of the list of ranked genes. The enrichment score is calculated by walking down the ranked gene list and increasing the running-sum statistic when a gene is in the gene set and decreasing when it is not. The maximum deviation from zero that is encountered through walking the list is the ES, where a positive ES indicates the gene set is enriched at the top of the ranked list and a negative ES indicates the gene set is enriched at the bottom of the ranked list. A leading-edge gene list is provided displaying a subset of the genes in the gene set that contribute the most to the ES."),
                                     h4("Gene Set Sources"),
                                     p("The Molecular Signatures Database (MSigDB) is a collection of over 32,000 annotated gene sets divided into 9 major collections as well as various sub-collections [PMID: 16199517, PMID: 12808457]. The MSigDB  gene sets were downloaded with the msigdbr package and processed through R using Tidyverse packages to combine the gene sets into a data frame."),
                            ),
                            tabPanel("Data Selection Tutorial",
                                     h3("Download Example Expression and Meta Data for Comparison"),
                                     downloadButton("ExampDownload3.R", "Download Example RNAseq Expression Matrix"),
                                     downloadButton("ExampDownload3.RM", "Download Example RNAseq Meta Data"),
                                     h3("Subsetting Expression Data:"),
                                     p("This is an optional selection depending on the large project data that is being analyzed. If there is a 'MetaSelector' file, which is described on the GitHub page, this option to subset the expression data will appear. This is useful if the expression data can be subset be a certain parameter, for instance in the example of a CCLE data set the samples come from various different diseases or Cancer Linages, so the selector file and subsetting the expression data allows the user to narrow down which subset is of interest."),
                                     h3("Selecting a Coniditon of Interest"),
                                     p("After subsetting the expression data, if that is an option, the user must select a condition/phentotype from the drop down box which is generated based on the meta data."),
                                     h3("Other Options:"),
                                     p("The user may choose to log2 transform the selected expression data as well as name the selected data being analyzed which will be used as a label throughout the analysis."),
                                     h3("Grouping User and Project Data:"),
                                     p("Once the data is uploaded the user may select which groups from each expression data to compare. The analysis groups together Group1 from both datasets and Group2 from both datasets, and thos become a new 'Group1' and 'Group2' which will be compared in the down stream analysis."),
                                     p(),
                                     p("You may find more information on our ",a("GitHub page.", href="https://github.com/shawlab-moffitt/DRPPM-EASY-LargeProject-Integration"))
                            )
                          )
                          
                        )
                      )
             ),
             
             ####----Data Input Tab----####
             
             tabPanel("Data Input",
                      fluidPage(
                        title = "Data Input",
                        sidebarLayout(
                          sidebarPanel(
                            width = 3,
                            
                            ##--User Data Input--##
                            
                            conditionalPanel(condition = "input.datainput == '1'",
                                             h4("User Matrix Input"),
                                             textInput("MAT1name","Matrix Name:", value = "UserMatrix"),
                                             fileInput("RNAexp_F", "Expression Matrix File Upload", accept = c(".tsv",".txt",".csv")),
                                             checkboxInput("RNAheader","Check if Header in Meta File", value = T),
                                             fileInput("RNAmeta_F", "Meta Matrix File Upload", accept = c(".tsv",".txt",".csv"))
                            ),
                            conditionalPanel(condition = "input.datainput == '2'",
                                             uiOutput("subexprLabel"),
                                             uiOutput("subselectexpr"),
                                             #h4("Subset Expression Data:"),
                                             #selectInput("SelecExpr", "",
                                             #            choices = names(subsetter_list),
                                             #            multiple = F),
                                             uiOutput("subsetselec"),
                                             h4("Condition Selection:"),
                                             selectInput("phenoChoice","",
                                                         choices = db_pheno_choice, selected = "-"),
                                             #h4("Log2 Transform Expression Data:"),
                                             checkboxInput("log2expr","Log2 Transform Expression Data:"),
                                             h4("Label for Project Data:"),
                                             textInput("MAT2name","", value = ProjName)
                            ),
                            conditionalPanel(condition = "input.datainput == '3'",
                                             h4("Comparison Groups From User Data:"),
                                             uiOutput("UserGroupAselec"),
                                             uiOutput("UserGroupBselec"),
                                             h4(paste("Comparison Groups From ",ProjName," Data:", sep = "")),
                                             uiOutput("DBGroupAselec"),
                                             uiOutput("DBGroupBselec"),
                                             h4("Label Group 1 and Group 2:"),
                                             textInput("g1label","Group 1 Label:","Group1"),
                                             textInput("g2label","Group 2 Label:","Group2"),
                                             h4("Download Compiled Meta Table:"),
                                             textInput("CompMetaName", "File Name for Download",
                                                       value = ""),
                                             downloadButton("UserDBmetaDL")
                            )
                          ),
                          mainPanel(
                            tabsetPanel(
                              id = "datainput",
                              tabPanel("Step1: User Data Input",
                                       uiOutput("exprTableLabel"),
                                       div(DT::dataTableOutput("UserExprTable"), style = "font-size:12px; height:450px"),
                                       uiOutput("metaTableLabel"),
                                       div(DT::dataTableOutput("UserMetaTable"), style = "font-size:12px; height:450px"),
                                       value = 1),
                              tabPanel(paste("Step2: ",ProjName," Data Selection", sep = ""),
                                       fluidRow(
                                         column(4,
                                                h4("Meta Data"),
                                                withSpinner(div(DT::dataTableOutput("testtable"), style = "font-size:10px; height:450px"), type = 6),
                                                uiOutput("downloadMETAbutton")),
                                         column(4,
                                                h4("Expression Data"),
                                                withSpinner(div(DT::dataTableOutput("testexprtable"), style = "font-size:10px; height:450px"), type = 6),
                                                uiOutput("downloadEXPRbutton")),
                                         column(4,
                                                uiOutput("NameMapLabel"),
                                                withSpinner(div(DT::dataTableOutput("NameMapTable"), style = "font-size:10px; height:450px"), type = 6),
                                                uiOutput("downloadNameMapbutton"))
                                       ),
                                       value = 2),
                              tabPanel("Step3: Designate Comparison Groups",
                                       div(DT::dataTableOutput("UserDBmeta"), style = "font-size:12px; height:500px"),
                                       value = 3)
                            )
                          )
                        )
                      )
             ),
             
             ####----Expression Scatter Plot----####
             
             tabPanel("Expression Scatter Plot",
                      fluidPage(
                        title = "Expression Scatter Plot",
                        sidebarPanel(
                          width = 3,
                          h4("Group Selection:"),
                          uiOutput("groupAselec"),
                          uiOutput("groupBselec"),
                          h4("Gene Selection:"),
                          uiOutput("GeneSelec"),
                          textInput("gsSelection2", "Text Input of Genes (space or tab delimited):", value = ""),
                          uiOutput("hover_info")
                        ),
                        mainPanel(
                          withSpinner(jqui_resizable(plotOutput("RNAvPROTscatter", height = "600px",
                                                                hover = hoverOpts("plot_hover", delay = 10, delayType = "debounce"))), type = 6),
                          DT::dataTableOutput("RNAvPROTtable"),
                          downloadButton("logFCtableDownload", "Download .tsv")
                        )
                      )
             ),
             
             ####----GSEA----####
             
             tabPanel("Reciprocal GSEA",
                      fluidPage(
                        title = "GSEA",
                        sidebarPanel(
                          width = 3,
                          h4("Group Selection:"),
                          uiOutput("groupAselec2"),
                          uiOutput("groupBselec2"),
                          h4("Gene Set Analysis Parameters:"),
                          numericInput("GSnumber","Number of Top Genes for Gene Set Creation",
                                       min = 5, value = 100),
                          numericInput("gseaPCutt", "P-Value Cutoff:", value = 1.0)
                        ),
                        mainPanel(
                          h3("GSEA Enrichment Plots"),
                          fluidRow(
                            column(6,
                                   #verbatimTextOutput("NESandPval1"),
                                   htmlOutput("NESandPval1", style = "font-size:12px;"),
                                   withSpinner(plotOutput("enrichplot1", width = "500px", height = "450px"), type = 6),
                                   h3("Leading Edge Genes (~Signal2Noise Ranking)"),
                                   div(DT::dataTableOutput("LeadingEdgeGenes1"), style = "font-size:12px; height:500px; overflow-y: scroll")
                            ),
                            column(6,
                                   #verbatimTextOutput("NESandPval2"),
                                   htmlOutput("NESandPval2", style = "font-size:12px;"),
                                   withSpinner(plotOutput("enrichplot2", width = "500px", height = "450px"), type = 6),
                                   h3("Leading Edge Genes (~Signal2Noise Ranking)"),
                                   div(DT::dataTableOutput("LeadingEdgeGenes2"), style = "font-size:12px; height:500px; overflow-y: scroll")
                            )
                          )
                        )
                      )
             ),
             
             ####----ssGSEA----####
             
             tabPanel("Reciprocal ssGSEA",
                      fluidPage(
                        title = "ssGSEA",
                        sidebarPanel(
                          width = 3,
                          h4("Group Selection:"),
                          uiOutput("groupAselec4"),
                          uiOutput("groupBselec4"),
                          h4("Gene Set Analysis Parameters:"),
                          numericInput("GSnumber2","Number of Top Genes for Gene Set Creation",
                                       min = 5, value = 100),
                          numericInput("gseaPCutt2", "P-Value Cutoff:", value = 1.0),
                          selectInput("ssGSEAtype","Choose Scoring Method",
                                      choices = c("ssgsea","gsva","zscore","plage"))
                        ),
                        mainPanel(
                          p(),
                          withSpinner(jqui_resizable(plotOutput('ssGSEAboxplot1', width = "100%", height = "400px")), type = 6),
                          div(DT::dataTableOutput("ssGSEAtable1"), style = "font-size:10px"),
                          downloadButton("ssGSEADown1","Download"),
                          withSpinner(jqui_resizable(plotOutput('ssGSEAboxplot2', width = "100%", height = "400px")), type = 6),
                          div(DT::dataTableOutput("ssGSEAtable2"), style = "font-size:10px"),
                          downloadButton("ssGSEADown2","Download"),
                          withSpinner(jqui_resizable(plotOutput('ssGSEAboxplot3', width = "100%", height = "400px")), type = 6),
                          div(DT::dataTableOutput("ssGSEAtable3"), style = "font-size:10px"),
                          downloadButton("ssGSEADown3","Download"),
                          withSpinner(jqui_resizable(plotOutput('ssGSEAboxplot4', width = "100%", height = "400px")), type = 6),
                          div(DT::dataTableOutput("ssGSEAtable4"), style = "font-size:10px"),
                          downloadButton("ssGSEADown4","Download")
                        )
                      )
             ),
             
             ####----VennDiagram and Stats----####
             
             tabPanel("Statistical Analysis",
                      fluidPage(
                        title = "Statistical Analysis",
                        sidebarPanel(
                          width = 3,
                          h4("Group Selection:"),
                          uiOutput("groupAselec3"),
                          uiOutput("groupBselec3"),
                          h4("Significance Thresholds:"),
                          numericInput("vennLog", "log2FC Threshold:",
                                       min = 0, value = 0, step = 0.1),
                          numericInput("adjpvalVenn", "Adj.P.Val Threshold:",
                                       min = 0, value = 0.05, step = 0.01),
                          tabsetPanel(
                            id = "genesets",
                            tabPanel("MSigDB Gene Sets",
                                     div(DT::dataTableOutput("msigdbtable"), style = "font-size:10px; height:500px; overflow-X: scroll"),
                                     h4("Download Statistics Table:"),
                                     textInput("msigdnldname","File Name for Download"),
                                     downloadButton("msigstatdnld"),
                                     value = 1),
                            tabPanel("User Provided GMT",
                                     fileInput("GMTcomp", "GMT File Upload", accept = ".gmt"),
                                     uiOutput("downloadheader"),
                                     uiOutput("StatNameTextIn"),
                                     uiOutput("StatButton"),
                                     value = 2)
                          )
                        ),
                        mainPanel(
                          fluidRow(
                            column(12,
                                   htmlOutput("VennExplan", style = "font-size:14px;"))
                          ),
                          fluidRow(
                            column(6,
                                   jqui_resizable(plotOutput("vennUP"))
                            ),
                            column(6,
                                   jqui_resizable(plotOutput("vennDN"))
                            )
                          ),
                          fluidRow(
                            column(12,
                                   p(),
                                   htmlOutput("StatExplan", style = "font-size:14px;"),
                                   p(),
                                   withSpinner(DT::dataTableOutput("VennStatTable"), type = 6),
                                   uiOutput("test")
                                   #plotOutput("statBPup"),
                                   #plotOutput("statBPdn")
                            )
                          )
                        )
                      )
             )
  )



####----Server----####


server <- function(input, output, session) {
  
  
  ####----Sample Selection and Database Data Initialization----####
  
  #observe user selection to choose based on lineage or disease, update UI based on selection
  observe({
    output$subsetselec <- renderUI({
      
      if (length(subsetter_list) > 1){
        
        firstChoice_options <- input$SelecExpr
        
        selectInput("firstChoice","Select Sample Subset:",choices = subsetter_list[[firstChoice_options]],
                    multiple = F)
        
      }
    })
  })
  
  #Reactive for db Meta table - TEMP
  meta_c_pre <- reactive({
    
    #assign user input values
    firstChoice <- input$firstChoice
    phenoChoice <- input$phenoChoice
    exprChoice <- gsub("_Choices","",input$SelecExpr)
    
    if (is.null(db_meta_selec) == F) {
      
      #get samples selected based on lineage choice
      samp_selec <- meta_subsetters[which(meta_subsetters[,2] == exprChoice & meta_subsetters[,3] == firstChoice),1]
      #subset phenotype data based on selected phenotype
      db_pheno.u <- meta_phenos[which(meta_phenos[,2] == phenoChoice),]
      #use samples selected to subset again phenotype meta
      db_meta.u <- db_pheno.u[which(db_pheno.u[,1] %in% samp_selec),]
      #reformat to proper meta table
      db_meta.u2 <- db_meta.u[,c(1,3)]
      colnames(db_meta.u2)[2] <- "Type"
      db_meta_sub <- db_meta.u2
      
    }
    if (is.null(db_meta_selec) == T) {
      
      #subset phenotype data based on selected phenotype
      db_pheno.u <- meta_phenos[which(meta_phenos[,2] == phenoChoice),]
      #reformat to proper meta table
      db_meta.u2 <- db_pheno.u[,c(1,3)]
      colnames(db_meta.u2)[2] <- "Type"
      db_meta_sub <- db_meta.u2
      
    }
    
    db_meta_sub
    
  })
  
  #Reactive for db expression data - TEMP
  expr_c_pre <- reactive({
    
    #assign user choices
    firstChoice <- input$firstChoice
    phenoChoice <- input$phenoChoice
    exprChoice <- gsub("_Choices","",input$SelecExpr)
    
    samp_selec <- unlist(meta_c_pre()[,1])
    db_expr_sub <- db_expr[,which(colnames(db_expr) %in% samp_selec), drop = F]
    db_expr_sub$gene <- rownames(db_expr_sub)
    db_expr_sub <- db_expr_sub %>%
      relocate(gene)
    
    #drop any rows with NA from the new subset data
    expr <- db_expr_sub %>%
      drop_na()
    #Assign the gene names to be row names, make unique if duplicates
    row.names(expr) <- make.names(expr[,1], unique = T)
    #remove first column of gene names
    expr <- expr[,-1]
    #make special characters uniform to .
    colnames(expr) <- gsub("[_.-]", ".", colnames(expr))
    #log transform if checked
    if (input$log2expr == TRUE) {
      
      expr_mat <- as.matrix(expr)
      expr <- as.data.frame(log2(expr + 1))
      
    }
    
    expr
    
  })
  
  #reactive for db expression and meta ensuring similar sample names
  expr_c <- reactive({
    
    #Non intersected db data
    expr <- expr_c_pre()
    meta <- meta_c_pre()
    
    sampsames <- intersect(colnames(expr),meta[,1])
    expr <- expr[,sampsames]
    
    expr
    
  })
  meta_c <- reactive({
    
    #Non intersected db data
    expr <- expr_c_pre()
    meta <- meta_c_pre()
    
    sampsames <- intersect(colnames(expr),meta[,1])
    meta <- meta[which(meta[,1] %in% sampsames),]
    
    meta
    
  })
  
  #db expression data table
  output$testexprtable <- DT::renderDataTable({
    DT::datatable(expr_c(),
                  options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"), scrollX = T))
    
  })
  
  #db meta data table
  output$testtable <- DT::renderDataTable({
    DT::datatable(meta_c(),
                  rownames = F,
                  options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"), scrollX = T))
    
  })
  
  output$NameMapLabel <- renderUI({
    
    if (is.null(db_namemap) == FALSE) {
      
      h4("Name Map Guide")
      
    }
    
  })
  
  # name map table
  output$NameMapTable <- DT::renderDataTable({
    
    if (is.null(db_namemap) == FALSE) {
      
      namemap <- db_namemap
      meta.new <- meta_c()
      samp_selected <- unlist(meta.new[,1])         
      namemap_selcted <- namemap[which(namemap[,1] %in% samp_selected),]
      DT::datatable(namemap_selcted,
                    rownames = F,
                    options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"), scrollX = T))
      
    }
    
  })
  
  output$downloadNameMapbutton <- renderUI({
    
    if (is.null(db_namemap) == FALSE) {
      
      downloadButton("downloadNameMap","Download Name Map")
      
    }
    
  })
  
  output$downloadNameMap <- downloadHandler(
    filename = function() {
      #meta.cond <- meta()
      firstChoice <- input$firstChoice
      phenoChoice <- input$phenoChoice
      exprChoice <- gsub("_Choices","",input$SelecExpr)
      paste("NameMap_",ProjName,"_",gsub(" ","",firstChoice),"_",phenoChoice,".tsv", sep = "")
    },
    content = function(file) {
      namemap <- db_namemap
      meta.new <- meta_c()
      samp_selected <- unlist(meta.new[,1])
      namemap_selcted <- namemap[which(namemap[,1] %in% samp_selected),]
      write_tsv(namemap_selcted,file)
      
    }
  )
  
  output$downloadMETAbutton <- renderUI({
    #req(input$updateMETA)
    downloadButton("downloadMETA","Download Meta Data")
  })
  
  output$downloadEXPRbutton <- renderUI({
    #req(input$updateMETA)
    downloadButton("downloadEXPR","Download Expression Data")
  })
  
  
  output$downloadMETA <- downloadHandler(
    filename = function() {
      #meta.cond <- meta()
      firstChoice <- input$firstChoice
      phenoChoice <- input$phenoChoice
      paste("meta_",ProjName,"_",phenoChoice,".tsv",sep = "")
    },
    content = function(file) {
      
      
      firstChoice <- input$firstChoice
      phenoChoice <- input$phenoChoice
      exprChoice <- gsub("_Choices","",input$SelecExpr)
      
      if (is.null(db_meta_selec) == F) {
        
        #get samples selected based on lineage choice
        samp_selec <- meta_subsetters[which(meta_subsetters[,2] == exprChoice & meta_subsetters[,3] == firstChoice),1]
        #subset phenotype data based on selected phenotype
        db_pheno.u <- meta_phenos[which(meta_phenos[,2] == phenoChoice),]
        #use samples selected to subset again phenotype meta
        db_meta.u <- db_pheno.u[which(db_pheno.u[,1] %in% samp_selec),]
        #reformat to proper meta table
        db_meta.u2 <- db_meta.u[,c(1,3)]
        colnames(db_meta.u2)[2] <- "Type"
        db_meta_sub <- db_meta.u2
        
      }
      if (is.null(db_meta_selec) == T) {
        
        #subset phenotype data based on selected phenotype
        db_pheno.u <- meta_phenos[which(meta_phenos[,2] == phenoChoice),]
        #reformat to proper meta table
        db_meta.u2 <- db_pheno.u[,c(1,3)]
        colnames(db_meta.u2)[2] <- "Type"
        db_meta_sub <- db_meta.u2
        
      }
      write_tsv(db_meta_sub,file)
    }
  )
  
  output$downloadEXPR <- downloadHandler(
    filename = function() {
      
      firstChoice <- input$firstChoice
      phenoChoice <- input$phenoChoice
      paste("expr_",ProjName,"_",phenoChoice,".tsv",sep = "")
    },
    content = function(file) {
      
      ##--update EXPR section--##
      #assign user lineage/disease choices
      firstChoice <- input$firstChoice
      phenoChoice <- input$phenoChoice
      exprChoice <- gsub("_Choices","",input$SelecExpr)
      
      if (is.null(db_meta_selec) == F) {
        
        #get samples selected based on lineage choice
        samp_selec <- meta_subsetters[which(meta_subsetters[,2] == exprChoice & meta_subsetters[,3] == firstChoice),1]
        #subset phenotype data based on selected phenotype
        db_pheno.u <- meta_phenos[which(meta_phenos[,2] == phenoChoice),]
        #use samples selected to subset again phenotype meta
        db_meta.u <- db_pheno.u[which(db_pheno.u[,1] %in% samp_selec),]
        #reformat to proper meta table
        db_meta.u2 <- db_meta.u[,c(1,3)]
        colnames(db_meta.u2)[2] <- "Type"
        db_meta_sub <- db_meta.u2
        
      }
      if (is.null(db_meta_selec) == T) {
        
        #subset phenotype data based on selected phenotype
        db_pheno.u <- meta_phenos[which(meta_phenos[,2] == phenoChoice),]
        #reformat to proper meta table
        db_meta.u2 <- db_pheno.u[,c(1,3)]
        colnames(db_meta.u2)[2] <- "Type"
        db_meta_sub <- db_meta.u2
        
      }
      
      samp_selec <- unlist(db_meta_sub[,1])
      db_expr_sub <- db_expr[,which(colnames(db_expr) %in% samp_selec), drop = F]
      db_expr_sub$gene <- rownames(db_expr_sub)
      db_expr_sub <- db_expr_sub %>%
        relocate(gene)
      write_tsv(db_expr_sub,file)
      
    }
  )
  
  
  ####----Render UI----####
  
  output$subexprLabel <- renderUI({
    
    if (length(subsetter_list) > 1) {
      
      h4("Subset Expression Data:")
      
    }
    
  })
  
  output$subselectexpr <- renderUI({
    
    if (length(subsetter_list) > 1) {
      
      selectInput("SelecExpr", "",
                  choices = names(subsetter_list),
                  multiple = F)
      
    }
    
  })
  
  
  output$exprTableLabel <- renderUI({
    
    req(input$RNAexp_F)
    h4("User Expression Data")
    
  })
  
  output$metaTableLabel <- renderUI({
    
    req(input$RNAmeta_F)
    h4("User Meta Data")
    
  })
  
  output$UserGroupAselec <- renderUI({
    
    metaR <- RNAmeta()
    metagroups <- as.vector(levels(factor(metaR[,2])))
    selectInput("UserGroup1", "Comparison Group 1:",
                choices = metagroups, selected = metagroups[1])
    
  })
  
  output$UserGroupBselec <- renderUI({
    
    metaR <- RNAmeta()
    metagroups <- as.vector(levels(factor(metaR[,2])))
    selectInput("UserGroup2", "Comparison Group 2:",
                choices = metagroups, selected = metagroups[2])
    
  })
  
  output$DBGroupAselec <- renderUI({
    
    metaR <- meta_c()
    metagroups <- as.vector(levels(factor(metaR[,2])))
    selectInput("DBGroup1", "Comparison Group 1:",
                choices = metagroups, selected = metagroups[1])
    
  })
  
  output$DBGroupBselec <- renderUI({
    
    metaR <- meta_c()
    metagroups <- as.vector(levels(factor(metaR[,2])))
    selectInput("DBGroup2", "Comparison Group 2:",
                choices = metagroups, selected = metagroups[2])
    
  })
  
  
  #user selection for group A based on expression matrices
  output$groupAselec <- renderUI({
    
    #metaR <- RNAmeta()
    #metaP <- meta_c()
    metaR <- meta_un()
    metaP <- meta_cn()
    #metasame <- merge(metaR,metaP)
    metasame <- intersect(metaR[,2],metaP[,2])
    #metasame <- union(metaR[,2],metaP[,2])
    metagroups <- as.vector(levels(factor(metasame)))
    selectInput("groupAcomp", "Comparison Group A:",
                choices = metagroups, selected = metagroups[1])
    
  })
  
  #user selection for group B based on expression matrices
  output$groupBselec <- renderUI({
    
    #metaR <- RNAmeta()
    #metaP <- meta_c()
    metaR <- meta_un()
    metaP <- meta_cn()
    #metasame <- merge(metaR,metaP)
    metasame <- intersect(metaR[,2],metaP[,2])
    #metasame <- union(metaR[,2],metaP[,2])
    metagroups <- as.vector(levels(factor(metasame)))
    selectInput("groupBcomp", "Comparison Group B:",
                choices = metagroups, selected = metagroups[2])
    
  })
  
  #user selection for group A based on expression matrices - GSEA
  output$groupAselec2 <- renderUI({
    
    #metaR <- RNAmeta()
    #metaP <- meta_c()
    metaR <- meta_un()
    metaP <- meta_cn()
    #metasame <- merge(metaR,metaP)
    metasame <- intersect(metaR[,2],metaP[,2])
    #metasame <- union(metaR[,2],metaP[,2])
    metagroups <- as.vector(levels(factor(metasame)))
    selectInput("groupAcomp2", "Comparison Group A:",
                choices = metagroups, selected = metagroups[1])
    
  })
  
  #user selection for group B based on expression matrices - GSEA
  output$groupBselec2 <- renderUI({
    
    #metaR <- RNAmeta()
    #metaP <- meta_c()
    metaR <- meta_un()
    metaP <- meta_cn()
    #metasame <- merge(metaR,metaP)
    metasame <- intersect(metaR[,2],metaP[,2])
    #metasame <- union(metaR[,2],metaP[,2])
    metagroups <- as.vector(levels(factor(metasame)))
    selectInput("groupBcomp2", "Comparison Group B:",
                choices = metagroups, selected = metagroups[2])
    
  })
  
  #user selection for group A based on expression matrices - GSEA
  output$groupAselec3 <- renderUI({
    
    #metaR <- RNAmeta()
    #metaP <- meta_c()
    metaR <- meta_un()
    metaP <- meta_cn()
    #metasame <- merge(metaR,metaP)
    metasame <- intersect(metaR[,2],metaP[,2])
    #metasame <- union(metaR[,2],metaP[,2])
    metagroups <- as.vector(levels(factor(metasame)))
    selectInput("groupAcomp3", "Comparison Group A:",
                choices = metagroups, selected = metagroups[1])
    
  })
  
  #user selection for group B based on expression matrices - GSEA
  output$groupBselec3 <- renderUI({
    
    #metaR <- RNAmeta()
    #metaP <- meta_c()
    metaR <- meta_un()
    metaP <- meta_cn()
    #metasame <- merge(metaR,metaP)
    metasame <- intersect(metaR[,2],metaP[,2])
    #metasame <- union(metaR[,2],metaP[,2])
    metagroups <- as.vector(levels(factor(metasame)))
    selectInput("groupBcomp3", "Comparison Group B:",
                choices = metagroups, selected = metagroups[2])
    
  })
  
  #user selection for group A based on expression matrices - GSEA
  output$groupAselec4 <- renderUI({
    
    #metaR <- RNAmeta()
    #metaP <- meta_c()
    metaR <- meta_un()
    metaP <- meta_cn()
    #metasame <- merge(metaR,metaP)
    metasame <- intersect(metaR[,2],metaP[,2])
    #metasame <- union(metaR[,2],metaP[,2])
    metagroups <- as.vector(levels(factor(metasame)))
    selectInput("groupAcomp4", "Comparison Group A:",
                choices = metagroups, selected = metagroups[1])
    
  })
  
  #user selection for group B based on expression matrices - GSEA
  output$groupBselec4 <- renderUI({
    
    #metaR <- RNAmeta()
    #metaP <- meta_c()
    metaR <- meta_un()
    metaP <- meta_cn()
    #metasame <- merge(metaR,metaP)
    metasame <- intersect(metaR[,2],metaP[,2])
    #metasame <- union(metaR[,2],metaP[,2])
    metagroups <- as.vector(levels(factor(metasame)))
    selectInput("groupBcomp4", "Comparison Group B:",
                choices = metagroups[2])
    
  })
  
  #render gene selection list for RNAvProt scatter
  output$GeneSelec <- renderUI({
    
    RNA <- RNAmat()
    prot <- expr_c()
    comgenes <- intersect(rownames(RNA),rownames(prot))
    selectInput("scatterGeneSelec", "Select Genes to Identify in Plot",
                choices = comgenes,
                multiple = T,
                selected = "-")
    
  })
  
  #render box plot stat compare method
  output$boxplotstat <- renderUI({
    
    metaR <- RNAmeta()
    metaP <- meta_c()
    #metasame <- merge(metaR,metaP)
    metasame <- intersect(metaR[,2],metaP[,2])
    #metasame <- union(metaR[,2],metaP[,2])
    metagroups <- as.vector(levels(factor(metasame)))
    #boxplot choices based on meta groups
    if (length(metagroups) == 2) {
      boxopt <- c("wilcox.test", "t.test", "none")
    }
    if (length(metagroups) >= 3) {
      boxopt <- c("kruskal.test", "anova", "none")
    }
    selectInput("boxplotcompare", "BoxPlot Stat Compare Method:",
                choices = boxopt)
    
  })
  
  #render download button when gmt file input
  output$StatButton <- renderUI({
    
    req(input$GMTcomp)
    downloadButton("StatTabDown","Download Table")
    
  })
  
  #render text input for download name if GMT file input
  output$StatNameTextIn <- renderUI({
    
    req(input$GMTcomp)
    textInput("StatTabName","File Name for Download:")
    
  })
  
  output$downloadheader <- renderUI({
    
    req(input$GMTcomp)
    h4("Download Statistics Table:")
    
  })
  
  output$hover_info <- renderUI({
    FCtable <- RNAvProtFC()
    M1 <- input$MAT1name
    M2 <- input$MAT2name
    #get delta
    FCtable$Delta_M1M2 <- FCtable[,2] - FCtable[,3]
    #edit column names for hover text to find
    colnames(FCtable) <- c("gene","log2FC_M1","log2FC_M2","DeltaLog2FC")
    hover <- input$plot_hover
    point <- nearPoints(FCtable, hover, threshold = 10, maxpoints = 1, addDist = FALSE)
    if (nrow(point) == 0) return(NULL)
    wellPanel(
      p(HTML(paste0("<b> Gene: </b>", point[1], "<br/>",
                    "<b> ",M1," Log2FC: </b>", round(point[2], digits = 4), "<br/>",
                    "<b> ",M2," Log2FC: </b>", round(point[3], digits = 4), "<br/>",
                    "<b> Delta Log2FC Between ",M1," and ",M2,": </b>", round(point[4], digits = 4), "<br/>",
                    NULL
      ))))
  })
  
  
  ####----Reactives----####
  
  #render RNA seq expression matrix from file inputs
  RNAmat <- reactive({
    
    RNA.u <- input$RNAexp_F
    ext <- tools::file_ext(RNA.u$datapath)
    req(RNA.u)
    validate(need(ext == c("csv","tsv","txt"), "Please upload tab delimited .tsv or .txt file"))
    if (ext == "tsv") {
      exprR <- as.data.frame(read_delim(RNA.u$datapath, delim = '\t'))
    }
    else if (ext == "txt") {
      exprR <- as.data.frame(read_delim(RNA.u$datapath, delim = '\t'))
    }
    else if (ext == "csv") {
      exprR <- as.data.frame(read_delim(RNA.u$datapath, delim = ','))
    }
    colnames(exprR)[1] <- "Gene"
    exprR <- exprR %>%
      drop_na()
    row.names(exprR) <- exprR[,1]
    exprR <- exprR[,-1]
    colnames(exprR) <- gsub("[_.-]", ".", colnames(exprR))
    exprR
    
  })
  
  #render RNA seq Meta table
  RNAmeta <- reactive({
    
    RNA.u.m <- input$RNAmeta_F
    ext <- tools::file_ext(RNA.u.m$datapath)
    req(RNA.u.m)
    validate(need(ext == c("csv","tsv","txt"), "Please upload tab delimited .tsv or .txt file"))
    if (ext == "tsv") {
      metaR <- read.delim(RNA.u.m$datapath, sep = '\t', header = input$RNAheader)
    }
    else if (ext == "txt") {
      metaR <- read.delim(RNA.u.m$datapath, sep = '\t', header = input$RNAheader)
    }
    else if (ext == "csv") {
      metaR <- read.delim(RNA.u.m$datapath, sep = ',', header = input$RNAheader)
    }
    metaR[,1] <- gsub("[_.-]", ".", metaR[,1])
    colnames(metaR) <- c("SampleName","Type")
    metaR
    
  })
  
  #render RNAseq and Protein FC table
  RNAvProtFC <- reactive({
    
    #Variables
    RNA <- RNAmat()
    prot <- expr_c()
    metaR <- meta_un()
    metaP <- meta_cn()
    #metasame <- merge(metaR,metaP)
    A1 <- metaR[which(metaR[,2] == input$groupAcomp),1]
    B1 <- metaR[which(metaR[,2] == input$groupBcomp),1]
    
    A2 <- metaP[which(metaP[,2] == input$groupAcomp),1]
    B2 <- metaP[which(metaP[,2] == input$groupBcomp),1]
    
    #RNAseq top table generation
    mat <- RNA[,c(A1,B1)]
    mat <- log2(mat + 1.0)
    groupAOther <- factor(c(rep("A", length(A1)), rep("B", length(B1))))
    designA <- model.matrix(~0 + groupAOther)
    fit <- lmFit(mat, design = designA)
    contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    options(digits = 4)
    top1_R <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
    
    #Protein top table generation
    mat <- prot[,c(A2,B2)]
    mat <- log2(mat + 1.0)
    groupAOther <- factor(c(rep("A", length(A2)), rep("B", length(B2))))
    designA <- model.matrix(~0 + groupAOther)
    fit <- lmFit(mat, design = designA)
    contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    options(digits = 4)
    top1_P <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
    
    #merge tables
    top1_R2 <- top1_R
    top1_R2$gene <-rownames(top1_R2)
    top1_R3 <- top1_R2 %>%
      select(gene,logFC)
    top1_P2 <- top1_P
    top1_P2$gene <-rownames(top1_P2)
    top1_P3 <- top1_P2 %>%
      select(gene,logFC)
    logFC_RP <- merge(top1_R3,top1_P3, by = "gene")
    colnames(logFC_RP) <- c("gene", input$MAT1name, input$MAT2name)
    logFC_RP
    
  })
  
  #generate GMT
  recGMT <- reactive({
    
    #Variables
    RNA <- RNAmat()
    prot <- expr_c()
    metaR <- meta_un()
    metaP <- meta_cn()
    #metasame <- merge(metaR,metaP)
    #metagroupsR <- as.vector(levels(factor(metaR[,2])))
    #metagroupsP <- as.vector(levels(factor(metaP[,2])))
    A1 <- metaR[which(metaR[,2] == input$groupAcomp2),1]
    B1 <- metaR[which(metaR[,2] == input$groupBcomp2),1]
    
    A2 <- metaP[which(metaP[,2] == input$groupAcomp2),1]
    B2 <- metaP[which(metaP[,2] == input$groupBcomp2),1]
    
    #RNAseq top table generation
    mat <- RNA[,c(A1,B1)]
    mat <- log2(mat + 1.0)
    groupAOther <- factor(c(rep("A", length(A1)), rep("B", length(B1))))
    designA <- model.matrix(~0 + groupAOther)
    fit <- lmFit(mat, design = designA)
    contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    options(digits = 4)
    top1_R <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
    
    #Protein top table generation
    mat <- prot[,c(A2,B2)]
    mat <- log2(mat + 1.0)
    groupAOther <- factor(c(rep("A", length(A2)), rep("B", length(B2))))
    designA <- model.matrix(~0 + groupAOther)
    fit <- lmFit(mat, design = designA)
    contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    options(digits = 4)
    top1_P <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
    
    #subset top 100 up and down regulated genes from RNA and Protein expression
    top1_RU <- top1_R[which(top1_R$logFC > 0),]
    topGenes_RU <- top1_RU[c(1:input$GSnumber),]
    top1_RD <- top1_R[which(top1_R$logFC < 0),]
    topGenes_RD <- top1_RD[c(1:input$GSnumber),]
    top1_PU <- top1_P[which(top1_P$logFC > 0),]
    topGenes_PU <- top1_PU[c(1:input$GSnumber),]
    top1_PD <- top1_P[which(top1_P$logFC < 0),]
    topGenes_PD <- top1_PD[c(1:input$GSnumber),]
    
    #convert genes into TSV gmt format for GSEA
    M1 <- input$MAT1name
    M2 <- input$MAT2name
    topGenes_RU_gmt <- data.frame(paste("UPreg_in_",gsub(" ","",input$groupAcomp2),"_",M1,sep = ""),rownames(topGenes_RU))
    colnames(topGenes_RU_gmt) <- c("term","gene")
    topGenes_RD_gmt <- data.frame(paste("DNreg_in_",gsub(" ","",input$groupAcomp2),"_",M1,sep = ""),rownames(topGenes_RD))
    colnames(topGenes_RD_gmt) <- c("term","gene")
    topGenes_PU_gmt <- data.frame(paste("UPreg_in_",gsub(" ","",input$groupAcomp2),"_",M2,sep = ""),rownames(topGenes_PU))
    colnames(topGenes_PU_gmt) <- c("term","gene")
    topGenes_PD_gmt <- data.frame(paste("DNreg_in_",gsub(" ","",input$groupAcomp2),"_",M2,sep = ""),rownames(topGenes_PD))
    colnames(topGenes_PD_gmt) <- c("term","gene")
    
    #combine to 4 gene set list
    RNAvProt_GS <- rbind(topGenes_RU_gmt,topGenes_RD_gmt,topGenes_PU_gmt,topGenes_PD_gmt)
    RNAvProt_GS
    
  })
  
  #perform GSEA function with made GMT - GROUP 1
  recGSEAg1 <- reactive({
    
    #Variables
    A <- as.matrix(RNAmat())
    #prot <- as.matrix(PROTmat())
    metaR <- meta_un()
    metaP <- meta_cn()
    #metasame <- merge(metaR,metaP)
    #metasame <- intersect(metaR[,2],metaP[,2])
    groupA <- metaR[which(metaR[,2] == input$groupAcomp2),1]
    groupB <- metaR[which(metaR[,2] == input$groupBcomp2),1]
    gmt <- recGMT()
    ##----Signal-to-Noise Calculation----##
    A <- A + 0.00000001
    P = as.matrix(as.numeric(colnames(A) %in% groupA))
    n1 <- sum(P[,1])
    M1 <- A %*% P
    M1 <- M1/n1
    A2 <- A*A
    S1 <- A2 %*% P
    S1 <- S1/n1 - M1*M1 
    S1 <- sqrt(abs((n1/(n1-1)) * S1))
    P = as.matrix(as.numeric(colnames(A) %in% groupB))
    n2 <- sum(P[,1])
    M2 <- A %*% P
    M2 <- M2/n2
    A2 <- A*A
    S2 <- A2 %*% P
    S2 <- S2/n2 - M2*M2
    S2 <- sqrt(abs((n2/(n2-1)) * S2))
    rm(A2)
    # small sigma "fix" as used in GeneCluster
    S2 <- ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
    S2 <- ifelse(S2 == 0, 0.2, S2)
    S1 <- ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
    S1 <- ifelse(S1 == 0, 0.2, S1)
    M1 <- M1 - M2
    rm(M2)
    S1 <- S1 + S2
    rm(S2)
    s2n.matrix <- M1/S1
    ##----Reformatting----##
    s2n.df <- as.data.frame(s2n.matrix)
    s2n.df$GeneID <- rownames(s2n.df)
    rownames(s2n.df) <- NULL
    data <- dplyr::select(s2n.df, GeneID, V1)
    data.gsea <- data$V1
    names(data.gsea) <- as.character(data$GeneID)
    s2n.matrix.s <- sort(data.gsea, decreasing = T)
    ##----GSEA Function----##
    M2 <- input$MAT2name
    sets <- unique(gmt[,1])
    GS <- sets[grep(paste('_',M2,sep = ''),unique(gmt[,1]))]
    GSEA(s2n.matrix.s, TERM2GENE = gmt[which(gmt[,1] == GS),],
         verbose = F, pvalueCutoff = input$gseaPCutt)
    
  })
  
  #perform GSEA function with made GMT - GROUP 2
  recGSEAg2 <- reactive({
    
    #Variables
    #RNA <- as.matrix(RNAmat())
    A <- as.matrix(expr_c())
    metaR <- meta_un()
    metaP <- meta_cn()
    #metasame <- merge(metaR,metaP)
    #metasame <- intersect(metaR[,2],metaP[,2])
    groupA <- metaP[which(metaP[,2] == input$groupAcomp2),1]
    groupB <- metaP[which(metaP[,2] == input$groupBcomp2),1]
    gmt <- recGMT()
    ##----Signal-to-Noise Calculation----##
    A <- A + 0.00000001
    P = as.matrix(as.numeric(colnames(A) %in% groupA))
    n1 <- sum(P[,1])
    M1 <- A %*% P
    M1 <- M1/n1
    A2 <- A*A
    S1 <- A2 %*% P
    S1 <- S1/n1 - M1*M1 
    S1 <- sqrt(abs((n1/(n1-1)) * S1))
    P = as.matrix(as.numeric(colnames(A) %in% groupB))
    n2 <- sum(P[,1])
    M2 <- A %*% P
    M2 <- M2/n2
    A2 <- A*A
    S2 <- A2 %*% P
    S2 <- S2/n2 - M2*M2
    S2 <- sqrt(abs((n2/(n2-1)) * S2))
    rm(A2)
    # small sigma "fix" as used in GeneCluster
    S2 <- ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
    S2 <- ifelse(S2 == 0, 0.2, S2)
    S1 <- ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
    S1 <- ifelse(S1 == 0, 0.2, S1)
    M1 <- M1 - M2
    rm(M2)
    S1 <- S1 + S2
    rm(S2)
    s2n.matrix <- M1/S1
    ##----Reformatting----##
    s2n.df <- as.data.frame(s2n.matrix)
    s2n.df$GeneID <- rownames(s2n.df)
    rownames(s2n.df) <- NULL
    data <- dplyr::select(s2n.df, GeneID, V1)
    data.gsea <- data$V1
    names(data.gsea) <- as.character(data$GeneID)
    s2n.matrix.s <- sort(data.gsea, decreasing = T)
    ##----GSEA Function----##
    M1 <- input$MAT1name
    sets <- unique(gmt[,1])
    GS <- sets[grep(paste('_',M1,sep = ''),unique(gmt[,1]))]
    GSEA(s2n.matrix.s, TERM2GENE = gmt[which(gmt[,1] == GS),],
         verbose = F, pvalueCutoff = input$gseaPCutt)
    
  })
  
  #generate GMT
  
  recGMTssGSEA <- reactive({
    
    #Variables
    RNA <- RNAmat()
    prot <- expr_c()
    metaR <- meta_un()
    metaP <- meta_cn()
    #metasame <- merge(metaR,metaP)
    #metagroupsR <- as.vector(levels(factor(metaR[,2])))
    #metagroupsP <- as.vector(levels(factor(metaP[,2])))
    A1 <- metaR[which(metaR[,2] == input$groupAcomp4),1]
    B1 <- metaR[which(metaR[,2] == input$groupBcomp4),1]
    
    A2 <- metaP[which(metaP[,2] == input$groupAcomp4),1]
    B2 <- metaP[which(metaP[,2] == input$groupBcomp4),1]
    
    #RNAseq top table generation
    mat <- RNA[,c(A1,B1)]
    mat <- log2(mat + 1.0)
    groupAOther <- factor(c(rep("A", length(A1)), rep("B", length(B1))))
    designA <- model.matrix(~0 + groupAOther)
    fit <- lmFit(mat, design = designA)
    contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    options(digits = 4)
    top1_R <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
    
    #Protein top table generation
    mat <- prot[,c(A2,B2)]
    mat <- log2(mat + 1.0)
    groupAOther <- factor(c(rep("A", length(A2)), rep("B", length(B2))))
    designA <- model.matrix(~0 + groupAOther)
    fit <- lmFit(mat, design = designA)
    contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    options(digits = 4)
    top1_P <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
    
    #subset top 100 up and down regulated genes from RNA and Protein expression
    top1_RU <- top1_R[which(top1_R$logFC > 0),]
    topGenes_RU <- top1_RU[c(1:input$GSnumber2),]
    top1_RD <- top1_R[which(top1_R$logFC < 0),]
    topGenes_RD <- top1_RD[c(1:input$GSnumber2),]
    top1_PU <- top1_P[which(top1_P$logFC > 0),]
    topGenes_PU <- top1_PU[c(1:input$GSnumber2),]
    top1_PD <- top1_P[which(top1_P$logFC < 0),]
    topGenes_PD <- top1_PD[c(1:input$GSnumber2),]
    
    #convert genes into TSV gmt format for GSEA
    M1 <- input$MAT1name
    M2 <- input$MAT2name
    topGenes_RU_gmt <- data.frame(paste("UPreg_in_",gsub(" ","",input$groupAcomp4),"_",M1,sep = ""),rownames(topGenes_RU))
    colnames(topGenes_RU_gmt) <- c("term","gene")
    topGenes_RD_gmt <- data.frame(paste("DNreg_in_",gsub(" ","",input$groupAcomp4),"_",M1,sep = ""),rownames(topGenes_RD))
    colnames(topGenes_RD_gmt) <- c("term","gene")
    topGenes_PU_gmt <- data.frame(paste("UPreg_in_",gsub(" ","",input$groupAcomp4),"_",M2,sep = ""),rownames(topGenes_PU))
    colnames(topGenes_PU_gmt) <- c("term","gene")
    topGenes_PD_gmt <- data.frame(paste("DNreg_in_",gsub(" ","",input$groupAcomp4),"_",M2,sep = ""),rownames(topGenes_PD))
    colnames(topGenes_PD_gmt) <- c("term","gene")
    
    #combine to 4 gene set list
    RNAvProt_GS <- rbind(topGenes_RU_gmt,topGenes_RD_gmt,topGenes_PU_gmt,topGenes_PD_gmt)
    RNAvProt_GS
    
  })
  
  #ssGSEA reactive - group1
  ssgseaG1 <- reactive({
    
    #make RData list with generated GMT - M1
    M2 <- input$MAT2name
    gmt <- recGMTssGSEA()
    A <- as.matrix(RNAmat())
    gsDataList <- list()
    for (i in unique(gmt[,1])){
      gsDataList[[i]] <- gmt[gmt[,1] == i,]$gene
    }
    gsva(A, gsDataList, method = input$ssGSEAtype, verbose = F)
    
  })
  
  #ssGSEA reactive - group2
  ssgseaG2 <- reactive({
    
    #make RData list with generated GMT - M1
    M1 <- input$MAT1name
    gmt <- recGMTssGSEA()
    A <- as.matrix(expr_c())
    gsDataList <- list()
    for (i in unique(gmt[,1])){
      gsDataList[[i]] <- gmt[gmt[,1] == i,]$gene
    }
    gsva(A, gsDataList, method = input$ssGSEAtype, verbose = F)
    
  })
  
  #ssGSEA master table
  ssGSEAmaster <- reactive({
    
    #meta info
    metaR <- meta_un()
    metaP <- meta_cn()
    #metasame <- merge(metaR,metaP)
    #metasame <- intersect(metaR[,2],metaP[,2])
    #meta types - no white space
    T1 <- gsub(" ","",input$groupAcomp4)
    T2 <- gsub(" ","",input$groupBcomp4)
    samporder1 <- metaR[,1]
    samporder2 <- metaP[,1]
    M1 <- input$MAT1name
    M2 <- input$MAT2name
    #ssgsea results
    ssgseaM1 <- ssgseaG1()
    ssgseaM2 <- ssgseaG2()
    
    #table from M1
    ssgseaM1_2 <- as.data.frame(t(ssgseaM1))
    ssgseaM1_3 <- as.data.frame(ssgseaM1_2[samporder1,])
    ssgseaM1_4 <- ssgseaM1_3 %>% 
      mutate(type = case_when(
        rownames(ssgseaM1_3) == metaR[,1] ~ metaR[,2],
      ))
    ssgseaM1_4 <- ssgseaM1_4 %>%
      relocate(type)
    ssgseaM1_4$Sample <- rownames(ssgseaM1_4)
    ssgseaM1_5 <- reshape2::melt(ssgseaM1_4)
    ssgseaM1_5$Matrix <- paste(M1,"_Matrix",sep = "")
    
    #table from M2
    ssgseaM2_2 <- as.data.frame(t(ssgseaM2))
    ssgseaM2_3 <- as.data.frame(ssgseaM2_2[samporder2,])
    ssgseaM2_4 <- ssgseaM2_3 %>% 
      mutate(type = case_when(
        rownames(ssgseaM2_3) == metaP[,1] ~ metaP[,2],
      ))
    ssgseaM2_4 <- ssgseaM2_4 %>%
      relocate(type)
    ssgseaM2_4$Sample <- rownames(ssgseaM2_4)
    ssgseaM2_5 <- reshape2::melt(ssgseaM2_4)
    ssgseaM2_5$Matrix <- paste(M2,"_Matrix",sep = "")
    #merge tables
    ssgseaM <- rbind(ssgseaM1_5,ssgseaM2_5)
    #change down regulated in type 1 to up regulated in type 2
    ssgseaM$variable <- gsub(paste("DNreg_in_",T1,sep = ""),paste("UPreg_in_",T2,sep = ""),ssgseaM$variable)
    ssgseaM
    
  })
  
  meta_cn <- reactive({
    
    #load in files
    meta.c <- meta_c()
    #Subset sample names for group 1 and 2 in each meta
    g1.c <- meta.c[which(meta.c[,2] == input$DBGroup1),1]
    g2.c <- meta.c[which(meta.c[,2] == input$DBGroup2),1]
    #group subsets between user and DB data
    g1 <- g1.c
    g2 <- g2.c
    #make new meta data frame
    meta.n <- data.frame(c(g1,g2))
    colnames(meta.n)[1] <- "SampleName"
    meta.n <- meta.n %>%
      mutate(Type = case_when(
        SampleName %in% g1 ~ input$g1label,
        SampleName %in% g2 ~ input$g2label
      )
      )
    meta.n
    
  })
  
  meta_un <- reactive({
    
    #load in files
    meta.u <- RNAmeta()
    #Subset sample names for group 1 and 2 in each meta
    g1.u <- meta.u[which(meta.u[,2] == input$UserGroup1),1]
    g2.u <- meta.u[which(meta.u[,2] == input$UserGroup2),1]
    #group subsets between user and DB data
    g1 <- g1.u
    g2 <- g2.u
    #make new meta data frame
    meta.n <- data.frame(c(g1,g2))
    colnames(meta.n)[1] <- "SampleName"
    meta.n <- meta.n %>%
      mutate(Type = case_when(
        SampleName %in% g1 ~ input$g1label,
        SampleName %in% g2 ~ input$g2label
      )
      )
    meta.n
    
  })
  
  
  ####----Data Table----####
  
  
  output$UserDBmeta <- DT::renderDataTable({
    
    #load in files
    meta.u <- RNAmeta()
    meta.c <- meta_c()
    #Subset sample names for group 1 and 2 in each meta
    g1.u <- meta.u[which(meta.u[,2] == input$UserGroup1),1]
    g2.u <- meta.u[which(meta.u[,2] == input$UserGroup2),1]
    g1.c <- meta.c[which(meta.c[,2] == input$DBGroup1),1]
    g2.c <- meta.c[which(meta.c[,2] == input$DBGroup2),1]
    #group subsets between user and DB data
    g1 <- c(g1.u,g1.c)
    g2 <- c(g2.u,g2.c)
    #make new meta data frame
    meta.n <- data.frame(c(g1,g2))
    colnames(meta.n)[1] <- "SampleName"
    meta.n <- meta.n %>%
      mutate(Type = case_when(
        SampleName %in% g1 ~ input$g1label,
        SampleName %in% g2 ~ input$g2label
      )
      )
    #display table
    DT::datatable(meta.n,
                  options = list(keys = TRUE,
                                 searchHighlight = TRUE,
                                 pageLength = 25,
                                 scrollX = T,
                                 scrollY = T,
                                 lengthMenu = c("10", "25", "50", "100")))
    
    
    
    
  })
  
  output$UserDBmetaDL <- downloadHandler(
    
    filename = function() {
      if (input$CompMetaName == ""){
        paste("user_",ProjName,"_compiled_meta.tsv", sep = "")
      }
      else if (input$CompMetaName != ""){
        paste(input$CompMetaName,".tsv", sep = "")
      }
    },
    content = function(file){
      #load in files
      meta.u <- RNAmeta()
      meta.c <- meta_c()
      #Subset sample names for group 1 and 2 in each meta
      g1.u <- meta.u[which(meta.u[,2] == input$UserGroup1),1]
      g2.u <- meta.u[which(meta.u[,2] == input$UserGroup2),1]
      g1.c <- meta.c[which(meta.c[,2] == input$DBGroup1),1]
      g2.c <- meta.c[which(meta.c[,2] == input$DBGroup2),1]
      #group subsets between user and DB data
      g1 <- c(g1.u,g1.c)
      g2 <- c(g2.u,g2.c)
      #make new meta data frame
      meta.n <- data.frame(c(g1,g2))
      colnames(meta.n)[1] <- "SampleName"
      meta.n <- meta.n %>%
        mutate(Type = case_when(
          SampleName %in% g1 ~ "Group1",
          SampleName %in% g2 ~ "Group2"
        )
        )
      write_tsv(meta.n,file)
    }
    
  )
  
  
  output$UserExprTable <- DT::renderDataTable({
    
    DT::datatable(RNAmat(),
                  options = list(keys = TRUE,
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 scrollX = T,
                                 scrollY = T,
                                 lengthMenu = c("10", "25", "50", "100")))
    
  })
  
  output$UserMetaTable <- DT::renderDataTable({
    
    DT::datatable(RNAmeta(),
                  options = list(keys = TRUE,
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 scrollX = T,
                                 scrollY = T,
                                 lengthMenu = c("10", "25", "50", "100")))
    
  })
  
  #table for MSigDB gmt selection
  output$msigdbtable <- DT::renderDataTable({
    
    DT::datatable(MSigDBcat3,
                  selection = list(mode = 'single', selected = 1),
                  options = list(keys = TRUE,
                                 searchHighlight = TRUE,
                                 pageLength = 20,
                                 lengthMenu = c("10", "25", "50", "100")))
    
  })
  
  output$RNAvPROTtable <- DT::renderDataTable({
    
    req(input$RNAexp_F,input$RNAmeta_F)
    FCtable <- RNAvProtFC()
    M1 <- input$MAT1name
    M2 <- input$MAT2name
    colnames(FCtable)[c(2,3)] <- c(paste("Log2FC_",M1,sep = ""),paste("Log2FC_",M2,sep = ""))
    #get delta
    FCtable$Delta_M1M2 <- FCtable[,2] - FCtable[,3]
    colnames(FCtable)[4] <- paste("Delta",M1,M2,sep = "_")
    #make table
    DT::datatable(FCtable,
                  rownames = F,
                  options = list(keys = TRUE,
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = c("10", "25", "50", "100")))
    
  })
  
  #Leading edge genes table for matrix 1
  output$LeadingEdgeGenes1 <- DT::renderDataTable({
    
    res <- recGSEAg1()
    gmt <- recGMT()
    gsea.df <- as.data.frame(res@result)
    #use GS devised from matrix 2
    M2 <- input$MAT2name
    sets <- unique(gmt[,1])
    GS <- sets[grep(paste('_',M2,sep = ''),unique(gmt[,1]))]
    #make list object of LEG
    gsLEGlist1 <- list()
    for (i in GS) {
      genes1 <- as.matrix(gsea.df[which(gsea.df$Description==i),"core_enrichment"])
      genes2 <- unlist(strsplit(genes1,"/"))
      gsLEGlist1[[i]] <- genes2
    }
    #transform to data frame
    LEGtable1 <- as.data.frame(stri_list2matrix(gsLEGlist1))
    colnames(LEGtable1) <- names(gsLEGlist1)
    #add rank
    LEGtable1$Rank <- rownames(LEGtable1)
    #move rank to front
    LEGtable1 <- LEGtable1 %>%
      select(Rank, everything())
    DT::datatable(LEGtable1, options = list(paging = F), rownames = F)
    
  })
  
  #Leading edge genes table for matrix 2
  output$LeadingEdgeGenes2 <- DT::renderDataTable({
    
    res <- recGSEAg2()
    gmt <- recGMT()
    gsea.df <- as.data.frame(res@result)
    #use GS derived from matrix 1
    M1 <- input$MAT1name
    sets <- unique(gmt[,1])
    GS <- sets[grep(paste('_',M1,sep = ''),unique(gmt[,1]))]
    #make list object of LEG
    gsLEGlist2 <- list()
    for (i in GS) {
      genes1 <- as.matrix(gsea.df[which(gsea.df$Description==i),"core_enrichment"])
      genes2 <- unlist(strsplit(genes1,"/"))
      gsLEGlist2[[i]] <- genes2
    }
    #transform to data frame
    LEGtable2 <- as.data.frame(stri_list2matrix(gsLEGlist2))
    colnames(LEGtable2) <- names(gsLEGlist2)
    #add rank
    LEGtable2$Rank <- rownames(LEGtable2)
    #move rank to front
    LEGtable2 <- LEGtable2 %>%
      select(Rank, everything())
    DT::datatable(LEGtable2, options = list(paging = F), rownames = F)
    
  })
  
  #ssGSEA table - 1
  output$ssGSEAtable1 <- DT::renderDataTable({
    
    #matrix names
    M1 <- input$MAT1name
    M2 <- input$MAT2name
    #master table
    ssgseaM <- ssGSEAmaster()
    #select GS
    sets <- as.vector(unique(ssgseaM[,3]))
    GS <- sets[1]
    #reformat
    sscores <- ssgseaM[which(ssgseaM$variable == GS),]
    sscores <- sscores %>%
      select(Sample,type,Matrix,value)
    colnames(sscores)[2] <- "Type"
    colnames(sscores)[4] <- GS
    #table output
    DT::datatable(sscores,
                  rownames = F,
                  options = list(keys = TRUE,
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = c("10", "25", "50", "100"),
                                 scrollX = T))
    
  })
  
  #ssGSEA table - 2
  output$ssGSEAtable2 <- DT::renderDataTable({
    
    #matrix names
    M1 <- input$MAT1name
    M2 <- input$MAT2name
    #master table
    ssgseaM <- ssGSEAmaster()
    #select GS
    sets <- as.vector(unique(ssgseaM[,3]))
    GS <- sets[2]
    #reformat
    sscores <- ssgseaM[which(ssgseaM$variable == GS),]
    sscores <- sscores %>%
      select(Sample,type,Matrix,value)
    colnames(sscores)[2] <- "Type"
    colnames(sscores)[4] <- GS
    #table output
    DT::datatable(sscores,
                  rownames = F,
                  options = list(keys = TRUE,
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = c("10", "25", "50", "100"),
                                 scrollX = T))
  })
  
  #ssGSEA table - 3
  output$ssGSEAtable3 <- DT::renderDataTable({
    
    #matrix names
    M1 <- input$MAT1name
    M2 <- input$MAT2name
    #master table
    ssgseaM <- ssGSEAmaster()
    #select GS
    sets <- as.vector(unique(ssgseaM[,3]))
    GS <- sets[3]
    #reformat
    sscores <- ssgseaM[which(ssgseaM$variable == GS),]
    sscores <- sscores %>%
      select(Sample,type,Matrix,value)
    colnames(sscores)[2] <- "Type"
    colnames(sscores)[4] <- GS
    #table output
    DT::datatable(sscores,
                  rownames = F,
                  options = list(keys = TRUE,
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = c("10", "25", "50", "100"),
                                 scrollX = T))
  })
  
  #ssGSEA table - 4
  output$ssGSEAtable4 <- DT::renderDataTable({
    
    #matrix names
    M1 <- input$MAT1name
    M2 <- input$MAT2name
    #master table
    ssgseaM <- ssGSEAmaster()
    #select GS
    sets <- as.vector(unique(ssgseaM[,3]))
    GS <- sets[4]
    #reformat
    sscores <- ssgseaM[which(ssgseaM$variable == GS),]
    sscores <- sscores %>%
      select(Sample,type,Matrix,value)
    colnames(sscores)[2] <- "Type"
    colnames(sscores)[4] <- GS
    #table output
    DT::datatable(sscores,
                  rownames = F,
                  options = list(keys = TRUE,
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = c("10", "25", "50", "100"),
                                 scrollX = T))
  })
  
  #Up and Down stats table
  output$VennStatTable <- DT::renderDataTable({
    
    if (input$genesets == 1) {
      
      #Assign matrices to variables and assign meta groups
      M1 <- RNAmat()
      M2 <- expr_c()
      meta1 <- meta_un()
      meta2 <- meta_cn()
      #metasame <- merge(meta1,meta2)
      A1 <- meta1[which(meta1[,2] == input$groupAcomp3),1]
      B1 <- meta1[which(meta1[,2] == input$groupBcomp3),1]
      A2 <- meta2[which(meta2[,2] == input$groupAcomp3),1]
      B2 <- meta2[which(meta2[,2] == input$groupBcomp3),1]
      
      #Assign user given matrix names for labeling
      M1Name <- input$MAT1name
      M2Name <- input$MAT2name
      
      #Matrix 1 top table generation
      mat <- M1[,c(A1,B1)]
      mat <- log2(mat + 1.0)
      groupAOther <- factor(c(rep("A", length(A1)), rep("B", length(B1))))
      designA <- model.matrix(~0 + groupAOther)
      fit <- lmFit(mat, design = designA)
      contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
      fit2 <- contrasts.fit(fit, contrast.matrix)
      fit2 <- eBayes(fit2)
      options(digits = 4)
      top1_1 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
      
      #Matrix 2 top table generation
      mat <- M2[,c(A2,B2)]
      mat <- log2(mat + 1.0)
      groupAOther <- factor(c(rep("A", length(A2)), rep("B", length(B2))))
      designA <- model.matrix(~0 + groupAOther)
      fit <- lmFit(mat, design = designA)
      contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
      fit2 <- contrasts.fit(fit, contrast.matrix)
      fit2 <- eBayes(fit2)
      options(digits = 4)
      top1_2 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
      
      #re-assign main files and variables
      M1top <- top1_1
      M2top <- top1_2
      
      #user selected GS category
      GScat <- MSigDBcat3[input$msigdbtable_rows_selected,2]
      #select GS names based on category selected
      GScat_sets <- MSigDBcat[which(MSigDBcat[,2] == GScat),3]
      #add GS to GSlist based off of what is selected
      gs_list <- gs[GScat_sets]
      gmt_unique <- unique(unlist(gs_list, use.names = F))
      
      #list of similar genes between 3 datasets
      SimGenes <- intersect(intersect(rownames(M1top),rownames(M2top)),gmt_unique)
      #subset data based on similar genes only
      M1top2 <- M1top[SimGenes,] 
      M2top2 <- M2top[SimGenes,]
      #extract genes from gslist that are not in matrices
      for (i in 1:length(gs_list)){
        gs_list[[i]] <- gs_list[[i]][which(gs_list[[i]] %in% SimGenes)] 
      }
      #list of unique but similar to matrices genes in GMT
      gmt2Uniq <- unique(unlist(gs_list, use.names = F))
      
      #user variables
      FC <- input$vennLog
      adjP <- input$adjpvalVenn
      
      ##get gene lists from top tables
      #subset up reg DEGs and non DEGs - M1
      DEG_M1U <- rownames(M1top2[which(M1top2$logFC > FC & M1top2$adj.P.Val <= adjP),])
      nonDEG_M1U <- rownames(M1top2[which(!(M1top2$logFC > FC & M1top2$adj.P.Val <= adjP)),])
      #subset up reg DEGs and non DEGs - M2
      DEG_M2U <- rownames(M2top2[which(M2top2$logFC > FC & M2top2$adj.P.Val <= adjP),])
      nonDEG_M2U <- rownames(M2top2[which(!(M2top2$logFC > FC & M2top2$adj.P.Val <= adjP)),])
      #subset down reg DEGs and non DEGs - M1
      DEG_M1D <- rownames(M1top2[which(M1top2$logFC < FC & M1top2$adj.P.Val <= adjP),])
      nonDEG_M1D <- rownames(M1top2[which(!(M1top2$logFC < FC & M1top2$adj.P.Val <= adjP)),])
      #subset down reg DEGs and non DEGs - M2
      DEG_M2D <- rownames(M2top2[which(M2top2$logFC < FC & M2top2$adj.P.Val <= adjP),])
      nonDEG_M2D <- rownames(M2top2[which(!(M2top2$logFC < FC & M2top2$adj.P.Val <= adjP)),])
      
      #Generate names for the DEG lists made above - will be used to compare with GMT
      gsName1 <- paste("UpRegDEG_in_",M1Name,sep = "")
      gsName2 <- paste("DnRegDEG_in_",M1Name,sep = "")
      gsName3 <- paste("UpRegDEG_in_",M2Name,sep = "")
      gsName4 <- paste("DnRegDEG_in_",M2Name,sep = "")
      
      #Add DEGs subset above to GMT list
      gsDataListGMT <- list() #start gene set list
      gsDataListGMT[[gsName1]] <- DEG_M1U
      gsDataListGMT[[gsName3]] <- DEG_M2U
      gsDataListGMT[[gsName2]] <- DEG_M1D
      gsDataListGMT[[gsName4]] <- DEG_M2D
      
      gsDataListGMT <- c(gsDataListGMT,gs_list)
      
      #world of all genes from matrices
      DEG_M1M2 <- intersect(rownames(M1top2), rownames(M2top2))
      
      #learn Jaccard index function
      jaccard <- function(a, b) {
        intersection = length(intersect(a, b))
        union = length(a) + length(b) - intersection
        return (intersection/union)
      }
      
      #list to store stat data
      StatDataList <- list()
      
      #loop through gmt file and get stats when comparing gene sets to expression matrices
      for (i in 1:length(gsDataListGMT)) {
        
        #Gene set name that is being compared
        GS <- names(gsDataListGMT)[[i]]
        
        #Integration of UP regulated genes with gmt gene set - M1
        M1DEG_UP_GS <- length(intersect(DEG_M1U,gsDataListGMT[[i]]))                        #DE up + in GS
        M1nDEG_UP_nGS <- length(intersect(nonDEG_M1U,setdiff(gmt2Uniq,gsDataListGMT[[i]]))) #non DE up + not in GS 
        M1DEG_UP_nGS <- length(intersect(DEG_M1U,setdiff(gmt2Uniq,gsDataListGMT[[i]])))     #DE up + not in GS
        M1nDEG_UP_GS <- length(intersect(nonDEG_M1U,gsDataListGMT[[i]]))                    #non DE up + in GS
        #Integration of DOWN regulated genes with gmt gene set - M1
        M1DEG_DN_GS <- length(intersect(DEG_M1D,gsDataListGMT[[i]]))                        #DE dn + in GS
        M1nDEG_DN_nGS <- length(intersect(nonDEG_M1D,setdiff(gmt2Uniq,gsDataListGMT[[i]]))) #non DE dn + not in GS
        M1DEG_DN_nGS <- length(intersect(DEG_M1D,setdiff(gmt2Uniq,gsDataListGMT[[i]])))     #DE dn + not in GS
        M1nDEG_DN_GS <- length(intersect(nonDEG_M1D,gsDataListGMT[[i]]))                    #non DE dn + in GS
        
        #Integration of UP regulated genes with gmt gene set - M2
        M2DEG_UP_GS <- length(intersect(DEG_M2U,gsDataListGMT[[i]]))
        M2nDEG_UP_nGS <- length(intersect(nonDEG_M2U,setdiff(gmt2Uniq,gsDataListGMT[[i]])))
        M2DEG_UP_nGS <- length(intersect(DEG_M2U,setdiff(gmt2Uniq,gsDataListGMT[[i]])))
        M2nDEG_UP_GS <- length(intersect(nonDEG_M2U,gsDataListGMT[[i]]))
        #Integration of DOWN regulated genes with gmt gene set - M2
        M2DEG_DN_GS <- length(intersect(DEG_M2D,gsDataListGMT[[i]]))
        M2nDEG_DN_nGS <- length(intersect(nonDEG_M2D,setdiff(gmt2Uniq,gsDataListGMT[[i]])))
        M2DEG_DN_nGS <- length(intersect(DEG_M2D,setdiff(gmt2Uniq,gsDataListGMT[[i]])))
        M2nDEG_DN_GS <- length(intersect(nonDEG_M2D,gsDataListGMT[[i]]))
        
        #number of genes in 'world' AND GS
        WorldNGS <- length(intersect(DEG_M1M2,gsDataListGMT[[i]]))
        
        ##Generate Matrices
        #up regulated - M1
        UpMat1GS <- matrix(data = c(M1DEG_UP_GS,M1DEG_UP_nGS,
                                    M1nDEG_UP_GS,M1nDEG_UP_nGS), nrow = 2)
        #write out numbers of matrix for last columns
        UpMat1GS.n <- paste(M1DEG_UP_GS,M1DEG_UP_nGS,
                            M1nDEG_UP_GS,M1nDEG_UP_nGS, sep = ",")
        
        #down regulated - M1
        DnMat1GS <- matrix(data = c(M1DEG_DN_GS,M1DEG_DN_nGS,
                                    M1nDEG_DN_GS,M1nDEG_DN_nGS), nrow = 2)
        #write out numbers of matrix for last columns
        DnMat1GS.n <- paste(M1DEG_DN_GS,M1DEG_DN_nGS,
                            M1nDEG_DN_GS,M1nDEG_DN_nGS, sep = ",")
        
        #up regulated - M2
        UpMat2GS <- matrix(data = c(M2DEG_UP_GS,M2DEG_UP_nGS,
                                    M2nDEG_UP_GS,M2nDEG_UP_nGS), nrow = 2)
        #write out numbers of matrix for last columns
        UpMat2GS.n <- paste(M2DEG_UP_GS,M2DEG_UP_nGS,
                            M2nDEG_UP_GS,M2nDEG_UP_nGS, sep = ",")
        
        #down regulated - M2
        DnMat2GS <- matrix(data = c(M2DEG_DN_GS,M2DEG_DN_nGS,
                                    M2nDEG_DN_GS,M2nDEG_DN_nGS), nrow = 2)
        #write out numbers of matrix for last columns
        DnMat2GS.n <- paste(M2DEG_DN_GS,M2DEG_DN_nGS,
                            M2nDEG_DN_GS,M2nDEG_DN_nGS, sep = ",")
        
        ##Fisher Exact Test
        #Up regulated - M1
        UpfishM1 <- fisher.test(UpMat1GS, alternative = "two.sided")
        UpFishPM1 <- UpfishM1[[1]]
        UpFishORM1 <- UpfishM1[[3]]
        #Down regulated - M1
        DnfishM1 <- fisher.test(DnMat1GS, alternative = "two.sided")
        DnFishPM1 <- DnfishM1[[1]]
        DnFishORM1 <- DnfishM1[[3]]
        #Up regulated - M2
        UpfishM2 <- fisher.test(UpMat2GS, alternative = "two.sided")
        UpFishPM2 <- UpfishM2[[1]]
        UpFishORM2 <- UpfishM2[[3]]
        #Down regulated - M2
        DnfishM2 <- fisher.test(DnMat2GS, alternative = "two.sided")
        DnFishPM2 <- DnfishM2[[1]]
        DnFishORM2 <- DnfishM2[[3]]
        
        ##Kappa Score
        #Mat1
        UpKappaM1 <- kappa(UpMat1GS)
        DnKappaM1 <- kappa(DnMat1GS)
        #Mat1
        UpKappaM2 <- kappa(UpMat2GS)
        DnKappaM2 <- kappa(DnMat2GS)
        
        #up reg genes in M1
        UpJaccM1 <- jaccard(DEG_M1U,gsDataListGMT[[i]])
        DnJaccM1 <- jaccard(DEG_M1D,gsDataListGMT[[i]])
        #up reg genes in M2
        UpJaccM2 <- jaccard(DEG_M2U,gsDataListGMT[[i]])
        DnJaccM2 <- jaccard(DEG_M2D,gsDataListGMT[[i]])
        
        #Assign stats to GS name and add to list
        StatDataList[[GS]] <- c(WorldNGS,
                                UpFishPM1,UpFishPM2,DnFishPM1,DnFishPM2,
                                UpFishORM1,UpFishORM2,DnFishORM1,DnFishORM2,
                                UpKappaM1,UpKappaM2,DnKappaM1,DnKappaM2,
                                UpJaccM1,UpJaccM2,DnJaccM1,DnJaccM2,
                                UpMat1GS.n,UpMat2GS.n,DnMat1GS.n,DnMat2GS.n)
        
      }
      
      #convert to data frame
      StatData.df <- do.call(rbind.data.frame, StatDataList)
      rownames(StatData.df) <- names(StatDataList)
      
      #colnames(StatData.df) <- c("Number_of_DEGenes_in_GeneSet",
      #                           paste(M1Name,"_UpRegulated_Pvalue",sep = ""),paste(M2Name,"_UpRegulated_Pvalue",sep = ""),
      #                           paste(M1Name,"_DownRegulated_Pvalue",sep = ""),paste(M2Name,"_DownRegulated_Pvalue",sep = ""),
      #                           paste(M1Name,"_UpRegulated_OddsRatio",sep = ""),paste(M2Name,"_UpRegulated_OddsRatio",sep = ""),
      #                           paste(M1Name,"_DownRegulated_OddsRatio",sep = ""),paste(M2Name,"_DownRegulated_OddsRatio",sep = ""),
      #                           paste(M1Name,"_UpRegulated_KappaScore",sep = ""),paste(M2Name,"_UpRegulated_KappaScore",sep = ""),
      #                           paste(M1Name,"_DownRegulated_KappaScore",sep = ""),paste(M2Name,"_DownRegulated_KappaScore",sep = ""),
      #                           paste(M1Name,"_UpRegulated_JaccardIndex",sep = ""),paste(M2Name,"_UpRegulated_JaccardIndex",sep = ""),
      #                           paste(M1Name,"_DownRegulated_JaccardIndex",sep = ""),paste(M2Name,"_DownRegulated_JaccardIndex",sep = ""),
      #                           paste(M1Name,"_UpRegulated_Matrix",sep = ""),paste(M2Name,"_UpRegulated_Matrix",sep = ""),
      #                           paste(M1Name,"_DownRegulated_Matrix",sep = ""),paste(M2Name,"_DownRegulated_Matrix",sep = ""))
      #spaces in the name allow for them to be wrapped in the column names
      colnames(StatData.df) <- c("Number of DEGs in GeneSet",
                                 paste(M1Name," UpRegulated Pvalue",sep = ""),paste(M2Name," UpRegulated Pvalue",sep = ""),
                                 paste(M1Name," DownRegulated Pvalue",sep = ""),paste(M2Name," DownRegulated Pvalue",sep = ""),
                                 paste(M1Name," UpRegulated OddsRatio",sep = ""),paste(M2Name," UpRegulated OddsRatio",sep = ""),
                                 paste(M1Name," DownRegulated OddsRatio",sep = ""),paste(M2Name," DownRegulated OddsRatio",sep = ""),
                                 paste(M1Name," UpRegulated KappaScore",sep = ""),paste(M2Name," UpRegulated KappaScore",sep = ""),
                                 paste(M1Name," DownRegulated KappaScore",sep = ""),paste(M2Name," DownRegulated KappaScore",sep = ""),
                                 paste(M1Name," UpRegulated JaccardIndex",sep = ""),paste(M2Name," UpRegulated JaccardIndex",sep = ""),
                                 paste(M1Name," DownRegulated JaccardIndex",sep = ""),paste(M2Name," DownRegulated JaccardIndex",sep = ""),
                                 paste(M1Name," UpRegulated Contingency Table",sep = ""),paste(M2Name," UpRegulated Contingency Table",sep = ""),
                                 paste(M1Name," DownRegulated Contingency Table",sep = ""),paste(M2Name," DownRegulated Contingency Table",sep = ""))
      
      #make labeled column for rownames
      StatData.df$GeneSets <- rownames(StatData.df)
      StatData.df <- StatData.df %>%
        relocate(GeneSets)
      #convert numeric column to as.numeric
      StatData.df[,c(2:18)] <- sapply(StatData.df[,c(2:18)], as.numeric)
      #table output
      DT::datatable(StatData.df,
                    rownames = F,
                    extensions = "FixedColumns",
                    options = list(keys = TRUE,
                                   searchHighlight = TRUE,
                                   pageLength = 10,
                                   lengthMenu = c("10", "25", "50", "100"),
                                   scrollX = T,
                                   autoWidth = TRUE,
                                   fixedColumns = list(leftColumns = 1))) %>%
        formatRound(columns=c(3:18), digits=4)
    }
    else if (input$genesets == 2) {
      
      #user upload of gmt file
      gmt.u <- input$GMTcomp
      ext <- tools::file_ext(gmt.u$datapath)
      req(gmt.u)
      validate(need(ext == c("gmt"), "Please upload a .gmt file"))
      gmt <- read.gmt(gmt.u$datapath)
      
      #Assign matrices to variables and assign meta groups
      M1 <- RNAmat()
      M2 <- expr_c()
      meta1 <- meta_un()
      meta2 <- meta_cn()
      #metasame <- merge(meta1,meta2)
      A1 <- meta1[which(meta1[,2] == input$groupAcomp3),1]
      B1 <- meta1[which(meta1[,2] == input$groupBcomp3),1]
      A2 <- meta2[which(meta2[,2] == input$groupAcomp3),1]
      B2 <- meta2[which(meta2[,2] == input$groupBcomp3),1]
      
      #Assign user given matrix names for labeling
      M1Name <- input$MAT1name
      M2Name <- input$MAT2name
      
      #Matrix 1 top table generation
      mat <- M1[,c(A1,B1)]
      mat <- log2(mat + 1.0)
      groupAOther <- factor(c(rep("A", length(A1)), rep("B", length(B1))))
      designA <- model.matrix(~0 + groupAOther)
      fit <- lmFit(mat, design = designA)
      contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
      fit2 <- contrasts.fit(fit, contrast.matrix)
      fit2 <- eBayes(fit2)
      options(digits = 4)
      top1_1 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
      
      #Matrix 2 top table generation
      mat <- M2[,c(A2,B2)]
      mat <- log2(mat + 1.0)
      groupAOther <- factor(c(rep("A", length(A2)), rep("B", length(B2))))
      designA <- model.matrix(~0 + groupAOther)
      fit <- lmFit(mat, design = designA)
      contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
      fit2 <- contrasts.fit(fit, contrast.matrix)
      fit2 <- eBayes(fit2)
      options(digits = 4)
      top1_2 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
      
      #re-assign main files and variables
      M1top <- top1_1
      M2top <- top1_2
      
      #list of similar genes between 3 datasets
      SimGenes <- intersect(intersect(rownames(M1top),rownames(M2top)),gmt[,2])
      #subset data based on similar genes only
      M1top2 <- M1top[SimGenes,] 
      M2top2 <- M2top[SimGenes,]
      gmt2 <- gmt[which(gmt[,2] %in% SimGenes),]
      #list of unique bu similar to matrices genes in GMT
      gmt2Uniq <- unique(gmt2[,2])
      
      #user variables
      FC <- input$vennLog
      adjP <- input$adjpvalVenn
      
      ##get gene lists from top tables
      #subset up reg DEGs and non DEGs - M1
      DEG_M1U <- rownames(M1top2[which(M1top2$logFC > FC & M1top2$adj.P.Val <= adjP),])
      nonDEG_M1U <- rownames(M1top2[which(!(M1top2$logFC > FC & M1top2$adj.P.Val <= adjP)),])
      #subset up reg DEGs and non DEGs - M2
      DEG_M2U <- rownames(M2top2[which(M2top2$logFC > FC & M2top2$adj.P.Val <= adjP),])
      nonDEG_M2U <- rownames(M2top2[which(!(M2top2$logFC > FC & M2top2$adj.P.Val <= adjP)),])
      #subset down reg DEGs and non DEGs - M1
      DEG_M1D <- rownames(M1top2[which(M1top2$logFC < FC & M1top2$adj.P.Val <= adjP),])
      nonDEG_M1D <- rownames(M1top2[which(!(M1top2$logFC < FC & M1top2$adj.P.Val <= adjP)),])
      #subset down reg DEGs and non DEGs - M2
      DEG_M2D <- rownames(M2top2[which(M2top2$logFC < FC & M2top2$adj.P.Val <= adjP),])
      nonDEG_M2D <- rownames(M2top2[which(!(M2top2$logFC < FC & M2top2$adj.P.Val <= adjP)),])
      
      #Generate names for the DEG lists made above - will be used to compare with GMT
      gsName1 <- paste("UpRegDEG_in_",M1Name,sep = "")
      gsName2 <- paste("DnRegDEG_in_",M1Name,sep = "")
      gsName3 <- paste("UpRegDEG_in_",M2Name,sep = "")
      gsName4 <- paste("DnRegDEG_in_",M2Name,sep = "")
      
      #Add DEGs subset above to GMT list
      gsDataListGMT <- list() #start gene set list
      gsDataListGMT[[gsName1]] <- DEG_M1U
      gsDataListGMT[[gsName3]] <- DEG_M2U
      gsDataListGMT[[gsName2]] <- DEG_M1D
      gsDataListGMT[[gsName4]] <- DEG_M2D
      
      #world of all genes from matrices
      DEG_M1M2 <- intersect(rownames(M1top2), rownames(M2top2))
      
      #make gene set list from GMT
      for (i in unique(gmt2[,1])){
        gsDataListGMT[[i]] <- gmt2[gmt2[,1] == i,]$gene
      }
      
      #learn Jaccard index function
      jaccard <- function(a, b) {
        intersection = length(intersect(a, b))
        union = length(a) + length(b) - intersection
        return (intersection/union)
      }
      
      #list to store stat data
      StatDataList <- list()
      
      #loop through gmt file and get stats when comparing gene sets to expression matrices
      for (i in 1:length(gsDataListGMT)) {
        
        #Gene set name that is being compared
        GS <- names(gsDataListGMT)[[i]]
        
        #Integration of UP regulated genes with gmt gene set - M1
        M1DEG_UP_GS <- length(intersect(DEG_M1U,gsDataListGMT[[i]]))                        #DE up + in GS
        M1nDEG_UP_nGS <- length(intersect(nonDEG_M1U,setdiff(gmt2Uniq,gsDataListGMT[[i]]))) #non DE up + not in GS 
        M1DEG_UP_nGS <- length(intersect(DEG_M1U,setdiff(gmt2Uniq,gsDataListGMT[[i]])))     #DE up + not in GS
        M1nDEG_UP_GS <- length(intersect(nonDEG_M1U,gsDataListGMT[[i]]))                    #non DE up + in GS
        #Integration of DOWN regulated genes with gmt gene set - M1
        M1DEG_DN_GS <- length(intersect(DEG_M1D,gsDataListGMT[[i]]))                        #DE dn + in GS
        M1nDEG_DN_nGS <- length(intersect(nonDEG_M1D,setdiff(gmt2Uniq,gsDataListGMT[[i]]))) #non DE dn + not in GS
        M1DEG_DN_nGS <- length(intersect(DEG_M1D,setdiff(gmt2Uniq,gsDataListGMT[[i]])))     #DE dn + not in GS
        M1nDEG_DN_GS <- length(intersect(nonDEG_M1D,gsDataListGMT[[i]]))                    #non DE dn + in GS
        
        #Integration of UP regulated genes with gmt gene set - M2
        M2DEG_UP_GS <- length(intersect(DEG_M2U,gsDataListGMT[[i]]))
        M2nDEG_UP_nGS <- length(intersect(nonDEG_M2U,setdiff(gmt2Uniq,gsDataListGMT[[i]])))
        M2DEG_UP_nGS <- length(intersect(DEG_M2U,setdiff(gmt2Uniq,gsDataListGMT[[i]])))
        M2nDEG_UP_GS <- length(intersect(nonDEG_M2U,gsDataListGMT[[i]]))
        #Integration of DOWN regulated genes with gmt gene set - M2
        M2DEG_DN_GS <- length(intersect(DEG_M2D,gsDataListGMT[[i]]))
        M2nDEG_DN_nGS <- length(intersect(nonDEG_M2D,setdiff(gmt2Uniq,gsDataListGMT[[i]])))
        M2DEG_DN_nGS <- length(intersect(DEG_M2D,setdiff(gmt2Uniq,gsDataListGMT[[i]])))
        M2nDEG_DN_GS <- length(intersect(nonDEG_M2D,gsDataListGMT[[i]]))
        
        #number of genes in 'world' AND GS
        WorldNGS <- length(intersect(DEG_M1M2,gsDataListGMT[[i]]))
        
        ##Generate Matrices
        #up regulated - M1
        UpMat1GS <- matrix(data = c(M1DEG_UP_GS,M1DEG_UP_nGS,
                                    M1nDEG_UP_GS,M1nDEG_UP_nGS), nrow = 2)
        #write out numbers of matrix for last columns
        UpMat1GS.n <- paste(M1DEG_UP_GS,M1DEG_UP_nGS,
                            M1nDEG_UP_GS,M1nDEG_UP_nGS, sep = ",")
        
        #down regulated - M1
        DnMat1GS <- matrix(data = c(M1DEG_DN_GS,M1DEG_DN_nGS,
                                    M1nDEG_DN_GS,M1nDEG_DN_nGS), nrow = 2)
        #write out numbers of matrix for last columns
        DnMat1GS.n <- paste(M1DEG_DN_GS,M1DEG_DN_nGS,
                            M1nDEG_DN_GS,M1nDEG_DN_nGS, sep = ",")
        
        #up regulated - M2
        UpMat2GS <- matrix(data = c(M2DEG_UP_GS,M2DEG_UP_nGS,
                                    M2nDEG_UP_GS,M2nDEG_UP_nGS), nrow = 2)
        #write out numbers of matrix for last columns
        UpMat2GS.n <- paste(M2DEG_UP_GS,M2DEG_UP_nGS,
                            M2nDEG_UP_GS,M2nDEG_UP_nGS, sep = ",")
        
        #down regulated - M2
        DnMat2GS <- matrix(data = c(M2DEG_DN_GS,M2DEG_DN_nGS,
                                    M2nDEG_DN_GS,M2nDEG_DN_nGS), nrow = 2)
        #write out numbers of matrix for last columns
        DnMat2GS.n <- paste(M2DEG_DN_GS,M2DEG_DN_nGS,
                            M2nDEG_DN_GS,M2nDEG_DN_nGS, sep = ",")
        
        ##Fisher Exact Test
        #Up regulated - M1
        UpfishM1 <- fisher.test(UpMat1GS, alternative = "two.sided")
        UpFishPM1 <- UpfishM1[[1]]
        UpFishORM1 <- UpfishM1[[3]]
        #Down regulated - M1
        DnfishM1 <- fisher.test(DnMat1GS, alternative = "two.sided")
        DnFishPM1 <- DnfishM1[[1]]
        DnFishORM1 <- DnfishM1[[3]]
        #Up regulated - M2
        UpfishM2 <- fisher.test(UpMat2GS, alternative = "two.sided")
        UpFishPM2 <- UpfishM2[[1]]
        UpFishORM2 <- UpfishM2[[3]]
        #Down regulated - M2
        DnfishM2 <- fisher.test(DnMat2GS, alternative = "two.sided")
        DnFishPM2 <- DnfishM2[[1]]
        DnFishORM2 <- DnfishM2[[3]]
        
        ##Kappa Score
        #Mat1
        UpKappaM1 <- kappa(UpMat1GS)
        DnKappaM1 <- kappa(DnMat1GS)
        #Mat1
        UpKappaM2 <- kappa(UpMat2GS)
        DnKappaM2 <- kappa(DnMat2GS)
        
        #up reg genes in M1
        UpJaccM1 <- jaccard(DEG_M1U,gsDataListGMT[[i]])
        DnJaccM1 <- jaccard(DEG_M1D,gsDataListGMT[[i]])
        #up reg genes in M2
        UpJaccM2 <- jaccard(DEG_M2U,gsDataListGMT[[i]])
        DnJaccM2 <- jaccard(DEG_M2D,gsDataListGMT[[i]])
        
        #Assign stats to GS name and add to list
        StatDataList[[GS]] <- c(WorldNGS,
                                UpFishPM1,UpFishPM2,DnFishPM1,DnFishPM2,
                                UpFishORM1,UpFishORM2,DnFishORM1,DnFishORM2,
                                UpKappaM1,UpKappaM2,DnKappaM1,DnKappaM2,
                                UpJaccM1,UpJaccM2,DnJaccM1,DnJaccM2,
                                UpMat1GS.n,UpMat2GS.n,DnMat1GS.n,DnMat2GS.n)
        
      }
      
      #convert to data frame
      StatData.df <- do.call(rbind.data.frame, StatDataList)
      rownames(StatData.df) <- names(StatDataList)
      #spaces in the name allow for them to be wrapped in the column names
      colnames(StatData.df) <- c("Number of DEGs in GeneSet",
                                 paste(M1Name," UpRegulated Pvalue",sep = ""),paste(M2Name," UpRegulated Pvalue",sep = ""),
                                 paste(M1Name," DownRegulated Pvalue",sep = ""),paste(M2Name," DownRegulated Pvalue",sep = ""),
                                 paste(M1Name," UpRegulated OddsRatio",sep = ""),paste(M2Name," UpRegulated OddsRatio",sep = ""),
                                 paste(M1Name," DownRegulated OddsRatio",sep = ""),paste(M2Name," DownRegulated OddsRatio",sep = ""),
                                 paste(M1Name," UpRegulated KappaScore",sep = ""),paste(M2Name," UpRegulated KappaScore",sep = ""),
                                 paste(M1Name," DownRegulated KappaScore",sep = ""),paste(M2Name," DownRegulated KappaScore",sep = ""),
                                 paste(M1Name," UpRegulated JaccardIndex",sep = ""),paste(M2Name," UpRegulated JaccardIndex",sep = ""),
                                 paste(M1Name," DownRegulated JaccardIndex",sep = ""),paste(M2Name," DownRegulated JaccardIndex",sep = ""),
                                 paste(M1Name," UpRegulated Contingency Table",sep = ""),paste(M2Name," UpRegulated Contingency Table",sep = ""),
                                 paste(M1Name," DownRegulated Contingency Table",sep = ""),paste(M2Name," DownRegulated Contingency Table",sep = ""))
      #make labeled column for rownames
      StatData.df$GeneSets <- rownames(StatData.df)
      StatData.df <- StatData.df %>%
        relocate(GeneSets)
      #convert numeric column to as.numeric
      StatData.df[,c(2:18)] <- sapply(StatData.df[,c(2:18)], as.numeric)
      #table output
      DT::datatable(StatData.df,
                    rownames = F,
                    extensions = "FixedColumns",
                    options = list(keys = TRUE,
                                   searchHighlight = TRUE,
                                   pageLength = 10,
                                   lengthMenu = c("10", "25", "50", "100"),
                                   scrollX = T,
                                   autoWidth = TRUE,
                                   fixedColumns = list(leftColumns = 1))) %>%
        formatRound(columns=c(3:18), digits=4)
      
    }
    
  })
  
  
  ####----Plots----####
  
  
  #render RNAseq vs Protein scatter plot
  output$RNAvPROTscatter <- renderPlot({
    
    req(input$RNAexp_F,input$RNAmeta_F)
    FCtable <- RNAvProtFC()
    if (length(input$scatterGeneSelec) >= 1 || input$gsSelection2 != "") {
      
      genesel.s <- NULL
      genesel.t <- NULL
      highlight2 <- NULL
      highlight1 <- NULL
      gene.selec <- NULL
      gene.selec <- input$scatterGeneSelec
      genesel.s <- unlist(strsplit(input$gsSelection2, " "))
      genesel.t <- unlist(strsplit(input$gsSelection2, "\t"))
      genesel.text <- c(genesel.s,genesel.t)
      highlight1 <- FCtable %>%
        filter(gene %in% gene.selec)
      highlight2 <- FCtable %>%
        filter(gene %in% genesel.text)
      highlight <- rbind(highlight1,highlight2)
      #edit column names for hover text to find
      colnames(FCtable) <- c("gene","log2FC_M1","log2FC_M2")
      #generate plot
      p <- ggplot(FCtable, aes(x = log2FC_M1, y = log2FC_M2)) +
        geom_point(color = "lightgray") +
        theme_minimal() +
        labs(x = paste("log2FC_", input$MAT1name, sep = ""), y = paste("log2FC_", input$MAT2name, sep = ""))
      #Color selected points
      p <- p +
        geom_point(data = highlight,
                   aes(x = highlight[,2], y = highlight[,3]),
                   color = "darkred",
                   size = 2)
      p <- p +
        geom_text_repel(data = highlight,
                        aes(label = highlight[,1], x = highlight[,2], y = highlight[,3]),
                        size = 6,
                        color="black",
                        nudge_x = 0.1,
                        nudge_y=0.1,
                        box.padding = unit(0.9, "lines"),
                        point.padding = unit(.3+4*0.1, "lines"),
                        max.overlaps = 50)
      #show plot
      p
      
    }
    else {
      #edit column names for hover text to find
      colnames(FCtable) <- c("gene","log2FC_M1","log2FC_M2")
      p <- ggplot(FCtable, aes(x = log2FC_M1, y = log2FC_M2)) +
        geom_point(color = "lightgray") +
        theme_minimal() +
        labs(x = paste("log2FC_", input$MAT1name, sep = ""), y = paste("log2FC_", input$MAT2name, sep = ""))
      #show plot
      p
    }
    
  })
  
  #Render GSEA Enrichment plot
  output$enrichplot1 <- renderPlot({
    
    res <- recGSEAg1()
    gmt <- recGMT()
    gsea.df <- as.data.frame(res@result)
    #use GS derived from matrix 2
    M2 <- input$MAT2name
    sets <- unique(gmt[,1])
    geneset <- sets[grep(paste('_',M2,sep = ''),unique(gmt[,1]))]
    gseaplot2(res,
              geneset,
              pvalue_table = F)
    
  })
  
  #Render GSEA Enrichment plot
  output$enrichplot2 <- renderPlot({
    
    res <- recGSEAg2()
    gmt <- recGMT()
    gsea.df <- as.data.frame(res@result)
    #use GS derived from matrix 1
    M1 <- input$MAT1name
    sets <- unique(gmt[,1])
    geneset <- sets[grep(paste('_',M1,sep = ''),unique(gmt[,1]))]
    gseaplot2(res,
              geneset,
              pvalue_table = F)
    
  })
  
  #render up reg venn diagram
  output$vennUP <- renderPlot({
    
    #Variables
    RNA <- RNAmat()
    prot <- expr_c()
    metaR <- meta_un()
    metaP <- meta_cn()
    metasame <- intersect(metaR[,2],metaP[,2])
    metagroups <- as.vector(levels(factor(metasame)))
    #metagroupsR <- as.vector(levels(factor(metaR[,2])))
    #metagroupsP <- as.vector(levels(factor(metaP[,2])))
    A1 <- metaR[which(metaR[,2] == input$groupAcomp3),1]
    B1 <- metaR[which(metaR[,2] == input$groupBcomp3),1]
    A2 <- metaP[which(metaP[,2] == input$groupAcomp3),1]
    B2 <- metaP[which(metaP[,2] == input$groupBcomp3),1]
    #matrix names
    M1 <- input$MAT1name
    M2 <- input$MAT2name
    
    #RNAseq top table generation
    mat <- RNA[,c(A1,B1)]
    mat <- log2(mat + 1.0)
    groupAOther <- factor(c(rep("A", length(A1)), rep("B", length(B1))))
    designA <- model.matrix(~0 + groupAOther)
    fit <- lmFit(mat, design = designA)
    contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    options(digits = 4)
    top1_R <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
    
    #Protein top table generation
    mat <- prot[,c(A2,B2)]
    mat <- log2(mat + 1.0)
    groupAOther <- factor(c(rep("A", length(A2)), rep("B", length(B2))))
    designA <- model.matrix(~0 + groupAOther)
    fit <- lmFit(mat, design = designA)
    contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    options(digits = 4)
    top1_P <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
    
    #user variables
    FC <- input$vennLog
    adjP <- input$adjpvalVenn
    
    #subset Up reg top tables
    topGenes_RU <- top1_R[top1_R$logFC > FC & top1_R$adj.P.Val <= adjP,]
    topGenes_PU <- top1_P[top1_P$logFC > FC & top1_P$adj.P.Val <= adjP,]
    
    #convert genes into gmt format
    topGenes_RU_gmt <- data.frame(paste("UPreg_in_",gsub(" ","",metagroups[1]),"_",M1,sep = ""),rownames(topGenes_RU))
    colnames(topGenes_RU_gmt) <- c("term","gene")
    topGenes_PU_gmt <- data.frame(paste("UPreg_in_",gsub(" ","",metagroups[1]),"_",M2,sep = ""),rownames(topGenes_PU))
    colnames(topGenes_PU_gmt) <- c("term","gene")
    
    #combine to 4 gene set list
    RNAvProt_GS <- rbind(topGenes_RU_gmt,topGenes_PU_gmt)
    
    #make R Data List
    gsDataList <- list()
    for (i in unique(RNAvProt_GS[,1])){
      gsDataList[[i]] <- RNAvProt_GS[RNAvProt_GS[,1] == i,]$gene
    }
    
    venn <- Venn(gsDataList)
    data <- process_data(venn)
    ggplot() +
      # 1. region count layer
      geom_sf(aes(fill = count), data = venn_region(data), show.legend = F) +
      # 2. set edge layer
      geom_sf(aes(color = name), data = venn_setedge(data), show.legend = F, size = 1) +
      # 3. set label layer
      geom_sf_text(aes(label = name), data = venn_setlabel(data), size = 4.5) +
      # 4. region label layer
      geom_sf_label(aes(label = paste0(count)), 
                    data = venn_region(data),
                    size = 4,
                    alpha = 0.5) +
      scale_fill_gradient(low = "#fef4f4", high = "#bf4949")+
      scale_color_manual(values = c("#fef4f4","#fef4f4"))+
      theme_void()
    
    
  })
  
  #render down reg venn diagram
  output$vennDN <- renderPlot({
    
    #Variables
    RNA <- RNAmat()
    prot <- expr_c()
    metaR <- meta_un()
    metaP <- meta_cn()
    metasame <- intersect(metaR[,2],metaP[,2])
    metagroups <- as.vector(levels(factor(metasame)))
    #metagroupsR <- as.vector(levels(factor(metaR[,2])))
    #metagroupsP <- as.vector(levels(factor(metaP[,2])))
    A1 <- metaR[which(metaR[,2] == input$groupAcomp3),1]
    B1 <- metaR[which(metaR[,2] == input$groupBcomp3),1]
    A2 <- metaP[which(metaP[,2] == input$groupAcomp3),1]
    B2 <- metaP[which(metaP[,2] == input$groupBcomp3),1]
    #matrix names
    M1 <- input$MAT1name
    M2 <- input$MAT2name
    
    #RNAseq top table generation
    mat <- RNA[,c(A1,B1)]
    mat <- log2(mat + 1.0)
    groupAOther <- factor(c(rep("A", length(A1)), rep("B", length(B1))))
    designA <- model.matrix(~0 + groupAOther)
    fit <- lmFit(mat, design = designA)
    contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    options(digits = 4)
    top1_R <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
    
    #Protein top table generation
    mat <- prot[,c(A2,B2)]
    mat <- log2(mat + 1.0)
    groupAOther <- factor(c(rep("A", length(A2)), rep("B", length(B2))))
    designA <- model.matrix(~0 + groupAOther)
    fit <- lmFit(mat, design = designA)
    contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    options(digits = 4)
    top1_P <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
    
    #user variables
    FC <- input$vennLog
    adjP <- input$adjpvalVenn
    
    #subset Down reg top tables
    topGenes_RD <- top1_R[top1_R$logFC < FC & top1_R$adj.P.Val <= adjP,]
    topGenes_PD <- top1_P[top1_P$logFC < FC & top1_P$adj.P.Val <= adjP,]
    
    #convert genes into gmt format
    topGenes_RD_gmt <- data.frame(paste("DNreg_in_",gsub(" ","",metagroups[1]),"_",M1,sep = ""),rownames(topGenes_RD))
    colnames(topGenes_RD_gmt) <- c("term","gene")
    topGenes_PD_gmt <- data.frame(paste("DNreg_in_",gsub(" ","",metagroups[1]),"_",M2,sep = ""),rownames(topGenes_PD))
    colnames(topGenes_PD_gmt) <- c("term","gene")
    
    #combine to 4 gene set list
    RNAvProt_GS <- rbind(topGenes_RD_gmt,topGenes_PD_gmt)
    
    #make R Data List
    gsDataList <- list()
    for (i in unique(RNAvProt_GS[,1])){
      gsDataList[[i]] <- RNAvProt_GS[RNAvProt_GS[,1] == i,]$gene
    }
    
    venn <- Venn(gsDataList)
    data <- process_data(venn)
    ggplot() +
      # 1. region count layer
      geom_sf(aes(fill = count), data = venn_region(data), show.legend = F) +
      # 2. set edge layer
      geom_sf(aes(color = name), data = venn_setedge(data), show.legend = F, size = 1) +
      # 3. set label layer
      geom_sf_text(aes(label = name), data = venn_setlabel(data), size = 4.5) +
      # 4. region label layer
      geom_sf_label(aes(label = paste0(count)), 
                    data = venn_region(data),
                    size = 4,
                    alpha = 0.5) +
      scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")+
      scale_color_manual(values = c("#F4FAFE","#F4FAFE"))+
      theme_void()
    
  })
  
  #ssGSEA boxplot1
  output$ssGSEAboxplot1 <- renderPlot({
    
    #matrix names
    M1 <- input$MAT1name
    M2 <- input$MAT2name
    #ssgsea master table
    ssgseaM <- ssGSEAmaster()
    #select GS
    sets <- unique(ssgseaM[,3])
    GS <- sets[1]
    x <- paste(GS," Gene Set",sep = "")
    y <- paste(input$ssGSEAtype," Score",sep = "")
    #plot
    ggplot(data = ssgseaM[which(ssgseaM$variable == GS),], aes(x = Matrix, y = value)) +
      geom_boxplot(aes(fill = type)) +
      theme_bw() +
      labs(x = x, y = y)
    
  })
  
  #ssGSEA boxplot 2
  output$ssGSEAboxplot2 <- renderPlot({
    
    #matrix names
    M1 <- input$MAT1name
    M2 <- input$MAT2name
    #ssgsea master table
    ssgseaM <- ssGSEAmaster()
    #select GS
    sets <- unique(ssgseaM[,3])
    GS <- sets[2]
    x <- paste(GS," Gene Set",sep = "")
    y <- paste(input$ssGSEAtype," Score",sep = "")
    #plot
    ggplot(data = ssgseaM[which(ssgseaM$variable == GS),], aes(x = Matrix, y = value)) +
      geom_boxplot(aes(fill = type)) +
      theme_bw() +
      labs(x = x, y = y)
    
  })
  
  #ssGSEA boxplot 3
  output$ssGSEAboxplot3 <- renderPlot({
    
    #matrix names
    M1 <- input$MAT1name
    M2 <- input$MAT2name
    #ssgsea master table
    ssgseaM <- ssGSEAmaster()
    #select GS
    sets <- unique(ssgseaM[,3])
    GS <- sets[3]
    x <- paste(GS," Gene Set",sep = "")
    y <- paste(input$ssGSEAtype," Score",sep = "")
    #plot
    ggplot(data = ssgseaM[which(ssgseaM$variable == GS),], aes(x = Matrix, y = value)) +
      geom_boxplot(aes(fill = type)) +
      theme_bw() +
      labs(x = x, y = y)
    
  })
  
  #ssGSEA boxplot 4
  output$ssGSEAboxplot4 <- renderPlot({
    
    #matrix names
    M1 <- input$MAT1name
    M2 <- input$MAT2name
    #ssgsea master table
    ssgseaM <- ssGSEAmaster()
    #select GS
    sets <- unique(ssgseaM[,3])
    GS <- sets[4]
    x <- paste(GS," Gene Set",sep = "")
    y <- paste(input$ssGSEAtype," Score",sep = "")
    #plot
    ggplot(data = ssgseaM[which(ssgseaM$variable == GS),], aes(x = Matrix, y = value)) +
      geom_boxplot(aes(fill = type)) +
      theme_bw() +
      labs(x = x, y = y)
    
  })
  
  
  ####----Text----####
  
  
  output$NESandPval1 <- renderUI({
    
    res <- recGSEAg1()
    gmt <- recGMT()
    M1 <- input$MAT1name
    M2 <- input$MAT2name
    gsea.df <- as.data.frame(res@result)
    #use GS derived from matrix 2
    M2 <- input$MAT2name
    sets <- unique(gmt[,1])
    GS <- sets[grep(paste('_',M2,sep = ''),unique(gmt[,1]))]
    ##make text output
    gs1 <- GS[1]
    gs2 <- GS[2]
    NES1 = gsea.df$NES[which(gsea.df[,1]==gs1)]
    Pval1 = gsea.df$pvalue[which(gsea.df[,1]==gs1)]
    NES.o1 <- paste0("<u><b>NES:</b></u> ", NES1)
    Pval.o1 <- paste0("<u><b>Pvalue:</b></u> ", Pval1)
    message1 <- paste("Gene sets derived from DEG analysis of <b>",M2, "</b> were used to perform GSEA on <b>",M1,"</b>.", sep = "")
    if (NES1 > 0){
      UpOrDown1 <- paste("Based on the normalized enrichment score above, the<b>", gs1, "</b>gene set is upregulated in<b>", input$groupAcomp2 , "</b>group.")
    }
    if (NES1 < 0) {
      UpOrDown1 <- paste("Based on the normalized enrichment score above, the<b>", gs1, "</b>gene set is downregulated in<b>", input$groupBcomp2 , "</b>group.")
    }
    NES2 = gsea.df$NES[which(gsea.df[,1]==gs2)]
    Pval2 = gsea.df$pvalue[which(gsea.df[,1]==gs2)]
    NES.o2 <- paste0("<u><b>NES:</b></u> ", NES2)
    Pval.o2 <- paste0("<u><b>Pvalue:</b></u> ", Pval2)
    if (NES2 > 0){
      UpOrDown2 <- paste("Based on the normalized enrichment score above, the<b>", gs2, "</b>gene set is upregulated in<b>", input$groupAcomp2 , "</b>group.")
    }
    if (NES2 < 0) {
      UpOrDown2 <- paste("Based on the normalized enrichment score above, the<b>", gs2, "</b>gene set is downregulated in<b>", input$groupBcomp2 , "</b>group.")
    }
    HTML(paste(message1,'<br/>',"<u><b>",gs1," Gene Set","</u></b>",'<br/>',NES.o1,'<br/>',Pval.o1,'<br/>',UpOrDown1,'<br/>',"<u><b>",gs2," Gene Set","</u></b>",'<br/>',NES.o2,'<br/>',Pval.o2,'<br/>',UpOrDown2, sep = ''))
    
  })
  
  output$NESandPval2 <- renderUI({
    
    res <- recGSEAg2()
    gmt <- recGMT()
    M1 <- input$MAT1name
    M2 <- input$MAT2name
    gsea.df <- as.data.frame(res@result)
    #use GS derived from matrix 2
    M1 <- input$MAT1name
    sets <- unique(gmt[,1])
    GS <- sets[grep(paste('_',M1,sep = ''),unique(gmt[,1]))]
    ##make text output
    gs1 <- GS[1]
    gs2 <- GS[2]
    NES1 = gsea.df$NES[which(gsea.df[,1]==gs1)]
    Pval1 = gsea.df$pvalue[which(gsea.df[,1]==gs1)]
    NES.o1 <- paste0("<u><b>NES:</b></u> ", NES1)
    Pval.o1 <- paste0("<u><b>Pvalue:</b></u> ", Pval1)
    message1 <- paste("Gene sets derived from DEG analysis of <b>",M1, "</b> were used to perform GSEA on <b>",M2,"</b>.", sep = "")
    if (NES1 > 0){
      UpOrDown1 <- paste("Based on the normalized enrichment score above, the<b>", gs1, "</b>gene set is upregulated in<b>", input$groupAcomp2 , "</b>group.")
    }
    if (NES1 < 0) {
      UpOrDown1 <- paste("Based on the normalized enrichment score above, the<b>", gs1, "</b>gene set is downregulated in<b>", input$groupBcomp2 , "</b>group.")
    }
    NES2 = gsea.df$NES[which(gsea.df[,1]==gs2)]
    Pval2 = gsea.df$pvalue[which(gsea.df[,1]==gs2)]
    NES.o2 <- paste0("<u><b>NES:</b></u> ", NES2)
    Pval.o2 <- paste0("<u><b>Pvalue:</b></u> ", Pval2)
    if (NES2 > 0){
      UpOrDown2 <- paste("Based on the normalized enrichment score above, the<b>", gs2, "</b>gene set is upregulated in<b>", input$groupAcomp2 , "</b>group.")
    }
    if (NES2 < 0) {
      UpOrDown2 <- paste("Based on the normalized enrichment score above, the<b>", gs2, "</b>gene set is downregulated in<b>", input$groupBcomp2 , "</b>group.")
    }
    HTML(paste(message1,'<br/>',"<u><b>",gs1," Gene Set","</u></b>",'<br/>',NES.o1,'<br/>',Pval.o1,'<br/>',UpOrDown1,'<br/>',"<u><b>",gs2," Gene Set","</u></b>",'<br/>',NES.o2,'<br/>',Pval.o2,'<br/>',UpOrDown2, sep = ''))
    
  })
  
  #written explanation for venn diagrams
  output$VennExplan <- renderUI({
    
    M1 <- input$MAT1name
    M2 <- input$MAT2name
    HTML(paste("The VennDiagrams below are comparing up-regulated and down-regulated differentially expressed genes between matrix 1: <u><b>",
               M1,"</u></b> and Matrix 2: <u><b>",M2,"</u></b>.",sep = ""))
    
  })
  
  #Written Explanation for Stat table
  output$StatExplan <- renderUI({
    
    #require user GMT input
    #req(input$GMTcomp)
    M1 <- input$MAT1name
    M2 <- input$MAT2name
    HTML(paste("The table generated below displays statistical analyses from comparing differentially expressed genes (DEGs) from the given matrices with gene sets from the user provided GMT file.<br/>
                   The first column represents the number of similar genes between all DEGs in both given matrices and the gene set being compared in that row.<br/>
                   The subsequent columns list the P-Value, Odds Ratio, Kappa Score, and Jaccard Index, which are derived from comparing the gene set of each row to the Up and Down regulated genes from each matrices.<br/>
                   Please note that the blank values denote Infinity."))
    
  })
  
  #render download button for stat table
  output$msigstatdnld <- downloadHandler(
    filename = function() {
      name <- input$msigdnldname
      paste(name,".tsv",sep = "")
    },
    content = function(file) {
      #Assign matrices to variables and assign meta groups
      M1 <- RNAmat()
      M2 <- expr_c()
      meta1 <- meta_un()
      meta2 <- meta_cn()
      #metasame <- merge(meta1,meta2)
      A1 <- meta1[which(meta1[,2] == input$groupAcomp3),1]
      B1 <- meta1[which(meta1[,2] == input$groupBcomp3),1]
      A2 <- meta2[which(meta2[,2] == input$groupAcomp3),1]
      B2 <- meta2[which(meta2[,2] == input$groupBcomp3),1]
      
      #Assign user given matrix names for labeling
      M1Name <- input$MAT1name
      M2Name <- input$MAT2name
      
      #Matrix 1 top table generation
      mat <- M1[,c(A1,B1)]
      mat <- log2(mat + 1.0)
      groupAOther <- factor(c(rep("A", length(A1)), rep("B", length(B1))))
      designA <- model.matrix(~0 + groupAOther)
      fit <- lmFit(mat, design = designA)
      contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
      fit2 <- contrasts.fit(fit, contrast.matrix)
      fit2 <- eBayes(fit2)
      options(digits = 4)
      top1_1 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
      
      #Matrix 2 top table generation
      mat <- M2[,c(A2,B2)]
      mat <- log2(mat + 1.0)
      groupAOther <- factor(c(rep("A", length(A2)), rep("B", length(B2))))
      designA <- model.matrix(~0 + groupAOther)
      fit <- lmFit(mat, design = designA)
      contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
      fit2 <- contrasts.fit(fit, contrast.matrix)
      fit2 <- eBayes(fit2)
      options(digits = 4)
      top1_2 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
      
      #re-assign main files and variables
      M1top <- top1_1
      M2top <- top1_2
      
      #user selected GS category
      GScat <- MSigDBcat3[input$msigdbtable_rows_selected,2]
      #select GS names based on category selected
      GScat_sets <- MSigDBcat[which(MSigDBcat[,2] == GScat),3]
      #add GS to GSlist based off of what is selected
      gs_list <- gs[GScat_sets]
      gmt_unique <- unique(unlist(gs_list, use.names = F))
      
      #list of similar genes between 3 datasets
      SimGenes <- intersect(intersect(rownames(M1top),rownames(M2top)),gmt_unique)
      #subset data based on similar genes only
      M1top2 <- M1top[SimGenes,] 
      M2top2 <- M2top[SimGenes,]
      #extract genes from gslist that are not in matrices
      for (i in 1:length(gs_list)){
        gs_list[[i]] <- gs_list[[i]][which(gs_list[[i]] %in% SimGenes)] 
      }
      #list of unique but similar to matrices genes in GMT
      gmt2Uniq <- unique(unlist(gs_list, use.names = F))
      
      #user variables
      FC <- input$vennLog
      adjP <- input$adjpvalVenn
      
      ##get gene lists from top tables
      #subset up reg DEGs and non DEGs - M1
      DEG_M1U <- rownames(M1top2[which(M1top2$logFC > FC & M1top2$adj.P.Val <= adjP),])
      nonDEG_M1U <- rownames(M1top2[which(!(M1top2$logFC > FC & M1top2$adj.P.Val <= adjP)),])
      #subset up reg DEGs and non DEGs - M2
      DEG_M2U <- rownames(M2top2[which(M2top2$logFC > FC & M2top2$adj.P.Val <= adjP),])
      nonDEG_M2U <- rownames(M2top2[which(!(M2top2$logFC > FC & M2top2$adj.P.Val <= adjP)),])
      #subset down reg DEGs and non DEGs - M1
      DEG_M1D <- rownames(M1top2[which(M1top2$logFC < FC & M1top2$adj.P.Val <= adjP),])
      nonDEG_M1D <- rownames(M1top2[which(!(M1top2$logFC < FC & M1top2$adj.P.Val <= adjP)),])
      #subset down reg DEGs and non DEGs - M2
      DEG_M2D <- rownames(M2top2[which(M2top2$logFC < FC & M2top2$adj.P.Val <= adjP),])
      nonDEG_M2D <- rownames(M2top2[which(!(M2top2$logFC < FC & M2top2$adj.P.Val <= adjP)),])
      
      #Generate names for the DEG lists made above - will be used to compare with GMT
      gsName1 <- paste("UpRegDEG_in_",M1Name,sep = "")
      gsName2 <- paste("DnRegDEG_in_",M1Name,sep = "")
      gsName3 <- paste("UpRegDEG_in_",M2Name,sep = "")
      gsName4 <- paste("DnRegDEG_in_",M2Name,sep = "")
      
      #Add DEGs subset above to GMT list
      gsDataListGMT <- list() #start gene set list
      gsDataListGMT[[gsName1]] <- DEG_M1U
      gsDataListGMT[[gsName3]] <- DEG_M2U
      gsDataListGMT[[gsName2]] <- DEG_M1D
      gsDataListGMT[[gsName4]] <- DEG_M2D
      
      gsDataListGMT <- c(gsDataListGMT,gs_list)
      
      #world of all genes from matrices
      DEG_M1M2 <- intersect(rownames(M1top2), rownames(M2top2))
      
      ##make gene set list from GMT
      #for (i in unique(gmt2[,1])){
      #    gsDataListGMT[[i]] <- gmt2[gmt2[,1] == i,]$gene
      #}
      
      #learn Jaccard index function
      jaccard <- function(a, b) {
        intersection = length(intersect(a, b))
        union = length(a) + length(b) - intersection
        return (intersection/union)
      }
      
      #list to store stat data
      StatDataList <- list()
      
      #loop through gmt file and get stats when comparing gene sets to expression matrices
      for (i in 1:length(gsDataListGMT)) {
        
        #Gene set name that is being compared
        GS <- names(gsDataListGMT)[[i]]
        
        #Integration of UP regulated genes with gmt gene set - M1
        M1DEG_UP_GS <- length(intersect(DEG_M1U,gsDataListGMT[[i]]))                        #DE up + in GS
        M1nDEG_UP_nGS <- length(intersect(nonDEG_M1U,setdiff(gmt2Uniq,gsDataListGMT[[i]]))) #non DE up + not in GS 
        M1DEG_UP_nGS <- length(intersect(DEG_M1U,setdiff(gmt2Uniq,gsDataListGMT[[i]])))     #DE up + not in GS
        M1nDEG_UP_GS <- length(intersect(nonDEG_M1U,gsDataListGMT[[i]]))                    #non DE up + in GS
        #Integration of DOWN regulated genes with gmt gene set - M1
        M1DEG_DN_GS <- length(intersect(DEG_M1D,gsDataListGMT[[i]]))                        #DE dn + in GS
        M1nDEG_DN_nGS <- length(intersect(nonDEG_M1D,setdiff(gmt2Uniq,gsDataListGMT[[i]]))) #non DE dn + not in GS
        M1DEG_DN_nGS <- length(intersect(DEG_M1D,setdiff(gmt2Uniq,gsDataListGMT[[i]])))     #DE dn + not in GS
        M1nDEG_DN_GS <- length(intersect(nonDEG_M1D,gsDataListGMT[[i]]))                    #non DE dn + in GS
        
        #Integration of UP regulated genes with gmt gene set - M2
        M2DEG_UP_GS <- length(intersect(DEG_M2U,gsDataListGMT[[i]]))
        M2nDEG_UP_nGS <- length(intersect(nonDEG_M2U,setdiff(gmt2Uniq,gsDataListGMT[[i]])))
        M2DEG_UP_nGS <- length(intersect(DEG_M2U,setdiff(gmt2Uniq,gsDataListGMT[[i]])))
        M2nDEG_UP_GS <- length(intersect(nonDEG_M2U,gsDataListGMT[[i]]))
        #Integration of DOWN regulated genes with gmt gene set - M2
        M2DEG_DN_GS <- length(intersect(DEG_M2D,gsDataListGMT[[i]]))
        M2nDEG_DN_nGS <- length(intersect(nonDEG_M2D,setdiff(gmt2Uniq,gsDataListGMT[[i]])))
        M2DEG_DN_nGS <- length(intersect(DEG_M2D,setdiff(gmt2Uniq,gsDataListGMT[[i]])))
        M2nDEG_DN_GS <- length(intersect(nonDEG_M2D,gsDataListGMT[[i]]))
        
        #number of genes in 'world' AND GS
        WorldNGS <- length(intersect(DEG_M1M2,gsDataListGMT[[i]]))
        
        ##Generate Matrices
        #up regulated - M1
        UpMat1GS <- matrix(data = c(M1DEG_UP_GS,M1DEG_UP_nGS,
                                    M1nDEG_UP_GS,M1nDEG_UP_nGS), nrow = 2)
        #write out numbers of matrix for last columns
        UpMat1GS.n <- paste(M1DEG_UP_GS,M1DEG_UP_nGS,
                            M1nDEG_UP_GS,M1nDEG_UP_nGS, sep = ",")
        
        #down regulated - M1
        DnMat1GS <- matrix(data = c(M1DEG_DN_GS,M1DEG_DN_nGS,
                                    M1nDEG_DN_GS,M1nDEG_DN_nGS), nrow = 2)
        #write out numbers of matrix for last columns
        DnMat1GS.n <- paste(M1DEG_DN_GS,M1DEG_DN_nGS,
                            M1nDEG_DN_GS,M1nDEG_DN_nGS, sep = ",")
        
        #up regulated - M2
        UpMat2GS <- matrix(data = c(M2DEG_UP_GS,M2DEG_UP_nGS,
                                    M2nDEG_UP_GS,M2nDEG_UP_nGS), nrow = 2)
        #write out numbers of matrix for last columns
        UpMat2GS.n <- paste(M2DEG_UP_GS,M2DEG_UP_nGS,
                            M2nDEG_UP_GS,M2nDEG_UP_nGS, sep = ",")
        
        #down regulated - M2
        DnMat2GS <- matrix(data = c(M2DEG_DN_GS,M2DEG_DN_nGS,
                                    M2nDEG_DN_GS,M2nDEG_DN_nGS), nrow = 2)
        #write out numbers of matrix for last columns
        DnMat2GS.n <- paste(M2DEG_DN_GS,M2DEG_DN_nGS,
                            M2nDEG_DN_GS,M2nDEG_DN_nGS, sep = ",")
        
        ##Fisher Exact Test
        #Up regulated - M1
        UpfishM1 <- fisher.test(UpMat1GS, alternative = "two.sided")
        UpFishPM1 <- UpfishM1[[1]]
        UpFishORM1 <- UpfishM1[[3]]
        #Down regulated - M1
        DnfishM1 <- fisher.test(DnMat1GS, alternative = "two.sided")
        DnFishPM1 <- DnfishM1[[1]]
        DnFishORM1 <- DnfishM1[[3]]
        #Up regulated - M2
        UpfishM2 <- fisher.test(UpMat2GS, alternative = "two.sided")
        UpFishPM2 <- UpfishM2[[1]]
        UpFishORM2 <- UpfishM2[[3]]
        #Down regulated - M2
        DnfishM2 <- fisher.test(DnMat2GS, alternative = "two.sided")
        DnFishPM2 <- DnfishM2[[1]]
        DnFishORM2 <- DnfishM2[[3]]
        
        ##Kappa Score
        #Mat1
        UpKappaM1 <- kappa(UpMat1GS)
        DnKappaM1 <- kappa(DnMat1GS)
        #Mat1
        UpKappaM2 <- kappa(UpMat2GS)
        DnKappaM2 <- kappa(DnMat2GS)
        
        #up reg genes in M1
        UpJaccM1 <- jaccard(DEG_M1U,gsDataListGMT[[i]])
        DnJaccM1 <- jaccard(DEG_M1D,gsDataListGMT[[i]])
        #up reg genes in M2
        UpJaccM2 <- jaccard(DEG_M2U,gsDataListGMT[[i]])
        DnJaccM2 <- jaccard(DEG_M2D,gsDataListGMT[[i]])
        
        #Assign stats to GS name and add to list
        StatDataList[[GS]] <- c(WorldNGS,
                                UpFishPM1,UpFishPM2,DnFishPM1,DnFishPM2,
                                UpFishORM1,UpFishORM2,DnFishORM1,DnFishORM2,
                                UpKappaM1,UpKappaM2,DnKappaM1,DnKappaM2,
                                UpJaccM1,UpJaccM2,DnJaccM1,DnJaccM2,
                                UpMat1GS.n,UpMat2GS.n,DnMat1GS.n,DnMat2GS.n)
        
      }
      
      #convert to data frame
      StatData.df <- do.call(rbind.data.frame, StatDataList)
      rownames(StatData.df) <- names(StatDataList)
      
      colnames(StatData.df) <- c("Number_of_DEGenes_in_GeneSet",
                                 paste(M1Name,"_UpRegulated_Pvalue",sep = ""),paste(M2Name,"_UpRegulated_Pvalue",sep = ""),
                                 paste(M1Name,"_DownRegulated_Pvalue",sep = ""),paste(M2Name,"_DownRegulated_Pvalue",sep = ""),
                                 paste(M1Name,"_UpRegulated_OddsRatio",sep = ""),paste(M2Name,"_UpRegulated_OddsRatio",sep = ""),
                                 paste(M1Name,"_DownRegulated_OddsRatio",sep = ""),paste(M2Name,"_DownRegulated_OddsRatio",sep = ""),
                                 paste(M1Name,"_UpRegulated_KappaScore",sep = ""),paste(M2Name,"_UpRegulated_KappaScore",sep = ""),
                                 paste(M1Name,"_DownRegulated_KappaScore",sep = ""),paste(M2Name,"_DownRegulated_KappaScore",sep = ""),
                                 paste(M1Name,"_UpRegulated_JaccardIndex",sep = ""),paste(M2Name,"_UpRegulated_JaccardIndex",sep = ""),
                                 paste(M1Name,"_DownRegulated_JaccardIndex",sep = ""),paste(M2Name,"_DownRegulated_JaccardIndex",sep = ""),
                                 paste(M1Name,"_UpRegulated_Matrix",sep = ""),paste(M2Name,"_UpRegulated_Matrix",sep = ""),
                                 paste(M1Name,"_DownRegulated_Matrix",sep = ""),paste(M2Name,"_DownRegulated_Matrix",sep = ""))
      
      #make labeled column for rownames
      StatData.df$GeneSets <- rownames(StatData.df)
      StatData.df <- StatData.df %>%
        relocate(GeneSets)
      #convert numeric column to as.numeric
      StatData.df[,c(2:18)] <- sapply(StatData.df[,c(2:18)], as.numeric)
      write_tsv(StatData.df,file)
      
    }
  )
  
  #render download button for stat table
  output$StatTabDown <- downloadHandler(
    filename = function() {
      name <- input$StatTabName
      paste(name,".tsv",sep = "")
    },
    content = function(file) {
      #user upload of gmt file
      gmt.u <- input$GMTcomp
      ext <- tools::file_ext(gmt.u$datapath)
      req(gmt.u)
      validate(need(ext == c("gmt"), "Please upload a .gmt file"))
      gmt <- read.gmt(gmt.u$datapath)
      
      #Assign matrices to variables and assign meta groups
      M1 <- RNAmat()
      M2 <- expr_c()
      meta1 <- meta_un()
      meta2 <- meta_cn()
      #metasame <- merge(meta1,meta2)
      A1 <- meta1[which(meta1[,2] == input$groupAcomp3),1]
      B1 <- meta1[which(meta1[,2] == input$groupBcomp3),1]
      A2 <- meta2[which(meta2[,2] == input$groupAcomp3),1]
      B2 <- meta2[which(meta2[,2] == input$groupBcomp3),1]
      
      #Assign user given matrix names for labeling
      M1Name <- input$MAT1name
      M2Name <- input$MAT2name
      
      #Matrix 1 top table generation
      mat <- M1[,c(A1,B1)]
      mat <- log2(mat + 1.0)
      groupAOther <- factor(c(rep("A", length(A1)), rep("B", length(B1))))
      designA <- model.matrix(~0 + groupAOther)
      fit <- lmFit(mat, design = designA)
      contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
      fit2 <- contrasts.fit(fit, contrast.matrix)
      fit2 <- eBayes(fit2)
      options(digits = 4)
      top1_1 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
      
      #Matrix 2 top table generation
      mat <- M2[,c(A2,B2)]
      mat <- log2(mat + 1.0)
      groupAOther <- factor(c(rep("A", length(A2)), rep("B", length(B2))))
      designA <- model.matrix(~0 + groupAOther)
      fit <- lmFit(mat, design = designA)
      contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
      fit2 <- contrasts.fit(fit, contrast.matrix)
      fit2 <- eBayes(fit2)
      options(digits = 4)
      top1_2 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
      
      #re-assign main files and variables
      M1top <- top1_1
      M2top <- top1_2
      
      #list of similar genes between 3 datasets
      SimGenes <- intersect(intersect(rownames(M1top),rownames(M2top)),gmt[,2])
      #subset data based on similar genes only
      M1top2 <- M1top[SimGenes,] 
      M2top2 <- M2top[SimGenes,]
      gmt2 <- gmt[which(gmt[,2] %in% SimGenes),]
      #list of unique bu similar to matrices genes in GMT
      gmt2Uniq <- unique(gmt2[,2])
      
      #user variables
      FC <- input$vennLog
      adjP <- input$adjpvalVenn
      
      ##get gene lists from top tables
      #subset up reg DEGs and non DEGs - M1
      DEG_M1U <- rownames(M1top2[which(M1top2$logFC > FC & M1top2$adj.P.Val <= adjP),])
      nonDEG_M1U <- rownames(M1top2[which(!(M1top2$logFC > FC & M1top2$adj.P.Val <= adjP)),])
      #subset up reg DEGs and non DEGs - M2
      DEG_M2U <- rownames(M2top2[which(M2top2$logFC > FC & M2top2$adj.P.Val <= adjP),])
      nonDEG_M2U <- rownames(M2top2[which(!(M2top2$logFC > FC & M2top2$adj.P.Val <= adjP)),])
      #subset down reg DEGs and non DEGs - M1
      DEG_M1D <- rownames(M1top2[which(M1top2$logFC < FC & M1top2$adj.P.Val <= adjP),])
      nonDEG_M1D <- rownames(M1top2[which(!(M1top2$logFC < FC & M1top2$adj.P.Val <= adjP)),])
      #subset down reg DEGs and non DEGs - M2
      DEG_M2D <- rownames(M2top2[which(M2top2$logFC < FC & M2top2$adj.P.Val <= adjP),])
      nonDEG_M2D <- rownames(M2top2[which(!(M2top2$logFC < FC & M2top2$adj.P.Val <= adjP)),])
      
      #Generate names for the DEG lists made above - will be used to compare with GMT
      gsName1 <- paste("UpRegDEG_in_",M1Name,sep = "")
      gsName2 <- paste("DnRegDEG_in_",M1Name,sep = "")
      gsName3 <- paste("UpRegDEG_in_",M2Name,sep = "")
      gsName4 <- paste("DnRegDEG_in_",M2Name,sep = "")
      
      #Add DEGs subset above to GMT list
      gsDataListGMT <- list() #start gene set list
      gsDataListGMT[[gsName1]] <- DEG_M1U
      gsDataListGMT[[gsName3]] <- DEG_M2U
      gsDataListGMT[[gsName2]] <- DEG_M1D
      gsDataListGMT[[gsName4]] <- DEG_M2D
      
      #world of all genes from matrices
      DEG_M1M2 <- intersect(rownames(M1top2), rownames(M2top2))
      
      #make gene set list from GMT
      for (i in unique(gmt2[,1])){
        gsDataListGMT[[i]] <- gmt2[gmt2[,1] == i,]$gene
      }
      
      #learn Jaccard index function
      jaccard <- function(a, b) {
        intersection = length(intersect(a, b))
        union = length(a) + length(b) - intersection
        return (intersection/union)
      }
      
      #list to store stat data
      StatDataList <- list()
      
      #loop through gmt file and get stats when comparing gene sets to expression matrices
      for (i in 1:length(gsDataListGMT)) {
        
        #Gene set name that is being compared
        GS <- names(gsDataListGMT)[[i]]
        
        #Integration of UP regulated genes with gmt gene set - M1
        M1DEG_UP_GS <- length(intersect(DEG_M1U,gsDataListGMT[[i]]))                        #DE up + in GS
        M1nDEG_UP_nGS <- length(intersect(nonDEG_M1U,setdiff(gmt2Uniq,gsDataListGMT[[i]]))) #non DE up + not in GS 
        M1DEG_UP_nGS <- length(intersect(DEG_M1U,setdiff(gmt2Uniq,gsDataListGMT[[i]])))     #DE up + not in GS
        M1nDEG_UP_GS <- length(intersect(nonDEG_M1U,gsDataListGMT[[i]]))                    #non DE up + in GS
        #Integration of DOWN regulated genes with gmt gene set - M1
        M1DEG_DN_GS <- length(intersect(DEG_M1D,gsDataListGMT[[i]]))                        #DE dn + in GS
        M1nDEG_DN_nGS <- length(intersect(nonDEG_M1D,setdiff(gmt2Uniq,gsDataListGMT[[i]]))) #non DE dn + not in GS
        M1DEG_DN_nGS <- length(intersect(DEG_M1D,setdiff(gmt2Uniq,gsDataListGMT[[i]])))     #DE dn + not in GS
        M1nDEG_DN_GS <- length(intersect(nonDEG_M1D,gsDataListGMT[[i]]))                    #non DE dn + in GS
        
        #Integration of UP regulated genes with gmt gene set - M2
        M2DEG_UP_GS <- length(intersect(DEG_M2U,gsDataListGMT[[i]]))
        M2nDEG_UP_nGS <- length(intersect(nonDEG_M2U,setdiff(gmt2Uniq,gsDataListGMT[[i]])))
        M2DEG_UP_nGS <- length(intersect(DEG_M2U,setdiff(gmt2Uniq,gsDataListGMT[[i]])))
        M2nDEG_UP_GS <- length(intersect(nonDEG_M2U,gsDataListGMT[[i]]))
        #Integration of DOWN regulated genes with gmt gene set - M2
        M2DEG_DN_GS <- length(intersect(DEG_M2D,gsDataListGMT[[i]]))
        M2nDEG_DN_nGS <- length(intersect(nonDEG_M2D,setdiff(gmt2Uniq,gsDataListGMT[[i]])))
        M2DEG_DN_nGS <- length(intersect(DEG_M2D,setdiff(gmt2Uniq,gsDataListGMT[[i]])))
        M2nDEG_DN_GS <- length(intersect(nonDEG_M2D,gsDataListGMT[[i]]))
        
        #number of genes in 'world' AND GS
        WorldNGS <- length(intersect(DEG_M1M2,gsDataListGMT[[i]]))
        
        ##Generate Matrices
        #up regulated - M1
        UpMat1GS <- matrix(data = c(M1DEG_UP_GS,M1DEG_UP_nGS,
                                    M1nDEG_UP_GS,M1nDEG_UP_nGS), nrow = 2)
        #write out numbers of matrix for last columns
        UpMat1GS.n <- paste(M1DEG_UP_GS,M1DEG_UP_nGS,
                            M1nDEG_UP_GS,M1nDEG_UP_nGS, sep = ",")
        
        #down regulated - M1
        DnMat1GS <- matrix(data = c(M1DEG_DN_GS,M1DEG_DN_nGS,
                                    M1nDEG_DN_GS,M1nDEG_DN_nGS), nrow = 2)
        #write out numbers of matrix for last columns
        DnMat1GS.n <- paste(M1DEG_DN_GS,M1DEG_DN_nGS,
                            M1nDEG_DN_GS,M1nDEG_DN_nGS, sep = ",")
        
        #up regulated - M2
        UpMat2GS <- matrix(data = c(M2DEG_UP_GS,M2DEG_UP_nGS,
                                    M2nDEG_UP_GS,M2nDEG_UP_nGS), nrow = 2)
        #write out numbers of matrix for last columns
        UpMat2GS.n <- paste(M2DEG_UP_GS,M2DEG_UP_nGS,
                            M2nDEG_UP_GS,M2nDEG_UP_nGS, sep = ",")
        
        #down regulated - M2
        DnMat2GS <- matrix(data = c(M2DEG_DN_GS,M2DEG_DN_nGS,
                                    M2nDEG_DN_GS,M2nDEG_DN_nGS), nrow = 2)
        #write out numbers of matrix for last columns
        DnMat2GS.n <- paste(M2DEG_DN_GS,M2DEG_DN_nGS,
                            M2nDEG_DN_GS,M2nDEG_DN_nGS, sep = ",")
        
        ##Fisher Exact Test
        #Up regulated - M1
        UpfishM1 <- fisher.test(UpMat1GS, alternative = "two.sided")
        UpFishPM1 <- UpfishM1[[1]]
        UpFishORM1 <- UpfishM1[[3]]
        #Down regulated - M1
        DnfishM1 <- fisher.test(DnMat1GS, alternative = "two.sided")
        DnFishPM1 <- DnfishM1[[1]]
        DnFishORM1 <- DnfishM1[[3]]
        #Up regulated - M2
        UpfishM2 <- fisher.test(UpMat2GS, alternative = "two.sided")
        UpFishPM2 <- UpfishM2[[1]]
        UpFishORM2 <- UpfishM2[[3]]
        #Down regulated - M2
        DnfishM2 <- fisher.test(DnMat2GS, alternative = "two.sided")
        DnFishPM2 <- DnfishM2[[1]]
        DnFishORM2 <- DnfishM2[[3]]
        
        ##Kappa Score
        #Mat1
        UpKappaM1 <- kappa(UpMat1GS)
        DnKappaM1 <- kappa(DnMat1GS)
        #Mat1
        UpKappaM2 <- kappa(UpMat2GS)
        DnKappaM2 <- kappa(DnMat2GS)
        
        #up reg genes in M1
        UpJaccM1 <- jaccard(DEG_M1U,gsDataListGMT[[i]])
        DnJaccM1 <- jaccard(DEG_M1D,gsDataListGMT[[i]])
        #up reg genes in M2
        UpJaccM2 <- jaccard(DEG_M2U,gsDataListGMT[[i]])
        DnJaccM2 <- jaccard(DEG_M2D,gsDataListGMT[[i]])
        
        #Assign stats to GS name and add to list
        StatDataList[[GS]] <- c(WorldNGS,
                                UpFishPM1,UpFishPM2,DnFishPM1,DnFishPM2,
                                UpFishORM1,UpFishORM2,DnFishORM1,DnFishORM2,
                                UpKappaM1,UpKappaM2,DnKappaM1,DnKappaM2,
                                UpJaccM1,UpJaccM2,DnJaccM1,DnJaccM2,
                                UpMat1GS.n,UpMat2GS.n,DnMat1GS.n,DnMat2GS.n)
        
      }
      
      #convert to data frame
      StatData.df <- do.call(rbind.data.frame, StatDataList)
      rownames(StatData.df) <- names(StatDataList)
      #assign column names
      colnames(StatData.df) <- c("Number_of_DEGenes_in_GeneSet",
                                 paste(M1Name,"_UpRegulated_Pvalue",sep = ""),paste(M2Name,"_UpRegulated_Pvalue",sep = ""),
                                 paste(M1Name,"_DownRegulated_Pvalue",sep = ""),paste(M2Name,"_DownRegulated_Pvalue",sep = ""),
                                 paste(M1Name,"_UpRegulated_OddsRatio",sep = ""),paste(M2Name,"_UpRegulated_OddsRatio",sep = ""),
                                 paste(M1Name,"_DownRegulated_OddsRatio",sep = ""),paste(M2Name,"_DownRegulated_OddsRatio",sep = ""),
                                 paste(M1Name,"_UpRegulated_KappaScore",sep = ""),paste(M2Name,"_UpRegulated_KappaScore",sep = ""),
                                 paste(M1Name,"_DownRegulated_KappaScore",sep = ""),paste(M2Name,"_DownRegulated_KappaScore",sep = ""),
                                 paste(M1Name,"_UpRegulated_JaccardIndex",sep = ""),paste(M2Name,"_UpRegulated_JaccardIndex",sep = ""),
                                 paste(M1Name,"_DownRegulated_JaccardIndex",sep = ""),paste(M2Name,"_DownRegulated_JaccardIndex",sep = ""),
                                 paste(M1Name,"_UpRegulated_Matrix",sep = ""),paste(M2Name,"_UpRegulated_Matrix",sep = ""),
                                 paste(M1Name,"_DownRegulated_Matrix",sep = ""),paste(M2Name,"_DownRegulated_Matrix",sep = ""))
      #make labeled column for rownames
      StatData.df$GeneSets <- rownames(StatData.df)
      StatData.df <- StatData.df %>%
        relocate(GeneSets)
      #convert numeric column to as.numeric
      StatData.df[,c(2:18)] <- sapply(StatData.df[,c(2:18)], as.numeric)
      write_tsv(StatData.df,file)
      
    }
  )
  
  #render download button for reciprocal ssGSEA table - 1
  output$ssGSEADown1 <- downloadHandler(
    filename = function() {
      M1 <- input$MAT1name
      M2 <- input$MAT2name
      ssType <- input$ssGSEAtype
      #master table
      ssgseaM <- ssGSEAmaster()
      #select GS
      sets <- as.vector(unique(ssgseaM[,3]))
      GS <- sets[1]
      paste(GS,ssType,"Score.tsv",sep = "_")
    },
    content = function(file) {
      #master table
      ssgseaM <- ssGSEAmaster()
      #select GS
      sets <- as.vector(unique(ssgseaM[,3]))
      GS <- sets[1]
      #reformat
      sscores <- ssgseaM[which(ssgseaM$variable == GS),]
      sscores <- sscores %>%
        select(Sample,type,Matrix,value)
      colnames(sscores)[2] <- "Type"
      colnames(sscores)[4] <- GS
      write_tsv(sscores, file)
    }
  )
  
  #render download button for reciprocal ssGSEA table - 2
  output$ssGSEADown2 <- downloadHandler(
    filename = function() {
      M1 <- input$MAT1name
      M2 <- input$MAT2name
      ssType <- input$ssGSEAtype
      #master table
      ssgseaM <- ssGSEAmaster()
      #select GS
      sets <- as.vector(unique(ssgseaM[,3]))
      GS <- sets[2]
      paste(GS,ssType,"Score.tsv",sep = "_")
    },
    content = function(file) {
      #master table
      ssgseaM <- ssGSEAmaster()
      #select GS
      sets <- as.vector(unique(ssgseaM[,3]))
      GS <- sets[2]
      #reformat
      sscores <- ssgseaM[which(ssgseaM$variable == GS),]
      sscores <- sscores %>%
        select(Sample,type,Matrix,value)
      colnames(sscores)[2] <- "Type"
      colnames(sscores)[4] <- GS
      write_tsv(sscores, file)
    }
  )
  
  #render download button for reciprocal ssGSEA table - 3
  output$ssGSEADown3 <- downloadHandler(
    filename = function() {
      M1 <- input$MAT1name
      M2 <- input$MAT2name
      ssType <- input$ssGSEAtype
      #master table
      ssgseaM <- ssGSEAmaster()
      #select GS
      sets <- as.vector(unique(ssgseaM[,3]))
      GS <- sets[3]
      paste(GS,ssType,"Score.tsv",sep = "_")
    },
    content = function(file) {
      #master table
      ssgseaM <- ssGSEAmaster()
      #select GS
      sets <- as.vector(unique(ssgseaM[,3]))
      GS <- sets[3]
      #reformat
      sscores <- ssgseaM[which(ssgseaM$variable == GS),]
      sscores <- sscores %>%
        select(Sample,type,Matrix,value)
      colnames(sscores)[2] <- "Type"
      colnames(sscores)[4] <- GS
      write_tsv(sscores, file)
    }
  )
  
  #render download button for reciprocal ssGSEA table - 4
  output$ssGSEADown4 <- downloadHandler(
    filename = function() {
      M1 <- input$MAT1name
      M2 <- input$MAT2name
      ssType <- input$ssGSEAtype
      #master table
      ssgseaM <- ssGSEAmaster()
      #select GS
      sets <- as.vector(unique(ssgseaM[,3]))
      GS <- sets[4]
      paste(GS,ssType,"Score.tsv",sep = "_")
    },
    content = function(file) {
      #master table
      ssgseaM <- ssGSEAmaster()
      #select GS
      sets <- as.vector(unique(ssgseaM[,3]))
      GS <- sets[4]
      #reformat
      sscores <- ssgseaM[which(ssgseaM$variable == GS),]
      sscores <- sscores %>%
        select(Sample,type,Matrix,value)
      colnames(sscores)[2] <- "Type"
      colnames(sscores)[4] <- GS
      write_tsv(sscores, file)
    }
  )
  
  #render download button for example file 1 - scatter
  output$ExampDownload1.s <- downloadHandler(
    
    filename = function() {
      paste(examp_file1.s, sep = '')
    },
    content = function(file) {
      examp1.s <- read.delim(examp_file1.s, sep = '\t', header = T)
      write_tsv(examp1.s, file)
    }
  )
  
  #render download button for example file  2 - scatter
  output$ExampDownload2.s <- downloadHandler(
    
    filename = function() {
      paste(examp_file2.s, sep = '')
    },
    content = function(file) {
      examp2.s <- read.delim(examp_file2.s, sep = '\t', header = T)
      write_tsv(examp2.s, file)
    }
  )
  
  #render download button for example file 1 - cor
  output$ExampDownload1.c <- downloadHandler(
    
    filename = function() {
      paste(examp_file1.c, sep = '')
    },
    content = function(file) {
      examp1.c <- read.delim(examp_file1.c, sep = '\t', header = T)
      write_tsv(examp1.c, file)
    }
  )
  
  #render download button for example file  2 - cor
  output$ExampDownload2.c <- downloadHandler(
    
    filename = function() {
      paste(examp_file2.c, sep = '')
    },
    content = function(file) {
      examp2.c <- read.delim(examp_file2.c, sep = '\t', header = T)
      write_tsv(examp2.c, file)
    }
  )
  
  #render download button for RNAvProt logFC table
  output$logFCtableDownload <- downloadHandler(
    
    filename = function() {
      M1 <- input$MAT1name
      M2 <- input$MAT2name
      paste(M1,"_",M2,"_ExprLogFC.tsv", sep = '')
    },
    content = function(file) {
      FCtable <- RNAvProtFC()
      M1 <- input$MAT1name
      M2 <- input$MAT2name
      ##add difference column in df
      FCtable$Delta_M1M2 <- FCtable[,2] - FCtable[,3]
      colnames(FCtable)[4] <- paste("Delta_",M1,"_",M2,sep = "")
      colnames(FCtable)[c(2,3)] <- c(paste("Log2FC_",M1,sep = ""),paste("Log2FC_",M2,sep = ""))
      write_tsv(FCtable, file)
    }
  )
  
  #render download button for example file  2 - cor
  output$ExampDownload3.R <- downloadHandler(
    
    filename = function() {
      paste(examp_file3.R, sep = '')
    },
    content = function(file) {
      examp3.R <- read.delim(examp_file3.R, sep = '\t', header = T)
      write_tsv(examp3.R, file)
    }
  )
  
  #render download button for example file  2 - cor
  output$ExampDownload3.RM <- downloadHandler(
    
    filename = function() {
      paste(examp_file3.RM, sep = '')
    },
    content = function(file) {
      examp3.RM <- read.delim(examp_file3.RM, sep = '\t', header = T)
      write_tsv(examp3.RM, file)
    }
  )
  
  #render download button for example file  3 - protein expression
  output$ExampDownload3.P <- downloadHandler(
    
    filename = function() {
      paste(examp_file3.P, sep = '')
    },
    content = function(file) {
      examp3.P <- read.delim(examp_file3.P, sep = '\t', header = T)
      write_tsv(examp3.P, file)
    }
  )
  
  #render download button for example file  3 - protein meta
  output$ExampDownload3.PM <- downloadHandler(
    
    filename = function() {
      paste(examp_file3.PM, sep = '')
    },
    content = function(file) {
      examp3.PM <- read.delim(examp_file3.PM, sep = '\t', header = F)
      write_tsv(examp3.PM, file)
    }
  )
  
}

# Run the application 
shinyApp(ui = ui, server = server)


