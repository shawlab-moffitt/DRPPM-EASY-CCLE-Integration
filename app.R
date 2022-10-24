


####----Install and load packages----####

packages <- c("shiny","shinythemes","shinyjqui","pheatmap","BiocManager",
              "dplyr","DT","ggplot2","ggpubr","reshape2","tibble",
              "plotly","readr","enrichR","ggrepel","tidyr","tools","shinycssloaders")

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



####----User Data Input----####

##--Project Name--##
ProjectName <- "CCLE"

##--Database Files--##
#Meta
db_meta_file <- 'CCLE_data/CCLE_meta.zip'
#Meta Selector File
db_meta_selec_file <- 'CCLE_data/CCLE_meta_selector.tsv'
#Expression Data
db_expr_file <- 'CCLE_data/CCLE_expr.zip'
#Name Map File
db_namemap_file <- 'CCLE_data/CCLE_NameMap.tsv'

##--User Gene Set Input--##
#write in the name of your gene set list for shiny UI
userGSlist_name <- 'CellMarker Gene Sets'
#path to your gene set file .gmt or .txt/.tsv
userGS_file <- 'CellMarker_gsNsym_HS.tsv'
#Does gene set file have header?
header.gs <- TRUE
#path to your R data list object for ssGSEA
userRData_file <- 'CellMarker_GS_HS.RData'



####----Backend Data Input----####

#MSigDB gene set
msigdb <- 'msigdb_gsNsym_HS.zip'
#MSigDB gene set FOR UI
msigdb2 <- 'msigdb_gsNcat_HS.tsv'
#gene set list for ssGSEA
load('msigdb_gs_HS.RData')
#Cytokine genes for human
CTKgenes <- c("IL2","IL12A","IL12B","IL17A","IFNA1","IFNB1","IFNG","IFNGR","CD11b",
              "ITGAM","CD33","ENTPD1","ICOSLG","CD275","CD278","TNFSF9","TNFRSF9",
              "CD40","CD40LG","CD70","CD27","TNFSF18","TNFRSF18","TNFSF14","TNFRSF14",
              "TNFSF4","TNFRSF4","HLA-A","CD3","CEACAM1","CD80","CD86","CTLA4","CD276",
              "VTCN1","PVR","CD226","TIGIT","CD96","LGALS3","LGALS3BP","LGALS9","LGALS9C",
              "HAVCR2","HHLA2","TMIGD2","CD274","PDCD1LG2","PDCD1","VSIR")

#Begin reactive values NULL
Gene <- NULL
geneList <- NULL
sampsames <- NULL
metagroups <- NULL
boxopt <- NULL


####----Read in Data and Reformat----####

##--MSigDB Gene Sets--##
#MSigDB gene sets
msigdb.gsea <- as.data.frame(read_delim(msigdb, delim = '\t'))
gmt <- msigdb.gsea
#MSigDB gene sets FOR UI
msigdb.gsea2 <- read.delim(msigdb2, header = T, sep = '\t')

##--Tab 2 'user' gene set--##
#tab2 User gene set
if (file_ext(userGS_file) == "gmt") {
  tab2 <- read.gmt(userGS_file)
}
if (file_ext(userGS_file) == "tsv" || file_ext(userGS_file) == "txt") {
  tab2 <- read.delim(userGS_file, header = header.gs, sep = '\t')
}
#tab2 back end
GeneSet2 <- as.data.frame(unique(tab2[,1]))
rownames(GeneSet2) <- 1:nrow(GeneSet2)
colnames(GeneSet2)[1] <- "Gene_Set"
#tab2 R Data list
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
gs2 <- loadRData(userRData_file)

##--Large Project Files--##
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



####----Other information and functions----####

#CV function for variance
cv <- function(x){
  (sd(x)/mean(x))*100
}

#enrichR pathway choices
enrichRtab <- listEnrichrDbs()
enrichRchoice <- enrichRtab[,"libraryName"]





####----Shiny UI----####


shinytheme("sandstone")

ui <-
  
  navbarPage(paste("{ ",ProjectName," Expression Analysis }", sep=""),
             
             ####----Intro Tab----####
             
             tabPanel("Intro/Methodology",
                      fluidPage(
                        mainPanel(
                          tabsetPanel(
                            id = 'intropage',
                            tabPanel("Intoduction",
                                     h3("Introduction"),
                                     p("This is an extention of the DRPPM Expression Analysis ShinY (EASY) App which integrates data from the Cancer Cell Line Encyclopedia (CCLE) for sample selection perform differential gene expression and gene set enrichment analyses. This R Shiny app is very similar in features to the main app, expect for the addition of the first tab which allows for CCLE sample selection. Based on this selection, the expression and meta data will be subset and imported in the back end to the app for further analysis and visualization. Below is a schematic for the DRPPM-EASY-CCLE pipeline."),
                                     h3("Unsupervised Clustering Methods"),
                                     p("Unsupervised clustering is performed by calculating the top variable genes through MAD, CV, or VAR as the variance measure. The choice of variance measure and number to top genes/probes can be determined by the user, as well as the number of clusters with k-means. The identified most variable genes can be downloaded as a table for the user. The â€œcompleteâ€ method was chosen as the default algorithm for clustering, but other methods such as Ward, average, and centroid could be chosen based on the â€˜hclust()â€™ function in R."),
                                     h3("Differential Gene Expression Methods"),
                                     p("Differential gene expression analysis is performed on the expression data between two groups defined by the provided metadata file. The samples chosen are log-transformed (log2 + 1). A model matrix is then designed for LIMMA linear modeling followed by empirical Bayes statistics to moderate gene variance and modeling of the global characteristic of the data. Of note, the current pipeline expects the matrix input quantification are pre-normalized, such as in CPM, TMM, or FPKM."),
                                     h3("Gene Set Enrichment Analysis (GSEA) Methods"),
                                     p("Gene Set Enrichment Analysis (GSEA) is performed through signal-to-noise calculations on the expression matrix. This calculation requires at least two types of sample groups and at least three samples for each grouping. The signal-to-noise calculation scales the difference of means by standard deviation and generates a ranked list of genes, effectively capturing the differences between phenotypes. The GSEA identifies the enrichment score (ES) for a gene set, which indicates the extent that the gene set is overrepresented at the beginning or end of the list of ranked genes. The enrichment score is calculated by walking down the ranked gene list and increasing the running-sum statistic when a gene is in the gene set and decreasing when it is not. The maximum deviation from zero that is encountered through walking the list is the ES, where a positive ES indicates the gene set is enriched at the top of the ranked list and a negative ES indicates the gene set is enriched at the bottom of the ranked list. A leading-edge gene list is provided displaying a subset of the genes in the gene set that contribute the most to the ES."),
                                     h4("Gene Set Sources"),
                                     p("The Molecular Signatures Database (MSigDB) is a collection of over 32,000 annotated gene sets divided into 9 major collections as well as various sub-collections [PMID: 16199517, PMID: 12808457]. The MSigDB  gene sets were downloaded with the msigdbr package and processed through R using Tidyverse packages to combine the gene sets into a data frame."),
                                     p("The Cell Marker gene sets were derived from a comprehensive list of cell markers from cell types in various tissues obtained through manually curating over 100,000 published papers [PMID: 30289549]. This data was obtained as an extensive data frame and parsed to make gene sets with the tidyverse packages."),
                                     p("The Library of Integrated Network-based Cellular Signatures (LINCS) L1000 gene sets were derived from the Connectivity Map (CMAP) projects data which treated 4 human cell lines with ~1300 drugs followed by profiling genome-wide mRNA expression following small molecule perturbation [PMID: 24906883].")
                            ),
                            tabPanel("Data Selection Tutorial",
                                     h3("Subsetting Expression Data:"),
                                     p("This is an optional selection depending on the large project data that is being analyzed. If there is a 'MetaSelector' file, which is described on the GitHub page, this option to subset the expression data will appear. This is useful if the expression data can be subset be a certain parameter, for instance in the example of a CCLE data set the samples come from various different diseases or Cancer Linages, so the selector file and subsetting the expression data allows the user to narrow down which subset is of interest."),
                                     h3("Selecting a Coniditon of Interest"),
                                     p("After subsetting the expression data, if that is an option, the user must select a condition/phentotype from the drop down box which is generated based on the meta data."),
                                     h3("Other Options:"),
                                     p("The user may choose to log2 transform the selected expression data as well as name the selected data being analyzed which will be used as a label throughout the analysis."),
                                     h3("Final Step:"),
                                     p("Lastly the user must press the 'Load/Update Selected Data' button to load the selected data into the program."),
                                     p(),
                                     p("You may find more information on our ",a("GitHub page.", href="https://github.com/shawlab-moffitt/DRPPM-EASY-LargeProject"))
                            )
                          ))
                      )
             ),
             
             ####----Sample Selection Tab----####
             tabPanel("Sample Selection",
                      fluidPage(
                        sidebarPanel(
                          width = 3,
                          uiOutput("subexprLabel"),
                          uiOutput("subselectexpr"),
                          uiOutput("subsetselec"),
                          h4("Condition Selection:"),
                          selectInput("phenoChoice","",
                                      choices = db_pheno_choice, selected = "-"),
                          checkboxInput("log2expr","Log2 Transform Expression Data:"),
                          h4("Label for Project Data:"),
                          textInput("MAT2name","", value = ProjectName),
                          actionButton("updateMETA", "Load/Update Selected Data")
                        ),
                        mainPanel(
                          fluidRow(
                            column(4,
                                   uiOutput("metadatalabel"),
                                   withSpinner(div(DT::dataTableOutput("testtable"), style = "font-size:10px; height:450px"), type = 6),
                                   uiOutput("downloadMETAbutton")),
                            column(4,
                                   uiOutput("exprdatalabel"),
                                   withSpinner(div(DT::dataTableOutput("testexprtable"), style = "font-size:10px; height:450px"), type = 6),
                                   uiOutput("downloadEXPRbutton")),
                            column(4,
                                   uiOutput("NameMapLabel"),
                                   withSpinner(div(DT::dataTableOutput("NameMapTable"), style = "font-size:10px; height:450px"), type = 6),
                                   uiOutput("downloadNameMapbutton"))
                          )
                        )
                      )
             ),
             
             ####----Data Exploration Tab----####
             
             tabPanel("Data Exploration",
                      fluidPage(
                        title = "Data Exploration",
                        sidebarLayout(
                          sidebarPanel(
                            
                            ##--Unsupervised Clustering Heatmap--##
                            
                            conditionalPanel(condition = "input.dataset == '22'",
                                             p(),
                                             numericInput("NumFeatures", step = 1, label = "Number of Genes", value = 100),
                                             h4("Clustering Parameters"),
                                             selectInput("VarianceMeasure", "Select Variance Measure",
                                                         choices = c("MAD","CV","VAR")),
                                             selectInput("ClusteringMethod",
                                                         "Select clustering Method",
                                                         choices = c("complete", "ward.D", "ward.D2", "single", "average", "mcquitty", "median", "centroid")),
                                             numericInput("NumClusters", step = 1, label = "Number of Clusters (Cut Tree with ~k)", value = 2),
                                             h4("Download Cluster Result:"),
                                             downloadButton("downloadClusters", "Download .tsv"),
                                             h4("Download Most Variable Genes List:"),
                                             downloadButton("MVGdownload", "Download MVG .tsv"),
                                             downloadButton("MVGdownloadgmt", "Download MVG .gmt"),
                                             p(),
                                             sliderInput("heatmapFont2.c", "Heatmap Font Size (columns):",
                                                         min = 5, max = 75,
                                                         value = 12, step = 1),
                                             sliderInput("heatmapFont2.r", "Heatmap Font Size (rows):",
                                                         min = 5, max = 75,
                                                         value = 9, step = 1)
                            ),
                            
                            ##--Custom Heatmap--##
                            
                            conditionalPanel(condition = "input.dataset == '33'",
                                             p(),
                                             selectizeInput("heatmapGeneSelec","Gene Selection:",
                                                            choices = sort(as.vector(geneList[,1])),
                                                            multiple = T, selected = CTKgenes),
                                             textInput("userheatgenes", "Text Input of Gene List (space or tab delimited):", value = ""),
                                             selectInput("ClusteringMethod2",
                                                         "Select clustering Method",
                                                         choices = c("complete", "ward.D", "ward.D2", "single", "average", "mcquitty", "median", "centroid")),
                                             selectizeInput("userheatsamp2", "Samples Selection:",
                                                            choices = sampsames, multiple = T, selected = sampsames),
                                             sliderInput("heatmapFont3.c", "Heatmap Font Size (columns):",
                                                         min = 5, max = 75,
                                                         value = 12, step = 1),
                                             sliderInput("heatmapFont3.r", "Heatmap Font Size (rows):",
                                                         min = 5, max = 75,
                                                         value = 12, step = 1)
                            ),
                            
                            ##--Gene scatter plot side panel--##
                            
                            conditionalPanel(condition = "input.dataset == '26'",
                                             p(),
                                             selectInput("scatterG1","Select Gene 1",
                                                         choices = Gene),
                                             selectInput("scatterG2","Select Gene 2",
                                                         choices = Gene, selected = Gene[2]),
                                             checkboxInput("logask","Log2 Transform Expression Data")
                            ),
                            
                            ##--Boxplot side panel--##
                            
                            conditionalPanel(condition = "input.dataset == '34'",
                                             p(),
                                             h3("Select A Gene:"),
                                             div(DT::dataTableOutput("GeneListTable2"), style = "font-size:10px"),
                                             p(),
                                             sliderInput("boxplot.font2", "Font Size:",
                                                         min = 5, max = 75,
                                                         value = 15, step = 1),
                                             sliderInput("boxplotDot2", "Box Plot Dot Size:",
                                                         min = 0, max = 5,
                                                         value = 0.75, step = .25)
                            )
                          ),
                          mainPanel(
                            tabsetPanel(
                              id = "dataset",
                              tabPanel("Unsupervised Clustering Heatmap",
                                       withSpinner(jqui_resizable(plotOutput("heatmap1", width = "100%", height = "1000px")), type = 6),
                                       value = 22),
                              tabPanel("Custom Heatmap",
                                       withSpinner(jqui_resizable(plotOutput("heatmap2", width = "100%", height = "1000px")), type = 6),
                                       value = 33),
                              tabPanel("Gene Scatter Plot",
                                       p(),
                                       withSpinner(jqui_resizable(plotlyOutput("geneScatter0", height = "500px")), type = 6),
                                       DT::dataTableOutput("geneScatterTable"),
                                       downloadButton("geneScatterDownload", "Download Non-log2 Transformed .tsv"),
                                       value = 26),
                              tabPanel("Box Plots",
                                       p(),
                                       withSpinner(jqui_resizable(plotOutput('boxplot3', width = "100%", height = "600px")), type = 6),
                                       value = 34)
                            )
                          )
                        )
                      )
             ),
             
             ####----Differential Expression Analysis----####
             
             tabPanel("Differential Expression Analysis",
                      fluidPage(
                        title = "Differential Expression Analysis",
                        sidebarLayout(
                          sidebarPanel(
                            
                            ##--Volcano and MA plot side panel--##
                            
                            conditionalPanel(condition = "input.datasettwo == '44' || input.datasettwo == '66'",
                                             p(),
                                             selectInput("comparisonA2", "Comparison: GroupA",
                                                         choices = metagroups, selected = metagroups[1]),
                                             selectInput("comparisonB2", "Comparison: GroupB",
                                                         choices = metagroups, selected = metagroups[2]),
                                             numericInput("fc_cutoff", "LogFC Threshold (Absolute Value)",
                                                          min = 0, max = 5, step = 0.1, value = 1),
                                             numericInput("p_cutoff", "Significance Threshold (-log10(P.Value)):",
                                                          min = 0, max = 10, step = 0.1, value = 0.05),
                                             numericInput("top_x", "Number of Top Hits:", value = 10,
                                                          min = 0),
                                             selectizeInput("userGeneSelec", "User Selected Hits:",
                                                            choices = sort(as.vector(geneList[,1])), multiple = T, selected = "-"),
                                             textInput("userGeneSelec2", "Text Input of Gene List (space or tab delimited):", value = "")
                            ),
                            
                            ##--Boxplot side panel--##
                            
                            conditionalPanel(condition = "input.datasettwo == '88'",
                                             p(),
                                             h3("Select A Gene:"),
                                             div(DT::dataTableOutput("GeneListTable"), style = "font-size:10px"),
                                             p(),
                                             sliderInput("boxplot.font", "Font Size:",
                                                         min = 5, max = 75,
                                                         value = 15, step = 1),
                                             sliderInput("boxplotDot", "Box Plot Dot Size:",
                                                         min = 0, max = 5,
                                                         value = 0.75, step = .25)
                            ),
                            
                            ##--EnrichR pathway side panel--##
                            
                            conditionalPanel(condition = "input.datasettwo == '28'",
                                             selectInput("comparisonA2.path", "Comparison: GroupA",
                                                         choices = metagroups, selected = metagroups[1]),
                                             selectInput("comparisonB2.path", "Comparison: GroupB",
                                                         choices = metagroups, selected = metagroups[2]),
                                             numericInput("pathpval", "Adjusted P-Value Cutoff:",
                                                          min = 0, value = 0.05),
                                             numericInput("pathFC", "Log2 Fold Change Cutoff:",
                                                          min = 0, value = 1),
                                             selectInput("SelectedPathway", "Select Pathway",
                                                         choices = enrichRchoice)
                            ),
                            
                            ##--DEG table side panel--##
                            
                            conditionalPanel(condition = "input.datasettwo == '24'",
                                             p(),
                                             selectInput("comparisonA2.DEG", "Comparison: GroupA",
                                                         choices = metagroups, selected = metagroups[1]),
                                             selectInput("comparisonB2.DEG", "Comparison: GroupB",
                                                         choices = metagroups, selected = metagroups[2]),
                                             h4("Download DEG Table as GMT File"),
                                             textInput("DEGfileName", "File Name for Download:",value = "DEGgeneSet"),
                                             numericInput("fc_cutoff2", "LogFC Threshold",
                                                          min = 0, max = 5, step = 0.1, value = 1),
                                             selectInput("UpDnChoice","Up-regulated or Down-regulated:",
                                                         choices = c("UpAndDown_Regulated","Up_Regulated","Down_Regulated")),
                                             numericInput("p_cutoff2", "Adj.P.Value Cutoff:",
                                                          min = 0, max = 10, step = 0.1, value = 0.05),
                                             numericInput("top_x2", "Number of Top Hits:", value = 100),
                                             downloadButton("DEGgmtDownload", "Download DEG .gmt")
                            )
                          ),
                          mainPanel(
                            tabsetPanel(
                              id = "datasettwo",
                              tabPanel("Volcano Plot",
                                       p(),
                                       verbatimTextOutput("VolGroupsText"),
                                       jqui_resizable(plotOutput('Volcano3', width = "800px", height = "550px",
                                                                 hover = hoverOpts("plot_hover", delay = 10, delayType = "debounce"))),
                                       uiOutput("hover_info"),
                                       value = 44),
                              tabPanel("MA Plot",
                                       p(),
                                       verbatimTextOutput("MAGroupsText"),
                                       jqui_resizable(plotOutput('MAPlot1', width = "800px", height = "550px",
                                                                 hover = hoverOpts("plot_hover2", delay = 10, delayType = "debounce"))),
                                       uiOutput("hover_info2"),
                                       value = 66),
                              tabPanel("Box Plots",
                                       p(),
                                       withSpinner(jqui_resizable(plotOutput('boxplot1', width = "100%", height = "600px")), type = 6),
                                       value = 88),
                              tabPanel("Pathway Analysis",
                                       uiOutput("UpRegPathLabel"),
                                       verbatimTextOutput("upregpath_text"),
                                       withSpinner(plotOutput('UpRegPathway1'), type = 6),
                                       div(DT::dataTableOutput("UpRegPathwayTable1"), style = "font-size:10px"),
                                       downloadButton("UpRegPathDownload", "Download Table .tsv"),
                                       downloadButton("UpRegPathDownloadgmt", "Download .gmt"),
                                       uiOutput("DnRegPathLabel"),
                                       verbatimTextOutput("downregpath_text"),
                                       withSpinner(plotOutput('DnRegPathway1'), type = 6),
                                       div(DT::dataTableOutput("DnRegPathwayTable1"), style = "font-size:10px"),
                                       downloadButton("DnRegPathDownload", "Download Table .tsv"),
                                       downloadButton("DnRegPathDownloadgmt", "Download .gmt"),
                                       p(),
                                       value = 28),
                              tabPanel("DEG Table",
                                       p(),
                                       verbatimTextOutput("degtext"),
                                       div(DT::dataTableOutput("DEGtable1"), style = "font-size:12px"),
                                       downloadButton("DEGtableDownload", "Download DEG table .tsv"),
                                       value = 24)
                            )
                          )
                        )
                      )
             ),
             
             ####----GSEA Tab----####
             
             tabPanel("GSEA Analysis",
                      fluidPage(
                        title = "GSEA Analysis",
                        sidebarLayout(
                          sidebarPanel(
                            tabsetPanel(id = "GSEA",
                                        tabPanel("GSEA Parameters",
                                                 p(),
                                                 selectInput("comparisonA", "Comparison: GroupA",
                                                             choices = metagroups, selected = metagroups[1]),
                                                 selectInput("comparisonB", "Comparison: GroupB",
                                                             choices = metagroups, selected = metagroups[2]),
                                                 numericInput("userPval", "Pvalue Cutoff", value = 1.0, width = '100%'),
                                                 h4("ssGSEA Parameters"),
                                                 selectInput("ssGSEAtype","Choose ssGSEA Method",
                                                             choices = c("ssgsea","gsva","zscore","plage")),
                                                 selectInput("boxplotcompare", "Boxplot Stat Compare Method:",
                                                             choices = boxopt),
                                                 tabsetPanel(
                                                   id = "tables",
                                                   tabPanel("MSigDB Gene Sets",
                                                            div(DT::dataTableOutput("msigdbTable"), style = "font-size:10px; height:500px; overflow-X: scroll"),
                                                            value = 1),
                                                   tabPanel(paste(userGSlist_name),
                                                            div(DT::dataTableOutput("tab2table"), style = "font-size:10px; height:500px; overflow-X: scroll"),
                                                            value = 3),
                                                   tabPanel("Use your own gene set",
                                                            p(),
                                                            uiOutput("user.gmt"),
                                                            uiOutput("user.RDataButton"),
                                                            uiOutput("RDataMessage"),
                                                            uiOutput("user.GStable"),
                                                            value = 5)
                                                 )
                                        ),
                                        tabPanel("Figure Parameters",
                                                 p(),
                                                 h4("Heatmap Parameters"),
                                                 sliderInput("heatmapFont1.c", "Heatmap Font Size (columns):",
                                                             min = 5, max = 75,
                                                             value = 12, step = 1),
                                                 sliderInput("heatmapFont1.r", "Heatmap Font Size (rows):",
                                                             min = 5, max = 75,
                                                             value = 10, step = 1),
                                                 h4("Box Plot Parameters"),
                                                 sliderInput("boxplotFontss", "Box Plot Font Size:",
                                                             min = 5, max = 75,
                                                             value = 15, step = 1),
                                                 sliderInput("boxplotDotss", "Box Plot Dot Size:",
                                                             min = 0, max = 5,
                                                             value = 0.75, step = .25)
                                        )
                            )
                          ),
                          mainPanel(
                            tabsetPanel(
                              id = "datasettwo",
                              tabPanel("Enrichment Plot",
                                       h3("GSEA Enrichment Plot"),
                                       verbatimTextOutput("NESandPval"),
                                       withSpinner(plotOutput("enrichplot0", width = "500px", height = "450px"), type = 6),
                                       h3("Leading Edge Genes (~Signal2Noise Ranking)"),
                                       downloadButton("LEGdownload", "Download .tsv"),
                                       div(DT::dataTableOutput("LeadingEdgeGenes"), style = "font-size:12px; height:500px; overflow-y: scroll"),
                                       value = 2),
                              tabPanel("Gene Expression Heatmap",
                                       withSpinner(jqui_resizable(plotOutput("heatmap0", width = "100%", height = "1000px")), type = 6),
                                       value = 4),
                              tabPanel("ssGSEA Boxplots",
                                       p(),
                                       withSpinner(jqui_resizable(plotOutput('boxplot2', width = "100%", height = "500px")), type = 6),
                                       DT::dataTableOutput("ssGSEAtable"),
                                       downloadButton("ssGSEAdownload", "Download .tsv"),
                                       value = 10)
                            )
                          )
                        )
                      )
             )
  )

####----Server----####


server <- function(input, output, session) {
  
  
  ####----Sample Selection and Data Initialization----####
  
  
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
  
  #observe user selection to choose subset expression data, update UI based on selection
  observe({
    output$subsetselec <- renderUI({
      
      if (length(subsetter_list) > 1){
        
        firstChoice_options <- input$SelecExpr
        
        selectInput("firstChoice","Select Sample Subset:",choices = subsetter_list[[firstChoice_options]],
                    multiple = F)
        
      }
    })
  })
  
  ###--Tim added the following function 1--###
  
  # initialize values
  updated_expr <- reactiveVal()
  meta <- reactiveVal()
  geneList <- reactiveVal()
  
  observeEvent(input$updateMETA, {
    
    ##assign user input values for meta selector choice
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
    
    #replace meta variable with new subset meta
    meta(as.data.frame(db_meta_sub))
    #assign the new meta to a variable
    meta = meta()
    #determine number of levels in meta sample type to provide specific options to the user
    metagroups <- as.vector(levels(factor(meta[,2])))
    #backend to determine type of stat compare method for box plot
    if (length(metagroups) == 2) {
      boxopt <- c("wilcox.test", "t.test", "none")
    }
    if (length(metagroups) >= 3) {
      boxopt <- c("kruskal.test", "anova", "none")
    }
    
    #update all user selections in app based on new data
    updateSelectInput(session, "boxplotcompare", "Boxplot Stat Compare Method:", choices = boxopt)
    updateSelectInput(session, "comparisonA2", "Comparison: GroupA", choices = metagroups, selected = metagroups[1])
    updateSelectInput(session, "comparisonA2.path", "Comparison: GroupA", choices = metagroups, selected = metagroups[1])
    updateSelectInput(session, "comparisonA2.DEG", "Comparison: GroupA", choices = metagroups, selected = metagroups[1])
    updateSelectInput(session, "comparisonA", "Comparison: GroupA", choices = metagroups, selected = metagroups[1])
    updateSelectInput(session, "comparisonB2", "Comparison: GroupB", choices = metagroups, selected = metagroups[2])
    updateSelectInput(session, "comparisonB2.path", "Comparison: GroupB", choices = metagroups, selected = metagroups[2])
    updateSelectInput(session, "comparisonB2.DEG", "Comparison: GroupB", choices = metagroups, selected = metagroups[2])
    updateSelectInput(session, "comparisonB", "Comparison: GroupB", choices = metagroups, selected = metagroups[2])
    
    
    ##--update EXPR section--##
    firstChoice <- input$firstChoice
    phenoChoice <- input$phenoChoice
    exprChoice <- gsub("_Choices","",input$SelecExpr)
    
    samp_selec <- unlist(meta()[,1])
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
    
    
    # gene list file from expression data
    Gene <- rownames(expr)
    geneList <- as.data.frame(Gene)
    #update gene list reactive value with genes from the subset data
    geneList(geneList)
    
    
    #meta
    meta <- meta();
    meta[,1] <- gsub("[_.-]", ".", meta[,1])
    colnames(meta) <- c("SampleName","Group")
    metagroups <- as.vector(levels(factor(meta[,2])))
    
    #boxplot choices based on meta groups
    if (length(metagroups) == 2) {
      boxopt <- c("wilcox.test", "t.test", "none")
    }
    if (length(metagroups) >= 3) {
      boxopt <- c("kruskal.test", "anova", "none")
    }
    
    #for heatmap sample selection
    sampsames <- intersect(colnames(expr),meta[,1])
    #ensure expression samples and meta are exact
    expr <- expr[,sampsames]
    meta = meta[which(meta[,1] %in% sampsames),]
    meta <- meta %>%
      distinct()
    #update expression and meta reactive values ensure both have the same samples
    updated_expr(expr)
    meta(meta)
    
    #update gene selection based on genes available
    updateSelectInput(session, "scatterG1","Select Gene 1", choices = Gene)
    updateSelectInput(session, "scatterG2","Select Gene 2", choices = Gene, selected = Gene[2])
    updateSelectInput(session, "heatmapGeneSelec","Gene Selection:",
                      choices = sort(as.vector(geneList[,1])),
                      selected = CTKgenes)
    updateSelectInput(session, "userheatsamp2", "Samples Selection:",
                      choices = sampsames, selected = sampsames)
    updateSelectInput(session, "userGeneSelec", "User Selected Hits:",
                      choices = sort(as.vector(geneList[,1])), selected = "-")
    
    
  })
  output$testexprtable <- DT::renderDataTable({
    expr = updated_expr();
    DT::datatable(expr,
                  options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"), scrollX = T))
    
  })
  output$testtable <- DT::renderDataTable({
    DT::datatable(meta(),
                  rownames = F,
                  options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"), scrollX = T))
    
  })
  
  output$metadatalabel <- renderUI({
    
    req(input$updateMETA)
    h4("Meta Data")
    
  })
  
  output$exprdatalabel <- renderUI({
    
    req(input$updateMETA)
    h4("Expression Data")
    
  })
  
  output$NameMapLabel <- renderUI({
    
    if (file.exists(db_namemap_file) == TRUE) {
      req(input$updateMETA)
      h4("Sample Name Guide")
    }
    
  })
  
  output$NameMapTable <- DT::renderDataTable(({
    req(input$updateMETA)
    namemap <- db_namemap
    meta.new <- meta()
    samp_selected <- unlist(meta.new[,1])
    namemap_selcted <- namemap[which(namemap[,1] %in% samp_selected),]
    DT::datatable(namemap_selcted,
                  rownames = F,
                  options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"), scrollX = T))
    
    
  }))
  
  output$downloadNameMapbutton <- renderUI({
    req(input$updateMETA)
    downloadButton("downloadNameMap","Download Cell Line Name Map")
  })
  
  output$downloadNameMap <- downloadHandler(
    filename = function() {
      firstChoice <- input$firstChoice
      phenoChoice <- input$phenoChoice
      exprChoice <- gsub("_Choices","",input$SelecExpr)
      paste(ProjectName,"_NameMap_",gsub(" ","",firstChoice),"_",phenoChoice,".tsv", sep = "")
    },
    content = function(file) {
      namemap <- db_namemap
      meta.new <- meta()
      samp_selected <- unlist(meta.new[,1])
      namemap_selcted <- namemap[which(namemap[,1] %in% samp_selected),]
      write_tsv(namemap_selcted,file)
      
    }
  )
  
  output$downloadMETAbutton <- renderUI({
    req(input$updateMETA)
    downloadButton("downloadMETA","Download Meta Data")
  })
  
  output$downloadEXPRbutton <- renderUI({
    req(input$updateMETA)
    downloadButton("downloadEXPR","Download Expression Data")
  })
  
  
  output$downloadMETA <- downloadHandler(
    filename = function() {
      firstChoice <- input$firstChoice
      phenoChoice <- input$phenoChoice
      exprChoice <- gsub("_Choices","",input$SelecExpr)
      paste("meta_",ProjectName,"_",phenoChoice,".tsv",sep = "")
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
      paste("expr_",ProjectName,"_",phenoChoice,".tsv",sep = "")
    },
    content = function(file) {
      
      ##--update EXPR section--##
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
      db_meta_sub <- db_expr[,which(colnames(db_expr) %in% samp_selec), drop = F]
      db_meta_sub$gene <- rownames(db_meta_sub)
      db_meta_sub <- db_meta_sub %>%
        relocate(gene)
      write_tsv(db_meta_sub,file)
      
    }
  )
  
  #render user gmt data upload if indicated
  output$user.gmt <- renderUI({
    fileInput("user.gmt.file", "Gene Set File (.gmt, .tsv, or .txt)", accept = c(".gmt",".tsv",".txt"))
  })
  
  #render gene set table based off gmt file given
  output$user.GStable <- renderUI({
    req(input$user.gmt.file)
    div(DT::dataTableOutput("GStable.u"), style = "font-size:10px; height:500px; overflow-X: scroll")
  })
  
  #Warning message for RData list generation
  output$RDataMessage <- renderUI({
    req(input$user.gmt.file)
    p("Generating RData list may take up to several minutes depending on size of GMT file.")
  })
  
  #render action button to create RData list for ssGSEA boxplots
  output$user.RDataButton <- renderUI({
    req(input$user.gmt.file)
    actionButton("user.RData.Gen", "Generate RData list for ssGSEA Boxplots")
  })
  
  #render UI for hover text in volcano plot
  output$hover_info <- renderUI({
    top2 <- topgenereact()
    df <- top2 %>%
      select(GeneName,logFC,P.Value,adj.P.Val)
    colnames(df)[3] <- "-log10(P.Value)"
    df$P.Value <- df$'-log10(P.Value)'
    df$`-log10(P.Value)` <- -log10(df$`-log10(P.Value)`)
    hover <- input$plot_hover
    point <- nearPoints(df, hover, threshold = 10, maxpoints = 1, addDist = FALSE)
    if (nrow(point) == 0) return(NULL)
    wellPanel(
      p(HTML(paste0("<b> Name: </b>", point[1], "<br/>",
                    "<b> Fold change: </b>", round(point[2], digits = 4), "<br/>",
                    "<b> P Value: </b>", point[5], "<br/>",
                    "<b> adj P Value: </b>", point[4], "<br/>",
                    NULL
      ))))
  })
  
  #render UI for hover text in MA plot
  output$hover_info2 <- renderUI({
    top2 <- topgenereact()
    df <- top2 %>%
      select(GeneName,AveExpr,logFC,P.Value,adj.P.Val)
    colnames(df)[4] <- "-log10(P.Value)"
    df$P.Value <- df$'-log10(P.Value)'
    df$`-log10(P.Value)` <- -log10(df$`-log10(P.Value)`)
    hover <- input$plot_hover2
    point <- nearPoints(df, hover, threshold = 10, maxpoints = 1, addDist = FALSE)
    if (nrow(point) == 0) return(NULL)
    wellPanel(
      p(HTML(paste0("<b> Name: </b>", point[1], "<br/>",
                    "<b> Fold change: </b>", round(point[3], digits = 4), "<br/>",
                    "<b> P Value: </b>", point[6], "<br/>",
                    "<b> adj P Value: </b>", point[5], "<br/>",
                    "<b> Avg. Expression: </b>", round(point[2], digits = 4), "<br/>",
                    NULL
      ))))
  })
  
  
  ####----Reactives----####
  
  
  #reactive to generate RData list
  RDataListGen <- eventReactive(input$user.RData.Gen, {
    gmt <- GStable.ubg()
    colnames(gmt) <- c("term","gene")
    RData.u <- list()
    for (i in unique(gmt[,1])){
      RData.u[[i]] <- gmt[gmt[,1] == i,]$gene
    }
    RData.u
  })
  
  #reactive for ssGSEA function
  ssGSEAfunc <- reactive({
    expr = updated_expr()
    meta = meta()
    A = as.matrix(expr)
    
    if (input$tables == 1) {
      GS <- gs[(msigdb.gsea2[input$msigdbTable_rows_selected,3])]
    }
    if (input$tables == 3) {
      GS <- gs2[(GeneSet2[input$tab2table_rows_selected,1])]
    }
    if (input$tables == 5) {
      GS <- RDataListGen()[(user_gs_mirror()[input$GStable.u_rows_selected,1])]
    }
    gsva(A, GS, method = input$ssGSEAtype, verbose = F)
  })
  
  #create background GMT from user input gene set table
  GStable.ubg <- reactive({
    gmt.u <- input$user.gmt.file
    ext <- tools::file_ext(gmt.u$datapath)
    req(gmt.u)
    validate(need(ext == c("gmt","tsv","txt"), "Please upload .gmt, .tsv, or .txt file"))
    if (ext == "gmt") {
      read.gmt(gmt.u$datapath)
    }
    else {
      as.data.frame(read_delim(gmt.u$datapath, delim = '\t'))
    }
  })
  
  #gs mirror from user input for selection help on back end
  user_gs_mirror <- reactive({
    GeneSet <- as.data.frame(unique(GStable.ubg()[1]))
    rownames(GeneSet) <- 1:nrow(GeneSet)
    colnames(GeneSet)[1] <- "Gene_Set"
    GeneSet
  })
  
  #perform sig2noise calculation and create GSEA result from user chosen gene set
  datasetInput <- reactive({
    expr = updated_expr();
    A <- as.matrix(expr)
    meta = as.matrix(meta());
    
    groupA <- meta[which(meta[,2] == input$comparisonA),1]
    groupB <- meta[which(meta[,2] == input$comparisonB),1]
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
    ##----GSEA----##
    if (input$tables == 1){
      GSEA(s2n.matrix.s, TERM2GENE = gmt[which(gmt$gs_name == msigdb.gsea2[input$msigdbTable_rows_selected,3]),],
           verbose = F, pvalueCutoff = input$userPval)
    }
    else if (input$tables == 3){
      GSEA(s2n.matrix.s, TERM2GENE = tab2[which(tab2$term == as.character(GeneSet2[input$tab2table_rows_selected,1])),],
           verbose = F, pvalueCutoff = input$userPval)
    }
    else if (input$tables == 5){
      GSEA(s2n.matrix.s, TERM2GENE = GStable.ubg()[which(GStable.ubg()$term == as.character(user_gs_mirror()[input$GStable.u_rows_selected,1])),],
           verbose = F, pvalueCutoff = input$userPval)
    }
  })
  
  #top genes data frame reactive
  topgenereact <- reactive({
    meta = as.matrix(meta())
    expr = as.matrix(updated_expr())
    
    #make group based on user input
    A <- meta[which(meta[,2] == input$comparisonA2),1]
    B <- meta[which(meta[,2] == input$comparisonB2),1]
    
    #A <- meta[,1][meta[,2] == input$comparisonA2]
    #B <- meta[,1][meta[,2] == input$comparisonB2]
    #make top table
    samples <- c(A,B)
    mat <- expr[,which(colnames(expr) %in% samples)]
    mat <- log2(mat + 1.0)
    groupAOther <- factor(c(rep("A", length(A)), rep("B", length(B))))
    designA <- model.matrix(~0 + groupAOther)
    fit <- lmFit(mat, design = designA)
    contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    options(digits = 4)
    top1 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
    top2 <- top1
    top2["GeneName"] <- rownames(top2)
    Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
    top2['group'] <- "NotSignificant"
    top2[which(top2$P.Value < input$p_cutoff & abs(top2$logFC) < abs(input$fc_cutoff)), "group"] <- "Significant"
    top2[which(top2$P.Value > input$p_cutoff & abs(top2$logFC) > abs(input$fc_cutoff)), "group"] <- "FoldChange"
    top2[which(top2$P.Value < input$p_cutoff & abs(top2$logFC) > abs(input$fc_cutoff)), "group"] <- "Significant&FoldChange"
    top2['FCgroup'] <- "NotSignificant"
    top2[which(abs(top2$logFC) > abs(input$fc_cutoff)), "group2"] <- "FoldChange"
    top2
  })
  
  
  ####----Data Tables----####
  
  
  #render MSigDB gene set table
  output$msigdbTable <- DT::renderDataTable({
    DT::datatable(msigdb.gsea2,
                  selection = list(mode = 'single', selected = 1),
                  options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100")))
    
  })
  
  #create user input gene set table
  output$GStable.u <- DT::renderDataTable({
    gmt.u <- input$user.gmt.file
    ext <- tools::file_ext(gmt.u$datapath)
    req(gmt.u)
    validate(need(ext == c("gmt","tsv","txt"), "Please upload .gmt, .tsv, or .txt file"))
    if (ext == "gmt") {
      gmt.us <- read.gmt(gmt.u$datapath)
    }
    else {
      gmt.us <- as.data.frame(read_delim(gmt.u$datapath, delim = '\t'))
    }
    GeneSet <- as.data.frame(unique(gmt.us[1]))
    rownames(GeneSet) <- 1:nrow(GeneSet)
    colnames(GeneSet)[1] <- "Gene_Set"
    DT::datatable(GeneSet,
                  selection = list(mode = 'single', selected = 1),
                  options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100")))
  })
  
  #render tab2 gene set table
  output$tab2table <- DT::renderDataTable({
    DT::datatable(GeneSet2,
                  selection = list(mode = 'single', selected = 1),
                  options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100")))
    
  })
  
  #render leading edge genes list
  output$LeadingEdgeGenes <- DT::renderDataTable({
    if (input$tables == 1){
      if (length(input$msigdbTable_rows_selected) > 0){
        res <- datasetInput()
        gsea.df <- as.data.frame(res@result)
        GS <- msigdb.gsea2[input$msigdbTable_rows_selected,3]
        ## Subset core enriched genes
        genes1 <- as.matrix(gsea.df[which(gsea.df$Description==GS),"core_enrichment"])
        genes2 <- strsplit(genes1,"/")
        GeneSymbol <- as.data.frame(genes2, col.names = "GeneSymbol")
        GeneSymbol$Rank <- rownames(GeneSymbol)
        GeneSymbol <- GeneSymbol[,c("Rank","GeneSymbol")]
        DT::datatable(GeneSymbol, options = list(paging = F), rownames = F)
      }
    }
    else if (input$tables == 3){
      if (length(input$tab2table_rows_selected) > 0){
        res <- datasetInput()
        gsea.df <- as.data.frame(res@result)
        GS <- as.character(GeneSet2[input$tab2table_rows_selected,1])
        ## Subset core enriched genes
        genes1 <- as.matrix(gsea.df[which(gsea.df$Description==GS),"core_enrichment"])
        genes2 <- strsplit(genes1,"/")
        GeneSymbol <- as.data.frame(genes2, col.names = "GeneSymbol")
        GeneSymbol$Rank <- rownames(GeneSymbol)
        GeneSymbol <- GeneSymbol[,c("Rank","GeneSymbol")]
        DT::datatable(GeneSymbol, options = list(paging = F), rownames = F)
      }
    }
    else if (input$tables == 5){
      if (length(input$GStable.u_rows_selected) > 0){
        res <- datasetInput()
        gsea.df <- as.data.frame(res@result)
        GS <- as.character(user_gs_mirror()[input$GStable.u_rows_selected,1])
        ## Subset core enriched genes
        genes1 <- as.matrix(gsea.df[which(gsea.df$Description==GS),"core_enrichment"])
        genes2 <- strsplit(genes1,"/")
        GeneSymbol <- as.data.frame(genes2, col.names = "GeneSymbol")
        GeneSymbol$Rank <- rownames(GeneSymbol)
        GeneSymbol <- GeneSymbol[,c("Rank","GeneSymbol")]
        DT::datatable(GeneSymbol, options = list(paging = F), rownames = F)
      }
    }
  })
  
  ##render pre-loaded enriched signatures table
  #output$enrich_sig_table <- DT::renderDataTable({
  #    gsea.df <- as_tibble(get(ES_Tab_List[which(ES_Tab_List == paste("ES_table",match(input$SigTableChoice, SigNames),sep = ""))]))
  #    DT::datatable(gsea.df,
  #                  extensions = c("KeyTable", "FixedHeader"),
  #                  caption = "Enriched Signatures",
  #                  options = list(keys = T, searchHighlight = T, pageLength = 10, lengthMenu = c("10", "25", "50", "100"),scrollX = T)) %>%
  #        formatRound(columns = c(2:10), digits = 2)
  #})
  
  ##render user generated enriched signature table based off other gmt
  #output$enrich_sig_table_gen <- DT::renderDataTable({
  #    expr = updated_expr()
  #    meta = meta()
  #    A = as.matrix(expr)
  #    
  #    if (input$tables == 3) {
  #        #groupA <- meta[,1][meta[,2] == input$comparisonA]
  #        #groupB <- meta[,1][meta[,2] == input$comparisonB]
  #        groupA <- meta[which(meta[,2] == input$comparisonA),1]
  #        groupB <- meta[which(meta[,2] == input$comparisonB),1]	    
  #        ##----Signal-to-Noise Calculation----##
  #        A <- A + 0.00000001
  #        P = as.matrix(as.numeric(colnames(A) %in% groupA))
  #        n1 <- sum(P[,1])
  #        M1 <- A %*% P
  #        M1 <- M1/n1
  #        A2 <- A*A
  #        S1 <- A2 %*% P
  #        S1 <- S1/n1 - M1*M1 
  #        S1 <- sqrt(abs((n1/(n1-1)) * S1))
  #        P = as.matrix(as.numeric(colnames(A) %in% groupB))
  #        n2 <- sum(P[,1])
  #        M2 <- A %*% P
  #        M2 <- M2/n2
  #        A2 <- A*A
  #        S2 <- A2 %*% P
  #        S2 <- S2/n2 - M2*M2
  #        S2 <- sqrt(abs((n2/(n2-1)) * S2))
  #        rm(A2)
  #        # small sigma "fix" as used in GeneCluster
  #        S2 <- ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
  #        S2 <- ifelse(S2 == 0, 0.2, S2)
  #        S1 <- ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
  #        S1 <- ifelse(S1 == 0, 0.2, S1)
  #        M1 <- M1 - M2
  #        rm(M2)
  #        S1 <- S1 + S2
  #        rm(S2)
  #        s2n.matrix <- M1/S1
  #        ##----Reformatting----##
  #        s2n.df <- as.data.frame(s2n.matrix)
  #        s2n.df$GeneID <- rownames(s2n.df)
  #        rownames(s2n.df) <- NULL
  #        data <- dplyr::select(s2n.df, GeneID, V1)
  #        data.gsea <- data$V1
  #        names(data.gsea) <- as.character(data$GeneID)
  #        s2n.matrix.s <- sort(data.gsea, decreasing = T)
  #        ##----GSEA----##
  #        gmt.i <- tab2
  #        gsea.res <- GSEA(s2n.matrix.s, TERM2GENE = gmt.i, verbose = F, pvalueCutoff = input$userPval)
  #        gsea.df <- as_tibble(gsea.res@result)
  #        ## displaying the GSEA results as interactive data table
  #        DT::datatable(gsea.df,
  #                      extensions = c("KeyTable", "FixedHeader"),
  #                      caption = "Enriched Signatures",
  #                      options = list(keys = T, searchHighlight = T, pageLength = 10, lengthMenu = c("10", "25", "50", "100"), scrollX = T)) %>%
  #            formatRound(columns = c(2:10), digits = 2)
  #    }
  #    else if (input$tables == 5) {
  #        #groupA <- meta[,1][meta[,2] == input$comparisonA]
  #        #groupB <- meta[,1][meta[,2] == input$comparisonB]
  #        groupA <- meta[which(meta[,2] == input$comparisonA),1]
  #        groupB <- meta[which(meta[,2] == input$comparisonB),1]	    
  #        ##----Signal-to-Noise Calculation----##
  #        A <- A + 0.00000001
  #        P = as.matrix(as.numeric(colnames(A) %in% groupA))
  #        n1 <- sum(P[,1])
  #        M1 <- A %*% P
  #        M1 <- M1/n1
  #        A2 <- A*A
  #        S1 <- A2 %*% P
  #        S1 <- S1/n1 - M1*M1 
  #        S1 <- sqrt(abs((n1/(n1-1)) * S1))
  #        P = as.matrix(as.numeric(colnames(A) %in% groupB))
  #        n2 <- sum(P[,1])
  #        M2 <- A %*% P
  #        M2 <- M2/n2
  #        A2 <- A*A
  #        S2 <- A2 %*% P
  #        S2 <- S2/n2 - M2*M2
  #        S2 <- sqrt(abs((n2/(n2-1)) * S2))
  #        rm(A2)
  #        # small sigma "fix" as used in GeneCluster
  #        S2 <- ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
  #        S2 <- ifelse(S2 == 0, 0.2, S2)
  #        S1 <- ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
  #        S1 <- ifelse(S1 == 0, 0.2, S1)
  #        M1 <- M1 - M2
  #        rm(M2)
  #        S1 <- S1 + S2
  #        rm(S2)
  #        s2n.matrix <- M1/S1
  #        ##----Reformatting----##
  #        s2n.df <- as.data.frame(s2n.matrix)
  #        s2n.df$GeneID <- rownames(s2n.df)
  #        rownames(s2n.df) <- NULL
  #        data <- dplyr::select(s2n.df, GeneID, V1)
  #        data.gsea <- data$V1
  #        names(data.gsea) <- as.character(data$GeneID)
  #        s2n.matrix.s <- sort(data.gsea, decreasing = T)
  #        ##----GSEA----##
  #        gmt.i <- GStable.ubg()
  #        gsea.res <- GSEA(s2n.matrix.s, TERM2GENE = gmt.i, verbose = F, pvalueCutoff = input$userPval)
  #        gsea.df <- as_tibble(gsea.res@result)
  #        ## displaying the GSEA results as interactive data table
  #        DT::datatable(gsea.df,
  #                      extensions = c("KeyTable", "FixedHeader"),
  #                      caption = "Enriched Signatures",
  #                      options = list(keys = T, searchHighlight = T, pageLength = 10, lengthMenu = c("10", "25", "50", "100"), scrollX = T)) %>%
  #            formatRound(columns = c(2:10), digits = 2)
  #    }
  #})
  
  #render variable genes list in sidebar from heatmap
  output$MostVariableGenesList <- DT::renderDataTable({
    expr = updated_expr()
    meta = meta()
    
    top_probes <- input$NumFeatures
    col_labels <- colnames(expr)
    isexpr <- rowSums(expr[,col_labels] > 1) >= length(col_labels)
    exp <- expr[isexpr,]
    mad <- NULL
    var <- NULL
    cv <- NULL
    var_type <- input$VarianceMeasure
    if (var_type == "MAD"){
      mad <- apply(log2(exp + 1), 1, mad)
      mad <- sort(mad, decreasing = T)
      mad <- head(mad, n = (top_probes +1))
      out <- cbind(names(mad), mad[names(mad)], exp[names(mad),])
      colnames(out) <- c("Gene", "MAD", colnames(exp))
      dataset <- exp[names(mad),]
      variable_gene_list <- names(mad)
    }
    if (var_type == "VAR"){
      var <- apply(log2(exp + 1), 1, var)
      var <- sort(var, decreasing = T)
      var <- head(var, n = (top_probes +1))
      out <- cbind(names(var), var[names(var)], exp[names(var),])
      colnames(out) <- c("Gene", "VAR", colnames(exp))
      dataset <- exp[names(var),]
      variable_gene_list <- names(var)
    }
    if (var_type == "CV"){
      cv <- apply(log2(exp + 1), 1, cv)
      cv <- sort(cv, decreasing = T)
      cv <- head(cv, n = (top_probes +1))
      out <- cbind(names(cv), cv[names(cv)], exp[names(cv),])
      colnames(out) <- c("Gene", "VAR", colnames(exp))
      dataset <- exp[names(cv),]
      variable_gene_list <- names(cv)
    }
    variable_gene_list <- as.data.frame(variable_gene_list)
    colnames(variable_gene_list)[1] <- "Genes"
    variable_gene_list$Rank <- rownames(variable_gene_list)
    variable_gene_list <- variable_gene_list %>%
      select(Rank, Genes)
    DT::datatable(variable_gene_list, options = list(paging = F), rownames = F)
  })
  
  #render DEG table
  output$DEGtable1 <- DT::renderDataTable({
    expr = updated_expr()
    meta = meta()
    A <- meta[which(meta[,2] == input$comparisonA2.DEG),1]
    B <- meta[which(meta[,2] == input$comparisonB2.DEG),1]
    
    #A <- meta[,1][meta[,2] == input$comparisonA2.DEG]
    #B <- meta[,1][meta[,2] == input$comparisonB2.DEG]
    mat <- expr[,c(A,B)]
    mat <- log2(mat + 1.0)
    groupAOther <- factor(c(rep("A", length(A)), rep("B", length(B))))
    designA <- model.matrix(~0 + groupAOther)
    fit <- lmFit(mat, design = designA)
    contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    options(digits = 4)
    top1 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
    DT::datatable(top1, options = list(lengthMenu = c(50,100,1000, 5000, 10000), pageLength = 100, scrollX = TRUE),
                  selection=list(mode = "multiple"))
  })
  
  #render up regulated pathway enrichment data table
  output$UpRegPathwayTable1 <- DT::renderDataTable({
    expr = updated_expr()
    meta = meta()
    adjp <- input$pathpval
    FC <- input$pathFC
    A <- meta[which(meta[,2] == input$comparisonA2.path),1]
    B <- meta[which(meta[,2] == input$comparisonB2.path),1]
    #A <- meta[,1][meta[,2] == input$comparisonA2.path]
    #B <- meta[,1][meta[,2] == input$comparisonB2.path]
    mat <- expr[,c(A,B)]
    mat <- log2(mat + 1.0)
    groupAOther <- factor(c(rep("A", length(A)), rep("B", length(B))))
    designA <- model.matrix(~0 + groupAOther)
    fit <- lmFit(mat, design = designA)
    contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    options(digits = 4)
    top1 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
    genes <- rownames(top1)[which(top1$adj.P.Val < adjp & top1$logFC > FC)]
    dbs <- listEnrichrDbs() 
    enrichRLive <- TRUE 
    if (is.null(dbs)) { 
      enrichRLive <- FALSE 
    }
    dbs <- input$SelectedPathway
    enriched <- enrichr(genes, dbs) # Plot top 20 GO-BP results ordered by P-value
    DT::datatable(enriched[[1]], options = list(lengthMenu = c(10,20,50,100,1000), pageLength = 10, scrollX = TRUE),
                  selection=list(mode = "multiple"))
  })
  
  #render down regulated pathway enrichment data table
  output$DnRegPathwayTable1 <- DT::renderDataTable({
    expr = updated_expr()
    meta = meta()
    adjp <- input$pathpval
    FC <- input$pathFC
    A <- meta[which(meta[,2] == input$comparisonA2.path),1]
    B <- meta[which(meta[,2] == input$comparisonB2.path),1]
    
    #A <- meta[,1][meta[,2] == input$comparisonA2.path]
    #B <- meta[,1][meta[,2] == input$comparisonB2.path]
    mat <- expr[,c(A,B)]
    mat <- log2(mat + 1.0)
    groupAOther <- factor(c(rep("A", length(A)), rep("B", length(B))))
    designA <- model.matrix(~0 + groupAOther)
    fit <- lmFit(mat, design = designA)
    contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    options(digits = 4)
    top1 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
    genes <- rownames(top1)[which(top1$adj.P.Val < adjp & top1$logFC < -FC)]
    dbs <- listEnrichrDbs() 
    enrichRLive <- TRUE 
    if (is.null(dbs)) { 
      enrichRLive <- FALSE 
    }
    dbs <- input$SelectedPathway
    enriched <- enrichr(genes, dbs) # Plot top 20 GO-BP results ordered by P-value
    DT::datatable(enriched[[1]], options = list(lengthMenu = c(10,20,50,100,1000), pageLength = 10, scrollX = TRUE),
                  selection=list(mode = "multiple"))
  })
  
  #render gene list table for boxplot selection
  output$GeneListTable <- DT::renderDataTable({
    geneList <- geneList()
    DT::datatable(geneList,
                  selection = 'single',
                  options = list(keys = TRUE,
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = c("10", "25", "50", "100")))
  })
  
  #render gene list table for boxplot selection
  output$GeneListTable2 <- DT::renderDataTable({
    geneList <- geneList()
    DT::datatable(geneList,
                  selection = 'single',
                  options = list(keys = TRUE,
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = c("10", "25", "50", "100")))
  })
  
  #render gene scatter plot data table
  output$geneScatterTable <- DT::renderDataTable({
    meta = meta();
    expr = updated_expr();
    
    #log if user designates
    if (input$logask == TRUE) {
      expr <- log2(expr + 1)
    }
    #transpose
    expr_t <- as.data.frame(t(expr))
    #reorder rowname to match meta for merging
    samporder <- meta[,1]
    expr_t2 <- as.data.frame(expr_t[samporder,])
    #add type
    expr_t3 <- expr_t2 %>% 
      mutate(type = case_when(
        rownames(expr_t2) == meta[,1] ~ meta[,2],
      ))
    expr_t3 <- expr_t3 %>%
      relocate(type)
    #user gene input
    gene1.u <- input$scatterG1
    gene2.u <- input$scatterG2
    #get columns and info based off user input
    gene1 <- expr_t3[,gene1.u]
    gene2 <- expr_t3[,gene2.u]
    if (input$logask == TRUE) {
      gene1.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene1.u)], " Expression (log2 +1)", sep = "")
      gene2.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene2.u)], " Expression (log2 +1)", sep = "")
    }
    else if (input$logask == FALSE) {
      gene2.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene2.u)], " Expression", sep = "")
      gene1.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene1.u)], " Expression", sep = "")
    }
    #make table
    Sample <- rownames(expr_t3)
    Type <- expr_t3[,'type']
    gene1col <- expr_t3[,gene1.u]
    gene2col <- expr_t3[,gene2.u]
    scatterTab <- data.frame(Sample, Type, gene1col, gene2col)
    colnames(scatterTab)[c(3,4)] <- c(gene1.n, gene2.n)
    #table output
    DT::datatable(scatterTab,
                  options = list(keys = TRUE,
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = c("10", "25", "50", "100")))
    
  })
  
  #render ssGSEA table
  output$ssGSEAtable <- DT::renderDataTable({
    meta = meta();
    expr = updated_expr();
    
    if (input$tables == 3) {
      if (length(input$tab2table_rows_selected) > 0){
        ssgsea <- ssGSEAfunc()
        ssgsea2 <- as.data.frame(t(ssgsea))
        samporder <- meta[,1]
        ssgsea3 <- as.data.frame(ssgsea2[samporder,])
        colnames(ssgsea3)[1] <- colnames(ssgsea2[1])
        rownames(ssgsea3) <- samporder
        ssgsea4 <- ssgsea3 %>% 
          mutate(type = case_when(
            rownames(ssgsea3) == meta[,1] ~ meta[,2],
          ))
        ssgsea4 <- ssgsea4 %>%
          relocate(type)
        #table output
        DT::datatable(ssgsea4,
                      options = list(keys = TRUE,
                                     searchHighlight = TRUE,
                                     pageLength = 10,
                                     lengthMenu = c("10", "25", "50", "100")))
      }
    }
    else if (input$tables == 1) {
      if (length(input$msigdbTable_rows_selected) > 0){
        ssgsea <- ssGSEAfunc()
        ssgsea2 <- as.data.frame(t(ssgsea))
        samporder <- meta[,1]
        ssgsea3 <- as.data.frame(ssgsea2[samporder,])
        colnames(ssgsea3)[1] <- colnames(ssgsea2[1])
        rownames(ssgsea3) <- samporder
        ssgsea4 <- ssgsea3 %>% 
          mutate(type = case_when(
            rownames(ssgsea3) == meta[,1] ~ meta[,2],
          ))
        ssgsea4 <- ssgsea4 %>%
          relocate(type)
        #table output
        DT::datatable(ssgsea4,
                      options = list(keys = TRUE,
                                     searchHighlight = TRUE,
                                     pageLength = 10,
                                     lengthMenu = c("10", "25", "50", "100")))
      }
    }
    else if (input$tables == 5) {
      if (length(input$GStable.u_rows_selected) > 0){
        ssgsea <- ssGSEAfunc()
        ssgsea2 <- as.data.frame(t(ssgsea))
        samporder <- meta[,1]
        ssgsea3 <- as.data.frame(ssgsea2[samporder,])
        colnames(ssgsea3)[1] <- colnames(ssgsea2[1])
        rownames(ssgsea3) <- samporder
        ssgsea4 <- ssgsea3 %>% 
          mutate(type = case_when(
            rownames(ssgsea3) == meta[,1] ~ meta[,2],
          ))
        ssgsea4 <- ssgsea4 %>%
          relocate(type)
        #table output
        DT::datatable(ssgsea4,
                      options = list(keys = TRUE,
                                     searchHighlight = TRUE,
                                     pageLength = 10,
                                     lengthMenu = c("10", "25", "50", "100")))
      }
    }
  })
  
  
  ####----Plots----####
  
  
  #render GSEA plot
  output$enrichplot0 <- renderPlot({
    meta = meta();
    expr = updated_expr();
    A = as.matrix(expr) 
    if (input$tables == 1){
      if (length(input$msigdbTable_rows_selected) > 0){
        res <- datasetInput()
        gsea.df <- as_tibble(res@result)
        geneset <- as.character(msigdb.gsea2[input$msigdbTable_rows_selected,3])
        title <- geneset
        gseaplot2(res,
                  geneset,
                  title,
                  pvalue_table = F)
      }
    }
    else if (input$tables == 3){
      if (length(input$tab2table_rows_selected) > 0){
        res <- datasetInput()
        gsea.df <- as_tibble(res@result)
        geneset <- as.character(GeneSet2[input$tab2table_rows_selected,1])
        title <- geneset
        gseaplot2(res,
                  geneset,
                  title,
                  pvalue_table = F)
      }
    }
    else if (input$tables == 5){
      if (length(input$GStable.u_rows_selected) > 0){
        res <- datasetInput()
        gsea.df <- as_tibble(res@result)
        geneset <- as.character(user_gs_mirror()[input$GStable.u_rows_selected,1])
        title <- geneset
        gseaplot2(res,
                  geneset,
                  title,
                  pvalue_table = F)
      }
    }
    
    #res <- datasetInput()
    #gsea.df <- as_tibble(res@result)
    #geneset <- as.matrix(msigdb.gsea2[input$msigdbTable_rows_selected,3])
    #geneset = "HALLMARK_ADIPOGENESIS"
    #title <- geneset
    #gseaplot2(res,
    #    geneset,
    #    title,
    #    pvalue_table = F)
  })
  
  #render heatmap
  output$heatmap0 <- renderPlot({
    expr = updated_expr()
    meta = meta()
    A = as.matrix(expr)
    
    if (input$tables == 1) {
      if (length(input$msigdbTable_rows_selected) > 0){
        groupA <- meta[which(meta[,2] == input$comparisonA),1]
        groupB <- meta[which(meta[,2] == input$comparisonB),1]		    
        #groupA <- meta[,1][meta[,2] == input$comparisonA]
        #groupB <- meta[,1][meta[,2] == input$comparisonB]
        res <- datasetInput()
        #gsea.df <- as.data.frame(res@result)
        gsea.df <- as_tibble(res@result)
        GS <- msigdb.gsea2[input$msigdbTable_rows_selected,3]
        genes1 <- as.matrix(gsea.df[which(gsea.df$Description==GS),"core_enrichment"])
        genes2 <- strsplit(genes1,"/")
        genes3 <- as.data.frame(genes2, col.names = "genes")
        gene_symbol <- genes3$genes
        ## convert expression matrix to numeric
        class(A) <- "numeric"
        ## Transforming data
        A <- A[,c(groupA,groupB)]
        exp.mat1 = log2(A + 1) # log
        exp.mat2 = apply(exp.mat1, 1, scale); # z score
        exp.mat3 = apply(exp.mat2, 1, rev); # transpose
        colnames(exp.mat3) = colnames(A) # set the column name
        exp.mat4 = exp.mat3[rowSums(is.na(exp.mat3[,])) == 0,] # remove NaN rows
        exp.mat5 = exp.mat4[rownames(exp.mat4) %in% gene_symbol,] # grab gene symbols
        # reassign data
        dataset <- exp.mat5
        ## generate color for pheatmap
        if (abs(min(dataset)) > abs(max(dataset))) {
          dataset[dataset < -abs(max(dataset))] = -abs(max(dataset))
        } else {
          dataset[dataset > abs(min(dataset))] = abs(min(dataset))
        }
        meta2 <- meta[which(meta[,1] %in% c(groupA,groupB)),]
        meta2 <- meta2[order(meta2[,2]),]
        type <- meta2[,2]
        meta3 <- as.data.frame(type)
        rownames(meta3) <- meta2[,1]
        zscore_range = 10;
        minimum = -zscore_range;
        maximum = zscore_range;
        bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
        hmcols<- colorRampPalette(c("dark blue","blue","white","red", "dark red"))(length(bk)-1)
        pheatmap(dataset,            #data
                 cluster_cols = F,    #cluster columns - NO
                 cluster_row = F,     #cluster rows - YES
                 fontsize_col = input$heatmapFont1.c,   #column fontsize
                 fontsize_row = input$heatmapFont1.r,
                 show_rownames = T,  
                 show_colnames = T,
                 color=hmcols,
                 annotation_col = meta3)
      }
    }
    else if (input$tables == 3) {
      if (length(input$tab2table_rows_selected) > 0){
        groupA <- meta[which(meta[,2] == input$comparisonA),1]
        groupB <- meta[which(meta[,2] == input$comparisonB),1]
        #groupA <- meta[,1][meta[,2] == input$comparisonA]
        #groupB <- meta[,1][meta[,2] == input$comparisonB]
        
        res <- datasetInput()
        gsea.df <- as.data.frame(res@result)
        GS <- as.character(GeneSet2[input$tab2table_rows_selected,1])
        genes1 <- as.matrix(gsea.df[which(gsea.df$Description==GS),"core_enrichment"])
        genes2 <- strsplit(genes1,"/")
        genes3 <- as.data.frame(genes2, col.names = "genes")
        gene_symbol <- genes3$genes
        ## convert expression matrix to numeric
        class(A) <- "numeric"
        ## Transforming data
        A <- A[,c(groupA,groupB)]
        exp.mat1 = log2(A + 1) # log
        exp.mat2 = apply(exp.mat1, 1, scale); # z score
        exp.mat3 = apply(exp.mat2, 1, rev); # transpose
        colnames(exp.mat3) = colnames(A) # set the column name
        exp.mat4 = exp.mat3[rowSums(is.na(exp.mat3[,])) == 0,] # remove NaN rows
        exp.mat5 = exp.mat4[rownames(exp.mat4) %in% gene_symbol,] # grab gene symbols
        # reassign data
        dataset <- exp.mat5
        ## generate color for pheatmap
        if (abs(min(dataset)) > abs(max(dataset))) {
          dataset[dataset < -abs(max(dataset))] = -abs(max(dataset))
        } else {
          dataset[dataset > abs(min(dataset))] = abs(min(dataset))
        }
        meta2 <- meta[which(meta[,1] %in% c(groupA,groupB)),]
        meta2 <- meta2[order(meta2[,2]),]
        type <- meta2[,2]
        meta3 <- as.data.frame(type)
        rownames(meta3) <- meta2[,1]
        zscore_range = 10;
        minimum = -zscore_range;
        maximum = zscore_range;
        bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
        hmcols<- colorRampPalette(c("dark blue","blue","white","red", "dark red"))(length(bk)-1)
        pheatmap(dataset,            #data
                 cluster_cols = F,    #cluster columns - NO
                 cluster_row = F,     #cluster rows - YES
                 fontsize_col = input$heatmapFont1.c,   #column fontsize
                 fontsize_row = input$heatmapFont1.r,
                 show_rownames = T,  
                 show_colnames = T,
                 color=hmcols,
                 annotation_col = meta3)
      }
    }
    else if (input$tables == 5) {
      if (length(input$GStable.u_rows_selected) > 0){
        groupA <- meta[which(meta[,2] == input$comparisonA),1]
        groupB <- meta[which(meta[,2] == input$comparisonB),1]
        #groupA <- meta[,1][meta[,2] == input$comparisonA]
        #groupB <- meta[,1][meta[,2] == input$comparisonB]
        
        res <- datasetInput()
        gsea.df <- as.data.frame(res@result)
        GS <- as.character(user_gs_mirror()[input$GStable.u_rows_selected,1])
        genes1 <- as.matrix(gsea.df[which(gsea.df$Description==GS),"core_enrichment"])
        genes2 <- strsplit(genes1,"/")
        genes3 <- as.data.frame(genes2, col.names = "genes")
        gene_symbol <- genes3$genes
        ## convert expression matrix to numeric
        class(A) <- "numeric"
        ## Transforming data
        A <- A[,c(groupA,groupB)]
        exp.mat1 = log2(A + 1) # log
        exp.mat2 = apply(exp.mat1, 1, scale); # z score
        exp.mat3 = apply(exp.mat2, 1, rev); # transpose
        colnames(exp.mat3) = colnames(A) # set the column name
        exp.mat4 = exp.mat3[rowSums(is.na(exp.mat3[,])) == 0,] # remove NaN rows
        exp.mat5 = exp.mat4[rownames(exp.mat4) %in% gene_symbol,] # grab gene symbols
        # reassign data
        dataset <- exp.mat5
        ## generate color for pheatmap
        if (abs(min(dataset)) > abs(max(dataset))) {
          dataset[dataset < -abs(max(dataset))] = -abs(max(dataset))
        } else {
          dataset[dataset > abs(min(dataset))] = abs(min(dataset))
        }
        meta2 <- meta[which(meta[,1] %in% c(groupA,groupB)),]
        meta2 <- meta2[order(meta2[,2]),]
        type <- meta2[,2]
        meta3 <- as.data.frame(type)
        rownames(meta3) <- meta2[,1]
        zscore_range = 10;
        minimum = -zscore_range;
        maximum = zscore_range;
        bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
        hmcols<- colorRampPalette(c("dark blue","blue","white","red", "dark red"))(length(bk)-1)
        pheatmap(dataset,            #data
                 cluster_cols = F,    #cluster columns - NO
                 cluster_row = F,     #cluster rows - YES
                 fontsize_col = input$heatmapFont1.c,   #column fontsize
                 fontsize_row = input$heatmapFont1.r,
                 show_rownames = T,  
                 show_colnames = T,
                 color=hmcols,
                 annotation_col = meta3)
      }
    }
  })
  
  #render custom DEG heatmap
  output$heatmap2 <- renderPlot({
    expr = updated_expr()
    meta = meta()
    
    if (length(input$heatmapGeneSelec) >= 2 || length(input$userheatgenes) >= 1) {
      genelist.uih <- NULL
      genelist.ush <- NULL
      genelist.uih2 <- NULL
      genelist.ush <- input$heatmapGeneSelec
      genelist.uih <- unlist(strsplit(input$userheatgenes, " "))
      genelist.uih2 <- unlist(strsplit(input$userheatgenes, "\t"))
      heatgenes <- c(genelist.ush,genelist.uih,genelist.uih2)
      usersamps <- input$userheatsamp2
      exp <- expr[heatgenes,usersamps]
      meta <- meta[which(meta[,1] %in% usersamps),]
      dataset <- exp
      dataset <- log2(dataset + 1)
      zdataset <- apply(dataset, 1, scale)
      zdataset <- apply(zdataset, 1, rev)
      colnames(zdataset) <- names(dataset)
      dataset <- as.matrix(zdataset)
      dataset[is.na(dataset)] <- 0
      
      dataset = dataset[apply(dataset[,-1], 1, function(x) !all(x==0)),]
      minimum = -5;
      maximum = 5;
      if (abs(min(dataset)) > abs(max(dataset))) {
        dataset[dataset < -abs(max(dataset))] = -abs(max(dataset))
      } else {
        dataset[dataset > abs(min(dataset))] = abs(min(dataset))
      }
      meta <- meta[order(meta[,2]),]
      type <- meta[,2]
      meta2 <- as.data.frame(type)
      rownames(meta2) <- meta[,1]
      bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
      hmcols<- colorRampPalette(c("dark blue","blue","white","red", "dark red"))(length(bk)-1)
      pheatmap(dataset,
               cluster_col = F,
               cluster_row = T,
               fontsize_row = input$heatmapFont3.r,
               fontsize_col = input$heatmapFont3.c,
               show_rownames = T ,
               show_colnames = T,
               annotation_col = meta2,
               clustering_method = input$ClusteringMethod2,
               color=hmcols,
               border_color = NA)
    }
    
  })
  
  #render DEG heatmap
  output$heatmap1 <- renderPlot({
    expr = updated_expr()
    meta = meta()
    top_probes <- input$NumFeatures
    exp <- expr
    mad <- NULL
    var <- NULL
    cv <- NULL
    var_type <- input$VarianceMeasure
    if (var_type == "MAD"){
      mad <- apply(log2(exp + 1), 1, mad)
      mad <- sort(mad, decreasing = T)
      mad <- head(mad, n = (top_probes +1))
      out <- cbind(names(mad), mad[names(mad)], exp[names(mad),])
      colnames(out) <- c("Gene", "MAD", colnames(exp))
      dataset <- exp[names(mad),]
      variable_gene_list <- names(mad)
    }
    else if (var_type == "VAR"){
      var <- apply(log2(exp + 1), 1, var)
      var <- sort(var, decreasing = T)
      var <- head(var, n = (top_probes +1))
      out <- cbind(names(var), var[names(var)], exp[names(var),])
      colnames(out) <- c("Gene", "VAR", colnames(exp))
      dataset <- exp[names(var),]
      variable_gene_list <- names(var)
    }
    else if (var_type == "CV"){
      cv <- apply(log2(exp + 1), 1, cv)
      cv <- sort(cv, decreasing = T)
      cv <- head(cv, n = (top_probes +1))
      out <- cbind(names(cv), cv[names(cv)], exp[names(cv),])
      colnames(out) <- c("Gene", "CV", colnames(exp))
      dataset <- exp[names(cv),]
      variable_gene_list <- names(cv)
    }
    dataset <- log2(dataset + 1)
    zdataset <- apply(dataset, 1, scale)
    zdataset <- apply(zdataset, 1, rev)
    colnames(zdataset) <- names(dataset)
    dataset <- as.matrix(zdataset)
    dataset[is.na(dataset)] <- 0
    
    dataset = dataset[apply(dataset[,-1], 1, function(x) !all(x==0)),]
    minimum = -5;
    maximum = 5;
    if (abs(min(dataset)) > abs(max(dataset))) {
      dataset[dataset < -abs(max(dataset))] = -abs(max(dataset))
    } else {
      dataset[dataset > abs(min(dataset))] = abs(min(dataset))
    }
    type <- meta[,2]
    meta2 <- as.data.frame(type)
    rownames(meta2) <- meta[,1]
    bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
    hmcols<- colorRampPalette(c("dark blue","blue","white","red", "dark red"))(length(bk)-1)
    pheatmap(dataset,
             cluster_col = T,
             cluster_row = T,
             fontsize_row = input$heatmapFont2.r,
             fontsize_col = input$heatmapFont2.c,
             show_rownames = T ,
             annotation_col = meta2,
             show_colnames = T,
             clustering_method = input$ClusteringMethod,
             color=hmcols,
             border_color = NA)
    
  })
  
  ##render RNAseq heatmap
  #output$heatmap1 <- renderPlot({
  #    expr = updated_expr()
  #    meta = meta()
  #    if (input$customs == 444){
  #        if (length(input$heatmapGeneSelec) >= 2 || length(input$userheatgenes) >= 1) {
  #            genelist.uih <- NULL
  #            genelist.ush <- NULL
  #            genelist.uih2 <- NULL
  #            genelist.ush <- input$heatmapGeneSelec
  #            genelist.uih <- unlist(strsplit(input$userheatgenes, " "))
  #            genelist.uih2 <- unlist(strsplit(input$userheatgenes, "\t"))
  #            heatgenes <- c(genelist.ush,genelist.uih,genelist.uih2)
  #            usersamps <- input$userheatsamp2
  #            exp <- expr[heatgenes,usersamps]
  #            meta <- meta[which(meta[,1] %in% usersamps),]
  #            dataset <- exp
  #            dataset <- log2(dataset + 1)
  #            zdataset <- apply(dataset, 1, scale)
  #            zdataset <- apply(zdataset, 1, rev)
  #            colnames(zdataset) <- names(dataset)
  #            dataset <- as.matrix(zdataset)
  #            dataset[is.na(dataset)] <- 0
  #            
  #            dataset = dataset[apply(dataset[,-1], 1, function(x) !all(x==0)),]
  #            minimum = -5;
  #            maximum = 5;
  #            if (abs(min(dataset)) > abs(max(dataset))) {
  #                dataset[dataset < -abs(max(dataset))] = -abs(max(dataset))
  #            } else {
  #                dataset[dataset > abs(min(dataset))] = abs(min(dataset))
  #            }
  #            meta <- meta[order(meta[,2]),]
  #            type <- meta[,2]
  #            meta2 <- as.data.frame(type)
  #            rownames(meta2) <- meta[,1]
  #            bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
  #            hmcols<- colorRampPalette(c("dark blue","blue","white","red", "dark red"))(length(bk)-1)
  #            pheatmap(dataset,
  #                     cluster_col = F,
  #                     cluster_row = T,
  #                     fontsize_row = input$heatmapFont3.r,
  #                     fontsize_col = input$heatmapFont3.c,
  #                     show_rownames = T ,
  #                     show_colnames = T,
  #                     annotation_col = meta2,
  #                     clustering_method = input$ClusteringMethod2,
  #                     color=hmcols,
  #                     border_color = NA)
  #        }
  #    }
  #    else {
  #        top_probes <- input$NumFeatures
  #        exp <- expr
  #        mad <- NULL
  #        var <- NULL
  #        cv <- NULL
  #        var_type <- input$VarianceMeasure
  #        if (var_type == "MAD"){
  #            mad <- apply(log2(exp + 1), 1, mad)
  #            mad <- sort(mad, decreasing = T)
  #            mad <- head(mad, n = (top_probes +1))
  #            out <- cbind(names(mad), mad[names(mad)], exp[names(mad),])
  #            colnames(out) <- c("Gene", "MAD", colnames(exp))
  #            dataset <- exp[names(mad),]
  #            variable_gene_list <- names(mad)
  #        }
  #        else if (var_type == "VAR"){
  #            var <- apply(log2(exp + 1), 1, var)
  #            var <- sort(var, decreasing = T)
  #            var <- head(var, n = (top_probes +1))
  #            out <- cbind(names(var), var[names(var)], exp[names(var),])
  #            colnames(out) <- c("Gene", "VAR", colnames(exp))
  #            dataset <- exp[names(var),]
  #            variable_gene_list <- names(var)
  #        }
  #        else if (var_type == "CV"){
  #            cv <- apply(log2(exp + 1), 1, cv)
  #            cv <- sort(cv, decreasing = T)
  #            cv <- head(cv, n = (top_probes +1))
  #            out <- cbind(names(cv), cv[names(cv)], exp[names(cv),])
  #            colnames(out) <- c("Gene", "VAR", colnames(exp))
  #            dataset <- exp[names(cv),]
  #            variable_gene_list <- names(cv)
  #        }
  #        dataset <- log2(dataset + 1)
  #        zdataset <- apply(dataset, 1, scale)
  #        zdataset <- apply(zdataset, 1, rev)
  #        colnames(zdataset) <- names(dataset)
  #        dataset <- as.matrix(zdataset)
  #        dataset[is.na(dataset)] <- 0
  #        
  #        dataset = dataset[apply(dataset[,-1], 1, function(x) !all(x==0)),]
  #        minimum = -5;
  #        maximum = 5;
  #        if (abs(min(dataset)) > abs(max(dataset))) {
  #            dataset[dataset < -abs(max(dataset))] = -abs(max(dataset))
  #        } else {
  #            dataset[dataset > abs(min(dataset))] = abs(min(dataset))
  #        }
  #        type <- meta[,2]
  #        meta2 <- as.data.frame(type)
  #        rownames(meta2) <- meta[,1]
  #        bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
  #        hmcols<- colorRampPalette(c("dark blue","blue","white","red", "dark red"))(length(bk)-1)
  #        pheatmap(dataset,
  #                 cluster_col = T,
  #                 cluster_row = T,
  #                 fontsize_row = input$heatmapFont2.r,
  #                 fontsize_col = input$heatmapFont2.c,
  #                 show_rownames = T ,
  #                 annotation_col = meta2,
  #                 show_colnames = T,
  #                 clustering_method = input$ClusteringMethod,
  #                 color=hmcols,
  #                 border_color = NA)
  #    }
  #    
  #})
  
  #render MA plot
  output$MAPlot1 <- renderPlot({
    top2 <- topgenereact()
    #add color categories based on FC and pval
    top2['threshold'] <- "none"
    top2[which(top2$logFC > abs(input$fc_cutoff)), "threshold"] <- "up"
    top2[which(top2$logFC < -abs(input$fc_cutoff)), "threshold"] <- "down"
    upRed <- "lightcoral"
    dnBlue <- "cadetblue3"
    mdGray <- "gray70"
    top2u <- top2[order(top2[,1], decreasing = TRUE),]
    top2d <- top2[order(top2[,1]),]
    top_hits_up <- top2u[head(which(top2u$logFC > abs(input$fc_cutoff)), n = input$top_x),]
    top_hits_dn <- top2d[head(which(top2d$logFC < -abs(input$fc_cutoff)), n = input$top_x),]
    #select genes to label based on user selection
    genesel.s <- NULL
    genesel.t <- NULL
    genesel.u <- NULL
    genesel.s <- unlist(strsplit(input$userGeneSelec2, " "))
    genesel.t <- unlist(strsplit(input$userGeneSelec2, "\t"))
    genesel.u <- input$userGeneSelec
    genesel.text <- c(genesel.s,genesel.t,genesel.u)
    top2_selec <- top2 %>%
      filter(GeneName %in% genesel.text)
    x <- ggplot(data = top2, aes(x = AveExpr, y = logFC)) +
      geom_point(size = 2, shape = 16) +
      theme_light(base_size = 16)
    #colors
    x <- x + aes(color = threshold) +
      scale_color_manual(values = c("up" = upRed,"down" = dnBlue, "none" = mdGray))
    x <- x + geom_hline(yintercept = c(-abs(input$fc_cutoff),abs(input$fc_cutoff)), linetype="dashed", color="gray20")
    x <- x + geom_text_repel(
      data =  top_hits_up,
      aes(label = rownames(top_hits_up)),
      color="gray20",
      size = 6,
      nudge_x = 0.2,
      nudge_y=0.2,
      box.padding = unit(0.9, "lines"),
      point.padding = unit(.3+4*0.1, "lines"),
      max.overlaps = 50)
    x <- x + geom_text_repel(
      data =  top_hits_dn,
      aes(label = rownames(top_hits_dn)),
      color="gray20",
      size = 6,
      nudge_x = 0.2,
      nudge_y=0.2,
      box.padding = unit(0.9, "lines"),
      point.padding = unit(.3+4*0.1, "lines"),
      max.overlaps = 50)
    x <- x + geom_text_repel(
      data =  top2_selec,
      aes(label = rownames(top2_selec)),
      color="gray20",
      size = 6,
      nudge_x = 0.2,
      nudge_y=0.2,
      box.padding = unit(0.9, "lines"),
      point.padding = unit(.3+4*0.1, "lines"),
      max.overlaps = 50)
    #coloring selected points
    x <- x + geom_point(data = top2_selec,
                        aes(x = AveExpr, y = logFC),
                        pch = 21,
                        color = "black",
                        size = 2)
    x <- x + geom_point(data = top_hits_dn,
                        aes(x = AveExpr, y = logFC),
                        pch = 21,
                        color = "black",
                        size = 2)
    x <- x + geom_point(data = top_hits_up,
                        aes(x = AveExpr, y = logFC),
                        pch = 21,
                        color = "black",
                        size = 2)
    x <- x + theme(legend.position="none")
    x <- x + labs(x = "Average Expression", y = "log2FC")
    x <- x + theme(axis.text = element_text(size=18))
    x <- x + theme(axis.title = element_text(size=24))
    x
  })
  
  #render boxplot
  output$boxplot1 <- renderPlot({
    geneList <- geneList()
    meta = meta()
    expr = updated_expr()
    if (length(input$GeneListTable_rows_selected) > 0){
      gene <- geneList[input$GeneListTable_rows_selected, 1]
      min <- min(log2(expr[gene,] + 1.0))
      max <- max(log2(expr[gene,] + 1.0))
      meta_temp <- meta
      rownames(meta_temp) <- meta[,1]
      meta_temp <- meta_temp %>%
        select(Group)
      data = merge(t(expr[gene,]), meta_temp, by=0)
      colnames(data) = c("SampleName", "GeneExpr", "Cluster")
      ggplot(data, aes(x=Cluster, y=log2(GeneExpr + 1.0))) +
        geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=1) +
        geom_dotplot(binaxis='y', stackdir='center', dotsize=input$boxplotDot) +
        stat_compare_means(label = "p.signif") + ylim(min * 0.9, max * 1.3) +
        theme_bw() +
        labs(title= paste(gene, "Expression (log2)")) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1),
              text = element_text(size = input$boxplot.font))
    }
  })
  
  #render boxplot
  output$boxplot3 <- renderPlot({
    geneList <- geneList()
    meta = meta()
    expr = updated_expr()
    if (length(input$GeneListTable2_rows_selected) > 0){
      gene <- geneList[input$GeneListTable2_rows_selected, 1]
      min <- min(log2(expr[gene,] + 1.0))
      max <- max(log2(expr[gene,] + 1.0))
      meta_temp <- meta
      rownames(meta_temp) <- meta[,1]
      meta_temp <- meta_temp %>%
        select(Group)
      data = merge(t(expr[gene,]), meta_temp, by=0)
      colnames(data) = c("SampleName", "GeneExpr", "Cluster")
      ggplot(data, aes(x=Cluster, y=log2(GeneExpr + 1.0))) +
        geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=1) +
        geom_dotplot(binaxis='y', stackdir='center', dotsize=input$boxplotDot2) +
        stat_compare_means(label = "p.signif") + ylim(min * 0.9, max * 1.3) +
        theme_bw() +
        labs(title= paste(gene, "Expression (log2)")) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1),
              text = element_text(size = input$boxplot.font2))
    }
  })
  
  #render up regulated pathway enrichment plot
  output$UpRegPathway1 <- renderPlot({
    meta = meta();
    expr = updated_expr();
    
    adjp <- input$pathpval
    FC <- input$pathFC
    A <- meta[which(meta[,2] == input$comparisonA2.path),1]
    B <- meta[which(meta[,2] == input$comparisonB2.path),1]
    
    #A <- meta[,1][meta[,2] == input$comparisonA2.path]
    #B <- meta[,1][meta[,2] == input$comparisonB2.path]
    mat <- expr[,c(A,B)]
    mat <- log2(mat + 1.0)
    groupAOther <- factor(c(rep("A", length(A)), rep("B", length(B))))
    designA <- model.matrix(~0 + groupAOther)
    fit <- lmFit(mat, design = designA)
    contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    options(digits = 4)
    top1 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
    genes <- rownames(top1)[which(top1$adj.P.Val < adjp & top1$logFC > FC)]
    dbs <- listEnrichrDbs() 
    enrichRLive <- TRUE 
    if (is.null(dbs)) { 
      enrichRLive <- FALSE 
    }
    dbs <- input$SelectedPathway
    enriched <- enrichr(genes, dbs) # Plot top 20 GO-BP results ordered by P-value 
    plotEnrich(enriched[[1]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value")
  })
  
  #render up regulated pathway enrichment plot
  output$DnRegPathway1 <- renderPlot({
    meta = meta();
    expr = updated_expr();
    
    
    adjp <- input$pathpval
    FC <- input$pathFC
    A <- meta[which(meta[,2] == input$comparisonA2.path),1]
    B <- meta[which(meta[,2] == input$comparisonB2.path),1]
    
    #A <- meta[,1][meta[,2] == input$comparisonA2.path]
    #B <- meta[,1][meta[,2] == input$comparisonB2.path]
    mat <- expr[,c(A,B)]
    mat <- log2(mat + 1.0)
    groupAOther <- factor(c(rep("A", length(A)), rep("B", length(B))))
    designA <- model.matrix(~0 + groupAOther)
    fit <- lmFit(mat, design = designA)
    contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    options(digits = 4)
    top1 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
    genes <- rownames(top1)[which(top1$adj.P.Val < adjp & top1$logFC < -FC)]
    dbs <- listEnrichrDbs() 
    enrichRLive <- TRUE 
    if (is.null(dbs)) { 
      enrichRLive <- FALSE 
    }
    dbs <- input$SelectedPathway
    enriched <- enrichr(genes, dbs) # Plot top 20 GO-BP results ordered by P-value 
    plotEnrich(enriched[[1]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value")
  })
  
  #render volcano plot
  output$Volcano3 <- renderPlot({
    meta = meta();
    expr = updated_expr();
    
    top2 <- topgenereact()
    #add color categories based on FC and pval
    top2['threshold'] <- "none"
    top2[which(top2$logFC > abs(input$fc_cutoff) & top2$P.Value < input$p_cutoff), "threshold"] <- "up"
    top2[which(top2$logFC < -abs(input$fc_cutoff) & top2$P.Value < input$p_cutoff), "threshold"] <- "down"
    upRed <- "lightcoral"
    dnBlue <- "cadetblue3"
    mdGray <- "gray70"
    #select number of top hits based on input
    top_hits_up <- top2[head(which(top2$logFC > abs(input$fc_cutoff) & top2$P.Value < input$p_cutoff), n = input$top_x),]
    top_hits_dn <- top2[head(which(top2$logFC < -abs(input$fc_cutoff) & top2$P.Value < input$p_cutoff), n = input$top_x),]
    #select genes to label based on user selection
    genesel.s <- NULL
    genesel.t <- NULL
    genesel.u <- NULL
    genesel.s <- unlist(strsplit(input$userGeneSelec2, " "))
    genesel.t <- unlist(strsplit(input$userGeneSelec2, "\t"))
    genesel.u <- input$userGeneSelec
    genesel.text <- c(genesel.s,genesel.t,genesel.u)
    top2_selec <- top2 %>%
      filter(GeneName %in% genesel.text)
    #create plot
    x <- ggplot(data = top2, aes(x = logFC, y = -log10(P.Value))) +
      geom_point(size = 2, shape = 16) +
      theme_light(base_size = 16)
    #colors
    x <- x + aes(color = threshold) +
      scale_color_manual(values = c("up" = upRed,"down" = dnBlue, "none" = mdGray))
    #FC and pval lines
    x <- x + geom_vline(xintercept = c(-abs(input$fc_cutoff),abs(input$fc_cutoff)), linetype="dashed", color="gray20")
    x <- x + geom_hline(yintercept = -log10(input$p_cutoff), linetype="dashed", color="gray20")
    #label top hits if needed
    if (input$top_x > 0) { 
      x <- x + geom_text_repel(
        data =  top_hits_up,
        aes(label = rownames(top_hits_up)),
        size = 6,
        color="gray20",
        nudge_x = 0.2,
        nudge_y=0.2,
        box.padding = unit(0.9, "lines"),
        point.padding = unit(.3+4*0.1, "lines"),
        max.overlaps = 50)
      x <- x + geom_text_repel(
        data =  top_hits_dn,
        aes(label = rownames(top_hits_dn)),
        size = 6,
        color="gray20",
        nudge_x = 0.2,
        nudge_y=0.2,
        box.padding = unit(0.9, "lines"),
        point.padding = unit(.3+4*0.1, "lines"),
        max.overlaps = 50)
    }
    x <- x + geom_text_repel(
      data =  top2_selec,
      aes(label = rownames(top2_selec)),
      size = 6,
      color="gray20",
      nudge_x = 0.2,
      nudge_y=0.2,
      box.padding = unit(0.9, "lines"),
      point.padding = unit(.3+4*0.1, "lines"),
      max.overlaps = 50)
    #coloring selected points
    x <- x + geom_point(data = top2_selec,
                        aes(x = logFC, y = -log10(P.Value)),
                        pch = 21,
                        color = "black",
                        size = 2)
    x <- x + geom_point(data = top_hits_dn,
                        aes(x = logFC, y = -log10(P.Value)),
                        pch = 21,
                        color = "black",
                        size = 2)
    x <- x + geom_point(data = top_hits_up,
                        aes(x = logFC, y = -log10(P.Value)),
                        pch = 21,
                        color = "black",
                        size = 2)
    #axis parameters
    x <- x + theme(legend.position="none")
    x <- x + theme(axis.text = element_text(size=18))
    x <- x + theme(axis.title = element_text(size=24))
    x <- x + labs(x = "log2FC", y = "-log10(P.Value)")
    x
  })
  
  #render ssGSEA boxplot
  output$boxplot2 <- renderPlot({
    meta = meta()
    expr = updated_expr()
    A = as.matrix(expr)
    
    #if tab 2 selected
    if (input$tables == 3) {
      if (length(input$tab2table_rows_selected) > 0){
        ssgsea <- ssGSEAfunc()
        ssgsea2 <- as.data.frame(t(ssgsea))
        samporder <- meta[,1]
        ssgsea3 <- as.data.frame(ssgsea2[samporder,])
        colnames(ssgsea3)[1] <- colnames(ssgsea2[1])
        rownames(ssgsea3) <- samporder
        ssgsea4 <- ssgsea3 %>% 
          mutate(type = case_when(
            rownames(ssgsea3) == meta[,1] ~ meta[,2],
          ))
        ssgsea4 <- ssgsea4 %>%
          relocate(type)
        if (input$boxplotcompare == "none") {
          ggplot(ssgsea4, aes(factor(type), ssgsea4[,2], fill = type)) +
            geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
            geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = input$boxplotDotss) +
            labs(x = "Group", y = paste(colnames(ssgsea4)[2], " expression", sep = ""),
                 title = paste("ssGSEA Expression Signature: ", colnames(ssgsea4)[2],sep = "")) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1),
                  text = element_text(size = input$boxplotFontss))
        }
        else {
          ggplot(ssgsea4, aes(factor(type), ssgsea4[,2], fill = type)) +
            geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
            geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = input$boxplotDotss) +
            labs(x = "Group", y = paste(colnames(ssgsea4)[2], " expression", sep = ""),
                 title = paste("ssGSEA Expression Signature: ", colnames(ssgsea4)[2],sep = "")) +
            theme_bw() +
            stat_compare_means(method = input$boxplotcompare) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1),
                  text = element_text(size = input$boxplotFontss))
        }
      }
    }
    #if MSigDB tab selected
    else if (input$tables == 1) {
      if (length(input$msigdbTable_rows_selected) > 0){
        ssgsea <- ssGSEAfunc()
        ssgsea2 <- as.data.frame(t(ssgsea))
        samporder <- meta[,1]
        ssgsea3 <- as.data.frame(ssgsea2[samporder,])
        colnames(ssgsea3)[1] <- colnames(ssgsea2[1])
        rownames(ssgsea3) <- samporder
        ssgsea4 <- ssgsea3 %>% 
          mutate(type = case_when(
            rownames(ssgsea3) == meta[,1] ~ meta[,2],
          ))
        ssgsea4 <- ssgsea4 %>%
          relocate(type)
        if (input$boxplotcompare == "none") {
          ggplot(ssgsea4, aes(factor(type), ssgsea4[,2], fill = type)) +
            geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
            geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = input$boxplotDotss) +
            labs(x = "Group", y = paste(colnames(ssgsea4)[2], " expression", sep = ""),
                 title = paste("ssGSEA Expression Signature: ",colnames(ssgsea4)[2],sep = "")) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1),
                  text = element_text(size = input$boxplotFontss))
        }
        else {
          ggplot(ssgsea4, aes(factor(type), ssgsea4[,2], fill = type)) +
            geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
            geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = input$boxplotDotss) +
            labs(x = "Group", y = paste(colnames(ssgsea4)[2], " expression", sep = ""),
                 title = paste("ssGSEA Expression Signature: ",colnames(ssgsea4)[2],sep = "")) +
            theme_bw() +
            stat_compare_means(method = input$boxplotcompare) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1),
                  text = element_text(size = input$boxplotFontss))
        }
      }
    }
    #if user generated GS table selected
    else if (input$tables == 5) {
      if (length(input$GStable.u_rows_selected) > 0){
        ssgsea <- ssGSEAfunc()
        ssgsea2 <- as.data.frame(t(ssgsea))
        samporder <- meta[,1]
        ssgsea3 <- as.data.frame(ssgsea2[samporder,])
        colnames(ssgsea3)[1] <- colnames(ssgsea2[1])
        rownames(ssgsea3) <- samporder
        ssgsea4 <- ssgsea3 %>% 
          mutate(type = case_when(
            rownames(ssgsea3) == meta[,1] ~ meta[,2],
          ))
        ssgsea4 <- ssgsea4 %>%
          relocate(type)
        if (input$boxplotcompare == "none") {
          ggplot(ssgsea4, aes(factor(type), ssgsea4[,2], fill = type)) +
            geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
            geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = input$boxplotDotss) +
            labs(x = "Group", y = paste(colnames(ssgsea4)[2], " expression", sep = ""),
                 title = paste("ssGSEA Expression Signature: ",colnames(ssgsea4)[2],sep = "")) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1),
                  text = element_text(size = input$boxplotFontss))
        }
        else {
          ggplot(ssgsea4, aes(factor(type), ssgsea4[,2], fill = type)) +
            geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
            geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = input$boxplotDotss) +
            labs(x = "Group", y = paste(colnames(ssgsea4)[2], " expression", sep = ""),
                 title = paste("ssGSEA Expression Signature: ",colnames(ssgsea4)[2],sep = "")) +
            theme_bw() +
            stat_compare_means(method = input$boxplotcompare) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1),
                  text = element_text(size = input$boxplotFontss))
        }
      }
    }
  })
  
  #render gene expression comparison scatter plot
  output$geneScatter0 <- renderPlotly({
    expr = updated_expr();
    meta = meta();
    
    #log if user designates
    if (input$logask == TRUE) {
      expr <- log2(expr + 1)
    }
    #transpose
    expr_t <- as.data.frame(t(expr))
    #reorder rowname to match meta for merging
    samporder <- meta[,1]
    expr_t2 <- as.data.frame(expr_t[samporder,])
    #add type
    expr_t3 <- expr_t2 %>% 
      mutate(type = case_when(
        rownames(expr_t2) == meta[,1] ~ meta[,2],
      ))
    expr_t3 <- expr_t3 %>%
      relocate(type)
    #user gene input
    gene1.u <- input$scatterG1
    gene2.u <- input$scatterG2
    #get columns and info based off user input
    gene1 <- expr_t3[,gene1.u]
    gene2 <- expr_t3[,gene2.u]
    if (input$logask == TRUE) {
      gene1.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene1.u)], " Expression (log2 +1)", sep = "")
      gene2.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene2.u)], " Expression (log2 +1)", sep = "")
    }
    else if (input$logask == FALSE) {
      gene2.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene2.u)], " Expression", sep = "")
      gene1.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene1.u)], " Expression", sep = "")
    }
    #plot
    p <- ggplot(expr_t3, aes(x = gene1, y = gene2,
                             color = type,
                             text = paste("</br> Sample: ", rownames(expr_t3),
                                          "</br> Type: ", type))) +
      geom_point() +
      theme_minimal() +
      labs(x = gene1.n, y = gene2.n,
           title = paste(gene1.n, " vs. ", gene2.n, sep = ''))
    
    ggplotly(p, tooltip = 'text')
    
  })
  
  
  ####----Download Handlers----####
  
  
  output$downloadClusters <- downloadHandler(
    filename = function() {
      paste(ProjectName,"_ClusterResults.tsv", sep = '')
    },
    content = function(file) {
      meta = meta();
      expr = updated_expr();
      A = as.matrix(expr)
      
      top_probes <- input$NumFeatures
      col_labels <- colnames(expr)
      isexpr <- rowSums(expr[,col_labels] > 1) >= length(col_labels)
      exp <- expr[isexpr,]
      mad <- NULL
      var <- NULL
      cv <- NULL
      var_type <- input$VarianceMeasure
      if (var_type == "MAD"){
        mad <- apply(log2(exp + 1), 1, mad)
        mad <- sort(mad, decreasing = T)
        mad <- head(mad, n = (top_probes +1))
        out <- cbind(names(mad), mad[names(mad)], exp[names(mad),])
        colnames(out) <- c("Gene", "MAD", colnames(exp))
        dataset <- exp[names(mad),]
        variable_gene_list <- names(mad)
      }
      if (var_type == "VAR"){
        var <- apply(log2(exp + 1), 1, var)
        var <- sort(var, decreasing = T)
        var <- head(var, n = (top_probes +1))
        out <- cbind(names(var), var[names(var)], exp[names(var),])
        colnames(out) <- c("Gene", "VAR", colnames(exp))
        dataset <- exp[names(var),]
        variable_gene_list <- names(var)
      }
      if (var_type == "CV"){
        cv <- apply(log2(exp + 1), 1, cv)
        cv <- sort(cv, decreasing = T)
        cv <- head(cv, n = (top_probes +1))
        out <- cbind(names(cv), cv[names(cv)], exp[names(cv),])
        colnames(out) <- c("Gene", "VAR", colnames(exp))
        dataset <- exp[names(cv),]
        variable_gene_list <- names(cv)
      }
      dataset <- log2(dataset + 1)
      zdataset <- apply(dataset, 1, scale)
      zdataset <- apply(zdataset, 1, rev)
      colnames(zdataset) <- names(dataset)
      dataset <- as.matrix(zdataset)
      dataset[is.na(dataset)] <- 0
      dataset = dataset[apply(dataset[,-1], 1, function(x) !all(x==0)),]
      minimum = -5;
      maximum = 5;
      if (abs(min(dataset)) > abs(max(dataset))) {
        dataset[dataset < -abs(max(dataset))] = -abs(max(dataset))
      } else {
        dataset[dataset > abs(min(dataset))] = abs(min(dataset))
      }
      results2 = hclust(dist(t(dataset)), method = input$ClusteringMethod)
      m = sort(cutree(results2, k=input$NumClusters))
      output = cbind(colnames(m), as.matrix(m))
      colnames(output) = c("Cluster")
      write_tsv(output,file)
    }
  )
  
  #download gene scatter plot expression data
  output$geneScatterDownload <- downloadHandler(
    filename = function() {
      #user gene input
      gene1.u <- input$scatterG1
      gene2.u <- input$scatterG2
      paste(gene1.u, "_vs_", gene2.u, "_Expression.tsv", sep = "")
    },
    content = function(file) {
      meta = meta();
      expr = updated_expr();
      A = as.matrix(expr)
      
      #transpose
      expr_t <- as.data.frame(t(expr))
      #reorder rowname to match meta for merging
      samporder <- meta[,1]
      expr_t2 <- as.data.frame(expr_t[samporder,])
      #add type
      expr_t3 <- expr_t2 %>% 
        mutate(type = case_when(
          rownames(expr_t2) == meta[,1] ~ meta[,2],
        ))
      expr_t3 <- expr_t3 %>%
        relocate(type)
      #user gene input
      gene1.u <- input$scatterG1
      gene2.u <- input$scatterG2
      #get columns and info based off user input
      gene1 <- expr_t3[,gene1.u]
      gene1.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene1.u)], " Expression", sep = "")
      gene2 <- expr_t3[,gene2.u]
      gene2.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene2.u)], " Expression", sep = "")
      #make table
      Sample <- rownames(expr_t3)
      Type <- expr_t3[,'type']
      gene1col <- expr_t3[,gene1.u]
      gene2col <- expr_t3[,gene2.u]
      scatterTab <- data.frame(Sample, Type, gene1col, gene2col)
      colnames(scatterTab)[c(3,4)] <- c(gene1.n, gene2.n)
      write_tsv(scatterTab, file)
    }
  )
  
  #download ssGSEA score table
  output$ssGSEAdownload <- downloadHandler(
    filename = function() {
      meta = meta();
      expr = updated_expr();
      A = as.matrix(expr)
      
      if (input$tables == 1) {
        paste('ssGSEAscore_',names(gs[(msigdb.gsea2[input$msigdbTable_rows_selected,3])]),'.tsv',sep = '')
      }
      else if (input$tables == 3) {
        paste('ssGSEAscore_',names(gs2[(GeneSet2[input$tab2table_rows_selected,1])]),'.tsv',sep = '')
      }
      else if (input$tables == 5) {
        paste('ssGSEAscore_',names(RDataListGen()[(GStable.ubg()[input$GStable.u_rows_selected,1])]),'.tsv',sep = '')
      }
    },
    content = function(file) {
      
      if (input$tables == 1) {
        GS <- gs[(msigdb.gsea2[input$msigdbTable_rows_selected,3])]
      }
      else if (input$tables == 3) {
        GS <- gs2[(GeneSet2[input$tab2table_rows_selected,1])]
      }
      else if (input$tables == 5) {
        GS <- RDataListGen()[(GStable.ubg()[input$GStable.u_rows_selected,1])]
      }
      ssgsea <- ssGSEAfunc()
      ssgsea2 <- as.data.frame(t(ssgsea))
      samporder <- meta[,1]
      ssgsea3 <- as.data.frame(ssgsea2[samporder,])
      colnames(ssgsea3)[1] <- colnames(ssgsea2[1])
      rownames(ssgsea3) <- samporder
      ssgsea4 <- ssgsea3 %>% 
        mutate(type = case_when(
          rownames(ssgsea3) == meta[,1] ~ meta[,2],
        ))
      ssgsea4 <- ssgsea4 %>%
        relocate(type)
      ssgsea4$sample <- rownames(ssgsea4)
      ssgsea5 <- ssgsea4[,c(3,1,2)]
      rownames(ssgsea5) <- 1:nrow(ssgsea5)
      write_tsv(ssgsea5, file)
    }
  )
  
  #render download button for DEG GMT
  output$DEGgmtDownload <- downloadHandler(
    filename = function() {
      paste(input$DEGfileName,".gmt",sep = "")
    },
    content = function(file) {
      meta = meta();
      expr = updated_expr();
      A = as.matrix(expr)
      
      top_probes <- input$NumFeatures
      col_labels <- colnames(expr)
      isexpr <- rowSums(expr[,col_labels] > 1) >= length(col_labels)
      exp <- expr[isexpr,]
      mad <- NULL
      var <- NULL
      cv <- NULL
      var_type <- input$VarianceMeasure
      if (var_type == "MAD"){
        mad <- apply(log2(exp + 1), 1, mad)
        mad <- sort(mad, decreasing = T)
        mad <- head(mad, n = (top_probes +1))
        out <- cbind(names(mad), mad[names(mad)], exp[names(mad),])
        colnames(out) <- c("Gene", "MAD", colnames(exp))
        dataset <- exp[names(mad),]
        variable_gene_list <- names(mad)
      }
      if (var_type == "VAR"){
        var <- apply(log2(exp + 1), 1, var)
        var <- sort(var, decreasing = T)
        var <- head(var, n = (top_probes +1))
        out <- cbind(names(var), var[names(var)], exp[names(var),])
        colnames(out) <- c("Gene", "VAR", colnames(exp))
        dataset <- exp[names(var),]
        variable_gene_list <- names(var)
      }
      if (var_type == "CV"){
        cv <- apply(log2(exp + 1), 1, cv)
        cv <- sort(cv, decreasing = T)
        cv <- head(cv, n = (top_probes +1))
        out <- cbind(names(cv), cv[names(cv)], exp[names(cv),])
        colnames(out) <- c("Gene", "VAR", colnames(exp))
        dataset <- exp[names(cv),]
        variable_gene_list <- names(cv)
      }
      dataset <- log2(dataset + 1)
      zdataset <- apply(dataset, 1, scale)
      zdataset <- apply(zdataset, 1, rev)
      colnames(zdataset) <- names(dataset)
      dataset <- as.matrix(zdataset)
      dataset[is.na(dataset)] <- 0
      dataset = dataset[apply(dataset[,-1], 1, function(x) !all(x==0)),]
      minimum = -5;
      maximum = 5;
      if (abs(min(dataset)) > abs(max(dataset))) {
        dataset[dataset < -abs(max(dataset))] = -abs(max(dataset))
      } else {
        dataset[dataset > abs(min(dataset))] = abs(min(dataset))
      }
      results2 = hclust(dist(t(dataset)), method = input$ClusteringMethod)
      m = sort(cutree(results2, k=input$NumClusters))
      output = cbind(colnames(m), as.matrix(m))
      colnames(output) = c("Cluster")
      A <- meta[which(meta[,2] == input$comparisonA2),1]
      B <- meta[which(meta[,2] == input$comparisonB2),1]
      
      #A <- meta[,1][meta[,2] == input$comparisonA2]
      #B <- meta[,1][meta[,2] == input$comparisonB2]
      mat <- expr[,c(A,B)]
      mat <- log2(mat + 1.0)
      groupAOther <- factor(c(rep("A", length(A)), rep("B", length(B))))
      designA <- model.matrix(~0 + groupAOther)
      fit <- lmFit(mat, design = designA)
      contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
      fit2 <- contrasts.fit(fit, contrast.matrix)
      fit2 <- eBayes(fit2)
      options(digits = 4)
      top1 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
      if (input$UpDnChoice == "UpAndDown_Regulated"){
        genes <- rownames(top1)[which(abs(top1$logFC) > abs(input$fc_cutoff2) & top1$adj.P.Val < input$p_cutoff2)]
        genes.h <- t(as.data.frame(head(genes, n=input$top_x2)))
        genes.h.gmt <- data.frame(paste(input$top_x2,input$UpDnChoice,"DEgenes", sep = ""),
                                  paste(input$UpDnChoice,"DEgenes",sep = ""),
                                  genes.h)
      }
      else if (input$UpDnChoice == "Up_Regulated"){
        genes <- rownames(top1)[which(top1$logFC > input$fc_cutoff2 & top1$adj.P.Val < input$p_cutoff2)]
        genes.h <- t(as.data.frame(head(genes, n=input$top_x2)))
        genes.h.gmt <- data.frame(paste(input$top_x2,input$UpDnChoice,"DEgenes", sep = ""),
                                  paste(input$UpDnChoice,"DEgenes",sep = ""),
                                  genes.h)
      }
      else if (input$UpDnChoice == "Down_Regulated"){
        genes <- rownames(top1)[which(top1$logFC < -abs(input$fc_cutoff2) & top1$adj.P.Val < input$p_cutoff2)]
        genes.h <- t(as.data.frame(head(genes, n=input$top_x2)))
        genes.h.gmt <- data.frame(paste(input$top_x2,input$UpDnChoice,"DEgenes", sep = ""),
                                  paste(input$UpDnChoice,"DEgenes",sep = ""),
                                  genes.h)
      }
      write_delim(genes.h.gmt, file, delim = '\t', col_names = F)
    }
  )
  
  #render download button for upreg gene pathways GMT
  output$UpRegPathDownloadgmt <- downloadHandler(
    filename = function() {
      paste("UpReg_",input$SelectedPathway,"_pathway.gmt",sep = "")
    },
    content = function(file) {
      meta = meta();
      expr = updated_expr();
      A = as.matrix(expr)
      
      top1 <- topgenereact()
      genes <- rownames(top1)[which(top1$adj.P.Val < 0.05 & top1$logFC > 1)]
      dbs <- listEnrichrDbs() 
      enrichRLive <- TRUE 
      if (is.null(dbs)) { 
        enrichRLive <- FALSE 
      }
      dbs <- input$SelectedPathway
      enriched <- enrichr(genes, dbs) # Plot top 20 GO-BP results ordered by P-value
      df <- enriched[[1]]
      path.genes <- as.data.frame(df[,'Genes'])
      genes.n <- c()
      for (i in path.genes[,1]) {
        g <- strsplit(i,";")
        for (j in g){
          genes.n <- c(genes.n, j)
        }
      }
      rm(g)
      genes.n <- t(as.data.frame(unique(genes.n)))
      gmt.df <- data.frame(paste("UpReg_",input$SelectedPathway,"_PathwayGenes",sep = ""),
                           input$SelectedPathway,genes.n)
      rownames(gmt.df) <- NULL
      write_delim(gmt.df, file, delim = '\t', col_names = F)
    }
  )
  
  #render download button for dnreg gene pathways GMT
  output$DnRegPathDownloadgmt <- downloadHandler(
    filename = function() {
      paste("DnReg_",input$SelectedPathway,"_pathway.gmt",sep = "")
    },
    content = function(file) {
      meta = meta();
      expr = updated_expr();
      A = as.matrix(expr)
      
      top1 <- topgenereact()
      genes <- rownames(top1)[which(top1$adj.P.Val < 0.05 & top1$logFC > 1)]
      dbs <- listEnrichrDbs() 
      enrichRLive <- TRUE 
      if (is.null(dbs)) { 
        enrichRLive <- FALSE 
      }
      dbs <- input$SelectedPathway
      enriched <- enrichr(genes, dbs) # Plot top 20 GO-BP results ordered by P-value
      df <- enriched[[1]]
      path.genes <- as.data.frame(df[,'Genes'])
      genes.n <- c()
      for (i in path.genes[,1]) {
        g <- strsplit(i,";")
        for (j in g){
          genes.n <- c(genes.n, j)
        }
      }
      rm(g)
      genes <- t(as.data.frame(unique(genes.n)))
      gmt.df <- data.frame(paste("DnReg_",input$SelectedPathway,"_PathwayGenes",sep = ""),
                           input$SelectedPathway,genes.n)
      rownames(gmt.df) <- NULL
      write_delim(gmt.df, file, delim = '\t', col_names = F)
    }
  )
  
  #render download button for upreg gene pathways
  output$UpRegPathDownload <- downloadHandler(
    filename = function() {
      paste("UpReg_",input$SelectedPathway,"_pathway.tsv",sep = "")
    },
    content = function(file) {
      meta = meta();
      expr = updated_expr();
      A = as.matrix(expr)
      
      top1 <- topgenereact()
      genes <- rownames(top1)[which(top1$adj.P.Val < 0.05 & top1$logFC > 1)]
      dbs <- listEnrichrDbs() 
      enrichRLive <- TRUE 
      if (is.null(dbs)) { 
        enrichRLive <- FALSE 
      }
      dbs <- input$SelectedPathway
      enriched <- enrichr(genes, dbs) # Plot top 20 GO-BP results ordered by P-value
      df <- enriched[[1]]
      write_delim(df, file, delim = '\t')
    }
  )
  
  #render download button for dnreg gene pathways
  output$DnRegPathDownload <- downloadHandler(
    filename = function() {
      paste("DnReg_",input$SelectedPathway,"_pathway.tsv",sep = "")
    },
    content = function(file) {
      meta = meta();
      expr = updated_expr();
      A = as.matrix(expr)
      
      top1 <- topgenereact()
      genes <- rownames(top1)[which(top1$adj.P.Val < 0.05 & top1$logFC < -1)]
      dbs <- listEnrichrDbs() 
      enrichRLive <- TRUE 
      if (is.null(dbs)) { 
        enrichRLive <- FALSE 
      }
      dbs <- input$SelectedPathway
      enriched <- enrichr(genes, dbs) # Plot top 20 GO-BP results ordered by P-value
      df <- enriched[[1]]
      write_delim(df, file, delim = '\t')
    }
  )
  
  #render Most Variable Gene download button
  output$MVGdownload <- downloadHandler(
    filename = function() {
      top_probes <- input$NumFeatures
      paste(top_probes,"_Most_Variable_Genes", ".tsv", sep = "")
    },
    content = function(file){
      meta = meta();
      expr = updated_expr();
      A = as.matrix(expr)
      
      top_probes <- input$NumFeatures
      col_labels <- colnames(expr)
      isexpr <- rowSums(expr[,col_labels] > 1) >= length(col_labels)
      exp <- expr[isexpr,]
      mad <- NULL
      var <- NULL
      cv <- NULL
      var_type <- input$VarianceMeasure
      if (var_type == "MAD"){
        mad <- apply(log2(exp + 1), 1, mad)
        mad <- sort(mad, decreasing = T)
        mad <- head(mad, n = (top_probes +1))
        out <- cbind(names(mad), mad[names(mad)], exp[names(mad),])
        colnames(out) <- c("Gene", "MAD", colnames(exp))
        dataset <- exp[names(mad),]
        variable_gene_list <- names(mad)
      }
      if (var_type == "VAR"){
        var <- apply(log2(exp + 1), 1, var)
        var <- sort(var, decreasing = T)
        var <- head(var, n = (top_probes +1))
        out <- cbind(names(var), var[names(var)], exp[names(var),])
        colnames(out) <- c("Gene", "VAR", colnames(exp))
        dataset <- exp[names(var),]
        variable_gene_list <- names(var)
      }
      if (var_type == "CV"){
        cv <- apply(log2(exp + 1), 1, cv)
        cv <- sort(cv, decreasing = T)
        cv <- head(cv, n = (top_probes +1))
        out <- cbind(names(cv), cv[names(cv)], exp[names(cv),])
        colnames(out) <- c("Gene", "VAR", colnames(exp))
        dataset <- exp[names(cv),]
        variable_gene_list <- names(cv)
      }
      variable_gene_list <- as.data.frame(variable_gene_list)
      colnames(variable_gene_list)[1] <- "Genes"
      variable_gene_list$Rank <- rownames(variable_gene_list)
      variable_gene_list <- variable_gene_list %>%
        select(Rank, Genes)
      write_delim(variable_gene_list, file, delim = '\t')
    }
  )
  
  #render Most Variable Gene GMT download button
  output$MVGdownloadgmt <- downloadHandler(
    filename = function() {
      top_probes <- input$NumFeatures
      paste(top_probes,"_Most_Variable_Genes", ".gmt", sep = "")
    },
    content = function(file){
      meta = meta();
      expr = updated_expr();
      A = as.matrix(expr)
      
      top_probes <- input$NumFeatures
      col_labels <- colnames(expr)
      isexpr <- rowSums(expr[,col_labels] > 1) >= length(col_labels)
      exp <- expr[isexpr,]
      mad <- NULL
      var <- NULL
      cv <- NULL
      var_type <- input$VarianceMeasure
      if (var_type == "MAD"){
        mad <- apply(log2(exp + 1), 1, mad)
        mad <- sort(mad, decreasing = T)
        mad <- head(mad, n = (top_probes +1))
        out <- cbind(names(mad), mad[names(mad)], exp[names(mad),])
        colnames(out) <- c("Gene", "MAD", colnames(exp))
        dataset <- exp[names(mad),]
        variable_gene_list <- names(mad)
      }
      if (var_type == "VAR"){
        var <- apply(log2(exp + 1), 1, var)
        var <- sort(var, decreasing = T)
        var <- head(var, n = (top_probes +1))
        out <- cbind(names(var), var[names(var)], exp[names(var),])
        colnames(out) <- c("Gene", "VAR", colnames(exp))
        dataset <- exp[names(var),]
        variable_gene_list <- names(var)
      }
      if (var_type == "CV"){
        cv <- apply(log2(exp + 1), 1, cv)
        cv <- sort(cv, decreasing = T)
        cv <- head(cv, n = (top_probes +1))
        out <- cbind(names(cv), cv[names(cv)], exp[names(cv),])
        colnames(out) <- c("Gene", "VAR", colnames(exp))
        dataset <- exp[names(cv),]
        variable_gene_list <- names(cv)
      }
      variable_gene_list <- as.data.frame(variable_gene_list)
      colnames(variable_gene_list)[1] <- "Genes"
      genes.n <- c()
      for (i in variable_gene_list[,1]) {
        g <- strsplit(i,";")
        for (j in g){
          genes.n <- c(genes.n, j)
        }
      }
      rm(g)
      genes <- t(as.data.frame(unique(genes.n)))
      gmt.df <- data.frame(paste(top_probes,"_Most_Variable_Genes", sep = ""),
                           "Most_Variable_Genes",genes)
      rownames(gmt.df) <- NULL
      write_delim(gmt.df, file, delim = '\t', col_names = F)
    }
  )
  
  #download button for leading edge genes
  output$LEGdownload <- downloadHandler(
    filename = function() {
      
      if (input$tables == 1){
        if (length(input$msigdbTable_rows_selected) > 0){
          GS <- msigdb.gsea2[input$msigdbTable_rows_selected,3]
        }
      }
      else if (input$tables == 3){
        if (length(input$tab2table_rows_selected) > 0){
          GS <- as.character(GeneSet2[input$tab2table_rows_selected,1])
        }
      }
      else if (input$tables == 5){
        if (length(input$GStable.u_rows_selected) > 0){
          GS <- as.character(user_gs_mirror()[input$GStable.u_rows_selected,1])
        }
      }
      paste(GS,"_leading_edge_genes", ".tsv", sep = "")
    },
    content = function(file){
      meta = meta();
      expr = updated_expr();
      A = as.matrix(expr)
      
      if (input$tables == 1){
        if (length(input$msigdbTable_rows_selected) > 0){
          GS <- msigdb.gsea2[input$msigdbTable_rows_selected,3]
        }
      }
      else if (input$tables == 3){
        if (length(input$tab2table_rows_selected) > 0){
          GS <- as.character(GeneSet2[input$tab2table_rows_selected,1])
        }
      }
      else if (input$tables == 5){
        if (length(input$GStable.u_rows_selected) > 0){
          GS <- as.character(user_gs_mirror()[input$GStable.u_rows_selected,1])
        }
      }
      res <- datasetInput()
      gsea.df <- as.data.frame(res@result)
      ## Subset core enriched genes
      genes1 <- as.matrix(gsea.df[which(gsea.df$Description==GS),"core_enrichment"])
      genes2 <- strsplit(genes1,"/")
      GeneSymbol <- as.data.frame(genes2, col.names = "GeneSymbol")
      GeneSymbol$Rank <- rownames(GeneSymbol)
      GeneSymbol <- GeneSymbol[,c("Rank","GeneSymbol")]
      write_delim(GeneSymbol, file, delim = '\t')
    })
  
  #render DEG table download button
  output$DEGtableDownload <- downloadHandler(
    filename = function() {
      paste("DEG_Table_",Sys.Date(), ".tsv", sep = "")
    },
    content = function(file){
      meta = meta();
      expr = updated_expr();
      A = as.matrix(expr)
      
      top_probes <- input$NumFeatures
      col_labels <- colnames(expr)
      isexpr <- rowSums(expr[,col_labels] > 1) >= length(col_labels)
      exp <- expr[isexpr,]
      mad <- NULL
      var <- NULL
      cv <- NULL
      var_type <- input$VarianceMeasure
      if (var_type == "MAD"){
        mad <- apply(log2(exp + 1), 1, mad)
        mad <- sort(mad, decreasing = T)
        mad <- head(mad, n = (top_probes +1))
        out <- cbind(names(mad), mad[names(mad)], exp[names(mad),])
        colnames(out) <- c("Gene", "MAD", colnames(exp))
        dataset <- exp[names(mad),]
        variable_gene_list <- names(mad)
      }
      if (var_type == "VAR"){
        var <- apply(log2(exp + 1), 1, var)
        var <- sort(var, decreasing = T)
        var <- head(var, n = (top_probes +1))
        out <- cbind(names(var), var[names(var)], exp[names(var),])
        colnames(out) <- c("Gene", "VAR", colnames(exp))
        dataset <- exp[names(var),]
        variable_gene_list <- names(var)
      }
      if (var_type == "CV"){
        cv <- apply(log2(exp + 1), 1, cv)
        cv <- sort(cv, decreasing = T)
        cv <- head(cv, n = (top_probes +1))
        out <- cbind(names(cv), cv[names(cv)], exp[names(cv),])
        colnames(out) <- c("Gene", "VAR", colnames(exp))
        dataset <- exp[names(cv),]
        variable_gene_list <- names(cv)
      }
      dataset <- log2(dataset + 1)
      zdataset <- apply(dataset, 1, scale)
      zdataset <- apply(zdataset, 1, rev)
      colnames(zdataset) <- names(dataset)
      dataset <- as.matrix(zdataset)
      dataset[is.na(dataset)] <- 0
      dataset = dataset[apply(dataset[,-1], 1, function(x) !all(x==0)),]
      minimum = -5;
      maximum = 5;
      if (abs(min(dataset)) > abs(max(dataset))) {
        dataset[dataset < -abs(max(dataset))] = -abs(max(dataset))
      } else {
        dataset[dataset > abs(min(dataset))] = abs(min(dataset))
      }
      results2 = hclust(dist(t(dataset)), method = input$ClusteringMethod)
      m = sort(cutree(results2, k=input$NumClusters))
      output = cbind(colnames(m), as.matrix(m))
      colnames(output) = c("Cluster")
      A <- meta[which(meta[,2] == input$comparisonA2),1]
      B <- meta[which(meta[,2] == input$comparisonB2),1]
      
      #A <- meta[,1][meta[,2] == input$comparisonA2]
      #B <- meta[,1][meta[,2] == input$comparisonB2]
      mat <- expr[,c(A,B)]
      mat <- log2(mat + 1.0)
      groupAOther <- factor(c(rep("A", length(A)), rep("B", length(B))))
      designA <- model.matrix(~0 + groupAOther)
      fit <- lmFit(mat, design = designA)
      contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
      fit2 <- contrasts.fit(fit, contrast.matrix)
      fit2 <- eBayes(fit2)
      options(digits = 4)
      top1 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
      top1$Gene <- rownames(top1)
      top1 <- top1 %>%
        relocate(Gene)
      write_delim(top1, file, delim = '\t')
    }
  )
  
  ##download button for Enrich Sig Table
  #output$enrich_sig_download <- downloadHandler(
  #    filename = function() {
  #        
  #        groupA <- input$comparisonA
  #        groupB <- input$comparisonB
  #        paste(input$SigTableChoice,".tsv", sep = "")
  #    },
  #    content = function(file){
  #        meta = meta();
  #        expr = updated_expr();
  #        A = as.matrix(expr)
  #        
  #        gsea.df <- as_tibble(get(ES_Tab_List[which(ES_Tab_List == paste("ES_table",match(input$SigTableChoice, SigNames),sep = ""))]))
  #        write_delim(gsea.df, file, delim = '\t')
  #    }
  #)
  
  #download button for user enriched signature table
  output$enrich_sig_download.u <- downloadHandler(
    
    filename = function() {
      groupA <- input$comparisonA
      groupB <- input$comparisonB
      paste("Enrich_Sig_Table_",groupA,"vs",groupB,".tsv", sep = "")
    },
    content = function(file) {
      expr = updated_expr()
      meta = meta()
      A = as.matrix(expr)
      
      if (input$tables == 3) {
        groupA <- meta[which(meta[,2] == input$comparisonA),1]
        groupB <- meta[which(meta[,2] == input$comparisonB),1]
        #groupA <- meta[,1][meta[,2] == input$comparisonA]
        #groupB <- meta[,1][meta[,2] == input$comparisonB]
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
        ##----GSEA----##
        gmt.i <- tab2
        gsea.res <- GSEA(s2n.matrix.s, TERM2GENE = gmt.i, verbose = F, pvalueCutoff = input$userPval)
        gsea.df <- as.data.frame(gsea.res@result)
        write_delim(gsea.df, file, delim = '\t')
      }
      else if (input$tables == 5) {
        groupA <- meta[which(meta[,2] == input$comparisonA),1]
        groupB <- meta[which(meta[,2] == input$comparisonB),1]
        #groupA <- meta[,1][meta[,2] == input$comparisonA]
        #groupB <- meta[,1][meta[,2] == input$comparisonB]
        
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
        ##----GSEA----##
        gmt.i <- GStable.ubg()
        gsea.res <- GSEA(s2n.matrix.s, TERM2GENE = gmt.i, verbose = F, pvalueCutoff = input$userPval)
        gsea.df <- as.data.frame(gsea.res@result)
        write_delim(gsea.df, file, delim = '\t')
      }
    }
  )
  
  
  ####----Text----####
  
  
  #NES and Pval output
  output$NESandPval <- renderText({
    meta = meta();
    expr = updated_expr();
    A = as.matrix(expr)
    
    if (input$tables == 1){
      if (length(input$msigdbTable_rows_selected) > 0){
        res <- datasetInput()
        gsea.df <- as.data.frame(res@result)
        GS = msigdb.gsea2[input$msigdbTable_rows_selected,3]
        NES = gsea.df$NES[which(gsea.df[,'ID']==GS)]
        Pval = gsea.df$pvalue[which(gsea.df[,'ID']==GS)]
        NES.o <- paste0("NES: ", NES)
        Pval.o <- paste0("Pvalue: ", Pval)
        if (NES > 0){
          UpOrDown <- paste("Based on the normalized enrichment score above, the", GS, "gene set is upregulated in", input$comparisonA , "group.")
        }
        else {
          UpOrDown <- paste("Based on the normalized enrichment score above, the", GS, "gene set is downregulated in", input$comparisonA , "group.")
        }
        paste(NES.o, Pval.o, UpOrDown, sep = '\n')
      }
      else if (length(input$msigdbTable_rows_selected) == 0){
        paste("Please select gene set from side panel table to begin.", sep = '')
      }
    }
    else if (input$tables == 3){
      if (length(input$tab2table_rows_selected) > 0){
        res <- datasetInput()
        gsea.df <- as.data.frame(res@result)
        GS = as.character(GeneSet2[input$tab2table_rows_selected,1])
        NES = gsea.df$NES[which(gsea.df[,'ID']==GS)]
        Pval = gsea.df$pvalue[which(gsea.df[,'ID']==GS)]
        NES.o <- paste0("NES: ", NES)
        Pval.o <- paste0("Pvalue: ", Pval)
        if (NES > 0){
          UpOrDown <- paste("Based on the normalized enrichment score above, the", GS, "gene set is upregulated in", input$comparisonA , "group.")
        }
        else {
          UpOrDown <- paste("Based on the normalized enrichment score above, the", GS, "gene set is downregulated in", input$comparisonA , "group.")
        }
        paste(NES.o, Pval.o, UpOrDown, sep = '\n')
      }
      else if (length(input$tab2table_rows_selected) == 0){
        paste("Please select gene set from side panel table to begin.", sep = '')
      }
    }
    else if (input$tables == 5){
      if (length(input$GStable.u_rows_selected) > 0){
        res <- datasetInput()
        gsea.df <- as.data.frame(res@result)
        GS <- as.character(user_gs_mirror()[input$GStable.u_rows_selected,1])
        NES = gsea.df$NES[which(gsea.df[,'ID']==GS)]
        Pval = gsea.df$pvalue[which(gsea.df[,'ID']==GS)]
        NES.o <- paste0("NES: ", NES)
        Pval.o <- paste0("Pvalue: ", Pval)
        if (NES > 0){
          UpOrDown <- paste("Based on the normalized enrichment score above, the", GS, "gene set is upregulated in", input$comparisonA , "group.")
        }
        else {
          UpOrDown <- paste("Based on the normalized enrichment score above, the", GS, "gene set is downregulated in", input$comparisonA , "group.")
        }
        paste(NES.o, Pval.o, UpOrDown, sep = '\n')
      }
      else if (length(input$GStable.u_rows_selected) == 0){
        paste("Please select gene set from side panel table to begin.", sep = '')
      }
    }
  })
  
  output$VolGroupsText <- renderText({
    paste("This volcano plot is comparing group A: ",input$comparisonA2, " and group B: ",input$comparisonB2, ".\nGenes with a positive log fold change are upregulated in the ",
          input$comparisonA2, " group.\nGenes with a negative log fold change are upregulated in the ",input$comparisonB2, " group.", sep = "")
  })
  
  output$MAGroupsText <- renderText({
    paste("This MA plot is comparing group A: ",input$comparisonA2, " and group B: ",input$comparisonB2, ".\nGenes with a positive log fold change are upregulated in the ",
          input$comparisonA2, " group.\nGenes with a negative log fold change are upregulated in the ",input$comparisonB2, " group.", sep = "")
  })
  
  output$upregpath_text <- renderText({
    paste("Genes in these enriched terms are upregulated in group A: ", input$comparisonA2.path, " group.", sep = "")
  })
  
  output$downregpath_text <- renderText({
    paste("Genes in these enriched terms are upregulated in group B: ", input$comparisonB2.path," group", sep = "")
  })
  
  output$degtext <- renderText({
    paste("This table represents differentially expressed genes when comparing group A: ",input$comparisonA2.DEG," and group B: ",input$comparisonB2.DEG,
          ".\nGenes with a positive logFC, indicate an upregulation in group A: ", input$comparisonA2.DEG,
          ".\nGenes with a negative logFC, indicate an upregulation in group B: ", input$comparisonB2.DEG,
          ".\nThe 'AveExpr' column represents the log transformed average expression between group A: ",input$comparisonA2.DEG," and group B: ",input$comparisonB2.DEG,".", sep = "")
  })
  
  output$UpRegPathLabel <- renderUI({
    
    FC <- input$pathFC
    h3(paste("(Up-regulated pathway (> ",FC," logFC)",sep = ""))
    
  })
  
  output$DnRegPathLabel <- renderUI({
    
    FC <- input$pathFC
    h3(paste("Down-regulated pathway (> -",FC," logFC)",sep = ""))
    
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
