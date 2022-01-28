# DRPPM-EASY-Database-Integration

This is an extention of the [DRPPM Expression Analysis ShinY (EASY) Integration App](https://github.com/shawlab-moffitt/DRPPM-EASY-Integration) which allows the user to integrate their data witht the data of a large project database, for example we use the [Cancer Cell Line Encyclopedia (CCLE)](https://sites.broadinstitute.org/ccle/) and a [Lung Squamous Cell Carcinoma](https://www.sciencedirect.com/science/article/pii/S0092867421008576?via%3Dihub) study from Clinical Proteomic Tumor Analysis Consortium (CPTAC). Through the integration of data sets, users may perform expression level differences and in-depth reciprocal Gene Set Enrichment Analysis (GSEA). This R Shiny app is very similar in features to the original Integration app, with the addition of a sample selection tab which allows the user to subset samples or a study of their choice. Based on this selection, the expression and meta data will be subset and imported in the back end to the app for further analysis and visualization. 


# Installation

* Download ZIP file from https://github.com/shawlab-moffitt/DRPPM-EASY-Database-Integration
* Unzip and load into directory as a project in R Studio
* Open the ‘App.R’ script and write in user input files and options as directed at the top of the script
  * ‘App.R’ script begins with example files loaded in from the ExampleData folder
* Press ‘Run App’ button in R Studio to run in application or browser window and enjoy!
  * The app script will install any missing packages that the user may not have locally

# Requirements

* `R` - https://cran.r-project.org/src/base/R-4/
* `R Studio` - https://www.rstudio.com/products/rstudio/download/

# R Dependencies

|  |  |  |  |  |
| --- | --- | --- | --- | --- |
| shiny_1.6.0 | shinythemes_1.2.0 | shinyjqui_0.4.0 | shinycssloaders_1.0.0 | tools_4.1.0 |
| dplyr_1.0.7 | tidyr_1.1.3 | readr_2.0.1 | tibble_3.1.3 | DT_0.18 |
| ggplot2_3.3.5 | plotly_4.9.4.1 | enrichplot_1.12.2 | pheatmap_1.0.12 | ggrepel_0.9.1 |
| enrichR_3.0 | limma_3.48.3 | clusterProfiler_4.0.5 | limma_3.48.3 | GSVA_1.40.1 |
| BiocManager_1.30.16 | reshape2_1.4.4 | ggpubr_0.4.0 |  |  |

# Required Files

* **MSigDB Gene Set Names:**
  * These [gene set files](https://github.com/shawlab-moffitt/DRPPM-EASY-Database-Integration/tree/main/GeneSet_data) were gathered from the [Molecular Signatures Database (MSigDB)](http://www.gsea-msigdb.org/gsea/msigdb/index.jsp) as separate collections and processed through R to generate a master gene set file with catagorical labels to use for GSEA and ssGSEA analysis.
  * This is used mainly for the UI for gmt category selection.
* **MSigDB Gene Set RData List:**
  * The RData gene set list is a more refined format of the gene set table.
  * This is a named list with over 32,000 gene sets from MSigDB paired with the genes they consist of.
  * This list is used for the back end analysis.
* **User Data Input Files:**
  * Expression Matrix:
    * Must be tab delimited with gene names as symbols located in the first column with subsequent columns consiting of the sample name as the header and expression data down the column
    * The current App expects lowly expressed genes filtered out and normalized data either to FPKM or TMM
      * Larger files might inflict memory issues for you local computer
  * Meta File:
    * Must be tab delimited with two columns. First column of sasmple names and second column as phenotype grouping of the samples
* **Project Database Files:**
  * Meta Data:
    * Three column, tab-delimited, format with columns in the order of Sample Name, Meta Group, Sample Type
    * This is used to group the expression data into comparison groups for differential expression analysis
  * Meta Selector Data (Optional):
    * This is used when the expression data is able to be subset for analysis
      * In the case of the CCLE example we can subset the expression data based on disease or lineage before grouping with the meta file
    * This is a two column, tab-delimited, file with the first column being the meta groups (as seen in the second column of the main meta data) and the second column is either "Phenotype" or "Selector"
      * "Selector" designates if the meta group is used to subset the expression data
      * "Phenotype" designates if the meta group is used to group the expression data
  * Expression Data:
    * This should be a tab-delimited metrix with the columns labeling each sample and the rows labeling the feature/gene and the cell values should be un-logged expression values
  * Name Map (Optional):
    * This should be a two column, tab-delimited, format with the first column being the sample name used in the expression and meta data and the second column annotating another possible name for the sample
    * This is usefull when samples may have longer names or alternate descriptive names that may identify them

# App Preparation

Below shows the only section of the script that needs to be updated, currently written for setting up the CCLE Analysis App. The user may enter the project name, the meta file, the meta selector file, the expression matrix file, and a name map file. For the optional files (meta selector and name map) you may leave the contents empty and it will run without them. Once these files are added the app may be run as is.

```
####----Project Name----####
ProjName <- 'CCLE'
####----File Names----####
##--Database Files--##
#Meta
db_meta_file <- '~/R/DRPPM-EASY-Database-Integration-main/CCLE_data/CCLE_meta_melt_nosub.zip'
#Meta Selector File
db_meta_selec_file <- '~/R/DRPPM-EASY-Database-Integration-main/CCLE_data/CCLE_meta_selector_nosub.tsv'
#Expression Data
db_expr_file <- '~/R/DRPPM-EASY-Database-Integration-main/CCLE_data/CCLE_expr_trim_NewName.zip'
#Name Map File
db_namemap_file <- '~/R/DRPPM-EASY-Database-Integration-main/CCLE_data/CCLE_NameMap.tsv'
```

# App Features

## Data Input and Sample Selection Steps

### Step 1: User Data Input

![alt text](https://github.com/shawlab-moffitt/DRPPM-EASY-Database-Integration/blob/main/App_Demo_Pictures/EASY_DB_INT_dataupload.PNG?raw=true)

1. Write in a name to identify your data
2. User upload of expression matrix and meta file

### Step 2: Project Data Selection

![alt text](https://github.com/shawlab-moffitt/DRPPM-EASY-Database-Integration/blob/main/App_Demo_Pictures/EASY_DB_INT_selectdata.PNG?raw=true)

1. This section to subset expression data will only show if the user includes a meta selector file that is able to subset the expression data based on a variable
   * In the case of the CCLE example the user may select to subset the expression data based on lineage or disease type
2. The condition selection designated which meta group to group the expression data with
3. The user may choose to log2 transform the expression data
4. The label is automatically filled with the Project Name given in the script but is able to be adjusted here

### Step 3: Designate Comparison Groups

![alt text](https://github.com/shawlab-moffitt/DRPPM-EASY-Database-Integration/blob/main/App_Demo_Pictures/EASY_DB_INT_step3.PNG?raw=true)

1. Meta group 1 and 2 of the user data may be designated here
2. Meta group 1 and 2 of the project data may be designated here
   * The samples of group 1 from the user and project data will be grouped together and the samples of group 2 from the user and project data will be grouped together
   * These new groups 1 and 2 will be used for the downstream analysis
3. The new groups can be labeled here
4. The new meta table that is generated and shown may be names and downloaded for further use
