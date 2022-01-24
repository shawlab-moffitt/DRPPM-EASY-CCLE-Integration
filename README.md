# DRPPM-EASY-Database-Integration

This is an extention of the [DRPPM Expression Analysis ShinY (EASY) Integration App](https://github.com/shawlab-moffitt/DRPPM-EASY-Integration) which allows the user to integrate their data witht the data of a large project database, for example we use the [Cancer Cell Line Encyclopedia (CCLE)](https://sites.broadinstitute.org/ccle/) and a [Lung Squamous Cell Carcinoma](https://www.sciencedirect.com/science/article/pii/S0092867421008576?via%3Dihub) study from Clinical Proteomic Tumor Analysis Consortium (CPTAC). Through the integration of data sets users may perform expression level differences and in-depth reciprocal Gene Set Enrichment Analysis (GSEA). This R Shiny app is very similar in features to the original Integration app, with the addition of a sample selection tab which allows the user to subset samples or a study of their choice. Based on this selection, the expression and meta data will be subset and imported in the back end to the app for further analysis and visualization. 


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
