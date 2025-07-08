# Omics Notebook Interactive [OmNI]

An interactive shiny application version of the original Omics Notebook, with increased functionalities.

## Dependencies

Dependencies for the application are installed and loaded via the `loads.R` file. The first time opening the application may take 10-15+ minutes if your system does not have many of the packages installed. 

## Downloading the ZIP File

All necessary packages should load on opening the app via the `loads.R` other than `shiny`, which needs to be installed before opening the application.

```
install.packages("shiny") # install shiny
```

Shiny applications can be opened via the R function `shiny::runApp(appDir = getwd())` (make sure your working directory is set to the folder all of the application files are in, or provide the path to the folder instead of using `getwd()`) or by opening the `app.R`, `ui.R`, or `server.R` files in RStudio and clicking on the "Run App" button that appears in the top right corner of the code window. Additionally the app contains two files which can be used to run the app from a shortcut via the command line, `run.R` and `runOmicsNotebook.bat`. These were created to work on a Windows machine and to output an Rout file by username. Individual use cases can adjust output and file paths.

## Opening in R via GitHub

Alternatively, shiny apps can be opening without downloading from GitHub using the `shiny::runGitHub()` function.

```
shiny::runGitHub("OmNI", "gracerhpotter")
```

## Development & Data Information

Additional information on the contents of each file, and the source of referenced databases for enrichment, PCSF, and S-Score is in the `info` folder in `devInfo.txt`. Additional reference scripts for the creation of some the files are also in that folder.

## Example Data & Annotation Files

Example data and annotation files are present in the `example` folder, but can also be found in a [public Google Drive folder](https://drive.google.com/drive/folders/1lyzmIhorrZy_CKuxabi1Bv1cLHIblJhk?usp=drive_link).

## Preload Packages

To pre-install/pre-load necessary packages the script below can be run in R.

```
# INSTALL CRAN/BIOCONDUCTOR PACKAGES
cran_packages <- c("Biobase", "BiocManager", "bsicons", "bslib", "colourpicker",
                   "colourvalues", "cowplot", "devtools", "DT", "ggcorrplot", 
                   "ggfortify", "ggpubr", "ggrepel", "ggridges", "ggupset", 
                   "ggvenn", "heatmaply", "htmltools", "igraph",
                   "influential", "kableExtra", "KEGGREST", "limma", "markdown", 
                   "openxlsx", "pacman", "pcaMethods", "pheatmap", 
                   "plotly", "psych", "readr", "rmarkdown", "shiny", 
                   "shinycssloaders", "tidyverse", "uwot", "VIM", "visNetwork")
for (package in cran_packages){
  if (!(package %in% installed.packages()[,"Package"])){install.packages(package, update = TRUE, ask = FALSE)}
}

# INSTALL BIOCONDUCTOR PACKAGES
bioconductor_packages <- c("clusterProfiler", "ComplexHeatmap", "enrichplot",
                           "EnsDb.Hsapiens.v86", "EnsDb.Mmusculus.v79",
                           "GeneTonic", "Glimma", "NormalyzerDE", "org.Ce.eg.db", 
                           "org.Dm.eg.db", "org.Dr.eg.db", "org.Hs.eg.db", "org.Mm.eg.db", 
                           "org.Rn.eg.db", "org.Sc.sgd.db", "pathview", "pcaMethods", 
                           "topGO", "vsn")
for (package in bioconductor_packages){
  if (!(package %in% installed.packages()[,"Package"])){BiocManager::install(package, update = TRUE, ask = FALSE)}
}

# INSTALL GITHUB PACKAGES
github_packages <- c('IOR-Bioinformatics/PCSF')
for (package in github_packages){
  if (!(sub(".*/", "", package) %in% installed.packages()[,"Package"])){devtools::install_github(package)}
}

# note that MAC needs quartz to be able to run cairo package which is needed for heatmaps. 
# XQuartz: https://www.xquartz.org/
```
