
# INSTALL CRAN/BIOCONDUCTOR PACKAGES
cran_packages <- c("Biobase", "BiocManager", "bsicons", "bslib", "colourpicker",
                   "colourvalues", "cowplot", "devtools", "DT", "ggcorrplot", 
                   "ggfortify", "ggpubr", "ggrepel", "ggridges", "ggupset", 
                   "ggvenn", "heatmaply", "htmltools", "igraph",
                   "influential", "kableExtra", "KEGGREST", "markdown", 
                   "openxlsx", "pacman", "pcaMethods", "pheatmap", 
                   "plotly", "psych", "readr", "rmarkdown", "shiny", 
                   "shinycssloaders", "tidyverse", "uwot", "VIM", "visNetwork")

for (package in cran_packages){
  if (!(package %in% installed.packages()[,"Package"])){install.packages(package, update = TRUE, ask = FALSE)}
}

# INSTALL BIOCONDUCTOR PACKAGES
bioconductor_packages <- c('clusterProfiler', 'ComplexHeatmap', 'enrichplot',
                           'EnsDb.Hsapiens.v86', 'EnsDb.Mmusculus.v79',
                           'GeneTonic', "Glimma", "NormalyzerDE", 'org.Ce.eg.db', 
                           'org.Dm.eg.db', 'org.Dr.eg.db', "org.Hs.eg.db", 'org.Mm.eg.db', 
                           'org.Rn.eg.db', 'org.Sc.sgd.db', "pathview", "pcaMethods", 
                           'topGO', "vsn")
for (package in bioconductor_packages){
  if (!(package %in% installed.packages()[,"Package"])){BiocManager::install(package, update = TRUE, ask = FALSE)}
}

# INSTALL GITHUB PACKAGES
github_packages <- c('IOR-Bioinformatics/PCSF')

for (package in github_packages){
  if (!(sub(".*/", "", package) %in% installed.packages()[,"Package"])){devtools::install_github(package)}
}

# LOAD ONLY NECESSARY PACKAGES
pacman::p_load(shiny, tidyverse, Biobase, NormalyzerDE, pathview)

rm(cran_packages, bioconductor_packages, github_packages, package)

# reinstall to fix interaction with the Matrix package since update to >= 1.6.2 or UMAP doesn't work
# pacman::p_install(irlba)

# note that MAC needs quartz to be able to run cairo package which is needed for heatmap. 
# XQuartz: https://www.xquartz.org/
