# USE RENV
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

message("Configuring OmNI environment. This may take a few minutes...")
renv::restore(prompt = FALSE)

pacman::p_load(shiny, tidyverse, Biobase, NormalyzerDE, pathview)
