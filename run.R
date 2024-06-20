options(repos = list(CRAN="http://cran.rstudio.com/"))

if (!("pacman" %in% installed.packages()[,"Package"])){install.packages('pacman', update = TRUE, ask = FALSE)}
pacman::p_load(shiny)

Sys.setenv(RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools")

folder_address = '/Users/pottegra/Documents/ShinyApp_V17'
runApp(folder_address, launch.browser = TRUE)