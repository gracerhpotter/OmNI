# Interactive Omics Notebook

An interactive shiny application version of the original Omics Notebook, with increased functionalities.

## Packages

All necessary packages should load on opening the app via the `loads.R` other than `shiny`, which needs to be installed before opening the application.


```
install.packages("shiny") # install shiny
```

## Opening the App

Shiny applications can be opened via the R function `shiny::runApp(appDir = getwd())` (make sure your working directory is set to the folder all of the application files are in, or provide the path to the folder instead of using `getwd()`) or by opening the `app.R`, `ui.R`, or `server.R` files in RStudio and clicking on the "Run App" button that appears in the top right corner of the code window. Additionally the app contains two files which can be used to run the app from a shortcut via the command line, `run.R` and `runOmicsNotebook.bat`. These were created to work on a Windows machine and to output an Rout file by username. Individual use cases can adjust output and file paths.

## Development & Data Information

Additional information on the contents of each file, and the source of referenced databases for enrichment, PCSF, and S-Score is in the `info` folder in `devInfo.txt`. Additional reference scripts for the creation of some the files are also in that folder.