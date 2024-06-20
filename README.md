# Interactive Omics Notebook

An interactive shiny application version of the original Omics Notebook, with increased functionalities.

## Packages

All necessary packages should load on opening the app via the `loads.R` other than `shiny`, which needs to be installed before opening the application.


```
install.packages("shiny") # install shiny
```

## Opening the App

Shiny applications can be opened via the R function `shiny::runApp(folder_address)` or by opening the `app.R`, `ui.R`, or `server.R` files in RStudio and clicking on the "Run App" button that appears in the top right corner of the code window.