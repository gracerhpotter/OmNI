
# RUN APP

source("ui.R")
source("server.R")

grDevices::graphics.off()

shinyApp(ui = "ui.R", server = "server.R")

################################################################################
# install.packages('C:/Users/pottegra/Downloads/DelayedArray_0.28.0.zip', repos=NULL, type='source')