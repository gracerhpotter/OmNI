
# RUN APP

source("ui.R")
source("server.R")

grDevices::graphics.off()

shinyApp(ui = "ui.R", server = "server.R")

################################################################################