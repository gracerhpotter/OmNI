#-------------------------------------------------
#' Draw Range vs Rank
#'
#' This function generates a plot of Rank versus abundance. Not currently in Shiny app.
#' 
#' @param eset Expression set object
#' @param type String name of the dataset
#' 
#' @output ggplot object

drawRange <- function(eset,
                      type){
  
  matrix <- Biobase::exprs(eset)
  matrix <- matrix[rowMeans(matrix)!=0,]
  
  plot_data <- data.frame(Abundance = rowMeans(matrix),
                          Rank = rank(-rowMeans(matrix)) )
                                
  plot <- ggplot(data = plot_data, aes(x = Rank, y = Abundance)) + 
            geom_point(size = 0.8) + 
            labs(title = paste(type, "Range", sep = "")) +
            scale_y_continuous(trans = "log10",
                               breaks = scales::breaks_log()) +
            theme_bw() + theme(plot.margin = margin(10, 20, 10, 10), axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(plot)
}

