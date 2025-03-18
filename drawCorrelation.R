
#' @title Correlation Plot
#' 
#' @description
#' Generates a correlation plot based on the expression set object.
#' 
#' @param eset expression set input - eset object from omics data
#' @param high_variance logical input - whether to plot only high variance subset of features
#' @param plot logical input - whether to return plot output
#' @param table logical input - whether to return correlation output
#' @param numbers logical input - whether to include correlation coefficients in plot
#' @param percent_var numeric input - for high variance plot, what percent of top variance to include
#' @param neg_color string input - hex code for color for negative correlation
#' @param pos_color string input - hex code for color for positive correlation
#' @param title string input - optional title addition
#' 
#' @return correlation plot object from ggcorrplot.
#'

variationPlot <- function(eset,
                          high_variance = FALSE,
                          plot = TRUE,
                          table = FALSE,
                          numbers = FALSE,
                          percent_var = 0.2,
                          neg_color = "",
                          pos_color = "",
                          title = "") {
  
  colors = c("white", pos_color)
  
  # CALCULATE CORRELATION COEFFICIENT MATRIX
  emat <- Biobase::exprs(eset)
  emat_sel <- emat
  emat_sel <- na.omit(scale(emat_sel))
  correlation_matrix <- stats::cor(emat_sel)
  
  # assign("correlation_matrix", correlation_matrix, envir = .GlobalEnv)
  
  if (any(correlation_matrix < 0)) {colors = c(neg_color, "white", pos_color)}
  
  # SET TEXT SIZE FOR NUMBER
  if (ncol(emat) <= 10) {number_size = 12
  } else if (ncol(emat) > 10 && ncol(emat) <= 16) {number_size = 10
  } else if (ncol(emat) > 16) {number_size = 8
  } else if (ncol(emat) > 30) {number_size = 4}
  
  paletteLength <- 100
  breaks <- c(seq(min(correlation_matrix), 0, length.out=ceiling(paletteLength/2) + 1), 
                seq(max(correlation_matrix)/paletteLength, max(correlation_matrix), length.out=floor(paletteLength/2)))
  
  if (high_variance == FALSE) {
    # GENERATE CORRELATION HEATMAP PLOT
    correlation_plot <- pheatmap::pheatmap(correlation_matrix,
                                           clustering_method = "average",
                                           scale = "none",
                                           color = colorRampPalette(colors)(100),
                                           main = paste("Correlation Plot: Total Sample", "\n", title, sep = ""),
                                           fontsize = 12,
                                           show_colnames = TRUE,
                                           show_rownames = TRUE,
                                           display_numbers = numbers,
                                           number_color = "black",
                                           border_color = "white",
                                           fontsize_number = number_size,
                                           breaks = breaks)
    
    message("** GENERATED BASIC CORRELATION PLOT **")
  
  } else {
    # CALCULATE CORRELATION COEFFICIENT MATRIX FOR SUBSET OF HIGHLY VARIABLE FEATURES
    subset_variable <- round((nrow(emat_sel) * percent_var))
    vars <- sort(apply(emat_sel, 1, var), decreasing = TRUE)[1:subset_variable]
    emat_sel <- emat_sel[names(vars),]
    high_var_correlation_matrix <- stats::cor(emat_sel)
    
    # GENERATE CORRELATION HEATMAP PLOT
    correlation_plot <- pheatmap::pheatmap(high_var_correlation_matrix,
                                           clustering_method = "average",  
                                           scale = "none",  
                                           color = colorRampPalette(colors)(100),  
                                           main = paste("Correlation Plot: Highest Variance Proteins (Top ", (percent_var * 100), "%)", "\n", title, sep = ""),
                                           fontsize = 12,
                                           show_colnames = TRUE,  
                                           show_rownames = TRUE,
                                           display_numbers = numbers,
                                           number_color = "black",
                                           border_color = "white",
                                           fontsize_number = number_size,
                                           breaks = breaks) 
    
    message("** GENERATED HIGH VARIANCE CORRELATION PLOT **")
  }
  
  return(correlation_plot);
}

#' #-------------------------------------------------
#' #' Draw XY Correlation Plot
#' #'
#' #' This function generates correlation plots
#' #' 
#' #' @param item_list a list of tibbles
#' #' @param item_name a name for the items
#' #' @param outputpath output file path for plots
#' #' @param file_name suffix for plot filename
#' #' 
#' #' @examples
#' #' 
#' #' @import ggplot2
#' #' @export
#' 
#' drawXYCorr <- function(item_list, item_name, outputpath=output_plots_path, file_name="", subset_genes=FALSE){
#'   
#'   plotXYCorr <- function(sublist = NULL,
#'                          expanded_title = NULL){
#'     
#'     for (j in 1:(length(item_list) - 1)){
#'       for (k in (j+1):length(item_list)){
#'         # make data to plot
#'         plot_data <- merge(item_list[[j]], item_list[[k]], by.x=item_name, by.y=item_name)
#'         colnames(plot_data) <- c(item_name, names(item_list)[j], names(item_list)[k])
#'         title=paste0("Avg. ",item_name, " ", file_name)
#'         if(!is.null(sublist)){
#'           plot_data <- plot_data[which(plot_data[,item_name] %in% sublist),]
#'           if(!is.null(expanded_title)) title=paste0(title," in ",expanded_title,": \n",
#'                                                     names(item_list)[j]," vs. ",names(item_list)[k])
#'         }
#'         if(nrow(plot_data)==0){ next; }
#'         # density colors
#'         x<-grDevices::densCols(plot_data[,2], plot_data[,3],colramp=grDevices::colorRampPalette(c("black","white")))
#'         plot_data$Density <- grDevices::col2rgb(x)[1,] + 1L
#'         # make plot
#'         corr_coef <- cor(plot_data[,2], plot_data[,3], method="pearson");
#'         label<-paste("italic(r) == ",round(corr_coef, digits=2), sep="");
#'         plot <- ggplot(data=plot_data, aes(x=plot_data[,2], y=plot_data[,3], color=Density)) + 
#'           geom_point() + viridis::scale_color_viridis(direction=-1) + labs(x=names(item_list)[j], y=names(item_list)[k]) +
#'           labs(title=title) +
#'           theme_bw() + annotate("text",x=-Inf, y=Inf, label=label, color="black", vjust=1.5, hjust=-0.4 , parse=TRUE) +
#'           geom_smooth(method=lm, se=FALSE, color="red", formula=y~x)
#'         print(plot+theme(legend.position="none"))
#'         gridExtra::grid.arrange(g_legend(plot))
#'       }
#'     }
#'   }
#'   
#'   output_filename<-file.path(outputpath,paste("Correlation_Plots_",item_name,"_",file_name,".pdf", sep=""))
#'   pdf(output_filename, width=3, height=3 )
#'   
#'   plotXYCorr()
#'   
#'   # Repeat for subset genes
#'   if( class(subset_genes)!="logical"){
#'     for( i in 1:length(subset_genes) ){ try({
#'       plotXYCorr(subset_genes[[i]],names(subset_genes)[i])
#'     }) }
#'   }
#'   dev.off()
#' }