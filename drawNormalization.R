#-------------------------------------------------------------------------------
#' @title Intensity Normalization Plots
#'
#' @description
#'   This function takes omics data set and performs normalization and QC plotting 
#'   to inspect data. Can be expanded with additional normalization methods.
#' 
#' @param eset_prenorm an ExpressionSet object with omics data before normalization
#' @param eset_postnorm an ExpressionSet object with omics data after normalization
#' @param norm string input - normalization method ("none", "quantile", "loess", "median", "z transform", "mad", "irs")
#' @param type string input - type of data to be processed
#' @param plottype string input - normalization plot to generate ("Boxplot", "Density", "MA", "RLE", "QQ")
#' @param zero_cutoff numeric input - filter out rows with more than this percentage of zeros/NAs, set to 0.3
#' @param norm_by_batches logical input - whether to do batch normalization
#' @param ma_array numeric input - if generating MA plot, which plot to make
#' @param qq_column string input - if generating QQ plot, which sample to use
#'
#' @return A plot object showing pre/post normalization values
#' 

intensityNormPlots <- function(eset_prenorm,
                               eset_postnorm,
                               norm, 
                               type, 
                               plottype = NULL,
                               zero_cutoff = 0.3,
                               norm_by_batches = F,
                               ma_array = 1,
                               qq_column = ""){
  
  # SET COLORS
  col_palette <- grDevices::rainbow(length(levels(as.factor(Biobase::pData(eset_prenorm)$SampleName))))
  
  # LOAD PRE-NORM ESET DATA
  eset_matrix <- Biobase::exprs(eset_prenorm)
  
  df <- data.frame(reshape2::melt(eset_matrix, id.vars = NULL));
  colnames(df) <- c("Feature","Sample", "Intensity");
  
  
  # GENERATE UN-NORMALIZED PLOTS WITH PRE-NORM DATA ----------------------------
  
  # BOX/VIOLIN PLOT
  plot1 <- ggplot(df, aes(x = Sample, y = Intensity, fill = Sample)) + 
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12), 
          plot.title = element_text(size = 18),
          legend.position = "none") +  
    labs(title = paste(type, " log2 Intensity: \nBefore Normalization", sep = ''))
  
  # DENSITY PLOT
  plot2 <- ggplot(df, aes(x = Intensity, colour = Sample)) + 
    geom_density() + 
    theme_bw() + 
    theme(plot.title = element_text(size = 18),
          legend.position = "none") +
    labs(title = paste(type, " log2 Intensity: \nBefore Normalization", sep = ''));
  
  # QQ PLOT
  QQ_column <- which(colnames(eset_matrix) == qq_column)
  plot5 <- ggpubr::ggqqplot(eset_matrix[, QQ_column]) +
    theme_bw() + 
    theme(plot.title = element_text(size = 18),
          legend.position = "none") +
    labs(title = paste(type, " ", qq_column, " QQ Plot ", toupper(norm)," Post-Normalization", sep = ''));
  
  # LOAD POST-NORM DATA
  eset_matrix_norm <- Biobase::exprs(eset_postnorm)
  
  # GENERATE NORMALIZED PLOTS WITH POST-NORM DATA ------------------------------
  if (!grepl("none", norm)) {
    df <- data.frame(reshape2::melt(eset_matrix_norm, id.vars = NULL));
    colnames(df) <- c("Feature","Sample", "Intensity");
    df[,"Sample"] <- as.character(df[,"Sample"])
    
    # BOX/VIOLIN PLOT
    plot3 <- ggplot(df, aes(x = Sample, y = Intensity, fill = Sample)) + 
      geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
      theme_minimal() + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12), 
            plot.title = element_text(size = 18),
            legend.position = "none") +
      labs(title = paste(type, " log2 Intensity: \nAfter ", toupper(norm)," Normalization", sep = ''))
    
    # DENSITY PLOT
    plot4 <- ggplot(df, aes(x = Intensity, colour = Sample)) + 
      geom_density()+ 
      theme_bw() + 
      theme(plot.title = element_text(size = 18),
            legend.position = "none") +
      labs(title = paste(type, " log2 Intensity: \nAfter ", toupper(norm)," Normalization", sep = ''));
    
    # QQ PLOT
    plot6 <- ggpubr::ggqqplot(eset_matrix_norm[, QQ_column]) +
      theme_bw() + 
      theme(plot.title = element_text(size = 18),
            legend.position = "none") +
      labs(title = paste(type, " ", qq_column, " QQ Plot ", toupper(norm)," Post-Normalization", sep = ''));
  }
  
  # ARRANGE PRE- AND POST- NORM PLOTS SIDE BY SIDE FOR COMPARISON
  # OR GENERATE MA/RLE PLOT
  if (norm != "none") { 
    if (plottype == "Boxplot"){
      plot <- gridExtra::grid.arrange(plot1 + labs(title = "Raw Data"), 
                                      plot3 + labs(title = "Normalized"),
                                      top = grid::textGrob(paste(type, ":", " ", toupper(norm), sep = ""), gp = grid::gpar(fontsize = 20)), 
                                      ncol = 2);
      
    } else if (plottype == "Density") {
      plot <- gridExtra::grid.arrange(plot2 + labs(title = "Raw Data"), 
                                      plot4 + labs(title = "Normalized"),
                                      top = grid::textGrob(paste(type, ":", " ", toupper(norm), sep = ""), gp = grid::gpar(fontsize = 20)), 
                                      ncol = 2);
      
    } else if (plottype == "MA") {
      plot <- limma::plotMA(eset_matrix_norm, array = ma_array, main = paste(type, toupper(norm), "Normalized MA Plot", ma_array, sep = " "))
      
    } else if (plottype == "RLE") {
        mn <- apply(eset_matrix_norm, 1, median);
        rle <- data.frame(sweep(eset_matrix_norm, MARGIN = 1, STATS = mn, FUN = '-'));
        
        graphics::par(mar = c(10, 4, 4, 2) + 0.1)
        
        plot <- graphics::boxplot(rle, main = "RLE (Relative Log Expression)\nShould be centered on 0 (blue line)",
                          names = colnames(eset_matrix_norm), las = 2)
        lines(x = c(0, ncol(eset_matrix_norm) + 1), y = rep(0, 2), col = "blue", lty = 2)
        lines(x = c(0,ncol(eset_matrix_norm) + 1), y = rep(0.1, 2), col = "red", lty = 2)
        
        graphics::par(mar = c(5, 4, 4, 2) + 0.1)
        
    } else if (plottype == "QQ"){
      plot <- gridExtra::grid.arrange(plot5 + labs(title = "Raw Data"), 
                                      plot6 + labs(title = "Normalized"),
                                      top = grid::textGrob(paste(type, " ", qq_column, ":", " ", toupper(norm), sep = ""), gp = grid::gpar(fontsize = 20)), 
                                      ncol = 2);
    }
    return(plot);
    
  } else {
    
    # IF NO NORMALIZATION SELECTED, RETURN ONLY PRE-NORM PLOTS
    if (plottype == "Boxplot") {
      plot <- plot1 + labs(title = paste(type, ": No Normalization", sep = ""))
      
    } else if (plottype == "Density") {
      plot <- plot2 + labs(title = paste(type, ": No Normalization", sep = ""))
      
    } else if (plottype == "MA") {
      for (i in 1:ncol(eset_postnorm)) {
        plot <- limma::plotMA(eset_matrix, array = i, main = paste(type, "MA Plot", i, sep = " "))
      }
    }
    
    return(plot); 
  }
  
}

