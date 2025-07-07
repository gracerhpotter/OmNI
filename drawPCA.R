#-------------------------------------------------------------------------------
#' @title Draw PCA
#'  
#' @description This function generates PCA plots based on data expression set 
#' and type. Requires shiny input of color scheme with optional inputs re: shape,
#' including density graphs, and customizing title.
#' 
#' @param eset an ExpressionSet object with omcis data
#' @param x_axis string input - which PC to use as x-axis
#' @param y_axis string input -which PC to use as y-axis
#' @param type string input - type of data to be processed, or name for the Omics set.
#' @param color string input - RcolorBrewer palette
#' @param include_densities logical input - whether to include density graphs on the x & y axis
#' @param shapes logical input - whether to make different groups different shapes
#' @param title_add string input - optional customization of title
#' @param add_labels logical input - whether to add labels to points
#' @param add_ellipse logical input - whether to add ellipses around points by group
#'
#' @return PCA ggplot2 object
#' 

drawPCA <- function(eset, 
                    x_axis = "PC1", 
                    y_axis = "PC2", 
                    type, 
                    color,
                    include_densities = FALSE,
                    shapes = FALSE,
                    title_add = "",
                    add_labels = FALSE,
                    add_ellipse = FALSE) {
  
  assign("eset", eset, envir = .GlobalEnv)
  
  # PREP DATA W/ PCA ANALYSIS
  data = t(Biobase::exprs(eset))
  
  if(any(dim(data) == 0)){
    return(NULL)
  }
  
  PC_data <- stats::prcomp(data)
  
  percent_variance <- summary(PC_data)$importance["Proportion of Variance",] * 100
  
  pdata <- Biobase::pData(eset)
  groups <- unique(pdata$Group)
  
  # PCA PLOT -------------------------------------------------------------------
  # SET THEME
  basic_theme <- theme(legend.position = "none", 
                       axis.text = element_blank(), 
                       axis.line = element_blank(), 
                       axis.ticks = element_blank(),
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), 
                       panel.border = element_blank(), 
                       axis.title = element_text(size = 8))
  
  pca_graph <- ggplot(data = as.data.frame(PC_data$x), aes(x = PC_data$x[,x_axis], y = PC_data$x[,y_axis])) + 
    scale_x_continuous(expand = expansion(mult = c(0.1, 0.2))) +
    scale_y_continuous(expand = expansion(mult = c(0.1))) +
    
    # ADD ELLIPSES
    {if(add_ellipse == TRUE){
      ggforce::geom_mark_ellipse(aes(fill = pData(eset)$Group, color = pData(eset)$Group))}} +
      # stat_ellipse(geom = 'polygon', aes(colour = pData(eset)$Group, fill = pData(eset)$Group), alpha = 0.3)}} +
    
    # ADD LABELS
    {if(add_labels == TRUE){
      ggrepel::geom_text_repel(point.padding = 0.5, label = pData(eset)$SampleName, max.overlaps = 20);}} +
    
    # CHANGE SHAPES
    {if(length(groups) > 1 && shapes == TRUE){
      geom_point(size = 3, aes(shape = pData(eset)$Group, colour = pData(eset)$Group));}} +
    
    {if(length(groups) > 1 && shapes == TRUE){ 
      scale_shape_manual(values = 1:length(groups))}
      else {geom_point(size = 3, aes(colour = pData(eset)$Group));}} + 
    
    theme_bw() + theme(legend.title=element_blank()) + 
    
    # CHANGE COLOR IF TOO MANY GROUPS
    {if(length(groups) > 12){
      scale_fill_discrete()}
      else {scale_color_brewer(palette = color)}} +
    
    # ADD AXES LABELS
    labs(title = paste("PCA: ", title_add, "\n", type, sep = ""),
         x = paste(x_axis, sprintf(" (%2.0f%%)", percent_variance[x_axis]), sep = ""),
         y = paste(y_axis, sprintf(" (%2.0f%%)", percent_variance[y_axis]), sep = ""));
  
  # DENSITY GRAPHS ALONG AXES --------------------------------------------------
  if (include_densities == TRUE){  
    gg_dist_1 <- ggplot(data = as.data.frame(PC_data$x), aes(x = PC_data$x[,x_axis], fill = pData(eset)$Group)) + 
      geom_density(alpha = 0.4, linewidth = 0.2) + 
      ylab(paste(x_axis, "Density", sep = "\n")) + 
      theme_bw() + basic_theme + scale_fill_brewer(palette = color)
    
    gg_dist_2 <- ggplot(data = as.data.frame(PC_data$x), aes(x = PC_data$x[,y_axis], fill = pData(eset)$Group)) + 
      geom_density(alpha = 0.4, linewidth = 0.2) + 
      ylab(paste(y_axis, "Density", sep = "\n") ) + 
      theme_bw() + basic_theme + scale_fill_brewer(palette = color)
    
    piece1 <- gg_dist_1 + theme(axis.title.x = element_blank(), plot.margin = unit(c(0.5, -0.3, 0, 0.6), "cm"))
    piece2 <- gg_dist_2 + theme(axis.title.y = element_blank(), plot.margin = unit(c(-0.3, 1.0, 2.5, 0.2), "cm")) + coord_flip()
    piece3 <- pca_graph + theme(legend.position = "bottom", legend.title = element_blank(), plot.title = element_blank(),
                                plot.caption = element_text(size = 12, hjust = 0), plot.margin = unit(c(0,0,0.5,0.5), "cm")) + 
              labs(caption = paste("PCA: ", title_add, "\n", type, sep = ""))
    
    pca_graph <- cowplot::plot_grid( cowplot::plot_grid( piece1, piece3, ncol = 1, rel_heights = c(1,4)),
                                     cowplot::plot_grid(NULL, piece2, ncol = 1, rel_heights = c(1,4)),
                                     ncol = 2, rel_widths = c(4,1))
  }
  
  return(pca_graph)
  
  # 3D PCA? --------------------------------------------------------------------
  # plotly::plot_ly(as.data.frame(PC_data$x), x = ~PC1, y = ~PC2, z = ~PC3,
  #                 color = pData(eset)$Group,
  #                 colors = c('#636EFA','#EF553B', 'green'),
  #                 marker = list(size = 7)) %>%
  #   plotly::add_markers() %>%
  #   plotly::layout(
  #     title = "PCA"
  #     # scene = list(bgcolor = "#e5ecf6")
  #   )
}

#-------------------------------------------------------------------------------
#' @title PCA Percentages Plot
#' 
#' @description
#' A bar plot that contains the percentage of variance accounted for by each principle
#' component. Helps user decide which PCs to look at in PCA analysis.
#' 
#' @param eset Expression set object
#' @param type string input - dataset name from annotation file
#'
#' @return ggplot barplot object
#'

drawPCApercentages <- function(eset,
                               type) {
  data = t(Biobase::exprs(eset))
  
  if(any(dim(data) == 0)){
    return(NULL)
  }
  
  # CALCULATE PC VARIANCE
  PC_data <- stats::prcomp(data)
  percent_variance <- summary(PC_data)$importance["Proportion of Variance",] * 100
  PCs <- rownames(data.frame(percent_variance))
  PCdf <- data.frame(
    PC = c(PCs),
    Variance_Proportion = c(percent_variance)
  )
  
  # BARPLOT OF PC VARIANCE
  PCV_plot <- ggplot(PCdf) + 
                geom_bar(aes(x = reorder(PC, -Variance_Proportion), y = Variance_Proportion, fill = PC), stat = "identity") + 
                geom_line(aes(x = reorder(PC, -Variance_Proportion), y = Variance_Proportion, group = 1), stat = "identity", linewidth = 1) +
                theme(legend.position = "none", 
                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10), 
                      panel.background = element_rect(fill = "white", colour = "white"),
                      panel.grid.major = element_line(colour = "lightgrey")) +
                labs(title = paste(type, "PC Variance Percentages"),
                     x = "PC",
                     y = "Variance Proportion")
  
  return(PCV_plot);
}

#-------------------------------------------------------------------------------
#' @title PCA Loadings Plot
#' 
#' @description
#' A bar plot that contains the percentage of variance accounted for by each principle
#' component. Helps user decide which PCs to look at in PCA analysis.
#' 
#' @param eset Expression set object
#' @param column string input - which PC to use to generate loadings plot, i.e. "PC1"
#'
#' @return ggplot barplot object
#'

drawPCAloadings <- function(eset,
                            column) {
  
  data = t(Biobase::exprs(eset))
  col <- gsub("PC", "", column)
  col <- as.numeric(col)
  
  if(any(dim(data) == 0)){
    return(NULL)
  }
  
  # CALCULATE PC LOADINGS
  PC_data <- stats::prcomp(data)
  
  ls <- PC_data$rotation[,col]
  abs_ls <- abs(ls)
  ls_ranked <- sort(abs_ls, decreasing = TRUE)
  loading_scores <- data.frame(PC_data$rotation[names(ls_ranked[1:50]), col])
  colnames(loading_scores) <- c("Score")
  
  # PLOT PC LOADINGS
  PCL_plot <- ggplot(loading_scores) + 
                geom_bar(aes(x = reorder(rownames(loading_scores), -abs(Score)), y = loading_scores[,1], fill = loading_scores[,1] < 0), stat = "identity") +
                scale_fill_manual(guide = FALSE, breaks = c(TRUE, FALSE), values=c("royalblue", "red")) +
                ggtitle("Top 50 Loading Scores") +
                xlab("Row") + 
                ylab("Loading Score") + 
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10), 
                      panel.background = element_rect(fill = "white", colour = "white"),
                      panel.grid.major = element_line(colour = "lightgrey"))
  
  return(PCL_plot);
}

