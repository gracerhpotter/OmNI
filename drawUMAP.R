#-------------------------------------------------------------------------------
#' @title Draw UMAP
#'
#' @description This function generates UMAP plots given data expression set, type,
#' and RcolorBrewer color.
#' 
#' @param eset ExpressionSet object with omcis data
#' @param type string input - Type of data to be processed, or name for the Omics set
#' @param color string input - RcolorBrewer palette
#' @param title_add string input - optional add title
#' @param add_labels logical input - whether to add labels to points
#' @param shapes logical input - whether to make different groups different shapes
#'
#' @return UMAP ggplot object
#' 
#' @import ggplot2
#' @import Biobase
#' 

drawUMAP <- function(eset, 
                     type, 
                     color, 
                     title_add = "",
                     add_labels = FALSE,
                     shapes = FALSE) {
  
  # PREP DATA ------------------------------------------------------------------
  data <- t(Biobase::exprs(eset))
  data.labels <- Biobase::pData(eset)$Group
  
  use_size <- 15
  if (length(data.labels) < 15){
    use_size <- (length(data.labels) / 2)
  }
  
  data.umap <- uwot::umap(data, n_neighbors = use_size) 
  
  pdata <- Biobase::pData(eset)
  groups <- unique(pdata$Group)
  
  # MAKE UMAP PLOT -------------------------------------------------------------
  umap_graph <- ggplot(data = as.data.frame(data.umap), aes(x = data.umap[,1], y = data.umap[,2])) +
    
    scale_x_continuous(expand = expansion(mult = c(0.1, 0.2))) +
    scale_y_continuous(expand = expansion(mult = c(0.1))) +
    
    # ADD LABELS TO POINTS
    {if(add_labels == TRUE){
      ggrepel::geom_text_repel(label = pdata$SampleName, max.overlaps = 100);}} +
    
    # CHANGE SHAPE BY GROUP
    {if(length(groups) > 1 && shapes == TRUE) { 
      geom_point(size = 3, aes(shape = pdata$Group, colour = pdata$Group));}} +
    
    {if(length(groups) > 1 && shapes == TRUE) { 
      scale_shape_manual(values = 1:length(groups))} 
      else {geom_point(size = 3, aes(colour = pdata$Group));}} +
    
    # APPLY THEME
    theme_bw() + theme(legend.title = element_blank()) + 
    
    # CHANGE COLOR IF TOO MANY GROUPS
    {if(length(groups) > 12){
      scale_fill_discrete()}
      else {scale_color_brewer(palette = color)}} +
    
    # APPLY AXIS LABELS
    labs(title = paste("UMAP: ", title_add, "\n", type, sep = "")) + xlab("UMAP1") + ylab("UMPA2") ;
  
  return(umap_graph);
}
