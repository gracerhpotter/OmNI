#-------------------------------------------------------------------------------
#' @title Draw Heatmap
#'
#' @description This function generates heatmap plots based on data expression set, 
#' type, annotation file, and color input. Additional optional inputs allow subsetting
#' data by gene name, highest row variance, and feature/column as well as clustering 
#' columns, k-clustering rows, taking row z-score, and taking log2 intensity.
#' 
#' @param eset ExpressionSet object
#' @param type string input - data input type
#' @param title_add string input - optional plot title customization
#' @param mapcolor string input - RColorBrewer palette
#' @param annot formatted annotation file input
#' @param subset string input - what kind of subsetting to do, options: none, subset_genes, subset_features, subset_variable 
#' @param subset_genes list of strings input - list of gene names generated from shiny selectInput
#' @param subset_genes_file list of strings input - list of gene names generated from gene names input file
#' @param show_row_names logical input - whether to show row names in heatmap
#' @param cluster_samples logical input - whether to cluster samples (columns) of heatmap
#' @param kclustering logical input - whether to k-cluster rows of heatmap
#' @param subset_features list of strings input - list of features (columns) to include in the heatmap
#' @param subset_variable numeric input - number of top variable genes to include in the heatmap
#' @param log2 logical input - whether to take log2 of all rows
#' @param zscore logical input - whether to take zscore of all rows
#' 
#' SPECIFIC TO DIFFERENTIAL HEATMAP USING LIMMA MODEL TOPTABLE OUTPUT
#' @param differential_heatmap logical input - whether heatmap is differential (input from limma model)
#' @param top_table table input - topTable object from limma linear modeling
#' @param dif_title string input - optional alteration to title of differential heatmap
#' @param dif_number_rows numeric input - number of topTable rows to include in differential heatmap
#'
#' @return ComplexHeatmap object
#' 
#' @examples
#' 
#' @import Biobase
#' @import ComplexHeatmap
#' 

drawHeatmaps <- function(eset, 
                         type, 
                         title_add = '', 
                         mapcolor,
                         group_color,
                         annot = annotation, 
                         subset = 'none',
                         subset_genes,
                         subset_genes_file,
                         show_row_names = NULL, 
                         cluster_samples = TRUE, 
                         kclustering = FALSE,
                         subset_features, 
                         subset_variable = 50, 
                         log2 = FALSE,
                         zscore = TRUE,
                         differential_heatmap = FALSE,
                         top_table,
                         dif_title,
                         dif_number_rows = 50){
  
  # FORMAT ESET
  eset <- eset[,order(pData(eset)$Group)]
  
  # GROUP COLORS
  annotLab <- data.frame(Group = factor(pData(eset)$Group, 
                                        levels = unique(pData(eset)$Group)));
  annotCol <- list(Group = grDevices::hcl.colors(length(levels(factor(pData(eset)$Group))), group_color))
  
  names(annotCol$Group) <- levels(factor(annotLab$Group))

  # TOP_ANNOTATION OBJECT
  ComplexHeatmap::ht_opt("message" = FALSE)
  ha_column <- ComplexHeatmap::HeatmapAnnotation(df = annotLab, 
                                                 col = annotCol)
  
  # SET MAPCOLOR
  if(mapcolor == "viridis"){mapcolor <- viridisLite::viridis(11); maponeway <- viridisLite::viridis(11);
  
  } else {maponeway <- rev(RColorBrewer::brewer.pal(11, mapcolor));
          mapcolor <- (rev(RColorBrewer::brewer.pal(11, mapcolor)))}

  # DEFINE BIG HEATMAP
  # allows for showing rownames only if heatmap is not too big
  isBigMap = function(emat){
    ncol(emat) > 50 || nrow(emat) > 75
  }

  # DEFINE SCALING FUNCTION
  # for z-score across rows
  safe_scale = function(m){
    stats::na.omit(t(scale(t(m))))
  }
  
  # GETHEATMAP FUNCTION
  getHeatmap = function(`matrix`, 
                        name = "Z-score", 
                        show_row_names. = show_row_names,
                        column_title,
                        cluster_columns = cluster_samples, 
                        split = NULL, 
                        title_names = "", 
                        col = mapcolor) {
    ComplexHeatmap::Heatmap (`matrix` = `matrix`, 
                             col = col, 
                             name = name, 
                             top_annotation = ha_column,
                             show_row_names = if(is.null(show_row_names.)) !isBigMap(`matrix`) else show_row_names.,
                             row_labels = substring(rownames(`matrix`), 0, 50),
                             cluster_columns = cluster_columns, 
                             use_raster = isBigMap(`matrix`), 
                             split = split, 
                             row_names_gp = grid::gpar(fontsize = 6),
                             column_title = paste(type, ": ", title_names, "\n", column_title, sep = '') )
  }

################################################################################
  
  # TITLE
  # Customize title based on contents of heatmap
  ctitle = ""
  
  if (differential_heatmap == TRUE) {ctitle = paste(dif_title, " Differential", sep = "")}
  
  if (subset == "none") {ctitle = paste(ctitle, "All features", sep = "")
  } else if (subset == "subset_genes") {ctitle = paste(ctitle, "Subset genes", sep = "")
  } else if (subset == "subset_features") {ctitle = paste(ctitle, "Subset features", sep = "")
  } else if (subset == "subset_variable") {ctitle = paste(ctitle, "Highest variation", sep = "")}
  
  if (kclustering == TRUE) {ctitle = paste(ctitle, ", K-clustered", sep = "")}
  if (log2 == TRUE) {ctitle = paste(ctitle, ", Log2 value", sep = "")
  } else {ctitle = paste(ctitle, ", Z-score", sep = "")}
  
################################################################################  
  
  # INITIALIZE EMAT_SEL
  emat_sel <- Biobase::exprs(eset)
 
  #----------------------------------------------------------------------------- 
  
  # STANDARD HEATMAP
  # no subsetting, optional z-score across rows
  if(subset == "none"){
  
    if(zscore == FALSE){
      name = "Range"
    
    } else{
      emat_sel <- safe_scale(emat_sel) # Z-score across rows
      emat_sel[emat_sel < -2] <- -2
      emat_sel[emat_sel > 2] <- 2
      name = "Z-Score"
    }
  }
  
  #-----------------------------------------------------------------------------
  
  # SUBSET GENES
  # subsetting heatmap by specific gene names (rows), either from selectInput in shiny or
  # a file input with column of gene names
  else if(subset == "subset_genes"){
    
    if(length(subset_genes) == 0){
      subset_genes <- subset_genes_file
    } 
    
    if (length(subset_genes) == 0) {stop("Please input gene names in the dropdown or file upload.")}
    
    exprs_eset <- emat_sel
    emat_sel <- matrix(ncol = ncol(eset))
    
    for (k in 1:length(subset_genes)) {
      if (length(subset_genes[[k]]) < 1) {next;}
      
      emat_sel <- rbind(emat_sel, subset(exprs_eset, grepl(subset_genes[k], row.names(exprs_eset))))

      if (zscore == FALSE){
        name = "Range"
        
      } else{
        emat_sel <- safe_scale(emat_sel) # Z-score across rows
        emat_sel[emat_sel < -2] <- -2
        emat_sel[emat_sel > 2] <- 2
        name = "Z-Score"
      }
    }
  }
    
  #-----------------------------------------------------------------------------
  
  # SUBSET FEATURES/GROUPS
  # subsetting heatmap by columns (groups) requires resetting the HeatmapAnnotation
  # to reflect number of columns being used
  else if(subset == "subset_features"){
    if (length(subset_features) == 0) {stop("Please select features to subset the heatmap.")}

    index = data.frame()

    for(k in 1:length(subset_features)){
      index <- rbind(index, which(na.omit(annot$SampleName) == subset_features[k]))
    }

    for(i in 1:nrow(annotLab)) {
      if(i %in% index[,1]){
      }
      else{
        annotLab[i,] <- NA
      }
    }

    annotLab <- na.omit(annotLab)
    
    ComplexHeatmap::ht_opt("message" = FALSE)
    ha_column <- ComplexHeatmap::HeatmapAnnotation(df = annotLab, col = annotCol)
    
    exprs_eset <- emat_sel
    emat_sel <- matrix(nrow = nrow(exprs_eset))
    
    for (k in 1:length(subset_features)) {
      if(length(subset_features[[k]]) < 1) next;
      
      emat_sel <- cbind(emat_sel, (exprs_eset[,colnames(exprs_eset) %in% subset_features[k], drop = F]))
    }
    
    emat_sel <- emat_sel[,-1]
    
    if (zscore == FALSE){
      name = "Range"
      
    } else{
      emat_sel <- safe_scale(emat_sel) # Z-score across rows
      emat_sel[emat_sel < -2] <- -2
      emat_sel[emat_sel > 2] <- 2
      name = "Z-Score"
    }
  }
    
  #-----------------------------------------------------------------------------
  
  # SUBSET HIGH VARIANCE ROWS
  # shiny input allows choice of how many of top variable rows to include
  else if(subset == "subset_variable"){
    exprs_eset <- emat_sel
    vars <- sort(apply(exprs_eset, 1, var), decreasing = TRUE)[1:subset_variable]
    emat_sel <- exprs_eset[names(vars),]
    
    if (zscore == FALSE){
      name = "Range"
      
    } else{
      emat_sel <- safe_scale(emat_sel) # Z-score across rows
      emat_sel[emat_sel < -2] <- -2
      emat_sel[emat_sel > 2] <- 2
      name = "Z-Score"
    }
  }
  
  #-----------------------------------------------------------------------------
  
  # IF DIFFERENTIAL, SUBSET BY TOP LOGFC
  # shiny input allows choice of how many of top rows to include
  else if(subset == "top_table_rows"){
    top <- top_table[order(abs(top_table$"logFC"), decreasing = TRUE),][1:dif_number_rows,]
    emat_sel <- exprs(eset[rownames(eset) %in% rownames(top)])
    
    if (zscore == FALSE){
      name = "Range"
      
    } else{
      emat_sel <- safe_scale(emat_sel) # Z-score across rows
      emat_sel[emat_sel < -2] <- -2
      emat_sel[emat_sel > 2] <- 2
      name = "Z-Score"
    }
  }
  
################################################################################
  
  # GENERATE HEATMAP
  # optional k-clustering of rows into 3 groups
  if(kclustering == TRUE){
    kclus <- stats::kmeans(emat_sel, 3);
    split <- paste0("Cluster ", kclus$cluster)
    
    htmp <- getHeatmap(matrix = emat_sel, 
                       cluster_columns = cluster_samples, 
                       name = name,
                       split = split, 
                       title_names = title_add, 
                       column_title = ctitle)
    
  } else if(kclustering == FALSE){
    htmp <- getHeatmap(matrix = emat_sel, 
                       cluster_columns = cluster_samples, 
                       name = name,
                       title_names = title_add, 
                       column_title = ctitle)
  }
  
  return(htmp);
}
