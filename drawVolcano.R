#-------------------------------------------------------------------------------
#' @title Label Function
#' 
#' @description
#' Get gene names and/or feature identifiers to use to label points in graphs.
#' 
#' @param dat data table input
#' 
#' @returns Label list for graphing
#' 

getLabelVars <- function(dat){
  
  # IF GENE COLUMN LABEL BY GENE NAME
  if( "Gene" %in% colnames(dat)){
    textsize <- c(2, 2)
    labname <- rep("Gene", 2)
    
  # IF METABOLOMICS LABEL BY MX & FEATURE IDENTIFIER
  } else if( "mz" %in% colnames(dat)) {
    textsize <- c(2, 1.5)
    labname <- c("round(mz, digits=2)", "feature_identifier")
    
  # OTHERWISE JUST USE FEATURE IDENTIFIER COLUMN
  } else {
    textsize <- c(1, 1.5)
    labname <- rep("feature_identifier", 2)
    
  }
  
  list(textsize = textsize, labname = labname)
}

#-------------------------------------------------------------------------------
#' @title Volcano/MD Plot Base Function
#' 
#' @description
#' Function to create volcano and/or MD plots.
#' 
#' @param dat table input - data table used for plotting (may be topTable output or top logFCs)
#' @param xvar string input - which variable is the x-axis
#' @param yvar string input - which variable is the y-axis
#' @param title string input - optional title add
#' @param label_names list of strings - names of which rows to use for labelling 
#' @param colorby string input - what variable to apply color by
#' @param add_labels logical input - whether to add labels to the plot
#' @param label_specific list of strings input - list of specific genes to be labelled
#' @param up_color string input - hex code of color for up-regulated rows
#' @param down_color string input - hex code of color for down-regulated rows
#' @param caption string input - caption containing cutoff values
#' @param top_values numeric input - p-value cutoff for volcano plots dashed line
#' @param top_fc numeric input - logFC cutoff for volcano plots dashed line
#' @param cutoff numeric input - logFC cutoff for MD plots dashed line
#' 
#' @returns ggplot object
#' 

drawVPlots <- function(dat, 
                       xvar, 
                       yvar, 
                       title, 
                       label_names, 
                       colorby, 
                       add_labels, 
                       label_specific,
                       up_color,
                       down_color,
                       caption,
                       top_values,
                       top_fc,
                       cutoff){
  
  lv <- getLabelVars(dat)
  
  suppressWarnings({
    
    plot1 <- ggplot(data = data.frame(dat), aes(text = paste("Name: ", feature_identifier))) + 
      geom_point(aes_string(x = xvar, y = yvar, colour = colorby), size = 2, pch = 16, alpha = 0.8) +
      scale_colour_manual(values = c(NS = "grey", Up = up_color, Down = down_color)) +
      labs(title = title,
           caption = caption) + 
      theme_bw() + 
      theme(legend.title = element_blank(), plot.title = element_text(size = 10)) +
    
    {if (xvar == "Mean") {geom_hline(yintercept = c(-cutoff, cutoff), linetype = "dashed")}} +
    {if (xvar == "logFC") {geom_hline(yintercept = -log10(top_values), linetype = "dashed")}} +
    {if (xvar == "logFC") {geom_vline(xintercept = c(-top_fc, top_fc), linetype = "dashed")}}
      
    #---------------------------------------------------------------------------
    
    if (add_labels == FALSE && label_specific == FALSE){
      plot <- plot1
      
    } else if ((add_labels == TRUE || label_specific == TRUE)){
      plot1 <- plot1 + 
        ggrepel::geom_text_repel(data = data.frame(dat[label_names,]),
                                 size = 4,
                                 max.overlaps = 100,
                                 aes_string(x = xvar, y = yvar, label = lv$labname[1])) +
        geom_point(data = data.frame(dat[label_names,]), size = 2, pch = 21,
                   aes_string(x = xvar, y = yvar))
      
      plot <- plot1
      
    }
  })
}

#-------------------------------------------------------------------------------
#' Draw Volcano Plots
#'
#' @description
#' This function makes a volcano plot based on the output of the limma model fit.
#' 
#' @param dat data table input - the result of topTable()
#' @param type string input - type of data to be processed, or a name for the Omics analysis
#' @param v string input - variable for volcano y axis, either p-value or adjusted p-value
#' @param up_color string input - hex code for color for up-regulated rows
#' @param down_color string input - hex code for color for down-regulated rows
#' @param title_add string input - optional add title
#' @param add_labels logical input - whether to add labels to points
#' @param label_specific logical input - whether to add label to specific gene(s)
#' @param label_specific_gene string input - what specific gene to label
#' @param top_values numeric input - p-value cutoff for significance
#' @param top_fc numeric input - logFC value cutoff for significance
#' @param label_number numeric input - how many points to label if add_labels is TRUE
#' 
#' @import ggplot2
#' 
#' @returns ggplot object volcano plot

drawVolcano <- function(dat, 
                        type,
                        v,
                        up_color = "red",
                        down_color = "royalblue",
                        title_add = "",
                        title_type,
                        add_labels = FALSE,
                        label_specific = FALSE,
                        label_specific_gene = NULL,
                        top_values = 0.05, 
                        top_fc = 0,
                        label_number = 20){ 
  
  # SPECIFY DOT LABELS
  
  if(label_specific == TRUE && add_labels == FALSE){ # if only labeling one specific gene by name
    label_names = c()
    
    for(gene in label_specific_gene){
      label_names <- c(label_names, rownames(subset(dat, gene == gsub("_.*", "", dat$feature_identifier))))
    }
    
  } else if (label_specific == FALSE && add_labels == TRUE) { # if not labeling specific gene but yes labeling top/bottom logFC genes/rows
    label_names <- c(rownames(dat[dat$logFC > 0,][order(dat$P.Value, decreasing = FALSE),])[1:(label_number / 4 + label_number %% 4)],
                     rownames(dat[dat$logFC < 0,][order(dat$P.Value, decreasing = FALSE),])[1:(label_number / 4)],
                     rownames(dat[order(dat$logFC, decreasing = FALSE),])[1:(label_number / 4)],
                     rownames(dat[order(dat$logFC, decreasing = TRUE),])[1:(label_number / 4)])
    
    label_names <- unique(label_names)
    
  } else if (label_specific == TRUE && add_labels == TRUE){ # if labeling both specific gene by name AND top/bottom logFC genes/rows
    label_names = c()
    
    for(gene in label_specific_gene){
      label_names <- c(label_names, rownames(subset(dat, gene == gsub("_.*", "", dat$feature_identifier))))
    }
    
    label_names <- c(label_names, rownames(dat[dat$logFC > 0,][order(dat$P.Value, decreasing = FALSE),])[1:(label_number / 4 + label_number %% 4)],
                     rownames(dat[dat$logFC < 0,][order(dat$P.Value, decreasing = FALSE),])[1:(label_number / 4)],
                     rownames(dat[order(dat$logFC, decreasing = FALSE),])[1:(label_number / 4)],
                     rownames(dat[order(dat$logFC, decreasing = TRUE),])[1:(label_number / 4)])
    
    label_names <- unique(label_names) 
  }
  
  # DOT COLORING UP VS. DOWN VS. NS
  s <- paste0(v, ".sig")
  dat[,s] <- "NS"
  dat[,s] [( dat[,v] <= top_values & dat$logFC >= top_fc )] <- "Up"
  dat[,s] [( dat[,v] <= top_values & dat$logFC <= (-1 * top_fc))] <- "Down"
  dat[,s] <- factor(dat[,s], levels = c("NS", "Up", "Down"))
  
  # PLOT VOLCANO PLOT
  plot <- drawVPlots(dat, 
                     xvar = "logFC", 
                     yvar = paste0("-log10(",v,")"), 
                     title = paste(type, ": ", title_add, "\n", "Volcano Plot ", title_type, sep = ""),
                     caption = paste("FC cutoff: ", top_fc, " & P-Val cutoff: ", top_values, sep = ""),
                     label_names = label_names, 
                     colorby = paste0(v, ".sig"),
                     add_labels = add_labels,
                     label_specific = label_specific,
                     up_color = up_color,
                     down_color = down_color,
                     top_values = top_values,
                     top_fc = top_fc,
                     cutoff = NULL) 
  
  return(plot);
}

################################################################################

#' 
#' @title Fold Change & MD Plot
#'
#' @description
#' Calculate log fold change by an internal method, in case limma fails, and create
#' MD plot. For the case where the differential analysis fails (e.g., due to too few 
#' samples per group), and we need to calculate an average value for each group,
#' and then calculate a simple difference between all groups.
#' 
#' @param annot formatted annotation file generated table
#' @param eset ExpressionSet object
#' @param type string input - name of dataset from annotation file
#' @param fc_cutoff numeric input - fold change cutoff for significance
#' @param up_color string input - hex code for color of up-regulated rows
#' @param down_color string input - hex code for color of down-regulated rows
#' @param title_add string input - optional add title
#' @param label_number numeric input - how many points to label if add_labels is TRUE
#' @param add_labels logical input - whether to add labels to top abs(FC) points
#' @param label_specific logical input - whether to add label to specific gene(s)
#' @param label_specific_gene string input - which gene(s) to label
#' @param return_logfc_index logical input - whether to return a list of logfc index names to choose which to use for MD plot render
#' @param logfc_index_choice string input - which logfc name to use to render MD plot
#'
#' @return ggplot object MD plot
#' 

drawMD = function(annot,
                  eset,
                  type,
                  fc_cutoff = 1,
                  up_color = "red",
                  down_color = "royalblue",
                  title_add = "",
                  label_number = 20,
                  add_labels = FALSE,
                  label_specific = FALSE,
                  label_specific_gene = NULL,
                  return_logfc_index = FALSE,
                  logfc_index_choice = ""){
  
  contrastgroups = unique(annot$Group)
  logfc_index <- c();
  
  # CALCULATE MEAN VALUE PER GROUP
  for(j in 1:length(contrastgroups)){
    if(sum(pData(eset)$Group == contrastgroups[j]) > 1){
      fData(eset)[, paste("mean_", contrastgroups[j], sep = "")] <-
        rowMeans(exprs(eset)[, pData(eset)$Group == contrastgroups[j]])
      
    } else {
      fData(eset)[, paste("mean_", contrastgroups[j], sep = "")] <-
        exprs(eset)[, pData(eset)$Group == contrastgroups[j]] 
      
    }
  }
  
  # CREATE LOG FOLD CHANGE B/W EACH PAIRING OF GROUPS
  for(j in 1:(length(contrastgroups) - 1)){
    for(k in 2:(length(contrastgroups))){
      col_name <- paste("logfc_", contrastgroups[k], "_", contrastgroups[j], sep = "");
      fData(eset)[, col_name] <- (fData(eset)[, paste("mean_", contrastgroups[k], sep = "")] -
                                    fData(eset)[, paste("mean_", contrastgroups[j], sep = "")]) 
      col_name_rev <- paste("logfc_", contrastgroups[j], "_", contrastgroups[k], sep = "")
      fData(eset)[, col_name_rev] <- (fData(eset)[, paste("mean_", contrastgroups[j], sep = "")] -
                                    fData(eset)[, paste("mean_", contrastgroups[k], sep = "")]) 
      
      logfc_index <- c(logfc_index, col_name, col_name_rev)
    }
  }
  
  # We calculate a maximum fold change value across all groups as a way to see which features are changing most across all groups.
  # This is analogous to the F-statistic from the differential analysis.
  fData(eset)[, "logfc_Overall"] <- apply(fData(eset)[, grep("mean", colnames(fData(eset)))], 1, function(x) max(x) - min(x))
  
  logfc_index <- unique(logfc_index);
  
  if (return_logfc_index == TRUE){
    return(logfc_index);
    
  } else {
  
  ## PLOTTING
  
    logfc_index <- logfc_index_choice
    
    for(ci in 1:length(logfc_index) ){ try({
      # For the data set and coefficient, make the summary table and name
      top_sum <- data.frame(logFC = fData(eset)[, logfc_index[ci]],
                            Mean = rowMeans(exprs(eset)), 
                            feature_identifier = fData(eset)[, "feature_identifier"]);
      
      if("Gene" %in% colnames(fData(eset))){ top_sum$Gene <- fData(eset)$Gene
      } else if("mz" %in% colnames(fData(eset)) ){ top_sum$mz <- fData(eset)$mz }
      
      type_name <- paste(type, " ", logfc_index[ci], sep = "");
      
      if(fc_cutoff == 0){ cutoff <- 1 } else { cutoff <- fc_cutoff }
      
      # Gene labels
      label_number = label_number / 2
      
      if(label_specific == TRUE & add_labels == FALSE){
        label_names = c()
        for(gene in label_specific_gene){
          label_names <- c(label_names, rownames(subset(top_sum, grepl(gene, top_sum$feature_identifier))))
        }
        
      } else if(label_specific == FALSE & add_labels == TRUE) {
        label_names <- c(rownames(top_sum[order(top_sum$logFC, decreasing = FALSE),])[1:label_number],
                         rownames(top_sum[order(top_sum$logFC, decreasing = TRUE),])[1:label_number])
      
      } else if(label_specific == TRUE & add_labels == TRUE){
        label_names = c()
        
        for(gene in label_specific_gene){
          label_names <- c(label_names, 
                           rownames(subset(top_sum, grepl(gene, top_sum$feature_identifier))))
        }
        
        label_names <- c(label_names, 
                         rownames(top_sum[order(top_sum$logFC, decreasing = FALSE),])[1:label_number],
                         rownames(top_sum[order(top_sum$logFC, decreasing = TRUE),])[1:label_number])
      }
      
      # Dot coloring
      top_sum$Significance <- "NoChange"
      top_sum$Significance[top_sum$logFC > cutoff ] <- "Up"
      top_sum$Significance[(top_sum$logFC < (-1 * cutoff)) ] <- "Down"
      top_sum$Significance <- factor(top_sum$Significance, levels = c("NoChange", "Up", "Down"))
      
      # Plot
      plot <- drawVPlots(top_sum, 
                         xvar = "Mean", 
                         yvar = "logFC", 
                         title = paste("MD Plot: ", title_add, "\n", type_name, sep = ""),
                         caption = paste("FC cutoff: ", cutoff, sep = ""),
                         label_names = label_names, 
                         colorby = "Significance", 
                         add_labels = add_labels,
                         label_specific = label_specific,
                         up_color = up_color,
                         down_color = down_color,
                         top_values = NULL,
                         top_fc = NULL,
                         cutoff = cutoff) 
      
      return(plot);
    })}
  }
}

#-------------------------------------------------------------------------------
#' 
#' @title Glimma Interactive Volcano Plot
#'
#' @description
#' Generate the arguments for the Glimma::glXYPlot, which has an interactive volcano 
#' plot and Normalization Intensity plot.
#' 
#' @param eset Expression Set object
#' @param fit limma lmFit output
#' @param top_table dataframe input - limma topTable output
#' @param type string input - name of dataset from annotation file
#' @param coef string input - contrast used for volcano plot
#'
#' @return list of arguments to input to Glimma::glXYPlot function
#'

interactiveVolcano <- function(eset, 
                               fit = NULL, 
                               top_table, 
                               type, 
                               coef){
  
  dt = limma::decideTests(fit)
  
  annot_columns <- "feature_identifier";
  if("Gene" %in% colnames(fData(eset))) { annot_columns <- c(annot_columns, "Gene"); }
  if("Protein" %in% colnames(fData(eset))) { annot_columns <- c(annot_columns, "Protein"); }
  if("Protein.names" %in% colnames(fData(eset))) { annot_columns <- c(annot_columns, "Protein.names"); }
  if("Uniprot" %in% colnames(fData(eset))) { annot_columns <- c(annot_columns, "Uniprot"); }
  if("mz" %in% colnames(fData(eset))) { annot_columns <- c(annot_columns, "mz"); }
  if(length(annot_columns)==1) { annot_columns <- c(annot_columns, "feature_identifier"); }
  
  args = list(counts = exprs(eset), 
              groups = pData(eset)$Group,
              side.ylab = "Normalized Intensity", 
              anno = fData(eset)[,annot_columns],
              side.main = 'feature_identifier',
              main = coef,
              html = paste(coef, "_", type, "_GlimmaVolcanoPlot_", Sys.Date(), sep = ""),
              launch = TRUE)
  
  if(class(fit)!="NULL"){
    args2 = list(x = fit$coef[,coef], 
                 y = fit$lod[,coef], 
                 xlab = "logFC", 
                 ylab = "logOdds", 
                 status = dt[,coef])
  } else {
    args2 = list(x = top_table$Mean, 
                 y = top_table$logFC, 
                 xlab = "Mean", 
                 ylab = "log2FC")
  }
  
  args = c(args, args2)
  return(args);
}
