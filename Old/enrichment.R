
safe_scale = function(m){
  stats::na.omit(t(scale(t(m))))
}

# geneList <- calculateGeneList(fit, coef = "Glucose-NoGlucose")
calculateGeneList <- function(fit,
                              coef){
  
  top_sum <- limma::topTable(fit, adjust.method = "BH", n = Inf, sort.by = 'p', coef = coef);
  ranked <- cbind(top_sum[,'Gene'], top_sum[,'logFC'])
  
  colnames(ranked)<-c("GeneName", "rank")
  ranked <- ranked[ranked[,"GeneName"]!="",]
  ranked <- ranked[order(abs(as.numeric(ranked[,"rank"])), decreasing=TRUE),]
  ranked <- ranked[!duplicated(ranked[,"GeneName"]),]
  ranked <- ranked[order(as.numeric(ranked[,"rank"]), decreasing=TRUE),]
  
  top_sum$fcSign = sign(top_sum$logFC)
  top_sum$logP = -log10(top_sum$P.Value)
  top_sum$metric = top_sum$logP / top_sum$fcSign
  ranked <- top_sum[,c("Gene", "metric")]
  
  # make gene list
  geneList <- ranked[,2]
  names(geneList) <- as.character(ranked[,1])
  geneList = sort(geneList, decreasing = TRUE)
  
  return(geneList)
  
}

#' @title clusterProfiler Enrichment
#'
#' @description
#' A short description...
#' 
#' @param fit description
#' @param coef description
#' @param gmt description
#' @param enrichment description
#' @param pval_cutoff description
#'

# setwd("~/ShinyApp_V17")
# enriched_GSEA <- clusterProfilerEnrichment(geneList, gmt = "HUMAN_KEGG", enrichment = "GSEA")
# enriched_ORA <- clusterProfilerEnrichment(geneList, gmt = "HUMAN_KEGG", enrichment = "enricher")
clusterProfilerEnrichment <- function(geneList,
                                      gmt = "",
                                      enrichment = "",
                                      pval_cutoff = 0.05){
  
  # selecting gmt file
  # setwd("/Users/pottegra/Documents/OmicsAppTesting/Enrichment")
  # gmt = "HUMAN_KEGG"
  
  GMT_file <- list.files("./GMTs", pattern = gmt)
  my_geneset = readr::read_delim(paste0("./GMTs/", GMT_file))

  
  my_geneset = dplyr::mutate(my_geneset, term = pathway)
  length(unique(my_geneset$pathway))
  term2features = my_geneset %>% dplyr::select(term, feature_ids)
  term2pathway = my_geneset %>% dplyr::select(term, pathway)
  rm(my_geneset)
  
  if (enrichment == 'enricher'){
    enriched <- clusterProfiler::enricher(gene = names(geneList),
                                          universe = term2features$feature_ids,
                                          pAdjustMethod = "fdr",
                                          pvalueCutoff = pval_cutoff,
                                          minGSSize = 5,
                                          TERM2GENE = term2features,
                                          TERM2NAME = term2pathway )
    
  } else if (enrichment == 'GSEA'){
    enriched <- GSEA(gene = geneList,
                     pAdjustMethod = "fdr",
                     pvalueCutoff = pval_cutoff,
                     minGSSize = 5,
                     eps = 0,
                     TERM2GENE = term2features,
                     TERM2NAME = term2pathway)
    
  } 
  
  return(enriched)
  
}

#' @title clusterProfiler Enrichment Results Table
#'
#' @description
#' A short description...
#' 
#' @param enriched description
#' @param enrichment 
#'

enrichedTable <- function(enriched, 
                          enrichment){
  if (enrichment == 'enricher'){
    enriched_table = enriched@result %>%
      as.data.frame() %>%
      tibble::remove_rownames() %>%
      dplyr::select(Description, pvalue, p.adjust, qvalue, Count, GeneRatio, intersection = geneID) %>%
      dplyr::mutate(intersection = stringr::str_replace_all(intersection, "/", ", "))
    
  } else {
    enriched_table = enriched@result %>%
      as.data.frame() %>%
      tibble::remove_rownames() %>%
      dplyr::select(Description, pvalue, p.adjust, qvalue, setSize, NES, intersection = core_enrichment) %>%
      dplyr::mutate(intersection = stringr::str_replace_all(intersection, "/", ", "))
  }
  
  return(enriched_table)
}
  
#' @title clusterProfiler Enrichment Plots
#'
#' @description
#' A short description...
#' 
#' @param enriched description
#' @param gmt description
#' @param plottype description
#' @param dotplot_categories description
#' @param snetplot_categories description
#' @param cnetplot_layout description
#' @param snetplot_label description
#'

enrichedPlots <- function(enriched,
                          enrichment = "",
                          gmt,
                          geneList,
                          contrast = "",
                          plottype = "dot",
                          dotplot_categories = 15,
                          dotplot_title = "",
                          cnetplot_categories = 5,
                          cnetplot_layout = "kk",
                          cnetplot_label = "category",
                          cnet_title = "",
                          emap_categories = 30,
                          emap_title = "",
                          emap_category_label = 1,
                          upset_categories = 10,
                          upsetplot_title = "",
                          tree_categories = 30,
                          tree_clusters = 5,
                          treeplot_title = "",
                          heatmap_color,
                          heatmap_title = "") {
  
  enriched_pairwise <- enrichplot::pairwise_termsim(enriched)
  
  if (plottype == "go"){
    if (grepl("GO", gmt)){
      plot <- goplot(enriched)
    }
    
  } else if (plottype == "dot"){
    plot <- dotplot(enriched, 
                    showCategory = dotplot_categories, 
                    font.size = 12) +
      ggtitle(paste0(contrast, " ", dotplot_title))
    
  } else if (plottype == "cnet"){
    plot <- cnetplot(enriched, 
                     categorySize = "pvalue", # pvalue or geneNum
                     color.params = list(foldChange = geneList, edge = TRUE, category = "grey"), 
                     node_label = cnetplot_label, # category, gene, all, or none
                     showCategory = cnetplot_categories,
                     cex.params = list(category_label = 1, gene_label = 0.7),
                     circular = FALSE,
                     layout = cnetplot_layout) +
      ggtitle(paste0(contrast, " ", cnet_title))
    
  } else if (plottype == "heat"){
    # plot <- heatplot(enriched, 
    #                  foldChange = geneList, 
    #                  showCategory = 30,
    #                  symbol = 'rect', # rect or dot
    #                  label_format = 100)

    colNames <- unique(names(geneList))
    
    top_enriched = enriched@result %>%
      as.data.frame() %>%
      tibble::remove_rownames() %>%
      dplyr::select(Description, p.adjust)
    
    if (length(top_enriched$Description) > 50){
      rowNames <- top_enriched$Description[1:50]
    } else {
      rowNames <- top_enriched$Description
    }
    
    heatmapDF <- data.frame(matrix(ncol = length(colNames), nrow = length(rowNames)))
    colnames(heatmapDF) <- colNames
    rownames(heatmapDF) <- rowNames
    
    for (row in 1:nrow(heatmapDF)){
      for (col in 1:ncol(heatmapDF)){
        pathway <- rownames(heatmapDF)[row]
        gene <- colnames(heatmapDF)[col]
        
        if (gene %in% enriched@geneSets[[pathway]]){
          heatmapDF[row, col] <- geneList[[gene]]
        }
      }
      print(row)
    }
    
    heatmapDF[is.na(heatmapDF)] <- 0
    heatmapDF <- heatmapDF[,colSums(heatmapDF != 0) > 0]
    
    heatmapDF <- safe_scale(heatmapDF) # Z-score across rows
    heatmapDF[heatmapDF < -2] <- -2
    heatmapDF[heatmapDF > 2] <- 2
    
    plot <- heatmaply(heatmapDF,
                      scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "royalblue", 
                                                                              high = "red", 
                                                                              midpoint = 0, 
                                                                              limits = c(-max(abs(heatmapDF)), max(abs(heatmapDF)))),
                      show_dendrogram = c(FALSE, FALSE), 
                      showticklabels = c(FALSE, TRUE), 
                      xlab = "Genes",
                      ylab = "Pathways",
                      main = paste(contrast, "Enrichment Heatmap", sep = " "),
                      height = 1000,
                      fontsize_row = 7)
    
  } else if (plottype == "emap"){
    plot <- emapplot(enriched_pairwise,
                     cex.params = list(category_label = emap_category_label),
                     showCategory = emap_categories,
                     repel = TRUE) +
      ggtitle(paste0(contrast, " ", emap_title))
    
  } else if (plottype == "tree"){
    plot <- treeplot(enriched_pairwise,
                     showCategory = tree_categories,
                     nCluster = tree_clusters) +
      ggtitle(paste0(contrast, " ", treeplot_title))
    
  } else if (plottype == "upset"){
    plot <- upsetplot(enriched, 
                      n = upset_categories) +
      ggtitle(paste0(contrast, " ", upsetplot_title))
    
  } else if (plottype == "ridge"){
    # must be GSEA
    plot <- ridgeplot(enriched,
              label_format = 100)
    
  } else if (plottype == "gseamulti"){
    # must be GSEA
    plot <- gseaplot2(enriched, 
              geneSetID = 1:5, 
              base_size = 11, 
              rel_heights = c(2.0, 0.5, 0.5))
    
  } else if (plottype == "gseasingle"){
    # must be GSEA
    for (i in 1:5){
      print(gseaplot2(enriched, 
                      geneSetID = i, 
                      title = enriched$Description[i]))
    }
  }
  
  rm(enriched_pairwise)
  return(plot)
}


