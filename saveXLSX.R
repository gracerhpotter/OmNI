
# getwd()

generateXLSX <- function(eset, 
                         limmaFit = NULL, 
                         data_format, 
                         mapcolor = map_color, 
                         type,
                         coef_index = NULL, 
                         time_index = NULL, 
                         contrast_strings = NULL) {
  
  wbOut <- openxlsx::createWorkbook()

  if(type != "met_combined"){
    eSet = eset[,order(pData(eset)$Group)]
    writeDataToSheets(wb = wbOut, 
                      eset = eSet,
                      limmaFit = fit,
                      type = type, 
                      data_format = data_format);
  }

  
  output_filename = file.path(getwd(), paste(gsub("\\.", "", make.names(type)), "_Summary_", Sys.Date(), ".xlsx", sep = ""))
  openxlsx::saveWorkbook(wbOut, file = output_filename, overwrite = TRUE)
}

#-------------------------------------------------
#' Write data to sheets
#'
#' This function outputs summary data to an excel sheet to share with collaborators
#' 
#' @param wb an openxlsx workbook object
#' @param eset an ExpressionSet object with omics data
#' @param limmaFit a lmFit object or FALSE if no differential analysis
#' @param data_format format of data to be processed. See Omics_Notebook example for options.
#' @param mapcolor specifies color scale to use for heatmap: "viridis" "RdBu" "RdYlBu"
#' @param type Type of data to be processed, or a name for the Omics data
#' @param coef_index
#' @param time_index
#' @param contrast_strings
#'
#' @return wb with summary data formatted
#' 

writeDataToSheets <- function(wb, 
                              eset, 
                              limmaFit = NULL, 
                              data_format, 
                              mapcolor = map_color, 
                              type,
                              coef_index = NULL, 
                              time_index = NULL, 
                              contrast_strings = NULL) {
  
  # limmaFit = fit
  # type = "Proteomics"
  # data_format = "Protein.Groups..MQ."
  
  eset <- eset[,order(pData(eset)$Group)]
  
  # Make colors for sample names
  annotLab <- data.frame(Group = factor(pData(eset)$Group));
  annotCol <- list(Group = grDevices::rainbow(length(levels(factor(pData(eset)$Group)))) )
  sampleCols <- annotCol$Group[1:length(levels(factor(pData(eset)$Group)))][factor(pData(eset)$Group)];
  mapcolor <- rev(RColorBrewer::brewer.pal(7, "RdYlBu"))
  
  if (class(limmaFit) != "NULL"){
    # Create table for each coefficient
    if (class(coef_index) != "NULL") {
      DiffList <- vector("list", length(coef_index))
    } else {
      DiffList <- vector("list", (ncol(limmaFit$coefficients)))
      coef_index <- 1:(ncol(limmaFit$coefficients))
    }
    
    if (class(coef_index) != "NULL" & class(time_index) != "NULL") {
      if (time_index > 0) {
        DiffList_T <- limma::topTable(limmaFit, adjust.method = "BH", n = Inf, coef = coef_index)
      } else {
        DiffList_T <- NULL
      }
    }
    
    if (length(DiffList) > 1) {
      try({
        DiffList_F <- limma::topTable(limmaFit, adjust.method = "BH", n = Inf)
      })
    }
    
    for (i in 1:length(coef_index)) {
      DiffList[[i]] = limma::topTable(limmaFit, adjust.method = "BH", n = Inf, sort.by = 'p', 
                                      coef = coef_index[i])[,c("P.Value","adj.P.Val","logFC")]
    }
    
    # Match the row order to the first contrast
    eset <- eset[match(rownames(DiffList[[1]]),rownames(eset)),]
  } else {DiffList <- vector("list", 0 )}
  
  # Create new sheet to add to the workbook
  stName <- as.character(type)
  if(nchar(stName) > 31) {stName = substring(stName, 1, 31)}
  openxlsx::addWorksheet(wb = wb, sheetName = stName)
  
  # Format Uniprot hyperlinks
  if ("Link" %in% colnames(fData(eset))) {
    links <- fData(eset)$Link
    names(links) <- fData(eset)$Protein
    class(links) <- "hyperlink"
  }
  
  # Make row z-score values for "heatmap"
  emat_sel <- t(scale(t(exprs(eset)))) # Z-score across rows
  emat_sel[emat_sel < -2] <- -2
  emat_sel[emat_sel > 2] <- 2
  
  # Create the data table
  d = list()
  d$formatted_table <- ""
  
  if ("feature_identifier" %in% colnames(fData(eset))) {  
    d$formatted_table <- data.frame(fData(eset)[,c("feature_identifier")], emat_sel)
    d$names <- c("Feature", colnames(exprs(eset)))
  } else if ("Gene" %in% colnames(fData(eset))) {  
    d$formatted_table <- data.frame(fData(eset)[,c("Gene")], emat_sel)
    d$names <- c("Gene", colnames(exprs(eset)))
  } else { 
    d$formatted_table <- data.frame(rownames(exprs(eset)), emat_sel )
    d$names <- c("Feature", colnames(exprs(eset)))
  }
  
  d$col_widths <- c(20, rep(4, ncol(eset)));
  
  fillOutputTable <- function (outlist, fieldmap_v) {
    
    fieldmap = matrix(byrow = T, ncol = 3, data = fieldmap_v)
    valid_inds = fieldmap[,1] %in% colnames(fData(eset))
    
    if (any(valid_inds)) {
      fieldmap = matrix(fieldmap[valid_inds,], ncol = 3)
      outlist$formatted_table <- data.frame(outlist$formatted_table, fData(eset)[,fieldmap[,1]]);
      outlist$names <- c(outlist$names, fieldmap[,2])
      outlist$col_widths <- c(outlist$col_widths, as.integer(fieldmap[,3]));
    }
    outlist
  }
  
  grepadd <- function (outlist, regex, colsize) {
    
    cname <- grep(regex, colnames(fData(eset)), value = T)
    if (length(cname) > 0) {
      fillOutputTable(outlist, matrix(t(c(cname,cname,rep(colsize,length(cname)))),ncol = 3))
    } else {outlist}
  }
  
  d <- fillOutputTable(d, c(
    "Gene", "Gene", 16,
    "Protein.names", "Protein", 50,
    "Protein", "Uniprot", 16,
    "mz", "MZ", 16,
    "rt", "RT", 16,
    "Adduct", "Adduct", 8,
    "Formula", "Formula", 8,
    "Metabolite.name", "Metabolite.name", 8,
    "MS.MS.assigned", "MS.MS.assigned", 8,
    "logfc_Overall", "MaxFC", 8
  ))
  
  col_index <- length(d$names)
  
  if (class(limmaFit) == "MArrayLM") { 
    
    for (i in 1:length(DiffList)) {
      DiffList[[i]] <- DiffList[[i]][match(rownames(DiffList[[1]]),rownames(DiffList[[i]])),];
      d$formatted_table <- data.frame(d$formatted_table, DiffList[[i]][,c("P.Value","adj.P.Val","logFC")] );
      # d$names <- c(d$names, "P.Value","adj.P.Val","logFC");
      contrast <- colnames(limmaFit$coefficients)[i]
      d$names <- c(d$names, paste0(contrast, "_Pval"),paste0(contrast, "_adjPval"),paste0(contrast, "_logFC"));
      d$col_widths <- c(d$col_widths, 8,8,8);
    } 
    
    if (length(DiffList) > 1) {try({
      DiffList_F <- DiffList_F[match(rownames(DiffList[[1]]),rownames(DiffList_F)),];
      d$formatted_table <- data.frame(d$formatted_table, DiffList_F[,c("F","adj.P.Val")])
      d$names <- c(d$names, "F-Statistic", "F-Stat adj.P.Val");
      d$col_widths <- c(d$col_widths, 8,8);
    })}
    
    try({ if (class(DiffList_T) != "NULL") {
      DiffList_T <- DiffList_T[match(rownames(DiffList[[1]]),rownames(DiffList_T)),];
      d$formatted_table <- data.frame(d$formatted_table, DiffList_T[,"F"])
      d$names <- c(d$names, "F-Statistic:Time Trajectory");
      d$col_widths <- c(d$col_widths, 8);
    } }, silent = T)
    
  } else{
    d <- grepadd(d, "logfc_", 8)
  }
  
  d <- grepadd(d, "mummichogID_", 8)
  d <- fillOutputTable(d, c(
    "Fasta.headers", "Fasta_Headers", 8,
    "Uniprot_Function", "Uniprot_Function", 50,
    "Uniprot_Cellular_Location", "Uniprot_Cellular_Location", 50,
    "Uniprot_Disease", "Uniprot_Disease", 50,
    "GO_biological_process", "GO_biological_process", 50,
    "GO_molecular_function", "GO_molecular_function", 50,
    "GO_cellular_component", "GO_cellular_component", 50,
    "GO_ID", "GO_ID", 8,
    "ReactomeID", "ReactomeID", 8,
    "KEGG_ID", "KEGG_ID", 8,
    "identifier", "identifier", 16,
    "KEGG", "KEGG", 8,
    "Sequence", "Sequence", 16,
    "Proteins", "Proteins", 16
  ))
  
  # Add phospho site probabilities and data
  try({ 
    if (grepl("PhosphoSites", data_format)) { 
      if (all(c("Amino.acid","Position") %in% colnames(fData(eset)))) {
        fData(eset)[,"AAPos"] = paste(fData(eset)[,"Amino.acid"], fData(eset)[,"Position"], sep = "")
      }
      
      d <- fillOutputTable(d, c(
        "Localization.prob", "Localization Probability", 16,
        "Probabilities","Site Probabilities", 16,
        "AAPos", "Amino Acid", 16,
        "Sequence.window", "Peptide Sequence", 16))
    }
  })
  
  formatted_table <- data.frame(d$formatted_table, exprs(eset));  
  names <- c(d$names,colnames(eset));
  col_widths <- c(d$col_widths, rep(4, ncol(eset)));
  
  formatted_table <- data.frame(formatted_table, stringsAsFactors = FALSE)
  colnames(formatted_table) <-  make.unique(toupper(names));
  
  # Write data table to sheet
  openxlsx::writeDataTable(wb = wb, sheet = stName, x = formatted_table, xy = c("A",2), keepNA = FALSE, tableStyle = "TableStyleLight1")
  if ("Link" %in% colnames(fData(eset))) {openxlsx::writeData(wb = wb, sheet = stName, x = links, xy = c(which(names == "Uniprot"), 3))}
  
  # Add heatmap color
  openxlsx::conditionalFormatting(wb = wb, sheet = stName, type = "colourScale", cols = 2:(1 + ncol(eset)), 
                                  rows = 3:(2 + nrow(eset)), style = mapcolor[c(1, 4, 7)])
  
  # Add color bars for fold change
  if (length(DiffList) > 0) { 
    i = 1
    openxlsx::conditionalFormatting(wb = wb, sheet = stName, cols = (3 + col_index + (3 * (i - 1))), rows = 3:(3 + nrow(limmaFit)), 
                                    type = "databar", style = c("royalblue", "red"), gradient = FALSE, showvalue = FALSE)
  }
  
  # Rotate text for sample names
  for (ci in 1:ncol(eset)) {
    fgfill = if(is.na(sampleCols[ci])) {NULL} else {sampleCols[ci]}
    
    for (c_offset in c(1, (length(names) - ncol(eset)))) {
      openxlsx::addStyle(wb = wb, sheet = stName, style = openxlsx::createStyle(fgFill = fgfill, textRotation = 90, halign = "center", valign = "top"),
                         rows = 2, cols = ci + c_offset)
    }
  }
  
  # Merge cells and add Intensity title
  openxlsx::mergeCells(wb = wb, sheet = stName, rows = 1, cols = 2:(1 + ncol(eset)))
  openxlsx::writeData(wb = wb, sheet = stName, x = "Normalized Log2 Intensity, Row Z Score", xy = c(2,1))
  openxlsx::addStyle(wb = wb, sheet = stName, style = openxlsx::createStyle(textDecoration = "bold", halign = "center"), 
                     rows = 1, cols = 2, stack = TRUE)
  
  # Add final intensity title
  openxlsx::mergeCells(wb = wb, sheet = stName, rows = 1, cols = (length(names):(length(names) - ncol(eset) + 1)))
  openxlsx::writeData(wb = wb, sheet = stName, x = "Normalized Log2 Intensity", xy = c((length(names) - ncol(eset) + 1), 1))
  openxlsx::addStyle(wb = wb, sheet = stName, style = openxlsx::createStyle(textDecoration = "bold", halign = "center"), rows = 1, 
                     cols = (length(names) - ncol(eset) + 1), stack = TRUE)
  
  # Add contrast titles
  if (class(contrast_strings) == "NULL") {
    contrast_strings <- colnames(limmaFit$coefficients)
  }
  
  if (length(DiffList) > 0) {
    for (i in 1:length(DiffList)) {
      startCol <- col_index + 1 + (3 * (i - 1))
      openxlsx::mergeCells(wb = wb, sheet = stName, rows = 1, cols = startCol:(startCol + 2))
      openxlsx::writeData(wb = wb, sheet = stName, x = contrast_strings[i], xy = c(startCol, 1))
      openxlsx::addStyle(wb = wb, sheet = stName, style = openxlsx::createStyle(textDecoration = "bold", halign = "center"), 
                         rows = 1, cols = startCol, stack = TRUE)
    } 
  }
  
  # Freeze columns/rows and bold column titles
  openxlsx::freezePane(wb = wb, sheet = stName, firstActiveRow = 3, firstActiveCol = 2)
  openxlsx::addStyle(wb = wb, sheet = stName, style = openxlsx::createStyle(textDecoration = "bold", fontColour = 'black'),
                     rows = 1:2, cols = 1:length(names), gridExpand = TRUE, stack = TRUE)   
  
  # Don't treat gene names as dates
  openxlsx::addStyle(wb = wb, sheet = stName, style = openxlsx::createStyle(numFmt = "TEXT"), 
                     rows = 3:(2 + nrow(eset)), cols = 1, gridExpand = TRUE)
  
  # Set row heights and column widths
  openxlsx::setRowHeights(wb = wb, sheet = stName, rows = 2, heights = 100)
  openxlsx::setColWidths(wb = wb, sheet = stName, cols = 1:ncol(formatted_table), widths = col_widths)
}

