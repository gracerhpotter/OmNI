#-------------------------------------------------------------------------------
#' @title Format Annotation File
#' 
#' @description
#' Takes the raw excel annotation file input and formats it into the format used by the
#' rest of the functions in Omics Notebook.
#' 
#' @param annot raw annotation file read in via readxlsx
#' @param group string input - which column is functioning as the group
#'
#' @return A formatted dataframe with annotation information in columns.
#'

formatAnnotation <- function(annot,
                             group) {
  
  # Make a group column that combines the chosen columns by row
  annot$Group <- do.call(paste0, c(annot[group]))
  
  # Make contrast groups based on similar entries in group column
  contrasts <- unique(gsub("\\.","", make.names(na.omit(annot$Group))))
  
  # Make sure Group names are valid
  if(length(annot$Group) == length(contrasts)){
  annot$Group = contrasts
  
  } else {
      annot[,"Group"] <- gsub("\\.", "", make.names(annot[,"Group"]))
      
  }
  
  # Make unique sample names
  annot[,"SampleName"] <- gsub("\\.", "", make.names(make.unique(as.character(annot[,"SampleName"]))))

  return(annot)
}

