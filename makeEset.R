# NCmisc::list.functions.in.file("makeEset.R", alphabetic = TRUE)

##-------------------------------------------------
#' @title Make Expression Set
#'
#' @description This function takes in a data frame that is a raw output from a 
#' proteomics search engine (MaxQuant) or other Omics data and processes into 
#' an ExpressionSet object. 
#'
#' @param data A data.frame that contains samples in rows, and intensity values and annotation data in columns
#' @param annotate An annotation data frame that specifies columns in data that are samples. See Omics_Notebook example for format.
#' @param type Type of data to be processed or a name for the Omics Data.
#' @param log_transform TRUE to log transform data
#' @param data_format Format of data. See Omics_Notebook example for options. 
#' @param uniprot_annotation Whether or not to query uniprot for annotation data. (Requires web access and time consuming).
#'
#' @return An Expression Set object
#' 
#' @examples
#' eset <- makeEset(proteinGroups, annot, type='phos', data_format = 'Protein.Groups..MQ.');
#' 
#' @export

makeEset <- function(data, 
                     annotate, 
                     type, 
                     log_transform = TRUE,
                     data_format, 
                     uniprot_annotation = FALSE) {
  
  # data <- read.delim("~/Documents/ShinyApp_DEV/example_data/example_proteingroups.txt", header = TRUE, sep = "\t")
  # annotate <- openxlsx::read.xlsx("~/Documents/ShinyApp_DEV/example_data/example_annotation.xlsx", 2, colNames = TRUE)
  # type <- "Proteomics"
  # data_format <- "ProteinGroups"
  # log_transform = TRUE
  # uniprot_annotation = FALSE
  
  # data <- read.delim("~/Downloads/example_preprocessed.txt", header = T, sep = "\t")
  # annotate <- openxlsx::read.xlsx("~/Downloads/preprocessed_annotation.xlsx", 2, colNames = TRUE)
  # type <- "PreFormattedData"
  # data_format <- "PreFormatted"
  
  if(data_format == "PreFormatted"){
    data.matrix <- data.frame(data)
    data.matrix[is.na(data.matrix)] <- 0 # NA values become zeros
    annotate[,"SampleName"] <- make.names(annotate[,"SampleName"]); #format sample names
    annotate[,type] <- make.names(annotate[,type] ); # data column headers
    samp_cols = annotate[ (annotate[,type] != "NA." & annotate[,"SampleName"] != "NA."), type]
    sample_column_names <- colnames(data[, which(colnames(data) %in% samp_cols)])
    rownames(data.matrix) <- data.matrix[,1]
    data.matrix <- data.matrix[,-1]
    colnames(data.matrix) <- annotate[which(annotate[,type] %in% sample_column_names), "SampleName"]
    feature_identifier <- data.frame(feature_identifier = data[,1])
    # make expression set object
    eset <- Biobase::ExpressionSet(assayData = as.matrix(data.matrix))
    Biobase::fData(eset) <- feature_identifier
    Biobase::pData(eset) <- cbind(annotate[ which(annotate[,type] %in% sample_column_names),],
                                  colnames(exprs(eset)));
    rownames(Biobase::pData(eset)) <- colnames(data.matrix)
    return(eset);
  }
  
  # For MaxQuant formats, remove potential contaminants and reverse hits
  if(("Protein.IDs" %in% colnames(data))){
    data <- data[!grepl("CON_", data$Protein.IDs),] # remove contaminants from MQ
    data <- data[!grepl("REV_", data$Protein.IDs),] # remove reverse from MQ
  }

  if(("Protein" %in% colnames(data))){
    data <- data[!grepl("contam_", data$Protein),] # remove contaminants from FragPipe
    data <- data[!grepl("CON_", data$Protein),] # remove contaminants from MQ Phospho
    data <- data[!grepl("REV_", data$Protein),] # remove reverse from MQ Phospho
  }
  
  # Get and format Sample Names from annotation
  annotate[,"SampleName"] <- make.names(annotate[,"SampleName"]); #format sample names
  annotate[,type] <- make.names(annotate[,type] ); # data column headers
  samp_cols = annotate[ (annotate[,type] != "NA." & annotate[,"SampleName"] != "NA."), type]
  missing_cols = setdiff(samp_cols, colnames(data))
  
  # Stop if annotation file referencing columns not present
  if(!grepl("PhosphoSites", data_format) && length(missing_cols) > 0) {
    stop(paste0("Annotation file requests columns not present in data file (", type, ", / ", 
                data_format, "): ", paste0(missing_cols, collapse = ", ")))
  }

  #-----------------------------------------------------------------------------
  
  #' @title Split Semi Function
  #'
  #' @description
  #' Copies portion of contents of one column to another column. Used to grab
  #' gene and protein IDs and put into 'Gene' and 'Protein' columns.
  #' 
  
  split_semi = function(data, dest, src){
    if(src %in% colnames(data)){
      data[,dest] = sub(';.*', '', data[,src])
    }
    
    else if(!(dest %in% colnames(data))){
      stop(paste0(data_format, "input file requires one of the following columns:", dest, src))
    }
    
    return(data);
  }

  #-----------------------------------------------------------------------------
  
  if(data_format == "ProteinGroups") {
    data.matrix <- as.matrix(data[, samp_cols]); # pull out sample intensity values based on annotation
    
    if ("Majority.protein.IDs" %in% colnames(data)){
      data = split_semi(data, "Protein", "Majority.protein.IDs") # parse out protein names in MQ
      data[,"Protein"] <- gsub(".*[|]([^.]+)[|].*", "\\1", data[,"Protein"])
      data[,"Gene"] <- gsub(".*[|]([^.]+)[_].*", "\\1", data[,"Protein"])
      
    } else {
      try({data = split_semi(data, "Protein", "Protein.ID")}) # parse out protein names in FP
      try({data = split_semi(data, "Protein", "Protein.Group")}) # parse out protein names in DIA-NN
      try({data = split_semi(data, "Protein", "Accession")}) # parse out protein names in PD
    }
    
  }
  
  if(data_format == 'Peptides') {
    data.matrix <- as.matrix(data[, samp_cols]); # pull out sample intensity values based on annotation
    try({data = split_semi(data, "Protein", "Proteins")}) # parse out protein names in MQ
    try({data = split_semi(data, "Protein", "Protein.ID")}) # parse out protein names in FP
    
  }
  
  if(data_format == "PhosphoSites") { # Get data for MaxQuant Sites format
    
    if("+" %in%  data[,"Diagnostic.peak"]){ 
      data <- data[data[,"Diagnostic.peak"] != '+',]; # remove diagnostic features
    } 
    
    data <- data[data[,"Localization.prob"] >= 0.70 ,]; #filter low probability features
    
    # pull out sample intensity values based on annotation and collapse multiple sites
    data1 <- data[, !(colnames(data) %in% samp_cols)];
    data1 <- data1[rep(rownames(data1), 3),];
    
    data2 = do.call("cbind",lapply(1:length(samp_cols), FUN = function(i){
      cols.keep <- grep(samp_cols[i], colnames(data), value = T);  
      reshape2::melt(data[, cols.keep], measure.vars = cols.keep,
                     variable.name = paste("Phospho.Site.", i - 1, sep = ''),
                     value.name = samp_cols[i]);
    }))
    
    data <- cbind(data1, data2)
    data.matrix <- as.matrix(data[, samp_cols]);
    data = split_semi(data, "Protein", "Protein") # parse out protein names
    
  }
  
  sample_column_names <- colnames(data[, which(colnames(data) %in% samp_cols)])
  
  if(data_format == "Generic"){
    data.matrix <- as.matrix(data[, sample_column_names]);
    
  }
  
  if(data_format == "MetabolitePos" | data_format == "MetaboliteNeg") { # Either OpenMS output or XCMS Online output
    data.matrix <- as.matrix(data[, samp_cols]);
    
    if("mzmed" %in% colnames(data)){data[,"mz"] <- data[,"mzmed"]} # if xmcs output, make mz column
    if("mz_cf" %in% colnames(data)){data[,"mz"] <- data[,"mz_cf"]} # if openms output, make mz column
    
    if("rtmed" %in% colnames(data)){data[,"rt"] <- data[,"rtmed"]} # if xmcs output, make mz column
    if("rt_cf" %in% colnames(data)){data[,"rt"] <- data[,"rt_cf"]} # if openms output, make mz column
    
    if("links" %in% colnames(data)){ # if there are CSIFingerID guesses, grab first one and list in column
      data[,"CSIFingerID_Pubmed"] <- str_match(data[,"links"], "PubChem:\\((.*?)[ )]")[,2]
      data[,"CSIFingerID_HMDB"] <- str_match(data[,"links"], "HMDB:\\((.*?)\\)")[,2]
      data[,"identifier"] <- data[,"CSIFingerID_HMDB"]
    }
    
  } 
  
  possible_formats <- c("MetabolitePos", "MetaboliteNeg", "Generic", "PhosphoSites", 'Peptides', "ProteinGroups")
  if (!(data_format %in% possible_formats)){
    stop(paste0("The provided data format is invalid. Please make sure that a valid data format was
                selected in the annotation file on the `Inputs` sheet. \n", "DATA FORMAT PROVIDED: '",
                data_format, "'", sep = ""),
         call. = FALSE)

  }
  
  #-----------------------------------------------------------------------------
  
  # log2 Intensity Values
  if(log_transform == TRUE){ 
    data.matrix <- log2(data.matrix + 1); 
  }

  merge_str_nonempty <- function(a, b){
    if(is.null(a)){return(b)}
    if(is.null(b)){return(a)}
    
    i <- is.na(a) | a == ""
    a[i] = b[i]
    
    return(a)
  }
  
  emptysub <- function(p, r, x){
    missing <- grep(p, x, invert = T)
    ret <- sub(p, r, x)
    ret[missing] <- ""
    
    return(ret)
  }
  
  #-----------------------------------------------------------------------------
  
  # UNIPROT ANNOTATION
  if("Protein" %in% colnames(data)){
    
    uprot = data$Protein
    
    if(uniprot_annotation == TRUE){ try({ # Add uniprot annotation
      
      udata <- getUniprotAnnotation(IDs = data[,"Protein"], 
                                    genes =! ("Gene.names" %in% colnames(data))) 
      missing <- udata$ENTRY == ""
      
      if(any(missing)){
        gn = Reduce(merge_str_nonempty,
                    list(data[missing, "Gene.names"],
                         emptysub(".*(gene|GN)=([^ ]*).*", "\\2", data[missing, "Fasta.headers"])))
        
        gn = paste0("gene=", gn)
        g_udata <- getUniprotAnnotation(IDs = gn, genes =! ("Gene.names" %in% colnames(data)))
        udata[missing,] <- merge_str_nonempty(udata[missing,],g_udata)
      }
      
      data <- cbind(data, udata)
      uprot = data$ENTRY
      
    }, silent = FALSE)}
    
    data[, "Link"] <- paste("https://www.uniprot.org/uniprot/", uprot, sep = "") # Make uniprot hyperlink
    data[,"Uniprot"] <- paste("<a href='https://www.uniprot.org/uniprot/",uprot,"'>",
                              uprot,"</a>", sep = "")
  }
  
  #-----------------------------------------------------------------------------
  
  # Parse MQ to get gene names
  if ("Gene.names" %in% colnames(data)){ 
    data = split_semi(data, "Gene", "Gene.names")
    
  } else if ("Gene.Symbol" %in% colnames(data)){ # Parse PD to get gene names
    data = split_semi(data, "Gene", "Gene.Symbol")
    
  } else if ("Genes" %in% colnames(data)){ # Parse PD to get gene names
    data = split_semi(data, "Gene", "Genes")
    
  } else if("Fasta.headers" %in% colnames(data)){ # For species that have Uniprot but not MQ gene annotations
    uni = grepl("^(sp|tr)\\|",data[,"Fasta.headers"])
    
    if(any(uni)){
      data[,"Gene"] = ""
      data[uni,"Gene"] = sub("^..\\|.*\\|([^_]*).*","\\1",data[uni,"Fasta.headers"])
    }
  }
  
  #-----------------------------------------------------------------------------
  
  # make feature identifiers/rownames
  if("Gene" %in% colnames(data)){
    if("Protein" %in% colnames(data)){
      if(("Amino.acid" %in% colnames(data))&("Position"%in%colnames(data))){
        data[,"feature_identifier"] <- make.unique( paste(data[,"Gene"],"_", data[,"Protein"],"_",
                                                          data[,"Amino.acid"],"", data[,"Position"], sep="" ) );
      } else {
        data[,"feature_identifier"] <- make.unique( paste(data[,"Gene"], data[,"Protein"], sep="_" ) );
        
      }
    } else if("Transcript" %in% colnames(data)){
      data[,"feature_identifier"] <- make.unique( paste(data[,"Gene"], data[,"Transcript"], sep="_" ) );
      
    } else {
      data[,"feature_identifier"] <- make.unique( data[,"Gene"] );
      
    }
    
  } else {
    data[,"feature_identifier"] <- make.unique( make.names(data[,1]) ); # take the first column as feature ID if nothing else is found.
  }
  
  # TODO: Add better cleaning for metabolite data?
  # else if ("identifier" %in% colnames(data)) {
  #   ids <- sub(',.*', '', data$identifier)
  #   ids <- sub('.*c.{1}', '', ids)
  #   ids[ids == ""] <- "noID"
  #   data[,"feature_identifier"] <- make.unique(paste0(ids, "_", round(data[, "mz_cf"]), "_", round(data[, "rt_cf"])))
  # } 
  
  #-----------------------------------------------------------------------------
  
  data.matrix[is.na(data.matrix)] <- 0 # NA values become zeros
  
  # Make column and row names for data matrix
  colnames(data.matrix) <- annotate[which(annotate[,type] %in% sample_column_names), "SampleName"]
  rownames(data.matrix) <- data[,"feature_identifier"]
  
  # make expression set object
  eset <- Biobase::ExpressionSet(assayData = data.matrix)

  Biobase::fData(eset) <- data[, which(!(colnames(data) %in% sample_column_names)) ];

  Biobase::pData(eset) <- cbind(annotate[ which(annotate[,type] %in% sample_column_names),],
                       colnames(exprs(eset)) );

  rownames(Biobase::pData(eset)) <- colnames(data.matrix)
  
  # exprs(eset)
  return(eset);

}