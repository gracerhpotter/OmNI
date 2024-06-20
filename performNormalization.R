#-------------------------------------------------------------------------------
#' @title Intensity Normalization
#'
#' @description
#'   This function takes omics data set and performs normalization and QC plotting 
#'   to inspect data. Can be expanded with additional normalization methods.
#' 
#' @param eset an ExpressionSet object with omics data
#' @param type string input - type of data to be processed
#' @param norm string input - normalization method ("none", "quantile", "loess", "median", "z transform")
#' @param zero_cutoff numeric input - filter out rows with more than this percentage of zeros/NAs, set to 0.3
#' @param min_feature numeric input - 
#' @param norm_by_batches logical input - whether to do batch normalization
#'
#' @return an expression set object, with normalized data
#' 
#' @examples
#'   eset <- intensityNorm(eset,
#'                         type = "Proteomics",
#'                         norm = "loess",
#'                         zero_cutoff = 0,
#'                         min_feature = 0.01,
#'                         norm_by_batches = F)
#' 
#' 


intensityNorm <- function(eset, 
                          type, 
                          norm,
                          zero_cutoff = 0.3,
                          min_feature = 0.01,
                          norm_by_batches = FALSE,
                          batch_column = '',
                          impute_missing = FALSE,
                          IRS_column = '',
                          outlier_cols = ''){
  
  message("** INITIALIZED NORMALIZATION")
  message(paste("-- DATA:", type))
  message(paste("-- TYPE:", norm))
  
  # CLEAN ESET DATA AND CREATE ESET MATRIX
  eset <- eset[apply(eset, 1, FUN = function(x){sum(x == 0)}) < (ncol(eset) * zero_cutoff),] # Filter out rows with >30% NAs
  eset <- eset[,colSums(exprs(eset) > 0) >= min_feature * nrow(exprs(eset))];
  
  for (i in outlier_cols){
    if(grepl(type, i)){
      column <- gsub(paste0(type, "_"), "", i)
      message(paste("-- REMOVED OUTLIER COLUMN:", column))
      column <- which(colnames(exprs(eset)) == column)
      eset <- eset[, -c(column)]
    }
  }
  
  eset_matrix <- Biobase::exprs(eset)
  
  eset_matrix_norm <- eset_matrix
  
  #-----------------------------------------------------------------------------
  # ACCOUNT FOR BATCHES IF APPLICABLE
  if (!norm_by_batches) {
    pData(eset)$Batch2 <- 1
    message("-- BATCH: FALSE")
    
  } else {
    pData(eset)$Batch2 <- pData(eset)[, batch_column]
    message("-- BATCH: TRUE")
    message(paste("---- Batches:", paste(unlist(unique(pData(eset)$Batch2)), collapse = ", ")))
    
  }
  
  #-----------------------------------------------------------------------------
  # PERFORM NORMALIZATION BASED ON INPUT OF TYPE
  for (i in 1:length(unique(pData(eset)$Batch2))) {
    index <- unique(pData(eset)$Batch2)[i]
    message("-- NORMALIZING")
    
    if (norm == 'quantile') {
      eset_matrix_norm[,grep(index,pData(eset)$Batch2)] <- as.matrix(preprocessCore::normalize.quantiles(eset_matrix[,grep(index,pData(eset)$Batch2)]))
      dimnames(eset_matrix_norm) <- dimnames(eset_matrix)
      
    } else if (norm == 'loess') {
      eset_matrix_norm[,grep(index,pData(eset)$Batch2)] <- as.matrix(limma::normalizeCyclicLoess(eset_matrix[,grep(index,pData(eset)$Batch2)], method = 'pairs'))
      
    } else if (norm == 'median') {
      eset_matrix[eset_matrix == 0] <- NA
      colMedians <- matrixStats::colMedians(eset_matrix, na.rm = T)
      meanColMedian <- mean(colMedians, na.rm = T)
      eset_matrix_norm <- matrix(nrow = nrow(eset_matrix), ncol = ncol(eset_matrix), byrow = T)
      
      normFunc <- function(colIndex){
        (eset_matrix[rowIndex, colIndex] / colMedians[colIndex]) * meanColMedian
        
      }
      
      for (rowIndex in seq_len(nrow(eset_matrix))) {
        eset_matrix_norm[rowIndex,] <- vapply(seq_len(ncol(eset_matrix)), normFunc, 0)
      }
      
      eset_matrix_norm[is.na(eset_matrix_norm)] <- 0
      colnames(eset_matrix_norm) <- colnames(eset_matrix)
      rownames(eset_matrix_norm) <- rownames(eset_matrix)
      
    } else if (norm == 'z transform') {
      eset_matrix_norm[,grep(index,pData(eset)$Batch2)] <- as.matrix(scale(eset_matrix[,grep(index,pData(eset)$Batch2)]));
      colnames(eset_matrix_norm) <- colnames(eset_matrix)
      rownames(eset_matrix_norm) <- rownames(eset_matrix)
      
    } else if (norm == "MAD") {
      eset_matrix_norm[,grep(index,pData(eset)$Batch2)] <- as.matrix(NormalyzerDE::performSMADNormalization(eset_matrix[,grep(index,pData(eset)$Batch2)], noLogTransform = TRUE))

      # eset_matrix[eset_matrix == 0] <- NA
      # # Calculate the median and MAD for each column
      # medians = apply(eset_matrix[,grep(index,pData(eset)$Batch2)], 2, median, na.rm = TRUE)
      # mad_values = apply(eset_matrix[,grep(index,pData(eset)$Batch2)], 2, mad, na.rm = TRUE)
      # 
      # # Perform medianMAD normalization
      # normalized_data = sweep(eset_matrix[,grep(index,pData(eset)$Batch2)], 2, medians, "-")
      # eset_matrix_norm[,grep(index,pData(eset)$Batch2)] = sweep(normalized_data, 2, mad_values, "/")
      
    } else if (norm == "IRS") {
      eset_matrix_norm[,grep(index,pData(eset)$Batch2)] <- IRS_normalize(eset_matrix[,grep(index,pData(eset)$Batch2)], pData(eset), IRS_column)
    }
    
    # MISSING VALUE IMPUTATION
    # library(mice)
    
    if (impute_missing == TRUE){
      message("-- IMPUTING MISSING VALUES")
      
      eset_matrix_norm[eset_matrix_norm == 0] <- NA
      message("---- NAs BEFORE IMPUTATION:", sum(is.na(eset_matrix_norm[,grep(index,pData(eset)$Batch2)])))
      
      tryCatch({
        # for (i in 1:ncol(eset_matrix_norm)){
        #   eset_matrix_norm[,i] <- complete(mice::mice(eset_matrix_norm, method = "lasso.norm", printFlag = FALSE))[,i]
        # }
        
        if (type == "Proteomics"){
          assign("matrix", eset_matrix_norm[,grep(index,pData(eset)$Batch2)], envir = .GlobalEnv)
        }
        
        eset_imputed <- pcaMethods::nni(t(eset_matrix_norm[,grep(index,pData(eset)$Batch2)]), method = c("llsImpute"), correlation = "pearson", verbose = TRUE)
        eset_matrix_norm[,grep(index,pData(eset)$Batch2)] <- as.matrix(t(pcaMethods::completeObs(eset_imputed)))
        # eset_matrix_norm[,grep(index,pData(eset)$Batch2)] <- as.matrix(complete(mice::mice(eset_matrix_norm[,grep(index,pData(eset)$Batch2)], method = "lasso.norm", printFlag = FALSE)))
        
        message("---- NAs AFTER IMPUTATION:", sum(is.na(eset_matrix_norm[,grep(index,pData(eset)$Batch2)])))
      },
      error = function(cond) {
        message("---- COULD NOT COMPLETE MISSING VALUE IMPUTATION")
        message(paste("*** ERROR MESSAGE:", conditionMessage(cond), "***"))
        # Choose a return value in case of error
        eset_matrix_norm
      })
      
      eset_matrix_norm[is.na(eset_matrix_norm)] <- 0
    }
    
    #---------------------------------------------------------------------------
    
    # RECREATE ESET OBJECT
    eset_norm <- ExpressionSet(assayData = eset_matrix_norm);
    pData(eset_norm) <- pData(eset);
    rownames(pData(eset_norm)) <- colnames(eset_matrix_norm);
    fData(eset_norm) <- fData(eset);
    
    message("COMPLETED NORMALIZATION **")
    
    if (!grepl("none", norm)) { 
      return (eset_norm);
      
    } else {
      return(eset);
      
    }
  }
}

#
#
#

IRS_normalize <- function(eset_matrix,
                          eset_pdata,
                          TMT_groups_column = ''){
  # eset <- eset[apply(eset, 1, FUN = function(x){sum(x == 0)}) < (ncol(eset) * 0.3),]
  # eset_matrix = exprs(eset)
  # eset_pdata = pData(eset)
  # TMT_groups_column = "Run"
  
  TMT_column <- which(colnames(eset_pdata) == TMT_groups_column)
  
  experiment_groups <- eset_pdata[,TMT_column]
  
  data_raw <- eset_matrix
  
  group_columns <- list()
  for (i in unique(experiment_groups)){
    group_columns[[i]] <- which(experiment_groups == i)
  }
  
  # separate the TMT data by experiment
  # we do not need to do this for the normalization factor calculation here,
  # but we will need these data frames for the IRS step below.
  TMT_raw <- list()
  for (i in 1:length(group_columns)){
    TMT_raw[[i]] <- data_raw[,c(group_columns[[i]])]
  }
  
  # first basic normalization is to adjust each TMT experiment to equal signal per channel
  # figure out the global scaling value
  
  col_sums <- c()
  for (i in 1:length(TMT_raw)){
    col_sums <- c(col_sums, colSums(TMT_raw[[i]]))
  }
  
  col_sums <- unname(col_sums)
  
  target <- mean(col_sums)
  
  # do the sample loading normalization before the IRS normalization
  # there is a different correction factor for each column
  TMT_SL <- list()
  
  for (i in 1:length(TMT_raw)){
    norm_facs <- target / colSums(TMT_raw[[i]])
    TMT_SL[[i]] <- sweep(TMT_raw[[i]], 2, norm_facs, FUN = "*")
  }
  
  # make a pre-IRS data frame after sample loading normalizations
  data_SL <- data.frame(TMT_SL[[1]])
  for (i in 2:length(TMT_SL)){
    data_SL <- cbind(data_SL, TMT_SL[[i]])
  }
  
  # boxplot(data_SL)
  
  # perform TMM on the SL-normed data and visualize resulting distributions
  SL_TMM <- edgeR::calcNormFactors(data_SL)
  data_SL_TMM <- sweep(data_SL, 2, SL_TMM, FUN = "/") # data after SL and TMM on original scale
  
  # boxplot(data_SL_TMM)
  
  # make new data frame with row sums from each frame
  IRS <- data.frame(rowSums(TMT_SL[[1]]))
  for (i in 2:length(TMT_SL)){
    IRS[, i] <- rowSums(TMT_SL[[i]])
  }
  
  # get the geometric average intensity for each protein
  IRS$average <- apply(IRS, 1, function(x) exp(mean(log(x))))
  
  # compute the scaling factor vectors
  for (i in 1:length(TMT_SL)){
    IRS <- cbind(IRS, (IRS$average / IRS[,i]))
  }
  
  IRS[IRS == "NaN"] <- NA
  
  ## make new data frame with IRS normalized data
  data_IRS <- TMT_SL[[1]] * IRS[,(length(TMT_SL) + 1 + 1)]
  for (i in 2:length(TMT_SL)){
    data_IRS <- cbind(data_IRS, (TMT_SL[[i]] * IRS[,(length(TMT_SL) + 1 + i)]))
  }
  
  # data_IRS <- na.omit(data.frame(data_IRS))
  # data_IRS  <- data_IRS[apply(data_IRS , 1, FUN = function(x){sum(x == 0)}) < (ncol(eset_matrix)),] 
  # 
  # boxplot((data_IRS))
  
  # this is data after SL, IRS, and TMM normalized on original scale
  IRS_TMM <- edgeR::calcNormFactors(data_IRS)
  data_IRS_TMM <- sweep(data_IRS, 2, IRS_TMM, FUN = "/")
  
  # boxplot(log2(data_IRS_TMM))
  
  return(data_IRS_TMM)
}

