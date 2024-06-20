# NCmisc::list.functions.in.file("limmaLinearModel.R", alphabetic = TRUE)

################################################################################

#' @title Limma Linear Model / Differential Analysis
#'
#' @description 
#' Run limma::lmFit for differential analysis by sample groups.
#' 
#' @param annot processed Annotation object from formatAnnotation()
#' @param eset processed Expression Set object from makeEset()
#' @param samples list of strings - which groups to include as factors in model
#' @param contrast_fit logical input - whether to perform fit on contrast matrix
#' @param time_series logical input - whether time series input is included in annotation object
#' @param pairs logical input - whether pairs input is included in annotation object
#' @param coef_names logical input - if TRUE return list of coeficient string options
#' @param int_term logical input - whether use opts to include interaction term
#' @param int string input - specify interaction term
#' @param int_base string input - specify what group interaction term is interacting with
#'
#' @return limma fit object: to be put into topTable() to generate table of top expressed rows
#' @return list of strings: coefficient options for topTable
#' 
#' @import limma 
#'

limmaLM = function(annot,
                   eset,
                   samples,
                   contrast_fit = FALSE,
                   time_series = FALSE,
                   time_col,
                   time_points_cont = 0,
                   time_points_disc,
                   time = "",
                   coef_names = FALSE,
                   return_contrasts = FALSE,
                   contrasts_subset,
                   covariate = FALSE,
                   covariate_col,
                   remove_batch_PCA = FALSE,
                   batch_column = "") {
  
  ##############################################################################
  
  # STANDARD CONTRASTS
  contrastgroups = samples
  contrast_strings = c()
  
  if(length(contrastgroups > 1)){ # compare all groups to each other, assume first group is control
    contrast_strings <- sapply(combn(rev(contrastgroups), 2, simplify = F), FUN = function(l)paste(l, collapse = "-"))
    contrast_strings <- c(contrast_strings, sapply(combn(contrastgroups, 2, simplify = F), FUN = function(l)paste(l, collapse = "-")))
    
  } else if(length(contrastgroups) == 1){ # if only one group listed
    contrast_strings <- contrastgroups
    
  }
  
  #-----------------------------------------------------------------------------
  
  # CONTRASTS SPECIFIC TO TIME SERIES
  # for continuous option of time variable
  
  if(time_series == TRUE && time == "continuous"){
    new_groups = c()

    for (k in 1:length(contrastgroups)){ # make contrast strings based on the TimeSeries variable
      for(l in 1:time_points_cont){
        new_groups <- c(new_groups, paste(contrastgroups[k], "_Time_", l, sep=""))
        
      }
    }
    new_contrastgroups <- c(new_groups) # combine strings
    contrast_strings <- sapply(combn(rev(new_contrastgroups), 2, simplify = F), FUN = function(l) paste(l, collapse = "-"))
    contrast_strings <- c(contrast_strings, sapply(combn(new_contrastgroups, 2, simplify = F), FUN = function(l)paste(l, collapse = "-")))
  }
  
  # for discrete option of time variable
  if(time_series == TRUE && time == "discrete"){
    new_groups = c()
    
    for (k in 1:length(contrastgroups)){ # make contrast strings based on the TimeSeries variable
      for(l in time_points_disc){
        new_groups <- c(new_groups, paste(contrastgroups[k], "_Time_", l, sep=""))
        
      }
    }
    new_contrastgroups <- c(new_groups) # combine strings
    contrast_strings <- sapply(combn(rev(new_contrastgroups), 2, simplify = F), FUN = function(l) paste(l, collapse = "-"))
    contrast_strings <- c(contrast_strings, sapply(combn(new_contrastgroups, 2, simplify = F), FUN = function(l)paste(l, collapse = "-")))
  }
  
  ##############################################################################
  
  if (return_contrasts == TRUE) {
    return(contrast_strings);
  }
  
  # contrasts_subset = c("Glucose-NoGlucose_Time_15", "Glucose-NoGlucose_Time_30")
  if (length(contrasts_subset) >= 1){
    contrast_strings = contrasts_subset
  }
  
  ##############################################################################
  
  eset <- eset[,which(pData(eset)$Group %in% contrastgroups)]
  factor <- factor(pData(eset)$Group)
  
  # MAKE DESIGN ----------------------------------------------------------------
  
  ## STANDARD/DEFAULT LINEAR MODEL
  design <- stats::model.matrix(~ 0 + factor);
  
  # REMOVE BATCH EFFECT VISUALIZATION ------------------------------------------
  if (remove_batch_PCA == TRUE){
    batch_corrected_matrix <- limma::removeBatchEffect(eset, pData(eset)[, batch_column], design = design)
    
    batch_corrected_eset <- ExpressionSet(assayData = batch_corrected_matrix)
    pData(batch_corrected_eset) <- pData(eset);
    rownames(pData(batch_corrected_eset)) <- colnames(batch_corrected_matrix);
    fData(batch_corrected_eset) <- fData(eset);
    
    return(batch_corrected_eset);
  }
  
  #-----------------------------------------------------------------------------
  
  # LINEAR MODEL + TIME SERIES
  if(time_series == TRUE){
    if(time == "continuous"){
      Time <- splines::ns(as.numeric(pData(eset)[,time_col]), df = time_points_cont)
      colnames(Time) <- paste0("_", colnames(Time))
    
    } else if(time == "discrete"){
      Time <- pData(eset)[,time_col]
      Time <- paste0("_", Time)
      Time <- factor(Time)
    }

    design <- stats::model.matrix(~ 0 + factor + factor:Time);
  }
  
  #-----------------------------------------------------------------------------
  
  # LINEAR MODEL + COVARIATE (NOT TIME)
  if(covariate == TRUE){
    if (length(covariate_col) >= 1){
      covariate_1 <- pData(eset)[,covariate_col[1]]
      covariate_1 <- paste0("_", covariate_1)
      covariate_1 <- factor(covariate_1)
      
      design <- stats::model.matrix(~ 0 + factor + covariate_1, data = eset);
      
      if (time_series == TRUE) {design <- stats::model.matrix(~ 0 + factor + factor:Time + covariate_1);}
    }
    
    if (length(covariate_col) >= 2){
      covariate_2 <- pData(eset)[,covariate_col[2]]
      covariate_2 <- paste0("_", covariate_2)
      covariate_2 <- factor(covariate_2)
      
      design <- stats::model.matrix(~ 0 + factor + covariate_1 + covariate_2, data = eset);
      
      if (time_series == TRUE) {design <- stats::model.matrix(~ 0 + factor + factor:Time + covariate_1 + covariate_2);}
    }
    
    if (length(covariate_col) == 3){
      covariate_3 <- pData(eset)[,covariate_col[3]]
      covariate_3 <- paste0("_", covariate_3)
      covariate_3 <- factor(covariate_3)
      
      design <- stats::model.matrix(~ 0 + factor + covariate_1 + covariate_2 + covariate_3, data = eset);
      
      if (time_series == TRUE) {design <- stats::model.matrix(~ 0 + factor + factor:Time + covariate_1 + covariate_2 + covariate_3);}
    }
  }
  
  # CLEAN DESIGN ---------------------------------------------------------------
  colnames(design) <- gsub("factor", "", colnames(design))
  colnames(design) <- gsub(":", "_", colnames(design))
  colnames(design) <- make.names(colnames(design));
  
  if (time_series == TRUE) {
    if(time == "discrete"){
      colnames = contrastgroups
      
      for(group in contrastgroups){
        colnames <- c(colnames, paste(group, "_Time_", time_points_disc, sep = ""))
      }
      
      design <- design[, c(colnames)]
    }
  }
  
  # FIT ------------------------------------------------------------------------
  
  # https://www.biostars.org/p/9509937/ 
  # In this case OP doesn't have the data required for input to DEqMS. I think that 
  # allowing the variance to depend on the number of PSMs when that data is available 
  # is a good idea, but limma can already do that as part of its native code. 
  # For example you could use:
  # 
  # fit <- lmFit(y, design)
  # fit$Amean <- log(PSM)
  # fit <- eBayes(fit, trend=TRUE)
  # fit$Amean <- rowMeans(y)
  #
  # Where PSM is the PSM count and then you would already have the same variance/PSM 
  # approach that DEqMS proposes without the need for an extra package and an extra 
  # layer of analysis. In the devel version of limma we've made it one step easier again:
  #     
  # fit <- lmFit(y, design)
  # fit <- eBayes(fit, trend=log(PSM))
  # 
  # The use of native limma code further allows the PSM variance trend to be combined 
  # with robust empirical Bayes, something that the DEqMS code does not allow.
  #   
  
  # According to MaxQuant manual https://www.cores.emory.edu/eipc/_includes/documents/definitions-mmaxquant.pdf
  # MS/MS count: Peptide spectrum matches
  
  fit <- limma::lmFit(eset, design)
  
  if(contrast_fit == FALSE){
    if(coef_names == TRUE){
      return(colnames(design));
    }
    
    fit <- limma::eBayes(fit);
    
  } else if(contrast_fit == TRUE){
    contrast_matrix <- limma::makeContrasts(contrasts = contrasts_subset, levels = design);
    
    if(coef_names == TRUE){
      return(colnames(contrast_matrix));
    }
    
    fit <- limma::contrasts.fit(fit, contrast_matrix);
    fit <- limma::eBayes(fit);
  }
  return(fit);
  
}
  
################################################################################  

#' @title Top Table
#'
#' @description
#' Generation of limma top table object using the fit object generated via limmaLM.
#' 
#' @param fit limma fit object output from limmaLM function
#' @param coef string input - which column/contrast to use as the coefficient in the topTable
#' @param coef_options list strings input - list of possible coefficients
#' @param number_genes numeric input - how many genes to include in topTable
#' @param logfc_th_add logical input - whether to add logFC threshold cutoff to topTable
#' @param logfc_th numeric input - if logfc_th_add TRUE, what logFC threshold is
#' 
#' @return limma topTable object, table of genes with logFC and p-val based on model fit
#'
#' @import limma
#' 
#' @example makeTopTable(lmFit(object), "Glucose-NoGlucose", c("Glucose-NoGlucose", "NoGlucose-Glucose"))
#'

makeTopTable <- function(fit,
                         coef,
                         coef_options,
                         logfc_th = 0){
  
  table <- limma::topTable(fit, 
                           adjust = "BH", 
                           number = nrow(fit$coefficients), 
                           coef = which(coef_options == coef), 
                           lfc = logfc_th)
  
  end_col <- ncol(table)
  start_col <- which((colnames(table)) == "feature_identifier")
  top_table <- table[,seq(start_col, end_col)]
  
  return(top_table);
  
}

################################################################################

#' @title P-Value Distribution Plot
#' 
#' @desccription
#' 
#' @param top_table dataframe
#' @param pval_type string P.Value or adj.P.Val
#' @param threshold numeric
#' 

plotPvalHistogram <- function(top_table,
                              pval_type = "P.Value",
                              threshold = 0.05){
  
  p <- top_table[, pval_type]
  
  title <- "P-Value"
  if (pval_type == "adj.P.Val") {title <- "BH Adjusted P-Value"}
  
  top_table$Significant <- ifelse(p <= threshold, "sig", "insig")
  
  ggplot(data = top_table, aes(x = top_table[, pval_type], fill = Significant)) +
    geom_histogram(colour = "black", bins = 35) +
    theme_minimal() +
    labs(title = paste0(title, " Distribution"),
         x = pval_type) +
    {if ("sig" %in% top_table$Significant) {scale_fill_manual(values = c("grey", "seagreen3"))}} +
    {if (!("sig" %in% top_table$Significant)) {scale_fill_manual(values = c("grey"))}} # tomato
  
  # Info on interpreting p-value distributions: http://varianceexplained.org/statistics/interpreting-pvalue-histogram/
}

################################################################################

#' @title Limma Model Equation
#' 
#' @description
#' Outputs a string object that shows what the linear model equation is based on 
#' the chosen inputs.
#' 
#' @param groups list of groups included in the model
#' @param covariates list of covariates included in the mocel
#' @param time_series logical input - whether time_series is included
#' 
#' @output print string object of equation
#' 

limmaEquation = function(groups = "",
                         include_covariates = FALSE,
                         covariate_col = NULL,
                         time_series = FALSE){
  
  group <- paste0(unlist(groups), collapse = "_")
  covariates <- ''
  time <- ''
  
  if (include_covariates == TRUE){
    covariates <- paste0(unlist(covariate_col), collapse = " + ")
    covariates <- paste0(" + ", covariates)
  }
  
  if (time_series == TRUE){
    time <- paste0( " + ", group, ":Time")
  }
  
  model = paste0("~ 0 + ", group, time, covariates)
    
  
  return(model);
  
}

################################################################################
