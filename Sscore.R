
# data_list <- makeDataList(top_table_list, data_format)
# weighted_zi_list <- computeWeightedZIList(data_list)
# sscore_dataframe <- sscoreIntegration(data_list, adj_pval_cutoff = 0.05)

# MAKE DATA LIST ###############################################################
#' @title Make Data List from TopTable Objects
#' 
#' @description
#' Use a list of limma::topTable dataframes for the same contrast but different 
#' datasets (i.e. proteomics, phosphoproteomics, metabolomics, etc.) to format a
#' data_list which contains reformatted versions of these tables named based on 
#' the data format of that dataset.
#'
#' @param top_table_list A list of dataframes from limma::topTable
#' @param data_format A list of strings which contain the data_format for each
#' uploaded dataset (i.e. ProteinGroups, PhosphoSites, etc.).
#'
#' @output A list of dataframes.
#'

makeDataList <- function(top_table_list,
                         data_format) {
  
  data_list <- list()
  
  for (i in 1:length(top_table_list)){
    # For non-metabolomics dataframes
    if(!grepl("Metabol", data_format[i])) {
      dataframe <- data.frame(matrix(nrow = nrow(top_table_list[[i]]), ncol = 3))
      dataframe[,1] <- sub(';.*', '', top_table_list[[i]]$Protein)
      dataframe[,2] <- rownames(top_table_list[[i]])
      dataframe[,3] <- top_table_list[[i]]$logFC
      colnames(dataframe) <- c("uniprot_id", "feature_id", "logfc")
      data_list[[i]] <- dataframe
      
    # For metabolomics dataframes
    } else if (grepl("Metabol", data_format[i])){
      dataframe <- data.frame(matrix(nrow = nrow(top_table_list[[i]]), ncol = 3))
      if ("identifier" %in% colnames(top_table_list[[i]])){
        dataframe[,1] <- sub(',.*', '', top_table_list[[i]]$identifier) # assuming HMDB column called "identifier" as in openMS output
        dataframe[,1] <- sub('.*c.{1}', '', dataframe[,1])
      } else if ("HMDBID" %in% colnames(top_table_list[[i]])){
        dataframe[,1] <- sub(';.*', '', top_table_list[[i]]$HMDBID) # assuming HMDB column called "HMDBID" as in MetaboAnalyst annotated_peaklist output
      } else if ("InChIKey" %in% colnames(top_table_list[[i]])){
        dataframe[,1] <- top_table_list[[i]]$InChIKey # assuming InChIKey column called "InChIKey" as in MS-DIAL output
      }
      
      dataframe[,2] <- rownames(top_table_list[[i]])
      dataframe[,3] <- top_table_list[[i]]$logFC
      
      colnames(dataframe) <- c("ID", "feature_id", "logfc")
      
      dataframe$ID <- gsub("^$|^ $", NA, dataframe$ID)
      dataframe <- na.omit(dataframe)
      metabolite_mapping <- readr::read_rds("PCSF_files/metaboliteMap.rds")
      
      if (grepl("^HMDB|^hmdb", dataframe$ID[1])) {
        dataframe <- merge(dataframe, metabolite_mapping[, c("hmdb_id", "chebi_id")], by.x = "ID", by.y = "hmdb_id")
      } else if (grepl("[A-Z]{14}-[A-Z]{10}-[A-Z]{1}", dataframe$ID[1])){
        dataframe <- merge(dataframe, metabolite_mapping[, c("inchikey_id", "chebi_id")], by.x = "ID", by.y = "inchikey_id")
      } else if (grepl("C[0-9]{5}", dataframe$ID[1])){
        dataframe <- merge(dataframe, metabolite_mapping[, c("kegg_id", "chebi_id")], by.x = "ID", by.y = "kegg_id")
      } else{
        dataframe$other_id <- dataframe$ID
      }
      
      # Other IDs are mostly numeric of different lengths. Would need to tell users to add identifier
      # before number to be able to tell the difference, i.e. "PUBCHEM:01234".
      
      dataframe$feature_id <- dataframe$ID
      dataframe$ID <- NULL
      dataframe <- na.omit(dataframe)
      rownames(dataframe) <- 1:nrow(dataframe)
      dataframe <- dataframe[, c("chebi_id", "feature_id", "logfc")]
      
      # assign("metabolite_sscoreInput", dataframe, envir = .GlobalEnv)
      
      data_list[[i]] <- dataframe
    }
  }
  
  names(data_list) <- data_format
  assign("data_list", data_list, envir = .GlobalEnv)
  return(data_list)
}

# UNIPROT ID TO GENE SYMBOL ####################################################
#' @title Convert Uniprot ID to Gene Symbol
#'
#' @description 
#' Using Uniprot protein IDs from dataframe, create a column of gene symbols.
#'
#' @param input_df S-score dataframe with `uniprot_id` column
#' @param key_column A string of the column name containing Uniprot IDs
#' @param keytype What type of ID is being mapped (Uniprot).
#'
#' @return A dataframe with a `gene_symbol` column.
#'

# function to get gene_symbol from uniprot_ids
convertUniprotSymbol <- function(input_df, 
                                 key_column = "uniprot_id", 
                                 keytype = "UNIPROT") {
  
  # Extract keys from the specified column
  keys <- input_df[[key_column]]
  
  input_df$gene_symbol <- rep("NA", nrow(input_df))
  
  # Perform mapping using AnnotationDbi::mapIds
  # TODO: EDIT SO MORE THAN HUMAN ACCEPTED 
  try({input_df$gene_symbol <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys, keytype, column = "SYMBOL")}) # human
  try({input_df$gene_symbol <- AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db, keys, keytype, column = "SYMBOL")}) # mouse
  
  # Rename the "SYMBOL" column to "gene_symbol" in the original dataframe
  try({input_df$gene_symbol <- stringr::str_squish(input_df$gene_symbol)})
  
  # Return the modified dataframe
  return(input_df)
}

# COMPUTE WEIGHTED Z-SCORE LIST ################################################
#' @title Compute Weighted ZI List
#' 
#' @description 
#' Compute weighted zi for each modified top_table and save the output as a list.
#' 
#' @param data_list A list of dataframes, each one containing a general ID column (Uniprot
#' or Chebi), a feature ID column (from topTable) and a log Fold Change column (from
#' topTable).
#'
#' @return A list of modified dataframes that contains additional columns for wk,
#' weighted zi, & zi
#'

# computation of weighted zi and save as a list
computeWeightedZIList <- function(data_list = data_list) {
  
  message(paste("There are", length(data_list), "elements in data_list") )
  message(" ")
  message("*** Input data format ***")
  message(" -- The input 'data_list' is a list of 'dataframes' with 3 columns: uniprot_id, feature_id and logfc")
  message(" -- The code will combine all elements of 'data_list' using 'uniprot_id'.")
  message(" -- If present, S-score for metabolite dataset will be computed separately and then combined with genes")
  message(" -- Use 'my_combine_genes_and_metabols_using_sscore' function if metabolomics data is present.")
  message(" -- Use 'my_combine_genes_using_sscore' function if metabolomics data is NOT present.")
  message(" -- Example 'feature_id' for proteomics: Q8BH50_ARK2N")
  message(" -- Example 'feature_id' for phosphoproteomics & other PTMs: Q8BH50_ARK2N_T74")
  message(" ")
  message("*** Citation ***")
  message("Citation for original publication: Nat Commun. 2013:4:2617. doi: 10.1038/ncomms3617" )
  
  # Define the biological importance weights
  bio_importance_weights <- list(ProteinGroups = 0.5, 
                                 PhosphoSites = 0.5, 
                                 RNA = 0.3,
                                 MetaboliteNeg = 0.5,
                                 MetabolitePos = 0.5,
                                 Peptides = 0.5,
                                 Generic = 0.5)
  
  assign("data_list", data_list, envir = .GlobalEnv)
  
  # Calculate weights based on the number of features and biological importance
  if ("ProteinGroups" %in% names(data_list)){
    ProteinGroups_weight <- bio_importance_weights$ProteinGroups / sqrt(nrow(data_list$ProteinGroups))
    
  }
  if ("PhosphoSites" %in% names(data_list)){
    PhosphoSites_weight <- bio_importance_weights$PhosphoSites / sqrt(nrow(data_list$PhosphoSites))
    
  }
  if ("RNA" %in% names(data_list)){
    RNA_weight <- bio_importance_weights$RNA / sqrt(nrow(data_list$RNA))
    
  }
  if ("MetaboliteNeg" %in% names(data_list)){
    MetaboliteNeg_weight <- bio_importance_weights$MetaboliteNeg / sqrt(nrow(data_list$MetaboliteNeg))
    
  }
  if ("MetabolitePos" %in% names(data_list)){
    MetabolitePos_weight <- bio_importance_weights$MetabolitePos / sqrt(nrow(data_list$MetabolitePos))
    
  }
  if ("Peptides" %in% names(data_list)){
    Peptides_weight <- bio_importance_weights$Peptides / sqrt(nrow(data_list$Peptides))
    
  }
  if ("Generic" %in% names(data_list)){
    Generic_weight <- bio_importance_weights$Generic / sqrt(nrow(data_list$Generic))
    
  }
  
  # COMPUTE WEIGHTED Z-SCORES ####################################################
  #' @title Compute Weighted Z-Scores
  #'
  #' @description
  #' Compute the weighted Z-score for each row using the logFC.
  #' 
  #' @param x dataframe from `data_list` described in `computeWeightedZIList` function
  #' 
  #' @return dataframe with additional columns
  #'
  
  # function for computing weighed z-scores
  computeWeightedZI <- function(x, dataset_type) {
    wk <- switch(dataset_type,
                 "ProteinGroups" = ProteinGroups_weight, 
                 "PhosphoSites" = PhosphoSites_weight, 
                 "RNA" = RNA_weight,
                 "MetaboliteNeg" = MetaboliteNeg_weight,
                 "MetabolitePos" = MetabolitePos_weight,
                 "Peptides" = Peptides_weight,
                 "Generic" = Generic_weight)
    
    x %>%
      dplyr::mutate(wk = wk,
                    zi = scale(logfc, center = TRUE, scale = TRUE)[,1],
                    weighted_zi = wk * zi) %>%
      dplyr::select(logfc, wk, weighted_zi, everything()) %>%
      dplyr::mutate(across(where(is.character), ~ na_if(., "NA"))) %>%
      as.data.frame()
  }
  
  for (i in 1:length(data_list)) {
    
    # Extract data frame from the list
    data_s <- data_list[[i]]
    
    # Extract element name from the list
    element_name <- names(data_list)[i]
    
    # Compute weighted Z-scores using the new approach
    data_s1 <- computeWeightedZI(data_s, element_name) %>%
      dplyr::rename_with(~ ifelse(.x %in% c("uniprot_id", "chebi_id"), .x, paste0(.x, "_", element_name)))
    
    # Save the computed data frame back into the list
    data_list[[i]] <- data_s1
    
  }
  
  # Return the modified list (optional)
  return(data_list)
}

# COMPUTE S-SCORES #############################################################
#' @title Compute S-Scores
#'
#' @description
#' Using previous calculations of weighted z-score, etc. calculate s-score for each
#' row of the combined dataframe and output this updated dataframe.
#' 
#' @param x Modified data list generated using `computeWeightedZIList` function.
#' @param ID_column The name of the column containing IDs consistent across data frames,
#' generally "uniprot_id" or "chebi_id"
#'
#' @return Modified dataframe with s-score calculations.
#'

# function for computing S-scores
computeSscore <- function(x, ID_column) {
  if (ID_column == "uniprot_id"){
    print("UniProt ID Column")
    
    x_sscore <- x %>%
      dplyr::mutate(
        uniprot_id = stringr::str_squish(uniprot_id),
        comb_weighted_zi = rowSums(dplyr::select(., dplyr::starts_with("weighted_zi_")), na.rm = TRUE),
        comb_wk = sqrt(rowSums(dplyr::select(., dplyr::starts_with("wk_"))^2, na.rm = TRUE))
      ) %>%
      dplyr::mutate(
        sscore = comb_weighted_zi / comb_wk,
        sscore_pval = stats::pnorm(abs(sscore), lower.tail = FALSE) * 2,
        sscore_adj_pval = stats::p.adjust(sscore_pval, method = 'BH')
      ) %>%
      dplyr::mutate(across(where(is.character), ~ na_if(., "NA"))) %>%
      as.data.frame() %>%
      dplyr::select(
        uniprot_id,
        dplyr::starts_with("feature_id"),
        dplyr::starts_with("gene_symbol"),
        sscore,
        sscore_pval,
        sscore_adj_pval,
        dplyr::starts_with("comb_"),
        dplyr::starts_with("logfc"),
        dplyr::starts_with("weighted"),
        dplyr::starts_with("wk")
      )
    
  } else if (ID_column == "chebi_id") {
    x_sscore <- x %>%
      dplyr::mutate(
        chebi_id = stringr::str_squish(chebi_id),
        comb_weighted_zi = rowSums(dplyr::select(., dplyr::starts_with("weighted_zi_")), na.rm = TRUE),
        comb_wk = sqrt(rowSums(dplyr::select(., dplyr::starts_with("wk_"))^2, na.rm = TRUE)),
        sscore = comb_weighted_zi / comb_wk,
        sscore_pval = stats::pnorm(abs(sscore), lower.tail = FALSE) * 2,
        sscore_adj_pval = stats::p.adjust(sscore_pval, method = 'BH'),
        feature_id_metabolome = chebi_id,
        uniprot_id = "NA",
        gene_symbol = "NA",
      ) %>%
      dplyr::mutate(across(where(is.character), ~ na_if(., "NA"))) %>%
      as.data.frame() %>%
      dplyr::select(
        uniprot_id,
        chebi_id,
        dplyr::starts_with("feature_id"),
        dplyr::starts_with("gene_symbol"),
        sscore,
        sscore_pval,
        sscore_adj_pval,
        dplyr::starts_with("comb_"),
        dplyr::starts_with("logfc"),
        dplyr::starts_with("weighted"),
        dplyr::starts_with("wk")
      )
  }
  
  message("S-SCORE CALCULATION COMPLETE")
  return(x_sscore)
}

# S-SCORE INTEGRATION ##########################################################
#' @title S-Score Integration
#' ** This is the MAIN FUNCTION which calls ALL ABOVE FUNCTIONS.
#'
#' @description
#' Using above functions, take the data list and generate the s-score integrated
#' dataframe for genes & metabolites.
#' 
#' @param data_list  A list of formatted dataframes from limma::topTable
#' @param adj_pval_cutoff A numerical cutoff for pval, usually 0.05
#'
#' @return A s-score integrated dataframe for all datasets
#'

sscoreIntegration <- function(data_list){
  
  message("--- BEGINNING S-SCORE INTEGRATION")
  ## LOADS ----------------------------------------------------------------------
  prefix__id_proteingroups <- "pr."
  prefix__id_rna <- "rna."
  prefix__id_rna_deseq <- "rna_deseq."
  prefix__id_rna_quant <- "rna_quant."
  prefix__id_phosphosites <- "ph."
  prefix__id_acetylome <- "ac."
  prefix__id_ubiquitylome <- "ub."
  prefix__id_metabolome <- "met."
  prefix__id_metabolitepos <- "metpos."
  prefix__id_metaboliteneg <- "metneg."
  
  message("--- computing weighted z-score list")
  ## COMPUTE WEIGHTED ZI LIST --------------------------------------------------
  data_list_mod <- computeWeightedZIList(data_list)
  
  message(colnames(data_list_mod))
  
  message("--- separating my metabolite vs. non-metabolite")
  # separate by non-metabolomics versus metabolomics
  elements_to_use <- list()
  
  # non-metabolomics dataset names
  elements_to_use[[1]] <- names(data_list_mod)[names(data_list_mod) != "MetaboliteNeg" & 
                                                 names(data_list_mod) != "MetabolitePos"]
  # metabolomics dataset names
  elements_to_use[[2]] <- names(data_list_mod)[names(data_list_mod) == "MetaboliteNeg" | 
                                                 names(data_list_mod) == "MetabolitePos"]
  message(str(elements_to_use))
  
  sscore_combined_list <- list()
  
  message("--- computing s-score")
  ## COMPUTE S-SCORE -----------------------------------------------------------
  # separately for metabolomics combined dataframe and non-metabolomics combined 
  # dataframes
  for (i in 1:length(elements_to_use)){
    ifelse(i == 1, ID_column <- "uniprot_id", ID_column <- "chebi_id")
    
    if(length(elements_to_use[[i]]) != 0){
      sscore_combined_list[[i]] <- data_list_mod[elements_to_use[[i]]] %>% 
        purrr::reduce(full_join, by = ID_column) %>% 
        computeSscore(ID_column = ID_column) %>%
        dplyr::arrange(sscore_adj_pval)
    }
  }
  
  message("SORTED S_SCORE DATAFRAME")
  message(str(sscore_combined_list))
  
  message("--- binding dataframes")
  ## BIND S-SCORE DF ROWS ------------------------------------------------------
  
  if (!is.null(sscore_combined_list[[1]])){
    sscore_combined <- sscore_combined_list[[1]] %>%
      dplyr::select(uniprot_id,
                    dplyr::starts_with("feature_id"), 
                    sscore, 
                    sscore_pval, 
                    sscore_adj_pval, 
                    dplyr::starts_with("comb_"), 
                    dplyr::starts_with("logfc"), 
                    dplyr::starts_with("weighted"), 
                    dplyr::starts_with("wk"))
    
    # if metabolomics included this step combines the non-metabolomics and metabolomics
    # s-score dataframes
    if (length(sscore_combined_list) == 2){
      sscore_combined <- dplyr::bind_rows(sscore_combined, sscore_combined_list[[2]]) %>%
        dplyr::select(uniprot_id,
                      chebi_id, 
                      dplyr::starts_with("feature_id"), 
                      sscore, 
                      sscore_pval, 
                      sscore_adj_pval, 
                      dplyr::starts_with("comb_"), 
                      dplyr::starts_with("logfc"), 
                      dplyr::starts_with("weighted"), 
                      dplyr::starts_with("wk"))
    }
  } else {
    sscore_combined <- sscore_combined_list[[2]] %>%
      dplyr::select(uniprot_id,
                    chebi_id,
                    dplyr::starts_with("feature_id"), 
                    sscore, 
                    sscore_pval, 
                    sscore_adj_pval, 
                    dplyr::starts_with("comb_"), 
                    dplyr::starts_with("logfc"), 
                    dplyr::starts_with("weighted"), 
                    dplyr::starts_with("wk")) %>%
      dplyr::mutate(gene_symbol = "NA")
  }
  
  message("COMBINED S_SCORE DATAFRAME")
  message(str(sscore_combined))
  
  message("--- cleaning & formatting")
  ## CLEAN & FORMAT DATAFRAME --------------------------------------------------
  columns_to_process <- sscore_combined %>% dplyr::select(starts_with("feature_id_")) %>% names()
  
  sscore_combined <- sscore_combined %>%
    convertUniprotSymbol() %>%
    as.data.frame() %>%
    dplyr::mutate(across(where(is.character), ~ na_if(., "NA"))) %>%
    dplyr::mutate(across(where(is.character), ~ na_if(., "NULL"))) %>%
    as.data.frame() %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(columns_to_process), ~stringr::str_replace(., ".*_", ""))) %>%
    dplyr::mutate(dplyr::across(dplyr::contains("feature_id_"), 
                                ~{
                                  column_name <- as.character(cur_column())
                                  matching_string <- str_extract(column_name, "_(.*)")
                                  
                                  if (!is.na(matching_string)) {
                                    prefix_variable <- paste0("prefix_", tolower(gsub("\\d", "", matching_string)))
                                    paste(get(prefix_variable), ., sep = "")
                                  } else {
                                    .  # For columns without matching conditions, keep the original value
                                  }
                                }
    )) %>%
    dplyr::mutate(across(where(is.character), ~ na_if(., "NA"))) %>%
    dplyr::mutate(across(where(is.character), ~ na_if(., "NULL"))) %>%
    as.data.frame() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(sscore_label = paste(c(gene_symbol, dplyr::across(dplyr::all_of(columns_to_process), as.character)), collapse = "_"),
                  feature_id = ifelse(length(sscore_combined_list) == 2, coalesce(uniprot_id, chebi_id), uniprot_id)) %>%
    dplyr::ungroup() %>%
    as.data.frame() %>%
    dplyr::select(sscore_label, feature_id, gene_symbol, sscore, sscore_pval, sscore_adj_pval, everything())
  
  # assign("sscore_dataframe", sscore_combined, envir = .GlobalEnv)
  message("--- finished s-score integration")
  return(sscore_combined);
}

# PCSF FORMATTING ##############################################################
#' @title PCSF Formatting
#' 
#' @description 
#' Take s-score integrated dataframe and format it to be input for PCSF.
#'
#' @param sscore_dataframe Dataframe output from `sscoreIntegration`
#'
#' @return PCSF input formatted dataframe
#'

PSCFFormatting <- function(sscore_dataframe){
  
  sscore_combined_genes_metabols <- sscore_dataframe
  
  # OUTPUT FORMATTED FOR PCSF
  out_for_pcsf_x <- sscore_combined_genes_metabols %>%
    dplyr::filter(sscore_adj_pval <= adj_pval_cutoff) %>%
    dplyr::mutate(uniprot_id = stringr::str_squish(uniprot_id)) %>%
    dplyr::mutate(logfc = sscore * -log10(sscore_adj_pval)) %>% # this is done to remove sscore ties
    dplyr::arrange(uniprot_id, desc(abs(logfc))) %>%
    dplyr::group_by(uniprot_id) %>%
    dplyr::slice_max(order_by = abs(logfc)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(abs(logfc))) %>%
    dplyr::distinct(gene_symbol, .keep_all = TRUE)  %>%
    dplyr::select(gene_symbol, logfc, sscore, sscore_adj_pval)
  
  out_for_pcsf_y <- sscore_combined_genes_metabols %>%
    dplyr::filter(sscore_adj_pval <= adj_pval_cutoff) %>%
    dplyr::mutate(chebi_id = stringr::str_squish(chebi_id)) %>%
    dplyr::mutate(logfc = sscore * -log10(sscore_adj_pval)) %>% # this is done to remove sscore ties
    dplyr::arrange(chebi_id, desc(abs(logfc))) %>%
    dplyr::group_by(chebi_id) %>%
    dplyr::slice_max(order_by = abs(logfc)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-gene_symbol) %>%
    dplyr::rename(gene_symbol = chebi_id) %>%
    dplyr::arrange(desc(abs(logfc))) %>%
    dplyr::distinct(gene_symbol, .keep_all = TRUE)  %>%
    dplyr::select(gene_symbol, logfc, sscore, sscore_adj_pval)
  
  out_for_pcsf <- dplyr::bind_rows(out_for_pcsf_x,
                                   out_for_pcsf_y ) %>%
    readr::write_delim(., paste0(folder_path, "sscore_pcsf_prize_", output_name, ".txt"), delim = "\t")
}

# GSVA FORMATTING ##############################################################
#' @title GSVA Formatting
#' 
#' @description 
#' Take s-score integrated dataframe and format it to be input for GSVA.
#'
#' @param sscore_dataframe Dataframe output from `sscoreIntegration`
#'
#' @return GSVA input formatted dataframe
#'

GSVAFormatting <- function(sscore_dataframe){
  
  sscore_combined_genes_metabols <- sscore_dataframe
  
  # OUTPUT FORMATTED FOR GSVA 
  out_for_gsva_x <- sscore_combined_genes_metabols %>%
    dplyr::mutate(uniprot_id = stringr::str_squish(uniprot_id)) %>%
    dplyr::arrange(uniprot_id, desc(abs(sscore))) %>%
    dplyr::group_by(uniprot_id) %>%
    dplyr::slice_max(order_by = abs(sscore)) %>%
    dplyr::ungroup() %>%
    dplyr::select(gene_symbol, sscore)
  
  out_for_gsva_y <- sscore_combined_genes_metabols %>%
    dplyr::mutate(chebi_id = stringr::str_squish(chebi_id)) %>%
    dplyr::arrange(chebi_id, desc(abs(sscore))) %>%
    dplyr::group_by(chebi_id) %>%
    dplyr::slice_max(order_by = abs(sscore)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-gene_symbol) %>%
    dplyr::rename(gene_symbol = chebi_id) %>%
    dplyr::select(gene_symbol, sscore)
  
  out_for_gsva <- dplyr::bind_rows(out_for_gsva_x,
                                   out_for_gsva_y ) %>%
    dplyr::arrange(desc(abs(sscore)))  %>%
    dplyr::distinct(gene_symbol, .keep_all = TRUE) %>%
    dplyr::rename_with(~ output_name, sscore) %>%
    readr::write_delim(., paste0(folder_path, "sscore_gsva_", output_name, ".txt"), delim = "\t")
  
}

# S-SCORE VOLCANO PLOT #########################################################
#' @title S-Score Volcano Plot
#' 
#' @description 
#' Take s-score integrated dataframe and create a volcano plot with integrated S-score
#' on the x-axis and -log10(adj. p-val) on the y-axis.
#'
#' @param sscore_dataframe Dataframe output from `sscoreIntegration`
#' @param adj_pval_cutoff a numeric input of the cutoff value for significance based
#' on the adjusted p-value
#'
#' @return ggplot volcano plot
#'

sscoreVolcanoPlot <- function(sscore_dataframe, 
                              adj_pval_cutoff = 0.05,
                              label_size = 1,
                              up_color = "red",
                              down_color = "blue",
                              add_labels = TRUE,
                              label_num = 10){
  
  message("--- STARTING S-SCORE VOLCANO PLOT GENERATION")
  sscore_combined_genes_metabols <- sscore_dataframe
  my_pval <- adj_pval_cutoff
  label_pval <- adj_pval_cutoff
  
  # VOLCANO PLOTS
  
  message("--- formatting data")
  data = sscore_combined_genes_metabols %>%
    
    dplyr::select(sscore_label, sscore, sscore_adj_pval) %>%
    
    dplyr::mutate(log10_score_adj_pval = -log10(sscore_adj_pval)) %>%
    
    dplyr::mutate(regulation = case_when((sscore_adj_pval <= my_pval & sscore > 0) ~ "Up",
                                         (sscore_adj_pval <= my_pval & sscore < 0) ~ "Down" )) %>%
    dplyr::mutate(regulation = replace_na(regulation, "NoChange")) %>%
    
    dplyr::mutate(color = case_when(regulation == "Up" ~ up_color,
                                    regulation == "Down" ~ down_color)) %>%
    dplyr::mutate(color = replace_na(color, "gray50")) %>%
    
    dplyr::mutate(size = case_when(regulation == "Up" ~ 4,
                                   regulation == "Down" ~ 4 )) %>%
    dplyr::mutate(size = replace_na(size, 0.2)) %>%
    
    dplyr::mutate(alpha = case_when(regulation == "Up" ~ 1,
                                    regulation == "Down" ~ 1 )) %>%
    dplyr::mutate(alpha = replace_na(alpha, 0.5))
  
  my_top_bott <- function(x, n, wt) {   
    x1 = x %>% dplyr::top_n(., n = n, wt = {{ wt }})   
    x2 = x %>% dplyr::top_n(., n = n, wt = -{{ wt }})   
    x = dplyr::bind_rows(x1, x2)   
    return(x) 
  }
  
  if (add_labels == TRUE){
    message("--- sorting data")
    topdown <- data %>% 
      dplyr::filter(sscore_adj_pval <= label_pval) %>% 
      my_top_bott(n = label_num, wt = sscore) %>%
      dplyr::select(sscore_label,
                    sscore,
                    sscore_adj_pval ) %>%
      dplyr::distinct(.)
  }
  
  # Note: A single dot in the volcano plot may be labelled with multiple features because its the combination of those features that led to the derived s-score.
  message("--- generating plot")
  plot <- ggplot(data, aes (x = sscore,
                            y = -log10(sscore_adj_pval),
                            text = paste("Identifier: ", sscore_label))) +
    
    geom_hline(yintercept = -log10(label_pval), linetype = "dashed") + # -log10(0.05) = 1.3
    
    geom_point(shape = 16, size = data$size, color = data$color, alpha = 0.5) +
    
    scale_color_manual(values = c("red" = "red", 
                                  "blue" = "blue",
                                  "black" = "black")) +
    
    {if(add_labels == TRUE){
      ggrepel::geom_text_repel(data = topdown, 
                               aes(label = sscore_label, size = label_size),
                               show.legend = FALSE, 
                               box.padding = 0.5, 
                               max.overlaps = Inf,
                               min.segment.length = 0)
    }} + 
    
    labs(title = "S-Score Volcano Plot",
         x = "Integrated S-score", 
         y = "-log10 adj. p-value",
         caption = paste("Num UP: ", nrow(data[data$regulation == "Up",]), 
                         " & Num DOWN: ", nrow(data[data$regulation == "Down",]),
                         sep = "")) +
    
    theme_bw() + 
    
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 13),
          axis.title = element_text(size = 16),
          plot.title = element_text(size = 18),
          plot.caption = element_text(size = 12))
  
  return(plot)
  
}

# S-SCORE VENN DIAGRAM #########################################################
#' @title S-Score Dataframe Venn Diagram
#'
#' @description
#' Generate a venn diagram showing the number of significantly expressed genes (adj_pval)
#' compared to the total number of genes.
#' 
#' @param sscore_dataframe Output of `sscoreIntegration` function.
#'
#' @return ggVennDiagram object

sscoreVennDiagram <- function(sscore_dataframe){
  message("** STARTED VENN DIAGRAM **")
  
  sscore_datasets <- sub("feature_id_", "", colnames(sscore_dataframe)[grepl("feature_id_", colnames(sscore_dataframe))])
  sscore_datasets <- sscore_datasets[sscore_datasets != "metabolome"]
  sscore_dataframe[sscore_dataframe == "NA"] <- NA
  
  datasets_dataframes <- list()
  for(i in 1:length(sscore_datasets)){
    col <- which(colnames(sscore_dataframe) == paste("logfc_", sscore_datasets[i], sep = ""))
    dataframe <- sscore_dataframe$sscore_label[!is.na(sscore_dataframe[, col]) & sscore_dataframe$sscore_adj_pval <= 0.05]
    dataframe <- c(dataframe, sscore_dataframe$chebi_id[!is.na(sscore_dataframe[, col]) & sscore_dataframe$sscore_adj_pval <= 0.05])
    dataframe <- dataframe[!is.na(dataframe)]
    dataframe <- dataframe[dataframe != ""]
    dataframe <- unique(dataframe)
    
    datasets_dataframes[[i]] <- as.list(dataframe)
  }
  
  names(datasets_dataframes) <- sscore_datasets
  
  # str(venn_datalists)
  # venn <- ggVennDiagram(datasets_dataframes)
  venn1 <- ggvenn::ggvenn(datasets_dataframes, 
                         stroke_size = 0.5, 
                         set_name_size = 6, 
                         text_size = 6,
                         fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF", "#44A849", "#9561CC", "#868686FF"),
                         stroke_color = "white")
  
  datasets_dataframes <- list()
  for(i in 1:length(sscore_datasets)){
    col <- which(colnames(sscore_dataframe) == paste("logfc_", sscore_datasets[i], sep = ""))
    dataframe <- sscore_dataframe$sscore_label[!is.na(sscore_dataframe[, col])]
    dataframe <- c(dataframe, sscore_dataframe$chebi_id[!is.na(sscore_dataframe[, col])])
    dataframe <- dataframe[!is.na(dataframe)]
    dataframe <- dataframe[dataframe != ""]
    dataframe <- unique(dataframe)
    
    datasets_dataframes[[i]] <- as.list(dataframe)
  }
  
  names(datasets_dataframes) <- sscore_datasets
  
  # str(venn_datalists)
  # venn <- ggVennDiagram(datasets_dataframes)
  venn2 <- ggvenn::ggvenn(datasets_dataframes, 
                          stroke_size = 0.5, 
                          set_name_size = 6, 
                          text_size = 6,
                          fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF", "#44A849", "#9561CC", "#868686FF"),
                          stroke_color = "white")
  
  # assign("venn_datalists", datasets_dataframes, envir = .GlobalEnv)
  
  message("COMPLETED VENN DIAGRAM")
  
  venn <- gridExtra::grid.arrange(venn1 + labs(title = "Significant Variables"), # + theme(plot.title = element_text(size = 16))
                                  venn2 + labs(title = "All Variables"),
                                  ncol = 2)
  return(venn);
}

# S-SCORE ENRICHMENT ###########################################################

#' @title S-Score Enrichment GeneList
#' 
#' @description
#' Format the s-score output into a two column dataframe for use in enrichment, and
#' conversion of that dataframe into a named list of genes.
#' 
#' @param sscore_output dataframe of S-score result output.
#' 
#' @return Named gene list with s-score rank.
#' 

formatSscoreGeneList <- function(sscore_output){
  # assign("sscore_enrichment_input", sscore_output, envir = .GlobalEnv)
  sscore_output[is.na(sscore_output)] <- ""
  
  enrich_input <- data.frame(matrix(ncol = 0, nrow = nrow(sscore_output)))
  
  if ("chebi_id" %in% colnames(sscore_output)){
    enrich_input$feature <- paste0(sscore_output$gene_symbol, sscore_output$chebi_id)
    enrich_input$value <- sscore_output$sscore
  } else {
    enrich_input$feature <- sscore_output$gene_symbol
    enrich_input$value <- sscore_output$sscore
  }
  
  enrich_input <- enrich_input[enrich_input[,1] != "", ]
  enrich_input <- enrich_input[order(abs(as.numeric(enrich_input[,2])), decreasing = TRUE),]
  enrich_input <- enrich_input[!duplicated(enrich_input[,1]),]
  
  geneList <- as.numeric(enrich_input[,2])
  names(geneList) <- as.character(enrich_input[,1])
  geneList = sort(geneList, decreasing = TRUE)
  return(geneList)
}

#' @title S-Score Enrichment
#' 
#' @description
#' Using the s-score gene list, enrich with the specified GMT dataset.
#' 
#' @param geneList named list input - s-score gene list
#' 
#' @return clusterProfiler enriched object
#' 

runSscoreEnrichment <- function(geneList,
                                gmt = "",
                                pval_cutoff = 0.05){
  # GENERATE GMT REFERENCES
  GMT_file <- list.files("./GMTs", pattern = gmt)
  my_geneset = readr::read_delim(paste0("./GMTs/", GMT_file))
  colnames(my_geneset) <- c("pathway", "feature_ids")
  
  my_geneset = dplyr::mutate(my_geneset, term = pathway)
  my_geneset$feature_ids <- gsub("\\-.*","",my_geneset$feature_ids)
  term2features = my_geneset %>% dplyr::select(term, feature_ids)
  term2pathway = my_geneset %>% dplyr::select(term, pathway)
  rm(my_geneset)
  
  enriched <- "NONE"
  
  try({
    # RANKED ENRICHMENT
    enriched <- clusterProfiler::GSEA(gene = geneList,
                                      pAdjustMethod = "fdr",
                                      pvalueCutoff = pval_cutoff,
                                      minGSSize = 2,
                                      eps = 0,
                                      TERM2GENE = term2features,
                                      TERM2NAME = term2pathway)
      
  })
  
  return(enriched)
}
