
# NCmisc::list.functions.in.file("server.R", alphabetic = TRUE)

source("loads.R")
source("makeEset.R")
source("formatAnnotation.R")
source("drawRange.R")
source("getUniprotAnnotation.R")
source("drawHeatmaps.R")
source("drawPCA.R")
source("drawUMAP.R")
source("limmaLinearModel.R")
source("drawVolcano.R")
source("performNormalization.R")
source("saveXLSX.R")
source("enrichment.R")
source("drawCorrelation.R")
source("drawNormalization.R")
source("Sscore.R")
source("PCSF.R")

# Track sessions outside the server function so the variable is permanent
# Alternative is to write to database or file
users = reactiveValues(
  logTable = data.frame(
    id = character(),
    user = character(),
    login = character(),
    logout = character()
  )
)

# SERVER #######################################################################

server <- function(input, output, session) {
  
  #Register session start. Isolate needed for reactive environment
  isolate({
    users$logTable = rbind(
      users$logTable,
      list(id = session$token, 
           user = Sys.getenv("RSTUDIO_USER_IDENTITY"),
           login = as.character(Sys.time()),
           logout = NA)
    )
  })
  
  # INCREASE MAX FILE SIZE TO 30KB
  options(shiny.maxRequestSize = 60 * 1024 ^ 2)
  
  ## SESSION INFO ##############################################################
  output$session_info <- renderPrint({
    sessionInfo()
  })
  
  ## DATA ######################################################################
  
  ### REACTIVE OBJECTS ---------------------------------------------------------
  data <- reactive({
    files <- list()
    if(input$use_example_data == FALSE) {
      for (i in 1:length(input$data_file[, 1])){
        files[[i]] <- utils::read.delim(input$data_file[[i, 'datapath']], header = TRUE, sep = "\t")
      }
      
    } else if (input$use_example_data == TRUE){
      files[[1]] <- utils::read.delim("example_data/example_proteingroups.txt", header = TRUE, sep = "\t")
      files[[2]] <- utils::read.delim("example_data/example_phosphosites.txt", header = TRUE, sep = "\t")
      files[[3]] <- utils::read.delim("example_data/example_openms_annotated_output_pos.txt", header = TRUE, sep = "\t")
      files[[4]] <- utils::read.delim("example_data/example_openms_annotated_output_neg.txt", header = TRUE, sep = "\t")
      
    }
    
    return(files)
  })
  
  ### RENDER UI ----------------------------------------------------------------
  output$dataOverview <- renderUI({
    req(isTruthy(input$data_file) || isTruthy(input$annotation_file) || isTruthy(input$use_example_data))
    h3("Data Overview")
  })
  
  output$dataOverview_help <- renderUI({
    req(isTruthy(input$data_file) || isTruthy(input$annotation_file) || isTruthy(input$use_example_data))
    helpText("Interpreting the data summary: ", "\n", "The top line identifies 
               the structure of the entry, the number of observations (or rows) and the number 
               of variables (or columns). Underneath there are column names on the left followed by a colon, 
               then the variable type (i.e. character, number, factor) followed by 
               the first four or so entries in the column.")
  })
  
  output$dataFile <- renderUI({
    req(isTruthy(input$data_file) || isTruthy(input$use_example_data))
    h5("Data File(s)")
  })
  
  output$filesHeader <- renderUI({
    req(isTruthy(input$data_file) || isTruthy(input$annotation_file))
    h4("File(s)")
  })
  
  ### PRINT/TABLES/PLOTS -------------------------------------------------------
  output$summary_data <- renderPrint({
    req(isTruthy(input$data_file) || isTruthy(input$use_example_data))
    cat("Loaded ", length(data()), " data file(s) for analysis.")
  })
  
  output$file_table <- renderTable({
    req(input$data_file)
    
    input$data_file
  })
  
  ## ANNOTATION ################################################################
  
  ### REACTIVE OBJECTS ---------------------------------------------------------
  annotation_files <- reactive({
    if (input$use_example_data == FALSE){
      annot <- data.frame(openxlsx::read.xlsx(input$annotation_file$datapath, 1, colNames = FALSE))
      
    } else if (input$use_example_data == TRUE){
      annot <- data.frame(openxlsx::read.xlsx("example_data/example_annotation.xlsx", 1, colNames = FALSE))
      
    }
    
    annot <- annot[-(1:10),]
    annot <- annot[,-c(4)]
    colnames(annot) <- c("File", "Name", "Format")
    annot <- na.omit(annot)
    
    return(annot)
  })
  
  annotation_initial <- reactive({
    if (input$use_example_data == FALSE){
      req(input$annotation_file)
      ext <- tools::file_ext(input$annotation_file$datapath)
      shiny::validate(need(ext == "xlsx", "Please upload a .xlsx for the annotation file."))
      annot <- data.frame(openxlsx::read.xlsx(input$annotation_file$datapath, 2, colNames = TRUE))
    } else if (input$use_example_data == TRUE){
      annot <- data.frame(openxlsx::read.xlsx("example_data/example_annotation.xlsx", 2, colNames = TRUE))
    }
  })
  
  annotation <- reactive({
    req(isTruthy(input$group) && (isTruthy(input$annotation_file) || isTruthy(input$use_example_data)))
    formatAnnotation(annotation_initial(),
                     group = input$group)
  })
  
  annot_columns_pre_group <- reactive ({
    annot_columns <- make.names(colnames(annotation_initial()))
    
    annot_names <- annotation_files()
    
    annot_names <- annot_names[,-(3:4)]
    annot_names <- as.list(make.names(annot_names$Name))
    
    annot_columns <- annot_columns[!annot_columns == "SampleName"]
    annot_columns <- annot_columns[!annot_columns == annot_names]
    
    return(annot_columns)
  })
  
  annot_columns <- reactive({
    annot_columns <- annot_columns_pre_group()
    annot_columns <- annot_columns[!annot_columns == input$group]
    
    return(annot_columns)
  })
  
  type <- reactive({
    files <- list()
    
    if (input$use_example_data == FALSE){
      for (i in 1:length(input$data_file[, 1])){
        files[[i]] <- input$data_file[[i, 'name']]
      }
      
    } else if (input$use_example_data == TRUE){
      files[[1]] <- "example_proteingroups"
      files[[2]] <- "example_phosphosites"
      files[[3]] <- "example_openms_annotated_output_pos"
      files[[4]] <- "example_openms_annotated_output_neg"
      
    }
    
    annot <- annotation_files()
    annot <- annot[,-c(3:4)]
    
    type <- c()
    for (i in 1:length(files)){
      row <- annot[grepl(make.names(files[[i]]), make.names(annot$File)), ]
      type[i] <- row$Name
    }
    
    type <- make.names(type)
    return(type)
  })
  
  data_format <- reactive({
    files <- list()
    
    if (input$use_example_data == FALSE){
      annot <- data.frame(openxlsx::read.xlsx(input$annotation_file$datapath, 1, colNames = FALSE))
      
      for (i in 1:length(input$data_file[, 1])){
        files[[i]] <- input$data_file[[i, 'name']]
      }
      
    } else if (input$use_example_data == TRUE){
      annot <- data.frame(openxlsx::read.xlsx("example_data/example_annotation.xlsx", 1, colNames = FALSE))
      files[[1]] <- "example_proteingroups"
      files[[2]] <- "example_phosphosites"
      files[[3]] <- "example_openms_annotated_output_pos"
      files[[4]] <- "example_openms_annotated_output_neg"
      
    }
    
    annot <- annot[-(1:10),]
    annot <- annot[,-c(2,4)]
    colnames(annot) <- c("File", "Format")
    annot <- na.omit(annot)
    
    data_format <- c()
    for (i in 1:length(files)){
      row <- annot[grepl(make.names(files[[i]]), make.names(annot$File)), ]
      data_format[i] <- row$Format
    }
    
    return(data_format)
  })
  
  ### RENDER UI ----------------------------------------------------------------
  output$select_group <- renderUI({
    req(isTruthy(input$annotation_file) || isTruthy(input$use_example_data))
    
    tagList(
      hr(),
      h6("Group"),
      selectInput("group",
                  label = "Select column(s) to form groups by",
                  choices = annot_columns_pre_group(),
                  multiple = TRUE,
                  selected = annot_columns_pre_group()[1]),
      helpText("Multiple options can be selected to group by a combination of factors."),
      br()
    )
  })
  
  output$annotationFile <- renderUI({
    req(isTruthy(input$annotation_file) || isTruthy(input$use_example_data))
    h5("Annotation File")
  })
  
  output$typeHeader <- renderUI({
    req((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data))
    h6("Dataset Name(s) & Format(s)")
  })
  
  ### PRINT/TABLES/PLOTS -------------------------------------------------------
  output$summary_annotation <- renderPrint({
    utils::str(annotation_initial())
  })
  
  output$annotation_table <- renderTable({
    req(input$annotation_file)
    input$annotation_file
  })
  
  output$type <- renderPrint({
    req((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data))
    
    if (!all(input$data_file$name %in% annotation_files()$File)) {
      cat("A file name provided does not match any uploaded file names. Please check the provided file names against the uploaded file names. \n")
      cat("\t Annotation File Names: ", annotation_files()$File, "\n",
          "\t Uploaded File Names: ", input$data_file$name)
    } else {
      for (i in 1:length(type())){
        
        cat(type()[i], "; ", data_format()[i], "\n", sep = "")
        
        annot_columns <- annotation_initial()[,type()[i]]
        
        missing <- c()
        for (j in 1:length(annot_columns)){
          if (!(make.names(annot_columns[j]) %in% colnames(data()[[i]]))){
            missing <- c(missing, annot_columns[j])
          }
        }
        
        if (data_format()[i] == "PhosphoSites" && length(missing) > 0){
          remove <- c()
          for (k in 1:length(missing)){
            if (grepl(make.names(missing[k]), paste0(colnames(data()[[i]]), collapse = "; "))){
              remove <- c(remove, missing[k])
            }
          }
          missing <- missing[!(missing %in% unique(remove))]
        }
        
        missing <- na.omit(missing)
        
        if (length(missing) > 0){
          cat("\t Columns from annotation file not found in data file: \n\t -", paste0(missing, collapse = "\n\t - "), "\n")
        } else {
          cat("\t All columns from annotation file matched to data file.\n")
        }
        cat("\n")
      }
    }

  })
  
  output$annotation_file_datatable <- renderTable({
    annotation()
  })
  
  output$view_annotation <- DT::renderDataTable({
    shiny::validate(need(input$group, message = "No group column selected. Please select a column under the 'group' header on the left."))
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$group))
    
    DT::datatable(annotation(), rownames = FALSE, options = list(scrollX = TRUE))
    
  })
  
  ## EXPRESSION SET ############################################################
  
  ### REACTIVE OBJECTS ---------------------------------------------------------
  eset_prenorm <- reactive({
    eset_obj <- list()
    
    for (i in 1:length(type())){
      eset_obj[[i]] <- makeEset(data()[[i]], 
                                annotation(), 
                                type = type()[i],
                                species = input$species,
                                data_format = data_format()[i],
                                log_transform = input$log_transform,
                                uniprot_annotation = input$uniprot_annotation)
    }
    
    assign("eset_prenorm", eset_obj, envir = .GlobalEnv)
    return(eset_obj)
  })
  
  eset <- reactive({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$group))
    
    eset <- list()
    
    withProgress(message = "Starting data cleaning...", value = 0, {
      for (i in 1:length(type())){
        incProgress(amount = (1/(length(type()))), detail = paste("Filtering, normalizing, and annotating dataset", i , "of", length(type())))
        
        eset[[i]] <- intensityNorm(eset = eset_prenorm()[[i]],
                                   norm = input$norm_eset,
                                   type = type()[i],
                                   norm_by_batches = input$batch_norm,
                                   batch_column = input$batch_column,
                                   zero_cutoff = input$zero_cutoff,
                                   impute_missing = input$impute_missing,
                                   IRS_column = input$TMT_group_column,
                                   outlier_cols = input$outlier_columns)
      }
    })
    # assign("table_eset", eset, envir = .GlobalEnv)
    return(eset)
  })
  
  ### RENDER UI ----------------------------------------------------------------
  
  output$eset_choose_dataset <- renderUI({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$group))
    
    options <- type()
    selectInput("eset_dataset",
                label = "Choose data to view expression set",
                choices = options)
  })
  
  output$IRS_column <- renderUI({
    if (input$norm_eset == "IRS"){
      shiny::validate(need(input$group, message = "No group column selected. Please select a column under the 'group' header, then you will be able to select the TMT grouping column for IRS normalization."))
      
      selectInput("TMT_group_column",
                  label = "Select which column contains TMT groups",
                  choices = annot_columns())
    }
  })
  
  output$IRS_explanation <- renderUI({
    req(input$norm_eset == "IRS")
    
    tagList(
      br(),
      helpText("Internal Reference Standard Normalization is meant to be applied when multiple TMT experiments are being combined
             in one analysis. This method is adapted from the techniques first described in this publication: Mol Cell Proteomics. 
             2017 May;16(5):873-890. doi: 10.1074/mcp.M116.065524."),
      br(),
      br()
    )
  })
  
  output$batch_column <- renderUI({
    req(input$batch_norm)
    shiny::validate(need(input$group, message = "No group column selected. Please select a column under the 'group' header, then you will be able to select the TMT grouping column for IRS normalization."))
    
    selectInput("batch_column",
                label = "Select which column contains batch information",
                choices = annot_columns())
  }) %>% bindEvent(input$batch_norm)
  
  output$outlier_column_selection <- renderUI({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$group) && isTruthy(input$remove_outlier))
    
    options <- c()
    for (i in 1:length(type())){
      type_options <- paste(type()[i], annotation()$SampleName[!is.na(annotation()[, type()[i]])], sep = "_")
      options <- c(options, type_options)
    }
    
    tagList(
      selectInput("outlier_columns",
                  label = "Select columns to remove",
                  choices = options,
                  multiple = TRUE),
      
      helpText("Options for removal are separated by dataset. The naming convention 
               is 'DatasetName_ColumnName'.")
    )
  })
  
  ### PRINT/TABLE/PLOTS --------------------------------------------------------
  output$view_raw_eset <- DT::renderDataTable({
    shiny::validate(need(input$group, message = "No group column selected. Please select a column under the 'group' header on the left."))
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$group) && isTruthy(input$eset_dataset))
    
    i <- which(type() == input$eset_dataset)
    exprs_eset = Biobase::exprs(eset_prenorm()[[i]])
    
    DT::datatable(exprs_eset, rownames = TRUE, options = list(pageLength = 15, scrollX = TRUE))  %>% 
      DT::formatRound(columns = c(1:ncol(exprs_eset)), digits = 3)
    
  })
  
  output$view_cleaned_eset <- DT::renderDataTable({
    shiny::validate(need(input$group, message = "No group column selected. Please select a column under the 'group' header on the left."))
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$group) && isTruthy(input$normplot_dataset))
    
    i <- which(type() == input$normplot_dataset)
    exprs_eset = Biobase::exprs(eset()[[i]])
    
    assign("table_eset_table", exprs_eset, envir = .GlobalEnv)
    
    DT::datatable(exprs_eset, rownames = TRUE, options = list(pageLength = 15, scrollX = TRUE))  %>% 
      DT::formatRound(columns = c(1:ncol(exprs_eset)), digits = 3)
    
  })
  
  ### DOWNLOADS ----------------------------------------------------------------
  output$eset_table_download <- downloadHandler(
    filename = function() {
      i <- which(type() == input$eset_dataset)
      paste(type()[i], "_RawEsetTable_", Sys.Date(), ".txt", sep = "")
    },
    content = function(file){
      i <- which(type() == input$eset_dataset)
      exprs_eset = Biobase::exprs(eset_prenorm()[[i]])
      write.table(exprs_eset, file, row.names = TRUE, sep = "\t", col.names = NA)
    }
  )
  
  output$processed_eset_download <- downloadHandler(
    filename = function() {
      i <- which(type() == input$normplot_dataset)
      paste(type()[i], "_ProcessedEsetTable_", Sys.Date(), ".txt", sep = "")
    },
    content = function(file){
      i <- which(type() == input$normplot_dataset)
      exprs_eset = Biobase::exprs(eset()[[i]])
      # assign("table_eset_download", exprs_eset, envir = .GlobalEnv)
      write.table(exprs_eset, file, row.names = TRUE, sep = "\t", col.names = NA)
    }
  )
  
  ## PRE-NORM DATA SUMMARY #####################################################
  
  ### REACTIVE OBJECTS ---------------------------------------------------------
  
  visualizing_missing_values_plot <- reactive({
    i <- which(type() == input$missingvalue_dataset)
    data <- data.frame(Biobase::exprs(eset_prenorm()[[i]]))
    data[data == 0] <- NA
    
    VIM_eset <- VIM::aggr(data, cex.axis = 0.7, numbers = TRUE, prop = FALSE)
    
    return(VIM_eset)
  })
  
  distribution_plot <- reactive({
    i <- which(type() == input$distribution_dataset)
    data <- Biobase::exprs(eset_prenorm()[[i]])
    
    df <- data.frame(reshape2::melt(data, id.vars = NULL));
    colnames(df) <- c("Feature","Sample", "Intensity");
    
    plot <- ggplot(df, aes(x = Sample, y = Intensity, fill = Sample)) + 
      geom_boxplot() +
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
      labs(title = paste(type()[i], "Data Distribution Boxplot"))
    
    return(plot)
  })
  
  prenorm_stat_summary <- reactive({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$group) && isTruthy(input$statsummary_dataset))
    
    i <- which(type() == input$statsummary_dataset)
    data <- data.frame(Biobase::exprs(eset_prenorm()[[i]]))
    table <- psych::describe(data, na.rm = FALSE)
    
    return(table)
  })
  
  fdr_calculation <- reactive({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$group) && isTruthy(input$statsummary_dataset))
    
    i <- which(type() == input$statsummary_dataset)
    data <- data()[[i]]
    assign("data", data, envir = .GlobalEnv)
    
    FDR <- "Insufficient information provided."
    
    if ("Reverse" %in% colnames(data)){
      decoy <- nrow(data[(data$Reverse == "+" & data$Potential.contaminant == ""),])
      non_decoy <- nrow(data[(data$Reverse == "" & data$Potential.contaminant == ""),])
      FDR <- decoy / non_decoy
      FDR <- paste0("FDR: ", round(FDR * 100, digits = 2), "%; Decoys: ", decoy, " & Non-Decoy Targets: ", non_decoy)
    } else {
      col_names <- c("Index", "Protein", "Accession", "Protein.IDs")
      column <- intersect(colnames(data), col_names)
      column <- ifelse(length(column) > 1, column[1], column)
      try({
        decoy <- nrow(data[grepl("REV_|reverse_", data[, column]) & !grepl("CON_|contam_", data[, column]),])
        non_decoy <- nrow(data[!grepl("REV_|reverse_", data[, column]) & !grepl("CON_|contam_", data[, column]),])
      })
      if(length(column) != 0){
        FDR <- decoy / non_decoy
        FDR <- paste0("FDR: ", round(FDR * 100, digits = 2), "%; Decoys: ", decoy, " & Non-Decoy Targets: ", non_decoy)
      }  
    }
    
    return(FDR)
  })
  
  ### RENDER UI ----------------------------------------------------------------
  output$statsummary_choose_dataset <- renderUI({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)))
    options <- type()
    selectInput("statsummary_dataset",
                label = "Choose dataset to visualize",
                choices = options)
  })
  
  output$missingvalue_choose_dataset <- renderUI({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)))
    options <- type()
    selectInput("missingvalue_dataset",
                label = "Choose dataset to visualize",
                choices = options)
  })
  
  output$distribution_choose_dataset <- renderUI({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)))
    options <- type()
    selectInput("distribution_dataset",
                label = "Choose dataset to visualize",
                choices = options)
  })
  
  ### PRINT/TABLE/PLOTS --------------------------------------------------------
  
  output$prenorm_describe <- DT::renderDataTable({
    table <- prenorm_stat_summary()
    
    DT::datatable(table, rownames = TRUE, options = list(scrollX = TRUE)) %>%
      DT::formatRound(columns = c(3:ncol(table)), digits = 3)
    
  })
  
  output$fdr <- renderPrint({
    fdr_calculation()
  })
  
  output$title_missing_value <- renderPrint({
    i <- which(type() == input$missingvalue_dataset)
    
    cat(type()[i], " Dataset Missing Value Plots")
  })
  
  output$visualizing_missing_values <- renderPlot({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$group) && isTruthy(input$missingvalue_dataset))
    
    visualizing_missing_values_plot()
  })
  
  output$distribution_plot <- renderPlot({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$group) && isTruthy(input$distribution_dataset))
    
    distribution_plot()
  })
  
  ### DOWNLOADS ----------------------------------------------------------------
  output$prenorm_table_download <- downloadHandler(
      filename = function() {
        i <- which(type() == input$statsummary_dataset)
        paste(type()[i], '_RawDataSummary_', Sys.Date(), '.txt', sep = '')
      },
      content = function(file) {
        i <- which(type() == input$statsummary_dataset)
        data <- data.frame(Biobase::exprs(eset_prenorm()[[i]]))
        table <- psych::describe(data, na.rm = FALSE)
        write.table(table, file, sep = "\t", col.names = NA)
      }
  )
  
  output$prenorm_missingvalue_download <- downloadHandler(
    filename = function() {
      i <- which(type() == input$missingvalue_dataset)
      paste(type()[i], '_RawDataMissingValues_', Sys.Date(), input$prenorm_missingvalue_ext, sep = '')
    },
    content = function(file) {
      i <- which(type() == input$missingvalue_dataset)
      if (input$prenorm_missingvalue_ext == ".png"){
        grDevices::png(file, width = 1200, height = 1000, units = "px", res = 150)
        plot(visualizing_missing_values_plot(), cex.axis = 0.7, numbers = TRUE, prop = FALSE)
        dev.off()
        
      } else if (input$prenorm_missingvalue_ext == ".pdf"){
        grDevices::pdf(file, width = 12, height = 10)
        plot(visualizing_missing_values_plot(), cex.axis = 0.7, numbers = TRUE, prop = FALSE)
        dev.off()
        
      } else if (input$prenorm_missingvalue_ext == ".svg"){
        svglite::svglite(file, width = 12, height = 10)
        plot(visualizing_missing_values_plot(), cex.axis = 0.7, numbers = TRUE, prop = FALSE)
        dev.off()
      }
    }
  )
  
  output$prenorm_distribution_download <- downloadHandler(
    filename = function() {
      i <- which(type() == input$distribution_dataset)
      paste(type()[i], '_RawDataBoxplot_', Sys.Date(), input$prenorm_distribution_ext, sep = '')
    },
    content = function(file) {
      ggsave(file, plot = distribution_plot(), width = 12, height = 7)
    }
  )
  
  ## NORMALIZATION #############################################################
  
  ### REACTIVE OBJECTS ---------------------------------------------------------
  normal_plot <- reactive({
    i <- which(type() == input$normplot_dataset)
    
    plot <- intensityNormPlots(eset_prenorm = eset_prenorm()[[i]],
                               eset_postnorm = eset()[[i]],
                               norm = input$norm_eset, 
                               type = type()[i], 
                               plottype = input$norm_plottype,
                               ma_array = input$ma_array,
                               qq_column = input$QQ_column,
                               zero_cutoff = input$zero_cutoff)
    plot
  }) %>% bindEvent(input$normal_button)
  
  ### RENDER UI ----------------------------------------------------------------
  
  output$normplots_choose_dataset <- renderUI({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)))
    options <- type()
    selectInput("normplot_dataset",
                label = "Choose dataset to visualize",
                choices = options)
  })
  
  output$ma_array_choose <- renderUI({
    req(input$norm_plottype == "MA")
    i <- which(type() == input$normplot_dataset)
    
    list <- 1:ncol(eset()[[i]])
    selectInput("ma_array",
                 label = "Choose column to generate array",
                 choices = list)
  }) %>% bindEvent(input$norm_plottype)
  
  output$qqplot_column <- renderUI({
    req(input$norm_plottype == "QQ")
    i <- which(type() == input$normplot_dataset)
    
    options <- colnames(Biobase::exprs(eset()[[i]]))
    selectInput("QQ_column",
                label = "Choose column to generate QQ plot",
                choices = options)
  }) %>% bindEvent(input$norm_plottype)
  
  output$shapiro_column <- renderUI({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$group) && isTruthy(input$normplot_dataset))
    i <- which(type() == input$normplot_dataset)
    
    tryCatch({options <- colnames(Biobase::exprs(eset()[[i]]))}, 
             error = function(cond){stop("COULD NOT PERFORM SHAPIRO TEST. PLEASE CHECK FOR ERRORS IN DATA INPUT.")})
    
    selectInput("shapiro_column",
                label = "Choose feature to check normality",
                choices = options)
  })
  
  output$shapiro_test <- renderPrint({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$group) && isTruthy(input$normplot_dataset) && isTruthy(input$shapiro_column))
    set.seed(1234) 
    
    i <- which(type() == input$normplot_dataset)
    if (nrow(data.frame(Biobase::exprs(eset()[[i]]))) > 5000){
      data <- dplyr::sample_n(data.frame(exprs(eset()[[i]])), 5000)
    } else {
      data <- data.frame(Biobase::exprs(eset()[[i]]))
    }
    
    shapiro.test(data[, input$shapiro_column])
    
    # ks.test(data[, input$shapiro_column], "pnorm")
  })
  
  ### PRINT/TABLE/PLOTS --------------------------------------------------------
  output$normal_plot <- renderPlot({
    normal_plot()
  })
  
  ### DOWNLOADS ----------------------------------------------------------------
  output$normal_download <- downloadHandler(
    filename = function() {
      i <- which(type() == input$normplot_dataset)
      paste(type()[i], "_", input$norm_plottype, "_NormPlot_", Sys.Date(), input$normal_extension , sep = "")
    },
    content = function(file){
      i <- which(type() == input$normplot_dataset)
      if(input$normal_extension == ".png"){
        png(file, width = 1200, height = 1000, units = "px", res = 150)
        
        intensityNormPlots(eset_prenorm = eset_prenorm()[[i]],
                           eset_postnorm = eset()[[i]],
                           norm = input$norm_eset, 
                           type = type()[i], 
                           plottype = input$norm_plottype,
                           ma_array = input$ma_array,
                           qq_column = input$QQ_column,
                           zero_cutoff = input$zero_cutoff)
        dev.off()
        
      } else if(input$normal_extension == ".pdf"){
        pdf(file, width = 12, height = 7)
        
        intensityNormPlots(eset_prenorm = eset_prenorm()[[i]],
                           eset_postnorm = eset()[[i]],
                           norm = input$norm_eset, 
                           type = type()[i], 
                           plottype = input$norm_plottype,
                           ma_array = input$ma_array,
                           qq_column = input$QQ_column,
                           zero_cutoff = input$zero_cutoff)
        dev.off()
        
      } else if(input$normal_extension == ".svg"){
        svglite::svglite(file, width = 12, height = 7)
        
        intensityNormPlots(eset_prenorm = eset_prenorm()[[i]],
                           eset_postnorm = eset()[[i]],
                           norm = input$norm_eset, 
                           type = type()[i], 
                           plottype = input$norm_plottype,
                           ma_array = input$ma_array,
                           qq_column = input$QQ_column,
                           zero_cutoff = input$zero_cutoff)
        dev.off()
        
      }
    }
  )
  
  ## CORRELATION ###############################################################
  
  ### REACTIVE OBJECTS ---------------------------------------------------------
  base_corr_plot <- reactive({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$group))
    i <- which(type() == input$corrplot_dataset)
    
    corr_plot <- variationPlot(eset()[[i]],
                               high_variance = FALSE,
                               plot = TRUE,
                               number = input$corr_number,
                               neg_color = input$corr_neg_color,
                               pos_color = input$corr_pos_color,
                               title = paste(type()[i], input$corr_title))
    return(corr_plot);
  }) %>% bindEvent(input$corr_button)
  
  var_corr_plot <- reactive({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$group))
    i <- which(type() == input$corrplot_dataset)
    
    corr_plot <- variationPlot(eset()[[i]],
                               high_variance = TRUE,
                               plot = TRUE,
                               number = input$corr_number,
                               percent_var = input$corr_variance,
                               neg_color = input$corr_neg_color,
                               pos_color = input$corr_pos_color,
                               title = paste(type()[i], input$corr_title))
    return(corr_plot);
  }) %>% bindEvent(input$corr_button)
  
  ### RENDER UI ----------------------------------------------------------------
  
  output$corrplot_choose_dataset <- renderUI({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$group))
    
    options <- type()
    selectInput("corrplot_dataset",
                label = "Choose data to view expression set",
                choices = options)
  })
  
  ### PRINT/TABLE/PLOTS --------------------------------------------------------
  output$correlation_plot <- renderPlot({
    base_corr_plot()
  })
  
  output$variance_correlation_plot <- renderPlot({
    var_corr_plot()
  })
  
  ### DOWNLOADS ----------------------------------------------------------------
  output$corr_plot_download <- downloadHandler(
    filename = function() {
      i <- which(type() == input$corrplot_dataset)
      paste(type()[i], "_CorrelationPlot_", Sys.Date(), input$corrplot_extension, sep = "")
    },
    content = function(file){
    ggsave(file, base_corr_plot(), height = 7, width = 7)
    }
  )
  
  output$var_corr_plot_download <- downloadHandler(
    filename = function() {
      i <- which(type() == input$corrplot_dataset)
      paste(type()[i], "_CorrelationPlot_", Sys.Date(), input$var_corrplot_extension, sep = "")
    },
    content = function(file){
      ggsave(file, var_corr_plot(), height = 7, width = 7)
    }
  )
  
  ## PCA #######################################################################
  
  ### REACTIVE OBJECTS ---------------------------------------------------------
  raw_PC_plot <- reactive({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$group) && isTruthy(input$PC_dataset))
    i <- which(type() == input$PC_dataset)
    drawPCA(eset_prenorm()[[i]], 
            x_axis = input$PCA_xaxis, 
            y_axis = input$PCA_yaxis, 
            type = type()[i], 
            color = input$PCA_color, 
            include_densities = input$include_densities, 
            shapes = input$PCA_shapes,
            title_add = paste("Pre-Normalization", input$PCA_title),
            add_labels = input$PCA_labels,
            add_ellipse = input$PCA_ellipse)
  }) %>% bindEvent(c(input$PCA_button))
  
  PC_plot <- reactive({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$group) && isTruthy(input$PC_dataset))
    i <- which(type() == input$PC_dataset)
    drawPCA(eset()[[i]], 
            x_axis = input$PCA_xaxis, 
            y_axis = input$PCA_yaxis, 
            type = type()[i], 
            color = input$PCA_color, 
            include_densities = input$include_densities, 
            shapes = input$PCA_shapes,
            title_add = paste("Post-Normalization", input$PCA_title),
            add_labels = input$PCA_labels,
            add_ellipse = input$PCA_ellipse)
    
  }) %>% bindEvent(c(input$PCA_button))
  
  PC_variance_plot <- reactive({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$group) && isTruthy(input$PC_dataset))
    i <- which(type() == input$PC_dataset)
    drawPCApercentages(eset()[[i]],
                       type()[i])
  })
  
  PC_loadings <- reactive({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$group) && isTruthy(input$PC_dataset))
    i <- which(type() == input$PC_dataset)
    
    data = t(Biobase::exprs(eset()[[i]]))
    PC_data <- stats::prcomp(data)
    loading_scores <- PC_data$rotation
  })
  
  PC_loadings_plot <- reactive({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$group) && isTruthy(input$PC_dataset) && isTruthy(input$PC_loadings_column))
    i <- which(type() == input$PC_dataset)
    drawPCAloadings(eset()[[i]],
                    column = input$PC_loadings_column)
  })
  
  PCs <- reactive ({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$group) && isTruthy(input$PC_dataset))
    
    i <- which(type() == input$PC_dataset)
    data = t(Biobase::exprs(eset()[[i]]))
    PC_data <- stats::prcomp(data)
    percent_variance <- summary(PC_data)$importance["Proportion of Variance",] * 100
    PCs <- rownames(data.frame(percent_variance))
  })
  
  top_PC_percent <- reactive({
    i <- which(type() == input$PC_dataset)
    
    data = t(Biobase::exprs(eset()[[i]]))
    PC_data <- stats::prcomp(data)
    percent_variance <- summary(PC_data)$importance["Proportion of Variance",] * 100
    unname(percent_variance)
    top_two <- percent_variance[1] + percent_variance[2]
    return(top_two)
  })
  
  ### RENDER UI ----------------------------------------------------------------
  output$PCplot_choose_dataset <- renderUI({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$group))
    options <- type()
    selectInput("PC_dataset",
                label = "Choose dataset to visualize",
                choices = options)
  })
  
  output$PC_xaxis <- renderUI({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$group) && isTruthy(input$PC_dataset))
    
    selectInput("PCA_xaxis",
                label = "X-Axis",
                choices = PCs(),
                selected = "PC1")
  })
  
  output$PC_yaxis <- renderUI({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$group) && isTruthy(input$PC_dataset))
    
    selectInput("PCA_yaxis",
                label = "Y-Axis",
                choices = PCs(),
                selected = "PC2")
  })
  
  output$PC_loadings_column <- renderUI({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$group) && isTruthy(input$PC_dataset))
    
    columns <- PCs()
    selectInput("PC_loadings_column",
                label = "Select principle component to view top loadings",
                choices = columns,
                multiple = FALSE)
  })
  
  ### PRINT/TABLE/PLOTS --------------------------------------------------------
  output$top_PC_print <- renderPrint({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$group) && isTruthy(input$PC_dataset))
        
    cat("The above plot suggests that dimensionality reduction using PC1 & PC2 can capture ≥", 
        top_PC_percent(), "% of the dataset’s variance.")
  })
  
  output$PC_plot <- renderPlot({
    i <- which(type() == input$PC_dataset)
    gridExtra::grid.arrange(raw_PC_plot() + theme(legend.position = "bottom"), 
                            PC_plot() + theme(legend.position = "bottom"),
                            top = paste0("PCA Plots ", type()[i]), 
                            ncol = 2);
  })
  
  output$PC_variance_plot <- renderPlot({
    PC_variance_plot()
  })
  
  output$PC_loadings_plot <- renderPlot({
    PC_loadings_plot()
  })
  
  ### DOWNLOADS ----------------------------------------------------------------
  
  output$PCA_download <- downloadHandler(
    filename = function() {
      i <- which(type() == input$PC_dataset)
      paste(type()[i], "_PCA_", Sys.Date(), input$PCA_extension , sep = "")
    },
    content = function(file){
      plots <- gridExtra::arrangeGrob(raw_PC_plot() + theme(legend.position = "bottom"), 
                                      PC_plot() + theme(legend.position = "bottom"), 
                                      nrow = 1) 
      ggplot2::ggsave(file, plots, width = 11, height = 6)
    }
  )
  
  output$PC_loadings_download <- downloadHandler(
    filename = function() {
      paste("PCloadings_", type(), "_", Sys.Date(), ".txt", sep = "")
    },
    content = function(file){
      write.table(PC_loadings(), file, sep = "\t", col.names = NA)
    }
  )
  
  ## UMAP ######################################################################
  
  ### REACTIVE OBJECTS ---------------------------------------------------------
  UMAP_plot <- reactive({
    req((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data) && isTruthy(input$UMAP_dataset))
    
    i <- which(type() == input$UMAP_dataset)
    drawUMAP(eset()[[i]], 
             type = type()[i], 
             color = input$UMAP_color,
             title_add = input$UMAP_title,
             shapes = input$UMAP_shapes,
             add_labels = input$UMAP_labels)
  }) %>% bindEvent(input$UMAP_button)
  
  ### RENDER UI ----------------------------------------------------------------
  
  output$UMAPplot_choose_dataset <- renderUI({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$group))
    options <- type()
    selectInput("UMAP_dataset",
                label = "Choose dataset to visualize",
                choices = options)
  })
  
  ### PRINT/TABLE/PLOTS --------------------------------------------------------
  output$UMAP_plot <- renderPlot({
    UMAP_plot()
  })
  
  ### DOWNLOADS ----------------------------------------------------------------
  output$UMAP_download <- downloadHandler(
    filename = function() {
      i <- which(type() == input$UMAP_dataset)
      paste(type()[i], "_UMAP_", Sys.Date(), input$UMAP_extension, sep = "")
    },
    content = function(file){
      ggplot2::ggsave(file, UMAP_plot(), width = 7, height = 7)
    }
  )
  
  ## HEATMAP ###################################################################
  
  ### REACTIVE OBJECTS ---------------------------------------------------------
  heatmap_plot <- reactive({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$subset))
    i <- which(type() == input$heatmap_dataset)
    
    drawHeatmaps(eset()[[i]], 
                 type = type()[i], 
                 mapcolor = input$heatmap_color, 
                 group_color = input$heatmap_group_color,
                 subset = input$subset,
                 subset_genes = input$gene_list,
                 subset_genes_file = gene_list_file(),
                 subset_features = input$features_list, 
                 subset_variable = input$variable_list,
                 kclustering = input$kclustering,
                 log2 = input$log_transform,
                 cluster_samples = input$cluster_samples, 
                 title_add = input$heatmap_title,
                 annot = annotation(),
                 zscore = input$zscore)
  }) %>% bindEvent(input$heatmap_button)
  
  gene_list_file <- reactive({
    req(input$gene_list_file)
    gene_list <- utils::read.table(input$gene_list_file$datapath, header = FALSE, sep = "\t")
    gene_list <- as.list(gene_list)
    gene_list$V1
  })
  
  ### RENDER UI ----------------------------------------------------------------
  
  output$heatmap_choose_dataset <- renderUI({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$group))
    options <- type()
    selectInput("heatmap_dataset",
                label = "Choose dataset to visualize",
                choices = options)
  })
  

  output$select_features <- renderUI({
    req((input$subset) == "subset_features")
    i <- which(type() == input$heatmap_dataset)
    
    type <- type()[i]
    annotation <- annotation()[,c("SampleName", type)]
    features <- annotation[stats::complete.cases(annotation),]
    features <- features$SampleName
    selectInput("features_list",
                label = "Select Features to Include",
                choices = features,
                multiple = TRUE)
  }) %>% bindEvent(input$subset, input$heatmap_dataset)

  

  output$gene_list_select <- renderUI({
    req((input$subset) == "subset_genes")
    i <- which(type() == input$heatmap_dataset)
    
    eset_genes <- sub("_.*", "", row.names(Biobase::exprs(eset()[[i]])))
    br()
    
    shiny::validate(need(try(eset_genes), "No gene names available."))
    
    selectizeInput("gene_list", 
                   label = "Select genes by typing their names into the text box or
                by selecting them from the dropdown",
                   eset_genes,
                   multiple = TRUE
    )
  }) %>% bindEvent(c(input$subset, input$heatmap_dataset))

  
  output$gene_list_file <- renderUI({
    req((input$subset) == "subset_genes")
    br()
    
    fileInput("gene_list_file",
              label = "OR Upload a tab-delimited file containing a column of gene names with no column header",
              accept = c("text",
                         "text/tab-separated-values,text/plain",
                         ".tsv")
    )
  }) %>% bindEvent(input$subset)

  

  output$select_variable <- renderUI({
    req((input$subset) == "subset_variable")
    numericInput("variable_list",
                 label = "Select number of genes with highest variance to include",
                 value = 50,
                 min = 1,
                 max = 1000,
                 step = 10)
  }) %>% bindEvent(input$subset)

  
  ### PRINT/TABLE/PLOTS --------------------------------------------------------
  output$heatmap_plot <- renderPlot({
    heatmap_plot()
  }, height = 800)
  
  ### DOWNLOADS ----------------------------------------------------------------
  output$heatmap_download <- downloadHandler(
    filename = function() {
      i <- which(type() == input$heatmap_dataset)
      paste(type()[i], "_Heatmap_", Sys.Date(), input$heatmap_extension, sep = "")
    },
    
    content = function(file){
      withProgress(message = "Downloading heatmap plot:", detail = "This may take a moment...",value = 0.3, {
        if(input$heatmap_extension == ".pdf"){
          grDevices::pdf(file, width = 10, height = 10)
          ComplexHeatmap::draw(heatmap_plot())
          grDevices::dev.off()
          setProgress(value = 0.9)
          
        } else if(input$heatmap_extension == ".png"){
          grDevices::png(file, width = 1000, height = 1000)
          ComplexHeatmap::draw(heatmap_plot())
          grDevices::dev.off()
          setProgress(value = 0.9)
            
        } else if(input$heatmap_extension == ".svg"){
          svglite::svglite(file, width = 10, height = 10)
          ComplexHeatmap::draw(heatmap_plot())
          grDevices::dev.off()
          setProgress(value = 0.9)
          
        }
      })
    }
  )
  
  ## MD PLOT ###################################################################
  
  ### REACTIVE OBJECTS ---------------------------------------------------------
  md_plot <- reactive({
    i <- which(type() == input$MD_dataset)
    drawMD(annot = annotation(),
           eset = eset()[[i]],
           type = type()[i],
           fc_cutoff = input$md_fc_cutoff,
           up_color = input$md_up_color,
           down_color = input$md_down_color,
           title_add = input$md_title,
           add_labels = input$md_add_labels,
           label_specific = input$md_label_specific,
           label_specific_gene = input$md_label_specific_gene,
           label_number = input$md_label_number,
           return_logfc_index = FALSE,
           logfc_index_choice = input$md_contrast)
  }) %>% bindEvent(input$md_button)
  
  md_contrasts <- reactive({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$group))
    i <- which(type() == input$MD_dataset)
    
    contrasts <- drawMD(annot = annotation(),
                        eset = eset()[[i]],
                        type = type()[i],
                        return_logfc_index = TRUE,
                        logfc_index_choice = NULL)
  })

  ### RENDER UI ----------------------------------------------------------------  
  output$MDplot_choose_dataset <- renderUI({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$group))
    options <- type()
    selectInput("MD_dataset",
                label = "Choose dataset to visualize",
                choices = options)
  })
  
  output$md_contrast <- renderUI({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$MD_dataset))
    selectInput("md_contrast",
                label = "Choose contrast to generate MD plot",
                choices = md_contrasts())
  })
  
  output$md_label_number <- renderUI({
    req(input$md_add_labels == TRUE)
    
    numericInput("md_label_number",
                 label = "Select number of genes to label",
                 value = 20)
  }) %>% bindEvent(input$md_add_labels)
  
  output$md_label_specific_gene <- renderUI({
    req(input$md_label_specific == TRUE)
    i <- which(type() == input$MD_dataset)
    
    genes <- sub("_.*", "", rownames(Biobase::exprs(eset()[[i]])))
    
    selectInput("md_label_specific_gene",
              label = "Type gene name",
              multiple = TRUE,
              choices = genes)
  }) %>% bindEvent(input$md_label_specific)
  
  ### PRINT/TABLE/PLOTS --------------------------------------------------------
  output$md_plot <- renderPlot({
    md_plot()
  }, height = 800)
  
  ### DOWNLOADS ----------------------------------------------------------------
  output$md_download <- downloadHandler(
    filename = function() {
      i <- which(type() == input$MD_dataset)
      paste(type()[i], "_MDPlot_", Sys.Date(), input$md_extension, sep = "")
    },
    content = function(file){
      ggsave(file, md_plot(), width = 10)
    }
  )
  
  ## LIMMA LINEAR MODELING #####################################################
  
  ### REACTIVE OBJECTS ---------------------------------------------------------
  linear_model <- reactive({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$group) && isTruthy(input$limma_button))
    
    linear_models <- list()
    
    for (i in 1:length(type())){
      linear_models[[i]] <- limmaLM(annot = annotation(),
                                    eset = eset()[[i]],
                                    samples = input$linear_factors,
                                    time_series = input$include_timeseries,
                                    time_col = input$time_col,
                                    time_points_cont = input$time_points_cont,
                                    time_points_disc = input$time_points_disc,
                                    time = input$time_type,
                                    contrast_fit = input$contrast_fit,
                                    contrasts_subset = input$limma_contrasts,
                                    covariate = input$add_covariate,
                                    covariate_col = input$covariate_col)
    }
    
    return(linear_models)
    
  }) %>% bindEvent(input$limma_button)
  
  coef_options <- reactive({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$group))
    
    coef <- input$limma_contrasts
  })
  
  linear_model_table <- reactive({
    req(input$coef_options)
    i <- which(type() == input$limma_dataset)
    
    makeTopTable(linear_model()[[i]],
                 coef = input$coef_options,
                 coef_options = input$limma_contrasts,
                 logfc_th = 0)
    
  })
  
  linear_model_results_table <- reactive({
    req(input$results_pval, input$results_logfc)
    
    i <- which(type() == input$limma_dataset)
    fit <- linear_model()[[i]]
    
    results <- data.frame(summary(limma::decideTests(fit,
                                                     p.value = input$results_pval,
                                                     lfc = input$results_logfc)))
    results$Var1 <- as.character(results$Var1)
    results$Var2 <- as.character(results$Var2)
    results$Freq <- as.numeric(results$Freq)
    
    table <- data.frame()
    
    for (row in 1:nrow(results)) {
      rowname <- results[row, 1]
      colname <- results[row, 2]
      table[rowname, colname] <- results[row , 3]
    }
    
    return(table)
    
  })
  
  pval_distribution_plot <- reactive({
    pval <- plotPvalHistogram(top_table = linear_model_table(),
                              pval_type = "P.Value",
                              threshold = 0.05)
    
    adj_pval <- plotPvalHistogram(top_table = linear_model_table(),
                                  pval_type = "adj.P.Val",
                                  threshold = 0.05)
    
    gridExtra::grid.arrange(pval + theme(legend.position = "bottom"), 
                            adj_pval + theme(legend.position = "bottom"), 
                            nrow = 1);
  })
  
  ### RENDER UI ----------------------------------------------------------------
  output$limma_choose_dataset <- renderUI({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$group))
    options <- type()
    selectInput("limma_dataset",
                label = "Choose dataset to visualize",
                choices = options)
  })
  
  output$linear_model_results_header <- renderUI({
    req(input$limma_button)
    tagList(
      h5("Model Results Overview"),
      
      helpText("To be considered significant in the model results table a row must have an adjusted 
                p-value less than or equal to the specified p-value cutoff and a logFC greater than or 
                equal to the specified logFC cutoff. The method for adjusting p-values is the 
                Benjamini-Hochberg procedure which uses the desired false discovery rate control."),
    )
  })
  
  output$linear_model_results_logfc_cutoff <- renderUI({
    req(input$limma_button)
    numericInput("results_logfc",
                 label = "LogFC cutoff",
                 value = 1)
  })
  
  output$linear_model_results_pval_cutoff <- renderUI({
    req(input$limma_button)
    numericInput("results_pval",
                 label = "Adj. P-Value cutoff",
                 value = 0.05)
  })
  
  output$linear_model_equation_header <- renderUI({
    req(input$limma_button)
    
    tagList(h5("Model Equation"))
  })
  
  output$linear_factors <- renderUI({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$group) && isTruthy(input$limma_dataset))
    
    i <- which(type() == input$limma_dataset)
    options <- unique(pData(eset()[[i]])$Group)
    
    selectInput("linear_factors",
                label = "Choose at least two groups for comparison, i.e. Treatment & Control",
                choices = options,
                selected = options,
                multiple = TRUE)
  })
  
  output$choose_contrasts <- renderUI({
    req((input$contrast_fit == TRUE), isTruthy(length(input$linear_factors) > 1))
    
    i <- which(type() == input$limma_dataset)
    contrasts <- limmaLM(annot = annotation(),
                         eset = eset()[[i]],
                         samples = input$linear_factors,
                         time_series = input$include_timeseries,
                         time_col = input$time_col,
                         time_points_cont = input$time_points_cont,
                         time_points_disc = input$time_points_disc,
                         time = input$time_type,
                         contrast_fit = input$contrast_fit,
                         covariate = input$add_covariate,
                         covariate_col = input$covariate_col,
                         return_contrasts = TRUE)
    
    tagList(
      selectInput("limma_contrasts",
                  label = "Choose contrasts to include in model",
                  choices = contrasts,
                  multiple = TRUE),
    )
  })
  
  output$choose_covariate <- renderUI({
    req(input$add_covariate == TRUE)
    
    selectizeInput("covariate_col",
                   label = "Select column to act as covariate",
                   choices = annot_columns(),
                   multiple = TRUE,
                   options = list(maxItems = 3))
    
  })
  
  output$choose_time <- renderUI({
    req(input$include_timeseries == TRUE)
    
    selectInput("time_col",
                 label = "Select column containting time series information",
                 choices = annot_columns(),
                 multiple = FALSE)
  })
  
  output$coef_options <- renderUI({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$limma_button))
    
    selectInput("coef_options",
                label = "Choose contrast for fit table",
                choices = coef_options(),
                multiple = FALSE)
  }) %>% bindEvent(input$limma_button)
  
  output$time_type <- renderUI({
    req(input$include_timeseries == TRUE)
    
    radioButtons("time_type",
                 label = "Select how to include the time variable in the model",
                 choices = c("Continuous" = "continuous",
                             "Discrete" = "discrete"))
  }) %>% bindEvent(input$include_timeseries)
  
  output$time_points <- renderUI({
    req(isTruthy(input$include_timeseries))
    
    i <- which(type() == input$limma_dataset)
    if(req(input$time_type) == "continuous"){
      number = length(unique(pData(eset()[[i]])[,input$time_col]))
      numericInput("time_points_cont",
                   label = "Choose number of time points to include in the model",
                   min = 0,
                   max = number,
                   value = 3) 
      
    } else if(req(input$time_type) == "discrete"){
      options <- unique(pData(eset()[[i]])[,input$time_col])
      options <- options[-order(options)[1]]
      options <- sort(options)
      selectInput("time_points_disc",
                  label = "Choose which time points to include in the model",
                  choices = options,
                  multiple = TRUE)
      
    }
    
  }) %>% bindEvent(input$time_type, input$include_timeseries, input$time_col)
  
  output$time_df_explain <- renderUI({
    req(input$include_timeseries == TRUE)
    
    tagList(helpText("When choosing the number of time points to include in the model, keep in mind that each added
             time point will decrease the overall residual degrees of freedom of the model. Once the residual
             DF reaches zero, the model will no longer run."),
            
            br(),
            br(),
            
            helpText("The first time point (usually 0) will be used as a base/reference automatically and so will not be an option for contrasts."),
            
            br(),
            br()
    )
    
  })
  
  output$limma_toptable_header <- renderUI({
    req(input$limma_button)
    tagList(
      h5("Limma Top Table"),
      helpText("To change the contrast used to generate the table or to threshold by logFC, use the selections on the right.")
    )
  })
  
  ### PRINT/TABLE/PLOTS --------------------------------------------------------
  
  output$model_equation <- renderPrint({
    req(input$limma_button)
    limmaEquation(groups = input$group,
                  include_covariates = input$add_covariate,
                  covariate_col = input$covariate_col,
                  time_series = input$include_timeseries)
  }) %>% bindEvent(input$limma_button)
  
  output$linear_model_table <- DT::renderDataTable({
    req(input$coef_options)
    
    table <- linear_model_table()
    rownames(table) <- 1:nrow(table)
    table <- table[order(-abs(table[['logFC']])),]
    
    DT::datatable(table, rownames = FALSE) %>% 
      DT::formatSignif(columns = 2:ncol(table), digits = 3) %>%
      DT::formatStyle('logFC', fontWeight = 'bold', background = GeneTonic::styleColorBar_divergent(table$logFC, 'lightskyblue', 'lightcoral')) %>%
      DT::formatStyle('adj.P.Val', backgroundColor = DT::styleInterval(0.05, c('darkseagreen', 'darksalmon')))
    
  }) %>% bindEvent(c(input$coef_options, input$limma_button))
  
  output$linear_model_results <- renderTable({
    req(input$results_pval, input$results_logfc)
    linear_model_results_table()
    
  }, rownames = TRUE, digits = 0, hover = TRUE, bordered = TRUE, striped = TRUE, width = '100%', align = 'l') %>% bindEvent(input$results_pval, input$results_logfc, input$limma_button)
  
  output$pval_distribution_plot <- renderPlot({
    req(input$coef_options)
    pval_distribution_plot()
  }) %>% bindEvent(c(input$coef_options, input$limma_button))
  
  ### DOWNLOADS ----------------------------------------------------------------
  
  output$linear_table_download <- downloadHandler(
    filename = function() {
      i <- which(type() == input$limma_dataset)
      paste(type()[i], "_LimmaTopTable_", Sys.Date(), ".txt", sep = "")
    },
    content = function(file){
      write.table(linear_model_table(), file, row.names = FALSE, sep = "\t")
    }
  )
  
  output$linear_results_table_download <- downloadHandler(
    filename = function() {
      i <- which(type() == input$limma_dataset)
      paste(type()[i],"_LimmaResults_", Sys.Date(), ".txt", sep = "")
    },
    content = function(file){
      write.table(linear_model_results_table(), file, row.names = TRUE, sep = "\t", col.names = NA)
    }
  )
  
  output$pvaldist_plot_download <- downloadHandler(
    filename = function() {
      i <- which(type() == input$limma_dataset)
      paste(type()[i], "_", input$coef_options, "_PValueDistribution_", Sys.Date(), input$pvaldist_extension, sep = "")
    },
    content = function(file){
      ggplot2::ggsave(file, pval_distribution_plot(), width = 11, height = 6)
    }
  )
  
  ## BATCH EFFECT PCA ##########################################################
  
  ### REACTIVE OBJECTS ---------------------------------------------------------
  
  batch_corrected_data <- reactive({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)))
    i <- which(type() == input$limma_dataset)
    data <- limmaLM(annot = annotation(),
                    eset = eset()[[i]],
                    samples = input$linear_factors,
                    time_series = input$include_timeseries,
                    time_col = input$time_col,
                    time_points_cont = input$time_points_cont,
                    time_points_disc = input$time_points_disc,
                    time = input$time_type,
                    contrast_fit = input$contrast_fit,
                    contrasts_subset = input$limma_contrasts,
                    covariate = input$add_covariate,
                    covariate_col = input$covariate_col,
                    remove_batch_PCA = TRUE,
                    batch_column = input$batch_column)

    message("---- BATCH CORRECTED DATA GENERATED")
    return(data)
  })
  
  batch_PCA_plot <- reactive({
    i <- which(type() == input$limma_dataset)
    drawPCA(eset = batch_corrected_data(),
            x_axis = "PC1",
            y_axis = "PC2",
            type = type()[[i]],
            color = input$batch_PCA_color,
            include_densities = input$batch_PCA_include_densities,
            shapes = input$batch_PCA_shapes,
            title_add = paste0("After Batch Correction on Column/Factor '", input$batch_column, "'"),
            add_labels = input$batch_PCA_labels,
            add_ellipse = input$batch_PCA_ellipse)
  })
  
  pre_batch_PCA_plot <- reactive({
    i <- which(type() == input$limma_dataset)
    drawPCA(eset = eset()[[i]],
            x_axis = "PC1",
            y_axis = "PC2",
            type = type()[[i]],
            color = input$batch_PCA_color,
            include_densities = input$batch_PCA_include_densities,
            shapes = input$batch_PCA_shapes,
            title_add = paste0("Before Batch Correction"),
            add_labels = input$batch_PCA_labels,
            add_ellipse = input$batch_PCA_ellipse)
  })
  
  ### RENDER UI ----------------------------------------------------------------
  
  output$batch_PCA_column <- renderUI({
    req((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data))
    selectInput("batch_column",
                label = "Select column containing batch information",
                choices = annot_columns())
  })
  
  ### PRINT/TABLE/PLOTS ---------------------------------------------------------
  
  output$PCA_batch <- renderPlot({
    shiny::validate(need(input$batch_column != "", "Invalid batch column selected. Please select a different column from the drop down or add a batch column in your annotation file."))
    i <- which(type() == input$limma_dataset)
    gridExtra::grid.arrange(pre_batch_PCA_plot() + theme(legend.position = "bottom"), 
                            batch_PCA_plot() + theme(legend.position = "bottom"),
                            top = paste0("Pre/Post Batch Correction PCA Plots ", type()[i]), 
                            ncol = 2);
  }) %>% bindEvent(input$batch_PCA_button)
  
  ### DOWNLOADS -----------------------------------------------------------------
  
  output$batch_PCA_download <- downloadHandler(
    filename = function() {
      i <- which(type() == input$limma_dataset)
      paste(type()[i], "_BatchCorrectedPCA_", Sys.Date(), input$batch_PCA_extension , sep = "")
    },
    content = function(file){
      ggplot2::ggsave(file, batch_PCA_plot(), width = 9, height = 9)
    }
  )
  
  ## VOLCANO PLOT ##############################################################
  
  ### REACTIVE OBJECTS ---------------------------------------------------------
  volcano_linear_model_table <- reactive({
    req(input$volcano_coef_options)
    i <- which(type() == input$limma_dataset)
    
    makeTopTable(linear_model()[[i]],
                 coef = input$volcano_coef_options,
                 coef_options = input$limma_contrasts)
    
  })
  
  volcano_plot <- reactive({
    req((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data))
    i <- which(type() == input$limma_dataset)
    
    drawVolcano(dat = volcano_linear_model_table(), 
                type = type()[i],
                add_labels = input$volcano_labels,
                v = input$volcano_yaxis,
                up_color = input$volcano_up_color,
                down_color = input$volcano_down_color,
                title_add = input$volcano_title,
                title_type = input$volcano_coef_options,
                top_values = input$volcano_top_values,
                top_fc = input$volcano_top_fc,
                label_specific = input$volcano_label_specific,
                label_specific_gene = input$volcano_label_specific_gene,
                label_number = input$volcano_label_number)
  })
  
  ### RENDER UI ----------------------------------------------------------------
  output$volcano_coef_options <- renderUI({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$limma_button))
    
    selectInput("volcano_coef_options",
                label = "Contrast to plot",
                choices = coef_options(),
                multiple = FALSE)
  }) %>% bindEvent(input$limma_button)
  
  output$volcano_label_specific_gene <- renderUI({
    req(input$volcano_label_specific == TRUE)
    
    genes <- sub("_.*", "", volcano_linear_model_table()$feature_identifier)
    
    selectInput("volcano_label_specific_gene",
                label = "Type gene name",
                multiple = TRUE,
                choices = genes)
  }) %>% bindEvent(input$volcano_label_specific)
  
  output$volcano_label_number <- renderUI({
    req(input$volcano_labels == TRUE)
    
    numericInput("volcano_label_number",
                 label = "Number of genes to label",
                 value = 20)
  }) %>% bindEvent(input$volcano_labels)
  
  ### PRINT/TABLE/PLOTS --------------------------------------------------------
  output$volcano_plot <- renderPlot({
    req(input$volcano_coef_options)
    volcano_plot()
  }) %>% bindEvent(c(input$volcano_button, input$limma_button))
  
  ### DOWNLOADS ----------------------------------------------------------------
  output$volcano_download <- downloadHandler(
    filename = function() {
      i <- which(type() == input$limma_dataset)
      paste(type()[i], "_", input$volcano_coef_options, "_VolcanoPlot_", Sys.Date(), input$volcano_extension, sep = "")
    },
    content = function(file){
      ggsave(file, volcano_plot(), width = 10, height = 8)
    }
  )
  
  ## INTERACTIVE VOLCANO PLOT ##################################################
  
  ### REACTIVE OBJECTS ---------------------------------------------------------
  int_volcano_linear_model_table <- reactive({
    req(input$int_volcano_coef_options)
    i <- which(type() == input$limma_dataset)
    
    linear_model <- limmaLM(annot = annotation(),
                            eset = eset()[[i]],
                            samples = input$linear_factors,
                            time_series = input$include_timeseries,
                            time_col = input$time_col,
                            time_points_cont = input$time_points_cont,
                            time_points_disc = input$time_points_disc,
                            time = input$time_type,
                            contrast_fit = input$contrast_fit,
                            contrasts_subset = input$limma_contrasts,
                            covariate = input$add_covariate,
                            covariate_col = input$covariate_col)
    
    table <- makeTopTable(linear_model,
                 coef = input$int_volcano_coef_options,
                 coef_options = input$limma_contrasts)
    
    return(table)
  })
  
  int_volcano_plot <- reactive({
    req((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data))
    i <- which(type() == input$limma_dataset)
    
    plot <- drawVolcano(dat = int_volcano_linear_model_table(), 
                type = type()[i],
                add_labels = FALSE,
                v = input$int_volcano_yaxis,
                up_color = input$int_volcano_up_color,
                down_color = input$int_volcano_down_color,
                title_add = NULL,
                title_type = input$int_volcano_coef_options,
                top_values = input$int_volcano_top_values,
                top_fc = input$int_volcano_top_fc,
                label_specific = FALSE,
                label_specific_gene = FALSE,
                label_number = FALSE)
    
    message("PLOT GENERATED")
    
    return(plot)
  }) %>% bindEvent(input$int_volcano_button)
  
  glimma_args <- reactive({
    i <- which(type() == input$limma_dataset)
    
    args <- interactiveVolcano(eset = eset()[[i]],
                               fit = linear_model()[[i]], 
                               top_table = int_volcano_linear_model_table(), 
                               type = type()[i], 
                               coef = input$int_volcano_coef_options)
    
  })
  
  ### RENDER UI ----------------------------------------------------------------
  output$int_volcano_coef_options <- renderUI({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$limma_button))
    
    selectInput("int_volcano_coef_options",
                label = "Choose the coefficient/contrast to use for the plot",
                choices = coef_options(),
                multiple = FALSE)
  }) %>% bindEvent(input$limma_button)
  
  output$glimma <- renderUI({
    i <- which(type() == input$limma_dataset)
    includeHTML(paste("glimma-plots/", type()[i], "_GlimmaVolcanoPlot_", Sys.Date(), ".html", sep = ""))
  })
  
  ### PRINT/TABLE/PLOTS --------------------------------------------------------
  output$volcano_plot_interactive <- plotly::renderPlotly({
    plot <- plotly::ggplotly(int_volcano_plot())
  })
  
  observeEvent(input$glimma_volcano_button, {
    do.call(Glimma::glXYPlot, c(glimma_args()))
  })
  
  ## DIFFERENTIAL HEATMAP ######################################################
  
  ### REACTIVE OBJECTS ---------------------------------------------------------
  dif_heatmap_linear_model_table <- reactive({
    req(input$dif_heatmap_coef_options)
    i <- which(type() == input$limma_dataset)
    
    coef_options <- limmaLM(annot = annotation(),
                            eset = eset()[[i]],
                            samples = input$linear_factors,
                            coef_names = TRUE,
                            time_series = input$include_timeseries,
                            time_col = input$time_col,
                            time_points_cont = input$time_points_cont,
                            time_points_disc = input$time_points_disc,
                            time = input$time_type,
                            contrast_fit = input$contrast_fit,
                            contrasts_subset = input$limma_contrasts,
                            covariate = input$add_covariate,
                            covariate_col = input$covariate_col)
    
    makeTopTable(linear_model()[[i]],
                 coef = input$dif_heatmap_coef_options,
                 coef_options)
    
  })
  
  dif_heatmap_plot <- reactive({
    req((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data))
    i <- which(type() == input$limma_dataset)
    
    drawHeatmaps(eset()[[i]], 
                 type = type()[i], 
                 mapcolor = input$dif_heatmap_color, 
                 group_color = input$dif_heatmap_group_color,
                 subset = "top_table_rows",
                 kclustering = input$dif_kclustering,
                 log2 = input$log_transform,
                 cluster_samples = input$dif_cluster_samples, 
                 title_add = input$dif_heatmap_title,
                 annot = annotation(),
                 zscore = input$dif_zscore,
                 differential_heatmap = TRUE,
                 top_table = dif_heatmap_linear_model_table(),
                 dif_title = input$dif_heatmap_coef_options,
                 dif_number_rows = input$dif_number_rows)
  }) %>% bindEvent(input$dif_heatmap_button)
  
  ### RENDER UI ----------------------------------------------------------------
  output$dif_heatmap_coef_options <- renderUI({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$limma_button))
    
    selectInput("dif_heatmap_coef_options",
                label = "Choose the coefficient/contrast to use for the plot",
                choices = coef_options(),
                multiple = FALSE)
  }) %>% bindEvent(input$limma_button)
  
  ### PRINT/TABLE/PLOTS --------------------------------------------------------
  output$dif_heatmap_plot <- renderPlot({
    dif_heatmap_plot()
  }, height = 800)
  
  ### DOWNLOADS ----------------------------------------------------------------
  output$dif_heatmap_download <- downloadHandler(
    filename = function() {
      i <- which(type() == input$limma_dataset)
      paste(type()[i], "_", input$dif_heatmap_coef_options, "_DifferentialHeatmap_", Sys.Date(), input$dif_heatmap_extension, sep = "")
    },
    
    content = function(file){
      withProgress(message = "Downloading heatmap plot:", detail = "This may take a moment...",value = 0.3, {
        if(input$dif_heatmap_extension == ".pdf"){
          grDevices::pdf(file, width = 10, height = 10)
          ComplexHeatmap::draw(dif_heatmap_plot())
          grDevices::dev.off()
          setProgress(value = 0.9)
          
        } else if(input$dif_heatmap_extension == ".png"){
          grDevices::png(file, width = 1000, height = 1000)
          ComplexHeatmap::draw(dif_heatmap_plot())
          grDevices::dev.off()
          setProgress(value = 0.9)
          
        } else if(input$dif_heatmap_extension == ".svg"){
          svglite::svglite(file, width = 10, height = 10)
          ComplexHeatmap::draw(dif_heatmap_plot())
          grDevices::dev.off()
          setProgress(value = 0.9)
          
        }
      })
    }
  )
  
  ## ENRICHMENT ################################################################
  
  ### REACTIVE OBJECTS ---------------------------------------------------------
  
  geneList <- reactive({
    req(input$gmt)
    i <- which(type() == input$limma_dataset)
    geneList <- calculateGeneList(fit = linear_model()[[i]],
                                  coef = input$enrichment_coef,
                                  gmt = input$gmt)
    # assign("geneList", geneList, envir = .GlobalEnv)
    return(geneList)
  })
  
  enriched <- reactive({
    req(input$gmt)
    enriched <- clusterProfilerEnrichment(geneList = geneList(),
                                          gmt = input$gmt,
                                          enrichment = input$enrichment_calculation,
                                          pval_cutoff = input$enrichment_pval)
    # assign("enrichmentResult", enriched, envir = .GlobalEnv)
    return(enriched)
  })

  enriched_table <- reactive({
    req(input$gmt)
    shiny::validate(need((dim(enriched()@result)[1] != 0), "There are no enriched pathways found with the selected conditions."))
    enrichedTable(enriched = enriched(),
                  enrichment = input$enrichment_calculation)
  })
  
  enriched_dotplot <- reactive({
    enrichedPlots(enriched = enriched(),
                  gmt = input$gmt,
                  geneList = geneList(),
                  contrast = input$enrichment_coef,
                  plottype = "dot",
                  dotplot_title = input$dotplot_title,
                  dotplot_categories = input$dotplot_categories)
  })
  
  enriched_cnetplot <- reactive({
    enrichedPlots(enriched = enriched(),
                  gmt = input$gmt,
                  geneList = geneList(),
                  contrast = input$enrichment_coef,
                  plottype = "cnet",
                  cnetplot_categories = input$cnet_categories,
                  cnetplot_layout = input$cnet_layout,
                  cnetplot_label = input$cnet_labels,
                  cnet_title = input$cnet_title)
  })
  
  enriched_emapplot <- reactive({
    enrichedPlots(enriched = enriched(),
                  gmt = input$gmt,
                  geneList = geneList(),
                  contrast = input$enrichment_coef,
                  plottype = "emap",
                  emap_categories = input$emap_categories,
                  emap_title = input$emap_title)
  })
  
  enriched_upsetplot <- reactive({
    enrichedPlots(enriched = enriched(),
                  gmt = input$gmt,
                  geneList = geneList(),
                  contrast = input$enrichment_coef,
                  plottype = "upset",
                  upset_categories = input$upset_categories,
                  upsetplot_title = input$upset_title)
  })
  
  enriched_treeplot <- reactive({
    enrichedPlots(enriched = enriched(),
                  gmt = input$gmt,
                  geneList = geneList(),
                  contrast = input$enrichment_coef,
                  plottype = "tree",
                  tree_categories = input$tree_categories,
                  tree_clusters = input$tree_clusters,
                  treeplot_title = input$tree_title)
  })
  
  enriched_heatplot <- reactive({
    enrichedPlots(enriched = enriched(),
                  enrichment = input$enrichment_calculation,
                  gmt = input$gmt,
                  geneList = geneList(),
                  contrast = input$enrichment_coef,
                  plottype = "heat",
                  heatmap_color = input$enriched_heatmap_color,
                  heatmap_title = input$enriched_heatmap_title)
  })
  
  enriched_GOplot <- reactive({
    enrichedPlots(enriched = enriched(),
                  enrichment = input$enrichment_calculation,
                  gmt = input$gmt,
                  geneList = geneList(),
                  contrast = input$enrichment_coef,
                  plottype = "GO")
  }) 
  
  enriched_KEGGplot <- reactive({
    enrichedPlots(enriched = enriched(),
                  enrichment = input$enrichment_calculation,
                  gmt = input$gmt,
                  geneList = geneList(),
                  contrast = input$enrichment_coef,
                  plottype = "KEGG",
                  KEGGID = KEGG_pathway_ID())
  })
  
  KEGG_pathway_ID <- reactive({
    getKEGGID(gmt = input$gmt,
              KEGG_pathway = input$KEGG_pathway)
  })

  ### RENDER UI ----------------------------------------------------------------
  
  output$enrichment_coef_options <- renderUI({
    shiny::validate(need(input$limma_button, "Please generate a model in the linear modeling tab."))
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$limma_button))
    
    selectInput("enrichment_coef",
                label = "Choose the coefficient/contrast to use for enrichment",
                choices = coef_options(),
                multiple = FALSE)
  })
  
  output$gmts <- renderUI({
    req(input$species)
    i <- which(type() == input$limma_dataset)
    format <- data_format()[i]
    
    if (format == "PhosphoSites"){
      if (input$species == "human"){
        selectInput("gmt",
                    label = "Choose database to use for enrichment",
                    choices = c("KEGG" = "HUMAN_KEGG",
                                "Reactome" = "HUMAN_Reactome",
                                "WikiPathways" = "HUMAN_WikiPathways",
                                "GO Biological Processes" = "HUMAN_GOBP",
                                "GO Cellular Components" = "HUMAN_GOCC",
                                "GO Molecular Functions" = "HUMAN_GOMF",
                                "GO All" = "HUMAN_GOALL",
                                "Hallmark MSigdb" = "HUMAN_Hallmark",
                                "MSigdb" = "HUMAN_MSigDB",
                                "MOMENTA BioCyc" = "HUMAN_MOMENTABiocycEXP",
                                "MOMENTA MFN" = "HUMAN_MOMENTAMFNEXP",
                                "PTM-SEA All" = "HUMAN_PHOSPHO_UNIPROT",
                                "PTM-SEA Diseases" = "HUMAN_PHOSPHO_DISEASE",
                                "PTM-SEA Kinases" = "HUMAN_PHOSPHO_KINASE",
                                "PTM-SEA Pathways" = "HUMAN_PHOSPHO_PATHWAY",
                                "PTM-SEA Perturbation Sites" = "HUMAN_PHOSPHO_PERT"),
                    multiple = FALSE)
        
      } else if (input$species == "mouse"){
        selectInput("gmt",
                    label = "Choose database to use for enrichment",
                    choices = c("Reactome" = "MOUSE_Reactome",
                                "GO All" = "MOUSE_GOALL",
                                "GO Biological Processes" = "MOUSE_GOBP",
                                "GO Cellular Components" = "MOUSE_GOCC",
                                "GO Molecular Functions" = "MOUSE_GOMF",
                                "WikiPathways" = "MOUSE_WikiPathways",
                                "MOMENTA BioCyc" = "MOUSE_MOMENTABiocycEXP",
                                "PTM-SEA All" = "MOUSE_PHOSPHO_UNIPROT",
                                "PTM-SEA Diseases" = "MOUSE_PHOSPHO_DISEASE",
                                "PTM-SEA Kinases" = "MOUSE_PHOSPHO_KINASE",
                                "PTM-SEA Pathways" = "MOUSE_PHOSPHO_PATHWAY",
                                "PTM-SEA Perturbation Sites" = "MOUSE_PHOSPHO_PERT"),
                    multiple = FALSE)
        
      } else if (input$species == "rat"){
        selectInput("gmt",
                    label = "Choose database to use for enrichment",
                    choices = c("Reactome" = "RAT_Reactome",
                                "GO All" = "RAT_GOALL",
                                "GO Biological Processes" = "RAT_GOBP",
                                "GO Cellular Components" = "RAT_GOCC",
                                "GO Molecular Functions" = "RAT_GOMF",
                                "WikiPathways" = "RAT_WikiPathways",
                                "PTM-SEA All" = "RAT_PHOSPHO_UNIPROT",
                                "PTM-SEA Diseases" = "RAT_PHOSPHO_DISEASE",
                                "PTM-SEA Kinases" = "RAT_PHOSPHO_KINASE",
                                "PTM-SEA Pathways" = "RAT_PHOSPHO_PATHWAY",
                                "PTM-SEA Perturbation Sites" = "RAT_PHOSPHO_PERT"),
                    multiple = FALSE)
        
      } else if (input$species == "celegans"){
        selectInput("gmt",
                    label = "Choose database to use for enrichment",
                    choices = c("KEGG" = "CELEGANS_KEGG",
                                "WikiPathways" = "CELEGANS_WikiPathways",
                                "GO Biological Processes" = "CELEGANS_GOBP",
                                "GO Cellular Components" = "CELEGANS_GOCC",
                                "GO Molecular Functions" = "CELEGANS_GOMF"),
                    multiple = FALSE)
        
      } else if (input$species == "fruitfly"){
        selectInput("gmt",
                    label = "Choose database to use for enrichment",
                    choices = c("KEGG" = "FRUITFLY_KEGG",
                                "WikiPathways" = "FRUITFLY_WikiPathways",
                                "GO Biological Processes" = "FRUITFLY_GOBP",
                                "GO Cellular Components" = "FRUITFLY_GOCC",
                                "GO Molecular Functions" = "FRUITFLY_GOMF",
                                "MOMENTA BioCyc" = "FRUITFLY_MOMENTABiocycEXP"),
                    multiple = FALSE)
        
      } else if (input$species == "zebrafish"){
        selectInput("gmt",
                    label = "Choose database to use for enrichment",
                    choices = c("KEGG" = "ZEBRAFISH_KEGG",
                                "WikiPathways" = "ZEBRAFISH_WikiPathways",
                                "GO Biological Processes" = "ZEBRAFISH_GOBP",
                                "GO Cellular Components" = "ZEBRAFISH_GOCC",
                                "GO Molecular Functions" = "ZEBRAFISH_GOMF"),
                    multiple = FALSE)
        
      } else if (input$species == "yeast"){
        selectInput("gmt",
                    label = "Choose database to use for enrichment",
                    choices = c("KEGG" = "YEAST_KEGG",
                                "WikiPathways" = "YEAST_WikiPathways",
                                "GO Biological Processes" = "YEAST_GOBP",
                                "GO Cellular Components" = "YEAST_GOCC",
                                "GO Molecular Functions" = "YEAST_GOMF",
                                "MOMENTA BioCyc" = "YEAST_MOMENTABiocycEXP"),
                    multiple = FALSE)
        
      }
    } else if (format == "ProteinGroups"){
      if (input$species == "human"){
        selectInput("gmt",
                    label = "Choose database to use for enrichment",
                    choices = c("KEGG" = "HUMAN_KEGG",
                                "Reactome" = "HUMAN_Reactome",
                                "WikiPathways" = "HUMAN_WikiPathways",
                                "GO Biological Processes" = "HUMAN_GOBP",
                                "GO Cellular Components" = "HUMAN_GOCC",
                                "GO Molecular Functions" = "HUMAN_GOMF",
                                "GO All" = "HUMAN_GOALL",
                                "Hallmark MSigdb" = "HUMAN_Hallmark",
                                "MSigdb" = "HUMAN_MSigDB",
                                "MOMENTA BioCyc" = "HUMAN_MOMENTABiocycEXP",
                                "MOMENTA MFN" = "HUMAN_MOMENTAMFNEXP"),
                    multiple = FALSE)
        
      } else if (input$species == "mouse"){
        selectInput("gmt",
                    label = "Choose database to use for enrichment",
                    choices = c("Reactome" = "MOUSE_Reactome",
                                "GO All" = "MOUSE_GOALL",
                                "GO Biological Processes" = "MOUSE_GOBP",
                                "GO Cellular Components" = "MOUSE_GOCC",
                                "GO Molecular Functions" = "MOUSE_GOMF",
                                "WikiPathways" = "MOUSE_WikiPathways",
                                "MOMENTA BioCyc" = "MOUSE_MOMENTABiocycEXP"),
                    multiple = FALSE)
        
      } else if (input$species == "rat"){
        selectInput("gmt",
                    label = "Choose database to use for enrichment",
                    choices = c("Reactome" = "RAT_Reactome",
                                "GO All" = "RAT_GOALL",
                                "GO Biological Processes" = "RAT_GOBP",
                                "GO Cellular Components" = "RAT_GOCC",
                                "GO Molecular Functions" = "RAT_GOMF",
                                "WikiPathways" = "RAT_WikiPathways"),
                    multiple = FALSE)
        
      } else if (input$species == "celegans"){
        selectInput("gmt",
                    label = "Choose database to use for enrichment",
                    choices = c("KEGG" = "CELEGANS_KEGG",
                                "WikiPathways" = "CELEGANS_WikiPathways",
                                "GO Biological Processes" = "CELEGANS_GOBP",
                                "GO Cellular Components" = "CELEGANS_GOCC",
                                "GO Molecular Functions" = "CELEGANS_GOMF"),
                    multiple = FALSE)
        
      } else if (input$species == "fruitfly"){
        selectInput("gmt",
                    label = "Choose database to use for enrichment",
                    choices = c("KEGG" = "FRUITFLY_KEGG",
                                "WikiPathways" = "FRUITFLY_WikiPathways",
                                "GO Biological Processes" = "FRUITFLY_GOBP",
                                "GO Cellular Components" = "FRUITFLY_GOCC",
                                "GO Molecular Functions" = "FRUITFLY_GOMF",
                                "MOMENTA BioCyc" = "FRUITFLY_MOMENTABiocycEXP"),
                    multiple = FALSE)
        
      } else if (input$species == "zebrafish"){
        selectInput("gmt",
                    label = "Choose database to use for enrichment",
                    choices = c("KEGG" = "ZEBRAFISH_KEGG",
                                "WikiPathways" = "ZEBRAFISH_WikiPathways",
                                "GO Biological Processes" = "ZEBRAFISH_GOBP",
                                "GO Cellular Components" = "ZEBRAFISH_GOCC",
                                "GO Molecular Functions" = "ZEBRAFISH_GOMF"),
                    multiple = FALSE)
        
      } else if (input$species == "yeast"){
        selectInput("gmt",
                    label = "Choose database to use for enrichment",
                    choices = c("KEGG" = "YEAST_KEGG",
                                "WikiPathways" = "YEAST_WikiPathways",
                                "GO Biological Processes" = "YEAST_GOBP",
                                "GO Cellular Components" = "YEAST_GOCC",
                                "GO Molecular Functions" = "YEAST_GOMF",
                                "MOMENTA BioCyc" = "YEAST_MOMENTABiocycEXP"),
                    multiple = FALSE)
        
      }
    } else {
      p("Enrichment is only available for ProteinGroups & PhosphoSites data formats.")
    }

  })
  
  output$KEGG_pathway_options <- renderUI({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$group) && isTruthy(grepl("KEGG", input$gmt)) && isTruthy(input$enrichment_button))
    
    tagList(
        h5("KEGG Pathway Overview"),
        helpText("Enriched genes will be highlighted green in the pathway diagram."),
        
        br(),
        br(),
        
        selectInput("KEGG_pathway",
                    label = "Choose pathway to generate KEGG pathway overview",
                    choices = enriched_table()[,1]),
        
        actionButton("enriched_KEGG_button",
                     label = "Generate Plot",
                     style = "color: #fff; background-color: #9BD79A; 
                              border-color: #9BD79A; border-radius: 10px; border-width: 2px"),
        
        hr(),
        helpText("To download this diagram, right click on the image and select \"Save Image\".")
    )

  })
  
  output$notKEGG_message <- renderPrint({
    if (!grepl("KEGG", input$gmt)){
      return(cat("This visualization is not available as the KEGG database was not selected."))
    }
  })
  
  ### PRINT/TABLE/PLOTS --------------------------------------------------------
  
  output$enriched_table <- DT::renderDataTable({
    req(input$enrichment_button, (!is.null(enriched_table())))
    
    shiny::validate(need(input$gmt, "Enrichment analysis is not available for this data type."))
    if(input$enrichment_calculation == "ranked"){
      DT::datatable(enriched_table(), rownames = FALSE, options = list(order = list(1, 'asc'), pageLength = 10)) %>% 
        DT::formatRound(columns = c(6), digits = 3) %>%
        DT::formatSignif(columns = c(2,3,4), digits = 3) %>%
        DT::formatStyle('NES', fontWeight = 'bold', background = GeneTonic::styleColorBar_divergent(enriched_table()$NES, 'lightskyblue', 'lightcoral'))
    } else{
      DT::datatable(enriched_table(), rownames = FALSE, options = list(order = list(1, 'asc'), pageLength = 10)) %>% 
        DT::formatSignif(columns = c(2,3,4), digits = 3)%>%
        DT::formatStyle('Count', background = DT::styleColorBar(enriched_table()$Count, "plum"))
    }
  }) %>% bindEvent(input$enrichment_button)
  
  output$enriched_dotplot <- renderPlot({
    req(input$enrichment_button)
    enriched_dotplot()
  }) %>% bindEvent(input$dotplot_button, input$enrichment_button)
  
  output$enriched_cnetplot <- renderPlot({
    req(input$enrichment_button)
    enriched_cnetplot()
  }) %>% bindEvent(input$cnetplot_button, input$enrichment_button)
  
  output$enriched_emapplot <- renderPlot({
    req(input$enrichment_button)
    enriched_emapplot()
  }) %>% bindEvent(input$emap_button, input$enrichment_button)
  
  output$enriched_upsetplot <- renderPlot({
    req(input$enrichment_button)
    enriched_upsetplot()
  }) %>% bindEvent(input$upset_button, input$enrichment_button)
  
  output$enriched_treeplot <- renderPlot({
    req(input$enrichment_button)
    enriched_treeplot()
  }) %>% bindEvent(input$tree_button, input$enrichment_button)
  
  output$enriched_heatplot <- plotly::renderPlotly({
    req(input$enrichment_button)
    enriched_heatplot()
  }) %>% bindEvent(input$enriched_heatmap_button, input$enrichment_button)
  
  output$enriched_KEGGplot <- renderImage({
    req(input$enrichment_button, grepl("KEGG", input$gmt))
    
    outfile <- tempfile(fileext = ".png")
    
    enriched_KEGGplot()
    
    list(src = paste0(KEGG_pathway_ID(), ".pathview.png", sep = ""),
         contentType = 'image/png',
         width = "100%")
    
  }, deleteFile = TRUE) %>% bindEvent(input$enriched_KEGG_button)
  
  output$enriched_GOplot <- renderPlot({
    req(input$enrichment_button, grepl("GO", input$gmt))
    
    enriched_GOplot()
  })
  
  ### DOWNLOADS ----------------------------------------------------------------
  output$enriched_table_download <- downloadHandler(
    filename = function() {
      i <- which(type() == input$limma_dataset)
      paste(type()[i], "_", input$enrichment_coef, "_", input$gmt, "_EnrichmentTable_", Sys.Date(), ".txt", sep = "")
    },
    
    content = function(file){
      write.table(enriched_table(), file, row.names = FALSE, sep = "\t")
    }
  )
  
  output$dotplot_download <- downloadHandler(
    filename = function() {
      i <- which(type() == input$limma_dataset)
      paste(type()[i], "_", input$enrichment_coef, "_", input$gmt, "_EnrichmentDotPlot_", Sys.Date(), input$dotplot_extension, sep = "")
    },
    
    content = function(file){
      if(input$dotplot_extension == ".pdf"){
        grDevices::pdf(file, height = 12, width = 10)
        plot(enriched_dotplot())
        grDevices::dev.off()
        
      } else if(input$dotplot_extension == ".png"){
        grDevices::png(file, height = 1200, width = 900)
        plot(enriched_dotplot())
        grDevices::dev.off()
        
      } else if(input$dotplot_extension == ".svg"){
        svglite::svglite(file, height = 12, width = 10)
        plot(enriched_dotplot())
        grDevices::dev.off()
        
      }
    }
  )
  
  output$cnetplot_download <- downloadHandler(
    filename = function() {
      i <- which(type() == input$limma_dataset)
      paste(type()[i], "_", input$enrichment_coef, "_", input$gmt, "_EnrichmentCNET_", Sys.Date(), input$cnet_extension, sep = "")
    },
    
    content = function(file){
      if(input$cnet_extension == ".pdf"){
        grDevices::pdf(file, height = 12, width = 12)
        plot(enriched_cnetplot())
        grDevices::dev.off()
        
      } else if(input$cnet_extension == ".png"){
        grDevices::png(file, height = 1200, width = 1200)
        plot(enriched_cnetplot())
        grDevices::dev.off()
        
      } else if(input$cnet_extension == ".svg"){
        svglite::svglite(file, height = 12, width = 12)
        plot(enriched_cnetplot())
        grDevices::dev.off()
        
      }
    }
  )
  
  output$emap_download <- downloadHandler(
    filename = function() {
      i <- which(type() == input$limma_dataset)
      paste(type()[i], "_", input$enrichment_coef, "_", input$gmt, "_EnrichmentEMAP_", Sys.Date(), input$emap_extension, sep = "")
    },
    
    content = function(file){
      if(input$emap_extension == ".pdf"){
        grDevices::pdf(file, height = 12, width = 12)
        plot(enriched_emapplot())
        grDevices::dev.off()
        
      } else if(input$emap_extension == ".png"){
        grDevices::png(file, height = 1200, width = 1200)
        plot(enriched_emapplot())
        grDevices::dev.off()
        
      } else if(input$emap_extension == ".svg"){
        svglite::svglite(file, height = 12, width = 12)
        plot(enriched_emapplot())
        grDevices::dev.off()
        
      }
    }
  )
  
  output$upset_download <- downloadHandler(
    filename = function() {
      i <- which(type() == input$limma_dataset)
      paste(type()[i], "_", input$enrichment_coef, "_", input$gmt, "_EnrichmentUpSetPlot_", Sys.Date(), input$upset_extension, sep = "")
    },
    
    content = function(file){
      if(input$upset_extension == ".pdf"){
        grDevices::pdf(file, height = 12, width = 10)
        plot(enriched_upsetplot())
        grDevices::dev.off()
        
      } else if(input$upset_extension == ".png"){
        grDevices::png(file, height = 1200, width = 900)
        plot(enriched_upsetplot())
        grDevices::dev.off()
        
      } else if(input$upset_extension == ".svg"){
        svglite::svglite(file, height = 12, width = 10)
        plot(enriched_upsetplot())
        grDevices::dev.off()
        
      }
    }
  )
  
  output$tree_download <- downloadHandler(
    filename = function() {
      i <- which(type() == input$limma_dataset)
      paste(type()[i], "_", input$enrichment_coef, "_", input$gmt, "_EnrichmentTreePlot_", Sys.Date(), input$tree_extension, sep = "")
    },
    
    content = function(file){
      if(input$tree_extension == ".pdf"){
        grDevices::pdf(file, height = 12, width = 10)
        plot(enriched_treeplot())
        grDevices::dev.off()
        
      } else if(input$tree_extension == ".png"){
        grDevices::png(file, height = 1200, width = 900)
        plot(enriched_treeplot())
        grDevices::dev.off()
        
      } else if(input$tree_extension == ".svg"){
        svglite::svglite(file, height = 12, width = 10)
        plot(enriched_treeplot())
        grDevices::dev.off()
        
      }
    }
  )
  
  output$KEGG_download <- downloadHandler(
    filename = function() {
      i <- which(type() == input$limma_dataset)
      paste(type()[i], "_", input$enrichment_coef, "_", KEGG_pathway_ID, "_", Sys.Date(),".png", sep = "")
    },
    
    content = function(file){
      enriched_KEGGplot()
    }
  )
  
  ## S-SCORE INTEGRATION #######################################################
  
  ### REACTIVE OBJECTS ---------------------------------------------------------
  sscore_contrasts <- reactive({
    contrasts <- limmaLM(annot = annotation(),
                         eset = eset()[[1]],
                         samples = input$sscore_linear_factors,
                         time_series = input$sscore_include_timeseries,
                         time_col = input$sscore_time_col,
                         time_points_cont = input$sscore_time_points_cont,
                         time_points_disc = input$sscore_time_points_disc,
                         time = input$sscore_time_type,
                         contrast_fit = TRUE,
                         covariate = input$sscore_add_covariate,
                         covariate_col = input$sscore_covariate_col,
                         return_contrasts = TRUE)
  })
  
  sscore_toptable_list <- reactive({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$sscore_limma_contrasts) && isTruthy(input$sscore_button))
    message("*** STARTED SSCORE TOPTABLE GENERATION ***")
    
    linear_models <- list()
    top_tables <- list()
    
    datasets <- c()
    for (i in input$sscore_datasets){
      datasets <- c(datasets, which(type() == i))
    }
    
    for (i in 1:length(datasets)){
      linear_models[[i]] <- limmaLM(annot = annotation(),
                                    eset = eset()[[datasets[i]]],
                                    samples = input$sscore_linear_factors,
                                    time_series = input$sscore_include_timeseries,
                                    time_col = input$sscore_time_col,
                                    time_points_cont = input$sscore_time_points_cont,
                                    time_points_disc = input$sscore_time_points_disc,
                                    time = input$sscore_time_type,
                                    contrast_fit = TRUE,
                                    contrasts_subset = input$sscore_limma_contrasts,
                                    covariate = input$sscore_add_covariate,
                                    covariate_col = input$sscore_covariate_col)
    }
    
    coefs <- input$sscore_limma_contrasts
    for (j in 1:length(coefs)){
      top_tables_subset <- list()
      
      for (i in 1:length(linear_models)){
        fit <- linear_models[[i]]
        n <- nrow(fit$coefficients)
        
        end_col <- ncol(limma::topTable(fit, adjust = "BH", coef = coefs[j]))
        start_col <- which((colnames(limma::topTable(fit, adjust = "BH", coef = coefs[j]))) == "feature_identifier")
        top_tables_subset[[i]] <- limma::topTable(fit, adjust = "BH", number = n, coef = coefs[j])
      }
      names(top_tables_subset) <- input$sscore_datasets
      
      top_tables[[j]] <- top_tables_subset
    }
    names(top_tables) <- coefs
    
    message("--- COMPLETED TOPTABLE GENERATION")
    
    return(top_tables);
  }) %>% bindEvent(input$sscore_button)
  
  sscore_data_list <- reactive({
    message("*** STARTED SSCORE DATA LIST GENERATION ***")
    
    toptable_list <- sscore_toptable_list()
    
    data_formats <- c()
    for (i in input$sscore_datasets){
      data_formats <- c(data_formats, data_format()[which(type() == i)])
    }
    
    data_list <- list()
    for (i in 1:length(input$sscore_limma_contrasts)){
      data_list[[i]] <- makeDataList(top_table_list = toptable_list[[i]],
                                     data_format = data_formats)
    }
    names(data_list) <- input$sscore_limma_contrasts
    
    # assign("Sscore_DataList", data_list, envir = .GlobalEnv)
    message("--- COMPLETED SSCORE DATA LIST")
    
    return(data_list);
  })
  
  sscore_dataframe <- reactive({
    message("*** STARTED SSCORE INTEGRATION ***")
    
    withProgress(message = 'Performing Integration:', detail = "This may take a while...", value = 0, {
      sscore <- list()
      for (i in 1:length(sscore_data_list())){
        setProgress(value = (i / length(sscore_data_list()) - ((1 / length(sscore_data_list())) / 2)), 
                    detail = paste("Integrating contrast ", i, " of ", length(sscore_data_list())))
        sscore[[i]] <- sscoreIntegration(data_list = sscore_data_list()[[i]])
      }
      names(sscore) <- input$sscore_limma_contrasts
    })
    
    message("--- COMPLETED SSCORE INTEGRATION")
    
    return(sscore);
  })
  
  sscore_volcano <- reactive({
    i <- which(names(sscore_dataframe()) == input$sscore_coef_options)
    sscoreVolcanoPlot(sscore_dataframe()[[i]], 
                      adj_pval_cutoff = input$sscore_volcano_pval,
                      up_color = input$sscore_volcano_up_color,
                      down_color = input$sscore_volcano_down_color,
                      add_labels = input$sscore_volcano_label,
                      label_num = input$sscore_volcano_label_num) 
  })
  
  sscore_venn <- reactive({
    i <- which(names(sscore_dataframe()) == input$sscore_coef_options)
    sscoreVennDiagram(sscore_dataframe()[[i]])
  })
  
  ### RENDER UI ----------------------------------------------------------------
  output$sscore_datasets <- renderUI({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)))
    
    options <- type()
    selectInput("sscore_datasets",
                label = "Choose at least 2 datasets to include",
                choices = options,
                multiple = TRUE)
  })
  
  output$sscore_coef_options <- renderUI({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$sscore_limma_contrasts))
    
    selectInput("sscore_coef_options",
                label = "Choose contrast to compare across datasets via S-score integration",
                choices = input$sscore_limma_contrasts,
                multiple = FALSE)
  })
  
  output$sscore_linear_factors <- renderUI({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$group))
    
    options <- unique(pData(eset()[[1]])$Group)
    
    selectInput("sscore_linear_factors",
                label = "Choose at least two groups for comparison",
                choices = options,
                multiple = TRUE,
                selected = options)
  })
  
  output$sscore_choose_covariate <- renderUI({
    req(input$sscore_add_covariate == TRUE)
    
    selectizeInput("sscore_covariate_col",
                   label = "Select column to act as covariate",
                   choices = annot_columns(),
                   multiple = TRUE,
                   options = list(maxItems = 3))
    
  })
  
  output$sscore_choose_time <- renderUI({
    req(input$sscore_include_timeseries == TRUE)
    
    selectInput("sscore_time_col",
                label = "Select column containting time series information",
                choices = annot_columns(),
                multiple = FALSE)
  })
  
  output$sscore_time_type <- renderUI({
    req(input$sscore_include_timeseries == TRUE)
    
    radioButtons("sscore_time_type",
                 label = "Select how to include the time variable in the model",
                 choices = c("Continuous" = "continuous",
                             "Discrete" = "discrete"))
  }) %>% bindEvent(input$sscore_include_timeseries)
  
  output$sscore_time_points <- renderUI({
    req(isTruthy(input$sscore_include_timeseries))
    
    if(req(input$sscore_time_type) == "continuous"){
      number = length(unique(pData(eset()[[1]])[,input$time_col]))
      numericInput("sscore_time_points_cont",
                   label = "Choose number of time points to include in the model",
                   min = 0,
                   max = number,
                   value = 3) 
      
    } else if(req(input$sscore_time_type) == "discrete"){
      options <- unique(pData(eset()[[1]])[,input$sscore_time_col])
      options <- options[-order(options)[1]]
      options <- sort(options)
      selectInput("sscore_time_points_disc",
                  label = "Choose which time points to include in the model",
                  choices = options,
                  multiple = TRUE)
      
    }
    
  }) %>% bindEvent(input$sscore_time_type, input$sscore_include_timeseries, input$sscore_time_col)
  
  output$sscore_time_df_explain <- renderUI({
    req(input$sscore_include_timeseries == TRUE)
    
    tagList(helpText("When choosing the number of time points to include in the model, keep in mind that each added
             time point will decrease the overall residual degrees of freedom of the model. Once the residual
             DF reaches zero, the model will no longer run."),
            
            br(),
            br(),
            
            helpText("The first time point (usually 0) will be used as a base/reference automatically and so will not be an option for contrasts."),
            
            br(),
            br()
    )
  })
  
  output$sscore_choose_contrast <- renderUI({
    req(isTruthy(length(input$sscore_linear_factors) > 1))
    
    tagList(
      selectInput("sscore_limma_contrasts",
                  label = "Choose contrasts to include in model",
                  choices = sscore_contrasts(),
                  multiple = TRUE),
    )
  })
  
  ### PRINT/TABLES/PLOTS -------------------------------------------------------
  
  output$sscore_toptable_summary <- renderPrint({
    req(input$sscore_button)
    str(sscore_data_list())
  }) %>% bindEvent(input$sscore_button)
  
  output$sscore_integration_summary <- DT::renderDataTable({
    req(input$sscore_button)
    i <- which(names(sscore_dataframe()) == input$sscore_coef_options)
    
    end_column <- which(colnames(sscore_dataframe()[[i]]) == "sscore_adj_pval")
    table <- sscore_dataframe()[[i]][,1:end_column]
    DT::datatable(table, options = list(order = list((end_column - 1), 'asc'), pageLength = 15, scrollX = TRUE), rownames = FALSE) %>% 
      DT::formatSignif(columns = (end_column - 2):end_column, digits = 3) %>%
      DT::formatStyle('sscore', fontWeight = 'bold', background = GeneTonic::styleColorBar_divergent(table$sscore, 'lightskyblue', 'lightcoral')) %>%
      DT::formatStyle('sscore_adj_pval', backgroundColor = DT::styleInterval(0.05, c('darkseagreen', 'darksalmon')))
  }) %>% bindEvent(c(input$sscore_button, input$sscore_coef_options))
  
  output$sscore_volcano <- renderPlot({
    req(input$sscore_button)
    sscore_volcano() 
  }) %>% bindEvent(c(input$sscore_button, input$sscore_coef_options, input$sscore_volcano_button))
  
  output$sscore_volcano_interactive <- plotly::renderPlotly({
    req(input$sscore_button)
    plot <- plotly::ggplotly(sscore_volcano())
  }) %>% bindEvent(c(input$sscore_button, input$sscore_coef_options, input$sscore_interactive_volcano_button))
  
  output$sscore_venn <- renderPlot({
    req(input$sscore_button)
    sscore_venn()
  }) %>% bindEvent(c(input$sscore_button, input$sscore_coef_options, input$sscore_venn_button))
  
  ### DOWNLOADS ----------------------------------------------------------------
  
  output$sscore_table_download <- downloadHandler(
    filename = function() {
      paste(input$sscore_coef_options, "_SscoreTable_", Sys.Date(), ".txt", sep = "")
    },
    
    content = function(file){
      i <- which(names(sscore_dataframe()) == input$sscore_coef_options)
      write.table(sscore_dataframe()[[i]], file, row.names = FALSE, sep = "\t")
    }
  )
  
  output$sscore_volcano_download <- downloadHandler(
    filename = function() {
      paste(input$sscore_datasets, '_VolcanoPlot_', Sys.Date(), input$sscore_volcano_extension, sep = '')
    },
    content = function(file) {
      ggsave(file, plot = sscore_volcano(), width = 12, height = 7)
    }
  )
  
  ## S-SCORE ENRICHMENT ########################################################
  
  ### REACTIVE ----------------------------------------------------------------
  sscore_geneList <- reactive({
    i <- which(names(sscore_dataframe()) == input$sscore_enrichment_coef_options)
    formatSscoreGeneList(sscore_output = sscore_dataframe()[[i]])
  })
  
  sscore_enriched <- reactive({
    runSscoreEnrichment(geneList = sscore_geneList(),
                        gmt = input$sscore_gmt,
                        pval_cutoff = input$sscore_enrichment_pval)
  })
  
  sscore_enriched_table <- reactive({
    shiny::validate(need((dim(sscore_enriched()@result)[1] != 0), "There are no enriched pathways found with the selected conditions."))
    enrichedTable(enriched = sscore_enriched(),
                  enrichment = "ranked")
  })
  
  sscore_enriched_dotplot <- reactive({
    i <- which(type() %in% input$sscore_datasets)
    enrichedPlots(enriched = sscore_enriched(),
                  gmt = input$sscore_gmt,
                  geneList = sscore_geneList(),
                  contrast = input$sscore_enrichment_coef_options,
                  plottype = "dot",
                  dotplot_title = paste0("S-SCORE: ", paste0(type()[i], collapse = "; "), "\n", input$sscore_dotplot_title),
                  dotplot_categories = input$sscore_dotplot_categories)
  })
  
  sscore_enriched_cnetplot <- reactive({
    i <- which(type() %in% input$sscore_datasets)
    enrichedPlots(enriched = sscore_enriched(),
                  gmt = input$sscore_gmt,
                  geneList = sscore_geneList(),
                  contrast = input$sscore_enrichment_coef_options,
                  plottype = "cnet",
                  cnetplot_categories = input$sscore_cnet_categories,
                  cnetplot_layout = input$sscore_cnet_layout,
                  cnetplot_label = input$sscore_cnet_labels,
                  cnet_title = paste0("S-SCORE: ", paste0(type()[i], collapse = "; "), "\n", input$sscore_cnet_title))
  })
  
  sscore_enriched_emapplot <- reactive({
    i <- which(type() %in% input$sscore_datasets)
    enrichedPlots(enriched = sscore_enriched(),
                  gmt = input$sscore_gmt,
                  geneList = sscore_geneList(),
                  contrast = input$sscore_enrichment_coef_options,
                  plottype = "emap",
                  emap_categories = input$sscore_emap_categories,
                  emap_title = paste0("S-SCORE: ", paste0(type()[i], collapse = "; "), "\n", input$sscore_emap_title))
  })

  sscore_enriched_upsetplot <- reactive({
    i <- which(type() %in% input$sscore_datasets)
    enrichedPlots(enriched = sscore_enriched(),
                  gmt = input$sscore_gmt,
                  geneList = sscore_geneList(),
                  contrast = input$sscore_enrichment_coef_options,
                  plottype = "upset",
                  upset_categories = input$sscore_upset_categories,
                  upsetplot_title = paste0("S-SCORE: ", paste0(type()[i], collapse = "; "), "\n", input$sscore_upset_title))
  })

  sscore_enriched_treeplot <- reactive({
    i <- which(type() %in% input$sscore_datasets)
    enrichedPlots(enriched = sscore_enriched(),
                  gmt = input$sscore_gmt,
                  geneList = sscore_geneList(),
                  contrast = input$sscore_enrichment_coef_options,
                  plottype = "tree",
                  tree_categories = input$sscore_tree_categories,
                  tree_clusters = input$sscore_tree_clusters,
                  treeplot_title = paste0("S-SCORE: ", paste0(type()[i], collapse = "; "), "\n", input$sscore_tree_title))
  })

  sscore_enriched_heatplot <- reactive({
    i <- which(type() %in% input$sscore_datasets)
    enrichedPlots(enriched = sscore_enriched(),
                  enrichment = "ranked",
                  gmt = input$sscore_gmt,
                  geneList = sscore_geneList(),
                  contrast = input$sscore_enrichment_coef_options,
                  plottype = "heat",
                  heatmap_color = input$sscore_enriched_heatmap_color,
                  heatmap_title = input$sscore_enriched_heatmap_title)
  })

  # sscore_enriched_GOplot <- reactive({
  #   enrichedPlots(enriched = enriched(),
  #                 enrichment = input$enrichment_calculation,
  #                 gmt = input$gmt,
  #                 geneList = geneList(),
  #                 contrast = input$enrichment_coef,
  #                 plottype = "GO")
  # }) 

  sscore_enriched_KEGGplot <- reactive({
    enrichedPlots(enriched = sscore_enriched(),
                  enrichment = "ranked",
                  gmt = input$sscore_gmt,
                  geneList = sscore_geneList(),
                  contrast = input$sscore_enrichment_coef_options,
                  plottype = "KEGG",
                  KEGGID = sscore_KEGG_pathway_ID())
  })

  sscore_KEGG_pathway_ID <- reactive({
    getKEGGID(gmt = input$sscore_gmt,
              KEGG_pathway = input$sscore_KEGG_pathway)
  })
  
  
  ### RENDER UI ----------------------------------------------------------------
  output$sscore_enrichment_coef_options <- renderUI({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$sscore_limma_contrasts))
    
    selectInput("sscore_enrichment_coef_options",
                label = "Choose contrast to compare across datasets via S-score integration",
                choices = input$sscore_limma_contrasts,
                multiple = FALSE)
  })
  
  output$sscore_gmts <- renderUI({
    req(isTruthy(input$species) & isTruthy(input$sscore_datasets))
    
    i <- which(type() %in% input$sscore_datasets)
    metabolites <- grepl("Metabolite", paste0(data_format()[i], collapse = ", "))
    
    if (!metabolites){
      if (input$species == "human"){
        selectInput("sscore_gmt",
                    label = "Choose database to use for enrichment",
                    choices = c("KEGG" = "HUMAN_KEGG",
                                "Reactome" = "HUMAN_Reactome",
                                "WikiPathways" = "HUMAN_WikiPathways",
                                "GO Biological Processes" = "HUMAN_GOBP",
                                "GO Cellular Components" = "HUMAN_GOCC",
                                "GO Molecular Functions" = "HUMAN_GOMF",
                                "GO All" = "HUMAN_GOALL",
                                "Hallmark MSigdb" = "HUMAN_Hallmark",
                                "MSigdb" = "HUMAN_MSigDB",
                                "MOMENTA BioCyc" = "HUMAN_MOMENTABiocycEXP",
                                "MOMENTA MFN" = "HUMAN_MOMENTAMFNEXP",
                                "PTM-SEA All" = "HUMAN_PHOSPHO_UNIPROT",
                                "PTM-SEA Diseases" = "HUMAN_PHOSPHO_DISEASE",
                                "PTM-SEA Kinases" = "HUMAN_PHOSPHO_KINASE",
                                "PTM-SEA Pathways" = "HUMAN_PHOSPHO_PATHWAY",
                                "PTM-SEA Perturbation Sites" = "HUMAN_PHOSPHO_PERT"),
                    multiple = FALSE)
        
      } else if (input$species == "mouse"){
        selectInput("sscore_gmt",
                    label = "Choose database to use for enrichment",
                    choices = c("Reactome" = "MOUSE_Reactome",
                                "GO All" = "MOUSE_GOALL",
                                "GO Biological Processes" = "MOUSE_GOBP",
                                "GO Cellular Components" = "MOUSE_GOCC",
                                "GO Molecular Functions" = "MOUSE_GOMF",
                                "WikiPathways" = "MOUSE_WikiPathways",
                                "MOMENTA BioCyc" = "MOUSE_MOMENTABiocycEXP",
                                "PTM-SEA All" = "MOUSE_PHOSPHO_UNIPROT",
                                "PTM-SEA Diseases" = "MOUSE_PHOSPHO_DISEASE",
                                "PTM-SEA Kinases" = "MOUSE_PHOSPHO_KINASE",
                                "PTM-SEA Pathways" = "MOUSE_PHOSPHO_PATHWAY",
                                "PTM-SEA Perturbation Sites" = "MOUSE_PHOSPHO_PERT"),
                    multiple = FALSE)
        
      } else if (input$species == "rat"){
        selectInput("sscore_gmt",
                    label = "Choose database to use for enrichment",
                    choices = c("Reactome" = "RAT_Reactome",
                                "GO All" = "RAT_GOALL",
                                "GO Biological Processes" = "RAT_GOBP",
                                "GO Cellular Components" = "RAT_GOCC",
                                "GO Molecular Functions" = "RAT_GOMF",
                                "WikiPathways" = "RAT_WikiPathways",
                                "PTM-SEA All" = "RAT_PHOSPHO_UNIPROT",
                                "PTM-SEA Diseases" = "RAT_PHOSPHO_DISEASE",
                                "PTM-SEA Kinases" = "RAT_PHOSPHO_KINASE",
                                "PTM-SEA Pathways" = "RAT_PHOSPHO_PATHWAY",
                                "PTM-SEA Perturbation Sites" = "RAT_PHOSPHO_PERT"),
                    multiple = FALSE)
        
      } else if (input$species == "celegans"){
        selectInput("sscore_gmt",
                    label = "Choose database to use for enrichment",
                    choices = c("KEGG" = "CELEGANS_KEGG",
                                "WikiPathways" = "CELEGANS_WikiPathways",
                                "GO Biological Processes" = "CELEGANS_GOBP",
                                "GO Cellular Components" = "CELEGANS_GOCC",
                                "GO Molecular Functions" = "CELEGANS_GOMF"),
                    multiple = FALSE)
        
      } else if (input$species == "fruitfly"){
        selectInput("sscore_gmt",
                    label = "Choose database to use for enrichment",
                    choices = c("KEGG" = "FRUITFLY_KEGG",
                                "WikiPathways" = "FRUITFLY_WikiPathways",
                                "GO Biological Processes" = "FRUITFLY_GOBP",
                                "GO Cellular Components" = "FRUITFLY_GOCC",
                                "GO Molecular Functions" = "FRUITFLY_GOMF",
                                "MOMENTA BioCyc" = "FRUITFLY_MOMENTABiocycEXP"),
                    multiple = FALSE)
        
      } else if (input$species == "zebrafish"){
        selectInput("sscore_gmt",
                    label = "Choose database to use for enrichment",
                    choices = c("KEGG" = "ZEBRAFISH_KEGG",
                                "WikiPathways" = "ZEBRAFISH_WikiPathways",
                                "GO Biological Processes" = "ZEBRAFISH_GOBP",
                                "GO Cellular Components" = "ZEBRAFISH_GOCC",
                                "GO Molecular Functions" = "ZEBRAFISH_GOMF"),
                    multiple = FALSE)
        
      } else if (input$species == "yeast"){
        selectInput("sscore_gmt",
                    label = "Choose database to use for enrichment",
                    choices = c("KEGG" = "YEAST_KEGG",
                                "WikiPathways" = "YEAST_WikiPathways",
                                "GO Biological Processes" = "YEAST_GOBP",
                                "GO Cellular Components" = "YEAST_GOCC",
                                "GO Molecular Functions" = "YEAST_GOMF",
                                "MOMENTA BioCyc" = "YEAST_MOMENTABiocycEXP"),
                    multiple = FALSE)
        
      }
    } else {
      if (input$species == "human"){
        selectInput("sscore_gmt",
                    label = "Choose database to use for enrichment",
                    choices = c("Reactome Genes & Metabolites" = "HUMAN_METABOLGENE_Reactome"),
                    multiple = FALSE)
        
      } else {
        stop("There are not currently any pathway databases available for non-human species 
             that includes metabolites.")
      }
    }
  })
  
  output$sscore_KEGG_pathway_options <- renderUI({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$group) && isTruthy(grepl("KEGG", input$sscore_gmt)) && isTruthy(input$sscore_enrichment_button))
    
    tagList(
      h5("KEGG Pathway Overview"),
      helpText("Enriched genes will be highlighted green in the pathway diagram."),
      
      br(),
      br(),
      
      selectInput("sscore_KEGG_pathway",
                  label = "Choose pathway to generate KEGG pathway overview",
                  choices = sscore_enriched_table()[,1]),
      
      actionButton("sscore_enriched_KEGG_button",
                   label = "Generate Plot",
                   style = "color: #fff; background-color: #9BD79A; 
                              border-color: #9BD79A; border-radius: 10px; border-width: 2px"),
      
      hr(),
      helpText("To download this diagram, right click on the image and select \"Save Image\".")
    )
    
  })
  
  ### PRINT/TABLES/PLOTS -------------------------------------------------------
  
  output$sscore_enriched_table <- DT::renderDataTable({
    req(input$sscore_enrichment_button, (!is.null(sscore_enriched_table())))
    
    DT::datatable(sscore_enriched_table(), rownames = FALSE, options = list(order = list(1, 'asc'), pageLength = 10)) %>% 
      DT::formatRound(columns = c(6), digits = 3) %>%
      DT::formatSignif(columns = c(2,3,4), digits = 3) %>%
      DT::formatStyle('NES', fontWeight = 'bold', background = GeneTonic::styleColorBar_divergent(sscore_enriched_table()$NES, 'lightskyblue', 'lightcoral'))
  }) %>% bindEvent(input$sscore_enrichment_button)
  
  output$sscore_enriched_dotplot <- renderPlot({
    req(input$sscore_enrichment_button)
    sscore_enriched_dotplot()
  }) %>% bindEvent(input$sscore_dotplot_button, input$sscore_enrichment_button)
  
  output$sscore_enriched_cnetplot <- renderPlot({
    req(input$sscore_enrichment_button)
    sscore_enriched_cnetplot()
  }) %>% bindEvent(input$sscore_cnetplot_button, input$sscore_enrichment_button)

  output$sscore_enriched_emapplot <- renderPlot({
    req(input$sscore_enrichment_button)
    sscore_enriched_emapplot()
  }) %>% bindEvent(input$sscore_emap_button, input$sscore_enrichment_button)

  output$sscore_enriched_upsetplot <- renderPlot({
    req(input$sscore_enrichment_button)
    sscore_enriched_upsetplot()
  }) %>% bindEvent(input$sscore_upset_button, input$sscore_enrichment_button)

  output$sscore_enriched_treeplot <- renderPlot({
    req(input$sscore_enrichment_button)
    sscore_enriched_treeplot()
  }) %>% bindEvent(input$sscore_tree_button, input$sscore_enrichment_button)

  output$sscore_enriched_heatplot <- plotly::renderPlotly({
    req(input$sscore_enrichment_button)
    sscore_enriched_heatplot()
  }) %>% bindEvent(input$sscore_enriched_heatmap_button, input$sscore_enrichment_button)

  output$sscore_enriched_KEGGplot <- renderImage({
    req(input$sscore_enrichment_button, grepl("KEGG", input$sscore_gmt))

    outfile <- tempfile(fileext = ".png")

    sscore_enriched_KEGGplot()

    list(src = paste0(sscore_KEGG_pathway_ID(), ".pathview.png", sep = ""),
         contentType = 'image/png',
         width = "100%",
         alt = "This is alternate text")

  }, deleteFile = TRUE) %>% bindEvent(input$sscore_enriched_KEGG_button)
  
  output$sscore_notKEGG_message <- renderPrint({
    if (!grepl("KEGG", input$sscore_gmt)){
      return(cat("This visualization is not available as the KEGG database was not selected."))
    }
  })
  
  ### DOWNLOADS ----------------------------------------------------------------
  output$sscore_enriched_table_download <- downloadHandler(
    filename = function() {
      i <- which(type() %in% input$sscore_datasets)
      paste("S-Score_", paste0(type()[i], collapse = "-"), "_", input$sscore_enrichment_coef_options, "_", input$sscore_gmt, "_EnrichmentTable_", Sys.Date(), ".txt", sep = "")
    },
    
    content = function(file){
      write.table(sscore_enriched_table(), file, row.names = FALSE, sep = "\t")
    }
  )
  
  output$sscore_dotplot_download <- downloadHandler(
    filename = function() {
      i <- which(type() %in% input$sscore_datasets)
      paste("S-Score_", paste0(type()[i], collapse = "-"), "_", input$sscore_enrichment_coef_options, "_", input$sscore_gmt, "_EnrichmentDotPlot_", Sys.Date(), input$sscore_dotplot_extension, sep = "")
    },
    
    content = function(file){
      if(input$sscore_dotplot_extension == ".pdf"){
        grDevices::pdf(file, height = 12, width = 10)
        plot(sscore_enriched_dotplot())
        grDevices::dev.off()
        
      } else if(input$sscore_dotplot_extension == ".png"){
        grDevices::png(file, height = 1200, width = 900)
        plot(sscore_enriched_dotplot())
        grDevices::dev.off()
        
      } else if(input$sscore_dotplot_extension == ".svg"){
        svglite::svglite(file, height = 12, width = 10)
        plot(sscore_enriched_dotplot())
        grDevices::dev.off()
        
      }
    }
  )
  
  output$sscore_cnetplot_download <- downloadHandler(
    filename = function() {
      i <- which(type() %in% input$sscore_datasets)
      paste("S-Score_", paste0(type()[i], collapse = "-"), "_", input$sscore_enrichment_coef_options, "_", input$sscore_gmt, "_EnrichmentCNET_", Sys.Date(), input$sscore_cnet_extension, sep = "")
    },

    content = function(file){
      if(input$sscore_cnet_extension == ".pdf"){
        grDevices::pdf(file, height = 12, width = 12)
        plot(sscore_enriched_cnetplot())
        grDevices::dev.off()

      } else if(input$sscore_cnet_extension == ".png"){
        grDevices::png(file, height = 1200, width = 1200)
        plot(sscore_enriched_cnetplot())
        grDevices::dev.off()

      } else if(input$sscore_cnet_extension == ".svg"){
        svglite::svglite(file, height = 12, width = 12)
        plot(sscore_enriched_cnetplot())
        grDevices::dev.off()

      }
    }
  )

  output$sscore_emap_download <- downloadHandler(
    filename = function() {
      i <- which(type() %in% input$sscore_datasets)
      paste("S-Score_", paste0(type()[i], collapse = "-"), "_", input$sscore_enrichment_coef_options, "_", input$sscore_gmt, "_EnrichmentEMAP_", Sys.Date(), input$sscore_emap_extension, sep = "")
    },

    content = function(file){
      if(input$sscore_emap_extension == ".pdf"){
        grDevices::pdf(file, height = 12, width = 12)
        plot(sscore_enriched_emapplot())
        grDevices::dev.off()

      } else if(input$sscore_emap_extension == ".png"){
        grDevices::png(file, height = 1200, width = 1200)
        plot(sscore_enriched_emapplot())
        grDevices::dev.off()

      } else if(input$sscore_emap_extension == ".svg"){
        svglite::svglite(file, height = 12, width = 12)
        plot(sscore_enriched_emapplot())
        grDevices::dev.off()

      }
    }
  )

  output$sscore_upset_download <- downloadHandler(
    filename = function() {
      i <- which(type() %in% input$sscore_datasets)
      paste("S-Score_", paste0(type()[i], collapse = "-"), "_", input$sscore_enrichment_coef_options, "_", input$sscore_gmt, "_EnrichmentUpSetPlot_", Sys.Date(), input$sscore_upset_extension, sep = "")
    },

    content = function(file){
      if(input$sscore_upset_extension == ".pdf"){
        grDevices::pdf(file, height = 12, width = 10)
        plot(sscore_enriched_upsetplot())
        grDevices::dev.off()

      } else if(input$sscore_upset_extension == ".png"){
        grDevices::png(file, height = 1200, width = 900)
        plot(sscore_enriched_upsetplot())
        grDevices::dev.off()

      } else if(input$sscore_upset_extension == ".svg"){
        svglite::svglite(file, height = 12, width = 10)
        plot(sscore_enriched_upsetplot())
        grDevices::dev.off()

      }
    }
  )

  output$sscore_tree_download <- downloadHandler(
    filename = function() {
      i <- which(type() %in% input$sscore_datasets)
      paste("S-Score_", paste0(type()[i], collapse = "-"), "_", input$sscore_enrichment_coef_options, "_", input$sscore_gmt, "_EnrichmentTreePlot_", Sys.Date(), input$sscore_tree_extension, sep = "")
    },

    content = function(file){
      if(input$sscore_tree_extension == ".pdf"){
        grDevices::pdf(file, height = 12, width = 10)
        plot(sscore_enriched_treeplot())
        grDevices::dev.off()

      } else if(input$sscore_tree_extension == ".png"){
        grDevices::png(file, height = 1200, width = 900)
        plot(sscore_enriched_treeplot())
        grDevices::dev.off()

      } else if(input$sscore_tree_extension == ".svg"){
        svglite::svglite(file, height = 12, width = 10)
        plot(sscore_enriched_treeplot())
        grDevices::dev.off()

      }
    }
  )

  output$sscore_KEGG_download <- downloadHandler(
    filename = function() {
      i <- which(type() %in% input$sscore_datasets)
      paste("S-Score_", paste0(type()[i], collapse = "-"), "_", input$sscore_enrichment_coef_options, "_KEGG_", KEGG_pathway_ID, "_", Sys.Date(),".png", sep = "")
    },

    content = function(file){
      sscore_enriched_KEGGplot()
    }
  )
  
  ## PCSF ######################################################################
  
  ### REACTIVE OBJECTS ---------------------------------------------------------
  
  PCSF_input <- reactive({
    if (input$network_input_type == "Single-omic linear model output") {
      i <- which(type() == input$network_SOLM_dataset)
      
      message("** STARTED PCSF INPUT GENERATION")
      message("-- DATASET: ", type()[i])
      
      fit <- limmaLM(annot = annotation(),
                     eset = eset()[[i]],
                     samples = input$linear_factors,
                     time_series = input$include_timeseries,
                     time_col = input$time_col,
                     time_points_cont = input$time_points_cont,
                     time_points_disc = input$time_points_disc,
                     time = input$time_type,
                     contrast_fit = input$contrast_fit,
                     contrasts_subset = input$limma_contrasts,
                     covariate = input$add_covariate,
                     covariate_col = input$covariate_col)
      
      coef <- which(input$limma_contrasts == input$network_SOLM_coef)
      top_table <- limma::topTable(fit, adjust = "BH", coef = coef, number = nrow(fit$coefficients))
      
      end_col <- ncol(top_table)
      
      # TODO: Add ability to pull from metabolite files. 
        # Issue: Inconsistent naming conventions for different metabolite input file types.
        # Example file uses "identifier", MetaboAnalyst has "HMDBID" column, and MS-DIAL has 
        # "InChIKey". Right now in Sscore.R the makeDataList() function references "identifier".
      
      # if ("identifier" %in% colnames(top_table_list[[i]])){
      #   dataframe[,1] <- sub(',.*', '', top_table_list[[i]]$identifier) # assuming HMDB column called "identifier" as in openMS output
      #   dataframe[,1] <- sub('.*c.{1}', '', dataframe[,1])
      # } else if ("HMDBID" %in% colnames(top_table_list[[i]])){
      #   dataframe[,1] <- sub(';.*', '', top_table_list[[i]]$HMDBID) # assuming HMDB column called "HMDBID" as in MetaboAnalyst annotated_peaklist output
      # } else if ("InChIKey" %in% colnames(top_table_list[[i]])){
      #   dataframe[,1] <- top_table_list[[i]]$InChIKey # assuming InChIKey column called "InChIKey" as in MS-DIAL output
      # }
      # 
      # dataframe[,2] <- rownames(top_table_list[[i]])
      # dataframe[,3] <- top_table_list[[i]]$logFC
      # 
      # colnames(dataframe) <- c("ID", "feature_id", "logfc")
      # 
      # dataframe$ID <- gsub("^$|^ $", NA, dataframe$ID)
      # dataframe <- na.omit(dataframe)
      
      shiny::validate(need(grepl("^Gene$", colnames(top_table)), "There needs to be a gene names column ('Gene') in the dataset."))
      
      start_col <- which((colnames(top_table)) == "Gene")
      dataframe <- top_table[,seq(start_col, end_col)]
      
      message("---- NUM ROWS: ", nrow(dataframe))
      message("---- NUM ROWS BELOW PVAL CUTOFF: ", nrow(dataframe[dataframe$adj.P.Val <= input$network_SOLM_pval_cutoff, ]))
      message("---- NUM ROWS ABOVE LOGFC CUTOFF: ", nrow(dataframe[abs(dataframe$logFC) >= input$network_SOLM_logfc_cutoff, ]))
      message("---- NUM ROWS MEETING BOTH: ", nrow(dataframe[(abs(dataframe$logFC) >= input$network_SOLM_logfc_cutoff & 
                                                               dataframe$adj.P.Val <= input$network_SOLM_pval_cutoff), ]))
      
      dataframe <- dataframe[(abs(dataframe$logFC) >= input$network_SOLM_logfc_cutoff & 
                                dataframe$adj.P.Val <= input$network_SOLM_pval_cutoff), ]
      
      PCSF_input <- dataframe[, c("Gene", "logFC", "adj.P.Val")]
      colnames(PCSF_input) <- c("gene_symbol", "logfc", "adj_pval")
      
    } else if (input$network_input_type == "Multi-omic S-score output"){
      coef <- which(input$sscore_limma_contrasts == input$network_MOLM_coef)
      dataframe <- sscore_dataframe()[[coef]]
      dataframe <- dataframe[(abs(dataframe$sscore) >= input$network_MOLM_sscore_cutoff &
                                dataframe$sscore_adj_pval <= input$network_MOLM_pval_cutoff), ]
      
      PCSF_input <- dataframe[, c("sscore", "sscore_adj_pval")]
      dataframe[is.na(dataframe)] <- ""
      PCSF_input$gene_symbol <- paste0(dataframe$gene_symbol, dataframe$chebi_id)
      PCSF_input <- PCSF_input[, c("gene_symbol", "sscore", "sscore_adj_pval")]
      colnames(PCSF_input) <- c("gene_symbol", "sscore", "adj_pval")
      
    }
    
    message("COMPLETED PCSF INPUT GENERATION **")
    # assign("PCSF_input", PCSF_input, envir = .GlobalEnv)
    
    return(PCSF_input)
  })
  
  PCSF_network <- reactive({
    pcsf_input <- PCSF_input()
    pcsf_net <- makePCSFNetwork(pcsf_input_data = pcsf_input,
                                pcsf_nval = 10,
                                mu = input$pcsf_mu)
    
    # assign("pcsf_net", pcsf_net, envir = .GlobalEnv)
    return(pcsf_net)
  }) %>% bindEvent(c(input$network_input_type, input$network_SOLM_dataset, input$network_SOLM_coef,
                     input$network_SOLM_pval_cutoff, input$network_SOLM_logfc_cutoff, input$network_MOLM_coef,
                     input$network_MOLM_pval_cutoff, input$network_MOLM_logfc_cutoff, input$pcsf_mu))
  
  PCSF_enriched <- reactive({
    PCSF_enriched <- pcsfRunEnrichment(pcsf_net = PCSF_network(),
                                       gmt = input$PCSF_gmt)
    
    return(PCSF_enriched)
  }) 
  
  PCSF_enriched_table <- reactive({
    # assign("enrichment_output", PCSF_enriched(), envir = .GlobalEnv)
    if (!is.null(PCSF_enriched())){
      table <- pcsfEnrichedTable(PCSF_enriched())
    } else {
      table <- NULL
    }
  })
  
  ### RENDER UI ----------------------------------------------------------------
  
  output$network_input_options <- renderUI({
    shiny::validate(need((((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$group)), 
                  "Select valid inputs in the Data tab to perform analyses."))
    
    options <- c()
    if (isTruthy(input$limma_button)) {options[1] <- "Single-omic linear model output"}
    if (isTruthy(input$sscore_button)) {options <- c(options, "Multi-omic S-score output")}
    
    shiny::validate(need((length(options) > 0), "Please run at least one analysis in the Single Omic Analysis or Multi-Omic Integration tab to begin Network Analysis."))
    
    radioButtons("network_input_type",
                 label = "Choose the data to use as input for network analysis",
                 choices = options)
  })
  
  output$network_SOLM_dataset <- renderUI({
    req(input$network_input_type == "Single-omic linear model output")
    
    options <- type()
    
    tagList(
      hr(),
      h6("Dataset"),
      selectInput("network_SOLM_dataset",
                  label = "Choose dataset to visualize",
                  choices = options)
    )
  })
  
  output$network_SOLM_contrast <- renderUI({
    req(input$network_input_type == "Single-omic linear model output")
    
    selectInput("network_SOLM_coef",
                label = "Choose contrast to analyze",
                choices = coef_options())
  })
  
  output$network_SOLM_pval_cutoff <- renderUI({
    req(input$network_input_type == "Single-omic linear model output")
    
    numericInput("network_SOLM_pval_cutoff",
                 label = "Adj. P-Val Cutoff",
                 value = 0.05,
                 min = 0,
                 max = 1,
                 step = 0.05)
  })
  
  output$network_SOLM_logfc_cutoff <- renderUI({
    req(input$network_input_type == "Single-omic linear model output")
    
    numericInput("network_SOLM_logfc_cutoff",
                 label = "LogFC Cutoff",
                 value = 1,
                 min = 0,
                 step = 0.5)
  })
  
  output$network_MOLM_contrast <- renderUI({
    req(input$network_input_type == "Multi-omic S-score output")
    
    selectInput("network_MOLM_coef",
                label = "Choose contrast to analyze",
                choices = input$sscore_limma_contrasts)
  })
  
  output$network_MOLM_pval_cutoff <- renderUI({
    req(input$network_input_type == "Multi-omic S-score output")
    
    numericInput("network_MOLM_pval_cutoff",
                 label = "Adj. P-Val Cutoff",
                 value = 0.05,
                 min = 0,
                 max = 1,
                 step = 0.05)
  })
  
  output$network_MOLM_sscore_cutoff <- renderUI({
    req(input$network_input_type == "Multi-omic S-score output")
    
    numericInput("network_MOLM_sscore_cutoff",
                 label = "S-Score Cutoff",
                 value = 0,
                 min = 0,
                 step = 1)
  })
  
  output$PCSF_gmts <- renderUI({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$species))
    
    if (!grepl("CHEBI:", paste0(PCSF_input()$gene_symbol, collapse = "; "))){
      if (input$species == "human"){
        selectInput("PCSF_gmt",
                    label = "Choose database to use for enrichment",
                    choices = c("KEGG" = "HUMAN_KEGG",
                                "Reactome" = "HUMAN_Reactome",
                                "WikiPathways" = "HUMAN_WikiPathways",
                                "GO Biological Processes" = "HUMAN_GOBP",
                                "GO Cellular Components" = "HUMAN_GOCC",
                                "GO Molecular Functions" = "HUMAN_GOMF",
                                "GO All" = "HUMAN_GOALL",
                                "Hallmark MSigdb" = "HUMAN_Hallmark",
                                "MSigdb" = "HUMAN_MSigDB",
                                "MOMENTA MFN" = "HUMAN_MOMENTAMFNEXP",
                                "MOMENTA BioCyc" = "HUMAN_MOMENTABiocycEXP",
                                "PTM-SEA All" = "HUMAN_PHOSPHO_UNIPROT",
                                "PTM-SEA Diseases" = "HUMAN_PHOSPHO_DISEASE",
                                "PTM-SEA Kinases" = "HUMAN_PHOSPHO_KINASE",
                                "PTM-SEA Pathways" = "HUMAN_PHOSPHO_PATHWAY",
                                "PTM-SEA Perturbation Sites" = "HUMAN_PHOSPHO_PERT"),
                    multiple = FALSE)
        
      } else if (input$species == "mouse"){
        selectInput("PCSF_gmt",
                    label = "Choose database to use for enrichment",
                    choices = c("Reactome" = "MOUSE_Reactome",
                                "GO All" = "MOUSE_GOALL",
                                "GO Biological Processes" = "MOUSE_GOBP",
                                "GO Cellular Components" = "MOUSE_GOCC",
                                "GO Molecular Functions" = "MOUSE_GOMF",
                                "WikiPathways" = "MOUSE_WikiPathways",
                                "MOMENTA BioCyc" = "MOUSE_MOMENTABiocycEXP",
                                "PTM-SEA All" = "MOUSE_PHOSPHO_UNIPROT",
                                "PTM-SEA Diseases" = "MOUSE_PHOSPHO_DISEASE",
                                "PTM-SEA Kinases" = "MOUSE_PHOSPHO_KINASE",
                                "PTM-SEA Pathways" = "MOUSE_PHOSPHO_PATHWAY",
                                "PTM-SEA Perturbation Sites" = "MOUSE_PHOSPHO_PERT"),
                    multiple = FALSE)
        
      } else if (input$species == "rat"){
        selectInput("PCSF_gmt",
                    label = "Choose database to use for enrichment",
                    choices = c("Reactome" = "RAT_Reactome",
                                "GO All" = "RAT_GOALL",
                                "GO Biological Processes" = "RAT_GOBP",
                                "GO Cellular Components" = "RAT_GOCC",
                                "GO Molecular Functions" = "RAT_GOMF",
                                "WikiPathways" = "RAT_WikiPathways",
                                "PTM-SEA All" = "RAT_PHOSPHO_UNIPROT",
                                "PTM-SEA Diseases" = "RAT_PHOSPHO_DISEASE",
                                "PTM-SEA Kinases" = "RAT_PHOSPHO_KINASE",
                                "PTM-SEA Pathways" = "RAT_PHOSPHO_PATHWAY",
                                "PTM-SEA Perturbation Sites" = "RAT_PHOSPHO_PERT"),
                    multiple = FALSE)
        
      } else if (input$species == "celegans"){
        selectInput("PCSF_gmt",
                    label = "Choose database to use for enrichment",
                    choices = c("KEGG" = "CELEGANS_KEGG",
                                "WikiPathways" = "CELEGANS_WikiPathways",
                                "GO Biological Processes" = "CELEGANS_GOBP",
                                "GO Cellular Components" = "CELEGANS_GOCC",
                                "GO Molecular Functions" = "CELEGANS_GOMF"),
                    multiple = FALSE)
        
      } else if (input$species == "fruitfly"){
        selectInput("PCSF_gmt",
                    label = "Choose database to use for enrichment",
                    choices = c("KEGG" = "FRUITFLY_KEGG",
                                "WikiPathways" = "FRUITFLY_WikiPathways",
                                "GO Biological Processes" = "FRUITFLY_GOBP",
                                "GO Cellular Components" = "FRUITFLY_GOCC",
                                "GO Molecular Functions" = "FRUITFLY_GOMF",
                                "MOMENTA BioCyc" = "FRUITFLY_MOMENTABiocycEXP"),
                    multiple = FALSE)
        
      } else if (input$species == "zebrafish"){
        selectInput("PCSF_gmt",
                    label = "Choose database to use for enrichment",
                    choices = c("KEGG" = "ZEBRAFISH_KEGG",
                                "WikiPathways" = "ZEBRAFISH_WikiPathways",
                                "GO Biological Processes" = "ZEBRAFISH_GOBP",
                                "GO Cellular Components" = "ZEBRAFISH_GOCC",
                                "GO Molecular Functions" = "ZEBRAFISH_GOMF"),
                    multiple = FALSE)
        
      } else if (input$species == "yeast"){
        selectInput("PCSF_gmt",
                    label = "Choose database to use for enrichment",
                    choices = c("KEGG" = "YEAST_KEGG",
                                "WikiPathways" = "YEAST_WikiPathways",
                                "GO Biological Processes" = "YEAST_GOBP",
                                "GO Cellular Components" = "YEAST_GOCC",
                                "GO Molecular Functions" = "YEAST_GOMF",
                                "MOMENTA BioCyc" = "YEAST_MOMENTABiocycEXP"),
                    multiple = FALSE)
        
      }
    } else {
      if (input$species == "human"){
        selectInput("PCSF_gmt",
                    label = "Choose database to use for enrichment",
                    choices = c("Reactome Genes & Metabolites" = "HUMAN_METABOLGENE_Reactome"),
                    multiple = FALSE)
        
      } else {
        stop("There are not currently any pathway databases available for non-human species 
             that includes metabolites.")
      }
    }
    
    
  })
  
  ### PRINT/TABLES/PLOTS -------------------------------------------------------
  
  output$PCSF_Nodes_Network <- visNetwork::renderVisNetwork({
    req(input$network_input_type, input$PCSF_Nodes_Network_button)
    shiny::validate(need(nrow(PCSF_input()) > 0, "There are no significant variables based on the applied cutoffs."))
    
    PCSFVisNodes(pcsf_net = PCSF_network(),
                 pcsf_input_data = PCSF_input()) 
  }) %>% bindEvent(input$PCSF_Nodes_Network_button)
  
  output$PCSF_Influential_Network <- visNetwork::renderVisNetwork({
    req(input$network_input_type, input$PCSF_Influential_Network_button)
    shiny::validate(need(nrow(PCSF_input()) > 0, "There are no significant variables based on the applied cutoffs."))
    
    PCSFVisInfluential(pcsf_net = PCSF_network()) 
  }) %>% bindEvent(input$PCSF_Influential_Network_button)
  
  output$PCSF_input_table <- DT::renderDataTable({
    req(input$network_input_type)
    DT::datatable(PCSF_input(), rownames = FALSE) %>% 
      DT::formatSignif(columns = c(2,3), digits = 3) %>%
      DT::formatStyle('adj_pval', backgroundColor = DT::styleInterval(0.05, c('darkseagreen', 'darksalmon')))
  })
  
  output$PCSF_enriched_table <- DT::renderDataTable({
    req(input$network_input_type, input$PCSF_enrichment_button)
    
    shiny::validate(need(!is.null(PCSF_enriched()), "There are no enriched clusters to display. Try choosing another database for enrichment, adjusting the cutoff values, or selecting another input."))
    
    DT::datatable(PCSF_enriched_table(), rownames = FALSE, options = list(order = list(3, 'asc'))) %>% 
      DT::formatSignif(columns = 4, digits = 3) %>%
      DT::formatStyle('Adj_PVal', backgroundColor = DT::styleInterval(0.05, c('darkseagreen', 'darksalmon')))
  }) %>% bindEvent(input$PCSF_enrichment_button)
  
  output$PCSF_Enrichment_Network <- visNetwork::renderVisNetwork({
    req(input$network_input_type, input$PCSF_enrichment_button)
    
    shiny::validate(need(!is.null(PCSF_enriched()), "There are no enriched clusters to display. Try choosing another database for enrichment, adjusting the cutoff values, or selecting another input."))
    
    pcsfEnrichedSubnet(pcsf_enrich_pathway = PCSF_enriched())
  }) %>% bindEvent(input$PCSF_enrichment_button)
  
  output$PCSF_Enrichment_Clusters_Network <- visNetwork::renderVisNetwork({
    req(input$network_input_type, input$PCSF_enrichment_button)
    
    shiny::validate(need(!is.null(PCSF_enriched()), "There are no enriched clusters to display. Try choosing another database for enrichment, adjusting the cutoff values, or selecting another input."))
    
    pcsfEnrichedContracted(pcsf_enrich_pathway = PCSF_enriched())
  }) %>% bindEvent(input$PCSF_enrichment_button)
  
  ### DOWNLOADS ----------------------------------------------------------------
  
  ## REPORT GENERATION #########################################################
  output$checkrender <- renderText({
    if (identical(rmarkdown::metadata$runtime, "shiny")) {
      TRUE
    } else {
      FALSE
    }
  })
  
  ### RMD REPORT DOWNLOAD/PARAMS -----------------------------------------------
  report_name <- reactive({
    
  })
  
  output$report <- downloadHandler(
    filename = function() {paste0(report_file_name(), ".html")},
    content = function(file) {
      withProgress(message = 'Rendering the HTML report:', detail = "Preparing Data...",
                   
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      {tempReport <- file.path(tempdir(), "report.Rmd")
      file.copy("report.Rmd", tempReport, overwrite = TRUE)
      
      # CHECKS
      message("STARTING HTML REPORT GENERATION")
      
      # Set up parameters to pass to Rmd document
      params <- list(file_name = report_file_name(),
                     project_name = input$report_name,
                     eset = eset(),
                     eset_prenorm = eset_prenorm(),
                     annot = annotation(),
                     type = type(),
                     log_transform = input$log_transform,
                     data_format = data_format(),
                     norm = input$norm_eset,
                     MD_contrasts = report_md_contrasts(),
                     has_linear_model = report_has_linear_model(),
                     linear_model = report_linear_model(),
                     contrasts = report_coef_options(),
                     model_equation = report_model_equation(),
                     include_batch_pca = include_batch_pca(), 
                     batch_pca_data = report_batch_corrected_data(),
                     rendered_by_shiny = TRUE,
                     include_enrichment = input$report_include_enrichment,
                     species = input$species,
                     enrichment = input$report_enrichment_calculation,
                     enrichplots = input$report_enrichplots,
                     enriched = report_enricheds(),
                     geneList = report_geneLists(),
                     group_column = input$group,
                     zero_cutoff = input$zero_cutoff,
                     missing_value_imputation = input$impute_missing,
                     uniprot_annotation = input$uniprot_annotation,
                     render_glimma = input$report_include_glimma,
                     glimma = report_glimma_args(),
                     include_sscore = input$report_include_sscore,
                     sscore = report_sscore_dataframe(),
                     sscore_plots = input$report_sscore_plots,
                     sscore_enriched = report_sscore_enriched(),
                     sscore_geneList = report_sscore_geneLists(),
                     include_pcsf = input$report_include_PCSF,
                     pcsf_input = report_PCSF_input(),
                     pcsf_net = report_PCSF_network(),
                     pcsf_plots = input$report_pcsf_visualizations,
                     pcsf_nodes_net = report_PCSF_nodes_network(),
                     pcsf_influential_net = report_PCSF_influential_network(),
                     pcsf_enriched = report_PCSF_enriched(),
                     pcsf_enriched_network = report_PCSF_enriched_network())
      
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      message("ALL REPORT INPUTS GENERATED")
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv()))}
      )
    }
  )
  
  ### REPORT NAME --------------------------------------------------------------
  report_file_name <- reactive({
    name <- paste0("OmNI_Report_", format(Sys.time(), "%y%m%d%H%M"))
    
    if (input$report_name != ""){
      name <- paste0(make.names(input$report_name), "_OmNI_Report_", format(Sys.time(), "%y%m%d%H%M"))
    } 
    
    return(name)
  })
  
  #### REPORT MD PLOT -----------------------------------------------------------
  report_md_contrasts <- reactive({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$group))
    
    contrasts <- drawMD(annot = annotation(),
                        eset = eset()[[1]],
                        type = type()[i],
                        return_logfc_index = TRUE,
                        logfc_index_choice = NULL)
  })
  
  #### REPORT LIMMA LINEAR MODEL -----------------------------------------------
  report_linear_model <- reactive({
    req((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data))
    linear_models <- list("NONE")
    
    message("INITIALIZED LINEAR MODEL(S) GENERATION")
    try({
      for (i in 1:length(type())){
        linear_models[[i]] <- limmaLM(annot = annotation(),
                                      eset = eset()[[i]],
                                      samples = input$report_linear_factors,
                                      time_series = input$report_include_timeseries,
                                      time_col = input$report_time_col,
                                      time_points_cont = input$report_time_points_cont,
                                      time_points_disc = input$report_time_points_disc,
                                      time = input$report_time_type,
                                      contrast_fit = input$report_contrast_fit,
                                      contrasts_subset = input$report_limma_contrasts,
                                      covariate = input$report_add_covariate,
                                      covariate_col = input$report_covariate_col)
      }
    })
    
    message("---- LINEAR MODEL(S) GENERATED")
    return(linear_models)
  })
  
  report_has_linear_model <- reactive({
    assign("report_linear_model", report_linear_model(), envir = .GlobalEnv)
    if (!(length(report_linear_model()[[1]]) == 1)){
      has_linear_model <- TRUE
    } else {
      has_linear_model <- FALSE
    }
    return(has_linear_model)
  })
  
  report_coef_options <- reactive({
    req((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data))
    
    message("INITIALIZED LM COEF OPTIONS GENERATION")
    coef <- limmaLM(annot = annotation(),
                    eset = eset()[[1]],
                    samples = input$report_linear_factors,
                    coef_names = TRUE,
                    time_series = input$report_include_timeseries,
                    time_col = input$report_time_col,
                    time_points_cont = input$report_time_points_cont,
                    time_points_disc = input$report_time_points_disc,
                    time = input$report_time_type,
                    contrast_fit = input$report_contrast_fit,
                    contrasts_subset = input$report_limma_contrasts,
                    covariate = input$report_add_covariate,
                    covariate_col = input$report_covariate_col)
    
    message("---- LM COEF OPTIONS GENERATED")
    return(coef)
  })
  
  report_model_equation <- reactive({
    message("INITIALIZED LM EQUATION GENERATION")
    equation <- limmaEquation(groups = input$group,
                              include_covariates = input$report_add_covariate,
                              covariate_col = input$report_covariate_col,
                              time_series = input$report_include_timeseries)
    message("---- LM EQUATION GENERATED")
    return(equation)
  })
  
  output$report_linear_factors <- renderUI({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$group))
    options <- unique(annotation()$Group)
    
    selectInput("report_linear_factors",
                label = "Choose at least two groups to include as factors in the model",
                choices = options,
                selected = options,
                multiple = TRUE)
  })
  
  output$report_covariate_col <- renderUI({
    req(input$report_add_covariate == TRUE)
    
    selectizeInput("report_covariate_col",
                   label = "Select up to 3 columns to act as covariates",
                   choices = annot_columns(),
                   multiple = TRUE,
                   options = list(maxItems = 3))
    
  }) %>% bindEvent(input$report_add_covariate)
  
  output$report_choose_contrasts <- renderUI({
    req((input$report_contrast_fit == TRUE), isTruthy(length(input$report_linear_factors) > 1))
    contrasts <- limmaLM(annot = annotation(),
                         eset = eset()[[1]],
                         samples = input$report_linear_factors,
                         time_series = input$report_include_timeseries,
                         time_col = input$report_time_col,
                         time_points_cont = input$report_time_points_cont,
                         time_points_disc = input$report_time_points_disc,
                         time = input$report_time_type,
                         contrast_fit = input$report_contrast_fit,
                         return_contrasts = TRUE,
                         covariate = input$report_add_covariate,
                         covariate_col = input$report_covariate_col)
    
    selectInput("report_limma_contrasts",
                label = "Choose contrasts to include in model",
                choices = contrasts,
                multiple = TRUE)
  })
  
  output$report_choose_time <- renderUI({
    req(input$report_include_timeseries == TRUE)
    
    selectInput("report_time_col",
                label = "Select column containting time series information",
                choices = annot_columns(),
                multiple = FALSE)
  })
  
  output$report_time_type <- renderUI({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$report_include_timeseries))
    
    radioButtons("report_time_type",
                 label = "Select how to include the time variable in the model",
                 choices = c("Continuous" = "continuous",
                             "Discrete" = "discrete"))
  }) %>% bindEvent(input$report_include_timeseries)
  
  output$report_time_points <- renderUI({
    req(isTruthy(input$report_include_timeseries))
    
    if(req(input$report_time_type) == "continuous"){
      number = length(unique(annotation()[,input$report_time_col]))
      numericInput("report_time_points_cont",
                   label = "Choose number of time points to include in the model",
                   min = 0,
                   max = number,
                   value = 3) 
      
    } else if(req(input$report_time_type) == "discrete"){
      options <- unique(annotation()[,input$report_time_col])
      options <- options[-order(options)[1]]
      options <- sort(options)
      selectInput("report_time_points_disc",
                  label = "Choose which time points to include in the model",
                  choices = options,
                  multiple = TRUE)
      
    }
    
  }) %>% bindEvent(input$report_time_type, input$report_include_timeseries, input$report_time_col)
  
  #### REPORT BATCH CORRECTED PCA ----------------------------------------------
  
  include_batch_pca <- reactive({
    include = FALSE
    if (length(input$report_covariate_col) >= 1){
      include = TRUE
    }
    return(include)
  })
  
  report_batch_corrected_data <- reactive({
    req((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data))
    
    if (include_batch_pca() == TRUE) {
      message("INITIALIZED BATCH CORRECTED DATA GENERATION")
      
      data <- list()
      for (i in 1:length(type())){
        type <- list()
        for (j in 1:length(input$report_covariate_col)){
          type[[j]] <- limmaLM(annot = annotation(),
                               eset = eset()[[i]],
                               samples = input$report_linear_factors,
                               time_series = input$report_include_timeseries,
                               time_col = input$report_time_col,
                               time_points_cont = input$report_time_points_cont,
                               time_points_disc = input$report_time_points_disc,
                               time = input$report_time_type,
                               contrast_fit = input$report_contrast_fit,
                               contrasts_subset = input$report_limma_contrasts,
                               covariate = input$report_add_covariate,
                               covariate_col = input$report_covariate_col[j],
                               remove_batch_PCA = TRUE,
                               batch_column = input$report_covariate_col[j])
        }
        data[[i]] <- type
        names(data[[i]]) <- make.names(input$report_covariate_col)
      }
      names(data) <- type()
      
      message("---- BATCH CORRECTED DATA GENERATED")
      return(data)
      
    } else { 
      return(NULL)
    }
  })
  
  #### REPORT GLIMMA ----------------------------------------------------------
  report_glimma_args <- reactive({
    message("INITIALIZED GLIMMA ARG GENERATION")
    
    args <- list()
    
    try({
      for (i in 1:length(type())){
        
        arg <- list()
        for (j in 1:length(report_coef_options())){
          
          table <- makeTopTable(report_linear_model()[[i]],
                                coef = report_coef_options()[[j]],
                                coef_options = report_coef_options())
          
          arg[[j]] <- interactiveVolcano(eset = eset()[[i]],
                                         fit = report_linear_model()[[i]], 
                                         top_table = table, 
                                         type = type()[i], 
                                         coef = report_coef_options()[[j]])
          
        }
        
        args[[i]] <- arg
      }
    })
    message("---- GLIMMA ARGS GENERATED")
    return(args)
  })
  
  #### REPORT ENRICHMENT --------------------------------------------------------
  output$report_database_header <- renderUI({
    req(input$report_include_enrichment)
    h6("Database(s)")
  })
  
  output$report_enrichment_calculation <- renderUI({
    req(input$report_include_enrichment)
    
    tagList(
      h6("Algorithm"),
      radioButtons("report_enrichment_calculation",
                   label = "Enrichment algorithm to perform",
                   choices = c("Ranked Enrichment Analysis" = "ranked",
                               "Unranked Over-Representation Analysis (ORA)" = "unranked")),
      hr(),
      h6("Plots")
    )
  })
  
  output$report_enrichment_plots <- renderUI({
    req(input$report_include_enrichment)
    
    checkboxGroupInput("report_enrichplots",
                       label = "Include enrichment plots",
                       choices = c("Dot Plot" = "dot",
                                   "Gene Set Network" = "cnet",
                                   "Enrichment Map" = "emap",
                                   "UpSet Plot" = "upset",
                                   "Tree Diagram" = "tree"),
                       selected = c("dot"))
  })
  
  report_geneLists <- reactive({
    if (input$report_include_enrichment == TRUE){
      message("INITIALIZED GENE LIST GENERATION")
      
      GSEA_geneList <- list()
      PKSEA_geneList <- list()
      
      for (j in 1:length(type())){
        sub_geneList <- list()
        
        if (data_format()[j] == "ProteinGroups" | data_format()[j] == "PhosphoSites"){
          for (i in 1:length(report_coef_options())){
            genes <- calculateGeneList(fit = report_linear_model()[[j]],
                                       coef = report_coef_options()[i],
                                       gmt = "GSEA")
            sub_geneList[[paste0(report_coef_options()[i])]] <- genes 
          }
          GSEA_geneList[[j]] <- sub_geneList
          names(GSEA_geneList)[j] <- paste(type()[j], "GSEA", sep = "_")
        }
        
        if (data_format()[j] == "PhosphoSites"){
          for (i in 1:length(report_coef_options())){
            genes <- calculateGeneList(fit = report_linear_model()[[j]],
                                       coef = report_coef_options()[i],
                                       gmt = "PHOSPHO")
            sub_geneList[[paste0(report_coef_options()[i])]] <- genes 
          }
          PKSEA_geneList[[j]] <- sub_geneList
          names(PKSEA_geneList)[j] <- paste(type()[j], "PKSEA", sep = "_")
        }
      }
      
      geneList <- c(GSEA_geneList, PKSEA_geneList)
      geneList <- geneList[!sapply(geneList, is.null)]
      
      # assign("enrich_geneList", geneList, envir = .GlobalEnv)
      # assign("report_coef_options", report_coef_options(), envir = .GlobalEnv)
      
      message("---- GENE LIST GENERATED")
      return(geneList)
    } else { 
      return(NULL)
    }
    
  })
  
  report_enricheds <- reactive({
    if (input$report_include_enrichment == TRUE){
      message("INITIALIZED ENRICHMENT")
      
      enriched <- list()
      
      if (input$species == "human"){
        Gene_GMTs <- c("HUMAN_KEGG", "HUMAN_Reactome", "HUMAN_WikiPathways", "HUMAN_GOBP",
                       "HUMAN_GOCC", "HUMAN_GOMF", "HUMAN_Hallmark", "HUMAN_MSigDB", "HUMAN_MOMENTABiocycEXP",
                       "HUMAN_MOMENTAMFNEXP")
        Phospho_GMT <- c("HUMAN_PHOSPHO_UNIPROT")
        
      } else if (input$species == "mouse"){
        Gene_GMTs <- c("MOUSE_Reactome", "MOUSE_GOBP", "MOUSE_GOCC", "MOUSE_GOMF", "MOUSE_WikiPathways",
                       "MOUSE_MOMENTABiocycEXP")
        Phospho_GMT <- c("MOUSE_PHOSPHO_UNIPROT")
        
      } else if (input$species == "rat"){
        Gene_GMTs <- c("RAT_Reactome", "RAT_GOBP", "RAT_GOCC", "RAT_GOMF", "RAT_WikiPathways")
        Phospho_GMT <- c("RAT_PHOSPHO_UNIPROT")
        
      } else if (input$species == "celegans") {
        Gene_GMTs <- c("CELEGANS_KEGG", "CELEGANS_WikiPathways", "CELEGANS_GOBP", "CELEGANS_GOCC", "CELEGANS_GOMF")
        Phospho_GMT <- "NONE"
        
      } else if (input$species == "fruitfly") {
        Gene_GMTs <- c("FRUITFLY_KEGG", "FRUITFLY_WikiPathways", "FRUITFLY_GOBP", "FRUITFLY_GOCC", "FRUITFLY_GOMF", 
                       "FRUITFLY_MOMENTABiocycEXP")
        Phospho_GMT <- "NONE"
        
      } else if (input$species == "zebrafish") {
        Gene_GMTs <- c("ZEBRAFISH_KEGG", "ZEBRAFISH_WikiPathways", "ZEBRAFISH_GOBP", "ZEBRAFISH_GOCC", "ZEBRAFISH_GOMF",)
        Phospho_GMT <- "NONE"
        
      } else if (input$species == "yeast") {
        Gene_GMTs <- c("YEAST_KEGG", "YEAST_WikiPathways", "YEAST_GOBP", "YEAST_GOCC", "YEAST_GOMF", "YEAST_MOMENTABiocycEXP")
        Phospho_GMT <- "NONE"
        
      }
      
      withProgress(message = 'Performing Enrichment:', detail = "This may take a while...", value = 0, {
      for (j in 1:length(report_geneLists())){
        setProgress(value = (j / length(report_geneLists()) - ((1 / length(report_geneLists())) / 2)), 
                    detail = paste("Enriching gene list ", j, " of ", length(report_geneLists())))
        contrast_enriched <- list()
        phospho_enriched <- list()
        for (i in 1:length(report_geneLists()[[j]])){
          if (grepl("_GSEA", names(report_geneLists())[j])){
            for(k in 1:length(Gene_GMTs)){
              message(names(report_geneLists())[j], Gene_GMTs[k])
              enrich <- clusterProfilerEnrichment(geneList = report_geneLists()[[j]][[i]],
                                                  gmt = Gene_GMTs[k],
                                                  enrichment = input$report_enrichment_calculation,
                                                  pval_cutoff = 1)
              
              if (!is.atomic(enrich)){
                enrich@result <- enrich[enrich@result$pvalue <= 0.05]
                contrast_enriched[[paste0(Gene_GMTs[k])]] <- enrich
              }
            }
            enriched[[paste0(sub("_GSEA", "", names(report_geneLists())[j]), "_", names(report_geneLists()[[j]])[i])]] <- contrast_enriched
          } else if (grepl("_PKSEA", names(report_geneLists())[j]) & Phospho_GMT != "NONE"){
            message(names(report_geneLists())[j])
            enrich <- clusterProfilerEnrichment(geneList = report_geneLists()[[j]][[i]],
                                                gmt = Phospho_GMT,
                                                enrichment = input$report_enrichment_calculation,
                                                pval_cutoff = 1)
            
            if (!is.atomic(enrich)){
              enrich@result <- enrich[enrich@result$pvalue <= 0.05]
              enriched[[paste0(sub("_PKSEA", "", names(report_geneLists())[j]), "_", names(report_geneLists()[[j]])[i])]][paste0("HUMAN_PHOSPHO_UNIPROT")] <- enrich
            }
          }
        }
      }
      })
      
      assign("enriched", enriched, envir = .GlobalEnv)
      message("---- ENRICHMENT COMPLETE")
      return(enriched)
    } else {
      return(NULL)
    }
  })
  
  #### REPORT S-SCORE ----------------------------------------------------------
  
  output$report_sscore_plots <- renderUI({
    req(input$report_include_sscore == TRUE)
    
    checkboxGroupInput("report_sscore_plots",
                       label = "Choose visualizations to include",
                       choices = c("Volcano Plot" = "volcano",
                                   "Venn Diagram" = "venn"),
                       selected = "volcano")
  })
  
  output$report_sscore_dataset <- renderUI({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$report_include_sscore))
    
    options <- type()
    selectInput("report_sscore_datasets",
                label = "Choose at least 2 datasets to include",
                choices = options,
                multiple = TRUE)
  })
  
  report_sscore_contrasts <- reactive({
    req(input$report_include_sscore)
    contrasts <- limmaLM(annot = annotation(),
                         eset = eset()[[1]],
                         samples = input$report_linear_factors,
                         time_series = input$report_include_timeseries,
                         time_col = input$report_time_col,
                         time_points_cont = input$report_time_points_cont,
                         time_points_disc = input$report_time_points_disc,
                         time = input$report_time_type,
                         contrast_fit = TRUE,
                         covariate = input$report_add_covariate,
                         covariate_col = input$report_covariate_col,
                         return_contrasts = TRUE)
    return(contrasts)
  })
  
  report_sscore_toptable_list <- reactive({
    req(((isTruthy(input$data_file) && isTruthy(input$annotation_file)) || isTruthy(input$use_example_data)) && isTruthy(input$report_limma_contrasts) && isTruthy(input$report_include_sscore))
    message("INITIALIZED SSCORE TOPTABLE GENERATION")
    
    linear_models <- list()
    top_tables <- list()
    
    datasets <- c()
    for (i in input$report_sscore_datasets){
      datasets <- c(datasets, which(type() == i))
    }
    
    for (i in 1:length(datasets)){
      linear_models[[i]] <- limmaLM(annot = annotation(),
                                    eset = eset()[[datasets[i]]],
                                    samples = input$report_linear_factors,
                                    time_series = input$report_include_timeseries,
                                    time_col = input$report_time_col,
                                    time_points_cont = input$report_time_points_cont,
                                    time_points_disc = input$report_time_points_disc,
                                    time = input$report_time_type,
                                    contrast_fit = TRUE,
                                    contrasts_subset = input$report_limma_contrasts,
                                    covariate = input$report_add_covariate,
                                    covariate_col = input$report_covariate_col)
    }
    
    coefs <- input$report_limma_contrasts
    for (j in 1:length(coefs)){
      top_tables_subset <- list()
      
      for (i in 1:length(linear_models)){
        fit <- linear_models[[i]]
        n <- nrow(fit$coefficients)
        
        end_col <- ncol(limma::topTable(fit, adjust = "BH", coef = coefs[j]))
        start_col <- which((colnames(limma::topTable(fit, adjust = "BH", coef = coefs[j]))) == "feature_identifier")
        top_tables_subset[[i]] <- limma::topTable(fit, adjust = "BH", number = n, coef = coefs[j])
      }
      names(top_tables_subset) <- input$report_sscore_datasets
      
      top_tables[[j]] <- top_tables_subset
    }
    names(top_tables) <- coefs
    
    message("---- S-SCORE TOPTABLE(S) GENERATED")
    return(top_tables);
  })
  
  report_sscore_data_list <- reactive({
    req(input$report_include_sscore)
    message("INITIALIZED SSCORE DATA LIST GENERATION")
    
    toptable_list <- report_sscore_toptable_list()
    
    data_formats <- c()
    for (i in input$report_sscore_datasets){
      data_formats <- c(data_formats, data_format()[which(type() == i)])
    }
    
    data_list <- list()
    for (i in 1:length(input$report_limma_contrasts)){
      data_list[[i]] <- makeDataList(top_table_list = toptable_list[[i]],
                                     data_format = data_formats)
    }
    names(data_list) <- input$report_limma_contrasts
    
    message("---- S-SCORE DATA LIST(S) GENERATED")
    
    return(data_list);
  })
  
  report_sscore_dataframe <- reactive({
    if (input$report_include_sscore == TRUE){
    message("INITIALIZED SSCORE INTEGRATION")
    
    withProgress(message = 'Integrating via S-Score:', detail = "This may take a while...", value = 0, {
      sscore <- list()
      for (i in 1:length(report_sscore_data_list())){
        setProgress(value = (i / length(report_sscore_data_list()) - ((1 / length(report_sscore_data_list())) / 2)), 
                    detail = paste("Integrating contrast ", i, " of ", length(report_sscore_data_list())))
        sscore[[i]] <- sscoreIntegration(data_list = report_sscore_data_list()[[i]])
      }
      names(sscore) <- paste0(paste0(input$report_sscore_datasets, collapse = "_"), "_", input$report_limma_contrasts)
    })
    
    # assign("report_sscore_dataframe", sscore, envir = .GlobalEnv)
    message("---- SSCORE INTEGRATION COMPLETE")
    return(sscore);
    } else {
      return(NULL)
    }
  })
  
  ##### ENRICHMENT -------------------------------------------------------------
  report_sscore_geneLists <- reactive({
    if ((input$report_include_sscore == TRUE) & (input$report_include_enrichment == TRUE)){
      i <- which(grepl(paste0(input$report_limma_contrasts, collapse = "|"), names(report_sscore_dataframe())))
      geneLists <- list()
      for (j in 1:length(i)){
        geneLists[[paste0(names(report_sscore_dataframe())[j])]] <- formatSscoreGeneList(sscore_output = report_sscore_dataframe()[[j]])
      }
      
      # assign("report_sscore_geneLists", geneLists, envir = .GlobalEnv)
    } else {
      geneLists <- NULL
    }
    
    return(geneLists)
  })

  report_sscore_enriched <- reactive({
    if ((input$report_include_sscore == TRUE) & (input$report_include_enrichment == TRUE)){
      message("INITIALIZING S-SCORE ENRICHMENT")
      GMTs <- c("HUMAN_GOBP", "MOUSE_GOBP", "RAT_GOBP", "YEAST_GOBP", "FRUITFLY_GOBP",
                "CELEGANS_GOBP", "ZEBRAFISH_GOBP")
      
      GMT <- GMTs[grep(toupper(input$species), GMTs)]
      
      enriched <- list()
      for (i in 1:length(report_sscore_geneLists())){
        if (!grepl("CHEBI:", paste0(names(report_sscore_geneLists()[[i]]), collapse = "; "))){
          enriched[[paste0(names(report_sscore_geneLists())[i])]] <- runSscoreEnrichment(geneList = report_sscore_geneLists()[[i]],
                                                                                         gmt = GMT,
                                                                                         pval_cutoff = 0.05)
          message("NON-METABOLITE ENRICHMENT PERFORMED: ", names(report_sscore_geneLists())[i])
        } else if (grepl("CHEBI:", paste0(names(report_sscore_geneLists()[[i]]), collapse = "; ")) && input$species == "human") {
          enriched[[paste0(names(report_sscore_geneLists())[i])]] <- runSscoreEnrichment(geneList = report_sscore_geneLists()[[i]],
                                                                                         gmt = "HUMAN_METABOLGENE",
                                                                                         pval_cutoff = 0.05)
          message("METABOLITE ENRICHMENT PERFORMED: ", names(report_sscore_geneLists())[i])
        } else {
          enriched[[paste0(names(report_sscore_geneLists())[i])]] <- NULL
          message("NO ENRICHED PATHWAYS: ", names(report_sscore_geneLists())[i])
        }
      }
      
      # assign("SSCORE_enriched", enriched, envir = .GlobalEnv)
      message("---- SSCORE ENRICHMENT COMPLETE")
    } else {
      enriched <- NULL
    }
    return(enriched)
  })
  
  #### REPORT PCSF -------------------------------------------------------------
  
  output$report_pcsf_plots <- renderUI({
    req(input$report_include_PCSF == TRUE)
    
    checkboxGroupInput("report_pcsf_visualizations",
                       label = "Choose networks to render",
                       choices = c("Node Network" = "node",
                                   "Influential Network" = "influential"),
                       selected = c("node", "influential"))
  })
  
  report_PCSF_input_linear_model <- reactive({
    message("INITIALIZED PCSF LM INPUT GENERATION")
    coefs <- input$report_limma_contrasts
    
    PCSF_input <- list()
    for (j in 1:length(coefs)){
      top_tables_subset <- list()
      
      for (i in 1:length(report_linear_model())){
        fit <- report_linear_model()[[i]]
        
        # TODO: Add ability to use metabolite feature ID, i.e. HMDB, InChIKey in PCSF
          # NOTE: Would need to be converted to CHEBI for enrichment with metabolite
          # enrichment file.
        top_tables_subset[[i]] <- limma::topTable(fit, adjust = "BH", number = nrow(fit$coefficients), coef = coefs[j])
        print("Gene" %in% colnames(top_tables_subset[[i]]))
        if (!("Gene" %in% colnames(top_tables_subset[[i]]))) {top_tables_subset[[i]] <- matrix(nrow = 0, ncol = 3)}
        try({
          end_col <- ncol(top_tables_subset[[i]])
          start_col <- which((colnames(top_tables_subset[[i]])) == "Gene")
          top_tables_subset[[i]] <- top_tables_subset[[i]][,seq(start_col, end_col)]
          top_tables_subset[[i]] <- top_tables_subset[[i]][(abs(top_tables_subset[[i]]$logFC) >= 0 & 
                                                              top_tables_subset[[i]]$adj.P.Val <= 0.05), ]
          
          top_tables_subset[[i]] <- top_tables_subset[[i]][, c("Gene", "logFC", "adj.P.Val")]
          if(nrow(top_tables_subset[[i]]) > 100){
            top_tables_subset[[i]] <- top_tables_subset[[i]][1:100,]
          }
          colnames(top_tables_subset[[i]]) <- c("gene_symbol", "logfc", "adj_pval")
        })
      }
      names(top_tables_subset) <- type()
      
      PCSF_input[[j]] <- top_tables_subset
    }
    names(PCSF_input) <- coefs
    
    message("-- PCSF LM INPUT GENERATED")
    return(PCSF_input)
  })
  
  report_PCSF_input_sscore <- reactive({
    message("INITIALIZED PCSF SS INPUT GENERATION")
    coefs <- input$report_limma_contrasts
    
    PCSF_input <- list()
    for (i in 1:length(coefs)){
      dataframe <- report_sscore_dataframe()[[i]]
      dataframe <- dataframe[(abs(dataframe$sscore) >= 0 &
                                dataframe$sscore_adj_pval <= 0.05), ]
      dataframe[is.na(dataframe)] <- ""
      
      PCSF_input[[i]] <- dataframe[, c("sscore", "sscore_adj_pval")]
      PCSF_input[[i]]$gene_symbol <- paste0(dataframe$gene_symbol, dataframe$chebi_id)
      PCSF_input[[i]] <- PCSF_input[[i]][, c("gene_symbol", "sscore", "sscore_adj_pval")]
      colnames(PCSF_input[[i]]) <- c("gene_symbol", "sscore", "adj_pval")
    }
    names(PCSF_input) <- paste0("Sscore_", paste0(input$report_sscore_datasets, collapse = "_"), "_", coefs)
    
    message("-- PCSF SS INPUT GENERATED")
    return(PCSF_input)
  })
  
  report_PCSF_input <- reactive({
    try(sscore <- input$report_include_sscore)
    if (input$report_include_PCSF == TRUE){
      message("FINALIZING PCSF INPUT")
      pcsf_input_1 <- list()
      k <- 1
      for (i in 1:length(report_PCSF_input_linear_model())){
        for (j in 1:length(report_PCSF_input_linear_model()[[i]])){
          if (dim(report_PCSF_input_linear_model()[[i]][[j]])[1] > 0){
            print(str(report_PCSF_input_linear_model()[[i]][[j]]))
            input <- report_PCSF_input_linear_model()[[i]][[j]]
            print(str(input))
            pcsf_input_1[[k]] <- input
            names(pcsf_input_1)[k] <- paste0(names(report_PCSF_input_linear_model()[[i]])[j], "_", names(report_PCSF_input_linear_model())[i], sep = "")
            k <- k + 1
          }
        }
      }
      
      if (sscore == TRUE){
        pcsf_input_2 <- list()
        for (i in 1:length(report_PCSF_input_sscore)){
          if (dim(report_PCSF_input_sscore()[[i]])[1] > 0){
            input <- report_PCSF_input_sscore()[[i]]
            print(str(input))
            pcsf_input_2[[i]] <- input
            names(pcsf_input_2)[i] <- names(report_PCSF_input_sscore())[i]
          }
        }
        
        pcsf_input <- c(pcsf_input_1, pcsf_input_2)
      } else {
        pcsf_input <- pcsf_input_1
      }
      
      # assign("PCSF_report_input", pcsf_input, envir = .GlobalEnv)
      message("---- PCSF INPUT FINALIZED")
      return(pcsf_input)
    } else {
      return(NULL)
    }
  })
  
  report_PCSF_network <- reactive({
    if (input$report_include_PCSF == TRUE){
      message("INITIALIZING PCSF NETWORK GENERATION")
      
      withProgress(message = 'Generating PCSF:', detail = "This may take a while...", value = 0, {
        pcsf_net <- list()
        for (i in 1:length(report_PCSF_input())){
          setProgress(value = (i / length(report_PCSF_input()) - ((1 / length(report_PCSF_input())) / 2)), 
                      detail = paste("Making PCSF Network ", i, " of ", length(report_PCSF_input())))
          pcsf_net[[i]] <- makePCSFNetwork(pcsf_input_data = report_PCSF_input()[[i]],
                                           pcsf_nval = 10)
        }
        names(pcsf_net) <- names(report_PCSF_input())
      })
      
      # assign("PCSF_report_network", pcsf_net, envir = .GlobalEnv)
      message("---- PCSF NETWORK(S) GENERATED")
      return(pcsf_net)
    } else {
      return(NULL)
    }
  })
  
  report_PCSF_nodes_network <- reactive({
    if (input$report_include_PCSF == TRUE){
      network <- list()
      for (i in 1:length(report_PCSF_network())){
        network[[i]] <- PCSFVisNodes(pcsf_net = report_PCSF_network()[[i]],
                                     pcsf_input_data = report_PCSF_input()[[i]]) 
      }
      names(network) <- names(report_PCSF_input())
      
      # assign("PCSF_nodes_network", network, envir = .GlobalEnv)
      message("---- PCSF NODE NETWORK VISGRAPH(S) GENERATED")
      return(network)
    }
  })
  
  report_PCSF_influential_network <- reactive({
    if (input$report_include_PCSF == TRUE){
      network <- list()
      for (i in 1:length(report_PCSF_network())){
        network[[i]] <- PCSFVisInfluential(pcsf_net = report_PCSF_network()[[i]])
      }
      names(network) <- names(report_PCSF_input())
      
      # assign("PCSF_influential_network", network, envir = .GlobalEnv)
      message("---- PCSF INFLUENTIAL NETWORK VISGRAPH(S) GENERATED")
      return(network)
    }
  })
  
  report_PCSF_enriched <- reactive({
    if (input$report_include_PCSF == TRUE){
      if (input$report_include_enrichment == TRUE){
        message("INITIALIZING PCSF ENRICHMENT")
        
        GMTs <- c("HUMAN_GOBP", "MOUSE_GOBP", "RAT_GOBP", "YEAST_GOBP", "FRUITFLY_GOBP",
                  "CELEGANS_GOBP", "ZEBRAFISH_GOBP")
        
        GMT <- GMTs[grep(toupper(input$species), GMTs)]
        
        enriched <- list()
        for (i in 1:length(report_PCSF_network())){
          if (!grepl("CHEBI:", paste0(report_PCSF_input()[[i]]$gene_symbol, collapse = "; "))){
            enriched[[paste0(names(report_PCSF_input())[i])]] <- pcsfRunEnrichment(pcsf_net = report_PCSF_network()[[i]],
                                                                                   gmt = GMT)
            message("NON-METABOLITE ENRICHMENT PERFORMED: ", names(report_PCSF_input())[i])
          } else if (grepl("CHEBI:", paste0(report_PCSF_input()[[i]]$gene_symbol, collapse = "; ")) && input$species == "human") {
            enriched[[paste0(names(report_PCSF_input())[i])]] <- pcsfRunEnrichment(pcsf_net = report_PCSF_network()[[i]],
                                                                                   gmt = "HUMAN_METABOLGENE")
            message("METABOLITE ENRICHMENT PERFORMED: ", names(report_PCSF_input())[i])
          } else {
            enriched[[paste0(names(report_PCSF_input())[i])]] <- NULL
            message("NO ENRICHED PATHWAYS: ", names(report_PCSF_input())[i])
          }
        }
        
        # assign("PCSF_enriched", enriched, envir = .GlobalEnv)
        message("---- PCSF ENRICHMENT COMPLETE")
        return(enriched)
      }
    }
  }) 
  
  report_PCSF_enriched_network <- reactive({
    if (input$report_include_PCSF == TRUE){
      if (input$report_include_enrichment == TRUE){
        enriched_network <- list()
        for (i in 1:length(report_PCSF_enriched())){
          enriched_network[[i]] <- pcsfEnrichedSubnet(pcsf_enrich_pathway = report_PCSF_enriched()[[i]])
        }
        names(enriched_network) <- names(report_PCSF_enriched())
        # assign("PCSF_enriched_network", enriched_network, envir = .GlobalEnv)
        message("---- PCSF ENRICHED NETWORK GENERATED")
        return(enriched_network)
      }
    }
  })
  
  ## CLOSE OUT #################################################################
  
  session$onSessionEnded(function() {
    # delete any KEGG hsa pathway files that were generated
    files <- list.files(pattern = "hsa")
    to_be_deleted <- grep("hsa", files, value = T)
    file.remove(to_be_deleted)
    
    # delete Glimma HTML/CS/JSS from glimma-plots folder
    glimma_files <- list.files("glimma-plots", full.names = TRUE)
    unlink(glimma_files, recursive = TRUE)
    
    # fill out row in user_logTable with use times
    isolate({
      users$logTable[users$logTable$id == session$token,
                     "logout"] = as.character(Sys.time())
      
      table <- read.csv("user_logTable.csv")
      table <- rbind(table, users$logTable)
      write.csv(table, "user_logTable.csv", row.names = FALSE)
    })
  
    # stop the app
    stopApp()
  })
}
