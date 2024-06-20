
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
source("intensityNorm.R")
source("normCorrPlots.R")
source("saveXLSX.R")
source("enrichment.R")

# SERVER #######################################################################

server <- function(input, output, session) {
  
  # INCREASE MAX FILE SIZE TO 30KB
  options(shiny.maxRequestSize = 30 * 1024 ^ 2)
  
  ## DATA ######################################################################
  
  ### REACTIVE OBJECTS ---------------------------------------------------------
  data <- reactive({
    req(input$data_file)
    
    utils::read.delim(input$data_file$datapath, header = TRUE, sep = "\t")
  })
  
  ### RENDER UI ----------------------------------------------------------------
  output$dataOverview <- renderUI({
    req(isTruthy(input$data_file) || isTruthy(input$annotation_file))
    h3("Data Overview")
  })
  
  output$dataOverview_help <- renderUI({
    req(isTruthy(input$data_file) || isTruthy(input$annotation_file))
    helpText("Interpreting the data summary: ", "\n", "The top line identifies 
               the structure of the entry, the number of observations (or rows) and the number 
               of variables (or columns). Underneath there are column names on the left followed by a colon, 
               then the variable type (i.e. character, number, factor) followed by 
               the first four or so entries in the column.")
  })
  
  output$dataFile <- renderUI({
    req(isTruthy(input$data_file))
    h5("Data File")
  })
  
  output$filesHeader <- renderUI({
    req(isTruthy(input$data_file) || isTruthy(input$annotation_file))
    h4("File(s)")
  })
  
  ### PRINT/TABLES/PLOTS -------------------------------------------------------
  output$summary_data <- renderPrint({
    utils::str(data())
  })
  
  output$file_table <- renderTable(input$data_file)
  
  ## ANNOTATION ################################################################
  
  ### REACTIVE OBJECTS ---------------------------------------------------------
  annotation_initial <- reactive({
    req(input$annotation_file)
    
    ext <- tools::file_ext(input$annotation_file$datapath)
    validate(need(ext == "xlsx", "Please upload a .xlsx for the annotation file."))
    
    annot <- data.frame(openxlsx::read.xlsx(input$annotation_file$datapath, 2, colNames = TRUE))
    
  })
  
  annotation <- reactive({
    req(input$group, input$annotation_file)
    formatAnnotation(annotation_initial(),
                     group = input$group)
  })
  
  type <- reactive({
    annot <- data.frame(openxlsx::read.xlsx(input$annotation_file$datapath, 1, colNames = FALSE))
    annot <- annot[-(1:10),]
    annot <- annot[,-3]
    colnames(annot) <- c("File", "Name")
    annot <- na.omit(annot)
    
    file <- input$data_file$name
    row <- annot[grep(file, annot$File), ]
    type <- row$Name
  })
  
  ### RENDER UI ----------------------------------------------------------------
  output$select_group <- renderUI({
    req(input$annotation_file)
    
    annot <- data.frame(openxlsx::read.xlsx(input$annotation_file$datapath, 2, colNames = TRUE))
    annot_columns <- colnames(annot)
    
    annot_names <- data.frame(openxlsx::read.xlsx(input$annotation_file$datapath, 1, colNames = FALSE))
    annot_names <- annot_names[-(1:10),]
    annot_names <- annot_names[,-3]
    colnames(annot_names) <- c("File", "Name")
    annot_names <- na.omit(annot_names)
    annot_names <- as.list(annot_names$Name)
    
    annot_columns <- annot_columns[!annot_columns == "SampleName"]
    annot_columns <- annot_columns[!annot_columns == annot_names]
    
    tagList(
      selectInput("group",
                  label = "Select column(s) to form groups by",
                  choices = annot_columns,
                  multiple = TRUE),
      helpText("Select multiple options to group by more than one factor."),
      br(),
      br()
    )
  })
  
  output$annotationFile <- renderUI({
    req(isTruthy(input$annotation_file))
    h5("Annotation File")
  })
  
  output$typeHeader <- renderUI({
    req(input$data_file, input$annotation_file)
    h6("Dataset Name")
  })
  
  ### PRINT/TABLES/PLOTS -------------------------------------------------------
  output$summary_annotation <- renderPrint({
    utils::str(annotation_initial())
  })
  
  output$annotation_table <- renderTable(input$annotation_file)
  
  output$type <- renderPrint({
    req(input$annotation_file, input$data_file)
    return(cat(type()))

  })
  
  ## INPUT CHECKS ##############################################################
  
  output$datainput_check <- renderPrint({
    validate(need(!is.null(input$data_file), "No data file input."))
    return(cat("Data file uploaded."))
  })
  
  output$annotationinput_check <- renderPrint({
    validate(need(!is.null(input$annotation_file), "No annotation file input"))
    return(cat("Annotation file uploaded."))
  })
  
  output$annotationformat_check <- renderPrint({
    req(input$annotation_file)
    validate(need(length(type()) > 0, "Issue with annotation formatting, unable to find dataset name."))
    return(cat(paste("Dataset name found: ", type(), ".", sep = "")))
  })
  
  output$groupinput_check <- renderPrint({
    req(input$annotation_file)
    validate(need(!is.null(input$group), "Group columns not selected."))
    return(cat(paste("Group column(s) selected: ", input$group, ".", sep = "")))
  })
  
  output$dataformat_check <- renderPrint({
    validate(need(input$data_format != "none", "No data format selected."))
    return(cat(paste("Data format selected: ", input$data_format, ".", sep = "")))
  })
  
  ## EXPRESSION SET ############################################################
  
  ### REACTIVE OBJECTS ---------------------------------------------------------
  eset <- reactive({
    req(input$data_file, input$annotation_file, (input$data_format != "none"), input$group)
    eset_obj = eset_prenorm()
    
    intensityNorm(eset = eset_obj,
                  norm = input$norm_eset,
                  type = type(),
                  zero_cutoff = input$zero_cutoff)
  })
  
  ### PRINT/TABLE/PLOTS --------------------------------------------------------
  output$view_eset <- DT::renderDataTable({
    req(input$data_file, input$annotation_file, (input$data_format != "none"), input$group)
    
    exprs_eset = exprs(eset())
    
    DT::datatable(exprs_eset, rownames = TRUE)  %>% 
      DT::formatRound(columns = c(1:ncol(exprs_eset)), digits = 3)
    
  })
  
  output$eset_table_download <- downloadHandler(
    filename = function() {
      paste("eset_table_", input$norm_eset, "_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file){
      write.csv(exprs_eset(), file, row.names = TRUE)
    }
  )
  
  ## NORMALIZATION #############################################################
  
  ### REACTIVE OBJECTS ---------------------------------------------------------
  eset_prenorm <- reactive({
    eset_obj <- makeEset(data(), 
                         annotation(), 
                         type = type(), 
                         data_format = input$data_format,
                         log_transform = input$log_transform,
                         uniprot_annotation = input$uniprot_annotation)
  })
  
  normal_plot <- reactive({
    eset_prenorm()
    
    plot <- intensityNormPlots(eset = eset_prenorm(),
                               norm = input$norm_eset, 
                               type = type(), 
                               plottype = input$norm_plottype,
                               ma_array = input$ma_array,
                               zero_cutoff = input$zero_cutoff)
    plot
  }) %>% bindEvent(input$normal_button)
  
  ### RENDER UI ----------------------------------------------------------------
  output$ma_array_choose <- renderUI({
    req(input$norm_plottype == "MA")
    
    list <- 1:ncol(eset())
    selectInput("ma_array",
                 label = "Choose column to generate array",
                 choices = list)
  }) %>% bindEvent(input$norm_plottype)
  
  ### PRINT/TABLE/PLOTS --------------------------------------------------------
  output$normal_plot <- renderPlot({
    normal_plot()
  })
  
  output$normal_download <- downloadHandler(
    filename = function() {
      paste("normplot_", type(), "_", Sys.Date(), ".", input$normal_extension , sep = "")
    },
    content = function(file){
      if(input$normal_extension == "png"){
        png(file, width = 1200, height = 700, units = "px", res = 150)
        eset_obj = makeEset(data(), 
                            annotation(), 
                            type = type(), 
                            data_format = input$data_format,
                            log_transform = input$log_transform,
                            uniprot_annotation = input$uniprot_annotation)
        
        intensityNormPlots(eset = eset_obj,
                           norm = input$norm_eset, 
                           type = type(), 
                           plottype = input$norm_plottype,
                           ma_array = input$ma_array,
                           zero_cutoff = input$zero_cutoff)
        dev.off()
      } else if(input$normal_extension == "pdf"){
        pdf(file, width = 12, height = 7)
        eset_obj = makeEset(data(), 
                            annotation(), 
                            type = type(), 
                            data_format = input$data_format,
                            log_transform = input$log_transform,
                            uniprot_annotation = input$uniprot_annotation)
        
        intensityNormPlots(eset = eset_obj,
                           norm = input$norm_eset, 
                           type = type(), 
                           plottype = input$norm_plottype,
                           ma_array = input$ma_array,
                           zero_cutoff = input$zero_cutoff)
        dev.off()
      }
      
      # ggplot2::ggsave(file, normal_plot(), device = input$normal_extension, width = 10, height = 4)
    }
  )
  
  ## CORRELATION ###############################################################
  
  ### REACTIVE OBJECTS ---------------------------------------------------------
  correlation_plot <- reactive({
    variationPlot(eset(),
                  plot = TRUE,
                  number = input$corr_number,
                  shape = input$corr_shape,
                  neg_color = input$corr_neg_color,
                  pos_color = input$corr_pos_color,
                  type = input$corr_type,
                  title = input$corr_title)
  }) %>% bindEvent(input$corr_button)
  
  ### PRINT/TABLE/PLOTS --------------------------------------------------------
  output$corr_plot <- renderPlot({
    correlation_plot()
  }, width = 800, height = 800)
  
  output$corr_plot_download <- downloadHandler(
    filename = function() {
      paste("corrplot_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file){
      ggsave(file, correlation_plot(), height = 11, width = 11)
    }
  )
  
  ## PCA #######################################################################
  
  ### REACTIVE OBJECTS ---------------------------------------------------------
  PC_plot <- reactive({
    req(input$data_file, input$annotation_file, (input$data_format != "none"), input$group)
    
    drawPCA(eset(), 
            x_axis = input$PCA_xaxis, 
            y_axis = input$PCA_yaxis, 
            type = type(), 
            color = input$PCA_color, 
            include_densities = input$include_densities, 
            shapes = input$PCA_shapes,
            title_add = input$PCA_title,
            add_labels = input$PCA_labels,
            add_ellipse = input$PCA_ellipse)
    
  }) %>% bindEvent(input$PCA_button)
  
  PC_variance_plot <- reactive({
    req(input$data_file, input$annotation_file, (input$data_format != "none"), input$group)
    
    drawPCApercentages(eset())
  })
  
  PCs <- reactive ({
    data = t(exprs(eset()))
    PC_data <- stats::prcomp(data)
    percent_variance <- summary(PC_data)$importance["Proportion of Variance",] * 100
    PCs <- rownames(data.frame(percent_variance))
  })
  
  ### RENDER UI ----------------------------------------------------------------
  output$PC_xaxis <- renderUI({
    req(input$data_file, input$annotation_file, (input$data_format != "none"), input$group)
    
    selectInput("PCA_xaxis",
                label = "X-Axis",
                choices = PCs(),
                selected = "PC1")
  })
  
  output$PC_yaxis <- renderUI({
    req(input$data_file, input$annotation_file, (input$data_format != "none"), input$group)
    
    selectInput("PCA_yaxis",
                label = "Y-Axis",
                choices = PCs(),
                selected = "PC2")
  })
  
  ### PRINT/TABLE/PLOTS --------------------------------------------------------
  output$PC_plot <- renderPlot({
    PC_plot()
  })
  
  output$PC_variance_plot <- renderPlot({
    PC_variance_plot()
  })
  
  output$PCA_download <- downloadHandler(
    filename = function() {
      paste("PCA_", type(), "_", Sys.Date(), ".", input$PCA_extension , sep = "")
    },
    content = function(file){
      ggplot2::ggsave(file, PC_plot(), device = input$PCA_extension, width = 9, height = 5)
    }
  )
  
  ## UMAP ######################################################################
  
  ### REACTIVE OBJECTS ---------------------------------------------------------
  UMAP_plot <- reactive({
    req(input$data_file, input$annotation_file)
    drawUMAP(eset(), 
             type = type(), 
             color = input$UMAP_color,
             title_add = input$UMAP_title,
             shapes = input$UMAP_shapes,
             add_labels = input$UMAP_labels)
  }) %>% bindEvent(input$UMAP_button)
  
  ### PRINT/TABLE/PLOTS --------------------------------------------------------
  output$UMAP_plot <- renderPlot({
    UMAP_plot()
  })
  
  output$UMAP_download <- downloadHandler(
    filename = function() {
      paste("UMAP_", type(), "_", Sys.Date(), ".", input$UMAP_extension, sep = "")
    },
    content = function(file){
      ggplot2::ggsave(file, UMAP_plot(), device = input$UMAP_extension, width = 9, height = 5)
    }
  )
  
  ## HEATMAP ###################################################################
  
  ### REACTIVE OBJECTS ---------------------------------------------------------
  heatmap_plot <- reactive({
    req(input$subset, input$data_file, input$annotation_file)
    
    drawHeatmaps(eset(), 
                 type = type(), 
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
  observeEvent(
    {input$subset}, {
    output$select_features <- renderUI({
      req((input$subset) == "subset_features")
      type <- type()
      annotation <- annotation()[,c("SampleName", type)]
      features <- annotation[stats::complete.cases(annotation),]
      features <- features$SampleName
      selectInput("features_list",
                  label = "Select Features to Include",
                  choices = features,
                  multiple = TRUE)
    })
  })
  
  observeEvent({input$data_file 
    input$subset},
    output$gene_list_select <- renderUI({
      req((input$subset) == "subset_genes")
      eset <- makeEset(data(), 
                       annotation(), 
                       type = type(), 
                       data_format = input$data_format,
                       log_transform = input$log_transform,
                       uniprot_annotation = FALSE)
      
      eset_genes <- sub("_.*", "", row.names(exprs(eset())))
      br()
      
      validate(need(try(eset_genes), "No gene names available."))
      
      selectizeInput("gene_list", 
                     label = "Select genes by typing their names into the text box or
                  by selecting them from the dropdown",
                     eset_genes,
                     multiple = TRUE
      )
    })
  )
  
  observeEvent({input$data_file
    input$subset},
    output$gene_list_file <- renderUI({
      req((input$subset) == "subset_genes")
      br()
      
      fileInput("gene_list_file",
                label = "OR Upload a tab-delimited file containing a column of gene names with no column header",
                accept = c("text",
                           "text/tab-separated-values,text/plain",
                           ".tsv")
      )
    })
  )
  
  observeEvent(
    {input$annotation_file 
     input$subset}, {
        output$select_variable <- renderUI({
          req((input$subset) == "subset_variable")
          numericInput("variable_list",
                       label = "Select number of genes with highest variance to include",
                       value = 50,
                       min = 1,
                       max = 1000,
                       step = 10)
        })
    }
  )
  
  ### PRINT/TABLE/PLOTS --------------------------------------------------------
  output$heatmap_plot <- renderPlot({
    heatmap_plot()
  }, height = 800)
  
  
  output$heatmap_download <- downloadHandler(
    filename = function() {
      paste("heatmap_", type(), "_", Sys.Date(), ".", input$heatmap_extension, sep = "")
    },
    
    content = function(file){
      if(input$heatmap_extension == "pdf"){
        
        grDevices::pdf(file, width = 10, height = 10)
        ComplexHeatmap::draw(heatmap_plot())
        grDevices::dev.off()
        
      } else if(input$heatmap_extension == "png"){
          grDevices::png(file, width = 1000, height = 1000)
          ComplexHeatmap::draw(heatmap_plot())
          grDevices::dev.off()
      }
    }
  )
  
  ## MD PLOT ###################################################################
  
  ### REACTIVE OBJECTS ---------------------------------------------------------
  md_plot <- reactive({
    drawMD(annot = annotation(),
           eset = eset(),
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
    req(input$data_file, input$annotation_file, (input$data_format != "none"), input$group)
    contrasts <- drawMD(annot = annotation(),
                        eset = eset(),
                        return_logfc_index = TRUE,
                        logfc_index_choice = NULL)
  })

  ### RENDER UI ----------------------------------------------------------------  
  output$md_contrast <- renderUI({
    req(input$annotation_file, input$data_file)
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
    
    genes <- sub("_.*", "", rownames(exprs(eset())))
    
    selectInput("md_label_specific_gene",
              label = "Type gene name",
              multiple = TRUE,
              choices = genes)
  }) %>% bindEvent(input$md_label_specific)
  
  ### PRINT/TABLE/PLOTS --------------------------------------------------------
  output$md_plot <- renderPlot({
    md_plot()
  }, height = 600)
  
  output$md_download <- downloadHandler(
    filename = function() {
      paste("mdplot_", type(), "_", Sys.Date(), ".", input$md_extension, sep = "")
    },
    content = function(file){
      ggsave(file, md_plot(), device = input$md_extension, width = 10)
    }
  )
  
  ## LIMMA LINEAR MODELING #####################################################
  
  ### REACTIVE OBJECTS ---------------------------------------------------------
  linear_model <- reactive({
    req(input$limma_button, input$data_file, input$annotation_file, (input$data_format != "none"), input$group)
    
    limmaLM(annot = annotation(),
            eset = eset(),
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
    
  }) %>% bindEvent(input$limma_button)
  
  coef_options <- reactive({
    req(input$limma_button, input$data_file, input$annotation_file, (input$data_format != "none"), input$group)
    
    coef <- limmaLM(annot = annotation(),
                    eset = eset(),
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
    
  })
  
  linear_model_table <- reactive({
    req(input$coef_options)
    
    coef_options <- limmaLM(annot = annotation(),
                            eset = eset(),
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
    
    makeTopTable(linear_model(),
                 coef = input$coef_options,
                 coef_options,
                 number_genes = input$number_genes,
                 logfc_th_add = input$logfc_threshold,
                 logfc_th = input$logfc_thresh)
    
  })
  
  linear_model_results_table <- reactive({
    req(input$limma_button, input$data_file, input$annotation_file)
    
    fit <- linear_model()
    
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
    
    table
    
  })
  
  annot_columns <- reactive ({
    validate(
      need(input$data_file, "Please input valid data file!"),
      need(input$annotation_file, "Please input valid annotation file!"))
    
    annot_columns <- colnames(annotation())
    
    annot_names <- data.frame(openxlsx::read.xlsx(input$annotation_file$datapath, 1, colNames = FALSE))
    annot_names <- annot_names[-(1:10),]
    annot_names <- annot_names[,-3]
    colnames(annot_names) <- c("File", "Name")
    annot_names <- na.omit(annot_names)
    annot_names <- as.list(annot_names$Name)
    
    annot_columns <- annot_columns[!annot_columns == "SampleName"]
    annot_columns <- annot_columns[!annot_columns == "Group"]
    annot_columns <- annot_columns[!annot_columns == annot_names]
  })
  
  ### RENDER UI ----------------------------------------------------------------
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
  
  output$linear_model_equation_header <- renderUI({
    req(input$limma_button)
    
    tagList(h5("Model Equation"))
  })
  
  output$linear_factors <- renderUI({
    
    req(input$data_file, input$annotation_file, (input$data_format != "none"), input$group)
    
    options <- unique(pData(eset())$Group)
    
    selectInput("linear_factors",
                label = "Choose at least two groups to include as factors in the model",
                choices = options,
                multiple = TRUE)
  })
  
  output$choose_contrasts <- renderUI({
    req((input$contrast_fit == TRUE), isTruthy(length(input$linear_factors) > 1))
    contrasts <- limmaLM(annot = annotation(),
                         eset = eset(),
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
      helpText("If you do not select any contrasts, the model will automatically generate with all possible contrasts.")
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
    req(input$limma_button, input$data_file, input$annotation_file)
    
    tagList(
      h6("Coefficient"),
      selectInput("coef_options",
                  label = "Choose the coefficient for logFC calculation and top table generation",
                  choices = coef_options(),
                  multiple = FALSE)
    )
  }) %>% bindEvent(input$limma_button)
  
  output$time_type <- renderUI({
    req(input$include_timeseries == TRUE)
    
    validate(
      need(input$data_file, "Please input valid data file!"),
      need(input$annotation_file, "Please input valid annotation file!"))
    
    radioButtons("time_type",
                 label = "Select how to include the time variable in the model",
                 choices = c("Continuous" = "continuous",
                             "Discrete" = "discrete"))
  }) %>% bindEvent(input$include_timeseries)
  
  output$time_points <- renderUI({
    req(isTruthy(input$include_timeseries))
    
    if(req(input$time_type) == "continuous"){
      number = length(unique(pData(eset())[,input$time_col]))
      numericInput("time_points_cont",
                   label = "Choose number of time points to include in the model",
                   min = 0,
                   max = number,
                   value = 3) 
      
    } else if(req(input$time_type) == "discrete"){
      options <- unique(pData(eset())[,input$time_col])
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
  
  output$number_genes <- renderUI({
    req(input$limma_button)
    
    num_rows <- nrow(eset())
    
    tagList(
      h6("Gene Number"),
      numericInput("number_genes",
                   label = "Number of top genes/rows to include in table",
                   value = num_rows,
                   min = 10,
                   max = NA,
                   step = 10
      ),
      hr()
    )
  })
  
  output$logfc_th <- renderUI({
    req(input$logfc_threshold == TRUE)
    
    numericInput("logfc_thresh",
                 label = "Choose logFC threshold",
                 step = 0.1,
                 value = 0.1)
  }) %>% bindEvent(input$logfc_threshold)
  
  ### PRINT/TABLE/PLOTS --------------------------------------------------------
  
  output$model_equation <- renderPrint({
    req(input$limma_button)
    limmaEquation(groups = input$group,
                  include_covariates = input$add_covariate,
                  covariate_col = input$covariate_col,
                  time_series = input$include_timeseries)
  }) %>% bindEvent(input$limma_button)
  
  output$linear_model_table <- DT::renderDataTable({
    
    table <- linear_model_table()
    rownames(table) <- 1:nrow(table)
    table <- table[order(-abs(table[['logFC']])),]
    
    DT::datatable(table, rownames = FALSE) %>% 
      DT::formatSignif(columns = 2:ncol(table), digits = 3) %>%
      DT::formatStyle('logFC', fontWeight = 'bold', color = styleInterval(0, c('blue', 'red'))) %>%
      DT::formatStyle('adj.P.Val', backgroundColor = styleInterval(0.05, c('darkseagreen', 'darksalmon')))
    
  }) %>% bindEvent(input$toptable_button)
  
  output$linear_model_results <- renderTable({
    linear_model_results_table()
    
  }, rownames = TRUE, digits = 0, hover = TRUE, bordered = TRUE, striped = TRUE, width = '100%', align = 'l') %>% bindEvent(input$limma_button)
  
  output$linear_table_download <- downloadHandler(
    filename = function() {
      paste("limma_topTable_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file){
      write.csv(linear_model_table(), file, row.names = FALSE)
    }
  )
  
  output$linear_results_table_download <- downloadHandler(
    filename = function() {
      paste("limma_results_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file){
      write.csv(linear_model_results_table(), file, row.names = TRUE)
    }
  )
  
  ## VOLCANO PLOT ##############################################################
  
  ### REACTIVE OBJECTS ---------------------------------------------------------
  volcano_linear_model_table <- reactive({
    req(input$volcano_coef_options)
    
    coef_options <- limmaLM(annot = annotation(),
                            eset = eset(),
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
    
    makeTopTable(linear_model(),
                 coef = input$volcano_coef_options,
                 coef_options,
                 number_genes = input$volcano_number_genes,
                 logfc_th_add = input$logfc_threshold,
                 logfc_th = input$logfc_thresh)
    
  })
  
  volcano_plot <- reactive({
    req(input$data_file, input$annotation_file)
    
    drawVolcano(dat = volcano_linear_model_table(), 
                type = type(),
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
  }) %>% bindEvent(input$volcano_button)
  
  ### RENDER UI ----------------------------------------------------------------
  output$volcano_coef_options <- renderUI({
    req(input$limma_button, input$data_file, input$annotation_file)
    
    selectInput("volcano_coef_options",
                label = "Choose the coefficient/contrast to use for the plot",
                choices = coef_options(),
                multiple = FALSE)
  }) %>% bindEvent(input$limma_button)
  
  output$volcano_number_genes <- renderUI({
    req(input$limma_button)
    
    num_rows <- nrow(eset())
    
    numericInput("volcano_number_genes",
                 label = "Number of top genes/rows to include in plot",
                 value = num_rows,
                 min = 10,
                 max = NA,
                 step = 10
    )
  })
  
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
                 label = "Select number of genes to label",
                 value = 20)
  }) %>% bindEvent(input$volcano_labels)
  
  ### PRINT/TABLE/PLOTS --------------------------------------------------------
  output$volcano_plot <- renderPlot({
    volcano_plot()
  }, height = 700, width = 1000)
  
  output$volcano_download <- downloadHandler(
    filename = function() {
      paste("volcanoplot_", type(), "_", Sys.Date(), ".", input$volcano_extension, sep = "")
    },
    content = function(file){
      ggsave(file, volcano_plot(), device = input$volcano_extension, width = 10, height = 8)
    }
  )
  
  ## INTERACTIVE VOLCANO PLOT ##################################################
  
  ### REACTIVE OBJECTS ---------------------------------------------------------
  int_volcano_linear_model_table <- reactive({
    req(input$int_volcano_coef_options)
    
    coef_options <- limmaLM(annot = annotation(),
                            eset = eset(),
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
    
    makeTopTable(linear_model(),
                 coef = input$int_volcano_coef_options,
                 coef_options,
                 number_genes = input$int_volcano_number_genes,
                 logfc_th_add = input$logfc_threshold,
                 logfc_th = input$logfc_thresh)
    
  })
  
  int_volcano_plot <- reactive({
    req(input$data_file, input$annotation_file)
    
    drawVolcano(dat = int_volcano_linear_model_table(), 
                type = type(),
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
  }) %>% bindEvent(input$int_volcano_button)
  
  glimma_args <- reactive({
    args <- interactiveVolcano(eset = eset(),
                               fit = linear_model(), 
                               top_table = int_volcano_linear_model_table(), 
                               type = type(), 
                               coef = input$int_volcano_coef_options)
    
  })
  
  ### RENDER UI ----------------------------------------------------------------
  output$int_volcano_coef_options <- renderUI({
    req(input$limma_button, input$data_file, input$annotation_file)
    
    selectInput("int_volcano_coef_options",
                label = "Choose the coefficient/contrast to use for the plot",
                choices = coef_options(),
                multiple = FALSE)
  }) %>% bindEvent(input$limma_button)
  
  output$int_volcano_number_genes <- renderUI({
    req(input$limma_button)
    
    num_rows <- nrow(eset())
    
    numericInput("int_volcano_number_genes",
                 label = "Number of top genes/rows to include in plot",
                 value = num_rows,
                 min = 10,
                 max = NA,
                 step = 10
    )
  })
  
  output$glimma <- renderUI({
    includeHTML(paste("glimma-plots/Volcano-Plot_", type(), ".html", sep = ""))
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
    
    coef_options <- limmaLM(annot = annotation(),
                            eset = eset(),
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
    
    makeTopTable(linear_model(),
                 coef = input$dif_heatmap_coef_options,
                 coef_options,
                 number_genes = nrow(eset()),
                 logfc_th_add = input$logfc_threshold,
                 logfc_th = input$logfc_thresh)
    
  })
  
  dif_heatmap_plot <- reactive({
    req(input$data_file, input$annotation_file)
    
    drawHeatmaps(eset(), 
                 type = type(), 
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
    req(input$limma_button, input$data_file, input$annotation_file)
    
    selectInput("dif_heatmap_coef_options",
                label = "Choose the coefficient/contrast to use for the plot",
                choices = coef_options(),
                multiple = FALSE)
  }) %>% bindEvent(input$limma_button)
  
  ### PRINT/TABLE/PLOTS --------------------------------------------------------
  output$dif_heatmap_plot <- renderPlot({
    dif_heatmap_plot()
  }, height = 800)
  
  output$dif_heatmap_download <- downloadHandler(
    filename = function() {
      paste("difheatmap_", type(), "_", Sys.Date(), ".", input$dif_heatmap_extension, sep = "")
    },
    
    content = function(file){
      if(input$dif_heatmap_extension == "pdf"){
        
        grDevices::pdf(file, height = 10, width = 10)
        ComplexHeatmap::draw(dif_heatmap_plot())
        grDevices::dev.off()
        
      } else if(input$dif_heatmap_extension == "png"){
        grDevices::png(file, height = 1000, width = 1000)
        ComplexHeatmap::draw(dif_heatmap_plot())
        grDevices::dev.off()
      }
    }
  )
  
  ## ENRICHMENT ################################################################
  
  ### REACTIVE OBJECTS ---------------------------------------------------------
  
  geneList <- reactive({
    calculateGeneList(fit = linear_model(),
                      coef = input$enrichment_coef)
  })
  
  enriched <- reactive({
    clusterProfilerEnrichment(geneList = geneList(),
                              gmt = input$gmt,
                              enrichment = input$enrichment_type,
                              pval_cutoff = 0.05)
  })

  enriched_table <- reactive({
    enrichedTable(enriched = enriched(),
                  enrichment = input$enrichment_type)
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
                  enrichment = input$enrichment_type,
                  gmt = input$gmt,
                  geneList = geneList(),
                  contrast = input$enrichment_coef,
                  plottype = "heat",
                  heatmap_color = input$enriched_heatmap_color,
                  heatmap_title = input$enriched_heatmap_title)
  })
  
  ### RENDER UI ----------------------------------------------------------------
  
  output$enrichment_coef_options <- renderUI({
    validate(need(input$limma_button, "Please generate a model in the data modeling tab."))
    req(input$limma_button, input$data_file, input$annotation_file)
    
    selectInput("enrichment_coef",
                label = "Choose the coefficient/contrast to use for enrichment",
                choices = coef_options(),
                multiple = FALSE)
  })
  
  output$gmts <- renderUI({
    req(input$species)
    
    if (input$species == "human"){
      selectInput("gmt",
                  label = "Choose database to use for enrichment",
                  choices = c("KEGG" = "HUMAN_KEGG",
                              "Reactome" = "HUMAN_MOUSE_Reactome",
                              "WikiPathways" = "HUMAN_WikiPathways",
                              "GO Biological Processes" = "HUMAN_GOBP",
                              "GO Cellular Components" = "HUMAN_GOCC",
                              "GO Molecular Functions" = "HUMAN_GOMF",
                              "GO All" = "HUMAN_GOALL",
                              "Hallmark MSigdb" = "HUMAN_Hallmark",
                              "MSigdb" = "HUMAN_MSigDB",
                              "All Available Databases" = "HUMAN_AllPathways"),
                  multiple = FALSE)
      
    } else if (input$species == "mouse"){
      selectInput("gmt",
                  label = "Choose database to use for enrichment",
                  choices = c("Reactome" = "HUMAN_MOUSE_Reactome",
                              "GO All" = "MOUSE_GOALL",
                              "GO Biological Processes" = "MOUSE_GOBP",
                              "GO Cellular Components" = "MOUSE_GOCC",
                              "GO Molecular Functions" = "MOUSE_GOMF",
                              "WikiPathways" = "MOUSE_WikiPathways",
                              "All Available Databases" = "MOUSE_AllPathways"),
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
                              "All Available Databases" = "RAT_AllPathways"),
                  multiple = FALSE)
      
    } else if (input$species == "celegans"){
      selectInput("gmt",
                  label = "Choose database to use for enrichment",
                  choices = c("KEGG" = "CELEGANS_KEGG",
                              "WikiPathways" = "CELEGANS_WikiPathways",
                              "GO Biological Processes" = "CELEGANS_GOBP",
                              "GO Cellular Components" = "CELEGANS_GOCC",
                              "GO Molecular Functions" = "CELEGANS_GOMF",
                              "All Available Databases" = "CELEGANS_AllPathways"),
                  multiple = FALSE)
      
    } else if (input$species == "fruitfly"){
      selectInput("gmt",
                  label = "Choose database to use for enrichment",
                  choices = c("KEGG" = "FRUITFLY_KEGG",
                              "WikiPathways" = "FRUITFLY_WikiPathways",
                              "GO Biological Processes" = "FRUITFLY_GOBP",
                              "GO Cellular Components" = "FRUITFLY_GOCC",
                              "GO Molecular Functions" = "FRUITFLY_GOMF",
                              "All Available Databases" = "FRUITFLY_AllPathways"),
                  multiple = FALSE)
      
    } else if (input$species == "zebrafish"){
      selectInput("gmt",
                  label = "Choose database to use for enrichment",
                  choices = c("KEGG" = "ZEBRAFISH_KEGG",
                              "WikiPathways" = "ZEBRAFISH_WikiPathways",
                              "GO Biological Processes" = "ZEBRAFISH_GOBP",
                              "GO Cellular Components" = "ZEBRAFISH_GOCC",
                              "GO Molecular Functions" = "ZEBRAFISH_GOMF",
                              "All Available Databases" = "ZEBRAFISH_AllPathways"),
                  multiple = FALSE)
      
    } else if (input$species == "yeast"){
      selectInput("gmt",
                  label = "Choose database to use for enrichment",
                  choices = c("KEGG" = "YEAST_KEGG",
                              "WikiPathways" = "YEAST_WikiPathways",
                              "GO Biological Processes" = "YEAST_GOBP",
                              "GO Cellular Components" = "YEAST_GOCC",
                              "GO Molecular Functions" = "YEAST_GOMF",
                              "All Available Databases" = "YEAST_AllPathways"),
                  multiple = FALSE)
      
    }

  })
  
  ### PRINT/TABLE/PLOTS --------------------------------------------------------
  
  output$enriched_table <- DT::renderDataTable({
    req(input$enrichment_button)
    
    if(input$enrichment_type == "GSEA"){
      DT::datatable(enriched_table()[,-7], rownames = FALSE, options = list(order = list(1, 'asc'), pageLength = 15)) %>% 
        DT::formatRound(columns = c(6), digits = 3) %>%
        DT::formatSignif(columns = c(2,3,4), digits = 3) %>%
        DT::formatStyle('setSize', background = styleColorBar(enriched_table()$setSize, "lightskyblue")) %>%
        DT::formatStyle('NES', fontWeight = 'bold', color = styleInterval(0, c('blue', 'red')))
      
    } else{
      DT::datatable(enriched_table()[,-7], rownames = FALSE, options = list(order = list(1, 'asc'), pageLength = 15)) %>% 
        DT::formatSignif(columns = c(2,3,4), digits = 3)%>%
        DT::formatStyle('Count', background = styleColorBar(enriched_table()$Count, "lightskyblue"))
    }
  }) %>% bindEvent(input$enrichment_button)
  
  output$enriched_table_download <- downloadHandler(
    filename = function() {
      paste("enrichment_", type(), "_", Sys.Date(), ".csv", sep = "")
    },
    
    content = function(file){
      write.csv(enriched_table(), file, row.names = FALSE)
    }
  )
  
  output$enriched_dotplot <- renderPlot({
    req(input$enrichment_button)
    enriched_dotplot()
  }, height = 800) %>% bindEvent(input$dotplot_button, input$enrichment_button)
  
  output$dotplot_download <- downloadHandler(
    filename = function() {
      paste("enriched_dotplot_", type(), "_", Sys.Date(), ".pdf", sep = "")
    },
    
    content = function(file){
      grDevices::pdf(file, height = 8, width = 10)
      plot(enriched_dotplot())
      grDevices::dev.off()
    }
  )
  
  output$enriched_cnetplot <- renderPlot({
    req(input$enrichment_button)
    enriched_cnetplot()
  }, height = 800) %>% bindEvent(input$cnetplot_button, input$enrichment_button)
  
  output$cnetplot_download <- downloadHandler(
    filename = function() {
      paste("enriched_cnetplot_", type(), "_", Sys.Date(), ".pdf", sep = "")
    },
    
    content = function(file){
      grDevices::pdf(file, height = 10, width = 12)
      plot(enriched_cnetplot())
      grDevices::dev.off()
    }
  )
  
  output$enriched_emapplot <- renderPlot({
    req(input$enrichment_button)
    enriched_emapplot()
  }, height = 800) %>% bindEvent(input$emap_button, input$enrichment_button)
  
  output$emap_download <- downloadHandler(
    filename = function() {
      paste("enriched_emap_", type(), "_", Sys.Date(), ".pdf", sep = "")
    },
    
    content = function(file){
      grDevices::pdf(file, height = 10, width = 12)
      plot(enriched_emapplot())
      grDevices::dev.off()
    }
  )
  
  output$enriched_upsetplot <- renderPlot({
    req(input$enrichment_button)
    enriched_upsetplot()
  }, height = 800) %>% bindEvent(input$upset_button, input$enrichment_button)
  
  output$upset_download <- downloadHandler(
    filename = function() {
      paste("enriched_upsetplot_", type(), "_", Sys.Date(), ".pdf", sep = "")
    },
    
    content = function(file){
      grDevices::pdf(file, height = 10, width = 12)
      plot(enriched_upsetplot())
      grDevices::dev.off()
    }
  )
  
  output$enriched_treeplot <- renderPlot({
    req(input$enrichment_button)
    enriched_treeplot()
  }, height = 800) %>% bindEvent(input$tree_button, input$enrichment_button)
  
  output$tree_download <- downloadHandler(
    filename = function() {
      paste("enriched_treeplot_", type(), "_", Sys.Date(), ".pdf", sep = "")
    },
    
    content = function(file){
      grDevices::pdf(file, height = 10, width = 12)
      plot(enriched_treeplot())
      grDevices::dev.off()
    }
  )
  
  output$enriched_heatplot <- renderPlotly({
    req(input$enrichment_button)
    enriched_heatplot()
  }) %>% bindEvent(input$enriched_heatmap_button, input$enrichment_button)
  
  ## REPORT GENERATION #########################################################
  
  output$excel_summary_download_button <- downloadHandler(
    filename = function() {
      paste(gsub("\\.", "", make.names(type())), "_Summary_", Sys.Date(), ".xlsx", sep = "")
    },
    content = function(file){
      wbOut <- openxlsx::createWorkbook()
      
      if(type() != "met_combined"){
        eSet = eset()[,order(pData(eset())$Group)]
        writeDataToSheets(wb = wbOut, 
                          eset = eSet,
                          limmaFit = report_linear_model(),
                          type = type(), 
                          data_format = input$data_format);
      }
      
      openxlsx::saveWorkbook(wbOut, file = file, overwrite = TRUE)
    }
  )
  
  output$checkrender <- renderText({
    if (identical(rmarkdown::metadata$runtime, "shiny")) {
      TRUE
    } else {
      FALSE
    }
  })
  
  output$report <- downloadHandler(
    
    filename = paste("Report_", Sys.Date(), ".html"),
    
    content = function(file) {
      withProgress(message = 'Rendering the HTML report:', detail = "Running Enrichment...",
                   
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      {tempReport <- file.path(tempdir(), "report.Rmd")
      file.copy("report.Rmd", tempReport, overwrite = TRUE)
      
      # Set up parameters to pass to Rmd document
      params <- list(eset = eset(),
                     eset_prenorm = eset_prenorm(),
                     annot = annotation(),
                     type = type(),
                     log_transform = input$log_transform,
                     data_format = input$data_format,
                     include_norm = input$report_norm,
                     norm = input$norm_eset,
                     norm_plots = input$report_norm_plot,
                     include_PCA = input$report_PCA,
                     PCA_densities = input$report_PCA_densities,
                     PCA_ellipse = input$report_PCA_ellipse,
                     PCA_label = input$report_PCA_labels,
                     PCA_shape = input$report_PCA_shapes,
                     include_UMAP = input$report_UMAP,
                     UMAP_shape = input$report_UMAP_shapes,
                     UMAP_label = input$report_UMAP_labels,
                     include_MD = input$report_MD,
                     MD_label = input$report_md_add_labels,
                     MD_cutoff = input$report_md_fc_cutoff,
                     MD_contrasts = md_contrasts(),
                     include_corr = input$report_corr,
                     corr_shape = input$report_corr_shape,
                     corr_type = input$report_corr_type,
                     corr_numbers = input$report_corr_number,
                     include_heatmap = input$report_heatmap,
                     heatmap_kcluster = input$report_kclustering,
                     heatmap_rowcluster = input$report_cluster_samples,
                     heatmap_zscore = input$report_zscore,
                     include_volcano = input$report_volcano,
                     volcano_label = input$report_volcano_labels,
                     volcano_yaxis = input$report_volcano_yaxis,
                     volcano_FCcutoff = input$report_volcano_top_fc,
                     volcano_Pcutoff = input$report_volcano_top_values,
                     linear_model = report_linear_model(),
                     contrasts = report_coef_options(),
                     model_equation = report_model_equation(),
                     include_dif_heatmap = input$report_dif_heatmap,
                     dif_heatmap_kcluster = input$report_dif_kclustering,
                     dif_heatmap_rowcluster = input$report_dif_cluster_samples,
                     dif_heatmap_zscore = input$report_dif_zscore,
                     rendered_by_shiny = TRUE,
                     species = input$report_species,
                     enrichment_database = input$report_gmt,
                     enrichment = input$report_enrichment_type,
                     enrichplots = input$report_enrichplots,
                     enriched = report_enricheds(),
                     geneList = report_geneLists())
      
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv()))}
      )
    }
  )
  
  ### REPORT NORM --------------------------------------------------------------
  output$report_norm_options <- renderUI({
    req(input$report_norm == TRUE)
    
    tagList(
      checkboxGroupInput("report_norm_plot",
                    label = "Choose which plots to include",
                    choices = c("Boxplot",
                                "Density",
                                "MA",
                                "RLE"),
                    selected = "Density")
    )
  }) %>% bindEvent(input$report_norm)
  
  ### REPORT PCA ---------------------------------------------------------------
  output$report_PCA_options <- renderUI({
    req(input$report_PCA == TRUE)
    
    tagList(
      checkboxInput("report_PCA_densities",
                    label = "Include PC density curves along axes"),
      
      checkboxInput("report_PCA_ellipse",
                    label = "Add ellipses"),
      
      checkboxInput("report_PCA_labels",
                    label = "Label points",
                    value = TRUE),
      
      checkboxInput("report_PCA_shapes",
                    label = "Make different groups different shapes")
    )
  }) %>% bindEvent(input$report_PCA)
  
  ### REPORT UMAP --------------------------------------------------------------
  output$report_UMAP_options <- renderUI({
    req(input$report_UMAP == TRUE)
    
    tagList(
      checkboxInput("report_UMAP_shapes",
                    label = "Make different groups different shapes"),
      
      checkboxInput("report_UMAP_labels",
                    label = "Label points",
                    value = TRUE)
    )
  }) %>% bindEvent(input$report_UMAP)
  
  ### REPORT MD PLOT -----------------------------------------------------------
  output$report_MD_options <- renderUI({
    req(input$report_MD == TRUE)
    
    tagList(
      checkboxInput("report_md_add_labels",
                    "Add labels to genes with highest/lowest logFC",
                    value = TRUE),
      
      numericInput("report_md_fc_cutoff",
                   label = "LogFC Cutoff",
                   value = 1,
                   step = .5)
    )
  }) %>% bindEvent(input$report_MD)
  
  ### REPORT CORRELATION -------------------------------------------------------
  output$report_corr_options <- renderUI({
    req(input$report_corr == TRUE)
    
    tagList(
      radioButtons("report_corr_shape",
                   label = "Shape of points",
                   choices = c("Square" = "square",
                               "Circle" = "circle")),
      
      radioButtons("report_corr_type",
                   label = "Select graph orientation",
                   choices = c("Full" = "full",
                               "Lower" = "lower",
                               "Upper" = "upper"),
                   selected = "lower"),
      
      checkboxInput("report_corr_number",
                    label = "Add correlation values",
                    value = TRUE)
    )
  }) %>% bindEvent(input$report_corr)
  
  ### REPORT HEATMAP -----------------------------------------------------------
  output$report_heatmap_options <- renderUI({
    req(input$report_heatmap == TRUE)
    
    tagList(
      checkboxInput("report_kclustering",
                    "K-cluster rows"),
      
      checkboxInput("report_cluster_samples",
                    label = "Cluster columns",
                    value = TRUE),
      
      checkboxInput("report_zscore",
                    label = "Z-score across rows",
                    value = TRUE)
    )
  }) %>% bindEvent(input$report_heatmap)
  
  ### REPORT LIMMA LINEAR MODEL ------------------------------------------------
  output$report_linear_factors <- renderUI({
    req(input$annotation_file, input$data_file)
    options <- unique(pData(eset())$Group)
    
    selectInput("report_linear_factors",
                label = "Choose at least two groups to include as factors in the model",
                choices = options,
                multiple = TRUE)
  })
  
  report_linear_model <- reactive({
    req(input$data_file, input$annotation_file)
    
    limmaLM(annot = annotation(),
            eset = eset(),
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
    
  })
  
  output$report_covariate_col <- renderUI({
    req(input$report_add_covariate == TRUE)
    
    annot_columns <- colnames(annotation())
    
    annot_names <- data.frame(openxlsx::read.xlsx(input$annotation_file$datapath, 1, colNames = FALSE))
    annot_names <- annot_names[-(1:10),]
    annot_names <- annot_names[,-3]
    colnames(annot_names) <- c("File", "Name")
    annot_names <- na.omit(annot_names)
    annot_names <- as.list(annot_names$Name)
    
    annot_columns <- annot_columns[!annot_columns == "SampleName"]
    annot_columns <- annot_columns[!annot_columns == "Group"]
    annot_columns <- annot_columns[!annot_columns == annot_names]
    
    selectizeInput("report_covariate_col",
                   label = "Select up to 3 columns to act as covariates",
                   choices = annot_columns,
                   multiple = TRUE,
                   options = list(maxItems = 3))
    
  }) %>% bindEvent(input$report_add_covariate)
  
  output$report_choose_contrasts <- renderUI({
    req((input$report_contrast_fit == TRUE), isTruthy(length(input$report_linear_factors) > 1))
    contrasts <- limmaLM(annot = annotation(),
                         eset = eset(),
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
  
  report_coef_options <- reactive({
    req(input$data_file, input$annotation_file)
    
    coef <- limmaLM(annot = annotation(),
                    eset = eset(),
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
    
  })
  
  report_model_equation <- reactive({
    limmaEquation(groups = input$group,
                  include_covariates = input$report_add_covariate,
                  covariate_col = input$report_covariate_col,
                  time_series = input$report_include_timeseries)
  })
  
  output$report_choose_time <- renderUI({
    req(input$report_include_timeseries == TRUE)
    
    selectInput("report_time_col",
                label = "Select column containting time series information",
                choices = annot_columns(),
                multiple = FALSE)
  })
  
  output$report_time_type <- renderUI({
    req(input$report_include_timeseries == TRUE)
    
    validate(
      need(input$data_file, "Please input valid data file!"),
      need(input$annotation_file, "Please input valid annotation file!"))
    
    radioButtons("report_time_type",
                 label = "Select how to include the time variable in the model",
                 choices = c("Continuous" = "continuous",
                             "Discrete" = "discrete"))
  }) %>% bindEvent(input$report_include_timeseries)
  
  output$report_time_points <- renderUI({
    req(isTruthy(input$report_include_timeseries))
    
    if(req(input$report_time_type) == "continuous"){
      number = length(unique(pData(eset())[,input$report_time_col]))
      numericInput("report_time_points_cont",
                   label = "Choose number of time points to include in the model",
                   min = 0,
                   max = number,
                   value = 3) 
      
    } else if(req(input$report_time_type) == "discrete"){
      options <- unique(pData(eset())[,input$report_time_col])
      options <- options[-order(options)[1]]
      options <- sort(options)
      selectInput("report_time_points_disc",
                  label = "Choose which time points to include in the model",
                  choices = options,
                  multiple = TRUE)
      
    }
    
  }) %>% bindEvent(input$report_time_type, input$report_include_timeseries, input$report_time_col)
  
  ### REPORT VOLCANO -----------------------------------------------------------
  output$report_volcano_options <- renderUI({
    req(input$report_volcano == TRUE)
    
    tagList(
      checkboxInput("report_volcano_labels",
                    label = "Add labels to genes with highest/lowest logFC",
                    value = TRUE),
      
      radioButtons("report_volcano_yaxis",
                   label = "Choose Y-Axis Variable",
                   choices = c("Adjusted P-Value" = "adj.P.Val",
                               "P-Value" = "P.Value")),
      
      numericInput("report_volcano_top_fc",
                   label = "LogFC Cut Off",
                   value = 1.0,
                   step = 0.5,
                   min = 0),
      
      numericInput("report_volcano_top_values",
                   label = "P-Value Cutoff",
                   value = 0.05,
                   step = 0.01,
                   min = 0,
                   max = 1)
    )
  }) %>% bindEvent(input$report_volcano)
  
  ### REPORT DIFFERENTIAL HEATMAP ----------------------------------------------
  output$report_dif_heatmap_options <- renderUI({
    req(input$report_dif_heatmap == TRUE)
    
    tagList(
      checkboxInput("report_dif_kclustering",
                    "K-cluster rows"),
      
      checkboxInput("report_dif_cluster_samples",
                    label = "Cluster columns",
                    value = TRUE),
      
      checkboxInput("report_dif_zscore",
                    label = "Z-score across rows",
                    value = TRUE)
    )
  }) %>% bindEvent(input$report_dif_heatmap)
  
  ### REPORT ENRICHMENT --------------------------------------------------------
  
  output$report_gmts <- renderUI({
    req(input$report_species)
    
    if (input$report_species == "human"){
      selectInput("report_gmt",
                  label = "Choose database to use for enrichment",
                  choices = c("KEGG" = "HUMAN_KEGG",
                              "Reactome" = "HUMAN_MOUSE_Reactome",
                              "WikiPathways" = "HUMAN_WikiPathways",
                              "GO Biological Processes" = "HUMAN_GOBP",
                              "GO Cellular Components" = "HUMAN_GOCC",
                              "GO Molecular Functions" = "HUMAN_GOMF",
                              "GO All" = "HUMAN_GOALL",
                              "Hallmark MSigdb" = "HUMAN_Hallmark",
                              "MSigdb" = "HUMAN_MSigDB",
                              "All Available Databases" = "HUMAN_AllPathways"),
                  multiple = FALSE)
      
    } else if (input$report_species == "mouse"){
      selectInput("report_gmt",
                  label = "Choose database to use for enrichment",
                  choices = c("Reactome" = "HUMAN_MOUSE_Reactome",
                              "GO All" = "MOUSE_GOALL",
                              "GO Biological Processes" = "MOUSE_GOBP",
                              "GO Cellular Components" = "MOUSE_GOCC",
                              "GO Molecular Functions" = "MOUSE_GOMF",
                              "WikiPathways" = "MOUSE_WikiPathways",
                              "All Available Databases" = "MOUSE_AllPathways"),
                  multiple = FALSE)
      
    } else if (input$report_species == "rat"){
      selectInput("report_gmt",
                  label = "Choose database to use for enrichment",
                  choices = c("Reactome" = "RAT_Reactome",
                              "GO All" = "RAT_GOALL",
                              "GO Biological Processes" = "RAT_GOBP",
                              "GO Cellular Components" = "RAT_GOCC",
                              "GO Molecular Functions" = "RAT_GOMF",
                              "WikiPathways" = "RAT_WikiPathways",
                              "All Available Databases" = "RAT_AllPathways"),
                  multiple = FALSE)
      
    } else if (input$report_species == "celegans"){
      selectInput("report_gmt",
                  label = "Choose database to use for enrichment",
                  choices = c("KEGG" = "CELEGANS_KEGG",
                              "WikiPathways" = "CELEGANS_WikiPathways",
                              "GO Biological Processes" = "CELEGANS_GOBP",
                              "GO Cellular Components" = "CELEGANS_GOCC",
                              "GO Molecular Functions" = "CELEGANS_GOMF",
                              "All Available Databases" = "CELEGANS_AllPathways"),
                  multiple = FALSE)
      
    } else if (input$report_species == "fruitfly"){
      selectInput("report_gmt",
                  label = "Choose database to use for enrichment",
                  choices = c("KEGG" = "FRUITFLY_KEGG",
                              "WikiPathways" = "FRUITFLY_WikiPathways",
                              "GO Biological Processes" = "FRUITFLY_GOBP",
                              "GO Cellular Components" = "FRUITFLY_GOCC",
                              "GO Molecular Functions" = "FRUITFLY_GOMF",
                              "All Available Databases" = "FRUITFLY_AllPathways"),
                  multiple = FALSE)
      
    } else if (input$report_species == "zebrafish"){
      selectInput("report_gmt",
                  label = "Choose database to use for enrichment",
                  choices = c("KEGG" = "ZEBRAFISH_KEGG",
                              "WikiPathways" = "ZEBRAFISH_WikiPathways",
                              "GO Biological Processes" = "ZEBRAFISH_GOBP",
                              "GO Cellular Components" = "ZEBRAFISH_GOCC",
                              "GO Molecular Functions" = "ZEBRAFISH_GOMF",
                              "All Available Databases" = "ZEBRAFISH_AllPathways"),
                  multiple = FALSE)
      
    } else if (input$report_species == "yeast"){
      selectInput("report_gmt",
                  label = "Choose database to use for enrichment",
                  choices = c("KEGG" = "YEAST_KEGG",
                              "WikiPathways" = "YEAST_WikiPathways",
                              "GO Biological Processes" = "YEAST_GOBP",
                              "GO Cellular Components" = "YEAST_GOCC",
                              "GO Molecular Functions" = "YEAST_GOMF",
                              "All Available Databases" = "YEAST_AllPathways"),
                  multiple = FALSE)
      
    }
  })
  
  report_geneLists <- reactive({
    geneList <- vector("list", length(report_coef_options()))

    for (i in 1:length(report_coef_options())){
      genes <- calculateGeneList(fit = report_linear_model(),
                                 coef = report_coef_options()[i])
      geneList[[i]] <- genes 
    }
    
    return(geneList)
  })
  
  report_enricheds <- reactive({
    enriched <- c()
    
    for (i in 1:length(report_geneLists())){
      enrich <- clusterProfilerEnrichment(geneList = report_geneLists()[[i]],
                                          gmt = input$report_gmt,
                                          enrichment = input$report_enrichment_type)
      enriched[[i]] <- enrich
    }
    
    return(enriched)
  })
  
  ## CLOSE OUT #################################################################
  
  session$onSessionEnded(function() {
    stopApp()
  })
}
