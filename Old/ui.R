
# NCmisc::list.functions.in.file("ui.R", alphabetic = TRUE)

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
source("enrichment.R")

ui <- fluidPage(
  
  shinyjs::useShinyjs(),
  
  tags$head(
    
    tags$style(".progress-bar{background-color:#61A6F9;}
               .tabbable > .nav > li > a {margin-bottom:50px;}",
               ".shiny-notification {
               height: 100px;
               width: 275px;
               position:fixed;
               top: calc(50% - 50px);;
               left: calc(50% - 125px);;}"),
  ),
  
  
  navbarPage("Omics Notebook",
    theme = bslib::bs_theme(bootswatch = "lux"),

    ## ABOUT ###################################################################
    navbarMenu("About",
      
      #### OVERVIEW ------------------------------------------------------------
      tabPanel("Overview",
        br(),
        fluidRow(
          column(2),
          column(8,
            includeMarkdown("Overview.md")
          ),
          column(2),
        ),
        fluidRow(
          column(12, align = "center",
          br(),
          br(),
          hr(),
          div(
            class = "footer",
            includeHTML("footer.Rhtml")
          )
          )
        )
      ),
      
      #### FILE INPUTS ---------------------------------------------------------
      tabPanel(title = "File Inputs", value = "fileinputspage",
        br(),
        fluidRow(
          column(2),
          column(8,
            includeMarkdown("fileInputs.md")
          ),
          column(2)
        ),
        fluidRow(
          column(12, align = "center",
                 br(),
                 br(),
                 hr(),
                 div(
                   class = "footer",
                   includeHTML("footer.Rhtml")
                 )
          )
        )
      )
    ),
    
    ## DATA ####################################################################
    tabPanel("1. Data",
      tabsetPanel(
        
        #### DATA UPLOAD -------------------------------------------------------
        tabPanel("Data Upload",
          sidebarLayout(
            sidebarPanel(width = 4,
              
              h4("Upload Data"),
              
              h6("The selections made in this tab will be used in all other analyses in Omics Notebook."),
              
              fileInput("data_file",
                        label = "Select Data File", 
                        buttonLabel = "Browse", 
                        multiple = FALSE,
                        accept = c("text",
                                   "text/tab-separated-values,text/plain",
                                   ".tsv",
                                   ".csv",
                                   "text/csv")),
                                                
              fileInput("annotation_file", 
                        label = "Select Annotation File", 
                        buttonLabel = "Browse", 
                        multiple = FALSE,
                        accept = c("application/vnd.ms-excel",
                                   "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                                   ".xlsx",
                                   ".xls")),
                                                
              helpText("For more information on formatting annotation files and accepted 
                 file types see the 'About' tab."),
              
              br(),
              br(),
              
              tags$div(
                tags$a(href = "https://drive.google.com/drive/folders/1lyzmIhorrZy_CKuxabi1Bv1cLHIblJhk?usp=sharing", 
                       "Click here"),
                " to download example data and annotation files. These files can be uploaded above to test Omics Notebook functionality."
              ),
              
              br(),
              
              uiOutput("typeHeader"),
              
              verbatimTextOutput("type"),
              
              uiOutput("select_group"),
              
              selectInput("data_format", 
                          label = "Select data format", 
                          choices = c("Choose One" = "none",
                                      "Protein Groups" = "Protein.Groups..MQ.", 
                                      "Peptides" = "Peptides..MQ.", 
                                      "Phosphosites" = "Sites..MQ.", 
                                      "Generic", 
                                      "Metabolites")),
                                                
              hr(),
              
              h5("Data Processing Options"),
              
              checkboxInput("log_transform",
                            label = "Log-transform",
                            value = TRUE),
              
              checkboxInput("uniprot_annotation",
                            label = "Uniprot annotation",
                            value = FALSE),
              
              selectInput("norm_eset",
                          label = "Normalization method",
                          choices = c("Quantile" = "quantile",
                                      "Median" = "median",
                                      "Loess" = "loess",
                                      "Z-Transform" = "z transform",
                                      "Median Absolute Deviation (MAD)" = "mad",
                                      "None" = "none")),
              
              numericInput("zero_cutoff",
                           label = "Choose zero cutoff",
                           value = 0.3,
                           min = 0,
                           max = 1),
              
              helpText("The cutoff value specifies the allowed proportion of zeros in a row before it is removed.
                       For example, a cutoff of 0.3 indicates that rows where more than 30% of values are zero will be removed.")
              
            ),
              
            mainPanel(
              fluidRow(
                column(8,
                       br(),
                       h3("View Data"),
                       
                       helpText("Summaries of the data that you upload will be available below."),
                       
                       br(),
                       br(),
                       
                       helpText("The data upload checks to the right can be used to verify that
                                your data has uploaded correctly and all necessary selections have been
                                made before moving on to further analysis."),
                       
                       br(),
                       br(),
                       
                       helpText("Another good way to check that your data inputs are being processed correctly is 
                                by loading the expression matrix in the ", code("View Expression Matrix"), " tab.")
                       ),
                column(4,
                       br(),
                       h5("Data Upload Checks"),
                       verbatimTextOutput("datainput_check"),
                       verbatimTextOutput("annotationinput_check"),
                       verbatimTextOutput("annotationformat_check"),
                       verbatimTextOutput("groupinput_check"),
                       verbatimTextOutput("dataformat_check"),
                       )
              ),
              
              hr(),
              
              uiOutput("filesHeader"),
               
              tableOutput("file_table"),
              tableOutput("annotation_table"),
               
              uiOutput("dataOverview"),
              
              uiOutput("dataOverview_help"),
              
              br(),
               
              uiOutput("annotationFile"),
              verbatimTextOutput("summary_annotation"),
               
              uiOutput("dataFile"),
              verbatimTextOutput("summary_data")
               
            )
          ),
          fluidRow(
            column(12, align = "center",
                   br(),
                   br(),
                   hr(),
                   div(
                     class = "footer",
                     includeHTML("footer.Rhtml")
                   )
            )
          )
        ),
        
        #### VIEW EXPRESSION MATRIX --------------------------------------------
        tabPanel("View Expression Matrix",
          mainPanel(
            fluidRow(
              column(1,
              ),
              
              column(10,
                br(),
                h3("Expression Matrix"),
                
                helpText("The table produced is the result of calling the ", 
                         code("exprs()"), 
                         " function from Biobase on the expression set object derived 
                         from the data. If a normalization method was chosen, it has been applied to
                         this expression set. LogFC is calculated and appended to this table."),
                
                br(),
                
                helpText("Note: If you have opted to include Uniprot annotation, this may 
                         take a couple of minutes to load."),
                
                hr(),
                
                verbatimTextOutput('pdata_eset'),
                
                shinycssloaders::withSpinner(DT::dataTableOutput("view_eset"), type = 8),
                
                br(),
                
                downloadButton("eset_table_download",
                               label = "Download Dataset",
                               style="color: #61A6F9; background-color: #D8EAFF; border-color: 
                          #D8EAFF; border-radius: 10px; border-width: 2px")
              ),
              
              column(1,
              )
            )
          ),
          
          fluidRow(
            column(12, align = "center",
                   br(),
                   br(),
                   hr(),
                   div(
                     class = "footer",
                     includeHTML("footer.Rhtml")
                   )
            )
          )
        ),
        
        #### NORMAZLIATION PLOTS -----------------------------------------------
        tabPanel("Normalization Plots",
          sidebarLayout(
            sidebarPanel(
              h4("Plot Options"),
              radioButtons("norm_plottype",
                           label = "Plot Type",
                           choices = c("Pre/Post Normalization Boxplot" = "Boxplot",
                                       "Pre/Post Normalization Density Plot" = "Density",
                                       "MA Plot" = "MA",
                                       "Relative Log Expression Plot" = "RLE")),
              
              uiOutput("ma_array_choose"),
              
              helpText("An MA plot visualizes the difference between one sample and others using the average log expression (M),
                       and the expression log ratio (A). The RLE plot helps visualize variation, and helps to show whether the
                       normalization procedure was successful."),
              
              br(),
              br(),
              
              actionButton("normal_button", "Generate Plot",
                           class = "btn-block",
                           style="color: #fff; background-color: #9BD79A; border-color: 
                          #9BD79A; border-radius: 10px; border-width: 2px"),
              
              br(),
              br(),
              
              h6("This analysis may take a few minutes to load. \n"),
              
              hr(),
              
              h5("Download Plot"),
              radioButtons("normal_extension", 
                           label = "Save As:",
                           choices = c("PDF" = "pdf", 
                                       "PNG" = "png"),
                           inline = TRUE),
              
              downloadButton("normal_download", 
                             label = "Download Plot",
                             style="color: #61A6F9; background-color: #D8EAFF; border-color: 
                          #D8EAFF; border-radius: 10px; border-width: 2px"),
            ),
            
            mainPanel(
              br(),
              h3("Normalization Plots"),
              hr(),
              
              shinycssloaders::withSpinner(plotOutput("normal_plot"), type = 8)
            )
          ),
          
          fluidRow(
            column(12, align = "center",
                   br(),
                   br(),
                   hr(),
                   div(
                     class = "footer",
                     includeHTML("footer.Rhtml")
                   )
            )
          )
        )
      )
    ),
    
    ## BASIC ANALYSIS ##########################################################
    tabPanel("2. Basic Analysis",
      tabsetPanel(
        
        #### PCA ---------------------------------------------------------------
        tabPanel("PCA",
          sidebarLayout(
            sidebarPanel(
              h4("Analysis Options"),
              selectInput("PCA_color", 
                          label = "Select Color Palette",
                          c("Set1",
                            "Set2",
                            "Set3",
                            "Pastel1",
                            "Pastel2",
                            "Paired",
                            "Dark2",
                            "Accent")),
              
              helpText("If there are more groups than the number of colors in the 
                       largest palette (12), the graph will default to a larger color palette with distinct hues."),
              
              br(),
              br(),
              
              textInput("PCA_title",
                        label = "Customize Plot Title",
                        placeholder = "Type title here"),
              
              uiOutput("PC_xaxis"),
              uiOutput("PC_yaxis"),
               
              checkboxInput("include_densities",
                            label = "Include PC density curves along axes"),
              
              checkboxInput("PCA_ellipse",
                            label = "Add ellipses"),
              
              checkboxInput("PCA_labels",
                            label = "Label points",
                            value = TRUE),
              
              checkboxInput("PCA_shapes",
                            label = "Make different groups different shapes"),
               
              hr(),
              
              actionButton("PCA_button", "Run Analysis",
                          class = "btn-block",
                          style="color: #fff; background-color: #9BD79A; border-color: 
                          #9BD79A; border-radius: 10px; border-width: 2px"),
               
              br(),
              br(),
               
              h6("This analysis may take a few minutes to load. \n"),
               
              hr(),
               
              h5("Download Plot"),
              radioButtons("PCA_extension", 
                          label = "Save As:",
                          choices = c("PDF" = "pdf", 
                                      "PNG" = "png"),
                          inline = TRUE),
               
              downloadButton("PCA_download", 
                            label = "Download PCA Plot",
                            style="color: #61A6F9; background-color: #D8EAFF; border-color: 
                          #D8EAFF; border-radius: 10px; border-width: 2px"),
            ),
                                   
            mainPanel(
              br(),
              h3("Principle Component Analysis"),
               
              helpText("PCA requires both ", em("Annotation "), "and ", em("Data "), 
                      "file inputs. It also requires selection of ", em("Data Type "), 
                      "in the Data Upload tab."),
              
              shinycssloaders::withSpinner(plotOutput("PC_variance_plot"), type = 8),
              
              helpText("The above plot shows which principle components account for the variance in the dataset."),
              
              hr(),
               
              shinycssloaders::withSpinner(plotOutput("PC_plot"), type = 8)
            ),
          ),
          
          fluidRow(
            column(12, align = "center",
                   br(),
                   br(),
                   hr(),
                   div(
                     class = "footer",
                     includeHTML("footer.Rhtml")
                   )
            )
          )
        ),
        
        #### UMAP --------------------------------------------------------------
        tabPanel("UMAP",
          sidebarLayout(
            sidebarPanel(
              h4("UMAP Analysis Options"),
              
              selectInput("UMAP_color", 
                          label = "Select Color Palette",
                          c("Set1",
                            "Set2",
                            "Set3",
                            "Pastel1",
                            "Pastel2",
                            "Paired",
                            "Dark2",
                            "Accent")),
              
              helpText("If there are more groups than the number of colors in the 
                       largest palette (12), the graph will default to a larger color palette with distinct hues."),
              
              br(),
              br(),
              
              textInput("UMAP_title",
                        label = "Customize Plot Title",
                        placeholder = "Type title here"),
              
              checkboxInput("UMAP_shapes",
                            label = "Make different groups different shapes"),
              
              checkboxInput("UMAP_labels",
                            label = "Label points",
                            value = TRUE),
              
              hr(),
              
              actionButton("UMAP_button", "Run Analysis",
                           class = "btn-block",
                           style="color: #fff; background-color: #9BD79A; border-color: 
                          #9BD79A; border-radius: 10px; border-width: 2px"),
              br(),
              br(),
              
              h6("This analysis may take a few minutes to load. \n"),
              
              hr(),
              
              h5("Download Plot"),
              
              radioButtons("UMAP_extension", 
                           label = "Save As:",
                           choices = c("PDF" = "pdf", 
                                       "PNG" = "png"),
                           inline = TRUE),
              
              downloadButton("UMAP_download", 
                             label = "Download UMAP Plot",
                             style="color: #61A6F9; background-color: #D8EAFF; border-color: 
                          #D8EAFF; border-radius: 10px; border-width: 2px")
            ),
            
            mainPanel(
              br(),
              h3("UMAP Plot"),
              
              helpText("UMAP requires both ", em("Annotation "), "and ", em("Data "), 
                       "file inputs. It also requires selection of ", em("Data Type "), 
                       "in the Data Upload tab."),
              hr(),
              
              shinycssloaders::withSpinner(plotOutput("UMAP_plot"), type = 8)
            )
          ),
          
          fluidRow(
            column(12, align = "center",
                   br(),
                   br(),
                   hr(),
                   div(
                     class = "footer",
                     includeHTML("footer.Rhtml")
                   )
            )
          )
        ),
        
        #### MD PLOT -----------------------------------------------------------
        tabPanel("MD Plot",
          sidebarLayout(
            sidebarPanel(
              h4("Analysis Options"),
              
              fluidRow(
                column(6, colourpicker::colourInput("md_down_color", "Select DOWN color", "royalblue", allowTransparent = T)),
                column(6, colourpicker::colourInput("md_up_color", "Select UP color", "red", allowTransparent = T))
              ),
              
              textInput("md_title",
                        label = "Customize Plot Title",
                        placeholder = "Type title here"),
              
              uiOutput("md_contrast"),
              
              checkboxInput("md_add_labels",
                            "Add labels to genes with highest/lowest logFC",
                            value = TRUE),
              
              uiOutput("md_label_number"),
              
              checkboxInput("md_label_specific",
                            "Add label to specific gene(s)"),
              
              uiOutput("md_label_specific_gene"),
              
              numericInput("md_fc_cutoff",
                           label = "Select FC cutoff",
                           value = 1,
                           step = .5),
              
              hr(),
              
              actionButton("md_button",
                           label = "Run Analysis",
                           style="color: #fff; background-color: #9BD79A; 
                          border-color: #9BD79A; border-radius: 10px; border-width: 2px"),
              
              br(),
              br(),
              
              h6("This analysis may take a few minutes to load. \n"),
              
              hr(),
              
              h5("Download Plot"),
              
              radioButtons("md_extension", 
                           label = "Save As:",
                           choices = c("PDF" = "pdf", 
                                       "PNG" = "png"),
                           inline = TRUE),
              
              downloadButton("md_download", 
                             label = "Download MD Plot",
                             style="color: #61A6F9; background-color: #D8EAFF; border-color: 
                          #D8EAFF; border-radius: 10px; border-width: 2px")
            ),
            
            mainPanel(
              
              br(),
              h3("Mean Difference Plot"),
              
              helpText("The MD plot requires both ", em("Annotation "), "and ", em("Data "), 
                       "file inputs. It also requires selection of ", em("Data Type "), 
                       "in the Data Upload tab."),
              hr(),
              
              shinycssloaders::withSpinner(plotOutput("md_plot"), type = 8)
            )
          ),
          
          fluidRow(
            column(12, align = "center",
                   br(),
                   br(),
                   hr(),
                   div(
                     class = "footer",
                     includeHTML("footer.Rhtml")
                   )
            )
          )
        ),
        
        #### CORRELATION -------------------------------------------------------
        tabPanel("Correlation",
          sidebarLayout(
            sidebarPanel(
              
              h4("Options"),
              
              fluidRow(
                column(6, colourpicker::colourInput("corr_neg_color", "Select NEG color", "royalblue", allowTransparent = T)),
                column(6, colourpicker::colourInput("corr_pos_color", "Select POS color", "red", allowTransparent = T))
              ),
              
              textInput("corr_title",
                        label = "Add title",
                        placeholder = "Type title here"),
              
              radioButtons("corr_shape",
                           label = "Shape of points",
                           choices = c("Square" = "square",
                                       "Circle" = "circle")),
              
              radioButtons("corr_type",
                          label = "Select graph orientation",
                          choices = c("Full" = "full",
                                      "Lower" = "lower",
                                      "Upper" = "upper")),
              
              helpText("'Full' will fill in the whole square, showing duplicate correlations, 
                       while 'Lower' and 'Upper' will show only the top or bottom half of correlations."),
              
              br(),
              br(),
              
              checkboxInput("corr_number",
                            label = "Add correlation values"),
              
              
              helpText("When this box is checked the calculated correlation values will be added to each square of the correlation plot."),
              
              br(),
              br(),
              
              actionButton("corr_button",
                           label = "Run Analysis",
                           style="color: #fff; background-color: #9BD79A; 
                          border-color: #9BD79A; border-radius: 10px; border-width: 2px"),
              br(),
              br(),
              
              h6("This analysis may take a few minutes to load. \n"),
              
              hr(),
              
              h5("Download Correlation Plot"),
              
              br(),
               
              downloadButton("corr_plot_download", 
                              label = "Download Plot",
                              style="color: #61A6F9; background-color: #D8EAFF; border-color: 
                  #D8EAFF; border-radius: 10px; border-width: 2px"),
              
              br(),
              br(),
              
              helpText("The correlation plot downloads as a PDF."),
              
              br(),
              
            ),
            
            mainPanel(
              br(),
              h3("Correlation"),
              
              helpText("The correlation plot uses a z-scaled version of the dataset which
                       has NA omission. The below correlation table can be downloaded as a CSV file."),
              
              hr(),
              
              shinycssloaders::withSpinner(plotOutput("corr_plot"), type = 8)
              
            )
          ),
          
          fluidRow(
            column(12, align = "center",
                   br(),
                   br(),
                   hr(),
                   div(
                     class = "footer",
                     includeHTML("footer.Rhtml")
                   )
            )
          )
        ),
        
        #### HEATMAP -----------------------------------------------------------
        tabPanel("Heatmap",
          sidebarLayout(
            sidebarPanel(
              h4("Heatmap Analysis Options"),
              selectInput("heatmap_color", 
                          label = "Select color palette for heatmap body",
                          c("Spectral",
                            "RdYlGn",
                            "RdYlBu",
                            "RdGy",
                            "RdBu",
                            "PuOr",
                            "PRGn",
                            "PiYG",
                            "BrBG")),
              
              selectInput("heatmap_group_color", 
                          label = "Select color palette for groups",
                          c("Emrld",
                            "Plasma",
                            "Heat",
                            "Sunset",
                            "Harmonic",
                            "Terrain",
                            "Purple-Yellow",
                            "Roma",
                            "Geyser",
                            "Viridis")),
              
              textInput("heatmap_title",
                        label = "Customize Plot Title",
                        placeholder = "Type title here"),
                     
              checkboxInput("kclustering",
                            "K-cluster rows"),
              
              checkboxInput("cluster_samples",
                            label = "Cluster columns",
                            value = TRUE),
              
              checkboxInput("zscore",
                            label = "Z-score across rows",
                            value = TRUE),
              
              selectInput("subset",
                          label = "Subset Data",
                          c("None" = "none",
                            "Subset by gene name" = "subset_genes",
                            "Subset by most variable genes/rows" = "subset_variable",
                            "Subset by feature" = "subset_features")),
              
              shinycssloaders::withSpinner(uiOutput("gene_list_select"), type = 7),
              
              uiOutput("gene_list_file"),
                     
              uiOutput("select_features"),
              
              uiOutput("select_variable"),
              
              hr(),
                     
              actionButton("heatmap_button",
                          label = "Run Analysis",
                          style = "color: #fff; background-color: #9BD79A; 
                          border-color: #9BD79A; border-radius: 10px; border-width: 2px"),
                     
              br(),
              br(),
                     
              h6("This analysis may take a few minutes to load. \n"),
                     
              hr(),
                     
              h5("Download Plot"),
                     
              radioButtons("heatmap_extension", 
                          label = "Save As:",
                          choices = c("PDF" = "pdf", 
                                      "PNG" = "png"),
                          inline = TRUE),
                     
              downloadButton("heatmap_download", 
                            label = "Download Heatmap Plot",
                            style="color: #61A6F9; background-color: #D8EAFF; border-color: 
                          #D8EAFF; border-radius: 10px; border-width: 2px"),
              
              br(),
              br(),
              
              h6("Please allow a few moments for downloading.")
              
              ),
                   
            mainPanel(
              br(),
              h3("Heatmap"),
                     
              helpText("Heatmap analysis requires both ", em("Annotation "), "and ", em("Data "), 
                      "file inputs. It also requires selection of ", em("Data Format "),
                      "and ", em("Data Type "), "in the Data Upload tab, and takes into 
                      account the selection of ", em("Log-Transform"), " for the expression set. "),
              
              hr(),
                     
              shinycssloaders::withSpinner(plotOutput("heatmap_plot"), type = 8)
              
            )
          ),
          
          fluidRow(
            column(12, align = "center",
                   br(),
                   br(),
                   hr(),
                   div(
                     class = "footer",
                     includeHTML("footer.Rhtml")
                   )
            )
          )
        ),
      )
    ),
    
    ## DATA MODELING ###########################################################
    tabPanel("3. Data Modeling",
      tabsetPanel(
        
        #### MODELING OVERVIEW -------------------------------------------------
        tabPanel("Modeling Overview",
                 br(),
                 fluidRow(
                   column(2),
                   column(8,
                    includeMarkdown("linearModeling.md")
                   ),
                   column(2)
                 ),
                 fluidRow(
                   column(12, align = "center",
                          br(),
                          br(),
                          hr(),
                          div(
                            class = "footer",
                            includeHTML("footer.Rhtml")
                          )
                   )
                 )
        ),
        
        #### LINEAR MODEL ------------------------------------------------------
        tabPanel("Linear Model",
          sidebarLayout(
           sidebarPanel(
             h4("Model Options"),
             h5("STEP ONE"),
             
             h6("Model factors"),
             
             uiOutput("linear_factors"),
             
             hr(),
             
             h6("Covariates"),
             
             checkboxInput("add_covariate",
                           label = "Add covariate(s) to the model"),
             
             uiOutput("choose_covariate"),
             
             helpText("Choose up to 3 columns to include as covariates in the model. This will not change the contrast choices."),
             
             hr(),
             
             h6("Time Series"),
             
             checkboxInput("include_timeseries",
                           label = "Include Time as an interaction term in the model"),
             
             uiOutput("choose_time"),
             
             uiOutput("time_type"),
             
             uiOutput("time_points"),
             
             uiOutput("time_df_explain"),
             
             hr(),
             
             h6("Contrasts"),
             
             checkboxInput("contrast_fit",
                           label = "Perform model fit on contrasts",
                           value = TRUE),
             
             uiOutput("choose_contrasts"),
             
             hr(),
             
             h6("Cutoff values for model results table overview"),
             
             fluidRow(
               column(6,
                      numericInput("results_logfc",
                                   label = "LogFC cutoff",
                                   value = 0.5)
                      ),
               column(6,
                      numericInput("results_pval",
                                   label = "Adj. P-Value cutoff",
                                   value = 0.05)
                      )
             ),
             
             actionButton("limma_button",
                          label = "1: Generate Fit",
                          style="color: #fff; background-color: #9BD79A; 
                          border-color: #9BD79A; border-radius: 10px; border-width: 2px"),
             
             br(),
             
             hr(),
             
             h5("STEP TWO"),
             
             uiOutput("number_genes"),
             
             uiOutput("coef_options"),
             
             hr(),
             
             h6("logFC Threshold"),
             
             checkboxInput("logfc_threshold",
                           label = "Subset by logFC threshold"),
             
             uiOutput("logfc_th"),
             
             actionButton("toptable_button",
                          label = "2: Generate Top Table",
                          style="color: #fff; background-color: #9BD79A; 
                          border-color: #9BD79A; border-radius: 10px; border-width: 2px"),
             
             br(),
             hr(),
             
             h5("Download Tables"),
             br(),
             
             fluidRow(
               column(6,
                      downloadButton("linear_results_table_download",
                                     label = "Download results",
                                     style="color: #61A6F9; background-color: #D8EAFF; border-color: 
                          #D8EAFF; border-radius: 10px; border-width: 2px"),
                      br(),
                      br(),
                      helpText("Generated in step one: run analysis.")
               ),
               column(6,
                      downloadButton("linear_table_download",
                                     label = "Download top table",
                                     style="color: #61A6F9; background-color: #D8EAFF; border-color: 
                          #D8EAFF; border-radius: 10px; border-width: 2px"),
                      br(),
                      br(),
                      helpText("Generated in step two: generate top table.")
                      )
             )
           ),
           
            mainPanel(
              column(9,
              br(),
              h3("Linear Model"),
             
              helpText("Linear Model Differential Analysis requires both ", em("Annotation "), "and ", em("Data "), 
                      "file inputs as well as Data Type and Data Format selections. If your input data does not have any replicates, 
                      the limma model will fail. In this case it is recommended to reference the MD plot under the ", 
                      code("Basic Analysis"), " header to assess logFC values."),
             
              br(),
              br(),
             
              helpText("Once calculated, the limma model will output a table with the following columns: ",
                      code("FEATURE_IDENTIFIER"), " is the gene/row name. ",
                      code("LOGFC"), " is the estimate of the log2-fold-change corresponding to the effect or contrast; ",
                      code("AVEEXPR"), " is the average log2-expression for the probe over all arrays and channels; ",
                      code("T"), " is the moderated t-statistic; ",
                      code("F"), " is the moderated F-statistic; ",
                      code("P.VALUE"), " is the raw p-value; ",
                      code("ADJ.P.VALUE"), " is the adjusted p-value or q-value; ",
                      code("B"), " is the log-odds that the gene is differentially expressed."),
             
              hr(),
              
              uiOutput("linear_model_results_header"),
              
              br(),
              
              tableOutput("linear_model_results"),
              
              uiOutput("linear_model_equation_header"),
              
              verbatimTextOutput("model_equation"),
              
              hr(),
              br(),
             
              shinycssloaders::withSpinner(DT::dataTableOutput("linear_model_table"), type = 8)
              ),
             
              column(3,
             
              )
            )
          ),
          
          fluidRow(
            column(12, align = "center",
                   br(),
                   br(),
                   hr(),
                   div(
                     class = "footer",
                     includeHTML("footer.Rhtml")
                   )
            )
          )
        ),
        
        #### VOLCANO PLOT ------------------------------------------------------
        tabPanel("Volcano Plot",
          sidebarLayout(
            sidebarPanel(
              h4("Plot Options"),
              
              fluidRow(
                column(6, colourpicker::colourInput("volcano_down_color", "Select DOWN color", "royalblue", allowTransparent = T)),
                column(6, colourpicker::colourInput("volcano_up_color", "Select UP color", "red", allowTransparent = T))
              ),
              
              textInput("volcano_title",
                        label = "Customize Plot Title",
                        placeholder = "Type title here"),
              
              uiOutput("volcano_coef_options"),
              
              uiOutput("volcano_number_genes"),
              
              checkboxInput("volcano_labels",
                            label = "Add labels to genes with highest/lowest logFC",
                            value = TRUE),
              
              uiOutput("volcano_label_number"),
              
              checkboxInput("volcano_label_specific",
                            label = "Add label to specific gene(s)"),
              
              uiOutput("volcano_label_specific_gene"),
              
              radioButtons("volcano_yaxis",
                          label = "Choose Y-Axis Variable",
                          choices = c("P-Value" = "P.Value",
                                      "Adjusted P-Value" = "adj.P.Val")),
              
              fluidRow(
                column(6,
                       numericInput("volcano_top_fc",
                                    label = "LogFC Cut Off",
                                    value = 0.5,
                                    step = 0.5,
                                    min = 0)
                ),
                
                column(6,
                       numericInput("volcano_top_values",
                                    label = "P-Value Cut Off",
                                    value = 0.05,
                                    step = 0.01,
                                    min = 0,
                                    max = 1)
                )
              ),
              
              hr(),
              
              helpText("The limma model in the ", code("Linear Model"), " tab must be
                       run prior to generating the volcano plot."),
              br(),
              br(),
              
              actionButton("volcano_button",
                           label = "Run Analysis",
                           style="color: #fff; background-color: #9BD79A; 
                          border-color: #9BD79A; border-radius: 10px; border-width: 2px"),
              
              br(),
              br(),
              
              h6("This analysis may take a few minutes to load. \n"),
              
              hr(),
              
              h5("Download Plot"),
              
              radioButtons("volcano_extension", 
                           label = "Save As:",
                           choices = c("PDF" = "pdf", 
                                       "PNG" = "png"),
                           inline = TRUE),
              
              downloadButton("volcano_download", 
                             label = "Download Volcano Plot",
                             style="color: #61A6F9; background-color: #D8EAFF; border-color: 
                          #D8EAFF; border-radius: 10px; border-width: 2px")
            ),
            
            mainPanel(
              br(),
              h3("Volcano Plot"),
              
              helpText("The Volcano Plot is generated based on the input in the ", em("Linear Model "), "tab. It requires ", em("Annotation "), "and ", em("Data "), 
                       "file inputs as well as selection of ", em("Data Type "), "in the ", em("Data Upload "), "tab."),
              
              br(),
              br(),
              
              helpText("The contrast/coefficient options depend on the model design in step one in the ", em("Linear Model"),
                       " tab."),
              
              hr(),
              
              shinycssloaders::withSpinner(plotOutput("volcano_plot"), type = 8)
              
            )
          ),
          
          fluidRow(
            column(12, align = "center",
                   br(),
                   br(),
                   hr(),
                   div(
                     class = "footer",
                     includeHTML("footer.Rhtml")
                   )
            )
          )
        ),
        
        #### INTERACTIVE VOLCANO PLOT ------------------------------------------
        tabPanel("Interactive Volcano Plot",
          sidebarLayout(
            sidebarPanel(
              h4("Plot Options"),
              
              fluidRow(
                column(6, colourpicker::colourInput("int_volcano_down_color", "Select DOWN color", "royalblue", allowTransparent = T)),
                column(6, colourpicker::colourInput("int_volcano_up_color", "Select UP color", "red", allowTransparent = T))
              ),
              
              uiOutput("int_volcano_coef_options"),
              
              uiOutput("int_volcano_number_genes"),
              
              radioButtons("int_volcano_yaxis",
                           label = "Choose Y-Axis Variable",
                           choices = c("P-Value" = "P.Value",
                                       "Adjusted P-Value" = "adj.P.Val")),
              fluidRow(
                column(6,
                       numericInput("int_volcano_top_fc",
                                    label = "LogFC Cut Off",
                                    value = 0.5,
                                    step = 0.5,
                                    min = 0)
                ),
                 
                column(6,
                       numericInput("int_volcano_top_values",
                                    label = "P-Value Cut Off",
                                    value = 0.05,
                                    step = 0.01,
                                    min = 0,
                                    max = 1)
                )
              ),
             
              hr(),
             
              helpText("The limma model in the ", code("Linear Model"), " tab must be
                       run prior to generating the volcano plot."),
              br(),
              br(),
             
              actionButton("int_volcano_button",
                           label = "Generate Plotly",
                           style="color: #fff; background-color: #9BD79A; 
                                  border-color: #9BD79A; border-radius: 10px; border-width: 2px"),
             
              br(),
              
              hr(),
              
              h5("Glimma XY Plot"),
              
              helpText("The glimma package allows for generation of an interactive volcano plot with an x-axis of logFC and a y-axis of logOdds. It also 
                       generates a corresponding plot showing Normalization Intensity by group and a table of associated values. If you click on the ", em("Generate Glimma Tab"),
                       " button below a new tab will open with these interactive plots."),
              
              br(),
              br(),
              
              helpText("Please note that the cutoff for significance in the Glimma volcano plot is not customizable and differs
                       from the cutoff for the other volcano plots. As such, the Glimma plot may show a different number of hits."),
              
              br(),
              br(),
              
              actionButton("glimma_volcano_button",
                           label = "Generate Glimma Tab",
                           style="color: #fff; background-color: #9BD79A; 
                                  border-color: #9BD79A; border-radius: 10px; border-width: 2px"),
              
              br(),
              br(),
              
              h6("Please allow a few moments for loading.")
              
            ),
            
            mainPanel(
              br(),
              h3("Interactive Volcano Plot"),
              
              helpText("The Volcano Plot is generated based on the input in the ", em("Linear Model "), "tab. It requires ", em("Annotation "), "and ", em("Data "), 
                       "file inputs as well as selection of ", em("Data Type "), "in the ", em("Data Upload "), "tab."),
              
              br(),
              br(),
              
              helpText("The plotly graph object generated below allows for interaction with the generated volcano plot via hovering over data points, clicking & dragging,
                       or use of the menu in the top right. PNG images can be downloaded of the graph, but will not include labels or information from the hovertabs."),
              
              hr(),
              br(),
              
              shinycssloaders::withSpinner(plotly::plotlyOutput("volcano_plot_interactive", height = "700px",), type = 8)
            )
          ),
          
          fluidRow(
            column(12, align = "center",
                   br(),
                   br(),
                   hr(),
                   div(
                     class = "footer",
                     includeHTML("footer.Rhtml")
                   )
            )
          )
        ),
        
        #### DIFFERENTIAL HEATMAP ----------------------------------------------
        tabPanel("Differential Heatmap",
                 sidebarLayout(
                   sidebarPanel(
                     h4("Heatmap Analysis Options"),
                     selectInput("dif_heatmap_color", 
                                 label = "Select color palette for heatmap body",
                                 c("Spectral",
                                   "RdYlGn",
                                   "RdYlBu",
                                   "RdGy",
                                   "RdBu",
                                   "PuOr",
                                   "PRGn",
                                   "PiYG",
                                   "BrBG")),
                     
                     selectInput("dif_heatmap_group_color", 
                                 label = "Select color palette for groups",
                                 c("Emrld",
                                   "Plasma",
                                   "Heat",
                                   "Sunset",
                                   "Harmonic",
                                   "Terrain",
                                   "Purple-Yellow",
                                   "Roma",
                                   "Geyser",
                                   "Viridis")),
                     
                     textInput("dif_heatmap_title",
                               label = "Customize Plot Title",
                               placeholder = "Type title here"),
                     
                     uiOutput("dif_heatmap_coef_options"),
                     
                     checkboxInput("dif_kclustering",
                                   "K-cluster rows"),
                     
                     checkboxInput("dif_cluster_samples",
                                   label = "Cluster columns",
                                   value = TRUE),
                     
                     checkboxInput("dif_zscore",
                                   label = "Z-score across rows",
                                   value = TRUE),
                     
                     numericInput("dif_number_rows",
                                  label = "Select number of rows to include in heatmap",
                                  value = 50,
                                  min = 1),
                     
                     hr(),
                     
                     actionButton("dif_heatmap_button",
                                  label = "Run Analysis",
                                  style = "color: #fff; background-color: #9BD79A; 
                          border-color: #9BD79A; border-radius: 10px; border-width: 2px"),
                     
                     br(),
                     br(),
                     
                     h6("This analysis may take a few minutes to load. \n"),
                     
                     hr(),
                     
                     h5("Download Plot"),
                     
                     radioButtons("dif_heatmap_extension", 
                                  label = "Save As:",
                                  choices = c("PDF" = "pdf", 
                                              "PNG" = "png"),
                                  inline = TRUE),
                     
                     downloadButton("dif_heatmap_download", 
                                    label = "Download Heatmap Plot",
                                    style="color: #61A6F9; background-color: #D8EAFF; border-color: 
                          #D8EAFF; border-radius: 10px; border-width: 2px"),
                     
                     br(),
                     br(),
                     
                     h6("Please allow a few moments for downloading.")
                     
                   ),
                   
                   mainPanel(
                     br(),
                     h3("Heatmap"),
                     
                     helpText("Heatmap analysis requires both ", em("Annotation "), "and ", em("Data "), 
                              "file inputs. It also requires selection of ", em("Data Format "),
                              "and ", em("Data Type "), "in the Data Upload tab, and takes into 
                              account the selection of ", em("Log-Transform"), " for the expression set. "),
                     
                     br(),
                     br(),
                     
                     helpText("This heatmap uses the data output from the linear model in the top table.
                              It uses the rows with the largest logFC absolute values."),
                     
                     hr(),
                     
                     shinycssloaders::withSpinner(plotOutput("dif_heatmap_plot"), type = 8)
                     
                   )
                 ),
                 
                 fluidRow(
                   column(12, align = "center",
                          br(),
                          br(),
                          hr(),
                          div(
                            class = "footer",
                            includeHTML("footer.Rhtml")
                          )
                   )
                 )
        )
      )
    ),
    
    ## ENRICHMENT ##############################################################
    tabPanel("4. Enrichment",
             tabsetPanel(
               
               ### OVERVIEW ----------------------------------------------------
               tabPanel("Enrichment Overview",
                        br(),
                        fluidRow(
                          column(2),
                          column(8,
                                 includeMarkdown("enrichmentOverview.md")
                          ),
                          column(2)
                        ),
                        fluidRow(
                          column(12, align = "center",
                                 br(),
                                 br(),
                                 hr(),
                                 div(
                                   class = "footer",
                                   includeHTML("footer.Rhtml")
                                 )
                          )
                        )
                 
               ),
               
               ### CLUSTER PROFILER --------------------------------------------
               tabPanel("Cluster Profiler",
                        sidebarLayout(
                          sidebarPanel(
                            h4("Enrichment Options"),
                            
                            helpText("The enrichment function uses the output of the data modeling linear model to
                          create a ranked gene list. It will not function if the linear model has not been run."),
                            
                            br(),
                            br(),
                            
                            h6("Species"),
                            
                            selectInput("species", 
                                        label = "Select species",
                                        c("Human (9606)" = "human",
                                          "Mouse (10090)" = "mouse",
                                          "Yeast (559292)" = "yeast",
                                          "Zebrafish (7955)" = "zebrafish",
                                          "C. elegans (6239)" = "celegans",
                                          "Fruit Fly (7227)" = "fruitfly",
                                          "Rat (10116)" = "rat",
                                          "Other")),
                            
                            hr(),
                            
                            h6("Database"),
                            
                            uiOutput("gmts"),
                            
                            hr(),
                            
                            h6("Enrichment"),
                            
                            radioButtons("enrichment_type",
                                         label = "Enrichment type to perform",
                                         choices = c("Gene Set Enrichment Analysis (GSEA)" = "GSEA",
                                                     "Over-Representation Analysis (ORA)" = "enricher")),
                            
                            hr(),
                            
                            h6("Contrast"),
                            
                            uiOutput("enrichment_coef_options"),
                            
                            actionButton("enrichment_button",
                                         label = "Run Enrichment Analysis",
                                         style = "color: #fff; background-color: #9BD79A; 
                          border-color: #9BD79A; border-radius: 10px; border-width: 2px"),
                            
                            br(),
                            br(),
                            
                            h6("Please allow a few moments for the enrichment analysis to run."),
                            
                            hr(),
                            
                            h5("Download Table"),
                            
                            downloadButton("enriched_table_download",
                                           label = "Download Enrichment Table",
                                           style = "color: #61A6F9; background-color: #D8EAFF; border-color: 
                                                            #D8EAFF; border-radius: 10px; border-width: 2px"),
                            
                            br(),
                            br(),
                            
                            helpText("The enrichment table includes the values shown in the table generated in the ", code("Overview"),
                                     " tab as well as a column containing a list of gene symbols that intersect with each pathway. This
                          table downloads as a CSV file."),
                            
                            width = 3),
                          
                          mainPanel(
                            tabsetPanel(
                              
                              #### OVERVIEW ------------------------------------
                              tabPanel("Overview",
                                       
                                       column(10,
                                              br(),
                                              
                                              h4("Enriched Pathways"),
                                              
                                              helpText("The ", code("clusterProfiler"), " package is used for pathway enrichment."),
                                              
                                              br(),
                                              br(),
                                              
                                              helpText("The table generated below when enrichment analysis is run contains the names of 
                                     the enriched pathways based on the species and database selections. For both GSEA and 
                                     general enrichment types there is an additional column for the adjusted p-value associated
                                     with each pathway. GSEA enrichment will have another column with normalized enrichment 
                                     score values, while general enrichment will have a column with the gene count ratio."),
                                              
                                              br(),
                                              br(),
                                              
                                              helpText("The plots in the additional tabs (except the heatmap) are generated using the ", code("enrichplot"), "package,
                                            available on Bioconductor and maintained by Guangchuang Yu."),
                                              hr(),
                                              
                                              shinycssloaders::withSpinner(DT::dataTableOutput("enriched_table"), type = 8)
                                       ),
                                       column(2)
                              ),
                              
                              #### DOT PLOT ------------------------------------
                              tabPanel("Dot Plot",
                                       
                                       br(),
                                       h4("Dot Plot"),
                                       
                                       helpText("The dot plot depicts the enrichment scores (p-values) via color-coded points,
                                     and gene ratio as position relative to the x-axis. Gene count is represented by the 
                                     size of the point."),
                                       hr(),
                                       
                                       fluidRow(
                                         column(9,
                                                shinycssloaders::withSpinner(plotOutput("enriched_dotplot"), type = 8)
                                         ),
                                         column(3,
                                                h5("Options"),
                                                
                                                textInput("dotplot_title",
                                                          label = "Customize plot title",
                                                          value = "Enriched Dot Plot"),
                                                
                                                numericInput("dotplot_categories",
                                                             label = "Select number of pathways to include",
                                                             value = 15),
                                                
                                                hr(),
                                                
                                                actionButton("dotplot_button",
                                                             label = "Generate Plot",
                                                             style = "color: #fff; background-color: #9BD79A; 
                                                          border-color: #9BD79A; border-radius: 10px; border-width: 2px"),
                                                
                                                br(),
                                                br(),
                                                
                                                downloadButton("dotplot_download",
                                                               label = "Download",
                                                               style = "color: #61A6F9; background-color: #D8EAFF; border-color: 
                                                            #D8EAFF; border-radius: 10px; border-width: 2px"),
                                         )
                                       )
                              ),
                              
                              #### CNET ----------------------------------------
                              tabPanel("Gene Set Network",
                                       br(),
                                       h4("Gene-Concept Network"),
                                       
                                       helpText("The circular network plot is constructed based on hierarchical clustering of enriched 
                                     terms, and connections represent similarities between terms. Nodes in the plot represent 
                                     individual terms or gene sets, and the connections between nodes indicate relationships 
                                     or similarities between them. The color of nodes can be customized based on adjusted 
                                     p-values or combined scores. The size of nodes is proportional to the number of genes 
                                     in each category."),
                                       hr(),
                                       
                                       fluidRow(
                                         column(9,
                                                shinycssloaders::withSpinner(plotOutput("enriched_cnetplot"), type = 8)
                                         ),
                                         column(3,
                                                h5("Options"),
                                                
                                                textInput("cnet_title",
                                                          label = "Customize plot title",
                                                          value = "Enriched Network Plot"),
                                                
                                                numericInput("cnet_categories",
                                                             label = "Select number of categories to include",
                                                             value = 5),
                                                
                                                selectInput("cnet_layout",
                                                            label = "Plot layout",
                                                            choices = c("Kamada-Kawai" = "kk",
                                                                        "Large Graph" = "lgl",
                                                                        "Circle" = "circle",
                                                                        "Fruchterman-Reingold" = "fr",
                                                                        "Grid" = "grid",
                                                                        "Star" = "star",
                                                                        "Gem" = "gem",
                                                                        "Multidimensional Scaling" = "mds",
                                                                        "Distributed Recursive" = "drl",
                                                                        "Random" = "randomly")),
                                                
                                                selectInput("cnet_labels",
                                                            label = "Label",
                                                            choices = c("Category Nodes Only" = "category",
                                                                        "Gene Nodes Only"= "gene",
                                                                        "Category & Gene Nodes" = "all",
                                                                        "None" = "none")),
                                                
                                                hr(),
                                                
                                                actionButton("cnetplot_button",
                                                             label = "Generate Plot",
                                                             style = "color: #fff; background-color: #9BD79A; 
                                                          border-color: #9BD79A; border-radius: 10px; border-width: 2px"),
                                                
                                                br(),
                                                br(),
                                                
                                                downloadButton("cnetplot_download",
                                                               label = "Download",
                                                               style = "color: #61A6F9; background-color: #D8EAFF; border-color: 
                                                            #D8EAFF; border-radius: 10px; border-width: 2px"),
                                         )
                                       )
                              ),
                              
                              #### EMAP ----------------------------------------
                              tabPanel("Enrichment Map",
                                       br(),
                                       h4("Enrichment Map"),
                                       
                                       helpText("According to the Yu lab, \"the enrichment map organizes enriched terms into 
                                     a network with edges connecting overlapping gene sets. In this way, mutually 
                                     overlapping gene sets are tend to cluster together, making it easy to identify 
                                     functional modules.\""),
                                       hr(),
                                       
                                       fluidRow(
                                         column(9,
                                                shinycssloaders::withSpinner(plotOutput("enriched_emapplot"), type = 8)
                                         ),
                                         column(3,
                                                h5("Options"),
                                                
                                                textInput("emap_title",
                                                          label = "Customize plot title",
                                                          value = "Enrichment Map"),
                                                
                                                numericInput("emap_categories",
                                                             label = "Select number of categories to include",
                                                             value = 30),
                                                
                                                hr(),
                                                
                                                actionButton("emap_button",
                                                             label = "Generate Plot",
                                                             style = "color: #fff; background-color: #9BD79A; 
                                                          border-color: #9BD79A; border-radius: 10px; border-width: 2px"),
                                                
                                                br(),
                                                br(),
                                                
                                                downloadButton("emap_download",
                                                               label = "Download",
                                                               style = "color: #61A6F9; background-color: #D8EAFF; border-color: 
                                                            #D8EAFF; border-radius: 10px; border-width: 2px"),
                                         )
                                       )
                              ),
                              
                              #### UPSET PLOT ----------------------------------
                              tabPanel("UpSet Plot",
                                       br(),
                                       h4("UpSet Plot"),
                                       
                                       helpText("An UpSet plot is a visual representation that displays the intersections 
                                     and overlaps between multiple sets in a dataset. It uses bars or boxes 
                                     to represent individual sets and vertical lines to indicate the presence 
                                     or absence of elements in specific intersections. UpSet plots are particularly 
                                     useful for understanding the combinations of sets and visualizing the 
                                     relationships between them."),
                                       hr(),
                                       
                                       fluidRow(
                                         column(9,
                                                shinycssloaders::withSpinner(plotOutput("enriched_upsetplot"), type = 8)
                                         ),
                                         column(3,
                                                h5("Options"),
                                                
                                                textInput("upset_title",
                                                          label = "Customize plot title",
                                                          value = "Enrichment UpSet Plot"),
                                                
                                                numericInput("upset_categories",
                                                             label = "Select number of categories to include",
                                                             value = 10),
                                                
                                                hr(),
                                                
                                                actionButton("upset_button",
                                                             label = "Generate Plot",
                                                             style = "color: #fff; background-color: #9BD79A; 
                                                          border-color: #9BD79A; border-radius: 10px; border-width: 2px"),
                                                
                                                br(),
                                                br(),
                                                
                                                downloadButton("upset_download",
                                                               label = "Download",
                                                               style = "color: #61A6F9; background-color: #D8EAFF; border-color: 
                                                            #D8EAFF; border-radius: 10px; border-width: 2px"),
                                         )
                                       )
                              ),
                              
                              #### TREE ----------------------------------------
                              tabPanel("Tree Diagram",
                                       br(),
                                       h4("Tree Diagram"),
                                       
                                       helpText("The tree plot is a funnctional grouping tree diagram for
                                     enrichment result of over-representation test or gene set 
                                     enrichment analysis."),
                                       hr(),
                                       
                                       fluidRow(
                                         column(9,
                                                shinycssloaders::withSpinner(plotOutput("enriched_treeplot"), type = 8)
                                         ),
                                         column(3,
                                                h5("Options"),
                                                
                                                textInput("tree_title",
                                                          label = "Customize plot title",
                                                          value = "Enriched Tree Plot"),
                                                
                                                numericInput("tree_categories",
                                                             label = "Select number of categories to include",
                                                             value = 30),
                                                
                                                numericInput("tree_clusters",
                                                             label = "Select number of clusters to include",
                                                             value = 5),
                                                
                                                hr(),
                                                
                                                actionButton("tree_button",
                                                             label = "Generate Plot",
                                                             style = "color: #fff; background-color: #9BD79A; 
                                                          border-color: #9BD79A; border-radius: 10px; border-width: 2px"),
                                                
                                                br(),
                                                br(),
                                                
                                                downloadButton("tree_download",
                                                               label = "Download",
                                                               style = "color: #61A6F9; background-color: #D8EAFF; border-color: 
                                                            #D8EAFF; border-radius: 10px; border-width: 2px"),
                                         )
                                       )
                              ),
                              
                              #### HEATMAP -------------------------------------
                              tabPanel("Heatmap",
                                       fluidRow(
                                         column(10,
                                                br(),
                                                h4("Enriched Heatmap"),
                                                
                                                helpText("The heatmap is a visual representation of gene expression patterns 
                                                         across different samples or conditions based on enrichment analysis
                                                         results. In the heatmap, rows correspond to enriched pathways, and 
                                                         columns represent genes. The color intensity in the heatmap reflects 
                                                         the expression levels of genes, highlighting how they are enriched or 
                                                         depleted in specific gene sets identified by enrichment. This type of 
                                                         heatmap helps visualize the coordinated changes in gene expression 
                                                         associated with biological pathways or gene sets, providing insights 
                                                         into the functional relevance of different experimental conditions.")
                                                ),
                                         column(2,
                                                br(),
                                                br(),
                                                br(),
                                                actionButton("enriched_heatmap_button",
                                                             label = "Generate Plot",
                                                             style = "color: #fff; background-color: #9BD79A; 
                                                          border-color: #9BD79A; border-radius: 10px; border-width: 2px")
                                                )
                                       ),
                                       hr(),
                                       br(),
                                       shinycssloaders::withSpinner(plotlyOutput("enriched_heatplot"), type = 8)
                              ),
                            ), width = 9
                          ),
                        )
                        ),
               
               ### PHOSPHO ENRICHMENT ------------------------------------------
               tabPanel("Phosphosite PSEA",
                        ),
               
               tabPanel("Kinase KSEA",
                        )
             )
             
    ),
    
    ## GENERATE REPORT #########################################################
    tabPanel("Generate Report",
             fluidRow(
               h3("RMarkdown Report"),
               
               br(),
               br(),
               
               column(6,
                      p("The following allows for the generation of an RMarkdown report based on customized inputs."),
                      
                      br(),
                      
                      p("All inputs and selections must be made in the ", em("Data Upload"), " tab in order for the 
                        report to generate. For relevant plots, such as MD, volcano, and differential heatmaps, a plot will
                        be generated for each possible contrast based on the groupings selected."),
                      
                      br(),
                      
                      p("The report is a great way to get an overall view of trends in your data. More detailed and customizable versions of the
                        plots and data in the report can be generated in their specific portion of Omics Notebook, i.e. if there is a specific volcano plot
                        that you want to explore further from the report, you can generate it in the ", em("Data Modeling"), " tab as an interactive plot,
                        or select specific gene names to display in a static volcano plot."),
                      
                      br(),
                      
                      p("Please allow a few minutes for download. Load time will increase when more contrasts are selected."),
                      
                      br(),
                      
                      downloadButton("report",
                                     label = "Download Report",
                                     style = "color: #61A6F9; background-color: #D8EAFF; border-color: 
                                             #D8EAFF; border-radius: 10px; border-width: 2px"),
                      br(),
                      br(),
                      
                      helpText("Once the model design has been completed and all desired output options have been adjusted, click the ",
                               em("Download Report"), " button to generate the HTML report. A loading bar will appear showing progression of
                               the download."),
                      br(),
                      br(),
                      hr(),
                      
                      h4("Excel Data Summary"),
                      
                      br(),
                      
                      p("Additionally an Excel data summary can be generated separately which contains a summary of information and calculations
                        for each feature, such as: intensity values for each group, uniprot links, log fold change values for each contrast, and 
                        heatmap values."),
                      
                      br(),
                      
                      downloadButton("excel_summary_download_button",
                                     label = "Download Excel Summary",
                                     style = "color: #61A6F9; background-color: #D8EAFF; border-color: 
                          #D8EAFF; border-radius: 10px; border-width: 2px"),
                      
                      br(),
                      br()
                      ),
               column(6,
                      h5("Model Design [REQUIRED]"),
                      
                      br(),
                      
                      helpText("Selection of groups to include in the model for the report is necessary for generation."),
                      
                      br(),
                      br(),
                      
                      fluidRow(
                        column(6, uiOutput("report_linear_factors")),
                        column(6, uiOutput("report_reference_factor"),
                               uiOutput("report_choose_contrasts"))
                      ),
                      
                      helpText("The reference factor is only used if the checkbox for comparing all groups to one another 
                                 (below) is not checked. In that case, the reference factor will be the group all other groups are compared to."),
                      
                      br(),
                      br(),
                      
                      fluidRow(
                        column(6,
                               checkboxInput("report_contrast_fit",
                                             label = "Perform model fit on contrasts",
                                             value = TRUE),
                               
                               checkboxInput("report_add_covariate",
                                             label = "Add covariate(s) to the model"),
                               
                               uiOutput("report_covariate_col")
                               ),
                        column(6,
                               checkboxInput("report_include_timeseries",
                                             label = "Include time as interaction term"),
                               
                               uiOutput("report_choose_time"),
                               
                               uiOutput("report_time_type"),
                               
                               uiOutput("report_time_points")
                               )
                      ),
                      
                      hr(),
                      
                      h5("Enrichment"),
                      
                      fluidRow(
                        column(6,
                               selectInput("report_species", 
                                           label = "Select species",
                                           c("Human (9606)" = "human",
                                             "Mouse (10090)" = "mouse",
                                             "Yeast (559292)" = "yeast",
                                             "Zebrafish (7955)" = "zebrafish",
                                             "C. elegans (6239)" = "celegans",
                                             "Fruit Fly (7227)" = "fruitfly",
                                             "Rat (10116)" = "rat",
                                             "Other")),
                               
                               uiOutput("report_gmts"),
                               
                               radioButtons("report_enrichment_type",
                                            label = "Enrichment type to perform",
                                            choices = c("Gene Set Enrichment Analysis (GSEA)" = "GSEA",
                                                        "Over-Representation Analysis (ORA)" = "enricher")),
                               
                               ),
                        column(6,
                               checkboxGroupInput("report_enrichplots",
                                                  label = "Include enrichment plots",
                                                  choices = c("Dot Plot" = "dot",
                                                              "Gene Set Network" = "cnet",
                                                              "Enrichment Map" = "emap",
                                                              "UpSet Plot" = "upset",
                                                              "Tree Diagram" = "tree",
                                                              "Heatmap" = "heatmap"),
                                                  selected = c("dot", "cnet", "emap", "upset", "tree", "heatmap"))
                               )
                      )
                      ),
               
               br(),
               hr(),
               
               h4("Output Options"),
               
               br(),
               
               helpText("All checked sections will be included in the report. The report is 
                        downloaded in HTML format with all plots embedded. Plots can be saved/opened 
                        individually via right-click. MD plots, volcano plots, and differential heatmaps
                        will only be generated for the contrasts selected in the model design, or if 
                        contrasts are not selected, on each group individually."),
               br(),
               br(),
               br(),
               
               hr(),
               br(),
             ),
             fluidRow(
               column(3,
                      h5("Normalization"),
                      checkboxInput("report_norm",
                                    label = "Include normalization plots",
                                    value = TRUE),
                      hr(),
                      uiOutput("report_norm_options")),
               column(3,
                      h5("PCA"),
                      checkboxInput("report_PCA",
                                    label = "Include PCA plot",
                                    value = TRUE),
                      hr(),
                      uiOutput("report_PCA_options")),
               column(3,
                      h5("UMAP"),
                      checkboxInput("report_UMAP",
                                    label = "Include UMAP",
                                    value = TRUE),
                      hr(),
                      uiOutput("report_UMAP_options")),
               column(3,
                      h5("MD Plot"),
                      checkboxInput("report_MD",
                                    label = "Include MD plot",
                                    value = TRUE),
                      hr(),
                      uiOutput("report_MD_options")),
               hr()
             ),
             fluidRow(
               column(3,
                      h5("Correlation"),
                      checkboxInput("report_corr",
                                    label = "Include correlation plot",
                                    value = TRUE),
                      hr(),
                      uiOutput("report_corr_options")),
               column(3,
                      h5("Heatmap"),
                      checkboxInput("report_heatmap",
                                    label = "Include heatmap",
                                    value = TRUE),
                      hr(),
                      uiOutput("report_heatmap_options")),
               column(3,
                      h5("Volcano Plots"),
                      checkboxInput("report_volcano",
                                    label = "Include volcano plots",
                                    value = TRUE),
                      hr(),
                      uiOutput("report_volcano_options")),
               column(3,
                      h5("Differential Heatmaps"),
                      checkboxInput("report_dif_heatmap",
                                    label = "Include differential heatmaps",
                                    value = TRUE),
                      hr(),
                      uiOutput("report_dif_heatmap_options")),
             ),
             fluidRow(
               column(12, align = "center",
                      br(),
                      br(),
                      hr(),
                      div(
                        class = "footer",
                        includeHTML("footer.Rhtml")
                      )
               )
             )
      
    )
  )
)
