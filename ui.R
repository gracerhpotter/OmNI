
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
source("performNormalization.R")
source("saveXLSX.R")
source("enrichment.R")
source("drawCorrelation.R")
source("drawNormalization.R")
source("Sscore.R")
source("PCSF.R")
source("pullExpected.R")

# inactivity <- "function idleTimer() {
#   var t = setTimeout(logout, 3600000);
#   window.onmousemove = resetTimer; // catches mouse movements
#   window.onmousedown = resetTimer; // catches mouse movements
#   window.onclick = resetTimer;     // catches mouse clicks
#   window.onscroll = resetTimer;    // catches scrolling
#   window.onkeypress = resetTimer;  //catches keyboard actions
# 
#   function logout() {
#     window.close();  //close the window
#   }
# 
#   function resetTimer() {
#     clearTimeout(t);
#     t = setTimeout(logout, 3600000);  // time is in milliseconds (1000 is 1 second)
#   }
# }
# idleTimer();"

ui <- fluidPage(
  
  shinyjs::useShinyjs(),
  
  tags$head(
    
    # format the loading bar pop-up
    tags$style(".progress-bar{background-color:#61A6F9;}
               .tabbable > .nav > li > a {margin-bottom:50px;}",
               ".shiny-notification {
               height: 100px;
               width: 275px;}"),
    # position:fixed;
    # top: calc(50% - 50px);;
    # left: calc(50% - 125px);;
    
    # make it so that bslib::card objects allow drop downs to overflow over edge
    tags$style(HTML("
    .bslib-card, .tab-content, .tab-pane, .card-body {
      overflow: visible !important;
    }")),
    
    # make it so sidebar goes fills page height
    # ISSUE: goes over footer and looks bad
    
    # tags$style(HTML("hr {border-top: 1px solid #000000;}
    #             .well {height: 1200px;}")),
    
    # make app close after idle period
    # tags$script(inactivity)
  ),
  
  
  navbarPage("OmNI",
    theme = bslib::bs_theme(bootswatch = "lux"),

    ## ABOUT ###################################################################
    navbarMenu("About",
      
      #### OVERVIEW ------------------------------------------------------------
      tabPanel("Overview",
        br(),
        fluidRow(
          column(2),
          column(8,
            includeMarkdown("gettingStarted.md")
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
      ),
      
      #### LIMMA STATISTICS INFO -----------------------------------------------
      tabPanel(title = "Linear Modeling",
               br(),
               fluidRow(
                 column(2),
                 column(8,
                        includeMarkdown("linearModeling.md"),
                        includeMarkdown("limmaStats.md")
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
      
      #### ENRICHMENT -----------------------------------------------
      tabPanel(title = "Enrichment",
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
      
      #### GLOSSARY ---------------------------------------------------------
      # tabPanel(title = "Glossary",
      #          br(),
      #          fluidRow(
      #            column(2),
      #            column(8,
      #                   includeMarkdown("glossary.md")
      #            ),
      #            column(2)
      #          ),
      #          fluidRow(
      #            column(12, align = "center",
      #                   br(),
      #                   br(),
      #                   hr(),
      #                   div(
      #                     class = "footer",
      #                     includeHTML("footer.Rhtml")
      #                   )
      #            )
      #          )
      # ),
      
      #### DOCUMENTATION -------------------------------------------------------
      tabPanel(title = "Documentation",
               br(),
               fluidRow(
                 column(2),
                 column(8,
                        br(),
                        
                        includeMarkdown("methodsReference.md"),
                        
                        h5("Session Information"),
                        
                        verbatimTextOutput("session_info")
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
        
        ### DATA UPLOAD -------------------------------------------------------
        tabPanel("Data Upload",
          sidebarLayout(
            sidebarPanel(width = 4,
              
              h4("Upload Data"),
              
              helpText("The selections made in this tab will be used in all other analyses in Omics Notebook."),
              
              hr(),
              
              h6("Data File",
                 bslib::tooltip(bsicons::bs_icon("info-circle"),
                                "Data files must be a TXT, TSV, or CSV filetype.",
                                placement = "right")),
              
              fileInput("data_file",
                        label = "Select file(s)", 
                        buttonLabel = "Browse", 
                        multiple = TRUE,
                        placeholder = "No data file(s) selected",
                        accept = c("text",
                                   "text/tab-separated-values,text/plain",
                                   ".tsv",
                                   ".csv",
                                   "text/csv")),
              
              h6("Annotation File",
                 bslib::tooltip(bsicons::bs_icon("info-circle"),
                                "The annotation file must be an XLSX or XLS filetype.",
                                placement = "right")),
                                                
              fileInput("annotation_file", 
                        label = "Select file", 
                        buttonLabel = "Browse", 
                        multiple = FALSE,
                        placeholder = "No annotation file selected",
                        accept = c("application/vnd.ms-excel",
                                   "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                                   ".xlsx",
                                   ".xls")),
                                                
              helpText("For more information on formatting annotation files and accepted 
                 file types see the 'About' tab."),
              
              hr(),
              
              h6("Example Data",
                 bslib::popover(
                   bsicons::bs_icon("info-circle"),
                   "The example data contains 4 datasets - proteomics, phosphoproteomics,
                   and metabolomics positive & negative - each with treatment (Glucose)
                   and control (NoGlucose) groups. There are 5-6 biological replicates of
                   each group depending on the dataset. The data is pulled from the original
                   Omics Notebook documentation at https://github.com/cnsb-boston/Omics_Notebook_Docs."
                 )),
              
              checkboxInput("use_example_data",
                            label = "Load example data"),
              
              tags$div(
                "OR",
                tags$a(href = "https://drive.google.com/drive/folders/1lyzmIhorrZy_CKuxabi1Bv1cLHIblJhk?usp=sharing", 
                       "Click here"),
                " to download example data and annotation files."
              ),
              
              uiOutput("select_group"),
              
              hr(),
              
              h6("Species",
                 bslib::tooltip(bsicons::bs_icon("info-circle"),
                                "This will be used for enrichment analysis.",
                                placement = "right")),
              
              selectInput("species", 
                          label = "Select species",
                          c("Human (9606)" = "human",
                            "Mouse (10090)" = "mouse",
                            "Yeast (559292)" = "yeast",
                            "Zebrafish (7955)" = "zebrafish",
                            "C. elegans (6239)" = "celegans",
                            "Fruit Fly (7227)" = "fruitfly",
                            "Rat (10116)" = "rat")),
              
              hr(),
              
              h5("Data Processing Options"),
              
              br(),
              
              h6("Normalization"),
              
              selectInput("norm_eset",
                          label = "Normalization method",
                          choices = c("Median" = "median",
                                      "Median Absolute Deviation (MAD)" = "MAD",
                                      "Quantile" = "quantile",
                                      "VSN" = "vsn",
                                      "Loess" = "loess",
                                      "Median MAD" = "medianMAD",
                                      "Z-Transform" = "z transform",
                                      "Internal Reference Scaling" = "IRS",
                                      "None" = "none")),
              
              uiOutput("IRS_column"),
              
              uiOutput("IRS_explanation"),
              
              checkboxInput("batch_norm",
                            label = "Normalize in batches"),
              
              uiOutput("batch_column"),
              
              helpText("If batch normalization is selected normalization will be applied
                       to each batch separately, rather than to the dataset as a whole."),
              
              hr(),
              
              h6("Data Transformation"),
              
              checkboxInput("log_transform",
                            label = "Log-transform",
                            value = TRUE),
              
              hr(),
              
              h6("Imputation",
                 bslib::tooltip(bsicons::bs_icon("info-circle"),
                                "The `mice` R package's lasso linear regression algorithm 
                                is used for imputation.",
                                placement = "right")),
              
              checkboxInput("impute_missing",
                            label = "Impute missing values",
                            value = FALSE),
              
              helpText("Missing value imputation can make initial data normalization
                       and cleaning take much longer. Please keep this in mind when
                       waiting for initial tables and figures to render. A load bar 
                       can be seen in the bottom right-hand corner of the window."),
              
              hr(),
              
              h6("Annotation",
                 bslib::tooltip(bsicons::bs_icon("info-circle"),
                                "The UniProt database is referenced to pull
                                protein IDs if they are missing.",
                                placement = "right")),
              
              checkboxInput("uniprot_annotation",
                            label = "Uniprot annotation",
                            value = TRUE),
              
              hr(),
              
              h6("Outliers"),
              
              checkboxInput("remove_outlier",
                            label = "Remove outlier column(s)"),
              
              uiOutput("outlier_column_selection"),
              
              hr(),
              
              h6("Zero Value Filter"),
              
              numericInput("zero_cutoff",
                           label = "Choose zero cutoff",
                           value = 0.3,
                           min = 0,
                           max = 1),
              
              helpText("The cutoff value specifies the allowed proportion of zeros in a row before it is removed.
                       For example, a cutoff of 0.3 indicates that rows where more than 30% of values are zero will be removed.")
              
            ),
              
            mainPanel(
              tabsetPanel(
                
                #### FILE INPUTS ---------------------------------------------------
                tabPanel("File Inputs",
                         
                         br(),
                         
                         h3("View Data"),
                         
                         helpText("Summaries of the data that you upload will be available below. A good way to check that your data inputs are being processed correctly is 
                                  by loading the expression matrix in the ", code("View Raw Expression Matrix"), " tab."),
                         
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
                         verbatimTextOutput("summary_data"),
                         uiOutput("typeHeader"),
                         verbatimTextOutput("type"),
                         
                         hr(),
                         h4("Troubleshooting"),
                         
                         br(),
                         h6("How do I know that my data has loaded in correctly?"),
                         
                         helpText("When data has loaded in correctly: "),
                         
                         br(),
                         
                         helpText("1. An input option for the column to group by will appear above ", 
                                  code("Data Processing Options"), " on the left."),
                         
                         br(),
                         
                         helpText("2. Headers ", code("Annotation File"), " and ", code("Data File(s)"),
                                  " will load in under ", code("Data Overview"), " with a summary of the annotation 
                                  file uploaded and a count of the number of data files uploaded."),
                         
                         br(),
                         
                         helpText("3. A list of dataset names and formats (formatted 'Name;Format'), pulled from the annotation file, 
                                  will load underneath the ", code("Data File(s)"), "."),
                         
                         br(),
                         br(),
                         
                         helpText("If any of these things have not loaded in correctly it means that there has been
                                  an issue loading and formatting your data."),
                         
                         br(),
                         br(),
                         br(),
                         
                         h6("What do I do if my data hasn't loaded correctly?"),
                         
                         helpText("The most common errors that lead to errors in data upload and formatting: "),
                         
                         br(),
                         
                         helpText("- The name of the dataset file in the annotation file is incorrect. On sheet one
                                  of the annotation file under ", code("File"), "the file name must match the name of
                                  the dataset file *exactly* AND include the file type i.e. '.txt'."),
                         br(),
                         helpText("- The analysis name provided on sheet one of the annotation file does not match the
                                  name of the corresponding column on sheet two. The name of the column containing the 
                                  intensity value column names on sheet two needs to *exactly* match the analysis name
                                  for that coressponding file so that Omics Notebook knows which dataset correlates with
                                  which intensity value column names."),
                         
                         br(),
                         helpText("- No, or an incorrect, data format was selected on sheet one of the annotation file.
                                  The data format changes how the data is read, formatted, and cleaned. Selecting the 
                                  incorrect format will lead to errors."),
                         br(),
                         helpText("- The wrong annotation file was uploaded. Make sure that the annotation file 
                                  uploaded corresponds with the dataset(s) you have uploaded!")
                         ),
                
                
                #### ANNOTATION FILE -------------------------------------------
                tabPanel("Annotation Table",
                         
                         br(),
                         
                         h3("Annotation Table"),
                         
                         helpText("The table is a formatted version of the annotation information provided in the
                                  annotation file."),
                         
                         hr(),
                         br(),
                         
                         shinycssloaders::withSpinner(DT::dataTableOutput("view_annotation"), type = 8)
                         ),
                
                #### RAW DATA STATISTICAL SUMMARY -----------------------
                tabPanel("Statistical Summary",
                         br(),
                         
                         fluidRow(
                           column(8,
                                  h4("Statistical Overview"),
                                  
                                  helpText("The table provides a comprehensive statistical overview for each
                                         variable/column of the selected dataset before normalization, filtering, and imputation.
                                         The only modification made to the data at this stage is optional log scaling for improved
                                         interpretability. ", code("SD"), " stands for Standard Deviation. ",
                                           code("TRIMMED"), " stands for trimmed means, a robust estimator that is calculated after
                                         dropping a percentage of the highest and lowest values to try and avoid skewing the mean
                                         due to outliers. ", code("MAD"), " stands for the Median Absolute Deviation, a measure used
                                         to assess variation, particularly in non-symmetric data distributions. ", code("SKEW"), "
                                         indicates the degree of asymmetry in a data set. Right skew is positive and left skew is
                                         negative. ", code("KURTOSIS"), " similarly assesses skew, but by quantifying the extremity
                                         of outliers. Finally, ", code("SE"), " stands for Standard Error which estimates how representative
                                         the calculated mean is expected to be.")
                                  ),
                           column(4,
                                  bslib::card(bslib::card_title("Dataset"),
                                              bslib::card_body(fillable = FALSE,
                                                               uiOutput("statsummary_choose_dataset"))),
                                  br(),
                                  downloadButton("prenorm_table_download", 
                                                 label = "Download Table",
                                                 style = "color: #61A6F9; background-color: #D8EAFF; border-color: 
                                                       #D8EAFF; border-radius: 10px; border-width: 2px")
                                  )
                         ),
                         
                         hr(),
                         br(),
                         
                         shinycssloaders::withSpinner(DT::dataTableOutput("prenorm_describe"), type = 8),
                         
                         br(),
                         verbatimTextOutput("fdr"),
                         
                         br()
                ),
                
                #### MISSING VALUES -------------------------------------
                tabPanel("Missing Values",
                         
                         br(),
                         
                         fluidRow(
                           column(8,
                                  h4("Missing Value Visualization"),
                                  
                                  helpText("These plots are generated using the ", code("VIM"), " R package. The selected dataset
                                           is assessed pre-normalization, filtering, and imputation. The barplot
                                           on the left shows the number of missing values by feature, while the plot on the
                                           right shows the number of missing values for specific groupings of features. Red 
                                           is using to indicate missing values, while blue indicates observed values.")
                                  
                                  ),
                           column(4,
                                  uiOutput("missingvalue_choose_dataset")
                                  )
                         ),
                         
                         hr(),
                         
                         br(),
                         
                         shinycssloaders::withSpinner(plotOutput("visualizing_missing_values", height = 550), type = 8),
                         
                         br(),
                         
                         radioButtons("prenorm_missingvalue_ext",
                                      label = "Choose filetype for download",
                                      choices = c("PNG" = ".png",
                                                  "PDF" = ".pdf",
                                                  "SVG" = ".svg"),
                                      inline = TRUE),
                         
                         downloadButton("prenorm_missingvalue_download", 
                                        label = "Download Plot",
                                        style = "color: #61A6F9; background-color: #D8EAFF; border-color: 
                                                       #D8EAFF; border-radius: 10px; border-width: 2px")
                ),
                
                #### OUTLIERS -------------------------------------------
                tabPanel("Distribution",

                         br(),

                         fluidRow(
                           column(8,
                                  h4("Data Distribution"),
                                  
                                  helpText("The boxplot generated below shows the data distribution by sample for the selected
                                  dataset. This boxplot can be used to assess differences among samples in the raw data
                                  and to identify outliers - which can be removed by selecting 'Remove outlier column(s)'
                                  on the left.")
                                  ),
                           
                           column(4,
                                  uiOutput("distribution_choose_dataset")
                                  )
                         ),

                         hr(),
                         br(),

                         shinycssloaders::withSpinner(plotOutput("distribution_plot", height = 550), type = 8),

                         br(),

                         radioButtons("prenorm_distribution_ext",
                                      label = "Choose filetype for download",
                                      choices = c("PNG" = ".png",
                                                  "PDF" = ".pdf",
                                                  "SVG" = ".svg"),
                                      inline = TRUE),

                         downloadButton("prenorm_distribution_download",
                                        label = "Download Plot",
                                        style = "color: #61A6F9; background-color: #D8EAFF; border-color:
                                                #D8EAFF; border-radius: 10px; border-width: 2px")
                         ),
                
                #### EXPRESSION MATRIX -----------------------------------------
                tabPanel("Raw Expression Matrix",
                         
                         br(),
                         
                         fluidRow(
                           column(8,
                                  h3("Raw Expression Matrix"),
                                  
                                  helpText("The table produced is the result of calling the ", 
                                           code("exprs()"), " function from Biobase on the expression set object derived 
                                           from the data. The only optional modification on this data matrix is log scaling
                                           for readability, otherwise it represents the raw dataset.")
                                  ),
                           
                           column(4,
                                  uiOutput("eset_choose_dataset")
                                  )
                         ),
                         
                         hr(),
                         
                         shinycssloaders::withSpinner(DT::dataTableOutput("view_raw_eset"), type = 8),
                         
                         br(),
                         br(),
                         
                         downloadButton("eset_table_download",
                                        label = "Download Dataset",
                                        style = "color: #61A6F9; background-color: #D8EAFF; border-color: 
                                                #D8EAFF; border-radius: 10px; border-width: 2px")
                         
                )
                
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
        
        ### PROCESSED DATA OVERVIEW --------------------------------------------
        tabPanel("Processed Data Overview",
          sidebarLayout(
            sidebarPanel(
              h4("Plot Options"),
              
              h6("Dataset"),
              
              uiOutput("normplots_choose_dataset"),
              
              hr(),
              
              radioButtons("norm_plottype",
                           label = "Plot Type",
                           choices = c("Pre/Post Normalization Violin/Box-Plot" = "Boxplot",
                                       "Pre/Post Normalization Density Plot" = "Density",
                                       "QQ (quantile-quantile) Plot" = "QQ",
                                       "MA Plot" = "MA",
                                       "Relative Log Expression Plot" = "RLE")),
              
              uiOutput("ma_array_choose"),
              
              uiOutput("qqplot_column"),
              
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
              
              h6("Shapiro-Wilks Test"),
              
              uiOutput("shapiro_column"),
              
              verbatimTextOutput("shapiro_test"),
              
              hr(),
              
              h5("Download Plot"),
              
              radioButtons("normal_extension", 
                           label = "Save As:",
                           choices = c("PDF" = ".pdf", 
                                       "PNG" = ".png",
                                       "SVG" = ".svg"),
                           inline = TRUE),
              
              downloadButton("normal_download", 
                             label = "Download Plot",
                             style = "color: #61A6F9; background-color: #D8EAFF; border-color: 
                                      #D8EAFF; border-radius: 10px; border-width: 2px"),
            ),
            
            mainPanel(
              tabsetPanel(
                tabPanel("Plots",
                         br(),
                         h3("Normalization Plots"),
                         hr(),
                         
                         shinycssloaders::withSpinner(plotOutput("normal_plot", height = 800), type = 8),
                         ),
                
                tabPanel("View Processed Expression Matrix",
                         br(),
                         h3("Processed Expression Matrix"),
                         
                         helpText("The table produced is the result of calling the ", 
                                  code("exprs()"), " function from Biobase on the expression set object derived 
                                  from the data. Unlike the expression matrix under the ", code("Data Upload"),
                                  " tab, this matrix has had all filtering, cleaning, imputation, and normalization
                                  steps applied."),
                         hr(),
                         
                         shinycssloaders::withSpinner(DT::dataTableOutput("view_cleaned_eset"), type = 8),
                         
                         br(),
                         
                         downloadButton("processed_eset_download", 
                                        label = "Download Table",
                                        style = "color: #61A6F9; background-color: #D8EAFF; border-color: 
                                                #D8EAFF; border-radius: 10px; border-width: 2px"),
                         )
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
        )
      )
    ),
    
    ## BASIC ANALYSIS ##########################################################
    tabPanel("2. Basic Analysis",
      tabsetPanel(
        
        #### CORRELATION -------------------------------------------------------
        tabPanel("Correlation",
                 sidebarLayout(
                   sidebarPanel(
                     
                     h4("Options"),
                     
                     h6("Dataset"),
                     
                     uiOutput("corrplot_choose_dataset"),
                     
                     hr(),
                     
                     fluidRow(
                       column(6, colourpicker::colourInput("corr_neg_color", "Select NEG color", "royalblue", allowTransparent = T)),
                       column(6, colourpicker::colourInput("corr_pos_color", "Select POS color", "red", allowTransparent = T))
                     ),
                     
                     textInput("corr_title",
                               label = "Customize title",
                               placeholder = "Type title here"),
                     
                     numericInput("corr_variance",
                                  label = "Proportion of highest variance rows to include",
                                  value = 0.10,
                                  min = 0.01,
                                  max = 1,
                                  step = 0.05),
                     
                     checkboxInput("corr_number",
                                   label = "Add correlation values",
                                   value = TRUE),
                     
                     helpText("When this box is checked the calculated correlation values will be added to each square of the correlation plot."),
                     
                     br(),
                     br(),
                     
                     actionButton("corr_button",
                                  label = "Run Analysis",
                                  style = "color: #fff; background-color: #9BD79A; 
                                   border-color: #9BD79A; border-radius: 10px; border-width: 2px"),
                     br(),
                     br(),
                     
                     h6("This analysis may take a few minutes to load. \n"),
                     
                     hr(),
                     
                     h5("Download Plots"),
                     
                     br(),
                     
                     fluidRow(
                       column(6,
                              h6("Total Sample Plot"),
                              
                              radioButtons("corrplot_extension", 
                                           label = "Save As:",
                                           choices = c("PDF" = ".pdf", 
                                                       "PNG" = ".png",
                                                       "SVG" = ".svg"),
                                           inline = TRUE),
                              
                              downloadButton("corr_plot_download", 
                                             label = "Download Plot",
                                             style = "color: #61A6F9; background-color: #D8EAFF; border-color: 
                                      #D8EAFF; border-radius: 10px; border-width: 2px")
                              ),
                       column(6,
                              h6("Highest Variance Plot"),
                              
                              radioButtons("var_corrplot_extension", 
                                           label = "Save As:",
                                           choices = c("PDF" = ".pdf", 
                                                       "PNG" = ".png",
                                                       "SVG" = ".svg"),
                                           inline = TRUE),
                              
                              downloadButton("var_corr_plot_download", 
                                             label = "Download Plot",
                                             style = "color: #61A6F9; background-color: #D8EAFF; border-color: 
                                                     #D8EAFF; border-radius: 10px; border-width: 2px")
                              )
                     ),
                     
                     br(),
                     
                   ),
                   
                   mainPanel(
                     br(),
                     h3("Correlation"),
                     
                     helpText("The correlation plot uses a z-scaled version of the dataset which
                       has NA omission. The below correlation table can be downloaded as a CSV file."),
                     
                     hr(),
                     
                     fluidRow(
                       column(6, shinycssloaders::withSpinner(plotOutput("correlation_plot", height = 750), type = 8)),
                       column(6, shinycssloaders::withSpinner(plotOutput("variance_correlation_plot", height = 750), type = 8))
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
        
        #### MD PLOT -----------------------------------------------------------
        tabPanel("MD Plot",
                 sidebarLayout(
                   sidebarPanel(
                     h4("Analysis Options"),
                     
                     h6("Dataset"),
                     
                     uiOutput("MDplot_choose_dataset"),
                     
                     hr(),
                     
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
                                  style = "color: #fff; background-color: #9BD79A; 
                                   border-color: #9BD79A; border-radius: 10px; border-width: 2px"),
                     
                     br(),
                     br(),
                     
                     h6("This analysis may take a few minutes to load. \n"),
                     
                     hr(),
                     
                     h5("Download Plot"),
                     
                     radioButtons("md_extension", 
                                  label = "Save As:",
                                  choices = c("PDF" = ".pdf", 
                                              "PNG" = ".png",
                                              "SVG" = ".svg"),
                                  inline = TRUE),
                     
                     downloadButton("md_download", 
                                    label = "Download MD Plot",
                                    style = "color: #61A6F9; background-color: #D8EAFF; border-color: 
                                     #D8EAFF; border-radius: 10px; border-width: 2px")
                   ),
                   
                   mainPanel(
                     
                     br(),
                     h3("Mean Difference Plot"),
                     
                     helpText("The MD plot requires both ", em("Annotation "), "and ", em("Data "), 
                              "file inputs. It also requires selection of a ", code("Group"), " column
                       in the Data Upload tab."),
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
        
        #### UMAP --------------------------------------------------------------
        tabPanel("UMAP",
                 sidebarLayout(
                   sidebarPanel(
                     h4("UMAP Analysis Options"),
                     
                     h6("Dataset"),
                     
                     uiOutput("UMAPplot_choose_dataset"),
                     
                     hr(),
                     
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
                                  style = "color: #fff; background-color: #9BD79A; border-color: 
                                   #9BD79A; border-radius: 10px; border-width: 2px"),
                     br(),
                     br(),
                     
                     h6("This analysis may take a few minutes to load. \n"),
                     
                     hr(),
                     
                     h5("Download Plot"),
                     
                     radioButtons("UMAP_extension", 
                                  label = "Save As:",
                                  choices = c("PDF" = ".pdf", 
                                              "PNG" = ".png",
                                              "SVG" = ".svg"),
                                  inline = TRUE),
                     
                     downloadButton("UMAP_download", 
                                    label = "Download UMAP Plot",
                                    style = "color: #61A6F9; background-color: #D8EAFF; border-color: 
                                     #D8EAFF; border-radius: 10px; border-width: 2px")
                   ),
                   
                   mainPanel(
                     br(),
                     h3("UMAP Plot"),
                     
                     helpText("UMAP requires both ", em("Annotation "), "and ", em("Data "), 
                              "file inputs. It also requires selection of a ", code("Group"), 
                              " column in the Data Upload tab."),
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
        
        #### HEATMAP -----------------------------------------------------------
        tabPanel("Heatmap",
                 sidebarLayout(
                   sidebarPanel(
                     h4("Heatmap Analysis Options"),
                     
                     h6("Dataset"),
                     
                     uiOutput("heatmap_choose_dataset"),
                     
                     hr(),
                     
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
                                   "Cluster rows",
                                   value = TRUE),
                     
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
                                   "Subset by feature" = "subset_features"),
                                 selected = "subset_variable"),
                     
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
                                  choices = c("PDF" = ".pdf", 
                                              "PNG" = ".png",
                                              "SVG" = ".svg"),
                                  inline = TRUE),
                     
                     downloadButton("heatmap_download", 
                                    label = "Download Heatmap Plot",
                                    style = "color: #61A6F9; background-color: #D8EAFF; border-color: 
                                    #D8EAFF; border-radius: 10px; border-width: 2px"),
                     
                     br(),
                     br(),
                     
                     h6("Please allow a few moments for downloading.")
                     
                   ),
                   
                   mainPanel(
                     br(),
                     h3("Heatmap"),
                     
                     helpText("Heatmap analysis requires both ", em("Annotation "), "and ", em("Data "), 
                              "file inputs. It also requires selection of a ", code("Group"),
                              " column in the Data Upload tab, and takes into 
                      account the selection of ", em("Log-Transform"), " for the expression set."),
                     
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
        
        #### PCA ---------------------------------------------------------------
        tabPanel("PCA",
          sidebarLayout(
            sidebarPanel(
              h4("Analysis Options"),
              
              h6("Dataset"),
              
              uiOutput("PCplot_choose_dataset"),
              
              hr(),
              
              h6("Color Palette"),
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
              
              hr(),
              h6("Title"),
              textInput("PCA_title",
                        label = "Customize Plot Title",
                        placeholder = "Type title here"),
              
              hr(),
              h6("Axes"),
              fluidRow(
                column(6, uiOutput("PC_xaxis")),
                column(6, uiOutput("PC_yaxis"))
              ),
               
              hr(),
              h6("Other Visualization Options"),
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
              
              actionButton("PCA_button", "Plot PCA",
                          class = "btn-block",
                          style = "color: #fff; background-color: #9BD79A; border-color: 
                                  #9BD79A; border-radius: 10px; border-width: 2px"),
               
              br(),
              br(),
               
              h6("This analysis may take a few minutes to load. \n"),
               
              hr(),
               
              h5("Download Plot"),
              
              radioButtons("PCA_extension", 
                          label = "Save As:",
                          choices = c("PDF" = ".pdf", 
                                      "PNG" = ".png",
                                      "SVG" = ".svg"),
                          inline = TRUE),
               
              downloadButton("PCA_download", 
                            label = "Download PCA Plot",
                            style = "color: #61A6F9; background-color: #D8EAFF; border-color: 
                                    #D8EAFF; border-radius: 10px; border-width: 2px"),
            ),
                                   
            mainPanel(
              br(),
              h3("Principle Component Analysis"),
              
              fluidRow(
                column(6, 
                       helpText("PCA requires both ", em("Annotation "), "and ", em("Data "), 
                                "file inputs. It also requires selection of a ", code("Group"), 
                                " column in the Data Upload tab.")
                       ),
                
                column(3, 
                       uiOutput("PC_loadings_column")
                       ),
                
                column(3, downloadButton("PC_loadings_download", 
                                         label = "Download PC Loadings Table",
                                         style = "color: #61A6F9; background-color: #D8EAFF; border-color: 
                                                 #D8EAFF; border-radius: 10px; border-width: 2px")
                       )
              ),
              
              
              fluidRow(
                column(6, 
                       shinycssloaders::withSpinner(plotOutput("PC_variance_plot"), type = 8),
                       helpText(textOutput("top_PC_print"))
                       ),
                
                column(6, 
                       shinycssloaders::withSpinner(plotOutput("PC_loadings_plot"), type = 8),
                       helpText("The PC Loadings plot shows which rows most highly impact the principle components.")
                       )
              ),
              
              hr(),
               
              shinycssloaders::withSpinner(plotOutput("PC_plot", height = 600), type = 8)
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
      )
    ),
    
    ## SINGLE OMIC ANALYSIS ####################################################
    
    tabPanel("3. Single Omic Analysis",
             tabsetPanel(
               ### LINEAR MODELING ---------------------------------------------
               tabPanel("Linear Modeling",
                        sidebarLayout(
                          sidebarPanel(width = 3,
                                       h4("Model Options"),
                                       
                                       helpText("There are many optional additions to linear modeling.
                                     At a minimum selections must be made under ", code("Model Factors"),
                                                " and ", code("Contrasts"), ". Important covariates
                                     to add may include variables associated with batch effects
                                     such as run groupings."),
                                       
                                       hr(),
                                       
                                       h6("Dataset"),
                                       
                                       uiOutput("limma_choose_dataset"),
                                       
                                       helpText("Only one uploaded dataset can be used at a time for 
                              linear modeling analysis."),
                                       
                                       hr(),
                                       
                                       h6("Model factors"),
                                       
                                       uiOutput("linear_factors"),
                                       
                                       helpText("These options are derived from the number of unique
                              groups in the uploaded annotation file."),
                                       
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
                                       
                                       helpText("Only contrasts chosen at this step will be available 
                              for volcano plots, enrichment, and other differential 
                              analysis."),
                                       
                                       hr(),
                                       
                                       actionButton("limma_button",
                                                    label = "Generate Fit",
                                                    style = "color: #fff; background-color: #9BD79A; 
                              border-color: #9BD79A; border-radius: 10px; border-width: 2px"),
                                       
                                       br(),
                                       hr(),
                                       
                                       helpText("See ", code("Linear Modeling"), " under the", code("About"), "tab for
                                     more information on linear modeling options and an
                                     overview of the process."),
                                       
                                       hr(),
                                       h5("Download Outputs"),
                                       
                                       helpText("Click ", code("Download Results"), " below to download the ", em("Model Results Table"),
                                                ", ", code("Download Top Table"), " to download the ", em("Limma Top Table"), ", or ",
                                                code("Download P-Val Plot"), " to download the ", em("P-Value distribution histograms.")),
                                       
                                       br(),
                                       br(),
                                       
                                       downloadButton("linear_results_table_download",
                                                      label = "Download Results",
                                                      style = "color: #61A6F9; background-color: #D8EAFF; border-color: 
                                         #D8EAFF; border-radius: 10px; border-width: 2px"),
                                       
                                       br(),
                                       br(),
                                       
                                       downloadButton("linear_table_download",
                                                      label = "Download top table",
                                                      style = "color: #61A6F9; background-color: #D8EAFF; border-color: 
                                         #D8EAFF; border-radius: 10px; border-width: 2px"),
                                       
                                       br(),
                                       br(),
                                       
                                       downloadButton("pvaldist_plot_download",
                                                      label = "Download P-val Plot",
                                                      style = "color: #61A6F9; background-color: #D8EAFF; border-color: 
                                         #D8EAFF; border-radius: 10px; border-width: 2px"),
                                       
                                       br(),
                                       br(),
                                       
                                       radioButtons("pvaldist_extension", 
                                                    label = "Save Plot As:",
                                                    choices = c("PDF" = ".pdf", 
                                                                "PNG" = ".png",
                                                                "SVG" = ".svg"),
                                                    inline = TRUE)
                                       
                          ),
                          
                          mainPanel(width = 9,
                                    tabsetPanel(
                                      #### FIT SUMMARY ----------------------------------------
                                      tabPanel("Fit Summary",
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
                                                       code("P.VALUE"), " is the raw p-value; ",
                                                       code("ADJ.P.VALUE"), " is the adjusted p-value or q-value; ",
                                                       code("B"), " is the log-odds that the gene is differentially expressed."),
                                               
                                               br(),
                                               br(),
                                               fluidRow(
                                                 column(8, uiOutput("linear_model_results_header")),
                                                 column(2, uiOutput("linear_model_results_logfc_cutoff")),
                                                 column(2, uiOutput("linear_model_results_pval_cutoff"))
                                               ),
                                               
                                               br(),
                                               
                                               tableOutput("linear_model_results"),
                                               
                                               uiOutput("linear_model_equation_header"),
                                               
                                               verbatimTextOutput("model_equation"),
                                               
                                               hr(),
                                               br(),
                                               
                                               fluidRow(
                                                 column(4, uiOutput("limma_toptable_header")),
                                                 column(4, uiOutput("coef_options")),
                                               ),
                                               
                                               shinycssloaders::withSpinner(DT::dataTableOutput("linear_model_table"), type = 8),
                                               br(),
                                               shinycssloaders::withSpinner(plotOutput("pval_distribution_plot"), type = 8)
                                                
                                      ),
                                      
                                      #### BATCH PCA ---------------------------
                                      tabPanel("Batch Correction PCA",
                                               
                                               fluidRow(
                                                 column(9,
                                                        br(),
                                                        h3("Batch Corrected PCA"),
                                                        
                                                        helpText("The batch corrected PCA plot is representative of how the addition of the
                                                                 selected batch column corrects for batch effect when added as a covariate
                                                                 to the linear model. It models the effect of adding the chosen batch factor 
                                                                 to a basic linear model design using the ", code("removeBatchEffect"), " function
                                                                 from the ", code("limma"), " package."),
                                                 ),
                                                 column(3,
                                                        h5("Download Plot"),
                                                        
                                                        radioButtons("batch_PCA_extension", 
                                                                     label = "Save As:",
                                                                     choices = c("PDF" = ".pdf", 
                                                                                 "PNG" = ".png",
                                                                                 "SVG" = ".svg"),
                                                                     inline = TRUE),
                                                        
                                                        downloadButton("batch_PCA_download", 
                                                                       label = "Download PCA Plot",
                                                                       style = "color: #61A6F9; background-color: #D8EAFF; border-color: 
                                                                                #D8EAFF; border-radius: 10px; border-width: 2px")
                                                 )
                                               ),
                                               
                                               hr(),
                                               
                                               fluidRow(
                                                 column(9,
                                                        shinycssloaders::withSpinner(plotOutput("PCA_batch", height = 800), type = 8)
                                                 ),
                                                 
                                                 column(3,
                                                        actionButton("batch_PCA_button",
                                                                     label = "Generate Plot",
                                                                     style = "color: #fff; background-color: #9BD79A; 
                                                                              border-color: #9BD79A; border-radius: 10px; border-width: 2px"),
                                                        
                                                        hr(),
                                                        h6("Plot Options"),
                                                        
                                                        uiOutput("batch_PCA_column"),
                                                        
                                                        selectInput("batch_PCA_color", 
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
                                                       
                                                        checkboxInput("batch_PCA_include_densities",
                                                                      label = "Include PC density curves along axes"),
                                                       
                                                        checkboxInput("batch_PCA_ellipse",
                                                                      label = "Add ellipses"),
                                                       
                                                        checkboxInput("batch_PCA_labels",
                                                                      label = "Label points",
                                                                      value = TRUE),
                                                       
                                                        checkboxInput("batch_PCA_shapes",
                                                                      label = "Make different groups different shapes"),
                                                        hr()
                                                 )
                                               )
                                        
                                      ),
                                      #### VOLCANO PLOT ---------------------------------------
                                      tabPanel("Volcano Plot",
                                               
                                               br(),
                                               fluidRow(
                                                 column(9,
                                                        h3("Volcano Plot"),
                                                        
                                                        helpText("The Volcano Plot is generated based on the input in the ", em("Linear Model "), "tab. It requires ", em("Annotation "), "and ", em("Data "), 
                                                                 "file inputs as well as selection of ", em("Data Type "), "in the ", em("Data Upload "), "tab. Make customizations and select ", code("Generate Plot"), " on the right
                                                  to generate the volcano plot. If changes are made to the linear model (on the left) it will change the volcano plot options and output, and the plot will
                                                  need to be re-generated."),
                                                 ),
                                                 column(3,
                                                        h5("Download Plot"),
                                                        
                                                        radioButtons("volcano_extension", 
                                                                     label = "Save As:",
                                                                     choices = c("PDF" = ".pdf", 
                                                                                 "PNG" = ".png",
                                                                                 "SVG" = ".svg"),
                                                                     inline = TRUE),
                                                        
                                                        downloadButton("volcano_download", 
                                                                       label = "Download Volcano Plot",
                                                                       style = "color: #61A6F9; background-color: #D8EAFF; border-color: 
                                                                    #D8EAFF; border-radius: 10px; border-width: 2px")
                                                 )
                                               ),
                                               
                                               hr(),
                                               
                                               fluidRow(
                                                 column(9,
                                                        shinycssloaders::withSpinner(plotOutput("volcano_plot", height = 800), type = 8)
                                                 ),
                                                 
                                                 column(3,
                                                        actionButton("volcano_button",
                                                                     label = "Generate Plot",
                                                                     style = "color: #fff; background-color: #9BD79A; 
                                                              border-color: #9BD79A; border-radius: 10px; border-width: 2px"),
                                                        
                                                        hr(),
                                                        h6("Plot Options"),
                                                        
                                                        fluidRow(
                                                          column(6, colourpicker::colourInput("volcano_down_color", "DOWN color", "royalblue", allowTransparent = T)),
                                                          column(6, colourpicker::colourInput("volcano_up_color", "UP color", "red", allowTransparent = T))
                                                        ),
                                                        
                                                        textInput("volcano_title",
                                                                  label = "Customize Plot Title",
                                                                  placeholder = "Type title here"),
                                                        
                                                        uiOutput("volcano_coef_options"),
                                                        
                                                        checkboxInput("volcano_labels",
                                                                      label = "Add labels to genes with highest/lowest logFC",
                                                                      value = TRUE),
                                                        
                                                        uiOutput("volcano_label_number"),
                                                        
                                                        checkboxInput("volcano_label_specific",
                                                                      label = "Add label to specific gene(s)"),
                                                        
                                                        uiOutput("volcano_label_specific_gene"),
                                                        
                                                        radioButtons("volcano_yaxis",
                                                                     label = "Choose Y-Axis Variable",
                                                                     choices = c("Adjusted P-Value" = "adj.P.Val",
                                                                                 "P-Value" = "P.Value")),
                                                        
                                                        fluidRow(
                                                          column(6,
                                                                 numericInput("volcano_top_fc",
                                                                              label = "LogFC Cut Off",
                                                                              value = 1,
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
                                                        
                                                        hr()
                                                 )
                                               )
                                      ),
                                      
                                      #### INTERACTIVE VOLCANO PLOTS --------------------------
                                      tabPanel("Interactive Volcano Plots",
                                               br(),
                                               
                                               fluidRow(
                                                 column(6,
                                                        h3("Interactive Volcano Plots"),
                                                        
                                                        helpText("The Volcano Plot is generated based on the input in the ", em("Linear Model "), "tab. It requires ", em("Annotation "), "and ", em("Data "), 
                                                                 "file inputs as well as selection of ", em("Data Type "), "in the ", em("Data Upload "), "tab."),
                                                        
                                                        br(),
                                                        br(),
                                                        
                                                        helpText("The plotly graph object generated below allows for interaction with the generated volcano plot via hovering over data points, clicking & dragging,
                                                  or use of the menu in the top right. PNG images can be downloaded of the graph, but will not include labels or information from the hovertabs. Once
                                                  the linear model fit has been run and customization options selected, click ", code("Generate Plotly"), " to view the volcano plotly."),
                                                        
                                                 ),
                                                 
                                                 column(6,
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
                                                                     style = "color: #fff; background-color: #9BD79A; 
                                                              border-color: #9BD79A; border-radius: 10px; border-width: 2px"),
                                                        
                                                        br(),
                                                        br(),
                                                        
                                                        h6("Please allow a few moments for loading.")
                                                 )
                                               ),
                                               
                                               hr(),
                                               
                                               fluidRow(
                                                 column(9,
                                                        shinycssloaders::withSpinner(plotly::plotlyOutput("volcano_plot_interactive", height = "700px",), type = 8)
                                                 ),
                                                 column(3,
                                                        actionButton("int_volcano_button",
                                                                     label = "Generate Plotly",
                                                                     style = "color: #fff; background-color: #9BD79A; 
                                                              border-color: #9BD79A; border-radius: 10px; border-width: 2px"),
                                                        
                                                        hr(),
                                                        
                                                        h6("Plot Options"),
                                                        
                                                        fluidRow(
                                                          column(6, colourpicker::colourInput("int_volcano_down_color", "DOWN color", "royalblue", allowTransparent = T)),
                                                          column(6, colourpicker::colourInput("int_volcano_up_color", "UP color", "red", allowTransparent = T))
                                                        ),
                                                        
                                                        uiOutput("int_volcano_coef_options"),
                                                        
                                                        radioButtons("int_volcano_yaxis",
                                                                     label = "Choose Y-Axis Variable",
                                                                     choices = c("Adjusted P-Value" = "adj.P.Val",
                                                                                 "P-Value" = "P.Value")),
                                                        fluidRow(
                                                          column(6,
                                                                 numericInput("int_volcano_top_fc",
                                                                              label = "LogFC Cut Off",
                                                                              value = 1,
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
                                                        
                                                 )
                                               )
                                      ),
                                      
                                      #### DIFFERENTIAL HEATMAP -------------------------------
                                      tabPanel("Differential Heatmap",
                                               br(),
                                               
                                               fluidRow(
                                                 column(9,
                                                        h3("Heatmap"),
                                                        
                                                        helpText("Heatmap analysis requires both ", em("Annotation "), "and ", em("Data "), 
                                                                 "file inputs. It also requires selection of ", em("Data Format "),
                                                                 "and ", em("Data Type "), "in the Data Upload tab, and takes into 
                                                  account the selection of ", em("Log-Transform"), " for the expression set. "),
                                                        
                                                        br(),
                                                        br(),
                                                        
                                                        helpText("This heatmap uses the data output from the linear model in the top table.
                                                  It uses the rows with the largest logFC absolute values."),
                                                 ),
                                                 
                                                 column(3,
                                                        h5("Download Plot"),
                                                        
                                                        radioButtons("dif_heatmap_extension", 
                                                                     label = "Save As:",
                                                                     choices = c("PDF" = ".pdf", 
                                                                                 "PNG" = ".png",
                                                                                 "SVG" = ".svg"),
                                                                     inline = TRUE),
                                                        
                                                        downloadButton("dif_heatmap_download", 
                                                                       label = "Download Heatmap Plot",
                                                                       style = "color: #61A6F9; background-color: #D8EAFF; border-color: 
                                                                #D8EAFF; border-radius: 10px; border-width: 2px")
                                                 )
                                               ),
                                               
                                               hr(),
                                               
                                               fluidRow(
                                                 column(9,
                                                        shinycssloaders::withSpinner(plotOutput("dif_heatmap_plot"), type = 8)
                                                 ),
                                                 column(3,
                                                        actionButton("dif_heatmap_button",
                                                                     label = "Generate Heatmap",
                                                                     style = "color: #fff; background-color: #9BD79A; 
                                                              border-color: #9BD79A; border-radius: 10px; border-width: 2px"),
                                                        hr(),
                                                        
                                                        h6("Plot Options"),
                                                        
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
                                                                      "Cluster rows",
                                                                      value = TRUE),
                                                        
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
                                                        
                                                        hr()
                                                        
                                                 )
                                               )
                                      ),
                                      
                                      ### LIMMA ENRICHMENT ----------------------
                                      tabPanel("Enrichment",
                                               
                                               br(),
                                               
                                               h4("Linear Model Enrichment"),
                                               
                                               fluidRow(
                                                 column(9, helpText("The enrichment function uses the output of the data modeling linear model to
                                                                    create a ranked gene list. It will not function if the linear model has not been run.")),
                                                 column(3, actionButton("enrichment_button",
                                                                        label = "Run Enrichment Analysis",
                                                                        style = "color: #fff; background-color: #9BD79A; 
                                                                                border-color: #9BD79A; border-radius: 10px; border-width: 2px"))
                                                 
                                               ),
                                               
                                               br(),
                                               br(),
                                               
                                               fluidRow(
                                                 column(3, 
                                                        bslib::card(bslib::card_title("Database"),
                                                                    bslib::card_body(fillable = FALSE,
                                                                                     uiOutput("gmts")))
                                                        
                                                        ),
                                                 column(3, 
                                                        bslib::card(bslib::card_title("Algorithm"),
                                                                    bslib::card_body(fillable = FALSE,
                                                                                     radioButtons("enrichment_calculation",
                                                                                                  label = "Enrichment algorithm to perform",
                                                                                                  choices = c("Ranked Enrichment Analysis" = "ranked",
                                                                                                              "Unranked Over-Representation Analysis (ORA)" = "unranked"))))
                                                        ),
                                                 column(3, 
                                                        bslib::card(bslib::card_title("Contrast"),
                                                                    bslib::card_body(fillable = FALSE,
                                                                                     uiOutput("enrichment_coef_options")))
                                                        
                                                        ),
                                                 column(3, 
                                                        bslib::card(bslib::card_title("Threshold"),
                                                                    bslib::card_body(fillable = FALSE,
                                                                                     numericInput("enrichment_pval",
                                                                                                  "Adjusted p-value cutoff",
                                                                                                  value = 0.05,
                                                                                                  min = 0,
                                                                                                  max = 1,
                                                                                                  step = 0.05)))
                                                        
                                                        )
                                               ),
                                               
                                               tabsetPanel(
                                                 
                                                 #### OVERVIEW ------------------------------------
                                                 tabPanel("Overview",
                                                          
                                                          br(),
                                                          
                                                          h4("Enriched Pathways"),
                                                           
                                                          helpText("The ", code("clusterProfiler"), " package is used for pathway enrichment. 
                                                                    The table generated below when enrichment analysis is run contains the names of 
                                                                    the enriched pathways based on the species and database selections. For both GSEA and 
                                                                    general enrichment types there is an additional column for the adjusted p-value associated
                                                                    with each pathway. GSEA enrichment will have another column with normalized enrichment 
                                                                    score values, while general enrichment will have a column with the gene count ratio."),
                                                           
                                                          br(),
                                                          br(),
                                                           
                                                          helpText("The plots in the additional tabs (except the heatmap) are generated using the ", code("enrichplot"), "package,
                                                                    available on Bioconductor and maintained by Guangchuang Yu."),
                                                           
                                                          br(),
                                                          br(),
                                                          
                                                          downloadButton("enriched_table_download",
                                                                          label = "Download Enrichment Table",
                                                                          style = "color: #61A6F9; background-color: #D8EAFF; border-color: 
                                                                         #D8EAFF; border-radius: 10px; border-width: 2px"),
                                                                   
                                                          
                                                          hr(),
                                                          
                                                          # verbatimTextOutput("enrichment_no_results"),
                                                          shinycssloaders::withSpinner(DT::dataTableOutput("enriched_table"), type = 8)
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
                                                                   shinycssloaders::withSpinner(plotOutput("enriched_dotplot", height = 1000), type = 8)
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
                                                                   
                                                                   radioButtons("dotplot_extension", 
                                                                                label = "Save As:",
                                                                                choices = c("PDF" = ".pdf", 
                                                                                            "PNG" = ".png",
                                                                                            "SVG" = ".svg"),
                                                                                inline = TRUE),
                                                                   
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
                                                                   shinycssloaders::withSpinner(plotOutput("enriched_cnetplot", height = 1000), type = 8)
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
                                                                   
                                                                   radioButtons("cnet_extension", 
                                                                                label = "Save As:",
                                                                                choices = c("PDF" = ".pdf", 
                                                                                            "PNG" = ".png",
                                                                                            "SVG" = ".svg"),
                                                                                inline = TRUE),
                                                                   
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
                                                                   shinycssloaders::withSpinner(plotOutput("enriched_emapplot", height = 1000), type = 8)
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
                                                                   
                                                                   radioButtons("emap_extension", 
                                                                                label = "Save As:",
                                                                                choices = c("PDF" = ".pdf", 
                                                                                            "PNG" = ".png",
                                                                                            "SVG" = ".svg"),
                                                                                inline = TRUE),
                                                                   
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
                                                                   shinycssloaders::withSpinner(plotOutput("enriched_upsetplot", height = 1000), type = 8)
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
                                                                   
                                                                   radioButtons("upset_extension", 
                                                                                label = "Save As:",
                                                                                choices = c("PDF" = ".pdf", 
                                                                                            "PNG" = ".png",
                                                                                            "SVG" = ".svg"),
                                                                                inline = TRUE),
                                                                   
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
                                                                   shinycssloaders::withSpinner(plotOutput("enriched_treeplot", height = 1000), type = 8)
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
                                                                   
                                                                   radioButtons("tree_extension", 
                                                                                label = "Save As:",
                                                                                choices = c("PDF" = ".pdf", 
                                                                                            "PNG" = ".png",
                                                                                            "SVG" = ".svg"),
                                                                                inline = TRUE),
                                                                   
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
                                                          shinycssloaders::withSpinner(plotly::plotlyOutput("enriched_heatplot"), type = 8)
                                                 ),
                                                 #### KEGG PATHWAY DIAGRAM -------------------------
                                                 
                                                 tabPanel("KEGG",
                                                          br(),
                                                          h4("KEGG Pathway Enrichment Diagram"),
                                                          helpText("Certain databases have associated visualization options. If you selected ",
                                                                   code("KEGG"), " enrichment pathway an additional plot can be rendered below."),
                                                          
                                                          hr(),
                                                          br(),
                                                          fluidRow(
                                                            column(9,
                                                                   verbatimTextOutput("notKEGG_message"),
                                                                   shinycssloaders::withSpinner(plotOutput("enriched_KEGGplot", height = "auto"), type = 8)
                                                            ),
                                                            column(3,
                                                                   uiOutput("KEGG_pathway_options")
                                                            )
                                                          ),
                                                 )
                                               )
                                               
                                               )
                        ),
                )
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
    
    ## MULTI OMIC INTEGRATION ##################################################
    
    tabPanel("4. Multi Omic Integration",
             sidebarLayout(
               sidebarPanel(width = 3,
                            h4("Options"),
                            
                            helpText("S-score integration begins with fitting linear models 
                                     for all selected datasets, therefore similar selections
                                     to those found in the ", code("Single Omic Analysis"),
                                     " linear modeling tab must be made."),
                            
                            hr(),
                            
                            h6("Datasets",
                               bslib::tooltip(bsicons::bs_icon("info-circle"),
                                              "At least one dataset MUST have a UniProt ID column.",
                                              placement = "right")),
                            
                            uiOutput("sscore_datasets"),
                            helpText("Choose at least two omics datasets to integrate. 
                                  Please note that the more datasets you select the more
                                  processing time as every potential row overlap is calculated
                                  between each dataset."),
                            
                            hr(),
                            
                            h6("Model factors"),
                            
                            uiOutput("sscore_linear_factors"),
                            
                            helpText("These options are derived from the number of unique
                                  groups in the uploaded annotation file."),
                            
                            hr(),
                            
                            h6("Covariates"),
                            
                            checkboxInput("sscore_add_covariate",
                                          label = "Add covariate(s) to the model"),
                            
                            uiOutput("sscore_choose_covariate"),
                            
                            helpText("Choose up to 3 columns to include as covariates in the model. This will not change the contrast choices."),
                            
                            hr(),
                            
                            h6("Time Series"),
                            
                            checkboxInput("sscore_include_timeseries",
                                          label = "Include Time as an interaction term in the model"),
                            
                            uiOutput("sscore_choose_time"),
                            
                            uiOutput("sscore_time_type"),
                            
                            uiOutput("sscore_time_points"),
                            
                            uiOutput("sscore_time_df_explain"),
                            
                            hr(),
                            
                            h6("Contrasts"),
                            
                            uiOutput("sscore_choose_contrast"),
                            
                            helpText("Only contrasts chosen at this step will be available
                              for volcano plots, enrichment, and other differential
                              analysis."),
                            
                            hr(),
                            
                            actionButton("sscore_button",
                                         label = "Run Integration",
                                         style = "color: #fff; background-color: #9BD79A; 
                                              border-color: #9BD79A; border-radius: 10px; border-width: 2px"),
                            
                            br(),
                            br(),
                            
                            h6("Please allow a few minutes for s-score integration to run."),
                            
                            br(),
                            
                            helpText("When ", code("Run Integration"), " is clicked integration
                                     will be run for all contrasts selected to include in the model.
                                     For this reason, the more datasets and contrasts that are included
                                     the longer the integration run will take.")
                            
               ),
               
               mainPanel(width = 9,
                        tabsetPanel(
                          ### SUMMARY TABLE ------------------------------------
                          tabPanel("Summary Table",
                                   br(),
                                   h3("Multi-Omic S-Score Integration"),
                                   
                                   br(),
                                   
                                   fluidRow(
                                     column(9,
                                            helpText("This method is meant to integrate across omics
                                          datasets to identify proteins that are abundantly
                                          differentially expressed by a combination of two
                                          or more datasets. Differential expression is determined
                                          in the context of each dataset and datasets are weighted
                                          based on size. Then the 'S-Score', or abundance score,
                                          is calculated using the z-score of the log fold change 
                                          of each protein, weighted by dataset."),
                                            
                                            br(),
                                            br(),
                                            
                                            helpText("Benefits of S-Score integration include integrated normalization,
                                          flexible weighting, and the relative linear relationship
                                          of scores for most values of S. For more information on this
                                          method please see the referenced publication under the ", code("About"),
                                                     "/", code("Documentation"), " tab."),
                                            
                                            br(),
                                     ),
                                     
                                     column(3,
                                            h6("Visualization",
                                               bslib::tooltip(bsicons::bs_icon("info-circle"),
                                                              "Choose a contrast to change the table and volcano plot visualizations.",
                                                              placement = "bottom")),
                                            
                                            uiOutput("sscore_coef_options")
                                     )
                                   ),
                                   
                                   br(),
                                   hr(),
                                     
                                   fluidRow(
                                     column(9, 
                                            h4("S-Score Integration Summary Table")
                                            ),
                                     column(3, 
                                            downloadButton("sscore_table_download", 
                                                           label = "Download Table",
                                                           style = "color: #61A6F9; background-color: #D8EAFF; border-color: 
                                                                   #D8EAFF; border-radius: 10px; border-width: 2px")
                                            )
                                   ),
                                   br(),
                                   shinycssloaders::withSpinner(DT::dataTableOutput("sscore_integration_summary"), type = 8)
                                   ),
                          
                          tabPanel("Volcano Plot",
                                   fluidRow(
                                     column(9, 
                                            br(),
                                            h4("S-Score Volcano Plot"),
                                            helpText("The volcano plot is generated using the calculated s-score values
                                                     which can be seen in the summary table plot.")
                                            ),
                                     column(3, 
                                            hr(),
                                            h5("Download Plot"),
                                            radioButtons("sscore_volcano_extension", 
                                                         label = "Save As:",
                                                         choices = c("PDF" = ".pdf", 
                                                                     "PNG" = ".png",
                                                                     "SVG" = ".svg"),
                                                         inline = TRUE),
                                            
                                            downloadButton("sscore_volcano_download", 
                                                           label = "Download Volcano Plot",
                                                           style = "color: #61A6F9; background-color: #D8EAFF; border-color: 
                                                                   #D8EAFF; border-radius: 10px; border-width: 2px")
                                            )
                                   ),
                                   
                                   hr(),
                                   
                                   fluidRow(
                                     column(9, shinycssloaders::withSpinner(plotOutput("sscore_volcano", height = 800), type = 8)),
                                     column(3, 
                                            actionButton("sscore_volcano_button",
                                                            label = "Generate Plot",
                                                            style = "color: #fff; background-color: #9BD79A; 
                                                                    border-color: #9BD79A; border-radius: 10px; border-width: 2px"),
                                            hr(),
                                            h6("Plot Options"),
                                            
                                            fluidRow(
                                              column(6, colourpicker::colourInput("sscore_volcano_down_color", "DOWN color", "royalblue", allowTransparent = T)),
                                              column(6, colourpicker::colourInput("sscore_volcano_up_color", "UP color", "red", allowTransparent = T))
                                            ),
                                            
                                            checkboxInput("sscore_volcano_label",
                                                          label = "Add labels to points with highest log adjusted p-values",
                                                          value = TRUE),
                                            
                                            numericInput("sscore_volcano_label_num",
                                                         label = "Number of up/down regulated points to label",
                                                         value = 10,
                                                         min = 1),
                                            
                                            numericInput("sscore_volcano_pval",
                                                         label = "P-Value Cut Off",
                                                         value = 0.05,
                                                         step = 0.5,
                                                         min = 0)
                                            )
                                   )
                                   
                          ),
                          
                          tabPanel("Interactive Volcano Plot",
                                   br(),
                                   h4("Interactive Volcano Plot"),
                                   
                                   fluidRow(
                                     column(9,
                                            helpText("The plotly graph object generated below allows for interaction with 
                                                     the generated S-score volcano plot via hovering over data points, clicking & 
                                                     dragging, or use of the menu in the top right. PNG images can be downloaded 
                                                     of the graph, but will not include labels or information from the hovertabs. 
                                                     Click ", code("Generate Plotly"), " to view the volcano plotly. Color and
                                                     threshold options will be the same as those selected in the ", em("Volcano Plot"),
                                                     " tab.")
                                            ),
                                     column(3,
                                            actionButton("sscore_interactive_volcano_button",
                                                         label = "Generate Plot",
                                                         style = "color: #fff; background-color: #9BD79A; 
                                                                    border-color: #9BD79A; border-radius: 10px; border-width: 2px")
                                            
                                            )
                                   ),
                                   
                                   hr(),
                                   
                                   fluidRow(
                                     column(9, 
                                            shinycssloaders::withSpinner(plotly::plotlyOutput("sscore_volcano_interactive", height = 800), type = 8)
                                            ),
                                     column(3,
                                            )
                                   )
                                   ),
                          
                          tabPanel("Venn Diagram",
                                   br(),
                                   h4("Venn Diagram"),
                                   
                                   fluidRow(
                                     column(9, 
                                            helpText("The Venn diagram contains the number of statistically significant differentially expressed rows
                                            in each dataframe based on S-score adjusted p-value calculations. Overlap is determined by identifiers - 
                                            since metabolite-based and protein-based datasets use different identifiers there will note be overlap
                                            between those datasets.")
                                            ),
                                     column(3,
                                            actionButton("sscore_venn_button",
                                                         label = "Generate Plot",
                                                         style = "color: #fff; background-color: #9BD79A; 
                                                                 border-color: #9BD79A; border-radius: 10px; border-width: 2px")
                                            )
                                   ),
                                   
                                   br(),
                                   hr(),
                                   shinycssloaders::withSpinner(plotOutput("sscore_venn", height = 800), type = 8)
                                   ),
                          
                          tabPanel("Enrichment",
                                   br(),
                                   
                                   h4("S-Score Enrichment"),
                                   
                                   fluidRow(
                                     column(9, helpText("The enrichment function uses the output of the s-score integration to
                                                        create a ranked gene list. It will not function if the integration has not been run.")
                                            ),
                                     column(3, actionButton("sscore_enrichment_button",
                                                            label = "Run Enrichment Analysis",
                                                            style = "color: #fff; background-color: #9BD79A; 
                                                                    border-color: #9BD79A; border-radius: 10px; border-width: 2px")
                                            )
                                   ),
                                   
                                   br(),
                                   br(),
                                   
                                   fluidRow(
                                     column(3, 
                                            bslib::card(bslib::card_title("Database"),
                                                        bslib::card_body(fillable = FALSE,
                                                                         uiOutput("sscore_gmts")))
                                            
                                     ),
                                     column(3, 
                                            bslib::card(bslib::card_title("Contrast"),
                                                        bslib::card_body(fillable = FALSE,
                                                                         uiOutput("sscore_enrichment_coef_options")))
                                            
                                     ),
                                     column(3, 
                                            bslib::card(bslib::card_title("Threshold"),
                                                        bslib::card_body(fillable = FALSE,
                                                                         numericInput("sscore_enrichment_pval",
                                                                                      "Adjusted p-value cutoff",
                                                                                      value = 0.05,
                                                                                      min = 0,
                                                                                      max = 1,
                                                                                      step = 0.05)))
                                            
                                     )
                                   ),
                                   
                                   tabsetPanel(
                                     
                                     #### OVERVIEW ------------------------------------
                                     tabPanel("Overview",
                                              
                                              br(),
                                              
                                              h4("Enriched Pathways"),
                                              
                                              helpText("The ", code("clusterProfiler"), " package is used for pathway enrichment."),
                                              
                                              br(),
                                              br(),
                                              
                                              helpText("The table generated below when enrichment analysis is run contains the names of 
                                                       the enriched pathways based on the species and database selections. GSEA enrichment will have another column with normalized enrichment 
                                                       score values, while general enrichment will have a column with the gene count ratio."),
                                              
                                              br(),
                                              br(),
                                              
                                              helpText("The plots in the additional tabs (except the heatmap) are generated using the ", code("enrichplot"), "package,
                                                       available on Bioconductor and maintained by Guangchuang Yu."),
                                              hr(),
                                              
                                              shinycssloaders::withSpinner(DT::dataTableOutput("sscore_enriched_table"), type = 8),
                                              
                                              downloadButton("sscore_enriched_table_download",
                                                             label = "Download Enrichment Table",
                                                             style = "color: #61A6F9; background-color: #D8EAFF; border-color: 
                                                                     #D8EAFF; border-radius: 10px; border-width: 2px")
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
                                                       shinycssloaders::withSpinner(plotOutput("sscore_enriched_dotplot", height = 1000), type = 8)
                                                ),
                                                column(3,
                                                       h5("Options"),
                                                       
                                                       textInput("sscore_dotplot_title",
                                                                 label = "Customize plot title",
                                                                 value = "Enriched Dot Plot"),
                                                       
                                                       numericInput("sscore_dotplot_categories",
                                                                    label = "Select number of pathways to include",
                                                                    value = 15),
                                                       
                                                       hr(),
                                                       
                                                       actionButton("sscore_dotplot_button",
                                                                    label = "Generate Plot",
                                                                    style = "color: #fff; background-color: #9BD79A; 
                                                                            border-color: #9BD79A; border-radius: 10px; border-width: 2px"),
                                                       
                                                       br(),
                                                       br(),
                                                       
                                                       radioButtons("sscore_dotplot_extension", 
                                                                    label = "Save As:",
                                                                    choices = c("PDF" = ".pdf", 
                                                                                "PNG" = ".png",
                                                                                "SVG" = ".svg"),
                                                                    inline = TRUE),
                                                       
                                                       downloadButton("sscore_dotplot_download",
                                                                      label = "Download",
                                                                      style = "color: #61A6F9; background-color: #D8EAFF; border-color: 
                                                                              #D8EAFF; border-radius: 10px; border-width: 2px")
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
                                                       shinycssloaders::withSpinner(plotOutput("sscore_enriched_cnetplot", height = 1000), type = 8)
                                                ),
                                                column(3,
                                                       h5("Options"),
                                                       
                                                       textInput("sscore_cnet_title",
                                                                 label = "Customize plot title",
                                                                 value = "Enriched Network Plot"),
                                                       
                                                       numericInput("sscore_cnet_categories",
                                                                    label = "Select number of categories to include",
                                                                    value = 5),
                                                       
                                                       selectInput("sscore_cnet_layout",
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
                                                       
                                                       selectInput("sscore_cnet_labels",
                                                                   label = "Label",
                                                                   choices = c("Category Nodes Only" = "category",
                                                                               "Gene Nodes Only"= "gene",
                                                                               "Category & Gene Nodes" = "all",
                                                                               "None" = "none")),
                                                       
                                                       hr(),
                                                       
                                                       actionButton("sscore_cnetplot_button",
                                                                    label = "Generate Plot",
                                                                    style = "color: #fff; background-color: #9BD79A; 
                                                                            border-color: #9BD79A; border-radius: 10px; border-width: 2px"),
                                                       
                                                       br(),
                                                       br(),
                                                       
                                                       radioButtons("sscore_cnet_extension", 
                                                                    label = "Save As:",
                                                                    choices = c("PDF" = ".pdf", 
                                                                                "PNG" = ".png",
                                                                                "SVG" = ".svg"),
                                                                    inline = TRUE),
                                                       
                                                       downloadButton("sscore_cnetplot_download",
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
                                                       shinycssloaders::withSpinner(plotOutput("sscore_enriched_emapplot", height = 1000), type = 8)
                                                ),
                                                column(3,
                                                       h5("Options"),
                                                       
                                                       textInput("sscore_emap_title",
                                                                 label = "Customize plot title",
                                                                 value = "Enrichment Map"),
                                                       
                                                       numericInput("sscore_emap_categories",
                                                                    label = "Select number of categories to include",
                                                                    value = 30),
                                                       
                                                       hr(),
                                                       
                                                       actionButton("sscore_emap_button",
                                                                    label = "Generate Plot",
                                                                    style = "color: #fff; background-color: #9BD79A; 
                                                                            border-color: #9BD79A; border-radius: 10px; border-width: 2px"),
                                                       
                                                       br(),
                                                       br(),
                                                       
                                                       radioButtons("sscore_emap_extension", 
                                                                    label = "Save As:",
                                                                    choices = c("PDF" = ".pdf", 
                                                                                "PNG" = ".png",
                                                                                "SVG" = ".svg"),
                                                                    inline = TRUE),
                                                       
                                                       downloadButton("sscore_emap_download",
                                                                      label = "Download",
                                                                      style = "color: #61A6F9; background-color: #D8EAFF; border-color: 
                                                                              #D8EAFF; border-radius: 10px; border-width: 2px")
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
                                                       shinycssloaders::withSpinner(plotOutput("sscore_enriched_upsetplot", height = 1000), type = 8)
                                                ),
                                                column(3,
                                                       h5("Options"),
                                                       
                                                       textInput("sscore_upset_title",
                                                                 label = "Customize plot title",
                                                                 value = "Enrichment UpSet Plot"),
                                                       
                                                       numericInput("sscore_upset_categories",
                                                                    label = "Select number of categories to include",
                                                                    value = 10),
                                                       
                                                       hr(),
                                                       
                                                       actionButton("sscore_upset_button",
                                                                    label = "Generate Plot",
                                                                    style = "color: #fff; background-color: #9BD79A; 
                                                          border-color: #9BD79A; border-radius: 10px; border-width: 2px"),
                                                       
                                                       br(),
                                                       br(),
                                                       
                                                       radioButtons("sscore_upset_extension", 
                                                                    label = "Save As:",
                                                                    choices = c("PDF" = ".pdf", 
                                                                                "PNG" = ".png",
                                                                                "SVG" = ".svg"),
                                                                    inline = TRUE),
                                                       
                                                       downloadButton("sscore_upset_download",
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
                                                       shinycssloaders::withSpinner(plotOutput("sscore_enriched_treeplot", height = 1000), type = 8)
                                                ),
                                                column(3,
                                                       h5("Options"),
                                                       
                                                       textInput("sscore_tree_title",
                                                                 label = "Customize plot title",
                                                                 value = "Enriched Tree Plot"),
                                                       
                                                       numericInput("sscore_tree_categories",
                                                                    label = "Select number of categories to include",
                                                                    value = 30),
                                                       
                                                       numericInput("sscore_tree_clusters",
                                                                    label = "Select number of clusters to include",
                                                                    value = 5),
                                                       
                                                       hr(),
                                                       
                                                       actionButton("sscore_tree_button",
                                                                    label = "Generate Plot",
                                                                    style = "color: #fff; background-color: #9BD79A; 
                                                                            border-color: #9BD79A; border-radius: 10px; border-width: 2px"),
                                                       
                                                       br(),
                                                       br(),
                                                       
                                                       radioButtons("sscore_tree_extension", 
                                                                    label = "Save As:",
                                                                    choices = c("PDF" = ".pdf", 
                                                                                "PNG" = ".png",
                                                                                "SVG" = ".svg"),
                                                                    inline = TRUE),
                                                       
                                                       downloadButton("sscore_tree_download",
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
                                                       actionButton("sscore_enriched_heatmap_button",
                                                                    label = "Generate Plot",
                                                                    style = "color: #fff; background-color: #9BD79A; 
                                                                            border-color: #9BD79A; border-radius: 10px; border-width: 2px")
                                                )
                                              ),
                                              hr(),
                                              br(),
                                              shinycssloaders::withSpinner(plotly::plotlyOutput("sscore_enriched_heatplot", height = "auto"), type = 8)
                                     ),
                                     #### KEGG PATHWAY DIAGRAM -------------------------
                                     
                                     tabPanel("KEGG",
                                              br(),
                                              h4("KEGG Pathway Enrichment Diagram"),
                                              helpText("Certain databases have associated visualization options. If you selected ",
                                                       code("KEGG"), " enrichment pathway an additional plot can be rendered below."),
                                              
                                              hr(),
                                              br(),
                                              fluidRow(
                                                column(9,
                                                       verbatimTextOutput("sscore_notKEGG_message"),
                                                       shinycssloaders::withSpinner(plotOutput("sscore_enriched_KEGGplot", height = "auto"), type = 8)
                                                ),
                                                column(3,
                                                       uiOutput("sscore_KEGG_pathway_options")
                                                )
                                              ),
                                     )
                                   )
                                   )
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
    
    ## NETWORK ANALYSIS ########################################################
    
    tabPanel("5. Network Analysis",
             sidebarLayout(
               sidebarPanel(width = 3,
                 h4("PCSF Options"),
                 
                 br(),
                 
                 h6("Input",
                    bslib::tooltip(bsicons::bs_icon("info-circle"),
                                   "PCSF analysis uses the data generated from the Single Omic Analysis
                                   or Multi Omics Integration tabs. If none of these have been run there 
                                   will be no options underneath this header.",
                                   placement = "right")),
                 
                 uiOutput("network_input_options"),
                 uiOutput("network_SOLM_dataset"),
                 
                 hr(),
                 h6("Contrast"),
                 uiOutput("network_SOLM_contrast"),
                 uiOutput("network_MOLM_contrast"),
                 
                 hr(),
                 h6("Cutoffs",
                    bslib::tooltip(bsicons::bs_icon("info-circle"),
                                   "Only those rows that meet the cutoff thresholds will be included in PCSF 
                                   analysis. If no rows meet the cutoffs the graphs will not render - adjust the
                                   cutoffs and re-run the analysis.",
                                   placement = "right")),
                 
                 fluidRow(
                   column(6, uiOutput("network_SOLM_pval_cutoff")),
                   column(6, uiOutput("network_SOLM_logfc_cutoff"))
                 ),
                 
                 fluidRow(
                   column(6, uiOutput("network_MOLM_pval_cutoff")),
                   column(6, uiOutput("network_MOLM_sscore_cutoff"))
                 ),
                 
                 hr(),
                 h6("Mu",
                    bslib::tooltip(bsicons::bs_icon("info-circle"),
                                   "A higher mu correlates to a smaller the number of Steiners, which correlates to higher stringency.",
                                   placement = "right")),
                 
                 numericInput("pcsf_mu",
                              label = "Mu Value",
                              value = 0.005,
                              step = 0.05,
                              min = 0,
                              max = 1),
                 
                 hr(),
                 h6("Input Overview"),
                 shinycssloaders::withSpinner(DT::dataTableOutput("PCSF_input_table"), type = 8)
               ),
               mainPanel(width = 9,
                 tabsetPanel(
                   ## PCSF -----------------------------------------------------
                   tabPanel("PCSF",
                            br(),
                            h4("PCSF"),
                            
                            helpText("Prize-Collecting Steiner Forest, or PCSF, is a graph optimization 
                                     approach for network-based interpretation of high-throughput data. For PCSF
                                     analysis a database of gene interactions is referenced which provides a
                                     score, type, and source for each interaction."),
                            
                            br(),
                            br(),
                            
                            tabsetPanel(
                              
                              ### GENE INTERACTION NETWORK ---------------------
                              tabPanel("Gene Interaction Network",
                                       fluidRow(
                                         column(8,
                                                br(),
                                                shinycssloaders::withSpinner(visNetwork::visNetworkOutput("PCSF_Nodes_Network", height = "auto"), type = 8)
                                                ),
                                         column(4,
                                                br(),
                                                h5("Interpretation"),
                                                br(),
                                                
                                                helpText("This network includes nodes based on the input file and shows the
                                                         associated/related genes and their interactions."),
                                                
                                                br(),
                                                br(),
                                                
                                                helpText("The edges in this graph represent an interaction between two genes. 
                                                         If the edge is not named, the interaction is 'interacts-with'. 
                                                         All other interaction types will be named."),
                                                
                                                br(),
                                                br(),
                                                
                                                helpText("Coloring is based on regulation. Up-regulated nodes are red, down-regulated
                                                         are blue, and all others are gray."),
                                                
                                                br(),
                                                br(),
                                                
                                                helpText("The size of the node is determined by the absoluate value of the log
                                                         fold change associated with that gene."),
                                                
                                                br(),
                                                br(),
                                                
                                                actionButton("PCSF_Nodes_Network_button",
                                                             label = "Generate Network",
                                                             style = "color: #fff; background-color: #9BD79A; 
                                                                      border-color: #9BD79A; border-radius: 10px; border-width: 2px")
                                                )
                                       )
                                       ),
                              
                              ### GENE INFLUENCE NETWORK -----------------------
                              tabPanel("Gene Influence Network",
                                       fluidRow(
                                         column(8,
                                                shinycssloaders::withSpinner(visNetwork::visNetworkOutput("PCSF_Influential_Network", height = "auto"), type = 8)
                                                ),
                                         column(4,
                                                br(),
                                                h5("Interpretation"),
                                                br(),
                                                
                                                helpText("This network organizes genes by their overall 'influence' using the package ",
                                                         code("influential"), ". This package has a function ", code("IVI"), ". Per
                                                         the package's documentation, 'IVI is the first integrative method for the 
                                                         identification of network most influential nodes in a way that captures all 
                                                         network topological dimensions. The IVI formula integrates the most important 
                                                         local (i.e. degree centrality and ClusterRank), semi-local (i.e. neighborhood 
                                                         connectivity and local H-index) and global (i.e. betweenness centrality and 
                                                         collective influence) centrality measures in such a way that both synergizes 
                                                         their effects and removes their biases.' The influence score assigned by this 
                                                         algorithm is the basis for this network visualization."),
                                                br(),
                                                br(),
                                                
                                                actionButton("PCSF_Influential_Network_button",
                                                             label = "Generate Network",
                                                             style = "color: #fff; background-color: #9BD79A; 
                                                                      border-color: #9BD79A; border-radius: 10px; border-width: 2px")
                                                )
                                       )
                                       ),
                              
                              ## PCSF ENRICHMENT -------------------------------
                              tabPanel("Enrichment",
                                       br(),
                                       h5("PCSF Enrichment"),
                                       br(),
                                       
                                       fluidRow(
                                         column(8, 
                                                helpText("PCSF enrichment is performed on a clustered PCSF network, clustered using the ",
                                                         code("cluster_louvain"), " function from the ", em("igraph"), " package. Enrichment
                                                         is done with the ", em("clusterProfiler"), " package."),
                                                br(),
                                                br(),
                                                actionButton("PCSF_enrichment_button",
                                                             label = "Run Enrichment",
                                                             style = "color: #fff; background-color: #9BD79A;
                                                                                                  border-color: #9BD79A; border-radius: 10px; border-width: 2px")
                                                ),
                                         column(4,
                                                bslib::card(bslib::card_title("Database"),
                                                            bslib::card_body(fillable = FALSE,
                                                                             uiOutput("PCSF_gmts"))
                                                            )
                                                )
                                       ),
                                       
                                       tabsetPanel(
                                         ### ENRICHMENT TABLE ------------------
                                         tabPanel("Table Overview",
                                                  br(),
                                                  shinycssloaders::withSpinner(DT::dataTableOutput("PCSF_enriched_table"), type = 8),
                                                  br(),
                                                  downloadButton("pcsf_enriched_table_download", 
                                                                 label = "Download Table",
                                                                 style = "color: #61A6F9; background-color: #D8EAFF; border-color: 
                                                                   #D8EAFF; border-radius: 10px; border-width: 2px")
                                         ),
                                         
                                         ### ENRICHMENT SUBNETWORK ----------------
                                         tabPanel("Enrichment SubNetwork",
                                                  br(),
                                                  shinycssloaders::withSpinner(visNetwork::visNetworkOutput("PCSF_Enrichment_Network", height = "100%"), type = 8)
                                                  ),
                                         
                                         ### ENRICHMENT NETWORK ----------------
                                         tabPanel("Enrichment Cluster Network",
                                                  br(),
                                                  shinycssloaders::withSpinner(visNetwork::visNetworkOutput("PCSF_Enrichment_Clusters_Network", height = "auto"), type = 8)
                                         )
                                       )
                                       )
                            )
                   )
                            
                   # tabPanel("CytoScape Download")
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
    
    ## GENERATE REPORT #########################################################
    tabPanel("Generate Report",
             fluidRow(
               h3("RMarkdown Report"),
               
               br(),
               br(),
               
               column(8,
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
                      
                      p("Please allow a few minutes for download. Load time will increase when more contrasts are selected.")
                      
               ),
               column(4,
                      bslib::card(bslib::card_title("Download Rmarkdown HTML Report"),
                                  bslib::card_body("Once the model design has been completed and all desired output options have been adjusted, click the
                                                   Download Report button to generate the HTML report. A loading bar will appear showing progression of
                                                   the download.",
                                                   textInput("report_name",
                                                             label = "Add a project name to customize outputs",
                                                             width = "100%",
                                                             placeholder = "Proteomics_Project_A1")),
                                  bslib::card_footer(downloadButton("report",
                                                                    label = "Download Report",
                                                                    style = "color: #61A6F9; background-color: #D8EAFF; border-color: 
                                                                   #D8EAFF; border-radius: 10px; border-width: 2px")))
               ),
               
               hr(),
               
               fluidRow(
                 column(4,
                        bslib::card(bslib::card_title("Model Design [REQUIRED]"),
                                    bslib::card_header("This section must be completed to generate the HTML report."),
                                    bslib::card_body(fillable = FALSE,
                                                     h6("Linear Model Factors"),
                                                     uiOutput("report_linear_factors"),
                                                     hr(),
                                                     
                                                     h6("Contrasts"),
                                                     checkboxInput("report_contrast_fit",
                                                                   label = "Perform model fit on contrasts",
                                                                   value = TRUE),
                                                     uiOutput("report_reference_factor"),
                                                     uiOutput("report_choose_contrasts"),
                                                     hr(),
                                                     
                                                     h6("Covariate(s)"),
                                                     checkboxInput("report_add_covariate",
                                                                   label = "Add covariate(s) to the model"),
                                                     uiOutput("report_covariate_col"),
                                                     hr(),
                                                     
                                                     h6("Time Series"),
                                                     checkboxInput("report_include_timeseries",
                                                                   label = "Include time as interaction term"),
                                                     uiOutput("report_choose_time"),
                                                     uiOutput("report_time_type"),
                                                     uiOutput("report_time_points")))
                 ),
                 column(4,
                        bslib::card(bslib::card_title("Enrichment"),
                                    bslib::card_header(checkboxInput("report_include_enrichment",
                                                                     label = "Include enrichment",
                                                                     value = FALSE)),
                                    bslib::card_body(fillable = FALSE,
                                                     # uiOutput("report_database_header"),
                                                     # uiOutput("report_GSEA_gmts"),
                                                     # uiOutput("report_PKSEA_gmts"),
                                                     uiOutput("report_enrichment_calculation"),
                                                     uiOutput("report_enrichment_plots")))
                 ),
                 column(4,
                        bslib::card(bslib::card_title("PCSF",
                                                      bslib::tooltip(bsicons::bs_icon("info-circle"),
                                                                     "PCSF networks are generated for limma model outputs and, if selected, S-score outputs.",
                                                                     placement = "right")),
                                    bslib::card_header(checkboxInput("report_include_PCSF",
                                                                     label = "Include PCSF",
                                                                     value = FALSE)),
                                    bslib::card_body(fillable = FALSE,
                                                     uiOutput("report_pcsf_plots"),
                                                     uiOutput("report_pcsf_enrichment"),
                                                     uiOutput("report_PCSF_gmts"))),
                        
                        bslib::card(bslib::card_title("S-Score Multiomic Integration"),
                                    bslib::card_header(checkboxInput("report_include_sscore",
                                                                     label = "Include S-score",
                                                                     value = FALSE)),
                                    bslib::card_body(fillable = FALSE,
                                                     uiOutput("report_sscore_dataset"),
                                                     uiOutput("report_sscore_plots")))
                        
                        
                 )
               ),
               
               br()
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
