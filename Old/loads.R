
if (!("pacman" %in% installed.packages()[,"Package"])){install.packages('pacman', update = TRUE, ask = FALSE)}

pacman::p_load(BiocManager, shiny, DT, ggplot2, markdown, colourpicker, uwot, ggcorrplot,
               shinycssloaders, plotly, openxlsx, ggrepel, cowplot, Glimma, kableExtra, Biobase,
               NormalyzerDE, ComplexHeatmap, dplyr, rmarkdown, clusterProfiler, enrichplot)