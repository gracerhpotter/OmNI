
# Prize-Collecting Steiner Forest (PCSF) Graph Optimization Approach
# PLoS Comput Biol. 2017 Jul; 13(7): e1005694. DOI: 10.1371/journal.pcbi.1005694.

# PCSF NETWORK ANALYSIS ########################################################

## MAKE PCSF NETWORK -----------------------------------------------------------
#' @title PCSF Network
#' 
#' @description
#' Create PCSF network using PCSF_rand function from PCSF package for use in createing PCSF
#' network visualization.
#' 
#' @param pcsf_input_data dataframe input - must contain gene_symbol column and value column (logFC or S-score)
#' @param pcsf_nval numeric input -  number of runs with random noise added edge costs
#' 
#' @return PCSF network, PCSF igraph list object
#'

makePCSFNetwork <- function(pcsf_input_data = pcsf_input_data,
                            pcsf_nval = 10,
                            mu = 0.005){
  
  # pcsf_input_data <- PCSF_input
  pcsf_ppi <- readr::read_rds("PCSF_files/pcsf_ppi.rds")
  
  names(pcsf_input_data)[2] <- 'value'
  
  # input prizes
  pcsf_input_prizes = pcsf_input_data %>% 
    dplyr::select(gene_symbol, value) %>% 
    dplyr::mutate(gene_symbol = stringr::str_squish(gene_symbol))
  
  pcsf_input_prizes$value = abs(pcsf_input_prizes$value) # PCSF takes only positive values
  pcsf_input_prizes = tibble::deframe(pcsf_input_prizes) 
  pcsf_input_prizes_df = tibble::enframe(pcsf_input_prizes)
  
  set.seed(123)
  pcsf_net = PCSF::PCSF_rand(ppi = pcsf_ppi, 
                             terminal = pcsf_input_prizes, 
                             n = pcsf_nval,
                             mu = mu) # higher mu = smaller the number of Steiners
  
  # assign("pcsf_net", pcsf_net, envir = .GlobalEnv)
  return(pcsf_net)
}

## INTERACTION NETWORK ---------------------------------------------------------
#' @title PCSF Interaction Network
#'
#' @description
#' Generates an interactive PCSF network where nodes are genes and edges are interactions.
#' 
#' @param pcsf_net PCSF network output from makePCSFNetwork()
#' @param pcsf_input_data dataframe input
#'
#' @return visNetwork::visIgraph object

### PERFORM PCSF NETWORK ANALYSIS
PCSFVisNodes = function(pcsf_net = pcsf_net,
                        pcsf_input_data = pcsf_input) {
  
  my_interactome <- readr::read_rds("PCSF_files/pcsf_interactome.rds")
  names(pcsf_input_data)[2] <- 'value'
  
  # EDGES
  my_pcsf_net_edges = as.data.frame(igraph::get.edgelist(pcsf_net)) %>%
    dplyr::mutate(COMBINATION = paste0(pmin(V1, V2), pmax(V1, V2))) %>%
    dplyr::left_join(., my_interactome, by = "COMBINATION") %>%
    dplyr::mutate(weights = WEIGHT) %>%
    dplyr::mutate(weights = tidyr::replace_na(weights, min(weights, na.rm = TRUE))) %>%
    dplyr::mutate(label = INTERACTION_TYPE) %>%
    dplyr::mutate(label = replace(label, label == "interacts-with", "")) %>%
    dplyr::select(V1, V2, weights, label)
  
  # assign("PCSF_edges", my_pcsf_net_edges, envir = .GlobalEnv)
  
  # NODES
  my_pcsf_net_nodes = data.frame(gene_symbol = c(my_pcsf_net_edges$V1, my_pcsf_net_edges$V2))
  my_pcsf_net_nodes = my_pcsf_net_nodes %>% 
    dplyr::distinct(gene_symbol) %>%
    dplyr::left_join(., pcsf_input_data, by = "gene_symbol") %>%
    dplyr::distinct(gene_symbol, .keep_all = TRUE) %>%
    dplyr::mutate(is_steiner = case_when(value == "NA" ~ "Yes", 
                                         value != "NA" ~ "No")) %>% 
    dplyr::mutate(is_steiner = replace_na(is_steiner, "Yes")) %>% # if no logfc then its a Steiner
    dplyr::mutate(regulation = case_when(value >= 0 ~ "Up",
                                         value <= -0 ~ "Down")) %>%
    dplyr::mutate(regulation = replace_na(regulation, "None")) %>%
    dplyr::mutate(value = replace_na(value, 0)) %>%
    dplyr::mutate(abs_value = log2(abs(value))) 
  
  steiner_value = min(my_pcsf_net_nodes$abs_value[is.finite(my_pcsf_net_nodes$abs_value)], na.rm = TRUE) 
  my_pcsf_net_nodes = my_pcsf_net_nodes %>% 
    dplyr::mutate_all(~ifelse(is.na(.x) | is.nan(.x) | .x == -Inf, steiner_value, .x))
  
  # Define igraph nodes and edge aesthetics. Use igraph-compatible names
  my_pcsf_net_nodes = my_pcsf_net_nodes %>%
    dplyr::mutate(color = case_when(regulation == "Up" ~ "#F14E2B",
                                    regulation == "Down" ~ "#3F5DF4",
                                    regulation == "None" ~ "#E7E7E8")) %>%
    dplyr::mutate(label = gene_symbol) %>%
    dplyr::mutate(label.cex = scales::rescale(abs_value, c(0.3, 0.9)))
  
  # assign("PCSF_nodes", my_pcsf_net_nodes, envir = .GlobalEnv)
  
  # create igraph object
  my_pcsf_net_vis = igraph::graph_from_data_frame(d = my_pcsf_net_edges, 
                                                  vertices = my_pcsf_net_nodes, 
                                                  directed = FALSE)
  igraph::write_graph(graph = my_pcsf_net_vis,
                      file = paste0("~/Downloads/pcsf_network_", Sys.Date(), ".graphml"),
                      format = "graphml")
  # plot
  set.seed(123)
  visIgraph_obj = visNetwork::visIgraph(igraph = my_pcsf_net_vis,
                                        layout = "layout_nicely",
                                        physics = TRUE,
                                        smooth = TRUE, 
                                        idToLabel = FALSE) %>%
    visNetwork::visNodes(shape = "circle")  %>%
    visNetwork::visEdges(smooth = TRUE,
                         value = 2) %>%
    visNetwork::visOptions(height = "1000px",
                           width = "1000px",
                           highlightNearest = TRUE,
                           nodesIdSelection = TRUE ) %>%
    visNetwork::visInteraction(navigationButtons = TRUE, 
                               dragNodes = TRUE,
                               dragView = TRUE, 
                               zoomView = TRUE) 
  
  return(visIgraph_obj)
}
  
## INFLUENTIAL NETWORK ---------------------------------------------------------
#' @title PCSF Influential Node Network
#' 
#' @description
#' Generates an interactive PCSF network where nodes are genes and they are colored
#' and sized based on measure of "influence" determined by the `influential` package.
#' 
#' @param pcsf_net PCSF network object
#' 
#' @return visNetwork::visIgraph object
#'

PCSFVisInfluential = function(pcsf_net = pcsf_net) {
  ### PERFORM INFLUENTIAL ANALYSIS ON PCSF NETWORK
  my_graph = pcsf_net
  my_graph_vertices = igraph::V(pcsf_net) 
  
  # Compute influential nodes
  set.seed(123)
  my_graph_ivi = influential::ivi(graph = my_graph, 
                                  vertices = my_graph_vertices, 
                                  weights = NULL, 
                                  directed = FALSE, 
                                  mode = "all",
                                  loops = TRUE, 
                                  d = 3, 
                                  scale = "range" )
  
  my_graph_ivi_df = data.frame(my_graph_ivi) %>% 
    dplyr::mutate(genes = rownames(.)) %>%
    dplyr::rename(influence_score = my_graph_ivi) %>%
    dplyr::select(genes, influence_score) %>%
    dplyr::arrange(desc(influence_score))
  
  # assign("PCSF_influential", my_graph_ivi_df, envir = .GlobalEnv)
  
  my_graph_ivi_edges = as.data.frame(igraph::get.edgelist(my_graph))
  my_graph_ivi_nodes = my_graph_ivi_df
  
  # create igraph object
  my_igraph_ivi_vis = igraph::graph_from_data_frame(d = my_graph_ivi_edges, 
                                                    vertices = my_graph_ivi_nodes, 
                                                    directed = FALSE) 
  
  # pass node and edge information
  igraph::V(my_igraph_ivi_vis)$label <- igraph::as_ids(igraph::V(my_igraph_ivi_vis))
  igraph::V(my_igraph_ivi_vis)$color <- colourvalues::colour_values(igraph::V(my_igraph_ivi_vis)$influence_score, palette = "spectral", include_alpha = FALSE)
  
  igraph::V(my_igraph_ivi_vis)$label.cex <- igraph::V(my_igraph_ivi_vis)$influence_score
  igraph::V(my_igraph_ivi_vis)$label.cex <- scales::rescale(igraph::V(my_igraph_ivi_vis)$label.cex, c(0.5, 1.5))
  igraph::V(my_igraph_ivi_vis)$size <- igraph::V(my_igraph_ivi_vis)$influence_score
  igraph::V(my_igraph_ivi_vis)$size <- scales::rescale(igraph::V(my_igraph_ivi_vis)$size, c(8, 24))
  
  igraph::write_graph(graph = my_igraph_ivi_vis,
                      file = paste0("~/Downloads/pcsf_influential_", Sys.Date(), ".graphml"),
                      format = "graphml")
  
  # plot
  visIgraph_obj = visNetwork::visIgraph(igraph = my_igraph_ivi_vis,
                                        layout = "layout_nicely",
                                        physics = TRUE,
                                        smooth = TRUE, 
                                        idToLabel = FALSE) %>%
    visNetwork::visNodes() %>%
    visNetwork::visEdges() %>%
    visNetwork::visOptions(height = "1000px",
                           width = "1000px",
                           highlightNearest = TRUE,
                           nodesIdSelection = TRUE ) %>%
    visNetwork::visInteraction(navigationButtons = TRUE, 
                               dragNodes = TRUE,
                               dragView = TRUE, 
                               zoomView = TRUE) 
  return(visIgraph_obj)
}

# PCSF ENRICHMENT ##############################################################

## CLUSTER PROFILER ------------------------------------------------------------
#' @title PCSF Run ClusterProfiler
#'
#' @description
#' Run clusterProfiler and format output.
#' 
#' @param clusters igraph::cluster_louvain(subnet)
#' @param subnet PCSF network
#' @param gmt name of GMT database file.
#' 
#' @return enrichment_result_complete

runClusterProfiler <- function(clusters,
                               subnet,
                               gmt) {
  
  # gmt = "HUMAN_KEGG"
  enrichment_result = as.list(1:length(clusters))
  enrichment_result_complete = as.list(1:length(clusters))
  enrichment_table = as.list(1:length(clusters))
  
  # GENERATE GMT REFERENCES
  GMT_file <- list.files("./GMTs", pattern = gmt)
  my_geneset = readr::read_delim(paste0("./GMTs/", GMT_file))
  
  my_geneset = dplyr::mutate(my_geneset, term = pathway) %>% dplyr::mutate(feature = stringr::str_squish(feature_ids))
  term2features = my_geneset %>% dplyr::select(term, feature) #TERM2GENE
  term2pathway = my_geneset %>% dplyr::select(term, pathway) #TERM2NAME
  length(unique(my_geneset$feature))
  
  set.seed(123)
  for (a in 1:length(clusters)) {
    
    # INDICATE CLUSTER NUMBER
    cluster_number <- a
    print(paste0("Cluster Number: ", a))
    
    assign("clusters", clusters, envir = .GlobalEnv)
    
    enrichment_input <- clusters[[a]]
    # background_genes <- unique(my_geneset$feature) # use geneset as background
    background_genes <- unique(igraph::V(subnet)$name) # use network as background
    
    
    tryCatch({
      enriched_pathways <- clusterProfiler::enricher(gene = enrichment_input,
                                                     universe = background_genes,
                                                     minGSSize = 4,
                                                     pAdjustMethod = "fdr",
                                                     pvalueCutoff = 0.05,
                                                     TERM2GENE = term2features,
                                                     TERM2NAME = term2pathway)
      
      # EMPTY DATAFRAME TO STORE OUTPUT
      enriched_pathways_df_a = data.frame(cluster = "NULL", 
                                          term = "NULL", 
                                          adj_pval = 1, 
                                          intersection = "NULL")
      # DATAFRAME WITH ENRICHER RESULTS
      enriched_pathways_df_b = base::data.frame(enriched_pathways@result) %>%
        tibble::remove_rownames() %>%
        dplyr::select(term = Description, 
                      adj_pval = p.adjust, 
                      intersection = geneID) %>%
        dplyr::mutate(intersection = stringr::str_replace_all(intersection, "/", ", "))
      
      # BIND DATAFRAMES
      enriched_pathways_df = dplyr::bind_rows(enriched_pathways_df_a, enriched_pathways_df_b) %>%
        dplyr::mutate(cluster = paste0("Cluster_", cluster_number) )
      
      # DATA TABLE OF TOP RESULTS
      res_table_top15 = enriched_pathways_df %>% 
        dplyr::top_n(n = 15, wt = -adj_pval) %>%
        dplyr::filter(adj_pval <= 0.05) %>%
        dplyr::select(term, adj_pval, intersection) %>%
        dplyr::arrange(adj_pval)
      
      
      enrich = "<!DOCTYPE html> <html> <head> <style>\n      table {font-family: arial, sans-serif; font-size: 10px; border-collapse: collapse;width: 100%;} td,\n      th { border: 1px solid #dddddd; text-align: center; padding: 5px;}\n      tr:nth-child(even) {background-color: #dddddd;}\n      </style> </head> <body>\n      <table> <tr>  <th>Term</th> <th>Adjusted P-value</th> <th>Intersection</th> </tr>"
      
      for (i in 1:nrow(res_table_top15)) {
        enrich = paste0(enrich, " <tr>")
        for (j in 1:ncol(res_table_top15)) {
          enrich = paste0(enrich, "<td>", res_table_top15[i, 
                                                          j], "</td>")
        }
        enrich = paste0(enrich, "</tr> ")
      }
      
      enrich = paste0(enrich, "</table> </body> </html>")
      enrichment_result[[a]] = enrich
      enrichment_result_complete[[a]] = enriched_pathways_df
      enrichment_table[[a]] = enriched_pathways_df_b
      
      # Further processing of enriched_pathways if needed
    }, error = function(e) {
      warning(paste("Error in enricher for cluster ", a, ": ", conditionMessage(e)))
      # You can choose to handle the error in any specific way or proceed to the next iteration
    })
  }
  
  # assign("enrich_output", list(enrichment_result, enrichment_result_complete), envir = .GlobalEnv)
  return(list(enrichment_result, enrichment_result_complete, enrichment_table))
}

## ENRICHMENT ------------------------------------------------------------------
#' @title PCSF Enrichment
#'
#' @description
#' Call gprofiler and map pathways to PCSF modules.
#' 
#' @param subnet pcsf network
#' @param gmt GMT database name
#' 
#' @return output = list(subnet, enrichment_tab)

PCSFModuleEnrichments <- function(subnet,
                                  gmt) {
  # subnet = pcsf_net
  # Perform clustering on the PCSF network
  set.seed(123)
  clusters = igraph::cluster_louvain(subnet)
  # clusters = igraph::cluster_fast_greedy(subnet) 
  # clusters = igraph::cluster_edge_betweenness(subnet) 
  
  # saveRDS(clusters, file = paste0(out_directory_rds, "pcsf_p005_cluster_memberships_", out_file_name_suffix, ".rds"))
  assign("subnet", subnet, envir = .GlobalEnv)
  
  enrich = runClusterProfiler(clusters,
                              subnet,
                              gmt)
  
  assign("enrich", enrich, envir = .GlobalEnv)
  
  enrichment = enrich[[1]]
  enrichment_complete = enrich[[2]]
  
  novals = which(unlist(sapply(enrich[[2]], function(x) is.null(dim(x)))))
  
  if (length(novals) > 0){
    enrichment_complete = enrichment_complete[-novals]
  } 
  
  if (length(enrichment_complete) == 0){
    return(NULL)
  }
  
  enrichment_tab = do.call(rbind, lapply(c(1:length(enrichment_complete)), function(x) data.frame(Cluster = x, enrichment_complete[[x]])))
  
  igraph::V(subnet)$group = clusters$membership
  igraph::V(subnet)$title = paste0("Cluster ", clusters$membership, ": Enrichment analysis")
  
  for (i in 1:length(igraph::V(subnet))) {
    igraph::V(subnet)$title[i] = paste0(igraph::V(subnet)$title[i], enrichment[[igraph::V(subnet)$group[i]]])
  }
  
  class(subnet) = c("PCSFe", "igraph")
  output = list(subnet, enrichment_tab)
  names(output) = c("subnet", "enrichment")
  # assign("ModuleEnrichmentOutput", output, envir = .GlobalEnv)
  return(output)
  
}

## ENRICHED NETWORK ------------------------------------------------------------
#' @title Run PCSF Enrichment
#' 
#' @description
#' Call PCSF Enrichment modules.
#' 
#' @param pcsf_net PCSF network
#' @param gmt GMT database name
#' 
#' @return Enrichment results, list enrichment and subnet.
#' 

pcsfRunEnrichment <- function(pcsf_net,
                              gmt){
  ### PERFORM PATHWAY ENRICHMENT ON PCSF NETWORK
  ## This step performs clustering and then calls enrichment on each cluster
  
  subnet = pcsf_net
  set.seed(123)
  pcsf_enrich_pathway = PCSFModuleEnrichments(subnet,
                                              gmt)
  
  assign("pcsf_enrich_pathway", pcsf_enrich_pathway, envir = .GlobalEnv)
  
  return(pcsf_enrich_pathway)
}

#-------------------------------------------------------------------------------
#' @title PCSF Enrichment Results Table
#' 
#' @description
#' Format enrichment results into table.
#' 
#' @param pcsf_enrich_pathway enrichment results output from pcsfRunEnrichment ()
#' 
#' @return Enrichment results dataframe
#' 

pcsfEnrichedTable <- function(pcsf_enrich_pathway){
  # Write enrichment results as a text file
  pcsf_enrich_pathway_df = data.frame(pcsf_enrich_pathway$enrichment)
  pcsf_enrich_pathway_df[pcsf_enrich_pathway_df == "NULL"] <- NA
  pcsf_enrich_pathway_df <- na.omit(pcsf_enrich_pathway_df)
  pcsf_enrich_pathway_df$cluster <- sub("Cluster_", "", pcsf_enrich_pathway_df$cluster)
  colnames(pcsf_enrich_pathway_df) <- c("Cluster_1", "Cluster_2", "Pathway", "Adj_PVal", "Intersection")
  
  return(pcsf_enrich_pathway_df)
}

#-------------------------------------------------------------------------------
#' @title PCSF Enrichment Contracted Network
#' 
#' @description
#' Visualize contracted network of PCSF enrichment results.
#' 
#' @param pcsf_enrich_pathway enrichment results output from pcsfRunEnrichment()
#' 
#' @return visNetwork::visIgraph object
#' 

pcsfEnrichedContracted <- function(pcsf_enrich_pathway){
  
  # Contract the modules
  my_g = pcsf_enrich_pathway$subnet
  contracted_graph = igraph::contract(my_g, igraph::V(my_g)$group)
  contracted_graph = igraph::simplify(contracted_graph)
  
  number_of_clusters = seq(1:length(contracted_graph))
  new_module_names = paste("Module-", number_of_clusters, sep = "")
  igraph::V(contracted_graph)$name = new_module_names
  igraph::V(contracted_graph)$group = new_module_names
  
  unique_clusters = unique(igraph::V(contracted_graph)$name)
  cluster_gene_counts = sapply(unique_clusters, length)
  igraph::V(contracted_graph)$size = cluster_gene_counts
  igraph::V(contracted_graph)$color = colourvalues::colour_values(igraph::V(contracted_graph)$name, palette = "spectral", include_alpha = FALSE)

  igraph::write_graph(graph = contracted_graph,
                      file = paste0("~/Downloads/pcsf_enriched_contracted_", Sys.Date(), ".graphml"),
                      format = "graphml")
  
  visIgraph_obj = visNetwork::visIgraph(igraph = contracted_graph,
                                        layout = "layout_nicely", #layout_nicely
                                        physics = TRUE,
                                        smooth = TRUE, 
                                        idToLabel = FALSE) %>%
    visNetwork::visNodes(font = list(size = 15),
                         shape = "circle") %>%
    visNetwork::visEdges(smooth = TRUE,
                         value = 3) %>%
    visNetwork::visOptions(height = "1000px",
                           width = "1000px",
                           highlightNearest = TRUE,
                           nodesIdSelection = TRUE,
                           selectedBy = "group") %>%
    visNetwork::visInteraction(navigationButtons = TRUE, 
                               dragNodes = TRUE,
                               dragView = TRUE, 
                               zoomView = TRUE) 
  
  return(visIgraph_obj)
}

#-------------------------------------------------------------------------------
#' @title PCSF Enriched Subnet Network
#' 
#' @description
#' Visualize subnet of PCSF enrichment results.
#'
#' @param pcsf_enrich_pathway enrichment results output from pcsfRunEnrichment()
#' 
#' @return visNetwork::visIgraph object
#'

pcsfEnrichedSubnet <- function(pcsf_enrich_pathway){
  # Create visNetwork object for nice visualization
  visIgraph_obj = pcsf_enrich_pathway$subnet
  
  igraph::V(visIgraph_obj)$cluster = gsub(":.*", "", igraph::V(visIgraph_obj)$title)
  igraph::V(visIgraph_obj)$color <- colourvalues::colour_values(igraph::V(visIgraph_obj)$cluster, palette = "spectral", include_alpha = FALSE)
  igraph::write_graph(graph = visIgraph_obj,
                      file = paste0("~/Downloads/pcsf_enriched_subnet", Sys.Date(), ".graphml"),
                      format = "graphml")
  
  visIgraph_obj = visNetwork::visIgraph(igraph = visIgraph_obj,
                                        layout = "layout_nicely", #layout_nicely
                                        physics = TRUE,
                                        smooth = TRUE, 
                                        idToLabel = FALSE) %>%
    visNetwork::visNodes(font = list(size = 30)) %>%
    visNetwork::visEdges(smooth = TRUE,
                         value = 3) %>%
    visNetwork::visOptions(height = "1000px",
                           width = "1000px",
                           highlightNearest = TRUE,
                           nodesIdSelection = TRUE,
                           selectedBy = "group") %>%
    visNetwork::visInteraction(navigationButtons = TRUE, 
                               dragNodes = TRUE,
                               dragView = TRUE, 
                               zoomView = TRUE)
  
  return(visIgraph_obj)
}
  
