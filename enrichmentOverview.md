## Enrichment Analysis

> Enrichment analysis is a computational method used in bioinformatics and functional genomics to identify whether a particular set of genes or proteins is overrepresented or underrepresented in a given experimental or biological context. The analysis aims to uncover biological themes, pathways, or functions associated with a set of genes, helping researchers interpret high-throughput genomic or proteomic data.

### Enrichment Databases

The options under the `Database` header will depend on the selection of species and the format of the data input. The databases available fall under three categories: gene symbol only & phosphosite/kinase, and gene/metabolite combined.

The gene set enrichment databases contain pathway information from a variety of common databases depending upon the species, such as KEGG, Reactome, & GO. Gene Set Enrichment Analysis (GSEA) analyzes which of the pathways in that database seem to be enriched based on the enriched genes output from the linear modeling step.

Phosphosite and kinase enrichment use data from ssGSEAs PTM-SEA (or Post Translational Modification Set Enrichment Analysis) developed at Broad Institute which is a data set of perturbation sites, kinases, diseases, and pathways associated with post translational modifications. The IDs for the associated sites are notated via the protein Uniprot ID along with the amino acid and position, i.e. `Q12345;S000` where `Q12345` is the protein ID, `S` is the amino acid, and `000` is the associated position on the protein.

The gene/metabolite database includes both gene symbols and metabolite CHEBI identifiers along with their associated pathways.

For all databases and enrichment type the `clusterProfiler` package described below is used to perform either ranked or unranked pathway enrichment.

### Using `clusterProfiler`

Omics Notebook enrichment analysis uses the `clusterProfiler` package. For detailed documentation on this package follow [this link](https://guangchuangyu.github.io/software/clusterProfiler/).

The following selections are available for this enrichment:

1.  `Database`: Which database correlating the the chosen species will be referenced for enrichment. The options vary slightly based upon which species was chosen. There is an option at the bottom which combines all available databases together and references that master list.

> Please be aware that the larger the database you select the longer enrichment will take. Go All (all Go databases combined) and the All Available Databases options will take the longest. For species that have less data available (yeast, zebrafish, c. elegans, & fruit fly) load times will be shorter.

2.  `Enrichment Algorithm`: There are two enrichment algorithms - ranked enrichment analysis (run via GSEA) and unranked Over-Representation Analysis (ORA). The main difference between ranked and unranked enrichment is that ranked enrichment uses a *ranked* list which is created using the output of linear modeling and ranks genes using a metric that relies on the calculated p-value and log fold change, or S-score. ORA, or unranked analysis, *does not* use a ranked list, but rather takes the list of genes represented in the data and compares that to a database of predefined gene sets to see if there are pathways that are overrepresented by genes in the data set.

> For enrichment of S-score output in the Multi Omic Integration tab, ranked enrichment analysis will automatically be used with
the s-score value acting as the rank. Conversely, network enrichment of PCSF clusters will always use unranked enrichment.

3.  `Contrast`: Contrast options are dependent on which groups are included in the generation of the model. The contrast that is chosen will determine the gene list that is compared to the database as the gene list is generated using the model output (either ranked or unranked).
4.  `Threshold`: The p-value threshold provides a significance cutoff that is applied to the enrichment function and output.

Once enrichment analysis has been run there is an option to download a table of the enrichment results. this table contains pathway names, p-values, and enrichment score or gene count ratio as well as a list of the genes from the gene list that are a match for that pathway. This output can be used as input in cytoscape if additional network exploration is desired.

Finally, the enrichment results can be explored in the additional plot tabs. For both `Single Omic Analysis` and `Multi Omics Integration` these include a dot plot, gene set network, enrichment map, upset plot, tree diagram, and interactive heatmap, while the PCSF `Network Analysis` enrichment will generate interactive networks based on the enriched clusters.

