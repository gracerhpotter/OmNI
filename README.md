
# OmNI Shiny App

## Overview

An interactive shiny application version of the original Omics Notebook, with increased functionalities.

**Omics Notebook Interactive (OmNI)** is an R-based, open-source, and modular framework engineered for streamlined multi-omics data integration and analysis across diverse data types, incorporating interactive visualizations at each processing step. OmNI performs differential expression analysis utilizing customizable linear models, accommodating various covariates and complex experimental designs. For cross-omic layer integration, OmNI employs a modified S-score statistic, ensuring sensitive detection of differential features. The framework also integrates network and metabolomics data, offering detailed insights into regulatory mechanisms through comprehensive enrichment analysis using multiple pathway databases. Outputs include interactive HTML reports, CSV/TSV files, and Cytoscape-compatible objects. OmNI is readily deployable in both local and high-performance computing environments, enabling scalable data processing.

An **example dataset** is included in the **Data Upload** tab so users can explore the workflow without uploading their own files.

`![Example Dataset Toggle](images/ExampleDataset.png)`

---

## 🚀 Getting Started

### Launch the App

Run locally:

```r
# Launch the app from GitHub
shiny::runGitHub("OmNI", "gracerhpotter")
```

```r
# Download ZIP and run from folder
shiny::runApp("path/to/OmNI")
```

Or access the deployed version:  
👉 [https://emili-laboratory.shinyapps.io/omni/](https://emili-laboratory.shinyapps.io/omni/)

---

## App Structure

The navigation bar includes an **About** dropdown and six primary workflow tabs.

---

### 🔍 ABOUT ▼

Contains introductory and reference information:
- **Overview:** Summary of OmNI’s capabilities and intended workflow.  
- **File Inputs:** Expected input formats and data layout examples.  
- **Linear Modeling:** Description of limma-based analysis steps.  
- **Enrichment:** Overview of supported enrichment databases.  
- **Documentation:** Links to R package references and methods.

---

### 1️⃣ DATA

Upload and explore your datasets or use the built-in **Example Data**.  

**Functionalities:**
- **Upload Data:** Import `.csv`, `.tsv`, or `.txt` expression file(s) and an annotation file (`.xlsx`).  
- **Missing Values:** Explore missingness using **VIM** plots.  
- **Statistical Summary:** Review descriptive statistics.  
- **Expression Matrices:** Compare **raw** and **normalized** matrices.  
- **Normalization Plots:** Visualize quality control via violin, QQ, RLE, and density plots.  

`![Data Upload](images/DataUploadOptions.png)`  

There are a number of options for cutomization of data processing. Most options have preselections, and there are many tooltips providing information on the parameters. Once data is uploaded be sure to check that the **Group** selection (pulled from the provided annotation file) aligns with the desired analysis.

---

### 2️⃣ BASIC ANALYSIS

Perform exploratory analysis and assess dataset structure.

**Tabs:**
- **Correlation Matrix:** Pairwise sample correlation heatmap.  
- **MD Plot:** Mean-difference distribution by contrast.  
- **UMAP:** Low-dimensional visualization of sample relationships.  
- **Heatmap:** Hierarchical clustering of samples/features.  
- **PCA:** Principal component visualization of group separation.  

---

### 3️⃣ SINGLE OMIC ANALYSIS

Conduct differential analysis within a single omic layer.

**Functionality:**
- **Limma Modeling:** Apply linear modeling to detect significant features.  
- **Fit Summary:** Interactive topTable of model statistics.  
- **Volcano Plot:** Fold-change vs. significance visualization.  
- **Differential Heatmap:** Cluster significant features.  
- **Enrichment:** Perform functional enrichment using **clusterProfiler**.  
   - Supports GO, KEGG, Reactome, WikiPathways, MSigDB, and others.  
   - Results shown as dot plots, enrichment tables, gene set networks, and other visualizations from `clusterPofiler`.

`![Linear Model](images/LinearModelOptions.png)`

---

### 4️⃣ MULTI-OMIC INTEGRATION

Integrate multiple omic datasets to identify shared biological patterns.

**Features:**
- **S-Score Integration:** Weighted combination of multi-omic significance.  
- **Summary Table:** Lists integrated statistics and combined directionality.  
- **Volcano Plot:** Display integrated differential results.  
- **Venn Diagram:** Highlight shared and unique omic hits.  
- **Enrichment:** Combined enrichment for integrated feature sets.

---

### 5️⃣ NETWORK ANALYSIS

Build molecular networks and perform cluster-level enrichment.

**Modules:**
- **PCSF Networks:** Generate Steiner-forest-based functional subnetworks.  
- **PPI Network:** Visualize protein–protein interactions.  
- **Influence Network:** Highlight highly connected or influential nodes.  
- **Clustered Enrichment:** Identify and annotate enriched subnetworks.  

---

### 🧾 GENERATE REPORT

Compile and export all results.

📸 *Screenshot placeholder:*  
`![Generate Report Tab](images/ReportPage.png)`

---

## Example Workflow

1. Click **Use Example Data** on the **DATA** tab.  
2. Proceed through **Tabs 1–5** sequentially.  
3. Review each output interactively.  
4. Generate and download the full HTML report.

---

## Notes
- Tabs are designed for sequential use but can be revisited at any time.  
- Hover over plots for tooltips and data details.  
- Large datasets may take several minutes; progress indicators show job status.


### Dependencies

Dependencies for the application are installed and loaded via the `loads.R` file. The first time opening the application may take 10-15+ minutes if your system does not have many of the packages installed. 

### Development & Data Information

Additional information on the contents of each file, and the source of referenced databases for enrichment, PCSF, and S-Score is in the `info` folder in `devInfo.txt`. Additional reference scripts for the creation of some the files are also in that folder.

### Example Data & Annotation Files

Example data and annotation files are present in the `example` folder, but can also be found in a [public Google Drive folder](https://drive.google.com/drive/folders/1lyzmIhorrZy_CKuxabi1Bv1cLHIblJhk?usp=drive_link).
