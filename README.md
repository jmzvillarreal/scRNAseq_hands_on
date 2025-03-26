# Single cell hands on Session 

```r
## Load the required packages

library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)

## Load the data 

## Load the NMU_D dataset
NMU_O_D.data <- Read10X(data.dir = "./NMU_O_D/")
### Initialize the Seurat object with the raw (non-normalized data).
NMU_O_D <- CreateSeuratObject(counts = NMU_O_D.data, project = "NMU_O_D", min.cells = 3, min.features = 200)

NMU_O_D

### How does a sparse matrix looks like ?

# Lets examine a few genes in the first thirty cells
NMU_O_D.data[c("Upk3a", "Krt5", "Psca"), 1:30]

# Sparse matrices are a much more efficient way of storing the data

dense.size <- object.size(as.matrix(NMU_O_D.data))
dense.size
sparse.size <- object.size(NMU_O_D.data)
sparse.size

# Load the NMU_P dataset
NMU_O_P.data <- Read10X(data.dir = "./NMU_O_P/")
# Initialize the Seurat object with the raw (non-normalized data).
NMU_O_P <- CreateSeuratObject(counts = NMU_O_P.data, project = "NMU_O_P", min.cells = 3, min.features = 200)

NMU_O_P

```
# Cell QC analysis 

```r
# Determine % of mitochondrial genes 

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
NMU_O_D[["percent.mt"]] <- PercentageFeatureSet(NMU_O_D, pattern = "^mt-")
NMU_O_P[["percent.mt"]] <- PercentageFeatureSet(NMU_O_P, pattern = "^mt-")


# Load housekeeping genes
HK_genes <- read.table("HK_genes_mouse.txt")
HK_genes <- as.vector(HK_genes$V1)

# Identify the housekeeping genes that are in the dataset
HK_genes_NMU_O_D <- HK_genes[HK_genes %in% rownames(NMU_O_D)]
HK_genes_NMU_O_P <- HK_genes[HK_genes %in% rownames(NMU_O_P)]

# Count the number of housekeeping genes expressed per cell
NMU_O_D$HK_genes <- colSums(GetAssayData(NMU_O_D, assay = "RNA", slot = "counts")[HK_genes_NMU_O_D, ] > 0)
NMU_O_P$HK_genes <- colSums(GetAssayData(NMU_O_P, assay = "RNA", slot = "counts")[HK_genes_NMU_O_P, ] > 0)

# Visualization of QC criteria
VlnPlot(NMU_O_D, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "HK_genes"), ncol = 4)
VlnPlot(NMU_O_P, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "HK_genes"), ncol = 4)
plot1 <- FeatureScatter(NMU_O_D, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2 <- FeatureScatter(NMU_O_D, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot3 <- FeatureScatter(NMU_O_D, feature1 = "percent.mt", feature2 = "HK_genes")
plot1 + plot2 + plot3

plot1 <- FeatureScatter(NMU_O_P, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2 <- FeatureScatter(NMU_O_P, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot3 <- FeatureScatter(NMU_O_P, feature1 = "percent.mt", feature2 = "HK_genes")
plot1 + plot2 + plot3
```

## Remove low quality cells

### We filter cells that have unique feature counts over 2,500 or less than 200
### We filter cells that have >20% mitochondrial counts
### We filter cells that have <40 expressed HK genes

```r

NMU_O_D
NMU_O_D <- subset(NMU_O_D, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 20 & HK_genes > 40)
NMU_O_D

NMU_O_P
NMU_O_P <- subset(NMU_O_P, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 20 & HK_genes > 40)
NMU_O_P
```

# Normalizing the data 
## “LogNormalize”: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor.

```r
# Log normalize

NMU_O_D <- NormalizeData(NMU_O_D, normalization.method = "LogNormalize", scale.factor = 10000)
NMU_O_P <- NormalizeData(NMU_O_P, normalization.method = "LogNormalize", scale.factor = 10000)
```

# Identification of highly variable features 

```r
NMU_O_D <- FindVariableFeatures(NMU_O_D, selection.method = "vst", nfeatures = 2000)
NMU_O_P <- FindVariableFeatures(NMU_O_P, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes

top10_D <- head(VariableFeatures(NMU_O_D), 10)
top10_P <- head(VariableFeatures(NMU_O_P), 10)

# plot variable features with and without labels

plot1 <- VariableFeaturePlot(NMU_O_D)
plot2 <- LabelPoints(plot = plot1, points = top10_D, repel = TRUE)
plot2

plot1 <- VariableFeaturePlot(NMU_O_P)
plot2 <- LabelPoints(plot = plot1, points = top10_P, repel = TRUE)
plot2
```

# Scaling the data 

```r
NMU_O_D <- ScaleData(NMU_O_D)
NMU_O_P <- ScaleData(NMU_O_P)

# Perform linear dimensional reduction 

NMU_O_D <- RunPCA(NMU_O_D, features = VariableFeatures(object = NMU_O_D))
NMU_O_P <- RunPCA(NMU_O_P, features = VariableFeatures(object = NMU_O_P))

DimPlot(NMU_O_D, reduction = "pca")
DimPlot(NMU_O_P, reduction = "pca")

# Visualize genes in first 2 PCAs 

VizDimLoadings(NMU_O_D, dims = 1:2, reduction = "pca")
VizDimLoadings(NMU_O_P, dims = 1:2, reduction = "pca")

# Heatmaps visualizing the first 15 PCs

DimHeatmap(NMU_O_D, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(NMU_O_P, dims = 1:15, cells = 500, balanced = TRUE)
```

# Determine the ‘dimensionality’ of the dataset 
## Randomly permutes a subset of data, and calculates projected PCA scores for these 'random' genes. Then compares the PCA scores for the 'random' genes with the observed PCA scores to determine statistical significance. End result is a p-value for each gene's association with each principal component.
```r

# JackStraw analysis can take a long time for large datasets

NMU_O_D <- JackStraw(NMU_O_D, num.replicate = 100,dims = 30)
NMU_O_D <- ScoreJackStraw(NMU_O_D, dims = 1:30)
JackStrawPlot(NMU_O_D, dims = 1:30)
ElbowPlot(NMU_O_D, ndims = 30)

NMU_O_P <- JackStraw(NMU_O_P, num.replicate = 100,dims = 30)
NMU_O_P <- ScoreJackStraw(NMU_O_P, dims = 1:30)
JackStrawPlot(NMU_O_P, dims = 1:30)
ElbowPlot(NMU_O_P, ndims = 30)
```

# Run non-linear dimensional reduction (UMAP/tSNE)

```r
NMU_O_D <- RunUMAP(NMU_O_D, dims = 1:20)
DimPlot(NMU_O_D, reduction = "umap")

NMU_O_P <- RunUMAP(NMU_O_P, dims = 1:20)
DimPlot(NMU_O_P, reduction = "umap")
```

# Cell cycle scoring 

```r
library(nichenetr)
# Preloaded seurat cycle genes (Human)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Convert Human genes to mouse homologues
s.genes <- convert_human_to_mouse_symbols(s.genes,version = 1)
g2m.genes <- convert_human_to_mouse_symbols(g2m.genes,version = 1)

NMU_O_D <- CellCycleScoring(NMU_O_D, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(NMU_O_D@meta.data)
DimPlot(NMU_O_D, group.by = "Phase", reduction = "umap")

NMU_O_P <- CellCycleScoring(NMU_O_P, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(NMU_O_P@meta.data)
DimPlot(NMU_O_P, group.by = "Phase", reduction = "umap")
```

# Doublet estimation 

## DoubletFinder can be broken up into 4 steps:
### 1.Generate artificial doublets from existing scRNA-seq data
### 2. Pre-process merged real-artificial data
### 3.Perform PCA and use the PC distance matrix to find each cell's proportion of artificial k nearest neighbors (pANN)
### 4.Rank order and threshold pANN values according to the expected number of doublets
```r
library(DoubletFinder)
```

## NMU_O_D
```r
# Estimate expected doublet rate (e.g., 5% of total cells)
nExp <- round(ncol(NMU_O_D) * 0.05)  

# Estimate best pK (This defines the PC neighborhood size used to compute pANN, expressed as a proportion of the merged real-artificial data.)
sweep.res <- paramSweep(NMU_O_D, PCs = 1:20, sct = FALSE)

sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

# Choose the best pK (highest bcmvn value)
best_pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

# Run DoubletFinder with best pK
NMU_O_D <- doubletFinder(NMU_O_D, PCs = 1:20, pN = 0.25, pK = best_pK, nExp = nExp, sct = FALSE)

# Add DoubletFinder results to metadata
NMU_O_D$DF_classification <- NMU_O_D$DF.classifications_0.25_0.27_239

# Visualize doublets using UMAP
DimPlot(NMU_O_D, group.by = "DF_classification") 

# Check percentage of doublets
table(NMU_O_D$DF_classification)
```

## NMU_O_P
```r
# Estimate expected doublet rate (e.g., 5% of total cells)
nExp <- round(ncol(NMU_O_P) * 0.05)  # Adjust percentage based on your experiment

# Estimate best pK (This defines the PC neighborhood size used to compute pANN, expressed as a proportion of the merged real-artificial data.)
sweep.res <- paramSweep(NMU_O_P, PCs = 1:20, sct = FALSE)

sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

# Choose the best pK (highest bcmvn value)
best_pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

# Run DoubletFinder with best pK
NMU_O_P <- doubletFinder(NMU_O_P, PCs = 1:20, pN = 0.25, pK = best_pK, nExp = nExp, sct = FALSE)

# Add DoubletFinder results to metadata
NMU_O_P$DF_classification <- NMU_O_P$DF.classifications_0.25_0.01_385

# Visualize doublets using UMAP
DimPlot(NMU_O_P, group.by = "DF_classification") 

# Check percentage of doublets
table(NMU_O_P$DF_classification)
```

# Clustering and markers analysis 
```r
library(clustree)
```

## NMU_O_D
```r
# Clustering and resolution selection

NMU_O_D <- FindNeighbors(NMU_O_D, dims = 1:20)
NMU_O_D <- FindClusters(NMU_O_D, resolution = c(0, 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
clustree(NMU_O_D)
Idents(NMU_O_D) <- NMU_O_D$RNA_snn_res.0.1
DimPlot(NMU_O_D, reduction = "umap")

# Cluster markers

NMU_O_D_0.1.markers <- FindAllMarkers(NMU_O_D, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
NMU_O_D_0.1.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5

# Visualization of markers of top5
DotPlot(NMU_O_D, features = top5$gene) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  scale_color_gradientn(colors = c("blue", "grey", "red"))

# Plot expression of canonical markers

FeaturePlot(NMU_O_D, features = c("Krt5", "Upk3a"), cols = c("grey", "red"), reduction = "umap")

# Annotate cell populations

new.cluster.ids <- c("Intermediate", "Basal", "Luminal")
names(new.cluster.ids) <- levels(NMU_O_D)
NMU_O_D <- RenameIdents(NMU_O_D, new.cluster.ids)
DimPlot(NMU_O_D, reduction = "umap", label = TRUE)
```

## NMU_O_P
```r
# Clustering and resolution selection

NMU_O_P <- FindNeighbors(NMU_O_P, dims = 1:20)
NMU_O_P <- FindClusters(NMU_O_P, resolution = c(0, 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
clustree(NMU_O_P)
Idents(NMU_O_P) <- NMU_O_P$RNA_snn_res.0.2
DimPlot(NMU_O_P, reduction = "umap")

# Cluster markers

NMU_O_P_0.2.markers <- FindAllMarkers(NMU_O_P, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
NMU_O_P_0.2.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5

# Visualization of markers of top5
DotPlot(NMU_O_P, features = top5$gene) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  scale_color_gradientn(colors = c("blue", "grey", "red"))

# Plot expression of canonical markers

FeaturePlot(NMU_O_P, features = c("Mki67","Krt14", "Psca"), cols = c("grey", "red"), reduction = "umap")

# Annotate cell populations

new.cluster.ids <- c("Basal", "Intermediate high", "Intermediate low", "Basal G2M")
names(new.cluster.ids) <- levels(NMU_O_P)
NMU_O_P <- RenameIdents(NMU_O_P, new.cluster.ids)
DimPlot(NMU_O_P, reduction = "umap", label = TRUE)
```

# Integrative analysis
```r

library(harmony)

# Merge two objects

NMU_O_merged <- merge(x = NMU_O_D, y = c(NMU_O_P),add.cell.ids = c("D","P"))
NMU_O_merged
NMU_O_merged <- JoinLayers(NMU_O_merged)

# Repeat preprocessing

NMU_O_merged <- NormalizeData(NMU_O_merged, normalization.method = "LogNormalize", scale.factor = 10000)
NMU_O_merged <- FindVariableFeatures(NMU_O_merged, selection.method = "vst", nfeatures = 2000)
NMU_O_merged <- ScaleData(NMU_O_merged)
NMU_O_merged <- RunPCA(NMU_O_merged, features = VariableFeatures(object = NMU_O_merged))

# Visualize both datsets
DimPlot(NMU_O_merged, reduction = "pca", group.by = "orig.ident")

# Correct batch effect by integration with harmony

NMU_O_integrated <- RunHarmony(NMU_O_merged, "orig.ident")
NMU_O_integrated <- NMU_O_integrated %>% RunUMAP(reduction = "harmony",  dims = 1:20)

# Cluster analysis of integrated dataset

NMU_O_integrated <- NMU_O_integrated %>%
  FindNeighbors(reduction = "harmony") %>%
  FindClusters(resolution = c(0, 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)) 
clustree(NMU_O_integrated)
Idents(NMU_O_integrated) <- NMU_O_integrated$RNA_snn_res.0.1

# Visualize both datsets
DimPlot(object = NMU_O_integrated, reduction = "umap",group.by = "orig.ident")

# Visualize clusters
DimPlot(object = NMU_O_integrated, reduction = "umap")

# Visualize the distribution of populations
plot_data <- as.data.frame(table(Idents(NMU_O_integrated), NMU_O_integrated$orig.ident))
colnames(plot_data) <- c("Cluster", "Sample", "Count")

ggplot(plot_data, aes(x = Sample, y = Count, fill = Cluster)) +
  geom_bar(stat = "identity", position = "fill") + 
  scale_y_continuous(labels = scales::percent) + 
  labs(fill = "Cell population"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 18),  
    axis.text.x = element_text(angle = 45, hjust = 1),  
    axis.title = element_blank(),  
    axis.text.y = element_blank()  
  )


# Annotate cell populations on the integrated dataset

new.cluster.ids <- c("Intermediate","Basal", "Luminal", "Basal G2M")
names(new.cluster.ids) <- levels(NMU_O_integrated)
NMU_O_integrated <- RenameIdents(NMU_O_integrated, new.cluster.ids)

# Visualize cell populations
DimPlot(NMU_O_integrated, reduction = "umap", label = TRUE)

# Calculate cell popuplation markers
NMU_O_integrated.markers <- FindAllMarkers(NMU_O_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(NMU_O_integrated.markers, "NMU_O_integrated_markers.txt", sep = '\t')

save.image(file = "scRNAseq.RData")
```
# Use single cell annotated populaions to create pseudobulk

```r
NMU_O_integrated$annotated_populations <- Idents(NMU_O_integrated)
NMU_O_pseudobulk <- AggregateExpression(NMU_O_integrated,group.by = "annotated_populations")
NMU_O_pseudobulk
```

# Leverage on the scRNAseq findings to project them on bulk study
# 1. To do so we will extract cell population markers and us those as gene signatures for single sample GSEA
## https://cloud.genepattern.org/gp/pages/index.jsf

```r
# Filter cell population markers by adjusted p-value (e.g., p_adj < 0.05)
significant_markers <- NMU_O_integrated.markers[NMU_O_integrated.markers$p_val_adj < 0.05, ]

# Create a list of gene sets (significant markers per cluster)
gene_sets <- lapply(unique(significant_markers$cluster), function(cluster_id) {
  cluster_markers <- significant_markers[significant_markers$cluster == cluster_id, ]
  significant_genes <- cluster_markers$gene  # All significant genes for this cluster
  return(significant_genes)
})

names(gene_sets) <- unique(NMU_O_integrated.markers$cluster)

# Convert mouse genes to human homologues
Basal_G2M_GS <- convert_mouse_to_human_symbols(gene_sets[["Basal G2M"]], version = 1)
Basal_GS <- convert_mouse_to_human_symbols(gene_sets[["Basal"]], version = 1)
Intermediate_GS <- convert_mouse_to_human_symbols(gene_sets[["Intermediate"]], version = 1)
Luminal_GS <- convert_mouse_to_human_symbols(gene_sets[["Luminal"]], version = 1)

Basal_G2M_GS <- as.data.frame(Basal_G2M_GS)
Basal_GS <- as.data.frame(Basal_GS)
Intermediate_GS <- as.data.frame(Intermediate_GS)
Luminal_GS <- as.data.frame(Luminal_GS)

gene_sets <- list(Basal_G2M_GS$Basal_G2M_GS,Basal_GS$Basal_GS,Intermediate_GS$Intermediate_GS,Luminal_GS$Luminal_GS)
names(gene_sets) <- unique(NMU_O_integrated.markers$cluster)

# Function to remove NA values from the gene set and save as GMT
save_as_gmt <- function(gene_sets, file_name) {

  file_conn <- file(file_name, open = "w")
  
  for (set_name in names(gene_sets)) {
        cleaned_genes <- na.omit(gene_sets[[set_name]])
       if (length(cleaned_genes) > 0) {
      line <- paste(set_name, "description", paste(cleaned_genes, collapse = "\t"), sep = "\t")
      writeLines(line, file_conn)
    }
  }
  
  close(file_conn)
}

# Save the cleaned gene sets as a GMT file
save_as_gmt(gene_sets, "NMU_O_gene_sets.gmt")
```
# 2. Use a subset of the scRNAseq data to deconvolute cell populations in the bulk data
## https://cibersortx.stanford.edu/
## Extract 10 cells per population to use it for Cibersort deconvolution

```r
# Get identities
idents <- unique(Idents(NMU_O_integrated))

# Subset 10 random cells per identity class
set.seed(123)  # For reproducibility
subset_cells <- unlist(lapply(idents, function(ident) {
  sample(Cells(NMU_O_integrated)[Idents(NMU_O_integrated) == ident], size = min(10, sum(Idents(NMU_O_integrated) == ident)), replace = FALSE)
}))

# Subset  object
NMU_O_subset <- subset(NMU_O_integrated, cells = subset_cells)

# Extract expression matrix
NMU_O_expr_matrix <- as.data.frame(GetAssayData(NMU_O_subset, slot = "counts")) 

# Label cell population
colnames(NMU_O_expr_matrix) <- Idents(NMU_O_subset)

# Convert mouse genes and save the subset as a matrix
human_converted <- as.data.frame(convert_mouse_to_human_symbols(rownames(NMU_O_expr_matrix), version = 1))
NMU_O_expr_matrix$human_converted <- human_converted$`convert_mouse_to_human_symbols(rownames(NMU_O_expr_matrix), version = 1)`
NMU_O_expr_matrix <- na.omit(NMU_O_expr_matrix)
NMU_O_expr_matrix <- NMU_O_expr_matrix[!duplicated(NMU_O_expr_matrix$human_converted), ]
rownames(NMU_O_expr_matrix) <- NMU_O_expr_matrix$human_converted
NMU_O_expr_matrix$human_converted <- NULL
write.table(NMU_O_expr_matrix, file = "NMU_O_expr_matrix.tsv", sep = '\t',col.names = NA)
```
