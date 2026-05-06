##CosMX Analysis


library(ggplot2)
library(ggpubr)
library(ggrepel)
library(gridExtra)
library(matrixStats)
library(patchwork)
library(pheatmap)
library(Seurat)
library(RColorBrewer)
library(reshape2)

if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")

if (!require("tiledb", quietly = TRUE))
  remotes::install_github("TileDB-Inc/TileDB-R", force = TRUE, 
                          ref = "0.17.0")

if (!require("tiledbsc", quietly = TRUE))
  remotes::install_github("tiledb-inc/tiledbsc", force = TRUE, 
                          ref = "8157b7d54398b1f957832f37fff0b173d355530e")

library(tiledb)
library(tiledbsc)


tiledbURI <- "/Users/aliamiryousefi/Desktop/a7dbe561-8db0-468e-8b3f-6df996bba91a_TileDB/"
tiledb_scdataset <- tiledbsc::SOMACollection$new(uri = tiledbURI, 
                                                 verbose = FALSE)
counts <- tiledb_scdataset$somas$RNA$X$members$counts$to_matrix(batch_mode = TRUE) 
dim(counts)
counts[1:4,1:4]
counts<- counts[, order(colnames(counts))]

#reading the first normalized prop sample 
norm<- tiledb_scdataset$somas$RNA_normalized_69a4bf1f.9099.422a.9996.a3d589249eda_1$X$members$data$to_matrix(batch_mode = TRUE)
dim(norm)
norm[1:4, 1:4]
norm<- norm[, order(colnames(norm))]

metadata <- tiledb_scdataset$somas$RNA$obs$to_dataframe()
dim(metadata)
colnames(metadata)
metadata<-metadata[order(metadata$cell_id), ]
metadata[1:4, 1:10]
#cellCoords <- tiledb_scdataset$somas$RNA$obs$to_dataframe(
#  attrs = colnames(metadata))
#head(cellCoords)

cellCoords <- tiledb_scdataset$somas$RNA$obs$to_dataframe(
  attrs = c("x_FOV_px", "y_FOV_px", "x_slide_mm", "y_slide_mm", 
            "slide_ID_numeric", "Run_Tissue_name", "fov"))
head(cellCoords)

length(unique(metadata$fov))

ggplot(cellCoords, aes(x=x_slide_mm, y=y_slide_mm))+
  geom_point(alpha = 0.05, size = 0.01)+
  facet_wrap(~Run_Tissue_name)+
  coord_equal()+
  labs(title = "Cell coordinates in XY space")

transcriptCoords <- tiledb::tiledb_array(
  tiledb_scdataset$somas$RNA$obsm$members$transcriptCoords$uri,
  return_as="data.frame")[]
head(transcriptCoords)

#ggplot(transcriptCoords, aes(x = as.factor(z_FOV_slice), fill = CellComp)) + geom_bar(stat = "count")

slide <- 1
fov <- 15

slideName <- unique(cellCoords$Run_Tissue_name[cellCoords$slide_ID_numeric == 
                                                 slide])

fovCoords <- cellCoords[cellCoords$slide_ID_numeric == slide & 
                          cellCoords$fov==fov,]
fovTranscriptCoords <- transcriptCoords[transcriptCoords$slideID == slide & 
                                          transcriptCoords$fov==fov,]

targetCounts <- table(fovTranscriptCoords$target)

targets <- names(targetCounts[which(targetCounts >= 2500 & 
                                      targetCounts <= 3000)])
fovTranscriptCoords <- fovTranscriptCoords[fovTranscriptCoords$target %in% 
                                             targets,]
#par(mfrow = c(1, 2))
ggplot(fovCoords, aes(x=x_FOV_px, y=y_FOV_px))+
  geom_point(alpha = 0.5, size = 0.4, color = "black")+
  geom_point(data = fovTranscriptCoords, 
             mapping = aes(x=x_FOV_px, 
                           y=y_FOV_px, 
                           color = target), 
             size = 0.1, alpha = 0.5)+
  theme_bw()+
  coord_equal()+
  guides(colour = guide_legend(override.aes = list(size=1,
                                                   alpha=1)))+
  labs(color = "RNA Target", title = paste0("RNA Transcripts in\n", 
                                            slideName, "\nFOV", fov))



pca <- tiledb_scdataset$somas$RNA$obsm$members$dimreduction_approximatepca_8f3e3e90.03f8.4df2.b7c3.885e74d30467_1$to_matrix()
pca[1:4, 1:4]

neighbs <- tiledb_scdataset$somas$RNA$obsp$members$graph_nn_dbd21e5c.b8e7.4cd7.b677.e923ca1cb8dd_1_nn$to_seurat_graph()


RNA_seurat <- tiledb_scdataset$to_seurat(somas = c("RNA"), batch_mode = TRUE)
RNA_seurat

RNA_seurat@meta.data <- metadata

Idents(RNA_seurat) <- RNA_seurat$spatialclust_f1f2b5bd.0503.46e6.9b94.56c2e6af8944_1_assignments

markers <- FindMarkers(RNA_seurat, ident.1 = unique(Idents(RNA_seurat))[1],
                       logfc.threshold = 0.01, test.use = "roc",
                       only.pos = TRUE)

VlnPlot(RNA_seurat, features = head(rownames(markers), 2),
        log = TRUE, pt.size = 0)

Idents(RNA_seurat) <- RNA_seurat$spatialclust_15ce7f34.a7a6.4702.884e.e1b66051b69d_1_assignments

markers <- FindMarkers(RNA_seurat, ident.1 = unique(Idents(RNA_seurat))[1],
                       logfc.threshold = 0.01, test.use = "roc",
                       only.pos = TRUE)

VlnPlot(RNA_seurat, features = head(rownames(markers), 2),
        log = TRUE, pt.size = 0)

Idents(RNA_seurat) <- RNA_seurat$spatialclust_edc5ef94.9d28.4675.85fa.501ad94b1dbe_1_assignments

markers <- FindMarkers(RNA_seurat, ident.1 = unique(Idents(RNA_seurat))[1],
                       logfc.threshold = 0.0001, test.use = "roc",
                       only.pos = TRUE)

VlnPlot(RNA_seurat, features = head(rownames(markers), 3),
        log = TRUE, pt.size = 0)

colorColumn = "spatialclust_f1f2b5bd.0503.46e6.9b94.56c2e6af8944_1_assignments"
#colorColumn = "Cell-type"
umap_TileDB <- as.data.frame(tiledb_scdataset$somas$RNA$obsm$members$dimreduction_approximateumap_5a181577.058c.4b05.a1ba.f3f7aaa4df67_1$to_matrix())
umap_TileDB$colorBy <- tiledb_scdataset$somas$RNA$obs$to_dataframe(attrs = colorColumn)[[1]]


#umap specifications

colrs <- brewer.pal.info[brewer.pal.info$colorblind == TRUE, ]
colorSchemes <- c("PuOr", "Dark2", "Set2", "BrBG")
colrs <- colrs[colorSchemes,]
col_vec <- unique(unlist(mapply(brewer.pal, colrs$maxcolors, colorSchemes)))

col_vec <- col_vec[-grep(col_vec, pattern = "#F|#E|#D")]

plotting <- function(data, title, Xcol, Ycol, Xname, Yname, color, 
                     size = 0.02, alpha = 0.05){
  gp <- ggplot(data, aes(x = data[[Xcol]], y = data[[Ycol]], 
                         color = data[[color]]))+
    geom_point(size = size, alpha = alpha)+
    coord_equal()+
    labs(title = title, color = colorColumn,
         x = Xname, y = Yname)+
    scale_color_manual(values = col_vec)+
    guides(colour = guide_legend(override.aes = list(size=1,
                                                     alpha=1)))
  
  return(gp)
}

umapPlot <- function(data, title, colorBy = "colorBy", ...){
  return(plotting(data, title, Xcol = "APPROXIMATEUMAP5A181577058C4B05A1BAF3F7AAA4DF671_1", 
                  Ycol = "APPROXIMATEUMAP5A181577058C4B05A1BAF3F7AAA4DF671_2", Xname = "UMAP1", 
                  Yname = "UMAP2", color = colorBy, ...))
}

pcaPlot <- function(data, title, colorBy = "colorBy", ...){
  return(plotting(data, title, Xcol = "APPROXIMATEPCA8F3E3E9003F84DF2B7C3885E74D304671_1", 
                  Ycol = "APPROXIMATEPCA8F3E3E9003F84DF2B7C3885E74D304671_2", Xname = "APPROXIMATEPCA8F3E3E9003F84DF2B7C3885E74D304671_1", 
                  Yname = "APPROXIMATEPCA8F3E3E9003F84DF2B7C3885E74D304671_2", color = colorBy, ...))
}

xySlidePlot <- function(data, title, colorBy = "colorBy", ...){
  return(plotting(data, title, Xcol = "x_slide_mm", 
                  Ycol = "y_slide_mm", Xname = "x_slide_mm", 
                  Yname = "y_slide_mm", color = colorBy, ...))
}


umapGP_TileDB <- umapPlot(umap_TileDB, "TileDB - R")

ggpubr::ggarrange(umapGP_TileDB,common.legend = TRUE, legend = "right", ncol = 1)



pca<- as.data.frame(pca)
pca$colorBy <- tiledb_scdataset$somas$RNA$obs$to_dataframe(attrs = colorColumn)[[1]]
pcaPlot(pca, "PCA")


slideMetadata <- metadata[metadata$Run_Tissue_name == slideName & 
                            metadata$fov %in% c(1:length(unique(metadata$fov))),]
colorColumn = "Spatial-clust1"
xySlidePlot(data = slideMetadata, title = "Cell Types in Space",
            colorBy = "spatialclust_f1f2b5bd.0503.46e6.9b94.56c2e6af8944_1_assignments", size = 0.3, alpha = 0.5)
colorColumn = "Spatial-clust2"
xySlidePlot(data = slideMetadata, title = "Cell Types in Space",
            colorBy = "spatialclust_15ce7f34.a7a6.4702.884e.e1b66051b69d_1_assignments", size = 0.3, alpha = 0.5)
colorColumn = "Spatial-clust3"
xySlidePlot(data = slideMetadata, title = "Cell Types in Space",
            colorBy = "spatialclust_edc5ef94.9d28.4675.85fa.501ad94b1dbe_1_assignments", size = 0.3, alpha = 0.5)




##writing a function that read the expressed transcripts for each cell from the transcripCoor 

colorColumn = "Phenotype"

No_genes<- numeric(dim(metadata)[1])
phenotype<- rep("other", dim(metadata)[1])
ind<-0
for (i in metadata$cell_id[1000:10000]){
  ind<- ind+1
  test<- subset(transcriptCoords, transcriptCoords$cell_id== i)
  if (sum(test$target %in% c("CD3D","CD3E","CD3G", "CD2", "CD4", "CD8A", "CD8B", "MS4A1", "ITGAX", "CD163", "SOX10", "CDH1", "KRT19", "PTPRC")) > 0){
    if(sum(test$target %in% c("CD3D","CD3E","CD3G", "CD2", "CD4", "CD8A", "CD8B")) > 0) {phenotype[ind]<- "T_cell"}
    else if (sum(test$target %in% c("MS4A1")) > 0) {phenotype[ind]<- "B_cell"}
    else if (sum(test$target %in% c("ITGAX", "CD163")) > 0) {phenotype[ind]<- "Myeloid_cell"}
    else if (sum(test$target %in% c("SOX10")) > 0) {phenotype[ind]<- "Tumor"}
    else if (sum(test$target %in% c("CDH1", "KRT19")) > 0) {phenotype[ind]<- "Epithelial_cell"}
    else {phenotype[ind]<- "Immune_cell"}
    }
}

metadata<- cbind(metadata, phenotype)

slideMetadata <- metadata[metadata$Run_Tissue_name == slideName & 
                            metadata$fov %in% c(1:6),]

xySlidePlot(data = slideMetadata, title = "Cell Types in Space",
            colorBy = "phenotype", size = 0.3, alpha = 0.5)

metadata <- metadata[, -which(colnames(metadata)=="phenotype")]




####parallelization


transcriptCoords<<- transcriptCoords
metadata<<- metadata

phen<- function(x){
  phenotype<- rep("other", length(x))
  ind<-x[1]-1
  for (i in metadata$cell_id[x]){
    ind<- ind+1
    test<- subset(transcriptCoords, transcriptCoords$cell_id== i)
    if (sum(test$target %in% c("CD3D","CD3E","CD3G", "CD2", "CD4", "CD8A", "CD8B", "MS4A1", "ITGAX", "CD163", "SOX10", "CDH1", "KRT19", "PTPRC")) > 0){
      if(sum(test$target %in% c("CD3D","CD3E","CD3G", "CD2", "CD4", "CD8A", "CD8B")) > 0) {phenotype[ind]<- "T_cell"}
      else if (sum(test$target %in% c("MS4A1")) > 0) {phenotype[ind]<- "B_cell"}
      else if (sum(test$target %in% c("ITGAX", "CD163")) > 0) {phenotype[ind]<- "Myeloid_cell"}
      else if (sum(test$target %in% c("SOX10")) > 0) {phenotype[ind]<- "Tumor"}
      else if (sum(test$target %in% c("CDH1", "KRT19")) > 0) {phenotype[ind]<- "Epithelial_cell"}
      else {phenotype[ind]<- "Immune_cell"}
    }
  }
  return(phenotype[x])
}

# Load the required packages
library(foreach)
library(doParallel)

# Number of elements and number of cores to use
num_elements <- dim(metadata)[1]
num_cores <- detectCores()

# Register a parallel backend using doParallel
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Generate a sequence of numbers
batch_size <- 2
elements<- 10
batch<- seq(1,  elements, batch_size)

# Perform parallel computation using foreach
result <- foreach(j = batch, .combine = c) %dopar% {
  j_batch <- phen(c(j:(j*batch_size)))
  return(j_batch)
}

# Stop the parallel backend
stopCluster(cl)
registerDoSEQ()  # Revert to the sequential backend




#T Cell Markers:
#  CD3
#CD2
#CD4
#CD8
#B Cell Markers:
#  CD20
#Myeloid Cell Markers:
#  CD11c
#CD163
#Stem Cell and Developmental Markers:
#  SOX10
#Epithelial Cell Markers:
#  ECAD (Epithelial Cadherin)
#PANCK (Pan-Cytokeratin)
#General Immune and Hematopoietic Markers:
#  CD45



#Cell Surface Marker,Gene Name
#CD45,PTPRC
#CD4,CD4
#CD3,CD3D,CD3E,CD3G
#CD2,CD2
#CD20,MS4A1
#CD11c,ITGAX
#SOX10,SOX10
#ECAD,CDH1
#PANCK,KRT19
#CD8,CD8A,CD8B
#CD163,CD163



##leiden clustering 

#xySlidePlot(data = slideMetadata, title = "Cell Clusters in Space",
#            colorBy = "RNA_pca_cluster_default",
#            size = 0.3, alpha = 0.5)+
# labs(color = "Cell Cluster")


hist(log(Matrix::rowSums(counts)), 100, main = "Total counts log distribution", freq = F, xlab = "Counts")
hist((Matrix::colSums(counts)), 100, main = "Gene frequency", freq = T, xlab = "Counts")

norm_d00<- tiledb_scdataset$somas$RNA_normalized_d00f6be4.8ef8.4e3c.96ad.80c53b285551_1$X$members$data$to_matrix(batch_mode = TRUE)
norm_69a<- tiledb_scdataset$somas$RNA_normalized_69a4bf1f.9099.422a.9996.a3d589249eda_1$X$members$data$to_matrix(batch_mode = TRUE)
norm_7b4<- tiledb_scdataset$somas$RNA_normalized_7b4c3ad3.ccc4.4ac5.a982.1178ef6e942a_1$X$members$data$to_matrix(batch_mode = TRUE)