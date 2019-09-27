# Assembled by Hanna F. Berg for Drug Farm based on and inspired by tutorials from Satijalabs Seurat. 

# A script for loading one single sample which doesn't need batch effect removal. Filter, normalize, cluster and analyse.

################################################# Load libraries ##############################################
library(Seurat)
library(dplyr)
library(ggplot2)
devtools::load_all("~/Tools/dftoolz")

########################################### Load DGE create Seurat #############################################

samp.name = "CP1"
DGE <-dftoolz::DGE_load("/home/drugfarm/LYH_test_data/test_R/CP1_dge.txt.gz")
seurat <- CreateSeuratObject(counts = DGE, project = CP1, min.cells = 3, min.features = 200)

########################################## Filter cells ########################################################
seurat[["percent.mt"]] <- PercentageFeatureSet(object = seurat, pattern = "^mt-")

VlnPlot(object = seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

########### Subset
#Subset the cells if needed. The extreme outliers should be removed.
seurat <- subset(x = seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 1350 & nCount_RNA>500 & nCount_RNA< 3500 & percent.mt < 2.7)

#plot again to confirm no mistakes were made
VlnPlot(object = seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# normalize if it's not previously done
seurat <- NormalizeData(object = seurat, normalization.method = "LogNormalize", scale.factor = 10000)
seurat <- FindVariableFeatures(object = seurat, selection.method = "vst", nfeatures = 2000)
seurat <- ScaleData(object = seurat, features = rownames(seurat@assays[["RNA"]]@data))
################### calculate principal components and pick dimensions so analyze ###############################

seurat <- RunPCA(object = seurat, features = VariableFeatures(object = seurat), verbose = F)
ElbowPlot(object = seurat)
PCA_dimensions<-1:8

######################################### Cluster and visualize #################################################
seurat <- FindNeighbors(object = seurat, dims = PCA_dimensions, k.param = 10)
seurat <- FindClusters(object = seurat, resolution = 0.3)

seurat <- RunTSNE(object = seurat, dims = PCA_dimensions)
DimPlot(object = seurat, reduction = "tsne")

##################################### Differential expression  ##################################################

##############If you already computed and saved the DE genes, load it now and skip to the Visualization part.
seurat.markers<-read.table("/home/drugfarm/Gene markers samples CP1 19.09.04.txt", head=T, sep="\t")

############### If this is the first time computing the DE genes, continue below:
seurat.markers <- FindAllMarkers(object = seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#Save DE genes
write.table(seurat.markers, file="/home/drugfarm/LYH_test_data/test_R/GeneMarkersCP1.txt", sep ="\t")

############################## Visualize the differentially expressed genes ####################################

#Select how many genes that will be plotted. 
top_n_genes <- seurat.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)


############# Heatmap

DoHeatmap(object = seurat, features = as.character(top_n_genes$gene), size=7) + NoLegend() + theme(axis.text.y = element_text(size = 5))+ggtitle(paste("Heatmap", samp.name))


############################################## Dot plot #######################################################
# Dot plot is a quick and useful tool to compare some specific genes between samples, clusters or whatever is in the active identity.

###### Select some genes y  ou want to plot and make them available to DotPlot.
features<-seurat.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
features<- as.character(features$gene)

DotPlot(object = seurat,features = features)+ theme(axis.text.x = element_text(size = 7, angle = 90))


############################################## save  #######################################################

#Save the filtered and normalized DGE
 DGE_save(DGE=seurat@assays[["RNA"]]@data, samp.name="CP1_19.09.04")

# Save the seurat object.
filen<-paste0(samp.name,".Seurat.rds")
saveRDS(seurat, file=filen)
 