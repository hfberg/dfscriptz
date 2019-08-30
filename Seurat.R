# Assembled by Hanna F. Berg for Drug Farm based on and inspired by tutorials from Satijalabs Seurat. 

# A script for loading one single sample which doesn't need batch effect removal. Filter, normalize, cluster and analyse.

################################################# Load libraries ##############################################
library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
devtools::load_all("~/dftoolz")

########################################### Load DGE create Seurat #############################################

samp.name = "u_CJL_pre_plasma"
DGE <-DGE_load(samp.name="/home/drugfarm/proj.zhqugen/BM/result/u_CJL_pre/subset plasma/u_CJL_B_cell:plasma_cell_filt_dge.txt", path=T,name.barc = F)
seurat <- CreateSeuratObject(counts = DGE, project = samp.name, min.cells = 3, min.features = 500)

########################################## Filter cells ########################################################
seurat[["percent.mt"]] <- PercentageFeatureSet(object = seurat, pattern = "^MT-")
VlnPlot(object = seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

########### Subset
#Subset the cells if needed. The extreme outliers should be removed.
seurat <- subset(x = seurat, subset = nFeature_RNA > 1000 & nFeature_RNA < 4000 & nCount_RNA>1000 & nCount_RNA< 3500 & percent.mt < 10)

#plot again to confirm no mistakes were made
VlnPlot(object = seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# normalize if it's not previously done
seurat <- NormalizeData(object = seurat, normalization.method = "LogNormalize", scale.factor = 10000)
seurat <- FindVariableFeatures(object = seurat, selection.method = "vst", nfeatures = 2000)

################### calculate principal components and pick dimensions so analyze ###############################
seurat <- ScaleData(object = seurat, features = rownames(seurat@assays[["RNA"]]@data))
seurat <- RunPCA(object = seurat, features = VariableFeatures(object = seurat), verbose = F)
ElbowPlot(object = seurat)
PCA_dimensions<-1:18

######################################### Cluster and visualize #################################################
seurat <- FindNeighbors(object = seurat, dims = PCA_dimensions)
seurat <- FindClusters(object = seurat, resolution = 1)

seurat <- RunTSNE(object = seurat, dims = PCA_dimensions)
DimPlot(object = seurat, reduction = "tsne")

##################################### Differential expression  ##################################################

##############If you already computed and saved the DE genes, load it now and skip to the Visualization part.
seurat.markers<-read.table("/location/for/Gene markers samples.txt", head=T, sep="\t")

############### If this is the first time computing the DE genes, continue below:
seurat.markers <- FindAllMarkers(object = seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#Save DE genes
write.table(seurat.markers, file="Gene markers samples.txt", sep ="\t")

############################## Visualize the differentially expressed genes ####################################

#Select how many genes that will be plotted. 
top10 <- seurat.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top70 <- seurat.markers %>% group_by(cluster) %>% top_n(n = 70, wt = avg_logFC)
top200 <- seurat.markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_logFC)

############# Heatmap
pdf(paste0("Heatmap ", samp.name, ".pdf"))
DoHeatmap(object = seurat, features = as.character(top10$gene), size=7) + NoLegend() + theme(axis.text.y = element_text(size = 5))+ggtitle(paste("Heatmap", samp.name))
dev.off()

############################################## Dot plot #######################################################
# Dot plot is a quick and useful tool to compare some specific genes between samples, clusters or whatever is in the active identity.

DotPlot(seurat,features = features,cex.use=4)+ theme(axis.text.x = element_text(size = 7, angle = 90))

############################################## save  #######################################################

#Save the filtered and normalized DGE
 DGE_save(DGE=seurat@assays[["RNA"]]@data, samp.name="Lgao_PBS_filt")

# Save the seurat object.
filen<-paste0(samp.name,".Seurat.rds")
saveRDS(seurat, file=filen)
 