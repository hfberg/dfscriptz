# Assembled by Hanna F. Berg for Drug Farm
# A script for extracting cells of a specific cell type annotated by SingleR. The gene-cell-matrix can be used for further analysis in Seurat.

############################################### Load libraries ###################################################
library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
devtools::load_all("~/dftoolz")

############################################## Load SingleR.S4 ####################################################
# load the SingleR.S4 object for the sample to be subsetted either from the sample folder or a nwely generated SingleR object.
singler.S4<-readRDS("/media/drugfarm/Elements/proj.zhqugen/BM/result/u_CHFH_pre.vs.u_CFH_post/u_CHFH_pre_vs_u_CFH_post.samp.SingleR.S4.rds")


################################### Subset DGE (gene-cell exoression matrix) #######################################
# extract all annotated cell types and their cell barcodes from the SingleR.S4 object
cell.type<- as.data.frame(singler.S4@labels[["HPCA"]])

#extract cells with a particular cell type, for example plasma cells.
subset_cell<-"B_cell:Plasma_cell"
cell.type <- subset(cell.type, grepl(subset_cell,cell.type[,1]))

# generate a gene-cell expression matrix containing only the cells of a specific cell type.
cell.type <- as.vector(t(rownames(cell.type)))
DGE<-as.data.frame(singler.S4@expr)
DGE <- DGE %>% select(cell.type)

# Extract the sample name
samp.name = singler.S4@project.name
print(samp.name)

########################################## Filter and normalize ###################################################
# Generate Seurat object and visualize outliers
seurat <- CreateSeuratObject(counts = DGE)

seurat[["percent.mt"]] <- PercentageFeatureSet(object = seurat, pattern = "^MT-")
VlnPlot(object = seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

## Apply filters if needed, and normalize data
seurat <- subset(x = seurat, subset = nFeature_RNA > 1000 & nFeature_RNA < 4500 & nCount_RNA > 1000 & nCount_RNA < 3500 & percent.mt < 4 )
VlnPlot(object = seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

seurat <- NormalizeData(object = seurat)


############################################ Save DGE #############################################################
## Save the subsetted, filtered DGE
filen=paste0(samp.name,"_",subset_cell,"_filt")
DGE_save(DGE, samp.name = filen)
