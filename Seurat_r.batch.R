# Assembled by Hanna F. Berg for Drug Farm based on and inspired by tutorials from Satijalabs Seurat. 

# A script to create seuratobject and remove batch effects from two samples already merged into one DGE (gene-cell expression matrix).


########################################## Load libraries ################################################
devtools::load_all("~/dftoolz")
library(Seurat)
library(dplyr)

############################## Load sample and create Seurat #############################################
# Load gene-cell expression matrix
DGE<- DGE_load("s_SH_SHWH+u_CHFH_CJL_KUXP_B_cell:Plasma_cell_filt", name.barc = F, path=T)

#load the meta data file with a sampe anem paired with a cell barcode.
annot<-read.table("/home/drugfarm/proj.zhqugen/PBMC/process (copy)/data/s_SH_SHWH+u_CHFH_CJL_KUXP_B_cell:Plasma_cell_filt_samp.txt", header = T, row.names = 1, sep ="\t")

# Create seurat object
seurat<-CreateSeuratObject(counts = DGE, meta.data = annot)

########################################### Split samples #########################################################

seurat.list <- SplitObject(object = seurat, split.by = "samp.name")

#log normalization of each individual sample
for (i in 1:length(x = seurat.list)) {
  seurat.list[[i]] <- FindVariableFeatures(object = seurat.list[[i]], selection.method = "vst", 
                                        nfeatures = 2000, verbose = FALSE)
}


############################################ Remove batch effects ###############################################

# Find common markers on which the batch effects can be evaluated
seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list, anchor.features = 2000, dims = 1:30)

# The slot "Assay" now holds a batch corrected matrix for all cells.
seurat.combined <- IntegrateData(anchorset = seurat.anchors, dims = 1:30)

########################################### Change script #######################################################

#Now the seurat object is done. Open the script called Seurat.ctrl.stim.R and continue with the same seurat object 
# at this part of the script:

###########
########### Start here for one sample instead of two 
###########
###########


################################### Cluster without batch removal ###############################################

#seurat_wbatch<- FindVariableFeatures(object = seurat, selection.method = "mean.var.plot", 
# nfeatures = 5000, verbose = FALSE)
#seurat_wbatch <- ScaleData(object = seurat_wbatch, verbose = FALSE)
#seurat_wbatch <- RunPCA(object = seurat_wbatch, npcs = 30, verbose = FALSE)
#seurat_wbatch <- RunTSNE(object = seurat_wbatch, reduction = "pca", dims = 1:30)
#DimPlot(object = seurat_wbatch, reduction = "tsne", group.by = "X.samp.name.")

#################################################################################################################

