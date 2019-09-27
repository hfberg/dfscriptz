# A basic script to extract certain types of cells from a seurat object and process with SingleR. THis will greatly speed up the cluster identifiaction process.

####################################### Load libraries ####################################

library(Seurat)
library(dplyr)
library(ggplot2)
devtools::load_all("~/Tools/dftoolz")
library(SingleR)

################################## Subset Seurat object #######################################

#this code assumes you have a seurat object from which you want to extract some cells. for example if a lot of cells are too bland to reveil the cell type.
# Fill in the name of the project in "samp.name" below. 
# The "idents" choose which clusters in the seurat object you want to use. To see what else you can subset, type in "?subset" or "?WhichCells" in the console.

samp.name= "CP1"
seurat_sub<-subset( seurat, cells=WhichCells(object=seurat, idents=c(4,5,6,7)))

################################### Prepare library ######################################
# Download the immgen library to the local computer and save it in the object Immgen.

ImmGen <- ImmGenData()

#reduce the extent of the library by evaluating only the genes found in your data set.

common <- intersect(rownames(seurat_sub@assays[["RNA"]]@data), rownames(ImmGen))
ImmGen<-ImmGen[common,]
################################### Process with SingleR ######################################

#Create a singler object with main lables

singler.main <- SingleR(test = as.matrix(seurat_sub@assays[["RNA"]]@data), ref = ImmGen, labels = ImmGen$label.main)

#Create a singler object with detailed lables

singler.fine <- SingleR(test = as.matrix(seurat_sub@assays[["RNA"]]@data), ref = ImmGen, labels = ImmGen$label.fine)

######################################### Visualize ##############################################

# Now the singler objects are done, you can start visualizing what it found.
# The table will summarize number of cells identified for each cell type.
table(singler.main$labels)

# We can also plot a heatmap of the scores for each cell and cell type. This can be used for annotation.
# The clusters can be fetched from the seurat object as long as the barcodes exist as well.

plotScoreHeatmap(results=singler.main,show.pruned = TRUE, clusters=seurat_sub@active.ident, order.by.clusters = T)

############################# Annotate clusters with cell type ####################################
#Now you have identified what clusters have which cell type. You can name the cluster by that cell type.

#give as many new cell identities as there are clusters in the seurat object you want to plot.
new.cluster.ids <- c("mixed/ macrophages",  "mixed/ macrophages", "mixed/ macrophages", "mixed/ macrophages", "t-cells/ NK-cells/ ILC", "stromal cells/ fibroblasts/ macrophages/ monocytes", "mast cells", "stromal cells/ fibroblasts")

# If the slot called "active.ident" has any other value than the cluster values 0,1,2,3,4 and so on, do this:

seurat_sub$customized_identity<-seurat_sub@active.ident
Idents(seurat_sub)<-seurat_sub$seurat_clusters


#####continue here

#Assign each cluster cell type with a cluster level.
names(new.cluster.ids)<-levels(seurat_sub)

#rename identities
seurat_sub <- RenameIdents(seurat_sub, new.cluster.ids)

#Plot tSNE
DimPlot(seurat_sub, reduction = "tsne", label = TRUE, label.size = 4, pt.size = 0.5) + NoLegend()

###################################### Save singler and environment ###############################

# If SingleRBRowser starts working with our SingleR objects again, then save the RDS for uploading
saveRDS(singler.S4, paste0(samp.name,'.singleR.S4.rds'))

# Before you save your working environment to a file, make sure to remove all test-data 
# and unnecessary data by useing rm(name).
save.image(file=paste0(samp.name,".RData"))
