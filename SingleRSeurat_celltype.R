#A script to subset and analyze seurat objects based on SingleR annotated cell types and cellsub types

############################################### Cell type lists ##############################################################
#create lists with detailed annotatioins and main annotations fro SingleR.

celltype.fine<-singler.fine@rownames
celltype_fine<-cbind(celltype.fine, singler.fine@listData[["labels"]])

#choose a cell type from the list, in this example, macrophages are chosen, but any cell type can be chosen.

macrophages <- subset(celltype.fine, grepl(pattern = "Macrophages", x = celltype.fine[,2]))

############################################### Subset seurat #######################################################
# subset the seurat object with the cell type. In this example, the subseted seurat is named "seurat_macrophages" and 
# the original seurat object is just named "seurat"

seurat_macrophages<-subset(x = seurat, cells=macrophages[,1])

####################################### Differential expression all ctrl stim ########################################
# Calculate the differential expression and plot a heatmap for ctrl stim all cell sub types

# Calculate DE
seurat_macrophages.markers <- FindAllMarkers(seurat_macrophages, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top_n_genes <- seurat_macrophages.markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_logFC)

# Plot heatmap
DoHeatmap(seurat_macrophages, features=as.character(top_n_genes$gene), size = 2) + theme(axis.text.y = element_text(size = 10))

################################################## New identities ####################################################
# If you want to save what you currently have stored in "active.ident" store these with a new name in meta data.
seurat_macrophages$custom.ident<-seurat_macrophages@active.ident

# Set the new active ident to the macrophages from the list you jsut extracted.
Idents(seurat_macrophages)<-as.factor(macrophages[,2])

#create new identities with a combineation of cell type and stim. The new identities will be stored in a vector named "new.idents"
CreateCombinedIdentity(seurat = seurat_macrophages, ident1 = "active.ident", ident2 = "stim")

# replace the active ident with the new ident, and add barcode names.
Idents(seurat_macrophages)<- new.ident

################################## Differential expression cell sub type ctrl stim ##################################
# The new identities are set, and we can calculate the difference in gene expression between ctrl and stim in every cell type

# If the function HeatpamIndividualClusters are going to work, the numbered clusters must be replaced with cell names.
seurat_macrophages@meta.data[["seurat_clusters"]]<-as.factor(macrophages[,2])

# Runt DE analysis
HeatmapIndividualClusters(seurat = seurat_macrophages, top_n_genes = 2, text_size = 5)
