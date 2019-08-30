# Interacive multimodal 3D tSNE plotting of scRNA sequencing datasets
# The following is a length of code generated to create nice 
# 3D tSNE plots of seurat v3.0.0 objects utilizing the visualization 
# package plot_ly

# R v3.5.3 (x64 bit) and RStudio v1.2.1335 (x64 bit) were used for running this code :)

# Seurat is a multimodal single Cell RNA seq analysis algorithm created by
# The Satija Lab. For more information please see: https://satijalab.org/seurat/

# Contributors (by their Github handles):
# @Dragonmasterx87 (Dept. of Cell Biology, UM)
# @msaadsadiq (Dept. of Electrical and Computer Engineering, UM)
# @andrewwbutler (Center for Genomics and Systems Biology, NYU)  
# Edited 2019-08-30 by Hanna F. Berg to suit df samples.

library(plotly)
library(Seurat)
library(dplyr)

# Construct a dataframe using data from your pre-clustered Seurat v3.0.0 object
# Here 'seurat_clusters' is list of numeric cluster identities, you can find it here: seurat[["seurat_cluster"]], 
# or seurat$seurat_clusters, where 'seurat' is a Seurat object created with Seurat v3.0.0
#seurat <- ThisIsWhateverseuratIs

# Re-run tSNE s that you have accurate calculations for all tSNE(s)
seurat <- RunTSNE(seurat,
                            reduction.use = "pca",
                            dims.use = 1:12,
                            dim.embed = 3)

# Extract tSNE information from Seurat Object
tsne_1 <- seurat[["tsne"]]@cell.embeddings[,1]
tsne_2 <- seurat[["tsne"]]@cell.embeddings[,2]
tsne_3 <- seurat[["tsne"]]@cell.embeddings[,3]


# Prepare a dataframe for cell plotting
plot.data <- Seurat::FetchData(object = seurat, vars = c("tSNE_1", "tSNE_2", "tSNE_3", 'ident'))

# Make a column of row name identities (these will be your cell/barcode names)
plot.data$label <- paste(rownames(plot.data))

# Plot your data, in this example my Seurat object had 21 clusters (0-20)
plot_ly(data = plot.data, 
        x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3, 
        color = ~ident, 
        colors = "Paired",
        type = "scatter3d", 
        mode = "markers", 
        marker = list(size = 2, width=2), # controls size of points
        text=~label, #This is that extra column we made earlier for which we will use
        hoverinfo="text") #When you visualize your plotly object, hovering your mouse pointer over a point shows cell names


########################################### Plot gene expression ########################################

# Say you wanto make a gene-expression 3D plot, where you can plot gene expression against a color scale
# Here using the same seurat object as above, we extract gene expression information for beta-actin 'ACTB'
# create a dataframe
plotting.data <- FetchData(object = seurat, vars = c("tSNE_1", "tSNE_2", "tSNE_3", "ACTB"))

# Say you want change the scale, so that every cell having an expression >1 will be one color
# Basically, you are re-adjusting the scale here, so that any cell having a certain expression 
# will light up on your 3D plot

# First make another column in your dataframe, where all values above 1 are re-assigned a value of 1
# This information is stored in the 'changed' column of your dataframe
plotting.data$changed <- ifelse(test = plotting.data$ACTB <1, yes = plotting.data$ACTB, no = 1)

# Add the label column, so that now the column has 'cellname-its expression value'
plotting.data$label <- paste(rownames(plotting.data)," - ", plotting.data$ACTB, sep="")

# Plot your data, in this example my Seurat object had 21 clusters (0-20), and cells express a gene called ACTB
plot_ly(data = plotting.data, 
        x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3, 
        color = ~changed, # you can just run this against the column for the gene as well using ~ACTB, the algorith will automatically scale in that case based on maximal and minimal values
        opacity = .5,
        colors = c('darkgreen', 'red'), 
        type = "scatter3d", 
        mode = "markers",
        marker = list(size = 5, width=2), 
        text=~label,
        hoverinfo="text"
)

# On running this code the HTML output should appear in RStudio. You can save the output as a
# HTML file. Once you have saved, just open the HTML file in any web browser (double click on the html- file
# and if asked select to open with any web browser like google chrome/safari/mozilla/explorer etc).
# It should be have all of the integrated features you saw in the RStudio output file.

# Thank you for reading and using this code to further your scRNAseq analysis!
# If you liked it, dont forget to acknowledge, fork and star!
# Have a wonderful day!!

