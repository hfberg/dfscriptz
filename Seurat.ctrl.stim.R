# Assembled by Hanna F. Berg for Drug Farm based on and inspired by tutorials from Satijalabs Seurat. 

# A script for loading two samples, merging tem, removing batch effects, clustering and analyzing the differences.

######################################### Load libraries ########################################
library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
devtools::load_all("~/dftoolz")

####################################### Load samples ##############################################

ctrl_samp.name="RA3"
stim_samp.name = "RS3"

ctrl.data<-DGE_load(samp.name = ctrl_samp.name)
stim.data<-DGE_load(samp.name = stim_samp.name)

################# Set up control object
ctrl <- CreateSeuratObject(counts = ctrl.data, project = "ctrl")
ctrl$stim <- "ctrl"

# Preform quality control, subset if needed
ctrl[["percent.mt"]] <- PercentageFeatureSet(object = ctrl, pattern = "^mt-")
VlnPlot(object = ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ctrl <- subset(x = ctrl, subset = nFeature_RNA > 700 & nFeature_RNA < 1200 & nCount_RNA>700 & nCount_RNA< 2000 & percent.mt < 3.2)
VlnPlot(object = ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Check that ctrl object is scaled and normalized and filtered. If not, do this.

ctrl <- NormalizeData(object = ctrl, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable features 
ctrl <- FindVariableFeatures(object = ctrl, selection.method = "vst", nfeatures = 2000)


################## Set up stimulated object
stim <- CreateSeuratObject(counts = stim.data, project = "stim", min.cells = 5)
stim$stim <- "stim"

# Preform quality control, subset if needed
stim[["percent.mt"]] <- PercentageFeatureSet(object = stim, pattern = "^mt-")
VlnPlot(object = stim, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

stim <- subset(x = stim, subset = nFeature_RNA > 500 & nFeature_RNA <900 & nCount_RNA>500 & nCount_RNA< 1800 & percent.mt < 3.5 )
VlnPlot(object = stim, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Check if stim object is scaled and normalized adn filtered. If not, do this.
stim <- NormalizeData(object = stim, normalization.method = "LogNormalize", scale.factor = 10000)

#Find variable features
stim <- FindVariableFeatures(object = stim, selection.method = "vst", nfeatures = 2000)




######################################## Integrate samples #######################################

#Find gene markers to base the integration on
seurat.anchors <- FindIntegrationAnchors(object.list = c(ctrl,stim), dims = 1:30)
seurat.combined <- IntegrateData(anchorset = seurat.anchors, dims = 1:30)

###########
########### Start here for one sample instead of two 
###########
###########

# Add metadata information on which samples are ctrl and stim or successful and unsuccessful.

inf<-ctrl@meta.data[["stim"]]
inf<- append(inf, stim@meta.data[["stim"]])
seurat.combined$stim<-inf

# Or do it like this

seurat.combined$stim <-read.table("path/to/meta-data/file.txt", head =T, row.names = 1)[,1]


############################### Analysis of the integrated sample ##############################

# Re-set the assay to be investigated. We will cluster all cells based on the newly found anchor genes rather than all genes.
DefaultAssay(object = seurat.combined) <- "integrated"

# Find variable genes
seurat.combined <- FindVariableFeatures(object = seurat.combined, selection.method = "vst", nfeatures = 2000)

############# scale data is necessary for every assya in the seurat object, but do it only ONCE per object!
seurat.combined<-ScaleData(seurat.combined)

############## Cluster and visualize

# PCA. npcs can be changed
seurat.combined <- RunPCA(object = seurat.combined, npcs = 30, verbose = FALSE)
ElbowPlot(object = seurat.combined)

PCA_dimensions<-1:5
# t-SNE and Clustering, dims and resolution can be changed
seurat.combined <- RunTSNE(object = seurat.combined, reduction = "pca", dims = PCA_dimensions)
seurat.combined <- FindNeighbors(object = seurat.combined, reduction = "pca", dims = PCA_dimensions)
seurat.combined <- FindClusters(seurat.combined, resolution = 0.5)

# visualize clusters. Select which meta data you want to visualize in "group.by" or "split.by".
DimPlot(object = seurat.combined, reduction = "tsne", group.by = "samp.name")
DimPlot(object = seurat.combined, reduction = "tsne", group.by = "stim")
DimPlot(object = seurat.combined, reduction = "tsne", label = TRUE)

DimPlot(object = seurat.combined, reduction = "tsne", split.by = "stim")


################################ Without batch effect removal ###################################

#Test what the samples would look like without the batch effect removal. This creates a new seurat object and does not overwrite all work that has been done.
# also, this is not necessary for the final seurat object, but can be a good comparison to see what the samples would look like with batch effects.

assay.RNA<-GetAssay(object=seurat.combined, assay="RNA")
DGE.wbatch<-as.matrix(GetAssayData(object=assay.RNA, slot = "data"))
seurat.combined.wbatch<-CreateSeuratObject(counts=DGE.wbatch, project="with batch effects")
seurat.combined.wbatch$stim<-seurat.combined@meta.data[["stim"]]

seurat.combined.wbatch <- FindVariableFeatures(object = seurat.combined.wbatch, selection.method = "vst")
seurat.combined.wbatch<-ScaleData(seurat.combined.wbatch)

#make sure these parameters are the same as in the object with the reduced batch effects, otherwise they will not be comparable.
seurat.combined.wbatch <- RunPCA(object = seurat.combined.wbatch, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
seurat.combined.wbatch <- RunTSNE(object = seurat.combined.wbatch, reduction = "pca", dims = 1:30)
seurat.combined.wbatch <- FindNeighbors(object = seurat.combined.wbatch, reduction = "pca", dims = 1:30)
seurat.combined.wbatch <- FindClusters(seurat.combined.wbatch, resolution = 1)

# visualize clusters
DimPlot(object = seurat.combined.wbatch, reduction = "tsne", group.by = "samp.name")
DimPlot(object = seurat.combined.wbatch, reduction = "tsne", group.by = "stim")
DimPlot(object = seurat.combined.wbatch, reduction = "tsne", label = TRUE)


DimPlot(object = seurat.combined.wbatch, reduction = "tsne", split.by = "stim")

# You are now done with the comparison of the samples with batch effects. The rest of the script analyses the integrated, batch effect removed seurat object.

##################################### Differential expression ctrl vs stim #################################

##############If you already computed and saved the DE genes, load it now and skip to the Visualization part.
seurat.combined.markers<-read.table("/location/for/Gene markers samples.txt", head=T, sep="\t")

############### If this is the first time computing the DE genes, continue below:

# Re-set the assay to be investigated. We will look at differential expression in all genes rather than the anchor genes.
DefaultAssay(object = seurat.combined) <- "RNA"

# Choose what differential genes you want to compute. Between clusters, between successful vs unsuccessful, between samples and so on.
Idents(object = seurat.combined)<- as.factor(seurat.combined@meta.data[["seurat_clusters"]]) 

# DO ONLY ONCE: scale ALL RNA data to make sure to include all genes when searching for DE genes. 
seurat.combined <- ScaleData(seurat.combined,features = rownames(GetAssay(seurat.combined, slot = "RNA")))


# Find all markers calculates the differential expression. This takes a long time for big data sets, make sure you save the DE genes after this.
seurat.combined.markers <- FindAllMarkers(seurat.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#Save DE genes
write.table(seurat.combined.markers, file="Gene markers samples.txt", sep ="\t")


############################## Visualize the differentially expressed genes ####################################

#Select how many genes that will be plotted. 
top10 <- seurat.combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
# or:
top50 <- seurat.combined.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)

############## Heatmap

#the heatmap can not plot more than about 18000 cells in one go. Create subsets if the seurat.combined object contains more cells than 18000

seurat_small<-subset(x=seurat.combined, downsample=8000 )

pdf("Heatmap.pdf")
DoHeatmap(seurat_small, features=as.character(top50$gene), size = 2) + theme(axis.text.y = element_text(size = 1))
dev.off()

# if seurat.combined dontains less than 18000 cells, plot all cells.
DoHeatmap(seurat.combined, features = as.character(top50$gene),draw.lines = F, slot = "scale.data", size = 3) + theme(axis.text.y = element_text(size = 3))+ggtitle(paste("Heatmap", ctrl_samp.name, "vs", stim_samp.name))

############### Create new meta.data
# To compare for example all successful and all unsuccessful within one cluster, this data has to be created. 
# The function below assigns a name such as "cluster 1, successful" "and cluster 1 unsuccessful" for each cell. Then re-run findAllMarkers and heatmap to visualize.

samp.clus.ident=c()

for (i in 1:length(seurat.combined@active.ident)){
  nw<-paste0(seurat.combined@meta.data[["seurat_clusters"]][i]," ",seurat.combined@meta.data[["stim"]][i])
  samp.clus.ident<-append(samp.clus.ident, nw, after = i)
}

# set the new meta data to the identity to be evaluated.
Idents(object = seurat.combined)<- as.factor(samp.clus.ident)

############### Heatmap of individual clusters
#Plot a heatmap of only cluster x to compare successful and unsuccessful.

# Create an empty vector only once.
top_DE_ctrl_stim<-data.frame()

# subset the clusters to compare
seurat_subset<-subset(x = seurat.combined, idents = c("1 ctrl", "1 stim"))

# Find differentially expressed genes.
markers_subset <- FindAllMarkers(seurat_subset, only.pos = TRUE, min.pct = 0.25)

# Plot heatmap
pdf("Heatmap 1", height = 8.50, width = 11)
DoHeatmap(seurat_subset, features = markers_subset$gene, slot = "scale.data", size = 3) + NoLegend() + theme(axis.text.y = element_text(size = 3))
dev.off()

# extract the most differentiated genes between successful and unsuccessful in this cluster.
top1_sub <- markers_subset %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)

# save those genes in a data frame.
top_DE_ctrl_stim<- rbind(top_DE_ctrl_stim, as.data.frame(top1_sub))


############################################### Dot plot #############################################
# Dot plot is a quick and useful tool to compare some specific genes between samples, clusters or whatever is in the active identity.

# Plot top differential genes
DotPlot(object = seurat.combined, features = rev(x = unique(top10$gene)), cols = c("blue", "red"),  dot.scale = 8, split.by = "stim") + theme(axis.text.x = element_text(size = 7, angle = 90))

# Plot the markers with the biggest difference between successful and unsuccessful within one cluster. 
DotPlot(object = seurat.combined, features = rev(x = unique(top_DE_ctrl_stim$gene)), cols = c("blue", "red"),  dot.scale = 8, split.by = "stim") + theme(axis.text.x = element_text(size = 7, angle = 90))

#plot selected features and title
features = c("GENE1", "GENE2", "GENE3")
DotPlot(object = seurat.combined, features = features, cols = c("blue", "red"),  dot.scale = 8, split.by = "stim") + theme(axis.text.x = element_text(size = 7, angle = 90))+ggtitle(paste("DotPlot", ctrl_samp.name, "vs", stim_samp.name))

# Add a vertical line
DotPlot(seurat.combined,features = features,cex.use=4)+ theme(axis.text.x = element_text(size = 7, angle = 90))+ geom_vline(xintercept=8.5, color="red")

DotPlot(seurat.combined,features = features,cex.use=4)+ theme(axis.text.x = element_text(size = 7, angle = 90))+ geom_vline(xintercept=45.5, color="red")+ annotate("text", x = 20, y=4.5, label = "Miscellaneous")+geom_vline(xintercept=113.5, color="red") + annotate("text", x = 79, y=4.5, label = "Enzymes")+geom_vline(xintercept=149.5, color="red") + annotate("text", x = 130, y=4.5, label = "Transcription factors and their modulators
")+ geom_vline(xintercept=154.5, color="red")+ annotate("text", x = 152, y=4.3, label = "Early response genes")+ geom_vline(xintercept=176.5, color="red")+ annotate("text", x = 166, y=2.4, label = "Growth factors, ligands, and their modulators
")+ geom_vline(xintercept=192.5, color="red")+ annotate("text", x = 185, y=4.3, label = "Stress response genes")+ geom_vline(xintercept=198.5, color="red")+ annotate("text", x = 195, y=4.5, label = "Acute phase proteins")+ geom_vline(xintercept=204.5, color="red")+ annotate("text", x = 200, y=2.5, label = "Cell adhesion molecules
")+ geom_vline(xintercept=210.5, color="red")+ annotate("text", x = 210, y=3.5, label = "Proteins involved in antigen presentation")+ geom_vline(xintercept=231.5, color="red")+ annotate("text", x = 221, y=4.5, label = "Immunoreceptors")+ annotate("text", x = 251, y=4.5, label = "Cytokines/Chemokines and their modulators
")+ggtitle(paste("DotPlot", ctrl_samp.name, "vs", stim_samp.name))


#features<-c("Tnfsf13b", "Prdm1", "Ccl5", "Ccl17", "Ccl19", "Ccl22", "Ccl28", "Cxcl1", "Cxcl10", "Cxcl3", "Tgfb1", "Icos", "Ifng", "Il1a", "Il1b", "Il1rn", "Il2", "Il6", "Il10", "Il11", "Il12b", "Il13", "Il15", "Il23a", "Il27", "Ebi3", "Ifnb1", "Cxcl5", "Iigp1", "Lta", "Ltb", "Ccl2", "Cxcl9", "Ccl3", "Ccl4", "Tnf", "Tnfsf10", "Tnfsf15", "Cd80","Cxcr5", "Ccr5", "Ccr7", "Cxcr2","Tnfrsf9", "Cd3g", "Cd38", "Cd40", "Cd48", "Cd83","Cd86", "Slc3a2", "Tnfrsf4", "F11r", "Fcgrt", "Il2ra", "Ighg2c", "Igkc", "Bdkrb1", "B2m", "Nod2", "Pglyrp1", "Tlr2", "Tnfrsf1b", "Trem1", "Cfb", "C3", "Psmb9", "Tap1", "Tapbp", "Cd44", "Eng", "Fn1", "Icam1", "Ncam1", "Tnc", "Vcam1",   "Lbp", "Ptx3", "Saa1",  "Saa3", "F3", "Plau", "Cyp2e1",  "Cyp7b1", "Ptgs2", "Fth1", "Gclc", "Gclm", "Hsp90aa1", "Alox5", "Nos2", "Map4k1", "Senp2", "Sod1", "Sod2", "Mx1", "Nqo1", "Inhba", "Angpt1", "Pik3ap1", "Bdnf", "Blnk","Bmp4", "Calcb",  "Fstl3", "Csf3", "Hgf", "Mdk",  "Spr", "Nrg1", "Spp1", "Pdgfb", "Pigf", "Penk", "Kitl", "Thbs1", "Thbs2", "Vegfc", "Wnt10b", "Tnfaip2","Egr1","Ier3", "Dctn4", "Klf10","Tnfaip3", "Tnip3",  "Bcl3", "Bmi1",  "Fos",  "Myc", "Rel", "E2f3", "Ahctf1", "Ier2", "Gata3", "Nr3c1", "Hif1a", "Hoxa9","Irf1", "Irf2", "Irf4", "Irf7", "Nfkbia", "Nfkbie", "Junb", "Kdm6b", "Creb3", "Nfkbiz", "Nfkb2", "Nfkb1", "Nr4a2", "Trp53","Spi1", "Relb", "Snai1", "Sox9", "Stat5a", "Tfec", "Twist1",  "Yy1", "Abcb9", "Gcnt1", "Adh1", "Aicda", "Amacr", "Arfrp1", "Ass1",   "Bace1", "Btk", "Ctsb", "Ctsl", "Cdk6", "Uggt1", "Chi3l1", "Rdh5", "Dpyd", "Dnase1l2","Lipg", "Gad1", "St8sia1", "Mmp9", "Gstp1", "Gzmb", "Hpse", "Hmox1", "Has1", "Dio2", "Ido1",   "Mthfr", "Dusp1", "Mmp3",   "Nos1", "Pde7a", "Pim1", "Plk3", "Pik3ca", "Ppp5c", "Prkaca", "Prkcd", "Plcd1", "Ptgis", "Ptges", "Ptpn1","Gnb2l1", "Rev3l", "Slfn2",  "Atr", "St6gal1", "Nuak2", "Sat1", "Supv3l1",  "Tgm1", "Tgm2", "Pafah2", "Upp1", "Xdh", "Asph", "Afp", "Mif", "Apobec2", "Apod", "Apoe", "Bgn", "Brca2", "Myoz1", "Cav1", "Cdkn1a", "Col1a2",  "Ccnd1", "Ccnd2", "Ccnd3", "Slc11a2", "Edn1", "F8", "Gadd45b", "Gnai2", "Mt3", "Lgals3", "Vps53", "Hmgn1","Ltf", "Lamb2", "Lcn2", "S100a4", "Mbp", "Slc16a1", "Tnip1", "Psme1", "Psme2", "Serpine1", "Pts", "Prf1", "Ppargc1b", "Pgk1", "Pomc", "Hk3", "Pten",  "Rbbp4", "Ripk2", "Serpine2", "S100a6", "Sh3bgrl",  "Skp2", "Sdc4", "Slc6a6", "Kcnk5", "Tfpi2", "Ticam1", "Trpc1", "Ube2m", "Ucp2",  "Vdr", "Vim")



############################################## Save seurat #########################################
samp.name= "The name of my sample"
filen<-paste0(samp.name,".Seurat.rds")
saveRDS(seurat.combined, file=filen)


