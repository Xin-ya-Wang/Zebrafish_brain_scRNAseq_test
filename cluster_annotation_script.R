#setwd("F:/R_projects/cluster_annotation_test")

library(Seurat)
library(ggplot2)
library(tidyverse)
library(GeneNMF)



#dirs <- list.files("data/" ,full.names = F)

#rp.genes = c("rpl18a","rplp2l","rpl34","rpl13","rplp0","rpl36a","rpl12","rpl7a","rpl19","rpl3","rpl27","rpl23","rpl5b",
             "rplp2","rpl5a","rpl7","rpl37","rpl24","rpl9","rpl8","rpl31","rpl18","rpl28","rpl7l1","rpl6","rpl10a",
             "rpl13a","rpl39","rpl26","rpl4","rpl35a","rpl38","rplp1","rpl30","rpl11","rpl14","rpl10","rpl37.1","rpl35",
             "rpl17","rpl23a","rpl29","rpl22","rpl21","rpl22l1","rpl36","rpl32","rps16","rps13","rps4x","rps17","rps6ka3a",
             "rps6ka4","rps2","rps15a","rps11","rps19bp1","rps27a","rpsa","rps26l","rps10","rps28","rps8a","rps3a","rps6",
             "rps27.2","rps19","rps9","rps6kc1","rps7","rps29","rps8b","rps6ka1","rps6ka5","rps6kl1","rps6kal","rps6ka2","rps24",
             "rps3","rps27.1","rps18","rps6kb1b","rps5","rps21","rps26","rps12","rps14","rps6kb1a","rps25","rps15","rps23",
             "rps6ka3b","RPS17","RPL41")


#for(x in dirs){
  #cts <- read.table(gzfile(paste0("data/",x)),row.names=1)
  #gene_indices_to_remove <- which(rownames(cts) %in% rp.genes)
  #cts <- cts[-gene_indices_to_remove,]
  #name <- gsub("GSM333.*_lab1474_","",x) 
  #name <- gsub("_CountTable.tab.gz","",name) 
  #assign(name, CreateSeuratObject(counts=cts,min.cells=5, min.features=200,project=name))
  #}




#merged_seurat <- merge(AB42,y=c(IL4,PBS),add.cell.ids = c("AB42","IL4","PBS"))
#merged_seurat$cell <- rownames(merged_seurat@meta.data)
#merged_seurat$mt_percent <- PercentageFeatureSet(merged_seurat,pattern="^mt-")

save(merged_seurat,file="merged_seurat.Rdata")

VlnPlot(merged_seurat,c("nFeature_RNA","nCount_RNA","mt_percent"),group.by = "orig.ident")
seurat_filtered <- subset(merged_seurat,subset=nFeature_RNA>500 & nFeature_RNA<2500 & nCount_RNA>1000 & nCount_RNA<15000
                                 &mt_percent<6)
VlnPlot(seurat_filtered,c("nFeature_RNA","nCount_RNA","mt_percent"),group.by = "orig.ident")

seurat_filtered1=seurat_filtered


seurat_filtered <- NormalizeData(seurat_filtered)
seurat_filtered <- ScaleData(seurat_filtered)
hist(colSums(seurat_filtered,slot="scale.data"),
     breaks=100,
     main="Total expression after normalization")
seurat_filtered  <- FindVariableFeatures(seurat_filtered)


seurat_filtered <- RunPCA(seurat_filtered)
ElbowPlot(seurat_filtered)
seurat_cca <- IntegrateLayers(object=seurat_filtered,method=CCAIntegration,orig.reduction = "pca",new.reduction="integrated.cca",verbose=FALSE)



seurat_cca <- JoinLayers(seurat_cca)
seurat_cca<- FindNeighbors(seurat_cca, dims = 1:30, verbose = FALSE,reduction="integrated.cca")
seurat_cca <- FindClusters(seurat_cca, resolution=1.0 ,verbose = FALSE)
seurat_cca <- RunTSNE(seurat_cca,dims=1:30,reduction="integrated.cca")
DimPlot(seurat_cca,reduction="tsne",group.by = c("seurat_clusters","orig.ident"),label=TRUE)

DefaultAssay(seurat_cca)

markers1 <- FindAllMarkers(seurat_cca,only.pos=TRUE,logfc.threshold=1)
write.csv(markers1,file="markers.1.csv")


#top5 marker genes and heatmap
markers1 %>% group_by(cluster) %>% top_n(n=5, wt=avg_log2FC)%>% dplyr::arrange(desc(avg_log2FC))%>% dplyr::arrange(cluster)->top5
DoHeatmap(seurat_cca, features=top5$gene) + NoLegend()
VlnPlot(seurat_cca, features=top5$gene, pt.size=0)


#main cell type
main_markers<-c("lcp1","pfn1", #immune cells
                "her4.1","fabp7a","id1", "s100b",#progenitor cells(PCs)
               "sv2a","nrgna", #neuron
                "olig1","olig2") #oligodenrocytes

p1<-FeaturePlot(seurat_cca,features=main_markers,ncol=3)
ggsave("1.pdf",p1,width=15,height=10)

p2<-DotPlot(seurat_cca,features=main_markers)
ggsave("2.pdf",p2,width=14,height=6)

p3<-VlnPlot(seurat_cca,features=main_markers,stack=T,flip=T)+NoLegend()
ggsave("3.pdf",p3,width=14,height=6)

#immune cells: 13
#progenitor cells(PCs):2,3,4,5,6,7,10,12
#neuron:0,1,8,9,11,14
#oligodenrocytes: 15,16

seurat_cca$celltype_main <- recode(seurat_cca$seurat_clusters,
                                   "0"="neuron",
                                   "1"="neuron",
                                   "2"="progenitor cells",
                                   "3"="progenitor cells",
                                   "4"="progenitor cells",
                                   "5"="progenitor cells",
                                   "6"="progenitor cells",
                                   "7"="progenitor cells",
                                   "8"="neuron",
                                   "9"="neuron",
                                   "10"="progenitor cells",
                                   "11"="neuron",
                                   "12"="progenitor cells",
                                   "13"="immune cells",
                                   "14"="neuron",
                                   "15"="oligodenrocytes",
                                   "16"="oligodenrocytes")
pdf(file="4.pdf",width=7,height=6)
DimPlot(seurat_cca,reduction="tsne",group.by="celltype_main",label=T)
dev.off()


#progenitor cells subtypes
seurat_cca_PC=subset(seurat_cca,celltype_main %in% c("progenitor cells"))
seurat_PC=subset(seurat_filtered1,cell %in% seurat_cca_PC$cell)
seurat_PC<-seurat_PC%>% NormalizeData(verbose=F) %>% FindVariableFeatures(selection.method="vst", nfeatures=2000,verbose=F)%>%ScaleData(verbose=F)

seurat_PC <- RunPCA(seurat_PC)
ElbowPlot(seurat_PC)
seurat_cca_PC <- IntegrateLayers(object=seurat_PC,method=CCAIntegration,orig.reduction = "pca",new.reduction="integrated.PC.cca",verbose=FALSE)

seurat_cca_PC <- JoinLayers(seurat_cca_PC)
seurat_cca_PC<- FindNeighbors(seurat_cca_PC, dims = 1:30, verbose = FALSE,reduction="integrated.PC.cca")
seurat_cca_PC <- FindClusters(seurat_cca_PC, resolution=1.0 ,verbose = FALSE)
seurat_cca_PC <- RunTSNE(seurat_cca_PC,dims=1:30,reduction="integrated.PC.cca")
DimPlot(seurat_cca_PC,reduction="tsne",group.by = c("seurat_clusters","orig.ident"),label=TRUE)

marker_PC=c("gsx2","zic1","foxj1a","dmrta2","foxp4","nr2f1b","pou3f1","iqgap2","pcna","aurkb","cdk1","ascl1a")
p5<-FeaturePlot(seurat_cca_PC,features=marker_PC,ncol=3)
ggsave("5.pdf",p5,width=15,height=10)




#GeneNMF
seu.list <- SplitObject(seurat_cca_PC, split.by = "orig.ident")
geneNMF.programs <- multiNMF(seu.list, 
                             assay="RNA", slot="data", 
                             k=4:15, 
                             nfeatures = 2000)

geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
                                        nMP=9,
                                        weight.explained = 0.7,
                                        max.genes=100)
ph <- plotMetaPrograms(geneNMF.metaprograms)
ph

geneNMF.metaprograms$metaprograms.metrics

t(as.data.frame(lapply(geneNMF.metaprograms$metaprograms.genes, head)))

