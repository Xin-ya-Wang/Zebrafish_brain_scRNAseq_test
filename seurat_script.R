setwd("F:/R_projects/Seurat_raw_script")

library(Seurat)#seurat v5.10.0
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(harmony)

dirs <- list.dirs("data/", recursive = F, full.names = F)

for(x in dirs){
  name <- gsub("_filtered_feature_bc_matrix","",x) 
  
  cts <- Read10X(paste0('data/',x))
  
  assign(name, CreateSeuratObject(counts=cts,min.cells=3, min.features=100))}

HB17_tumor <- subset(HB17_tumor,downsample=1000)
HB17_background<- subset(HB17_background,downsample=1000)
HB53_background<- subset(HB53_background,downsample=1000)
HB53_tumor<- subset(HB53_tumor,downsample=1000)

merged_seurat <- merge(HB17_background,y=c(HB17_tumor,HB53_background,HB53_tumor),add.cell.ids = ls()[2:5])
#re-run from Rdata change to ls()[2:5]

merged_seurat$sample_barcode <- rownames(merged_seurat@meta.data)
merged_seurat@meta.data <- separate(merged_seurat@meta.data,col="sample_barcode",into=c('Patient','Type','Barcode'),sep='_')
merged_seurat$sample <- paste0(merged_seurat$Patient,"_",merged_seurat$Type)
merged_seurat$mt_percent <- PercentageFeatureSet(merged_seurat,pattern="^MT-")

ncount_q <- quantile(merged_seurat@meta.data$nCount_RNA,c(0.025,0.975))
nfeature_q <-quantile(merged_seurat@meta.data$nFeature_RNA,c(0.025,0.975))
VlnPlot(merged_seurat,c("nFeature_RNA","nCount_RNA","mt_percent"),group.by = "sample")

merged_seurat_filtered <- subset(merged_seurat,subset=nFeature_RNA>as.numeric(nfeature_q)[1]&nFeature_RNA<as.numeric(nfeature_q)[2]&nCount_RNA>as.numeric(ncount_q)[1]&nCount_RNA<as.numeric(ncount_q[2])
                                 &mt_percent<20)
VlnPlot(merged_seurat_filtered ,c("nFeature_RNA","nCount_RNA","mt_percent"),group.by = "sample")



####normalize+scale+findvariable
obj.nsf <- NormalizeData(merged_seurat_filtered)
obj.nsf <- ScaleData(obj.nsf)
hist(colSums(obj.nsf,slot="scale.data"),
     breaks=100,
     main="Total expression after normalization")
obj.nsf <- FindVariableFeatures(obj.nsf)

obj.nsf <- RunPCA(obj.nsf)
ElbowPlot(obj.nsf)

obj.nsf <- RunUMAP(obj.nsf,reduction = "pca",dims=1:20,reduction.name = "umap.unintegrated")
p1 <- DimPlot(obj.nsf,reduction = "umap.unintegrated",group.by="sample")+ggtitle("logNormalized_unintegrated")
####ScTransform
obj.sct <- SCTransform(object=merged_seurat_filtered,vars.to.regress = c("mt_percent","nFeature_RNA","nCount_RNA"),verbose = TRUE)
obj.sct <- RunPCA(obj.sct,verbose=FALSE)
ElbowPlot(obj.sct)
obj.sct <- RunUMAP(obj.sct,reduction = "pca",dims=1:20,reduction.name = "umap.unintegrated")
p2 <- DimPlot(obj.sct,reduction = "umap.unintegrated",group.by="sample")+ggtitle("ScTransform unintegrated")


####nsf+CCA####
#没有使用目标的方法
obj.nsf.cca <- IntegrateLayers(object=obj.nsf,method=CCAIntegration,orig.reduction = "pca",new.reduction="integrated.cca",verbose=FALSE)
obj.nsf.cca <- RunUMAP(obj.nsf.cca,dims=1:30,reduction="integrated.cca",reduction.name="umap.cca")
p3 <- DimPlot(obj.nsf.cca,reduction = "umap.cca",group.by="sample")+ggtitle("logNormalized CCA")


#SCT+CCA
obj.sct.cca <- IntegrateLayers(object=obj.sct,normalization.method='SCT',method=CCAIntegration,orig.reduction = "pca",new.reduction="integrated.cca",verbose=FALSE)
#必须要指定前面normalization.method,不然会报错！
obj.sct.cca <- RunUMAP(obj.sct.cca,dims=1:30,reduction="integrated.cca")
p4 <- DimPlot(obj.sct.cca,reduction = "umap",group.by="sample")+ggtitle("SCtransform CCA")



#SCT+harmony
#obj.sct.harmony <- IntegrateLayers(object=obj.sct,assay='SCT',method=HarmonyIntegration,orig.reduction = "pca",new.reduction="sct.harmony",verbose=FALSE)

obj.sct.harmony <- obj.sct %>%
  RunHarmony(group.by.vars='sample')
obj.sct.harmony <- RunUMAP(obj.sct.harmony,dims=1:30,reduction="harmony")
p5 <- DimPlot(obj.sct.harmony,reduction = "umap",group.by="sample")+ggtitle("SCtransform harmony")

#Nsf+harmony
obj.nsf.harmony <- obj.nsf %>% RunHarmony(group.by.vars='sample')
obj.nsf.harmony <- RunUMAP(obj.nsf.harmony,dims=1:30,reduction="harmony")
p6 <- DimPlot(obj.nsf.harmony,reduction = "umap",group.by="sample")+ggtitle("logNormalized harmony")



####cluster###
#use obj,nsf.cca, for example
obj.nsf.cca1 <- JoinLayers(obj.nsf.cca)
obj.nsf.cca1 <- FindNeighbors(obj.nsf.cca1, dims = 1:30, verbose = FALSE,reduction = "integrated.cca")
obj.nsf.cca1 <- FindClusters(obj.nsf.cca1, resolution=0.3,verbose = FALSE)

p7 <- DimPlot(obj.nsf.cca1,reduction="umap.cca",group.by = c("seurat_clusters","sample"),label=TRUE)





###DEG
#a.FindAllMarkers for Annotation
obj_markers<-FindAllMarkers(
  object=obj.nsf.cca1,
  only.pos=TRUE, 
  test.use="wilcox", 
  slot="data", 
  min.pct=0.01, 
  logfc.threshold=0.25)
  

#b.FindConservedMarkers
for (i in unique(obj.nsf.cca1$seurat_clusters)){
  name <- paste0(print(i),"_","Markers")
  
  marker <- FindConservedMarkers(obj.nsf.cca1, ident.1=print(i),grouping.var="sample",slot="data")

  assign(name,marker)}

#c.FindMarkers#compare same cluster between samples
degs_cluster1 <- subset(obj.nsf.cca1,subset=seurat_clusters=="1")
unique(degs_cluster1@meta.data$sample)
degs_cluster1_HB17_tumor<-FindMarkers(degs_cluster1,
                          logfc.threshold=0.25,
                          ident.1="HB17_tumor", ident.2=c("HB17_background","HB53_background","HB53_tumor"),
                          group.by="sample")

#d.pseudo bulk for comparing same cluster between samples/different conditions
#compare DEG expression results by FindMarkers and pseudobulk
#for example: DEGs in cluster 2 between background and tumor

#(1)Find DE features between "background_cluster0" and "tumor_cluster0" by FindMarkers

#obj.nsf.cca1$cluster_type <- paste0(obj.nsf.cca1$seurat_clusters,"_",obj.nsf.cca1$Type)
#Idents(obj.nsf.cca1) <- "sample"
#cluster0_de_markers <- FindMarkers(obj.nsf.cca1,ident.1="0_tumor",ident.2 = "0_background",verbose=FALSE)

#(2)Find DE features between "background_cluster0" and "tumor_cluster0" after pseudobulking
#pseudobulk the counts based on type-cluster
#obj.nsf.cca2=obj.nsf.cca1
#pseudo_obj <- AggregateExpression(obj.nsf.cca2,assays="RNA",slot="counts",return.seurat=FALSE,group.by = c("seurat_clusters","sample"))
#cts <- pseudo_obj$RNA
#cts <- as.matrix(cts)

#cts.t <- t(cts)

#cts.t <- as.data.frame(cts.t )

#SplitRows<- gsub("_.*","",rownames(cts.t))

#cts.split <- split.data.frame(cts.t,
                 #f=factor(SplitRows))

#cts.split.modified <- lapply(cts.split, function(x){
  #rownames(x) <-  gsub(".*_","",rownames(x))
  #t(x)
#})

#counts_cluster0 <- cts.split.modified$'g0'

#colData <- data.frame(samples=colnames(counts_cluster0))
#colData<- colData %>%
  #mutate(condition=ifelse(grepl('tumor',samples),'tumor','control'))%>%
  #column_to_rownames(var='samples')


#library(DESeq2)
#dds <- DESeqDataSetFromMatrix(countData = counts_cluster0,
                              #colData = colData,
                              #design = ~condition)

#keep <- rowSums(count(dds)) >=10
#dds <- dds[keep,]

#dds <- DESeq(dds)

#dds <- resultsNames(dds)
#res <- results(dds,names="cluster0_tumor_vs_background")



#pseudo_obj <- AggregateExpression(obj.nsf.cca2,assays="RNA",slot="counts",return.seurat=FALSE,group.by = c("seurat_clusters","sample"))
#tail(Cells(pseudo_obj))

#pseudo_obj$celltype_type <- paste0(pseudo_obj$Type,"_",pseudo_obj$seurat_clusters)
#Idents(pseudo_obj) <- "celltype_type"
#bulk_DEG_cluster0 <- FindMarkers(object = pseudo_obj,
                       # ident.1="tumor_0",ident.2 ="background_0", 
                        #test.use="DESeq2")
#not work：ValidateCellGroups(object = object, cells.1 = cells.1, cells.2 = cells.2, : Cell group 1 has fewer than 3 cells



CreatePseudoBulkData <- function(raw.data,normalized.data,sample.id,fun="sum")
{
  # fun can be "sum" or "mean"
  library(muscat)
  library(SingleCellExperiment)
  
  if (fun=="sum")
  {
    sce <- SingleCellExperiment(assays=list(counts=raw.data),colData=DataFrame(sample_id=sample.id))
  } else if (fun=="mean") {
    sce <- SingleCellExperiment(assays=list(counts=normalized.data),colData=DataFrame(sample_id=sample.id))
  }
  pb <- aggregateData(sce, assay = "counts", fun=fun, by=c("sample_id"))
  return(pb)
}


RunPseudobulkMethod <- function(raw.data,normalized.data,individual,group,test="ROTS",sum.or.mean="sum")
{
  pb <- CreatePseudoBulkData(raw.data = raw.data,normalized.data = normalized.data, 
                             sample.id = individual,fun = sum.or.mean)
  
  if (test=="ROTS")
  {
    
    group <- apply(as.matrix(table(individual,group)),1,function(x) x[1]==0)
    group_names <- names(group)
    group <- as.numeric(group)
    names(group) <- group_names
    group <- group[colnames(pb)]
    
    library(SingleCellExperiment)
    library(edgeR)
    library(ROTS)
    
    if (sum.or.mean=="sum")
    {
      y <- DGEList(assay(pb), remove.zeros = T)
      y <- calcNormFactors(y)
      logcpm <- edgeR::cpm(y, normalized.lib.sizes=T, prior.count=1, log=T)
      resrots <- ROTS(data = logcpm, groups = group, seed = 1234)
    } 
    else if (sum.or.mean == "mean") 
    {
      y <- assay(pb)
      print(dim(y))
      y <- y[apply(y,1,sum)!=0,]
      print(dim(y))
      resrots <- ROTS(data = y, groups = group, seed = 1234)
    }
    rotsdf <- as.data.frame(resrots$logfc)
    rotsdf$gene <- rownames(rotsdf)
    rotsdf <- cbind(rotsdf, resrots$pvalue)
    rotsdf <- cbind(rotsdf, resrots$FDR)
    rownames(rotsdf) <- 1:nrow(rotsdf)
    colnames(rotsdf) <- c("logFC", "gene", "pvalue", "FDR")
    rotsdf <- rotsdf[,c("logFC", "pvalue", "FDR","gene")]
    colnames(rotsdf) <- c("logFC","pvalue","padj","gene")
    return(rotsdf)
    
  }
  
  else if (test=="Limma")
  {
    library(edgeR)
    library(limma)
    
    group <- apply(as.matrix(table(individual,group)),1,function(x) x[1]==0)
    group_names <- names(group)
    group <- as.numeric(group)
    names(group) <- group_names
    group <- group[colnames(pb)]
    group <- plyr::mapvalues(group,from = c(0,1),to = c("A","B"))
    
    
    #Make the model and contrast for statistical testing
    mm_highcells <- model.matrix(~0 + factor(as.character(group)))
    dimnames(mm_highcells) <- list(names(group), levels(factor(as.character(group))))
    mm_highcells <- mm_highcells[colnames(pb),]
    contrast_highcells <- makeContrasts("B-A", levels = mm_highcells)
    
    if (sum.or.mean == "sum")
    {
      dge <- DGEList(counts = assay(pb), remove.zeros = T)
      dge <- calcNormFactors(dge)
      v <- voom(dge, design = mm_highcells, plot = F)
      fit <- lmFit(v, design = mm_highcells)
      fit2 <- contrasts.fit(fit, contrasts = contrast_highcells)
      fit2 <- eBayes(fit2)
      limma_out <- topTable(fit2, number = nrow(dge))
      
    } else {
      v <- assay(pb)
      v <- v[apply(v,1,sum)!=0,]
      fit <- lmFit(v, design = mm_highcells)
      fit2 <- contrasts.fit(fit, contrasts = contrast_highcells)
      fit2 <- eBayes(fit2)
      limma_out <- topTable(fit2, number = nrow(v))
      
    }
    df <- limma_out[,c("logFC","P.Value","adj.P.Val")]
    df$gene <- rownames(limma_out)
    colnames(df) <- c("logFC","pvalue","padj","gene")
    
    return(df)
  }
  
  
  else if (test=="edgeR")
  {
    library(edgeR)
    library(limma)
    
    group <- apply(as.matrix(table(individual,group)),1,function(x) x[1]==0)
    group_names <- names(group)
    group <- as.numeric(group)
    names(group) <- group_names
    group <- group[colnames(pb)]
    group <- plyr::mapvalues(group,from = c(0,1),to = c("A","B"))
    
    
    #Make the model and contrast for statistical testing
    mm_highcells <- model.matrix(~0 + factor(as.character(group)))
    dimnames(mm_highcells) <- list(names(group), levels(factor(as.character(group))))
    mm_highcells <- mm_highcells[colnames(pb),]
    contrast_highcells <- makeContrasts("B-A", levels = mm_highcells)
    
    
    dge <- DGEList(counts = assay(pb), group = group, remove.zeros = T)
    dge <- calcNormFactors(dge)
    dge <- estimateDisp(dge, design = mm_highcells)
    fit <- glmQLFit(dge, design = mm_highcells)
    fit2 <- glmQLFTest(fit, contrast = contrast_highcells)
    tt <- topTags(fit2, n = nrow(dge))
    edger_out <- tt$table
    
    df <- edger_out[,c("logFC","PValue","FDR")]
    df$gene <- rownames(edger_out)
    colnames(df) <- c("logFC","pvalue","padj","gene")
    
    
    return(df)
  }
  
  
  else if (test=="DESeq2")
  {
    library(limma)
    library(edgeR)
    library(DESeq2)
    
    group <- apply(as.matrix(table(individual,group)),1,function(x) x[1]==0)
    group_names <- names(group)
    group <- as.numeric(group)
    names(group) <- group_names
    group <- group[colnames(pb)]
    group <- plyr::mapvalues(group,from = c(0,1),to = c("A","B"))
    
    
    #Make the model and contrast for statistical testing
    mm_highcells <- model.matrix(~0 + factor(as.character(group)))
    dimnames(mm_highcells) <- list(names(group), levels(factor(as.character(group))))
    mm_highcells <- mm_highcells[colnames(pb),]
    contrast_highcells <- makeContrasts("B-A", levels = mm_highcells)
    
    pb_counts <- assay(pb)
    mode(pb_counts) <- "integer"
    dge <- DESeqDataSetFromMatrix(pb_counts, colData = colData(pb), design = mm_highcells)
    dds <- DESeq(dge)
    res <- results(dds, contrast = contrast_highcells)
    deseq_out <- as.data.frame(res@listData, row.names = res@rownames)
    
    df <- deseq_out[,c("log2FoldChange","pvalue","padj")]
    df$gene <- rownames(deseq_out)
    colnames(df) <- c("logFC","pvalue","padj","gene")
    return(df)
  }
  
  
}



# replace names with the ones you have in your Seurat object
obj.nsf.cca2=obj.nsf.cca1
individual <- obj.nsf.cca2@meta.data$sample
group <- obj.nsf.cca2@meta.data$Type
raw_data <- obj.nsf.cca2@assays$RNA@layers$counts
normalized_data <- obj.nsf.cca2@assays$RNA@layers$data
clustering <- obj.nsf.cca2@meta.data$seurat_clusters

# Run pseudo-bulk sum aggregation and testing with ROTS for each cluster separately

pb_rots_out_list <- list()
for (cluster_name in names(table(clustering)))
{
  pb_rots_out <- RunPseudobulkMethod(raw.data = raw_data[,clustering==cluster_name],
                                     normalized.data = normalized_data[,clustering==cluster_name],
                                     individual = individual[clustering==cluster_name],
                                     group = group[clustering==cluster_name],
                                     test = "ROTS",  ###Limma,DESeq2,edgeR
                                     sum.or.mean = "sum")
  pb_rots_out_list[[cluster_name]] <- pb_rots_out
  
}

# Get results with FDR <= 0.05 as a separate list

pb_rots_signif_list <- lapply(pb_rots_out_list,function(x) x[x$padj<=0.05,])

save(pb_rots_out_list,pb_rots_signif_list,file="results_pseudobulk_ROTS_sum.RData")




