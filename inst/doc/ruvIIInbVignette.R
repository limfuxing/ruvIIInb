## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
knitr::opts_chunk$set(fig.width=9, fig.height=6)


## ---- include=FALSE, resuls='hide', warning=FALSE, message=FALSE--------------
rm(list=ls())

## ---- include=TRUE, resuls='hide', warning=FALSE, message=FALSE,eval=FALSE----
#  #Install ruvIIInb package
#  devtools::install_github("limfuxing/ruvIIInb",build_vignettes = TRUE)
#  library(ruvIIInb)
#  

## ---- include=TRUE, resuls='hide', warning=FALSE, message=FALSE---------------

library(SingleCellExperiment)
library(scater)
library(scran)
library(scuttle)
library(edgeR)
library(SingleR)
library(celldex)
library(hrbrthemes)
library(tidyverse)
library(ggplot2)
library(uwot)
library(scMerge)
library(Seurat)
library(randomcoloR)
library(dittoSeq)



## ---- include=FALSE, echo=FALSE,results='hide',warning=FALSE,eval=FALSE-------
#  
#  ###########################################################
#  # read the data
#  data <- edgeR::read10X(path='/stornext/Home/data/allstaff/d/delivera.a/NSCLC_ADL_wehi/filtered_gene_bc_matrices')
#  dim(data)
#  
#  ###########################################################
#  # change gene names from ENSEMBL to SYMBOL
#  genenames <- unlist(data$genes)
#  head(genenames)
#  length(genenames)
#  dim(data$counts)
#  #Only keeping non-duplicated genes
#  data$counts <- data$counts[!duplicated(genenames),]
#  genenames <- genenames[!duplicated(genenames)]
#  names(genenames) <- NULL
#  # create SCE object
#  nsclc_obs <- data$counts
#  dim(nsclc_obs)
#  rownames(nsclc_obs) <- genenames
#  nsclc_samples<-data$samples
#  head(nsclc_samples)
#  
#  dim(nsclc_samples)
#  nsclc_samples<-nsclc_samples[,c("lib.size","Barcode")]
#  head(nsclc_samples)
#  
#  saveRDS(nsclc_obs,'nsclc_obs.rds')
#  saveRDS(nsclc_samples,'nsclc_samples.rds')
#  

## ---- include=FALSE, echo=FALSE, results='hide',warning=FALSE-----------------
nsclc_obs<-readRDS('../inst/extdata/nsclc_obs.rds')
nsclc_samples<-readRDS('../inst/extdata/nsclc_samples.rds')


## -----------------------------------------------------------------------------
sce <- SingleCellExperiment(assays=list(counts=nsclc_obs),
                            colData=as.data.frame(nsclc_samples))


## -----------------------------------------------------------------------------
sce <- addPerCellQCMetrics(x = sce,subsets=list(Mito=grep("MT-",rownames(sce))))

## ----warning=FALSE------------------------------------------------------------
libsize_drop <- isOutlier(
  metric = sce$total, 
  nmads = 2,
  type = "lower", 
  log = TRUE)
colData(sce)$libsize_drop<-libsize_drop


mito_drop <- isOutlier(
  metric = colData(sce)$subsets_Mito_percent, 
  nmads = 3, 
  type = "higher")
colData(sce)$mito_drop<-mito_drop

## ----LibsizeFig, fig.cap='A histogram showing the distribution of log cell-wise total counts flagging those with low library size', echo=FALSE,warning=FALSE----
plot_df <- data.frame(logtotal=log(sce$total), 
                      libsize_drop=factor(libsize_drop),
                      mito_drop=factor(mito_drop),
                      logdetected=log(sce$detected))
plot_df$mito_drop<-relevel(plot_df$mito_drop, "TRUE")
plot_df$libsize_drop<-relevel(plot_df$libsize_drop, "TRUE")

#A histogram showing the distribution of flagged cells
p <- plot_df %>%
  ggplot( aes(x=logtotal, fill=libsize_drop)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity',bins=30) +
  scale_fill_manual(values=c( "#69b3a2","#404080")) +
  theme_ipsum() +
  labs(fill="") +xlab("Cell-level log total count")+ylab("Frequency")
p+ggtitle('The distribution of cell-level log total counts\n  flagging cells with low library size')


## ----MitoFig, fig.cap='A scatterplot showing cell-level log number of detected genes against log total count', echo=FALSE, warning=FALSE----

#A plot showing the number of detected genes vs total cell count
p2<- plot_df %>%
  arrange(desc(mito_drop)) %>%
  ggplot( aes(x=logtotal,y=logdetected, fill=mito_drop)) +
  geom_point(
  mapping = aes(colour = mito_drop, shape = libsize_drop),
size = 2,alpha = 5 / 6) + scale_color_manual(values=c("#69b3a2","#404080"))+
  xlab("Log total cell-level count")+ylab("Log number of cell-level detected genes")
p2+ggtitle('The cell-level number of detected genes vs total cell-level count (on log scale)')


## -----------------------------------------------------------------------------
sce <- addPerFeatureQCMetrics(x = sce)

#Remove genes with zero counts for each gene
sce<-sce[-which(rowData(sce)$mean==0),]

lowcount_drop <- isOutlier(
  metric = rowData(sce)$mean,
  nmads = 1.5, 
  type = "lower", 
  log = TRUE)


## ----LowCellCountFig, fig.cap='A histogram showing the distribution of\n gene-wise log mean count', echo=FALSE,warning=FALSE----
#mean count of genes across all cells
plot_df2 <- data.frame(mean_genecount=log(rowData(sce)$mean), lowcount_drop=factor(lowcount_drop))
p <- plot_df2 %>%
  arrange(desc(lowcount_drop)) %>%
  ggplot( aes(x=mean_genecount, fill=lowcount_drop)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity',bins=50) +
  scale_fill_manual(values=c("#404080", "#69b3a2")) +
  theme_ipsum() +
  labs(fill="") +xlab("Log mean count across all cells")+ylab("Frequency")
p+ggtitle('The distribution of log mean count across all cells,\n flagging those with low count')



## -----------------------------------------------------------------------------
sce <- sce[!(lowcount_drop), !(libsize_drop | mito_drop)]

## ----warning=FALSE,message=FALSE,include=FALSE--------------------------------
sce<-readRDS('../inst/extdata/NSCLC_pre_27_09_2021.rds') 

## ----warning=FALSE,message=FALSE,eval=FALSE-----------------------------------
#  sce <- computeSumFactors(sce,assay.type="counts")
#  data_norm_pre <- sweep(assays(sce)$counts,2,sce$sizeFactor,'/')
#  assays(sce, withDimnames=FALSE)$lognormcounts<- log(data_norm_pre+1)

## ----message=FALSE,warning=FALSE,, include=FALSE, eval=FALSE------------------
#  #Saving the pre-NSCLC data file as computeSumFactors() takes a long time
#  #saveRDS(sce,file='NSCLC_pre_27_09_2021.rds')

## ----warning=FALSE,message=FALSE----------------------------------------------
hvg_df <- modelGeneVar(sce, assay.type = "lognormcounts")
hvg_df<- hvg_df[rownames(sce),]
hvg_genes   <- hvg_df$bio>quantile(hvg_df$bio,prob=0.75)
rowData(sce)$hvg_genes<-hvg_genes
sce<-sce[hvg_genes,]

## ----VarMeanFig, fig.cap='Gene-level variance vs Mean log-expression', echo=FALSE,warning=FALSE----
plot(hvg_df$mean, hvg_df$total,
     xlab="Mean log-expression", 
     ylab="Variance",pch=16,col="grey", main="Gene-level variance vs mean log-expression\n with the HVGs shown in red")
curve(metadata(hvg_df)$trend(x), col="blue", add=TRUE)
points(hvg_df$mean[hvg_genes], hvg_df$total[hvg_genes],col="red",pch=16)

## ----message=FALSE------------------------------------------------------------
# Reading in the control genes
ctl <- read.csv('../data/scHK_human.csv', header=T)[,1]

## ----LSFig, fig.cap='Comparing the overall log library size vs the log library size using the control genes', echo=FALSE,warning=FALSE----
#compare overall log Library sizes vs log library sizes using control genes only
sce$logLS <- log(colSums(assays(sce)$counts))
rowData(sce)$ctlLogical<-rownames(assays(sce)$counts) %in% ctl
sce$logLS_ctl <- log(colSums(assays(sce)$counts[rowData(sce)$ctlLogical,]))

plot(sce$logLS,sce$logLS_ctl,col="grey",
     xlab="Log library size using all genes",
     ylab="Log library size using the control genes")

## ----message=FALSE,warning=FALSE,  include=FALSE,eval=FALSE-------------------
#  assays(sce, withDimnames=FALSE)

## ----message=FALSE,warning=FALSE, eval=FALSE----------------------------------
#  
#  # Perform initial clustering to identify pseudo-replicates
#  snn_gr_init <- buildSNNGraph(sce, assay.type = "lognormcounts")
#  clusters_init <- igraph::cluster_louvain(snn_gr_init)
#  sce$cluster_init <- factor(clusters_init$membership)
#  
#  # Identify cell-type for each cluster
#  ref_se <- HumanPrimaryCellAtlasData()
#  
#  # Annotate initial cluster
#  pred10xnsclc_init  <- SingleR(test = sce, ref = ref_se,
#                                labels = ref_se$label.fine,
#                                assay.type.test='lognormcounts',
#                                clusters=factor(sce$cluster_init))
#  sce$ctype_init <- pred10xnsclc_init$labels[sce$cluster_init]
#  
#  # Anotate cells
#  pred10xnsclc_init2  <- SingleR(test = sce, ref = ref_se,
#                                labels=ref_se$label.fine,
#                                assay.type.test='lognormcounts')
#  sce$ctype_init2 <- pred10xnsclc_init2$labels
#  
#  
#  sce <- runUMAP(sce,exprs_values = "lognormcounts")
#  sce <- runPCA(sce,exprs_values = "lognormcounts")
#  
#  

## ----message=FALSE,warning=FALSE,, include=FALSE, eval=FALSE------------------
#  #Undo this to save this file for future use
#  #Saving the preproccessed NSCLC data file
#  #saveRDS(sce,file='NSCLC_preprocessed_init_annonated_27_09_2021.rds')

## ----message=FALSE,warning=FALSE, include=FALSE-------------------------------
#Saving the preproccessed NSCLC data file 
rm(list=ls())
library(ruvIIInb)
library(SingleCellExperiment)
library(scater)
library(scran)
library(scuttle)
library(edgeR)
require(SingleR)
library(celldex)
library(hrbrthemes)
library(tidyverse)
library(ggplot2)
library(uwot)
library(pheatmap)
library(SingleR)
library(dittoSeq)



sce<-readRDS('../inst/extdata/NSCLC_preprocessed_init_annonated_27_09_2021.rds') 

## ----ClFig, fig.cap='Visualising the initial clusters for identifying pseudo-replicates', echo=FALSE,warning=FALSE----
p1 <- plotUMAP(sce,colour_by='logLS')
p3 <- dittoDimPlot(sce, "ctype_init",colors=c(1,2,3,4,5,6,7,8))
gridExtra::grid.arrange(p1,p3,ncol=2)


## ----message=FALSE,eval=FALSE-------------------------------------------------
#  # Construct the replicate matrix M using pseudo-replicates identified using initial clustering
#  M <- matrix(0,ncol(assays(sce)$counts),length(unique(sce$cluster_init)))
#  cl<- sort(unique(as.numeric(unique(sce$cluster_init))))
#  for(CL in cl)
#    M[which(as.numeric(sce$cluster_init)==CL),CL] <- 1
#  
#  #RUV-III-NB code
#  ruv3nb_out<-ruvIII.nb(Y=as.matrix(assays(sce)$counts), # count matrix with genes as rows and cells as columns
#                      M=M, #Replicate matrix constructed as above
#                      ctl=rowData(sce)$ctlLogical, #A vector denoting control genes
#                      k=2, # dimension of unwanted variation factors
#                      strata=sce$cluster_init)
#  

## ----message=FALSE,eval=FALSE,include=FALSE-----------------------------------
#  which(rowData(sce)$mean==0)
#  which(rowSums(as.matrix(assays(sce)$counts))==0)
#  which(colSums(as.matrix(assays(sce)$counts))==0)

## ----message=FALSE,eval=FALSE-------------------------------------------------
#  
#  ruv3nb_out_ps<-ruvIII.nb(Y=as.matrix(assays(sce)$counts),
#                      M=diag(dim(sce)[2]),
#                      ctl=rowData(sce)$ctlLogical,
#                      k=2,
#                      use.pseudosample=TRUE, #Indicating whether pseudo-cells are to be used to to define pseudo-replicates
#                      batch=rep(1,dim(sce)[2]),# NSCLC data comes from a single batch
#                      strata=sce$cluster_init #The cells which belong to initial clusters are not pooled together
#                      )

## ----message=FALSE,include=FALSE----------------------------------------------
ruv3nb_out<- readRDS('../inst/extdata/NSCLC_ruvIIInb_ADL3_new0.1_27092021_K2_clM.rds')
# ruv3nb_out<- readRDS('../inst/extdata/NSCLC_ruvIIInb_ADL3_new0.1_18092021_K2_pcM.rds')

## ----WFig, fig.cap='Scatterplots showing the log library size vs components of unwanted variation', echo=FALSE,warning=FALSE----
W1data<-data.frame(W1=ruv3nb_out$W[,1],
                   W2=ruv3nb_out$W[,2],LogLS=colData(sce)$logLS)
p1<-ggplot(W1data, aes(x=LogLS, y=W1, color=LogLS)) + geom_point()+ theme(legend.position = "none") 
p2<-ggplot(W1data, aes(x=W2, y=W1, color=LogLS)) + geom_point()
gridExtra::grid.arrange(p1,p2,ncol=2)

## ----message=FALSE,include=FALSE----------------------------------------------
sce_ruv3nb<-readRDS('../inst/extdata/NSCLC_ruvIIInb_ADL3_new0.1_27092021_K2_clM_annotated.rds')

## ----message=FALSE,eval=FALSE-------------------------------------------------
#  
#  #Creating a SingleCellExperiment object
#  sce_ruv3nb <- makeSCE(ruv3nb_out,cData=colData(sce))

## ----message=FALSE,warning=FALSE----------------------------------------------
seurat_ruv3nb <- as.Seurat(sce_ruv3nb, counts = "counts", data ="logcorrected")

## ----message=FALSE,eval=FALSE,include=FALSE,eval=FALSE------------------------
#  # DimPlot(seurat_ruv3nb, reduction = "PCA", group.by = "ctype_init2") + NoLegend()
#  # p1 <- DimPlot(seurat_ruv3nb, reduction = "PCA", group.by = "Source") + NoLegend()
#  # p2 <- RidgePlot(manno.seurat, features = "ACTB", group.by = "Source")
#  # CombinePlots(plots = list(p1, p2))
#  #
#  

## ----eval=FALSE---------------------------------------------------------------
#  # PCA of RUV-III-NB normalized data
#  sce_ruv3nb <- runPCA(sce_ruv3nb, exprs_values = "logcorrected")
#  

## ----PCARUVFig, fig.cap='Visualising PCA following RUV-III-NB normalisation', echo=FALSE,warning=FALSE----
# PCA visualization
p1 <- plotPCA(sce, ncomponents=2,colour_by = "logLS",point_alpha=0.5, point_size=1) + ggtitle('Scran normalized')+ theme(legend.position = "none") 

p2 <- plotPCA(sce_ruv3nb, ncomponents=2,
    colour_by = "logLS",point_alpha=0.5, point_size=1) +  ggtitle('RUV-III-NB normalized')

gridExtra::grid.arrange(p1,p2,ncol=2)


## ----message=FALSE,eval=FALSE-------------------------------------------------
#  
#  # Perform clustering
#  snn_gr <- buildSNNGraph(sce_ruv3nb, assay.type = "logcorrected")
#  clusters <- igraph::cluster_louvain(snn_gr)
#  sce_ruv3nb$cluster_ruv <- factor(clusters$membership)
#  
#  # Identify cell-type for each cluster
#  ref_se <- HumanPrimaryCellAtlasData()
#  
#  # annotate cluster
#  pred10xnsclc_ruv  <- SingleR(test = sce_ruv3nb, ref = ref_se,
#                               labels = ref_se$label.fine,
#                               assay.type.test='logcorrected',
#                               clusters=factor(sce_ruv3nb$cluster_ruv))
#  sce_ruv3nb$ctype_ruv <- pred10xnsclc_ruv$labels[sce_ruv3nb$cluster_ruv]
#  
#  # anotate cells
#  pred10xnsclc_ruv2  <- SingleR(test = sce_ruv3nb, ref = ref_se,
#                                labels=ref_se$label.fine,
#                                assay.type.test='logcorrected')
#  sce_ruv3nb$ctype_ruv2 <- pred10xnsclc_ruv2$labels
#  
#  #UMAP
#  sce_ruv3nb <- runUMAP(sce_ruv3nb, exprs_values = "logcorrected")

## ----message=FALSE,include=FALSE,eval=FALSE-----------------------------------
#  #Saving the preproccessed AND ANNONATED NSCLC data file
#  saveRDS(sce_ruv3nb,file='NSCLC_ruvIIInb_ADL3_new0.1_27092021_K2_clM_annotated.rds')
#  

## ----ClRUVFig, fig.cap='Visualising the clusters following RUV-III-NB normalisation', echo=FALSE,warning=FALSE----
p1 <- plotUMAP(sce_ruv3nb,colour_by='logLS')+ggtitle('RUV-III-NB normalized')
p3<-dittoDimPlot(sce_ruv3nb, "ctype_ruv",colors=c(1,3,4,10,5,9, 6,11,7))
p1s <- plotUMAP(sce,colour_by='logLS')+ggtitle('Scran normalized')
p3s<-dittoDimPlot(sce, "ctype_init",colors=c(1,2,3,4,5,6,7,8))
gridExtra::grid.arrange(p1s,p1,p3s,p3,ncol=2)


## ----message=FALSE,warnings=TRUE----------------------------------------------
markers_ruv3nb <- findMarkers(x=as.matrix(assays(sce_ruv3nb)$logcorrected),
                              groups=sce_ruv3nb$ctype_ruv)#identify a combination of marker genes that together defines one cluster against the rest. 
cluster3_markers <- markers_ruv3nb[[3]] #Cluster 3 ("CMP" cell-type) is used here for demonstration purposes.
cluster3_top5 <- cluster3_markers[cluster3_markers $Top <= 5,] #Collects the top DE genes from each pairwise comparison involving cluster 3 
logFCs_cluster3_top5 <- getMarkerEffects(cluster3_top5) #Extracts log fold changes

## ----heatmapRUVFig, fig.cap='A heatmap of the DE log fold changes', echo=FALSE,warning=FALSE----
pheatmap(logFCs_cluster3_top5, breaks=seq(-5, 5, length.out=101))


## ----include=FALSE,eval=FALSE-------------------------------------------------
#  #Downloading after saving the output to save time
#  rm(list=ls())
#  sce_ruv3nb_celline<-readRDS('../inst/extdata/celline_ruvIIInb_ADL0.1_18092021_K2_added.rds')
#  #new.labels is misspelt
#  names_ruv3nb_celline<-names(colData(sce_ruv3nb_celline))
#  names(colData(sce_ruv3nb_celline))[which(names_ruv3nb_celline=="new.lables")]<-"new.labels"
#  names(colData(sce_ruv3nb_celline))
#  sce_celline<-readRDS('../inst/extdata/sce_celline_added.rds')
#  names_celline<-names(colData(sce_celline))
#  names(colData(sce_celline))[which(names_celline=="new.lables")]<-"new.labels"
#  names(colData(sce_celline))

## ----include=FALSE------------------------------------------------------------
rm(list=ls())
sce_celline<- readRDS('../inst/extdata/celline_preproc_ADL30.1_18092021.rds')

#new.labels is misspelt
names_celline<-names(colData(sce_celline))
names(colData(sce_celline))[which(names_celline=="new.lables")]<-"new.labels"
names(colData(sce_celline))

celline_obs<-as.matrix(assays(sce_celline)$counts)
celline_samples<-colData(sce_celline)
celline_rows<-rowData(sce_celline)


## -----------------------------------------------------------------------------
sce_celline <- SingleCellExperiment(assays=list(counts=celline_obs),
                            colData=as.data.frame(celline_samples),
                            rowData=as.data.frame(celline_rows))


## ----include=FALSE------------------------------------------------------------
names_celline<-names(colData(sce_celline))
names(colData(sce_celline))[which(names_celline=="new.lables")]<-"new.labels"
names(colData(sce_celline))

## ----message=FALSE,eval=FALSE-------------------------------------------------
#  
#  # Construct the replicate matrix M using the known cell-types
#  unique_ctype<- unique(sce_celline$new.lables)
#  M <- matrix(0,ncol(assays(sce_celline)$counts),length(unique_ctype))
#  dim(M)
#  for(i in 1:length(unique_ctype))
#   M[sce_celline$new.lables==unique_ctype[i],i] <- 1
#  
#  
#  #RUV-III-NB code
#  ruv3nb_out_celline<-ruvIII.nb(Y=as.matrix(assays(sce_celline)$counts),
#                             M=M,
#                             ctl=ctl,
#                             k=2)
#  

## ----message=FALSE,include=FALSE----------------------------------------------
#Loading the pre-saved RUV-III-NB sce output 
ruv3nb_out_celline<- readRDS('../inst/extdata/celline_ruvIIInb_ADL0.1_18092021_K2.rds')

## ----message=FALSE,warning=FALSE----------------------------------------------

#Creating a SingleCellExperiment object
sce_ruv3nb_celline <- makeSCE(ruv3nb_out_celline,cData=colData(sce_celline))

#Converting to a Seurat object if necessary
seurat_ruv3nb_celline <- as.Seurat(sce_ruv3nb_celline, counts = "counts", data ="logcorrected")


## ----WFig2, fig.cap='Scatterplots showing batch, log library size vs the components of unwanted variation', echo=FALSE,warning=FALSE----
Wdata_celline<-data.frame(W1=ruv3nb_out_celline$W[,1],
                   W2=ruv3nb_out_celline$W[,2],
                   batch=colData(sce_celline)$batch,
                   LogLS=log(colSums(assays(sce_celline)$counts)))
p1<-ggplot(Wdata_celline, aes(x=LogLS, y=W1, color=batch)) + geom_point()+ theme(legend.position = "none") 
p2<-ggplot(Wdata_celline, aes(x=W2, y=W1, color=batch)) + geom_point()
gridExtra::grid.arrange(p1,p2,ncol=2)

## -----------------------------------------------------------------------------
# PCA and UMAP of RUV-III-NB normalized data 
sce_ruv3nb_celline <- runPCA(sce_ruv3nb_celline, exprs_values = "logcorrected")
sce_ruv3nb_celline <- runUMAP(sce_ruv3nb_celline, exprs_values = "logcorrected")

## -----------------------------------------------------------------------------
# PCA and UMAP of Scran normalized data 
sce_celline <- computeSumFactors(sce_celline,assay.type="counts")
data_norm_celline_pre <- sweep(assays(sce_celline)$counts,2,sce_celline$sizeFactor,'/')
assays(sce_celline, withDimnames=FALSE)$lognormcounts<- log(data_norm_celline_pre+1)
sce_celline <- runPCA(sce_celline,exprs_values = "lognormcounts")
sce_celline <- runUMAP(sce_celline,exprs_values = "lognormcounts")

## ----message=FALSE,warning=FALSE, include=FALSE, eval=FALSE-------------------
#  #save the files here for download
#  saveRDS(sce_ruv3nb_celline,file='celline_ruvIIInb_ADL0.1_18092021_K2_added.rds')
#  saveRDS(sce_celline,file='sce_celline_added.rds')
#  

## ----PCARUVFigcelline, fig.cap='Visualising PCA following RUV-III-NB normalisation of the Celline data', echo=FALSE,warning=FALSE----

# PCA visualization
p1 <- plotPCA(sce_celline, ncomponents=2,colour_by = "new.labels",point_alpha=0.5, point_size=1) + ggtitle('Scran normalized')

p2 <- plotPCA(sce_ruv3nb_celline, ncomponents=2,
              colour_by = "new.labels",point_alpha=0.5, point_size=1) +  ggtitle('RUV-III-NB normalized')

gridExtra::grid.arrange(p1,p2,ncol=2)


## ----ClRUVFigcelline, fig.cap='Visualising the clusters following RUV-III-NB normalisation of the Celline data', echo=FALSE,warning=FALSE----

p1 <- plotUMAP(sce_ruv3nb_celline,colour_by='new.labels')+ggtitle('RUV-III-NB normalized')
#p2 <- plotUMAP(sce_ruv3nb,colour_by='cluster10X')
p3 <- plotUMAP(sce_ruv3nb_celline,colour_by='new.labels')+ggtitle('RUV-III-NB normalized')
p1s <- plotUMAP(sce_celline,colour_by='batch')+ggtitle('Scran normalized')
p3s <- plotUMAP(sce_celline,colour_by='batch')+ggtitle('Scran normalized')


gridExtra::grid.arrange(p1s,p1,p3s,p3,ncol=2)


## ----message=FALSE,warnings=TRUE----------------------------------------------
markers_ruv3nb_celline <- findMarkers(x=as.matrix(assays(sce_ruv3nb_celline)$logcorrected),
                              groups=sce_ruv3nb_celline$new.labels)
markers_celline293t <- markers_ruv3nb_celline[["293t"]]
top5_celline293t <- markers_celline293t[markers_celline293t$Top <= 5,]
top5_celline293t

## ----message=FALSE,eval=FALSE-------------------------------------------------
#  # Construct the replicate matrix M using the `scReplicate` algorithm
#  sce_celline <- computeSumFactors(sce_celline,assay.type="counts")
#  data_norm_pre <- sweep(assays(sce_celline)$counts,2,sce_celline$sizeFactor,'/')
#  assays(sce_celline, withDimnames=FALSE)$lognormcounts<- log(data_norm_pre+1)
#  scm.M<-scReplicate(sce_combine=sce_celline,
#                     exprs = "lognormcounts",
#                     batch=as.numeric(factor(sce_celline$batch)),
#                     kmeansK=c(1,1,2))
#  mode(scm.M) <- 'logical'
#  
#  
#  ruv3nb_out_celline_scMerge<- ruvIIInb::ruvIII.nb(as.matrix(assays(sce_celline)$counts),
#                              M=scm.M,
#                              ctl=ctl,
#                              k=2)

## ----message=FALSE,include=FALSE,eval=FALSE-----------------------------------
#  #SAVED after adding pca, umap etc to save time
#  sce_ruv3nb_celline_scMerge<- readRDS('../inst/extdata/celline_ruvIIInb_ADL0.1_18092021_scMerge_K2_added.rds')

## ----message=FALSE,include=FALSE----------------------------------------------
ruv3nb_out_celline_scMerge<- readRDS('../inst/extdata/celline_ruvIIInb_ADL0.1_18092021_scMerge_K2.rds')

## ----message=FALSE------------------------------------------------------------
#Creating a SingleCellExperiment object
sce_ruv3nb_celline_scMerge <- makeSCE(ruv3nb_out_celline_scMerge,
                                         cData=colData(sce_celline))

#Converting to a Seurat object if necessary
seurat_ruv3nb_celline_scMerge <- as.Seurat(sce_ruv3nb_celline_scMerge, counts = "counts", data ="logcorrected")



## -----------------------------------------------------------------------------

# PCA and UMAP of RUV-III-NB normalized data 
sce_ruv3nb_celline_scMerge <- runPCA(sce_ruv3nb_celline_scMerge, 
                                     exprs_values = "logcorrected")
sce_ruv3nb_celline_scMerge <- runUMAP(sce_ruv3nb_celline_scMerge, 
                                      exprs_values = "logcorrected")

## ----message=FALSE,warning=FALSE, include=FALSE, eval=FALSE-------------------
#  #save the files here for download
#  saveRDS(sce_ruv3nb_celline_scMerge,file='celline_ruvIIInb_ADL0.1_18092021_scMerge_K2_added.rds')

## ----scMergeClRUVFigcelline, fig.cap='Visualising the clusters following RUV-III-NB normalisation of the Celline data using scReplicate', echo=FALSE,warning=FALSE----

p1 <- plotUMAP(sce_ruv3nb_celline_scMerge,colour_by='new.labels')+ggtitle('RUV-III-NB normalized')
p3 <- plotUMAP(sce_ruv3nb_celline_scMerge,colour_by='batch')+ggtitle('RUV-III-NB normalized')
p1s <- plotUMAP(sce_celline,colour_by='new.labels')+ggtitle('Scran normalized')
p3s <- plotUMAP(sce_celline,colour_by='batch')+ggtitle('Scran normalized')


gridExtra::grid.arrange(p1s,p1,p3s,p3,ncol=2)


## ----message=FALSE,warning=FALSE----------------------------------------------
markers_ruv3nb_celline_scMerge <- findMarkers(x=as.matrix(assays(sce_ruv3nb_celline_scMerge)$logcorrected),
                                      groups=sce_ruv3nb_celline_scMerge$new.labels)
markers_celline293t_scMerge <- markers_ruv3nb_celline_scMerge[["293t"]]
top5_celline293t_scMerge <- markers_celline293t_scMerge[markers_celline293t_scMerge$Top <= 5,]
top5_celline293t_scMerge

## -----------------------------------------------------------------------------
sessionInfo()

