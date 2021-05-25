#' Add metrics of normalized data to SingleCellExperiment object
#'
#' This function takes ruvIIInb output and creates SingleCellExperiment object that contains various metrics of normalized data.
#' @param obj An ruvIII.nb output. 
#' @param cData A data frame containing cell-level metadata.
#' @param batch numeric vector containing batch information for each cell. Must corresponds to columns of count matrix. Only needed if batch-specific dispersion parameter is fitted.
#' @return A SingleCellExperiment object with metrics of normalized data in the assay slot. The \code{pearson} component of the assay slot contains the pearson residuals and the \code{logcorrected} component contains the log of percentile-invariant adjusted data + 1.
#' @export
makeSCE <- function(obj,cData,batch=NULL) {
 if(class(obj$counts)!='matrix') {
   obj$counts <- as.matrix(obj$counts)
 }
 cData$logLS<- log10(colSums(obj$counts))
 data.norm <- get.res(obj,type='logcounts',batch=batch)
 sce.obj   <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=obj$counts,logcounts=data.norm),colData=cData)
 assays(sce.obj)$pearson      <-  get.res(obj,type='pearson',batch=batch)
 assays(sce.obj)$logcorrected <- log(get.res(obj,type='quantile',batch=batch)+1)
 sce.obj
}


