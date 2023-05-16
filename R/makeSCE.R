#' Create SingleCellExperiment object containing normalized data
#'
#' This function takes ruvIII.nb or fastruvIII.nb output as input and creates a SingleCellExperiment (SCE) object containing various metrics of normalized data.  
#'
#' @param obj  object containing output of call to ruvIII.nb or fastruvIII.nb function.
#' @param cData A data frame containing cell-level metadata. This data frame will be used as 'colData' in the SingleCellExperiment object and must be specified.
#' @param batch numeric vector containing batch information for each sample. Must correspond to columns of count matrix. Only needed if batch-specific dispersion parameter is fitted.

#' @return A SingleCellExperiment object with normalized data added to the 'assays' slot. The log percentile-adjusted count is stored in the 'logPAC' component of the 'assays' slot and the Pearson residuals is in the 'pearson' component.

makeSCE<-function (obj, cData=NULL, batch = NULL)
{
  batch<-obj$batch
  if(is.null(cData))
    stop("'cData' must be a user-specified data frame.")

  data.norm <- get.res(obj,type='logcounts',batch = batch)
  gc(verbose=FALSE)
  sce.obj <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = obj$counts,
                                                                      logcounts = data.norm), 
                                                        colData = cData)
  assays(sce.obj, withDimnames=FALSE)$pearson <- get.res(obj,type='pearson',batch = batch)
  assays(sce.obj, withDimnames=FALSE)$logPAC <- log(get.res(obj,type='quantile',batch = batch)+1)

  sce.obj
  }
