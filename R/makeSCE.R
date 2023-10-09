#' Create SingleCellExperiment object containing normalized data
#'
#' This function takes ruvIII.nb or fastruvIII.nb output as input and creates a SingleCellExperiment (SCE) object containing various metrics of normalized data.  
#'
#' @param obj  object containing output of call to ruvIII.nb or fastruvIII.nb function.
#' @param cData A data frame containing cell-level metadata. This data frame will be used as 'colData' in the SingleCellExperiment object and must be specified.
#' @param batch numeric vector containing batch information for each sample. Must correspond to columns of count matrix. Only needed if batch-specific dispersion parameter is fitted.
#' @param assays type of corrected data to return. Available types are "logcounts","logPAC" and "pearson".

#' @return A SingleCellExperiment object with normalized data added to the 'assays' slot. The log percentile-adjusted count is stored in the 'logPAC' component of the 'assays' slot and the Pearson residuals is in the 'pearson' component.

makeSCE<-function (obj, cData=NULL, batch = NULL, assays = c('pearson','logPAC'))
{
  batch<-obj$batch
  if(is.null(cData))
    stop("'cData' must be a user-specified data frame.")
  if(any(! assays %in% c('logcounts','logPAC','pearson')))
    stop("'Wrong names of corrected assays data was specified. Available corrected data assays are logcounts,logPAC and pearson")

  sce.obj <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = as(obj$counts,"sparseMatrix")), colData = cData)
  if(any(assays=='logcounts'))
   assays(sce.obj, withDimnames=FALSE)$logcounts <- tryCatch({as( get.res(obj,type='logcounts',batch = batch),"sparseMatrix" )}, error = function(e) {NULL})
  if(any(assays=='pearson'))
   assays(sce.obj, withDimnames=FALSE)$pearson <- tryCatch({as( get.res(obj,type='pearson',batch = batch), "sparseMatrix")}, error = function(e) {NULL})
  if(any(assays=='logPAC'))
   assays(sce.obj, withDimnames=FALSE)$logPAC <- tryCatch({as(  log(get.res(obj,type='quantile',batch = batch)+1), "sparseMatrix")}, error = function(e) {NULL})

  sce.obj
}
