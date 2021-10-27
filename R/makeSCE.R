makeSCE<-function (obj, cData=NULL, batch = NULL) #pseudo.cells=FALSE 
{
  if (class(obj$counts)[1] != "matrix") {
    obj$counts <- as.matrix(obj$counts)
  }
#  cData$logLS <- log10(colSums(obj$counts))
  batch<-obj$batch
  data.norm <- get.res(obj, type = "logcounts", batch = batch)
  sce.obj <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = obj$counts,
                                                                      logcounts = data.norm), 
                                                        colData = cData)
  assays(sce.obj, withDimnames=FALSE)$pearson <- get.res(obj, type = "pearson", 
                                                         batch = batch)
  
  assays(sce.obj, withDimnames=FALSE)$logcorrected <- log(get.res(obj, type = "quantile", 
                                                                  batch = batch) + 1)
  
  sce.obj
  }
