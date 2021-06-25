#' Produce normalized count data after adjusting for unwanted variations
#'
#' This function takes ruvIII.nb output as input and produce various metrics of normalized data. 
#' Option for metrics of normalized data includes log normalized count, Pearson and Deviance residuals. 
#' @param out output of call to ruvIII.nb function.
#' @param type type of normalized data metrics. Options include pearson residuals and percentile-invariant log adjusted count.
#' @param batch numeric vector containing batch information for each sample.Must correspond to columns of count matrix. Only needed if batch-specific dispersion parameter is fitted.
#' @return The normalized count data.
#' @export
get.res <- function(out,type='pearson',batch=NULL) {
    Y <- out$counts
    if(class(Y)!='matrix') {
       Y <- as.matrix(Y)
    }

    if(class(out$pi0)!='matrix') {
       out$pi0 <- as.matrix(out$pi0)
    }
    
    mu <- exp(out$a%*%t(out$W) + out$gmean) * (1-out$pi0)
    mu.full <- exp(out$a%*%t(out$W) + out$gmean + out$Mb[,apply(out$M,1,which)]) 
    Wa.mean = out$a %*% as.matrix(colMeans(out$W))

    # smooth outlier estimates
    idx.gene <- abs(log(rowMeans(Y+1)/rowMeans(mu+1)))>log(5)
    if(sum(idx.gene)>1) 
      mu[idx.gene,] <- rowMeans(Y[idx.gene,])
    if(sum(idx.gene)==1) 
      mu[idx.gene,] <- mean(Y[idx.gene,])

    idx.gene <- abs(log(rowMeans(Y+1)/rowMeans(mu.full+1)))>log(5)
    mu.full[idx.gene,] <- Y[idx.gene,] * (1-out$pi0[idx.gene,])

    mu <- mu/(1-out$pi0)
    if(type=='pearson') {
       if(is.null(dim(out$psi)))
        dev <- (Y-mu*(1-out$pi0))/sqrt(mu*(1-out$pi0)*(1+mu*(out$psi.DE+out$pi0)))
       if(!is.null(dim(out$psi)))
        dev <- (Y-mu*(1-out$pi0))/sqrt(mu*(1-out$pi0)*(1+mu*(out$psi.DE[,batch]+out$pi0)))
	lims<- quantile(c(dev),prob=c(0.0005,0.9995),na.rm=T)
	dev[dev< lims[1]] <- lims[1] ; dev[dev>lims[2]] <- lims[2]
    }
    if(type=='logcounts') {
        mu <- exp( out$a%*%t(out$W)) * (1-out$pi0)
        dev <- log(Y/mu+1)
    }
    if(type=='quantile') {
        mu.noUV <- exp(out$Mb[,apply(out$M,1,which)] + out$gmean + c(Wa.mean))
        idx.gene <- abs(log(rowMeans(Y+1)/rowMeans(mu.noUV+1)))>log(5)
        if(sum(idx.gene)>1) 
          mu.noUV[idx.gene,] <- rowMeans(Y[idx.gene,]) 
        if(sum(idx.gene)==1) 
          mu.noUV[idx.gene,] <- mean(Y[idx.gene,]) 

        if(is.null(dim(out$psi))) {
          a <- ZIM::pzinb(Y-1,lambda=mu.full,k=1/out$psi.DE,omega=out$pi0)
          b <- ZIM::pzinb(Y,lambda=mu.full,k=1/out$psi.DE,omega=out$pi0)
          a[a>0.999] <- 0.9985
          b[b>0.999] <- 0.999
          avg.pi0 <- rowMeans(out$pi0)
          p <- dev <- Y
          for(i in 1:ncol(Y)) {
            p[,i]   <- runif(nrow(Y),a[,i],b[,i])
            dev[,i] <- ZIM::qzinb(p[,i],lambda=mu.noUV[,i],k=1/out$psi.DE,omega=avg.pi0)
          }
        }      
        if(!is.null(dim(out$psi))) {
          a <- ZIM::pzinb(Y-1,lambda=mu.full,k=1/out$psi.DE[,batch],omega=out$pi0)
          b <- ZIM::pzinb(Y,lambda=mu.full,k=1/out$psi.DE[,batch],omega=out$pi0)
          a[a>0.999] <- 0.9985
          b[b>0.999] <- 0.999
          avg.pi0 <- rowMeans(out$pi0)
          p <- dev <- Y
          for(i in 1:ncol(Y)) {
            p[,i]   <- runif(nrow(Y),a[,i],b[,i])
            dev[,i] <- ZIM::qzinb(p[,i],lambda=mu.noUV[,i],k=1/out$psi.DE[,batch[i]],omega=avg.pi0)
          }
        }      
    }
  dev
}

