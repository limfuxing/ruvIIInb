#' Produce normalized count data after adjusting for unwanted variations
#'
#' This internal function takes ruvIII.nb or fastruvIII.nb output as input and produce various metrics of normalized data.  Option for metrics of normalized data includes percentile-adjusted count (PAC) and Pearson residuals.
#'
#' @param out  output of call to ruvIII.nb or fastruvIII.nb function.
#' @param type type of normalized data. Supported options are 'quantile' (percentile-adjusted count), 'pearson' (pearson residuals) and logcounts (log normalized count).
#' @param batch numeric vector containing batch information for each sample.Must correspond to columns of count matrix. Only needed if batch-specific dispersion parameter is fitted.
#' @param block.size the maximum number of cells for block processing when returning the corrected expression matrix. Larger block size can be quicker but requires higher RAM. Default = 5000.


#' @return A matrix containing the normalized data
get.res<-function (out, type = c("logcounts","pearson","quantile"),batch=NULL,block.size=5000) 
{
 Y <- out$counts

 if (! any(class(out$pi0) != "matrix") & !is.null(dim(out$pi0))) {
  out$pi0 <- as.matrix(out$pi0)
 }
 require(DelayedArray)
 ns    <- ncol(Y)
 block.size <- min(block.size,ncol(Y))
 setAutoRealizationBackend("HDF5Array")
 sink <- AutoRealizationSink(c(nrow(Y), ncol(Y)))
 sink_grid <- RegularArrayGrid(dim(sink), spacings=c(nrow(sink),block.size))

 nb    <- ceiling(ncol(out$counts)/block.size)
 for (bid in seq_along(sink_grid)) {
  viewport <- sink_grid[[bid]]
  start.idx<- (bid-1)*block.size+1
  end.idx  <- min(bid*block.size,ns)
  curr.batch <- batch[start.idx:end.idx] 
  Ysub  <- as.matrix(Y[,start.idx:end.idx])
  Wa <- Matrix::tcrossprod(out$a,out$W[start.idx:end.idx,])
  mu <- exp(Wa + out$gmean) 
  mu.full <- exp(Wa + out$gmean + as.matrix(out$Mb[,start.idx:end.idx,drop=FALSE]))
  Wa.mean = out$a %*% as.matrix(colMeans(out$W))

  # winsorize mean (to avoid slow process for extremely large mean)
  max.val  <- exp(matrixStats::rowMedians(log(mu)) + 4*matrixStats::rowMads(log(mu)))
  winsorize<- mu>max.val
  mu <- mu*(1-winsorize) + winsorize*max.val

  max.val  <- exp(matrixStats::rowMedians(log(mu.full)) + 4*matrixStats::rowMads(log(mu.full)))
  winsorize<- mu.full>max.val
  mu.full <- mu.full*(1-winsorize) + winsorize*max.val


  # get Pearson residual
  if (any(type == "pearson")) {
    if (is.null(dim(out$psi))) {
      if(!is.null(dim(out$pi0))) 
        dev <- (Ysub - mu * (1 - out$pi0[,curr.batch]))/sqrt(mu.full * (1 - out$pi0[,curr.batch]) *(1 + mu.full * (out$psi + out$pi0[,curr.batch])))
      if(is.null(dim(out$pi0))) 
        dev <- (Ysub - mu)/sqrt(mu.full *(1 + mu.full * out$psi))
    }

    if (!is.null(dim(out$psi))) {
      if(!is.null(dim(out$pi0)))
        dev <- (Ysub - mu * (1 - out$pi0[,curr.batch]))/sqrt(mu.full * (1 - out$pi0[,curr.batch]) *(1 + mu.full * (out$psi[,curr.batch] + out$pi0[,curr.batch])))
      if(is.null(dim(out$pi0)))
        dev <- (Ysub - mu)/sqrt(mu.full *(1 + mu.full * (out$psi[,curr.batch])))
    }
    lims <- quantile(c(dev), prob = c(0.005, 0.995), na.rm = T)
    # clip residual
    dev[dev < lims[1]] <- lims[1]
    dev[dev > lims[2]] <- lims[2]
  }
  if (any(type == "logcounts")) {
    if(!is.null(dim(out$pi0)))
      mu <- exp(Wa) * (1 - out$pi0[,curr.batch])
    if(is.null(dim(out$pi0)))
      mu <- exp(Wa) 
    dev <- log(Ysub/mu + 1)
  }
  if (any(type == "quantile")) {
    mu.noUV <- exp(as.matrix(out$Mb[,start.idx:end.idx,drop=FALSE]) + out$gmean + c(Wa.mean))
    # winsorize
    max.val  <- exp(matrixStats::rowMedians(log(mu.noUV)) + 4*matrixStats::rowMads(log(mu.noUV)))
    winsorize<- mu.noUV>max.val
    mu.noUV <- mu.noUV*(1-winsorize) + winsorize*max.val

    if (is.null(dim(out$psi))) {
      if(!is.null(dim(out$pi0)))  {
        a <- ZIM::pzinb(Ysub - 1, lambda = mu.full, k = 1/out$psi,omega = out$pi0[,curr.batch])
        b <- a + ZIM::dzinb(Ysub, lambda = mu.full, k = 1/out$psi,omega = out$pi0[,curr.batch])
      }
      if(is.null(dim(out$pi0)))  {
        a <- pnbinom(Ysub - 1, mu = mu.full, size = 1/out$psi)
        b <- a + dnbinom(Ysub, mu = mu.full, size = 1/out$psi)
      }
      p <- (a+b)/2
      p[p > 0.995] <- 0.995 ; p[p < 0.005] <- 0.005

      if(!is.null(dim(out$pi0))) {
       avg.pi0 <- rowMeans(out$pi0)
       dev <- matrix(ZIM::qzinb(p,lambda = mu.noUV,k = 1/out$psi, omega = avg.pi0),nrow(Ysub),ncol(Ysub))
      }
      if(is.null(dim(out$pi0))) {
       dev <- qnbinom(p,mu = mu.noUV,size = 1/out$psi)
      }
    }
    if (!is.null(dim(out$psi))) {
      if(!is.null(dim(out$pi0)))  {
       a <- ZIM::pzinb(Ysub - 1, lambda = mu.full, k = 1/out$psi[,curr.batch], omega = out$pi0[,curr.batch])
       b <- a + ZIM::dzinb(Ysub, lambda = mu.full, k = 1/out$psi[,curr.batch], omega = out$pi0[,curr.batch])
      }
      if(is.null(dim(out$pi0)))  {
       a <- pnbinom(Ysub - 1, mu = mu.full, size = 1/out$psi[,curr.batch])
       b <- a + dnbinom(Ysub, mu = mu.full, size = 1/out$psi[,curr.batch])
      }
      p <- (a+b)/2
      p[p > 0.995] <- 0.995 ; p[p < 0.005] <- 0.005

      if(!is.null(dim(out$pi0))) {
       avg.pi0 <- rowMeans(out$pi0)
       ref.batch <- which.min(colMeans(out$psi))
       dev <- matrix(ZIM::qzinb(p,lambda = mu.noUV,k = 1/out$psi[,ref.batch], omega = out$pi0[,ref.batch]),nrow(Ysub),ncol(Ysub))
      }
      if(is.null(dim(out$pi0))) {
       ref.batch <- which.min(colMeans(out$psi))
       dev <- qnbinom(p,mu = mu.noUV,size = 1/out$psi[,ref.batch])
      }
    }
  }
  sink <- write_block(sink, viewport, dev)
 } # end of block iteration
 tmp <- as(sink, "DelayedArray")
 tmp
}

