#' Produce normalized count data after adjusting for unwanted variations
#'
#' This internal function takes ruvIII.nb or fastruvIII.nb output as input and produce various metrics of normalized data.  Option for metrics of normalized data includes percentile-adjusted count (PAC) and Pearson residuals.
#'
#' @param out  output of call to ruvIII.nb or fastruvIII.nb function.
#' @param type type of normalized data. Supported options are 'quantile' (percentile-adjusted count), 'pearson' (pearson residuals) and logcounts (log normalized count).
#' @param batch numeric vector containing batch information for each sample.Must correspond to columns of count matrix. Only needed if batch-specific dispersion parameter is fitted.

#' @return A matrix containing the normalized data
get.res<-function (out, type = "pearson", batch = NULL) 
{
  Y <- out$counts
    if (class(Y)[1] != "matrix") {
      
    Y <- as.matrix(Y)
  }
  if (class(out$pi0) != "matrix" & !is.null(dim(out$pi0))) {
    out$pi0 <- as.matrix(out$pi0)
  }

  Wa <- Matrix::tcrossprod(out$a,out$W)
  if(!is.null(dim(out$pi0))) 
   mu <- exp(Wa + out$gmean) * (1 - out$pi0)
  if(is.null(dim(out$pi0))) 
   mu <- exp(Wa + out$gmean) 

  mu.full <- exp(Wa + out$gmean + out$Mb[,apply(out$M, 1, which)])
  Wa.mean = out$a %*% as.matrix(colMeans(out$W))
  # fixed outlier genes
  idx.gene <- abs(log(rowMeans(Y + 1)/rowMeans(mu + 1))) > log(5)
  if (sum(idx.gene) > 1) 
    mu[idx.gene, ] <- rowMeans(Y[idx.gene, ])
  if (sum(idx.gene) == 1) 
    mu[idx.gene, ] <- mean(Y[idx.gene, ])

  idx.gene <- abs(log(rowMeans(Y + 1)/rowMeans(mu.full + 1))) > log(5)
  if(!is.null(dim(out$pi0))) {
   mu.full[idx.gene, ] <- Y[idx.gene, ] * (1 - out$pi0[idx.gene,])
   mu <- mu/(1 - out$pi0)
  }
  if(is.null(dim(out$pi0))) {
   mu.full[idx.gene, ] <- Y[idx.gene, ] 
  }

  # get Pearson residual
  if (type == "pearson") {
    if (is.null(dim(out$psi))) {
      if(!is.null(dim(out$pi0))) 
        dev <- (Y - mu * (1 - out$pi0))/sqrt(mu * (1 - out$pi0) *(1 + mu * (out$psi + out$pi0)))
      if(is.null(dim(out$pi0))) 
        dev <- (Y - mu)/sqrt(mu *(1 + mu * out$psi))
    }

    if (!is.null(dim(out$psi))) {
      if(!is.null(dim(out$pi0)))
        dev <- (Y - mu * (1 - out$pi0))/sqrt(mu * (1 - out$pi0) *(1 + mu * (out$psi[, batch] + out$pi0)))
      if(is.null(dim(out$pi0)))
        dev <- (Y - mu)/sqrt(mu *(1 + mu * (out$psi[, batch])))
    }
    lims <- quantile(c(dev), prob = c(5e-04, 0.9995), na.rm = T)
    dev[dev < lims[1]] <- lims[1]
    dev[dev > lims[2]] <- lims[2]
  }
  if (type == "logcounts") {
    if(!is.null(dim(out$pi0)))
      mu <- exp(Wa) * (1 - out$pi0)
    if(is.null(dim(out$pi0)))
      mu <- exp(Wa) 
    dev <- log(Y/mu + 1)
  }
  if (type == "quantile") {
    mu.noUV <- exp(out$Mb[, apply(out$M, 1, which)] + out$gmean + c(Wa.mean))
    idx.gene <- abs(log(rowMeans(Y + 1)/rowMeans(mu.noUV + 1))) > log(5)
    if (sum(idx.gene) > 1) 
      mu.noUV[idx.gene, ] <- rowMeans(Y[idx.gene, ])
    if (sum(idx.gene) == 1) 
      mu.noUV[idx.gene, ] <- mean(Y[idx.gene, ])
    if (is.null(dim(out$psi))) {
      if(!is.null(dim(out$pi0)))  {
        a <- ZIM::pzinb(Y - 1, lambda = mu.full, k = 1/out$psi,omega = out$pi0)
        b <- ZIM::pzinb(Y, lambda = mu.full, k = 1/out$psi,omega = out$pi0)
      }
      if(is.null(dim(out$pi0)))  {
        a <- pnbinom(Y - 1, mu = mu.full, size = 1/out$psi)
        b <- pnbinom(Y, mu = mu.full, size = 1/out$psi)
      }
      a[a > 0.999] <- 0.9985
      b[b > 0.999] <- 0.999
      p <- (a+b)/2
      if(!is.null(dim(out$pi0))) {
       avg.pi0 <- rowMeans(out$pi0)
       dev <- matrix(ZIM::qzinb(p,lambda = mu.noUV,k = 1/out$psi, omega = avg.pi0),nrow(out$counts),ncol(out$counts))
      }
      if(is.null(dim(out$pi0))) {
       dev <- qnbinom(p,mu = mu.noUV,size = 1/out$psi)
      }
    }
    if (!is.null(dim(out$psi))) {
      if(!is.null(dim(out$pi0)))  {
       a <- ZIM::pzinb(Y - 1, lambda = mu.full, k = 1/out$psi[,batch], omega = out$pi0)
       b <- ZIM::pzinb(Y, lambda = mu.full, k = 1/out$psi[,batch], omega = out$pi0)
      }
      if(is.null(dim(out$pi0)))  {
       a <- pnbinom(Y - 1, mu = mu.full, size = 1/out$psi[,batch])
       b <- pnbinom(Y, mu = mu.full, size = 1/out$psi[,batch])
      }
      a[a > 0.999] <- 0.9985
      b[b > 0.999] <- 0.999
      p <- (a+b)/2
      if(!is.null(dim(out$pi0))) {
       avg.pi0 <- rowMeans(out$pi0)
       dev <- matrix(ZIM::qzinb(p,lambda = mu.noUV,k = 1/out$psi[,batch], omega = avg.pi0),nrow(out$counts),ncol(out$counts))
      }
      if(is.null(dim(out$pi0))) {
       dev <- qnbinom(p,mu = mu.noUV,size = 1/out$psi[,batch])
      }
    }
  }
  dev
}
