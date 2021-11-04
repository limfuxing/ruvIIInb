#' Removing unwanted variation from sequencing count data 
#'
#' This function performs ruvIIInb method to remove unwanted variation from sequencing count data. The method has been applied to single-cell RNA-seq data but in principle it is applicable to count data from any sequencing platforms. It takes raw count matrix with features (genes/transcript in scRNA-seq) as rows and samples (cells in scRNA-seq) as columns. 
#' Users need to specify the set of negative control genes (i.e genes where the between cells variation is assumed to be solely due to the unwanted variation) 
#' and the (pseudo)replicate matrix that define sets of cells that can be considered technical replicates.
#'
#' @param Y raw count matrix with genes as rows and cells as columns
#' @param M replicate matrix with number of rows equal to the number of cells and number of columns equal to the number of distinct sub-population of cells: M(i,j) = 1 if cell i is part of replicate j and 0 otherwise. If cell-type is known, this can be used to form the replicate matrix, with each cell-type corresponds to one sub-population of cells. If cell-type is unknown, the scMerge::scReplicate command can be used to identify the pseudo-replicates.
#' @param ctl  either logical vector of length n with TRUE = control genes OR a character vector of names of control genes
#' @param k dimension of unwanted variation against which the data is to be normalized against. Default is 2.
#' @param robust logical value. Whether to use Huber's robust weight in the iterative reweighted least squares (IRLS) algorithm. Default is FALSE.
#' @param ortho.W logical value. Whether the unwanted factors (W) are to be orthogonalized. The default is FALSE.
#' @param lambda.a smoothing parameter for regularizing regression coefficients associated with W. The default is 0.01
#' @param lambda.b smoothing parameter for regularizing the gene-level intercepts. The default is 16.
#' @param batch numeric vector containing batch information for each cell. Must correspond to columns of count matrix. Only needed if batch-specific dispersion parameter is
#'        fitted AND/OR use.pseudosamples=TRUE.
#' @param step.fac multiplicative factor to decrease IRLS step by when log-likelihood diverges. Default is 0.5
#' @param inner.maxit the maximum number of IRLS iteration for estimating mean parameters for a given dispersion parameter. Default is 50
#' @param outer.maxit the maximum number of IRLS iteration. Default is 25 
#' @param ncores The number of cores used for parallel processing. The default number of workers (cores) is 2.
#' @param use.pseudosample whether to use pseudocell to define replicates (default is FALSE). Note that the replicates defined by the pseudocells will be used in addition to any replicates defined by the M matrix above. If no replicate is defined by the M matrix (ie M matrix is an identity matrix) then the unwanted variation will be estimated using the pseudoreplicates defined by the pseudocells only. We recommend that use.pseudosample=TRUE be used only for data with UMI.
#' @param nc.pool the number of cells per pool (used for defining pseudocells). The default is 20. Only relevant when use.pseudosample=TRUE.
#' @param strata By default the creation of pseudo-samples are stratified by batch and library size ONLY. This default strategy should work well if the biological factor of 
#'        interest is not correlated with batch. If they are correlated, this strategy will risk removing too much biology. Specifying the biological factor as 'strata' variable helps reducing the amount of biology removed.
#' @param batch.disp whether to fit batch-speficic dispersion parameter for each gene. The default is FALSE.
#' @param zeroinf logical vector indicating whether to fit zero-inflated negative binomial (ZINB) instead of NB model (the default) for some/all cells. The length of the vector #'        must equal the number of cells. We recommend that non-UMI data be fitted using ZINB and UMI data be fitted using NB model.
  
#' @return A list containing the raw data (as sparse matrix), the unwanted factors and regression coefficients associated with the unwanted factors.
#' @export
ruvIII.nb <- function(Y,M,ctl,k=2,robust=FALSE,ortho.W=FALSE,lambda.a=0.01,lambda.b=16,batch=NULL,step.fac=0.5,inner.maxit=50,outer.maxit=25,
        ncores=2,use.pseudosample=FALSE,nc.pool=20,strata=NULL,batch.disp=FALSE,zeroinf=NULL) {

# register parallel backend
#register(BPPARAM)

# setup cluster for doParallel
doParallel::registerDoParallel(ncores)
BPPARAM=BiocParallel::DoparParam()

parallel <- as.logical(ncores>1)

if(is.null(zeroinf))
  zeroinf <- rep(FALSE,ncol(Y))

if(!is.logical(ctl)) 
  ctl <- rownames(Y) %in% ctl
idx <- which(ctl)
nctl<- sum(ctl)

# force M matrix into logical
mode(M) <- 'logical'

# Winsorize gene-by-gene and calculate size-factor
Y <- ceiling(t(apply(Y,1,DescTools::Winsorize,probs=c(0,0.995),na.rm=TRUE)))
sf <- colSums(Y)
sf <- sf/mean(sf) 
# remove genes with all zero counts (if any)
zero.genes <- rowSums(Y>0)==0
if(any(zero.genes)) {
   Y <- Y[!zero.genes,]
   ctl <- ctl[!zero.genes]
   print(paste0(sum(zero.genes), 'genes with all zero counts are removed following a Winsorization step.'))
}


# if no batch info is supplied but batch.disp=TRUE, forced batch.disp=FALSE
if(is.null(batch) & batch.disp) {
  print(paste0('Batch variable not supplied...cannot fit batch-specific dispersion parameter.'))
  batch.disp <- FALSE
}

# adding pseudosamples
ncells <- rep(1,length(sf))
sample.type <- rep('sc',length(sf))
nsc <- ncol(Y)
if(use.pseudosample) {
 if(is.null(strata)) {
   batch.ext <- batch
   for(btc in unique(batch)) {
    nq  <- max(round(sum(batch==btc)/nc.pool),1)
    qsf <-  as.numeric(gtools::quantcut(sf[batch==btc],q=nq))
    for(q in 1:max(qsf)) {
       pool <- Y[,1:nsc][,batch==btc][,qsf==q]
       dns.pool <- rbinom(nrow(pool),size=rowSums(pool),prob=1/ncol(pool))
       Y <- cbind(Y,dns.pool)
       #ncells <- c(ncells, sum(sf[batch==btc][qsf==q])^2/sum(sf[batch==btc][qsf==q]^2))
       sample.type <- c(sample.type,'sc')
       batch.ext  <- c(batch.ext,btc+max(batch))
       M      <- rbind(M,FALSE)
    }
   }
 ps.mcol <- rep(FALSE,nrow(M))
 ps.mcol[seq(nsc+1,ncol(Y),1)] <- TRUE
 M <- cbind(M,ps.mcol)
 batch <- batch.ext
 }

 if(!is.null(strata)) {
  batch.ext <- batch
  for(st in unique(na.omit(strata))) {
    temp.M <- M
    for(btc in unique(batch)) {
     nq  <- max(round(sum(batch==btc & strata==st,na.rm=TRUE)/nc.pool),1)
     if(sum(batch==btc & strata==st, na.rm=TRUE)>1) {
      qsf <-  na.omit(as.numeric(gtools::quantcut(sf[batch==btc & strata==st],q=nq)))
      for(q in 1:max(qsf)) {
       pool <- Y[,1:nsc][,batch==btc & strata==st & !is.na(strata)][,qsf==q]
       dns.pool <- rbinom(nrow(pool),size=rowSums(pool),prob=1/ncol(pool))
       Y <- cbind(Y,dns.pool)
       #ncells <- c(ncells, sum(sf[batch==btc][qsf==q])^2/sum(sf[batch==btc][qsf==q]^2))
       sample.type <- c(sample.type,'sc')
       batch.ext  <- c(batch.ext,btc+max(batch))
       temp.M      <- rbind(temp.M,FALSE)
      }
     } #if
    } #for batch
   ps.mcol <- rep(FALSE,nrow(temp.M))
   if(sum(strata==st,na.rm=TRUE)>1) {
    ps.mcol[seq(nrow(M)+1,nrow(temp.M),1)] <- TRUE
    temp.M <- cbind(temp.M,ps.mcol)
    M      <- temp.M
   }
  } # for st
 if(!batch.disp) 
   batch.ext <- ifelse(batch.ext>max(batch),2,1)
   
 batch <- batch.ext

 } # if
}

# if using pseudosamples, only NB models can be fitted
if(use.pseudosample)
  zeroinf <- rep(FALSE,ncol(Y))


# get Huber's k
k.huber<- ifelse(robust,1.345,100)

ncells <- rep(1,ncol(Y))
sf <- rep(1,ncol(Y))
#sf <- colSums(Y)
sf <- sf/mean(sf) 
Y.norm <- sweep(as.matrix(Y),2,sf,'/')
off.g  <- log(rowMeans(Y.norm))
Z.res  <- log(Y.norm+1)

 
# initiate W by performing SVD on the control log(Yc+1)
Y.norm <- sweep(Y.norm,1,exp(off.g),'/')
svd.out <- irlba::irlba(log(Y.norm[ctl,]+1),nv=k)
W  <- svd.out$v

#format and convert M matrix into list of indices to save space
rep.ind <- apply(M,2,which)


# initiate beta and alpha: by performing poisson GLM gene-by-gene
nW   <- ncol(W)
nM   <- ncol(M)
ngene<- nrow(Y)
ns   <- ncol(Y)
coef <- matrix(0,ngene,nW+1)
wt   <- matrix(1,ngene,ns)

Ylist <- lapply(1:ngene,FUN=subsetMat,mat=Z.res,MARGIN=1)
wtlist<- lapply(1:ngene,FUN=subsetMat,mat=Y+1,MARGIN=1)
coef  <- drop(BiocParallel::bpmapply(FUN=get.coef.wt,y=Ylist,wt=wtlist,MoreArgs=list(W=W),BPPARAM=BPPARAM))

# alpha(n x k) matrix
alpha<- as.matrix(t(coef)[,-1])
alpha.mean <- colMeans(alpha)
alpha.dev <- sweep(alpha,2,alpha.mean,'-')

# Mb (n x m_1) matrix
Mb   <- matrix(0,nrow=ngene,ncol=nM,byrow=FALSE)

#handling NA/Inf
Mb[which(is.na(Mb) | is.infinite(Mb))] <- 0
alpha[which(is.na(alpha) | is.infinite(alpha))] <- 0

#gmean
gmean<- c(coef[1,])

# initiate zinb weight
wt.zinb<- matrix(1,ngene,ns)
pi0    <- rep(0,ngene)

# estimate initial psi 
Wa   <- alpha %*% t(W) 
offs <- Wa + Mb[,apply(M,1,which)] 
nsc  <- sum(sample.type=='sc')
if(parallel) {
 if(is.null(batch)) 
  psi <- estimateDisp.par(Y[,sample.type=='sc'],as.matrix(rep(1,nsc)),offset=offs[,sample.type=='sc'],tagwise=TRUE,robust=TRUE,BPPARAM=BPPARAM,
				weights=wt.zinb[,sample.type=='sc'])$tagwise.dispersion
 if(!is.null(batch)) {
  psi <- estimateDisp.par(Y[,batch==1 & sample.type=='sc'],as.matrix(rep(1,sum(batch==1 & sample.type=='sc'))),
				offset=offs[,batch==1 & sample.type=='sc'],tagwise=TRUE,robust=TRUE,BPPARAM=BPPARAM,
				weights=wt.zinb[,batch==1 & sample.type=='sc'])$tagwise.dispersion
  for(B in 2:max(batch)) 
      psi <- cbind(psi,estimateDisp.par(Y[,batch==B & sample.type=='sc'],as.matrix(rep(1,sum(batch==B & sample.type=='sc'))),
				offset=offs[,batch==B & sample.type=='sc'],tagwise=TRUE,robust=TRUE,BPPARAM=BPPARAM,
				weights=wt.zinb[,batch==B & sample.type=='sc'])$tagwise.dispersion)
 }
} 
if(!parallel) {
 if(is.null(batch)) 
  psi <- edgeR::estimateDisp(Y,as.matrix(rep(1,nsc)),offset=offs[,sample.type=='sc'],tagwise=TRUE,robust=TRUE,weights=wt.zinb[,sample.type=='sc'])$tagwise.dispersion
 if(!is.null(batch)) {
  psi <- edgeR::estimateDisp(Y[,batch==1 & sample.type=='sc'],as.matrix(rep(1,sum(batch==1 & sample.type=='sc'))),
				offset=offs[,batch==1 & sample.type=='sc'],tagwise=TRUE,robust=TRUE,
				weights=wt.zinb[,batch==1 & sample.type=='sc'])$tagwise.dispersion
  for(B in 2:max(batch)) 
   psi <- cbind(psi,edgeR::estimateDisp(Y[,batch==B & sample.type=='sc'],as.matrix(rep(1,sum(batch==B & sample.type=='sc'))),
				offset=offs[,batch==B & sample.type=='sc'],tagwise=TRUE,robust=TRUE,
				weights=wt.zinb[,batch==B & sample.type=='sc'])$tagwise.dispersion)
 }
}

# estimate pi0 (if needed)
if(any(zeroinf)) {
 mu  <- exp(offs+gmean)
 tmp <- foreach(i=1:ngene, .combine=c, .packages="ZIM") %dopar% {
   y     <- Y[i,zeroinf]
   muvec <- mu[i,zeroinf]
   if(is.null(batch))
    psivec<- rep(psi[i],sum(zeroinf,na.rm=T))
   if(!is.null(batch))
    psivec<- psi[i,batch][zeroinf]

   out <- tryCatch({ optim(p=-2,fn=logl.zinb,y=y,mu=muvec,size=1/psivec)$p}, error = function(e) {NA})
   return(1/(1+exp(-out)))
 }
 pi0[!is.na(tmp)] <- tmp[!is.na(tmp)]
}

# update zinb weight
if(any(zeroinf)) {
 for(i in which(zeroinf)) {
  zero.Y <- which(Y[,i]==0)
  if(is.null(batch))
    psivec<- psi[zero.Y]
  if(!is.null(batch))
    psivec<- psi[zero.Y,batch[i]]
  
  wt.zinb[zero.Y,i] <- 1 - pi0[zero.Y]/(pi0[zero.Y] + (1-pi0[zero.Y])*dnbinom(0,mu=exp(offs[zero.Y,i]+gmean[zero.Y]),size=1/psivec))
 }
}

# initiate
conv <- FALSE
iter.outer <- 0
logl.outer <- NULL
step <- rep(1,ngene)
while(!conv) {
 iter.outer <- iter.outer + 1
 logl.beta <- NULL
 
 ### calculate initial loglik
 lmu.hat    <- gmean + Mb[,apply(M,1,which)] + alpha %*% t(W) 

 if(!is.null(batch))
  temp <- ZIM::dzinb(Y,lambda=exp(lmu.hat),k=sweep(1/psi[,batch],2,ncells,'*'),omega=outer(pi0,zeroinf),log=TRUE)
 if(is.null(batch))
  temp <- ZIM::dzinb(Y,lambda=exp(lmu.hat),k=outer(1/psi,ncells,'*'),omega=outer(pi0,zeroinf),log=TRUE)
 
 temp[is.infinite(temp) | is.na(temp)] <- log(.Machine$double.xmin)
 loglik <- rowSums(temp) #- lambda.a*rowSums(scale(alpha,scale=FALSE)^2) - lambda.b*rowSums(Mb^2)
 #############################
 step <- rep(1,ngene)
 conv.beta <- FALSE
 # saving temporary best estimate
 best.W <- W ; best.psi <- psi ; best.a <- alpha; best.Mb <- Mb; best.adev <- alpha.dev; best.gmean <- gmean; best.pi0 <- pi0
 
 iter <- 1
 halving <- 0
 while(!conv.beta) {
 print(paste0('Inner iter:', iter))
 # calc working vector Z
 lmu.hat    <- gmean + Mb[,apply(M,1,which)] + alpha %*% t(W) 

 # weight based on NB
 if(!is.null(batch))
  sig.inv<- 1/(exp(-lmu.hat) + sweep(psi[,batch],2,ncells,'/'))
 if(is.null(batch))
  sig.inv<- 1/(exp(-lmu.hat) +outer(psi,ncells,'/'))
 
 # calculate signed deviance
 if(!is.null(batch)) {
   signdev <- sign(Y-exp(lmu.hat))* sqrt(2*(dnbinom(Y,mu=Y,size=sweep(1/psi[,batch],2,ncells,'*'),log=TRUE) - 
						dnbinom(Y,mu=exp(lmu.hat),size=sweep(1/psi[,batch],2,ncells,'*'),log=TRUE) ))
 }
 if(is.null(batch)) {
   signdev <- sign(Y-exp(lmu.hat))* sqrt(2*(dnbinom(Y,mu=Y,size=outer(1/psi,ncells,'*'),log=TRUE) - dnbinom(Y,mu=exp(lmu.hat),size=outer(1/psi,ncells,'*'),log=TRUE) ))
 }  
 
 
 Z      <- lmu.hat + step*((Y+0.01)/(exp(lmu.hat)+0.01) - 1)

 # project Z and W into M
 
 # step1: obtain projecion of Z and W onto the range of M and get the residual
 # RM_Z is (n x m) matrix, m=nsample,n=ngene
 # RM_W is (k x m) matrix
 proj.Z <- projection(rep.ind,Zmat=Z)[,apply(M,1,which)]
 RM_Z <- Z - proj.Z
 
 if(nW>1) 
   proj.tW <- projection(rep.ind,Zmat=t(W))
 if(nW==1) 
  proj.tW <- t(as.matrix(projection.K1(rep.ind,Zmat=t(W))))
 
 if(nW>1) 
  RM_W <- W - t(proj.tW[,apply(M,1,which)])

 if(nW==1) 
  RM_W <- W - as.matrix(proj.tW[,apply(M,1,which)])

 # store current alpha est
 alpha.old <- alpha

 # get alpha mean
 RM_Z.list <- lapply(1:ngene,FUN=subsetMat,mat=RM_Z,MARGIN=1)
 adev.list <- lapply(1:ngene,FUN=subsetMat,mat=alpha.dev,MARGIN=1)
 wtlist    <- lapply(1:ngene,FUN=subsetMat,mat=sig.inv*wt.zinb,MARGIN=1)
 signdev.list <- lapply(1:ngene,FUN=subsetMat,mat=signdev,MARGIN=1)
 results <- BiocParallel::bpmapply(FUN=get.coef2,RM_Z=RM_Z.list,signdev=signdev.list,alpha.dev=adev.list,wt=wtlist,MoreArgs=list(RM_W=RM_W,k.huber=k.huber),BPPARAM=BPPARAM)
 if(k>1)
  alpha.mean <- solve(as.matrix(apply(results[,-(nW+1),],c(1,2),sum)),as.matrix(apply(results[,nW+1,],1,sum)))
 if(k==1)
  alpha.mean <- c(sum(results[,(nW+1),])/sum(results[,-(nW+1),]))
 # get alpha dev 
 alpha.dev <- drop(BiocParallel::bpmapply(FUN=get.adev,RM_Z=RM_Z.list,signdev=signdev.list,wt=wtlist,
		MoreArgs=list(RM_W=RM_W,lambda.a=lambda.a,alpha.mean=alpha.mean,k.huber=k.huber),BPPARAM=BPPARAM))

 if(k>1)
  alpha.dev <- t(alpha.dev)
 if(k==1)
  alpha.dev <- as.matrix(alpha.dev)

 # get alpha
 alpha <- sweep(alpha.dev,2,alpha.mean,'+')
 # if new alpha est has missing/inf values, revert to previous estimate
 print(paste0('Problematic alpha entries=',sum(is.na(alpha) | is.infinite(alpha))))
 if(any(is.na(alpha) | is.infinite(alpha) )) alpha <- alpha.old

 # reduce outliers
 a.med <- apply(alpha,2,median)
 a.mad <- apply(alpha,2,mad)
 for(i in 1:nW) {
   alpha[which(alpha[,i]> (a.med[i]+4*a.mad[i])),i] <- a.med[i]+4*a.mad[i]
   alpha[which(alpha[,i]< (a.med[i]-4*a.mad[i])),i] <- a.med[i]-4*a.mad[i]
 }


 # step2: update W using WLS of Z_c (Z for control genes) on alpha.c (alpha for control genes)
 W.old <- W 
 Z.list <- lapply(1:ns,FUN=subsetMat,mat=Z-gmean,MARGIN=2)
 wtlist    <- lapply(1:ns,FUN=subsetMat,mat=sig.inv*wt.zinb,MARGIN=2)
 signdev.list <- lapply(1:ns,FUN=subsetMat,mat=signdev,MARGIN=2)
 zeroinf.list <- lapply(1:ns,FUN=subsetMat,mat=as.matrix(zeroinf),MARGIN=1)
 W     <- BiocParallel::bpmapply(FUN=getW,Z=Z.list,signdev=signdev.list,wt=wtlist,zeroinf=zeroinf.list,
			MoreArgs=list(alpha=alpha,ctl=ctl,k.huber=k.huber,k=k),BPPARAM=BPPARAM)
 if(k>1) 
  W <- t(W)
 if(k==1)
  W <- as.matrix(W)

 colNorm.W <- colSums(W^2)
 # normalize W to have unit length for each col
 W  <- sweep(W,2,sqrt(colNorm.W),'/')
 
 # orthogonalize W (optional)
 if(nW>1 & ortho.W)  
   W <- matlib::GramSchmidt(W)

 # if new W has missing/inf values, revert to previous estimate
 print(paste0('Problematic W entries=',sum(is.na(W) | is.infinite(W))))
 if(any(is.na(W) | is.infinite(W) )) W <- W.old
 
 #step3a: update gmean
 Z.res  <- Z - alpha %*% t(W) - Mb[,apply(M,1,which)]

 temp <- tryCatch({
 for(i in 1:ngene) {
    sigma.est <- mad(signdev[i,])/0.6745
    wt[i,]    <- MASS::psi.huber(signdev[i,],k=k.huber*sigma.est)
    wtvec     <- wt[i,] * sig.inv[i,] * wt.zinb[i,]
    wtvec     <- wtvec/sum(wtvec)
    gmean[i]  <- sum(Z.res[i,]*wtvec)
 }
 }, error = function(e) {NA})

 #step 3b: update Mb
 Z.res  <- Z - alpha %*% t(W) - gmean
 Mb.old <- Mb
 Mb <- tryCatch({ projection(rep.ind,Zmat=Z.res,Wmat=wt*sig.inv*wt.zinb,lambda=lambda.b)}, error = function(e) {NA})

 # if new Mb has missing/inf values, revert to previous estimate
 print(paste0('Problematic Mb entries=',sum(is.na(Mb) | is.infinite(Mb))))
 if(any(is.na(Mb))) Mb <- Mb.old
 Mb[which(is.na(Mb) | is.infinite(Mb))] <- 0

 # reduce outliers
 Mb.med <- apply(Mb,1,median)
 Mb.mad <- apply(Mb,1,mad)
 for(i in 1:ngene) {
   Mb[i,which(Mb[i,]> (Mb.med[i]+4*Mb.mad[i]))] <- Mb.med[i]+4*Mb.mad[i]
   Mb[i,which(Mb[i,]< (Mb.med[i]-4*Mb.mad[i]))] <- Mb.med[i]-4*Mb.mad[i]
 }

 # step4a: update pi0 (if needed)
 pi0.old <- pi0
 if(any(zeroinf)) {
  mu  <- exp(alpha %*% t(W) + Mb[,apply(M,1,which)] + gmean)
  tmp <- foreach(i=1:ngene, .combine=c, .packages="ZIM") %dopar% {
   y     <- Y[i,zeroinf]
   muvec <- mu[i,zeroinf]
   if(is.null(batch))
    psivec<- rep(psi[i],sum(zeroinf,na.rm=T))
   if(!is.null(batch))
    psivec<- psi[i,batch][zeroinf]

   out <- tryCatch({ optim(p=-2,fn=logl.zinb,y=y,mu=muvec,size=1/psivec)$p}, error = function(e) {NA})
   return(1/(1+exp(-out)))
  }
  pi0[!is.na(tmp)] <- tmp[!is.na(tmp)]
 }
 # if new pi0 has missing/inf values, revert to previous estimate
 print(paste0('Problematic pi0 entries',sum(is.na(pi0) | is.infinite(pi0))))
 if(any(is.na(pi0) | is.infinite(pi0) )) pi0 <- pi0.old

 # update zinb weight
 if(any(zeroinf)) {
  for(i in which(zeroinf)) {
   zero.Y <- which(Y[,i]==0)
   if(is.null(batch))
    psivec<- psi[zero.Y]
   if(!is.null(batch))
    psivec<- psi[zero.Y,batch[i]]
  
   wt.zinb[zero.Y,i] <- 1 - pi0[zero.Y]/(pi0[zero.Y] + (1-pi0[zero.Y])*dnbinom(0,mu=mu[zero.Y,i],size=1/psivec))
  }
 }

 # calculate current logl
 lmu.hat    <- gmean + Mb[,apply(M,1,which)] + alpha %*% t(W) 
 if(!is.null(batch))
  temp <- ZIM::dzinb(Y,lambda=exp(lmu.hat),k=sweep(1/psi[,batch],2,ncells,'*'),omega=outer(pi0,zeroinf),log=TRUE)
 if(is.null(batch))
  temp <- ZIM::dzinb(Y,lambda=exp(lmu.hat),k=outer(1/psi,ncells,'*'),omega=outer(pi0,zeroinf),log=TRUE)
 
 temp[is.infinite(temp) | is.na(temp)] <- log(.Machine$double.xmin)
 loglik.tmp <- rowSums(temp) #- lambda.a*rowSums(scale(alpha,scale=FALSE)^2) - lambda.b*rowSums(Mb^2)
 # check degenerate case
 if(iter>=1) {
  degener<- sum(loglik)>sum(loglik.tmp)
  degener[is.na(degener)] <- TRUE
  if(degener) {
    check.gene <- loglik>loglik.tmp
    step[check.gene] <- step[check.gene]*step.fac
    W    <- best.W; alpha <- best.a ; Mb <- best.Mb ; psi <- best.psi ; alpha.dev <- best.adev ; gmean <- best.gmean ; pi0 <- best.pi0
    # revert zinb weight
    if(any(zeroinf)) {
     mu  <- exp(alpha %*% t(W) + Mb[,apply(M,1,which)] + gmean)
     for(i in which(zeroinf)) {
      zero.Y <- which(Y[,i]==0)
      if(is.null(batch))
       psivec<- psi[zero.Y]
      if(!is.null(batch))
       psivec<- psi[zero.Y,batch[i]]
      wt.zinb[zero.Y,i] <- 1 - pi0[zero.Y]/(pi0[zero.Y] + (1-pi0[zero.Y])*dnbinom(0,mu=mu[zero.Y,i],size=1/psivec))
     }
    }
    halving <- halving + 1
    if(halving>=3) {
     loglik.tmp <- loglik
     degener <- !degener
    }
  }
  if(!degener) {
   loglik     <- loglik.tmp
   logl.beta  <- c(logl.beta,sum(loglik))
   print(paste0('Outer Iter ',iter.outer, ', Inner iter ', iter, ' logl-likelihood:', logl.beta[length(logl.beta)]))
   # saving temporary best estimate
   best.W <- W ; best.psi <- psi ; best.a <- alpha; best.Mb <- Mb ; best.adev <- alpha.dev; best.gmean <- gmean ; best.pi0 <- pi0
   iter <- iter + 1
   halving <- 0
  }
 }

conv.logl <- FALSE
if(iter>2) {
 conv.logl <- ifelse( (logl.beta[length(logl.beta)] - logl.beta[length(logl.beta)-1])/abs(logl.beta[length(logl.beta)-1]) < 1e-04,TRUE,FALSE)
}
conv.beta <- iter>=inner.maxit  | conv.logl
} # end of IRLS inner loop

#update psi
logl.outer.tmp <- sum(loglik,na.rm=T)
updt.psi <- TRUE
if(iter.outer>1)
   updt.psi <- abs(logl.outer[length(logl.outer)]-logl.outer.tmp)/abs(logl.outer[length(logl.outer)]) > 1e-08  

Wa = best.a %*% t(best.W)
offs <- Wa + best.Mb[,apply(M,1,which)] 
if(updt.psi) {
if(parallel) {
 if(is.null(batch)) 
  psi.new <- tryCatch({ estimateDisp.par(Y[,sample.type=='sc'],as.matrix(rep(1,nsc)),offset=offs[,sample.type=='sc'],
				tagwise=TRUE,robust=TRUE,BPPARAM=BPPARAM,weights=wt.zinb[,sample.type=='sc'])$tagwise.dispersion}, error = function(e) {NA})
 if(!is.null(batch)) {
  psi.new <- psi
  for(B in 1:max(batch)) 
   psi.new[,B]<-tryCatch({estimateDisp.par(Y[,batch==B & sample.type=='sc'],as.matrix(rep(1,sum(batch==B & sample.type=='sc'))),
				offset=offs[,batch==B & sample.type=='sc'],tagwise=TRUE,robust=TRUE,BPPARAM=BPPARAM,
				weights=wt.zinb[,batch==B & sample.type=='sc'])$tagwise.dispersion},error= function(e) {NA})
 }
}

if(!parallel) {
 if(is.null(batch)) 
  psi.new <- tryCatch({ edgeR::estimateDisp(Y[,sample.type=='sc'],as.matrix(rep(1,nsc)),offset=offs[,sample.type=='sc'],tagwise=TRUE,robust=TRUE,
				weights=wt.zinb[,sample.type=='sc'])$tagwise.dispersion},error = function(e) {NA})
 if(!is.null(batch)) {
  psi.new <- psi
  for(B in 1:max(batch)) 
   psi.new[,B] <- tryCatch({ edgeR::estimateDisp(Y[,batch==B & sample.type=='sc'],as.matrix(rep(1,sum(batch==B & sample.type=='sc'))),
				offset=offs[,batch==B & sample.type=='sc'],tagwise=TRUE,robust=TRUE,
				weights=wt.zinb[,batch==B & sample.type=='sc'])$tagwise.dispersion},error = function(e) {NA})
 }
}
# update psi
psi[!is.na(psi.new) & !is.infinite(psi.new)] <- psi.new[!is.na(psi.new) & !is.infinite(psi.new)]
}

# record outer logl
W <- best.W
alpha <- best.a
Mb    <- best.Mb
alpha.dev <- best.adev
gmean <- best.gmean
pi0   <- best.pi0
# new changes
best.psi <- psi

# recalculate working vector and weight
lmu.hat    <- gmean + Mb[,apply(M,1,which)] + alpha %*% t(W) 
 if(!is.null(batch))
  temp <- ZIM::dzinb(Y,lambda=exp(lmu.hat),k=sweep(1/psi[,batch],2,ncells,'*'),omega=outer(pi0,zeroinf),log=TRUE)
 if(is.null(batch))
  temp <- ZIM::dzinb(Y,lambda=exp(lmu.hat),k=outer(1/psi,ncells,'*'),omega=outer(pi0,zeroinf),log=TRUE)

temp[is.infinite(temp) | is.na(temp)] <- log(.Machine$double.xmin)
logl.outer <- c(logl.outer, sum(temp)) #-lambda.a*sum(scale(alpha,scale=FALSE)^2)-lambda.b*sum(Mb^2))
print(paste0('Outer Iter ',iter.outer, ' logl-likelihood:', logl.outer[length(logl.outer)]))

if(iter.outer>1) {
  conv <- ifelse( (logl.outer[iter.outer] - logl.outer[iter.outer-1])/abs(logl.outer[iter.outer-1]) < 1e-04 | iter.outer>=outer.maxit,TRUE,FALSE)
}

if(conv) {
 W <- best.W
 alpha <- best.a
 gmean <- best.gmean
 Mb  <- best.Mb
 psi   <- best.psi
 pi0 <- best.pi0
 alpha.c <- as.matrix(alpha)[ctl,]
}
} # end of outer IRLS loop

# calculate psi for DE (where Mb is left out of the offset)
psi.DE <- psi
# save raw count as sparse matrix
Y  <- Matrix::Matrix(Y,sparse=TRUE)
# prepare and save zero-inflation parameter as parse matrix
pi0<- outer(pi0,as.numeric(zeroinf))
pi0<- Matrix::Matrix(pi0,sparse=TRUE)
#output
return( list("counts"=Y,"W"=W, "M"=M, "ctl"=ctl, "logl"=logl.outer, "a"=alpha,"Mb"=Mb, "gmean"=gmean,"pi0"=pi0,
		"psi"=psi,'psi.DE'=psi.DE,'L.a'=lambda.a,'L.b'=lambda.b,batch=batch) )
}

#### Helpers functions ######################################################################################
rob.wt <- function(Y,mu,psi) {
 sign.dev <- sign(Y-mu)* sqrt(2*(dnbinom(Y,mu=Y,size=1/psi,log=TRUE) - dnbinom(Y,mu=mu,size=1/psi,log=TRUE) ))
 mad.est  <- mad(sign.dev)
 sign.dev <- sign.dev/mad.est
 matrix(MASS::psi.huber(sign.dev),byrow=T,nrow=nrow(sign.dev),ncol=ncol(sign.dev))
}

projection <- function(ind,Zmat,Wmat=NULL,lambda=0) {
 if(is.null(Wmat)) 
   out <- sapply(ind,rowAvg,Z=Zmat)
 if(!is.null(Wmat)) 
   out <- sapply(ind,rowWtAvg,Z=Zmat,Wt=Wmat,lambda=lambda)
   
out
}

projection.K1 <- function(ind,Zmat,Wmat=NULL) {
 if(is.null(Wmat)) 
   out <- sapply(ind,rowAvg.K1,Z=Zmat)
 if(!is.null(Wmat)) 
   out <- sapply(ind,rowWtAvg.K1,Z=Zmat,Wt=Wmat)
   
out
}
   
rowAvg <- function(ind,Z) {
  if(length(ind)>1) {
   out <- rowMeans(as.matrix(Z[,ind]), na.rm=T)
  }
  if(length(ind)==1) {
    out <- Z[,ind]
  }
out
}

rowWtAvg <- function(ind,Z,Wt,lambda=0) {
 if(length(ind)>1) {
  out <- rowSums(as.matrix(Z[,ind]*Wt[,ind]), na.rm=T)/(rowSums(as.matrix(Wt[,ind]), na.rm=T) + lambda)  
 }
 if(length(ind)==1) {
  out <- c(Z[,ind])*c(Wt[,ind])/c(Wt[,ind]+lambda)
 }
out
}
 
rowAvg.K1 <- function(ind,Z) {
  if(length(ind)>1) {
   out <- colMeans(as.matrix(Z[,ind]), na.rm=T)
  }
  if(length(ind)==1) {
    out <- Z[,ind]
  }
out
}

rowWtAvg.K1 <- function(ind,Z,Wt) {
 if(length(ind)>1) {
  out <- colSums(as.matrix(Z[,ind]*Wt[,ind]), na.rm=T)/colSums(as.matrix(Wt[,ind]), na.rm=T)   
 }
 if(length(ind)==1) {
  out <- Z[,ind]
 }
out
}

subsetMat <- function(index,MARGIN,mat) {
    if(MARGIN==1)
         return(mat[index,])
    if(MARGIN==2)
         return(mat[,index])
}

get.coef.wt <- function(y,wt,W) {
 Rfast::lmfit(y=y,x=cbind(1,W),w=wt)$be
}

get.coef <- function(y,W) {
 Rfast::lmfit(y=y,x=cbind(1,W),w=wt)$be
}

get.coef2 <- function(RM_Z,signdev,alpha.dev,wt,RM_W,k.huber=1.345) {
 sigma.est<- mad(signdev)/0.6745
 wtvec    <- MASS::psi.huber(signdev,k=k.huber*sigma.est) * wt 
 wtvec    <- wtvec/mean(wtvec)
 TXW      <- sweep(t(RM_W),2,wtvec,'*')
 A        <- TXW %*% RM_W
 b        <- TXW %*% as.matrix(RM_Z - RM_W %*% as.matrix(alpha.dev))
return(A=cbind(A,b))
}

get.adev <- function(RM_Z,signdev,wt,RM_W,lambda.a,alpha.mean,k.huber=1.345) {
  sigma.est <- mad(signdev)/0.6745
  wtvec   <- MASS::psi.huber(signdev,k=k.huber*sigma.est) * wt 
  wtvec <- wtvec/mean(wtvec)
  TXW<- sweep(t(RM_W),2,wtvec,'*')
  A <- TXW %*% RM_W
  b <- TXW %*% as.matrix(RM_Z - RM_W %*% alpha.mean)
return(solve(A+lambda.a*diag(ncol(A)),b))
}

getW  <- function(Z,signdev,wt,alpha,ctl,k.huber=1.345,k=2,zeroinf=FALSE) {
  W <- rep(0,k)
  sigma.est <- mad(signdev[ctl])/0.6745
  wtvec     <- MASS::psi.huber(signdev[ctl],k=k.huber*sigma.est) * wt[ctl] 
  wtvec     <- wtvec/mean(wtvec)
  W[1]      <- Rfast::lmfit(y=Z[ctl],x=alpha[ctl,1],w=wtvec)$be
  if(k>1) {
    for(j in 2:k) {
      W[j] <- Rfast::lmfit(y=Z[ctl]-rowSums(as.matrix(alpha[ctl,1:(j-1)]) %*% as.matrix(W[1:(j-1)])),x=alpha[ctl,j],w=wtvec)$be
    } 
  }	   
 W
}

# estimate lambda.a parameter using RIDGM approach (Dempster et al 1977).
est.lambda.a <- function(alpha.dev,signdev,wt,RM_W,k.huber=1.345) {
 k <- length(alpha.dev)
 lambda.a.root <- function(lambda.a,alpha.dev,signdev,wt,RM_W,k.huber=k.huber) {
  sigma.est <- mad(signdev)/0.6745
  wtvec   <- MASS::psi.huber(signdev,k=k.huber*sigma.est) * wt 
  wtvec <- wtvec/mean(wtvec)
  TXW<- sweep(t(RM_W),2,wtvec,'*')
  A <- TXW %*% RM_W
  if(k>1)
   diff <- sum(alpha.dev^2/diag(solve(A+lambda.a*diag(k)))) - k
  if(k==1)
   diff <- sum(alpha.dev^2 * c(A+lambda.a)) - k

 diff
 }

 lambda.a.est <- tryCatch({uniroot(lambda.a.root,interval=c(0.001,10),extendInt="yes",alpha.dev=alpha.dev,signdev=signdev,wt=wt,RM_W=RM_W,k.huber=k.huber)$root},
 			    error = function(e) {NA})
 lambda.a.est
}  

     
# ZINB logl to estimate pi0
logl.zinb <- function(p,y,mu,size) {
  -sum(ZIM::dzinb(y,lambda=mu,k=size,omega=1/(1+exp(-p)),log=TRUE))
}



