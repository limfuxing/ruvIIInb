#' Removing unwanted variation from sequencing count data 
#'
#' This function performs fast version of ruvIIInb method to remove unwanted variation from sequencing count data. Currently, only Negative Binomial model is implemented. It takes raw count matrix with features (genes/transcript in scRNA-seq) as rows and samples (cells in scRNA-seq) as columns. 
#' Users need to specify the set of negative control genes (i.e genes where the between cells variation is assumed to be solely due to the unwanted variation) 
#' and the (pseudo)replicate matrix that define sets of cells that can be considered technical replicates.
#'
#' @param Y raw count as DelayedMatrix object with genes as rows and cells as columns
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
#'        interest is not associated with batch. If they are associated, this strategy will risk removing too much biology. Specifying the biological factor as 'strata' variable helps reducing the amount of biology removed.
#' @param batch.disp whether to fit batch-speficic dispersion parameter for each gene. The default is FALSE.
#' @param pCells.touse the proportion of cells used to estimate alpha and dispersion parameter (for speed-up). Default=20%.
  
#' @return A list containing the raw data (as sparse matrix), the unwanted factors and regression coefficients associated with the unwanted factors.
#' @export
fastruvIII.nb <- function(Y,M,ctl,k=2,robust=FALSE,ortho.W=FALSE,lambda.a=0.01,lambda.b=16,batch=NULL,step.fac=0.5,inner.maxit=50,outer.maxit=25,
        ncores=2,use.pseudosample=FALSE,nc.pool=20,strata=NULL,batch.disp=FALSE,pCells.touse=0.2) {

# register parallel backend
#register(BPPARAM)

# setup cluster for doParallel
doParallel::registerDoParallel(ncores)
BPPARAM=BiocParallel::DoparParam()
require(doParallel)
require(foreach)
parallel <- as.logical(ncores>1)

if(!is.logical(ctl)) 
  ctl <- rownames(Y) %in% ctl
idx <- which(ctl)
nctl<- sum(ctl)

# force M matrix into logical
mode(M) <- 'logical'

#format and convert M matrix into list of indices to save space
rep.ind <- apply(M,2,which)

# force Y to DelayedMatrix object
if(class(Y)!='DelayedMatrix')
  Y <- DelayedArray::DelayedArray(Y)

# Winsorize gene-by-gene and calculate size-factor
bef <- Sys.time()
max.val  <- ceiling(DelayedMatrixStats::rowQuantiles(Y,prob=0.995))
winsorize<- Y>max.val
Y <- Y*(1-winsorize) + winsorize*max.val
aft <- Sys.time()
print(paste0('time to Winsorize count matrix:', difftime(aft,bef,units='secs')))

# remove genes with all zero counts (if any)
zero.genes <- DelayedMatrixStats::rowSums2(Y>0)==0
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

# use.pseudosample not available for HDF5Matrix input
if(use.pseudosample) {
  if(class(Y)=='DelayedMatrix') {
   print('Pseudo-sample (cell) feature is not available for DelayedMatrix input. Proceeding with use.pseudosample=FALSE')
   use.pseudosample <- FALSE
  }
}

# get Huber's k
k.huber<- ifelse(robust,1.345,100)

# select a subset of cells representative of sub-pops defined by M matrix
subsamples <- sample(rep.ind[[1]],size=ceiling(pCells.touse*length(rep.ind[[1]])))
if(length(rep.ind)>1) {
  for(i in 2:length(rep.ind)) {
    subsamples <- c(subsamples,sample(rep.ind[[i]],size=ceiling(pCells.touse*length(rep.ind[[i]]))))
  }
}
subsamples <- sort(subsamples)
Ysub  <- as.matrix(Y[,subsamples])

off.g  <- log(DelayedMatrixStats::rowMeans2(Ysub))

# initiate W by performing SVD on the control log(Yc+1)
Y.norm <- sweep(Ysub,1,exp(off.g),'/')
# determine the number of blocks
#nblock <- ceiling(ns*ngene/1e+08)
#cell.blocks <- round(seq(0,ns,length=nblock+1))
#gene.blocks <- round(seq(0,ngene,length=nblock+1))
Wsub <- irlba::irlba(log(as.matrix(Y.norm[ctl,])+1),nv=k)$v
Msub <- as.matrix(M[subsamples,])
rep.sub <- apply(Msub,2,which)
sub.batch <- batch[subsamples]
nsub <- ncol(Ysub)

# initiate beta and alpha: by performing poisson GLM gene-by-gene
nW   <- ncol(Wsub)
nM   <- ncol(M)
ngene<- nrow(Y)
ns   <- ncol(Y)
ncells <- rep(1,ncol(Y))

bef   <- Sys.time()
coef <- foreach(i=1:ngene, .combine=cbind, .packages="Rfast") %dopar% {
          out <- Rfast::lmfit(y=log(Ysub[i,]+1),x=cbind(1,Wsub),w=Ysub[i,]+1)$be
          return(as.matrix(out))
}
#print(dim(coef))
aft   <- Sys.time()
#print(paste0('time to compute initial est:', difftime(aft,bef,units='secs')))
coef <- as.matrix(coef)

# alpha(n x k) matrix
alpha<- as.matrix(t(coef)[,-1])
alpha.mean <- matrixStats::colMeans2(alpha)
alpha.dev <- sweep(alpha,2,alpha.mean,'-')

#print(dim(alpha))
#print(dim(Wsub))

alpha.mean <- matrixStats::colMeans2(alpha)
alpha.dev <- sweep(alpha,2,alpha.mean,'-')

# Mb (n x m_1) matrix
Mb   <- matrix(0,nrow=ngene,ncol=nM,byrow=FALSE)

#handling NA/Inf
Mb[which(is.na(Mb) | is.infinite(Mb))] <- 0
alpha[which(is.na(alpha) | is.infinite(alpha))] <- 0

#gmean
gmean<- c(coef[1,])

# estimate initial psi 
Wa   <- Matrix::tcrossprod(alpha,Wsub) 
offs <- Wa + Mb[,apply(Msub,1,which)] 
bef  <- Sys.time()
if(parallel) {
 if(is.null(batch)) 
  psi <- estimateDisp.par(Ysub,as.matrix(rep(1,nsub)),offset=offs,tagwise=TRUE,robust=TRUE,BPPARAM=BPPARAM)$tagwise.dispersion
 if(!is.null(batch)) {
  psi <- estimateDisp.par(Ysub[,sub.batch==1],as.matrix(rep(1,sum(sub.batch==1))),
				offset=offs[,sub.batch==1],tagwise=TRUE,robust=TRUE,BPPARAM=BPPARAM)$tagwise.dispersion
  for(B in 2:max(batch)) 
      psi <- cbind(psi,estimateDisp.par(Ysub[,sub.batch==B],as.matrix(rep(1,sum(sub.batch==B))),
				offset=offs[,sub.batch==B],tagwise=TRUE,robust=TRUE,BPPARAM=BPPARAM)$tagwise.dispersion)
 }
} 
if(!parallel) {
 if(is.null(batch)) 
  psi <- edgeR::estimateDisp(Ysub,as.matrix(rep(1,nsub)),offset=offs,tagwise=TRUE,robust=TRUE)$tagwise.dispersion
 if(!is.null(batch)) {
  psi <- edgeR::estimateDisp(Ysub[,sub.batch==1],as.matrix(rep(1,sum(sub.batch==1))),
				offset=offs[,sub.batch==1],tagwise=TRUE,robust=TRUE)$tagwise.dispersion
  for(B in 2:max(batch)) 
   psi <- cbind(psi,edgeR::estimateDisp(Ysub[,sub.batch==B],as.matrix(rep(1,sum(sub.batch==B))),
				offset=offs[,sub.batch==B],tagwise=TRUE,robust=TRUE)$tagwise.dispersion)
 }
}
aft <- Sys.time()
#print(paste0('time to calculate psi:',  difftime(aft,bef,units='secs')))

# initiate
conv <- FALSE
iter.outer <- 0
logl.outer <- NULL
step <- rep(1,ngene)
print('Start...')
while(!conv) {
 iter.outer <- iter.outer + 1
 logl.beta <- NULL
 
 ### calculate initial loglik
 lmu.hat    <- gmean + Mb[,apply(Msub,1,which)] +  Matrix::tcrossprod(alpha,Wsub) 
 if(!is.null(batch)) {
  temp <- foreach(i=1:ngene, .combine=c) %dopar% {
    tmp <- dnbinom(Ysub[i,],mu=exp(lmu.hat[i,]),size=1/psi[i,sub.batch],log=TRUE)
    tmp[is.infinite(tmp) | is.na(tmp)] <- log(.Machine$double.xmin)
    return(sum(tmp))
  }
 } 
 if(is.null(batch)) {
  temp <- foreach(i=1:ngene, .combine=c) %dopar% {
    tmp <- dnbinom(Ysub[i,],mu=exp(lmu.hat[i,]),size=1/psi[i],log=TRUE)
    tmp[is.infinite(tmp) | is.na(tmp)] <- log(.Machine$double.xmin)
    return(sum(tmp))
  }
 } 
 loglik <- temp

 #############################
 step <- rep(1,ngene)
 conv.beta <- FALSE
 # saving temporary best estimate
 best.W <- Wsub ; best.psi <- psi ; best.a <- alpha; best.Mb <- Mb; best.adev <- alpha.dev; best.gmean <- gmean
 
 iter <- 1
 halving <- 0
 while(!conv.beta) {
 print(paste0('Inner iter:', iter))
 # calc working vector Z for a subset of cells 
 lmu.hat    <- gmean + Mb[,apply(Msub,1,which)] + Matrix::tcrossprod(alpha,Wsub) 

 # weight based on NB
 if(!is.null(batch))
  sig.inv<- 1/(exp(-lmu.hat) + sweep(psi[,sub.batch],2,rep(1,nsub),'/'))
 if(is.null(batch))
  sig.inv<- 1/(exp(-lmu.hat) +outer(psi,rep(1,nsub),'/'))
 
 # calculate signed deviance
 if(!is.null(batch) & robust) {
   signdev <- sign(Ysub-exp(lmu.hat))* sqrt(2*(dnbinom(Ysub,mu=Ysub,size=sweep(1/psi[,sub.batch],2,rep(1,nsub),'*'),log=TRUE) - 
						dnbinom(Ysub,mu=exp(lmu.hat),size=sweep(1/psi[,sub.batch],2,rep(1,nsub),'*'),log=TRUE) ))
 }
 if(is.null(batch) & robust) {
   signdev <- sign(Ysub-exp(lmu.hat))* sqrt(2*(dnbinom(Ysub,mu=Ysub,size=outer(1/psi,rep(1,nsub),'*'),log=TRUE) - 
						dnbinom(Ysub,mu=exp(lmu.hat),size=outer(1/psi,rep(1,nsub),'*'),log=TRUE) ))
 }  
 
 
 Z      <- lmu.hat + step*((Ysub+0.01)/(exp(lmu.hat)+0.01) - 1)

 # project Z and W into M
 
 # step1: obtain projecion of Z and W onto the range of M and get the residual
 # RM_Z is (n x m) matrix, m=nsample,n=ngene
 # RM_W is (m x k) matrix
 proj.Z <- projection(rep.sub,Zmat=Z)[,apply(Msub,1,which)]
 RM_Z <- Z - proj.Z
 
 if(nW>1) 
   proj.tW <- projection(rep.sub,Zmat=t(Wsub))
 if(nW==1) 
  proj.tW <- t(as.matrix(projection.K1(rep.sub,Zmat=t(Wsub))))
 
 if(nW>1) 
  RM_W <- Wsub - t(proj.tW[,apply(Msub,1,which)])

 if(nW==1) 
  RM_W <- Wsub - as.matrix(proj.tW[,apply(Msub,1,which)])

 # store current alpha est
 alpha.old <- alpha

 # get alpha mean
 acomb <- function(...) abind::abind(..., along=3)
 bef   <- Sys.time()
 results <- foreach(i=1:ngene, .combine=acomb, .packages="Matrix") %dopar% {
  wtvec    <- sig.inv[i,] 
  wtvec    <- wtvec/mean(wtvec)
  A        <- Matrix::crossprod(RM_W*wtvec,RM_W)
  b1       <- Matrix::crossprod(RM_W*wtvec, as.matrix(RM_Z[i,]) - RM_W %*% as.matrix(alpha.dev[i,]))
  return(cbind(A,b1))
 }
 aft   <- Sys.time()
 #print(paste0('time to calculate ameans:',  difftime(aft,bef,units='secs')))
 #print(dim(results))

 if(k>1) {
  alpha.mean <- solve(as.matrix(apply(results[,1:nW,],c(1,2),sum)),as.matrix(apply(results[,nW+1,],1,sum)))
  A  <- results[,1:nW,,drop=FALSE]
 }
 if(k==1) {
  alpha.mean <- c(sum(results[,(nW+1),])/sum(results[,1:nW,]))
  A <- results[,1:nW,,drop=FALSE]
 }
 #print(dim(A))
 #print(dim(b))
 # get alpha dev
 acomb2 <- function(...) abind::abind(..., along=2)
 bef=Sys.time()
 RM_W_amean<- RM_W %*% as.matrix(alpha.mean)
 alpha.dev <- foreach(i=1:ngene, .combine=acomb2, .packages="Matrix") %dopar% { 
      wtvec    <- sig.inv[i,] 
      wtvec    <- wtvec/mean(wtvec)
      b        <- Matrix::crossprod(RM_W*wtvec, as.matrix(RM_Z[i,]) - RM_W_amean)
      solve(A[,,i]+lambda.a*diag(nW),b)
 }
 aft <- Sys.time()
 #print(dim(alpha.dev))
 #print(paste0('time to calculate adev:',  difftime(aft,bef,units='secs')))
 alpha.dev <- as.matrix(alpha.dev)
 if(ncol(alpha.dev)!=k)
  alpha.dev <- t(alpha.dev)

 # get alpha
 alpha <- sweep(alpha.dev,2,alpha.mean,'+')
 # if new alpha est has missing/inf values, revert to previous estimate
 #print(paste0('Problematic alpha entries=',sum(is.na(alpha) | is.infinite(alpha))))
 if(any(is.na(alpha) | is.infinite(alpha) )) alpha <- alpha.old

 # reduce outliers
 a.med <- matrixStats::colMedians(alpha)
 a.mad <- matrixStats::colMads(alpha)
 for(i in 1:nW) {
   alpha[which(alpha[,i]> (a.med[i]+4*a.mad[i])),i] <- a.med[i]+4*a.mad[i]
   alpha[which(alpha[,i]< (a.med[i]-4*a.mad[i])),i] <- a.med[i]-4*a.mad[i]
 }


 # step2: update W using WLS of Z_c (Z for control genes) on alpha.c (alpha for control genes)
 Wsub.old <- Wsub 

 bef <- Sys.time()
 #W   <- BiocParallel::bpmapply(FUN=getW,Z=Z.list,signdev=signdev.bycol,wt=wtlist,zeroinf=zeroinf.list,
 #			MoreArgs=list(alpha=alpha,ctl=ctl,k.huber=k.huber,k=k),BPPARAM=BPPARAM)


 Wsub <- foreach(i=1:nsub, .combine=cbind, .packages="Matrix") %dopar% { 
  out <- rep(0,k)
  wtvec     <- sig.inv[ctl,i] 
  wtvec     <- wtvec/mean(wtvec)
  out[1]      <- Rfast::lmfit(y=Z[ctl,i]-gmean[ctl],x=alpha[ctl,1],w=wtvec)$be
  if(k>1) {
    for(j in 2:k) {
      out[j] <- Rfast::lmfit(y=Z[ctl,i]-gmean[ctl]-matrixStats::rowSums2(as.matrix(alpha[ctl,1:(j-1)]) %*% as.matrix(out[1:(j-1)])),x=alpha[ctl,j],w=wtvec)$be
    } 
  }
  as.matrix(out)
 }
	   
 aft <- Sys.time()
 #print(paste0('time to calculate W:',  difftime(aft,bef,units='secs')))
 Wsub <- as.matrix(Wsub)
 if(ncol(Wsub)!=k) 
  Wsub <- t(Wsub)

 colNorm.Wsub <- matrixStats::colSums2(Wsub^2)
 # normalize W to have unit length for each col
 Wsub  <- sweep(Wsub,2,sqrt(colNorm.Wsub),'/')
 
 # orthogonalize W (optional)
 if(nW>1 & ortho.W)  
   Wsub <- matlib::GramSchmidt(Wsub)

 # if new W has missing/inf values, revert to previous estimate
 #print(paste0('Problematic W entries=',sum(is.na(Wsub) | is.infinite(Wsub))))
 if(any(is.na(Wsub) | is.infinite(Wsub) )) Wsub <- W.old
 
 #step3a: update gmean
 Z.res  <- Z -  Matrix::tcrossprod(alpha,Wsub)  - Mb[,apply(Msub,1,which)]

 bef <- Sys.time()
 temp <- tryCatch({
 for(i in 1:ngene) {
    wtvec     <- sig.inv[i,]
    wtvec     <- wtvec/sum(wtvec)
    gmean[i]  <- sum(Z.res[i,]*wtvec)
 }
 }, error = function(e) {NA})
 aft <- Sys.time()
 #print(paste0('time to calculate gmean:',  difftime(aft,bef,units='secs')))

 #step 3b: update Mb
 Z.res  <- Z -  Matrix::tcrossprod(alpha,Wsub)  - gmean
 Mb.old <- Mb
 bef <- Sys.time()
 Mb  <- tryCatch({ projection(rep.sub,Zmat=Z.res,Wmat=sig.inv,lambda=lambda.b)}, error = function(e) {NA})
 aft <- Sys.time()
 #print(paste0('time to calculate Mb:',  difftime(aft,bef,units='secs')))

 # if new Mb has missing/inf values, revert to previous estimate
 #print(paste0('Problematic Mb entries=',sum(is.na(Mb) | is.infinite(Mb))))
 if(any(is.na(Mb))) Mb <- Mb.old
 Mb[which(is.na(Mb) | is.infinite(Mb))] <- 0

 # reduce outliers
 Mb.med <- matrixStats::rowMedians(Mb)
 Mb.mad <- matrixStats::rowMads(Mb)
 for(i in 1:ngene) {
   Mb[i,which(Mb[i,]> (Mb.med[i]+4*Mb.mad[i]))] <- Mb.med[i]+4*Mb.mad[i]
   Mb[i,which(Mb[i,]< (Mb.med[i]-4*Mb.mad[i]))] <- Mb.med[i]-4*Mb.mad[i]
 }

 # calculate current logl
 lmu.hat    <- gmean + Mb[,apply(Msub,1,which)] +  Matrix::tcrossprod(alpha,Wsub) 
 if(!is.null(batch)) {
  temp <- foreach(i=1:ngene, .combine=c) %dopar% {
    tmp <- dnbinom(Ysub[i,],mu=exp(lmu.hat[i,]),size=1/psi[i,sub.batch],log=TRUE)
    tmp[is.infinite(tmp) | is.na(tmp)] <- log(.Machine$double.xmin)
    return(sum(tmp))
  }
 } 
 if(is.null(batch)) {
  temp <- foreach(i=1:ngene, .combine=c) %dopar% {
    tmp <- dnbinom(Ysub[i,],mu=exp(lmu.hat[i,]),size=1/psi[i],log=TRUE)
    tmp[is.infinite(tmp) | is.na(tmp)] <- log(.Machine$double.xmin)
    return(sum(tmp))
  }
 } 
 loglik.tmp <- temp 

 # check degenerate case
 if(iter>=1) {
  degener<- sum(loglik)>sum(loglik.tmp)
  degener[is.na(degener)] <- TRUE
  if(degener) {
    check.gene <- loglik>loglik.tmp
    step[check.gene] <- step[check.gene]*step.fac
    Wsub    <- best.W; alpha <- best.a ; Mb <- best.Mb ; psi <- best.psi ; alpha.dev <- best.adev ; gmean <- best.gmean 
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
   best.W <- Wsub ; best.psi <- psi ; best.a <- alpha; best.Mb <- Mb ; best.adev <- alpha.dev; best.gmean <- gmean 
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

Wa =  Matrix::tcrossprod(best.a,best.W) 
offs <- Wa + best.Mb[,apply(Msub,1,which)] 
if(updt.psi) {
if(parallel) {
 if(is.null(batch)) 
  psi.new <- tryCatch({ estimateDisp.par(Ysub,as.matrix(rep(1,nsub)),offset=offs,
				tagwise=TRUE,robust=TRUE,BPPARAM=BPPARAM)$tagwise.dispersion}, error = function(e) {NA})
 if(!is.null(batch)) {
  psi.new <- psi
  for(B in 1:max(batch)) 
   psi.new[,B]<-tryCatch({estimateDisp.par(Ysub[,sub.batch==B],as.matrix(rep(1,sum(sub.batch==B))),
				offset=offs[,sub.batch==B],tagwise=TRUE,robust=TRUE,BPPARAM=BPPARAM)$tagwise.dispersion},error= function(e) {NA})
 }
}

if(!parallel) {
 if(is.null(batch)) 
  psi.new <- tryCatch({ edgeR::estimateDisp(Ysub,as.matrix(rep(1,nsub)),offset=offs,tagwise=TRUE,robust=TRUE)$tagwise.dispersion},error = function(e) {NA})
 if(!is.null(batch)) {
  psi.new <- psi
  for(B in 1:max(batch)) 
   psi.new[,B] <- tryCatch({ edgeR::estimateDisp(Ysub[,sub.batch==B],as.matrix(rep(1,sum(sub.batch==B))),
				offset=offs[,sub.batch==B],tagwise=TRUE,robust=TRUE)$tagwise.dispersion},error = function(e) {NA})
 }
}
# update psi
psi[!is.na(psi.new) & !is.infinite(psi.new)] <- psi.new[!is.na(psi.new) & !is.infinite(psi.new)]
}

# record outer logl
Wsub <- best.W
alpha <- best.a
Mb    <- best.Mb
alpha.dev <- best.adev
gmean <- best.gmean
# new changes
best.psi <- psi

# recalculate logl
lmu.hat    <- gmean + Mb[,apply(Msub,1,which)] +  Matrix::tcrossprod(alpha,Wsub) 
if(!is.null(batch)) {
  temp <- foreach(i=1:ngene, .combine=c) %dopar% {
    tmp <- dnbinom(Ysub[i,],mu=exp(lmu.hat[i,]),size=1/psi[i,sub.batch],log=TRUE)
    tmp[is.infinite(tmp) | is.na(tmp)] <- log(.Machine$double.xmin)
    return(sum(tmp))
  }
} 
if(is.null(batch)) {
  temp <- foreach(i=1:ngene, .combine=c) %dopar% {
    tmp <- dnbinom(Ysub[i,],mu=exp(lmu.hat[i,]),size=1/psi[i],log=TRUE)
    tmp[is.infinite(tmp) | is.na(tmp)] <- log(.Machine$double.xmin)
    return(sum(tmp))
  }
} 
logl.outer <- c(logl.outer, sum(temp)) #-lambda.a*sum(scale(alpha,scale=FALSE)^2)-lambda.b*sum(Mb^2))
print(paste0('Outer Iter ',iter.outer, ' logl-likelihood:', logl.outer[length(logl.outer)]))

if(iter.outer>1) {
  conv <- ifelse( (logl.outer[iter.outer] - logl.outer[iter.outer-1])/abs(logl.outer[iter.outer-1]) < 1e-04 | iter.outer>=outer.maxit,TRUE,FALSE)
}

if(conv) {
 alpha <- best.a
 gmean <- best.gmean
 Mb  <- best.Mb
 psi   <- best.psi
 # estimate W for all samples (initial)
 print('Estimating W for all samples...')
 bef=Sys.time()
 W <- matrix(0,ncol(Y),nW)
 W[subsamples,] <- Wsub
 W <- foreach(i=1:ns, .combine=cbind, .packages="Matrix") %dopar% { 
  out <- rep(0,k)
  lmu  <- gmean + Mb[,which(M[i,])] + alpha %*% as.matrix(W[i,])
  if(!is.null(batch))
   wtvec<- c(1/(exp(-lmu) + psi[,batch[i]]))
  if(is.null(batch))
   wtvec<- c(1/(exp(-lmu) + psi))
  wtvec <- wtvec/mean(wtvec)

  Zctl   <- lmu[ctl] + ( (Y[ctl,i]+0.01)/(exp(lmu[ctl])+0.01) - 1)
  out[1] <- Rfast::lmfit(y=Zctl-gmean[ctl],x=alpha[ctl,1],w=wtvec[ctl])$be
  if(k>1) {
    for(j in 2:k) {
      out[j] <- Rfast::lmfit(y=Zctl-gmean[ctl]-matrixStats::rowSums2(as.matrix(alpha[ctl,1:(j-1)]) %*% as.matrix(out[1:(j-1)])),x=alpha[ctl,j],w=wtvec[ctl])$be
    } 
  }
  as.matrix(out)
 }
 W <- as.matrix(W)
 if(ncol(W)!=k) 
  W <- t(W)

 # normalize W and alpha
 calnorm.W <- sqrt(matrixStats::colSums2(W^2))
 #W     <- sweep(W,2,calnorm.W,'/')
 #alpha <- sweep(alpha,2,calnorm.W,'*')

 # orthogonalize W (optional)
 if(nW>1 & ortho.W)  
   W <- matlib::GramSchmidt(W)

 # update W until convergence
 conv.W <- FALSE
 while(!conv.W) {
  W.old <- W
  W <- foreach(i=1:ns, .combine=cbind, .packages="Matrix") %dopar% { 
   out <- rep(0,k)
   lmu  <- gmean + Mb[,which(M[i,])] + alpha %*% as.matrix(W[i,])
   Zctl<-  lmu[ctl] + ( (Y[ctl,i]+0.01)/(exp(lmu[ctl])+0.01) - 1)
   if(!is.null(batch))
    wtvec<- c(1/(exp(-lmu) + psi[,batch[i]]))
   if(is.null(batch))
    wtvec<- c(1/(exp(-lmu) + psi))
   wtvec <- wtvec/mean(wtvec)
   out[1] <- Rfast::lmfit(y=Zctl-gmean[ctl],x=alpha[ctl,1],w=wtvec[ctl])$be
   if(k>1) {
     for(j in 2:k) {
       out[j] <- Rfast::lmfit(y=Zctl-gmean[ctl]-matrixStats::rowSums2(as.matrix(alpha[ctl,1:(j-1)]) %*% as.matrix(out[1:(j-1)])),x=alpha[ctl,j],w=wtvec[ctl])$be
     } 
   }
   as.matrix(out)
  }
  W <- as.matrix(W)
  if(ncol(W)!=k) 
   W <- t(W)

  # normalize W and alpha
  calnorm.W <- sqrt(matrixStats::colSums2(W^2))
  #W     <- sweep(W,2,calnorm.W,'/')
  #alpha <- sweep(alpha,2,calnorm.W,'*')

  # orthogonalize W (optional)
  if(nW>1 & ortho.W)  
    W <- matlib::GramSchmidt(W) 

  crit.W <- mean( (abs(W-W.old)/abs(W.old))^2)
  #print(round(crit.W,10))
  conv.W <- crit.W< 1e-6
 }
 aft=Sys.time()
 #print(paste0('Time to estimate W for all samples:',difftime(aft,bef,units='secs')))
}
} # end of outer IRLS loop

# calibrate alpha and W
for(k in 1:nW) {
 ab.W     <- Rfast::lmfit(y=Wsub[,k],x=cbind(1,W[subsamples,k]))$be
 ab.alpha <- Rfast::lmfit(y=best.a[,k],x=cbind(1,alpha[,k]))$be
 W[,k]    <- ab.W[1] + ab.W[2]*W[,k]
 alpha[,k]<- ab.alpha[1] + ab.alpha[2]*alpha[,k]
}
#alpha <- sweep(alpha,2,scfac.alpha,'/')
cor.check <- diag(cor(Wsub,W[subsamples,]))
#output
return( list("counts"=Y,"W"=W, "Wsub"=Wsub,"M"=M, "ctl"=ctl, "logl"=logl.outer, "a"=alpha,"asub"=best.a,"Mb"=Mb, "gmean"=gmean,
		"psi"=psi,'L.a'=lambda.a,'L.b'=lambda.b,batch=batch,sub=subsamples,corW=cor.check) )
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
   out <- matrixStats::rowMeans2(as.matrix(Z[,ind]), na.rm=T)
  }
  if(length(ind)==1) {
    out <- Z[,ind]
  }
out
}

rowWtAvg <- function(ind,Z,Wt,lambda=0) {
 if(length(ind)>1) {
  out <- rowSums(as.matrix(Z[,ind]*Wt[,ind]), na.rm=T)/(matrixStats::rowSums2(as.matrix(Wt[,ind]), na.rm=T) + lambda)  
 }
 if(length(ind)==1) {
  out <- c(Z[,ind])*c(Wt[,ind])/c(Wt[,ind]+lambda)
 }
out
}
 
rowAvg.K1 <- function(ind,Z) {
  if(length(ind)>1) {
   out <- matrixStats::colMeans2(as.matrix(Z[,ind]), na.rm=T)
  }
  if(length(ind)==1) {
    out <- Z[,ind]
  }
out
}

rowWtAvg.K1 <- function(ind,Z,Wt) {
 if(length(ind)>1) {
  out <- matrixStats::colSums2(as.matrix(Z[,ind]*Wt[,ind]), na.rm=T)/matrixStats::colSums2(as.matrix(Wt[,ind]), na.rm=T)   
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

get.coef2 <- function(RM_Z,signdev,alpha.dev,wt,RM_W,k.huber=1.345) {
 sigma.est<- mad(signdev)/0.6745
 if(sigma.est!=0)
  wtvec    <- MASS::psi.huber(signdev,k=k.huber*sigma.est) * wt 
 if(sigma.est==0)
  wtvec    <- wt 
 wtvec    <- wtvec/mean(wtvec)
 A        <- Matrix::crossprod(RM_W*wtvec,RM_W)
 b        <- Matrix::crossprod(RM_W*wtvec, as.matrix(RM_Z - RM_W %*% as.matrix(alpha.dev)))
return(A=cbind(A,b))
}

get.adev <- function(RM_Z,signdev,wt,RM_W,lambda.a,alpha.mean,k.huber=1.345) {
  sigma.est <- mad(signdev)/0.6745
  if(sigma.est==0)
    wtvec   <- wt 
  if(sigma.est!=0)
    wtvec   <- MASS::psi.huber(signdev,k=k.huber*sigma.est) * wt 
  wtvec <- wtvec/mean(wtvec)
  A <- Matrix::crossprod(RM_W*wtvec,RM_W)
  b <- Matrix::crossprod(RM_W*wtvec,as.matrix(RM_Z - RM_W %*% alpha.mean))
return(solve(A+lambda.a*diag(ncol(A)),b))
}

getW  <- function(Z,signdev,wt,alpha,ctl,k.huber=1.345,k=2,zeroinf=FALSE) {
  W <- rep(0,k)
  sigma.est <- mad(signdev[ctl])/0.6745
  if(sigma.est!=0)
   wtvec     <- MASS::psi.huber(signdev[ctl],k=k.huber*sigma.est) * wt[ctl] 
  if(sigma.est==0)
   wtvec     <- wt[ctl] 
  wtvec     <- wtvec/mean(wtvec)
  W[1]      <- Rfast::lmfit(y=Z[ctl],x=alpha[ctl,1],w=wtvec)$be
  if(k>1) {
    for(j in 2:k) {
      W[j] <- Rfast::lmfit(y=Z[ctl]-matrixStats::rowSums2(as.matrix(alpha[ctl,1:(j-1)]) %*% as.matrix(W[1:(j-1)])),x=alpha[ctl,j],w=wtvec)$be
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



