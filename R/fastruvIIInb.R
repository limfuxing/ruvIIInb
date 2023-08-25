#' Removing unwanted variation from sequencing count data 
#'
#' This function performs fast version of ruvIIInb method to remove unwanted variation from sequencing count data. Currently, only Negative Binomial model for UMI data is supported. It takes raw count matrix with features (genes/transcript in scRNA-seq) as rows and samples (cells in scRNA-seq) as columns. 
#' Users need to specify the set of negative control genes (i.e genes where the between cells variation is assumed to be solely due to the unwanted variation) 
#' and the (pseudo)replicate matrix that define sets of cells that can be considered technical replicates.
#'
#' @param Y raw count as DelayedMatrix object with genes as rows and cells as columns
#' @param M replicate matrix with number of rows equal to the number of cells and number of columns equal to the number of distinct sub-population of cells: M(i,j) = 1 if cell i is part of sub-population j and 0 otherwise. If a row has all zero entries, it means the corresponding cell is not assumed to belong to any of the known sub-populations (a priori unannotated cells). 
#' @param ctl  either logical vector of length n with TRUE = control genes OR a character vector of names of control genes
#' @param k dimension of unwanted variation against which the data is to be normalized against. Default is 2.
#' @param robust logical value. Whether to use Huber's robust weight in the iterative reweighted least squares (IRLS) algorithm. Default is FALSE.
#' @param ortho.W logical value. Whether the unwanted factors (W) are to be orthogonalized. The default is FALSE.
#' @param lambda.a smoothing parameter for regularizing regression coefficients associated with W. The default is 0.01
#' @param lambda.b smoothing parameter for regularizing the gene-level intercepts. The default is 16.
#' @param batch numeric vector containing batch information for each cell. Must correspond to columns of count matrix. If not supplied, cells are assumed to come from one batch.
#' @param step.fac multiplicative factor to decrease IRLS step by when log-likelihood diverges. Default is 0.5
#' @param inner.maxit the maximum number of IRLS iteration for estimating mean parameters for a given dispersion parameter. Default is 50
#' @param outer.maxit the maximum number of IRLS iteration. Default is 25 
#' @param ncores The number of cores used for parallel processing. The default number of workers (cores) is 2.
#' @param use.pseudosample whether to use pseudocells to define replicates (default is FALSE). Note that the replicates defined by the pseudocells will be used in addition to any replicates defined by the M matrix above. We recommend that use.pseudosample=TRUE be used only for data with UMI.
#' @param nc.pool the number of cells per pool (used for defining pseudocells). The default is 20. Only relevant when use.pseudosample=TRUE.
#' @param batch.disp whether to fit batch-specific dispersion parameter for each gene. The default is FALSE.
#' @param pCells.touse the proportion of a priori annotated cells used to estimate alpha and dispersion parameter (Default=20%). When pseudocells are used (use.pseudosample=TRUE), ultimately only 10% of the original subset of cells will be used to estimate alpha.
  
#' @return A list containing the raw data (as sparse matrix), the unwanted factors, regression coefficients associated with the unwanted factors, the M matrix and the estimates of dispersion parameters
#' @export
fastruvIII.nb <- function(Y,M,ctl,k=2,robust=FALSE,ortho.W=FALSE,lambda.a=0.01,lambda.b=16,batch=NULL,step.fac=0.5,inner.maxit=50,outer.maxit=25,
        ncores=2,use.pseudosample=FALSE,nc.pool=20,batch.disp=FALSE,pCells.touse=0.2) {

# register parallel backend
#register(BPPARAM)

# setup cluster for doParallel
doParallel::registerDoParallel(ncores)
BPPARAM=BiocParallel::DoparParam()
require(doParallel,quietly=TRUE)
require(foreach,quietly=TRUE)
parallel <- as.logical(ncores>1)

if(!is.logical(ctl)) 
  ctl <- rownames(Y) %in% ctl
idx <- which(ctl)
nctl<- sum(ctl)

# check if any rows of M matrix has zero sum 
zero.Mrows <- any(rowSums(M)==0)
if(zero.Mrows) {
 n.unannot <- sum(rowSums(M)==0)
 print(paste0(n.unannot, ' cells have all zero entries in the M matrix. They will be considered as un-annotated cells'))
 if(n.unannot < (0.01*ncol(Y))) 
   stop('Less than 1% of cells have known annotation. Please increase the pct of cells with known annotation to 1%')
}


# check if any cols of M matrix has zero sum 
zero.Mcols <- any(colSums(M)==0)
if(zero.Mcols) {
  stop(paste0('Columns ', which(colSums(M)==0),' of the M matrix has zero sum. Pls remove these columns and re-run'))
}

# force M matrix into logical matrix
mode(M) <- 'logical'
M       <- as.matrix(M)

#format and convert M matrix into list of indices to save space
rep.ind <- apply(M,2,which,simplify=FALSE)
# define strata as columns of M matrix for which a cell is part of
strata <- apply(M,1,which)

# force Y to DelayedMatrix object
if(!any(class(Y)=='DelayedMatrix'))
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

# if no batch info is supplied, assumes cells come from one batch
if(is.null(batch)) {
  print(paste0('Batch variable not supplied...assuming cells come from one batch'))
  batch <- rep(1,ncol(Y))
}

# if batch is not numeric vector, stop and let user knows
if(!is.null(batch) & !is.numeric(batch)) {
  stop(paste0('Batch variable is not a numeric vector. Please supply batch variable as numeric vector with 1,2...B values'))
}


# get Huber's k
k.huber<- ifelse(robust,1.345,100)

# adjust pCells.touse to allow max 3000 cells
max.pCells <- 3000/sum(rowSums(M)>0)
if(max.pCells < pCells.touse) {
  pCells.touse <- max.pCells
  print(paste0('pCells.touse parameter is too large and has now been set to ', round(pCells.touse ,6)))
}
# select a subset of cells representative of sub-pops defined by M matrix
subsamples <- sample(rep.ind[[1]],size=ceiling(pCells.touse*length(rep.ind[[1]])))
if(length(rep.ind)>1) {
  for(i in 2:length(rep.ind)) {
    subsamples <- c(subsamples,sample(rep.ind[[i]],size=ceiling(pCells.touse*length(rep.ind[[i]]))))
  }
}
subsamples <- sort(subsamples)
Ysub  <- as.matrix(Y[,subsamples])
Msub <- as.matrix(M[subsamples,])
rep.sub <- apply(Msub,2,which,simplify=FALSE)
sub.batch <- batch[subsamples]
if(!is.null(strata)) 
 strata <- strata[subsamples]

# adding pseudosamples
sf <- colSums(Ysub)
sf <- sf/mean(sf) 
sample.type <- rep('sc',ncol(Ysub))
# number of batches before pseudocells addition
nbatch.org <- ifelse(is.null(sub.batch),1,length(unique(sub.batch)))
# if batch.disp=FALSE,force nbatch.org==1
nbatch.org <- ifelse(!batch.disp,1,nbatch.org)

# number of single cells before pseudocells addition
nsub <- ncol(Ysub)
if(use.pseudosample) {
 if(is.null(strata)) {
   batch.ext <- sub.batch
   for(btc in unique(sub.batch)) {
    nq  <- max(round(sum(batch==btc)/nc.pool),1)
    qsf <-  as.numeric(gtools::quantcut(sf[sub.batch==btc],q=nq))
    for(q in 1:max(qsf)) {
       pool <- Ysub[,1:nsub][,sub.batch==btc][,qsf==q]
       dns.pool <- rbinom(nrow(pool),size=rowSums(pool),prob=1/ncol(pool))
       Ysub <- cbind(Ysub,dns.pool)
       sample.type <- c(sample.type,'sc')
       batch.ext  <- c(batch.ext,btc+max(sub.batch))
       Msub      <- rbind(Msub,FALSE)
    }
   }
 ps.mcol <- rep(FALSE,nrow(Msub))
 ps.mcol[seq(nsub+1,ncol(Ysub),1)] <- TRUE
 Msub <- cbind(Msub,ps.mcol)
 sub.batch <- batch.ext
 }

 if(!is.null(strata)) {
  batch.ext <- sub.batch
  for(st in unique(na.omit(strata))) {
    for(btc in unique(sub.batch)) {
     nq  <- max(round(sum(sub.batch==btc & strata==st,na.rm=TRUE)/nc.pool),1)
     if(sum(sub.batch==btc & strata==st, na.rm=TRUE)>1) {
      qsf <-  na.omit(as.numeric(gtools::quantcut(sf[sub.batch==btc & strata==st],q=nq)))
      for(q in 1:max(qsf)) {
       pool <- Ysub[,1:nsub][,sub.batch==btc & strata==st & !is.na(strata)][,qsf==q]
       dns.pool <- rbinom(nrow(pool),size=rowSums(pool),prob=1/ncol(pool))
       Ysub <- cbind(Ysub,dns.pool)
       sample.type <- c(sample.type,'sc')
       batch.ext  <- c(batch.ext,btc+max(sub.batch))
       ps.mrow <- rep(FALSE,ncol(Msub)) 
       ps.mrow[st] <- TRUE
       Msub    <- rbind(Msub,ps.mrow)
      }
     } #if
    } #for batch
  } # for st
 if(!batch.disp) 
   batch.ext <- ifelse(batch.ext>max(sub.batch),2,1)
   
 sub.batch <- batch.ext

 } # if
}
nbatch <- ifelse(!use.pseudosample,nbatch.org,length(unique(sub.batch)))

subsamples.org     <- subsamples
subsubsamples.org  <- 1:nsub
if(use.pseudosample) {
 # here, we take only 10% of the original subsampled cells + all pseudocells
 subsamples.org     <- subsamples
 subsubsamples.org  <- sort(sample(nsub,size=round(0.1*nsub)))
 subsamples         <- c(subsubsamples.org,(nsub+1):ncol(Ysub))

 Ysub <- Ysub[,subsamples]
 Msub <- as.matrix(Msub[subsamples,])
 rep.sub <- apply(Msub,2,which,simplify=FALSE)
 sub.batch <- sub.batch[subsamples]
 nsub <- ncol(Ysub)
}

# select a maximum of 100 cells per batch for estimating psi
psi.idx   <- NULL
for(B in sort(unique(sub.batch))) 
  psi.idx  <- c(psi.idx,sample(which(sub.batch==B),size=min(100,round(0.5*sum(sub.batch==B)))))
psi.batch  <- sub.batch[psi.idx]
nsub.psi   <- length(psi.idx)
Ysub.psi   <- Ysub[,psi.idx]

# initial estimate
lmu.hat<- Z <- log(Ysub+1) 
gmean  <- rowMeans(Z)
proj.Z <- projection(rep.sub,Zmat=Z)[,apply(Msub,1,which)]
RM_Z   <- Z - proj.Z
U0     <- irlba::irlba(RM_Z,nv=k)$v
# alpha(n x k) matrix
alpha  <- Z %*% U0

# reduce outliers
a.med <- matrixStats::colMedians(alpha)
a.mad <- matrixStats::colMads(alpha)
for(i in 1:ncol(alpha)) {
 alpha[which(alpha[,i]> (a.med[i]+4*a.mad[i])),i] <- a.med[i]+4*a.mad[i]
 alpha[which(alpha[,i]< (a.med[i]-4*a.mad[i])),i] <- a.med[i]-4*a.mad[i]
}
alpha.c <- as.matrix(alpha[ctl,])
Z.c    <- Z[ctl,] - gmean[ctl]

# initiate W (m x k)
#Wsub <- Matrix::crossprod(Z.c,alpha.c) %*% solve(Matrix::crossprod(alpha.c)) 
Wsub      <- matrix(0,nsub,k)
Wsub[,1]  <- Matrix::crossprod(Z.c,as.matrix(alpha.c[,1])) %*% solve(lambda.a + Matrix::crossprod(as.matrix(alpha.c[,1]),as.matrix(alpha.c[,1]))) 
if(k>1) {
 for(j in 2:k) {
   Z.c <- Z.c - outer(alpha.c[,j-1],Wsub[,j-1])
   Wsub[,j] <- Matrix::crossprod(Z.c,as.matrix(alpha.c[,j])) %*% solve(lambda.a + Matrix::crossprod(as.matrix(alpha.c[,j]),as.matrix(alpha.c[,j]))) 
 }
}
# reduce outliers
w.med <- matrixStats::colMedians(Wsub)
w.mad <- matrixStats::colMads(Wsub)
for(i in 1:ncol(Wsub)) {
  Wsub[which(Wsub[,i]> (w.med[i]+4*w.mad[i])),i] <- w.med[i]+4*w.mad[i]
  Wsub[which(Wsub[,i]< (w.med[i]-4*w.mad[i])),i] <- w.med[i]-4*w.mad[i]
}

# initiate Mb and gmean (intercept): by performing poisson GLM gene-by-gene
nW   <- ncol(Wsub)
nM   <- ncol(M)
ngene<- nrow(Y)
ns   <- ncol(Y)

bef   <- Sys.time()
# update Mb
Wa     <- Matrix::tcrossprod(alpha,Wsub)
Z.res  <- Z -  Wa  - gmean
Mb  <- tryCatch({ projection(rep.sub,Zmat=Z.res,Wmat=NULL,lambda=lambda.b)}, error = function(e) {NA})

# estimate initial psi
bef <- Sys.time() 
Wa   <- Matrix::tcrossprod(alpha,Wsub) 
offs <- Wa + Mb[,apply(Msub,1,which)] 
offs.psi <- offs[,psi.idx]
bef  <- Sys.time()
if(parallel) {
 if(nbatch==1) 
  psi <- estimateDisp.par(Ysub.psi,as.matrix(rep(1,nsub.psi)),offset=offs.psi,tagwise=TRUE,robust=TRUE,BPPARAM=BPPARAM)$tagwise.dispersion
 if(nbatch>1) {
  psi <- estimateDisp.par(Ysub.psi[,psi.batch==1],as.matrix(rep(1,sum(psi.batch==1))),
				offset=offs.psi[,psi.batch==1],tagwise=TRUE,robust=TRUE,BPPARAM=BPPARAM)$tagwise.dispersion
  for(B in 2:max(psi.batch)) {
      psi <- cbind(psi,estimateDisp.par(Ysub.psi[,psi.batch==B],as.matrix(rep(1,sum(psi.batch==B))),
				offset=offs.psi[,psi.batch==B],tagwise=TRUE,robust=TRUE,BPPARAM=BPPARAM)$tagwise.dispersion)
  }
 }
} 
if(!parallel) {
 if(nbatch==1) 
  psi <- edgeR::estimateDisp(Ysub.psi,as.matrix(rep(1,nsub.psi)),offset=offs.psi,tagwise=TRUE,robust=TRUE)$tagwise.dispersion
 if(nbatch>1) {
  psi <- edgeR::estimateDisp(Ysub.psi[,psi.batch==1],as.matrix(rep(1,sum(psi.batch==1))),
				offset=offs.psi[,psi.batch==1],tagwise=TRUE,robust=TRUE)$tagwise.dispersion
  for(B in 2:max(sub.batch)) 
   psi <- cbind(psi,edgeR::estimateDisp(Ysub.psi[,psi.batch==B],as.matrix(rep(1,sum(psi.batch==B))),
				offset=offs.psi[,psi.batch==B],tagwise=TRUE,robust=TRUE)$tagwise.dispersion)
 }
}
aft <- Sys.time()
#print(paste0('time to calculate psi:',  difftime(aft,bef,units='secs')))

# initiate
conv <- FALSE
iter.outer <- 0
logl.outer <- NULL
step <- rep(1,nsub)
print('Start...')
while(!conv) {
 iter.outer <- iter.outer + 1
 logl.beta <- NULL
 
 ### calculate initial loglik
 lmu.hat    <- gmean + Mb[,apply(Msub,1,which)] +  Matrix::tcrossprod(alpha,Wsub) 
 if(nbatch>1) {
  temp <- foreach(i=1:nsub, .combine=c) %dopar% {
    tmp <- dnbinom(Ysub[,i],mu=exp(lmu.hat[,i]),size=1/psi[,sub.batch[i]],log=TRUE)
    tmp[is.infinite(tmp) | is.na(tmp)] <- log(.Machine$double.xmin)
    return(sum(tmp))
  }
 } 
 if(nbatch==1) {
  temp <- foreach(i=1:nsub, .combine=c) %dopar% {
    tmp <- dnbinom(Ysub[,i],mu=exp(lmu.hat[,i]),size=1/psi,log=TRUE)
    tmp[is.infinite(tmp) | is.na(tmp)] <- log(.Machine$double.xmin)
    return(sum(tmp))
  }
 } 
 loglik <- temp

 #############################
 step <- rep(1,nsub)
 conv.beta <- FALSE
 # saving temporary best estimate
 best.W <- Wsub ; best.psi <- psi ; best.a <- alpha; best.Mb <- Mb; best.gmean <- gmean
 
 iter <- 1
 halving <- 0
 while(!conv.beta) {
 print(paste0('Inner iter:', iter))
 # calc working vector Z for a subset of cells 
 lmu.hat    <- gmean + Mb[,apply(Msub,1,which)] + Matrix::tcrossprod(alpha,Wsub) 

 # weight based on NB
 if(nbatch>1)
  sig.inv<- 1/(exp(-lmu.hat) + sweep(psi[,sub.batch],2,rep(1,nsub),'/'))
 if(nbatch==1)
  sig.inv<- 1/(exp(-lmu.hat) +outer(psi,rep(1,nsub),'/'))
 
 # calculate signed deviance
 if(nbatch>1 & robust) {
   signdev <- sign(Ysub-exp(lmu.hat))* sqrt(2*(dnbinom(Ysub,mu=Ysub,size=sweep(1/psi[,sub.batch],2,rep(1,nsub),'*'),log=TRUE) - 
						dnbinom(Ysub,mu=exp(lmu.hat),size=sweep(1/psi[,sub.batch],2,rep(1,nsub),'*'),log=TRUE) ))
 }
 if(nbatch==1 & robust) {
   signdev <- sign(Ysub-exp(lmu.hat))* sqrt(2*(dnbinom(Ysub,mu=Ysub,size=outer(1/psi,rep(1,nsub),'*'),log=TRUE) - 
						dnbinom(Ysub,mu=exp(lmu.hat),size=outer(1/psi,rep(1,nsub),'*'),log=TRUE) ))
 }  
 
 
 Z      <- lmu.hat + sweep( ((Ysub+0.01)/(exp(lmu.hat)+0.01) - 1),2,step,'*')

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

 # get alpha 
 bef=Sys.time()
 wt.cell  <- matrixStats::colMeans2(sig.inv)
 # prevent outlier large weight
 p98      <- quantile(wt.cell,probs=0.98)
 wt.cell[which(wt.cell>=p98)] <- p98
 #wt.cell  <- rep(1,ncol(sig.inv))
 b        <- RM_Z %*% diag(wt.cell) %*% RM_W
 alpha    <- b %*% Matrix::chol2inv(chol(Matrix::crossprod(RM_W*wt.cell,RM_W)+lambda.a*diag(nW)))
 aft <- Sys.time()
 #print(dim(alpha))
 #print(paste0('time to calculate alpha:',  difftime(aft,bef,units='secs')))

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

 #Wsub <- foreach(i=1:nsub, .combine=cbind, .packages="Matrix") %dopar% { 
 # out <- rep(0,k)
 # wtvec     <- sig.inv[ctl,i] 
 # wtvec     <- wtvec/mean(wtvec)
 # out[1]      <- Rfast::lmfit(y=Z[ctl,i]-gmean[ctl],x=alpha[ctl,1],w=wtvec)$be
 # if(k>1) {
 #   for(j in 2:k) {
 #     out[j] <- Rfast::lmfit(y=Z[ctl,i]-gmean[ctl]-matrixStats::rowSums2(as.matrix(alpha[ctl,1:(j-1)]) %*% as.matrix(out[1:(j-1)])),x=alpha[ctl,j],w=wtvec)$be
 #   } 
 # }
 # as.matrix(out)
 #}
 Z.c     <- Z[ctl,]-gmean[ctl]
 alpha.c <- alpha[ctl,,drop=FALSE]
 wt.ctl <- rowMeans(sig.inv[ctl,])

 Wsub[,1]  <- Matrix::crossprod(Z.c,as.matrix(alpha.c[,1]*wt.ctl)) %*% solve(lambda.a + Matrix::crossprod(as.matrix(alpha.c[,1]*wt.ctl),as.matrix(alpha.c[,1]))) 
 if(k>1) {
  for(j in 2:k) {
   Z.c <- Z.c - outer(alpha.c[,j-1],Wsub[,j-1])
   Wsub[,j] <- Matrix::crossprod(Z.c,as.matrix(alpha.c[,j]*wt.ctl)) %*% solve(lambda.a +Matrix::crossprod(as.matrix(alpha.c[,j]*wt.ctl),as.matrix(alpha.c[,j]))) 
  }
 }
	   
 aft <- Sys.time()
 #print(paste0('time to calculate W:',  difftime(aft,bef,units='secs')))
 Wsub <- as.matrix(Wsub)
 if(ncol(Wsub)!=k) 
  Wsub <- t(Wsub)

 colNorm.Wsub <- matrixStats::colSums2(Wsub^2)
 # normalize W to have unit length for each col
 #Wsub  <- sweep(Wsub,2,sqrt(colNorm.Wsub),'/')
 
 # orthogonalize W (optional)
 if(nW>1 & ortho.W)  
   Wsub <- matlib::GramSchmidt(Wsub)

 # if new W has missing/inf values, revert to previous estimate
 #print(paste0('Problematic W entries=',sum(is.na(Wsub) | is.infinite(Wsub))))
 if(any(is.na(Wsub) | is.infinite(Wsub) )) Wsub <- W.old
 
 #step3a: update gmean
 Z.res  <- Z -  Matrix::tcrossprod(alpha,Wsub)  - Mb[,apply(Msub,1,which)]
 bef <- Sys.time()
 gmean<- matrixStats::rowSums2(Z.res*sig.inv)/matrixStats::rowSums2(sig.inv)
 aft <- Sys.time()
 #print(paste0('time to calculate gmean:',  difftime(aft,bef,units='secs')))

 #step 3b: update Mb
 Z.res  <- Z -  Matrix::tcrossprod(alpha,Wsub)  - gmean
 Mb.old <- Mb
 bef <- Sys.time()
 Mb  <- tryCatch({ projection(rep.sub,Zmat=Z.res,Wmat=sig.inv,lambda=lambda.b)}, error = function(e) {NA})
 # set Mb=0 for control genes
 Mb[ctl,] <- 0
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
 if(nbatch>1) {
  temp <- foreach(i=1:nsub, .combine=c) %dopar% {
    tmp <- dnbinom(Ysub[,i],mu=exp(lmu.hat[,i]),size=1/psi[,sub.batch[i]],log=TRUE)
    tmp[is.infinite(tmp) | is.na(tmp)] <- log(.Machine$double.xmin)
    return(sum(tmp))
  }
 } 
 if(nbatch==1) {
  temp <- foreach(i=1:nsub, .combine=c) %dopar% {
    tmp <- dnbinom(Ysub[,i],mu=exp(lmu.hat[,i]),size=1/psi,log=TRUE)
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
    Wsub    <- best.W; alpha <- best.a ; Mb <- best.Mb ; psi <- best.psi ; gmean <- best.gmean 
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
   best.W <- Wsub ; best.psi <- psi ; best.a <- alpha; best.Mb <- Mb ; best.gmean <- gmean 
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
offs.psi <- offs[,psi.idx]
if(updt.psi) {
if(parallel) {
 if(nbatch==1) 
  psi.new <- tryCatch({ estimateDisp.par(Ysub.psi,as.matrix(rep(1,nsub.psi)),offset=offs.psi,
				tagwise=TRUE,robust=TRUE,BPPARAM=BPPARAM)$tagwise.dispersion}, error = function(e) {NA})
 if(nbatch>1) {
  psi.new <- psi
  for(B in 1:max(psi.batch)) 
   psi.new[,B]<-tryCatch({estimateDisp.par(Ysub.psi[,psi.batch==B],as.matrix(rep(1,sum(psi.batch==B))),
				offset=offs.psi[,psi.batch==B],tagwise=TRUE,robust=TRUE,BPPARAM=BPPARAM)$tagwise.dispersion},error= function(e) {NA})
 }
}

if(!parallel) {
 if(nbatch==1) 
  psi.new <- tryCatch({ edgeR::estimateDisp(Ysub.psi,as.matrix(rep(1,nsub.psi)),offset=offs.psi,tagwise=TRUE,robust=TRUE)$tagwise.dispersion},error = function(e) {NA})
 if(nbatch>1) {
  psi.new <- psi
  for(B in 1:max(psi.batch)) 
   psi.new[,B] <- tryCatch({ edgeR::estimateDisp(Ysub.psi[,psi.batch==B],as.matrix(rep(1,sum(psi.batch==B))),
				offset=offs.psi[,psi.batch==B],tagwise=TRUE,robust=TRUE)$tagwise.dispersion},
				error = function(e) {NA})
 }
}
# update psi
psi[!is.na(psi.new) & !is.infinite(psi.new)] <- psi.new[!is.na(psi.new) & !is.infinite(psi.new)]
}

# record outer logl
Wsub <- best.W
alpha <- best.a
Mb    <- best.Mb
gmean <- best.gmean
# new changes
best.psi <- psi

# recalculate logl
lmu.hat    <- gmean + Mb[,apply(Msub,1,which)] +  Matrix::tcrossprod(alpha,Wsub) 
if(nbatch>1) {
  temp <- foreach(i=1:nsub, .combine=c) %dopar% {
    tmp <- dnbinom(Ysub[,i],mu=exp(lmu.hat[,i]),size=1/psi[,sub.batch[i]],log=TRUE)
    tmp[is.infinite(tmp) | is.na(tmp)] <- log(.Machine$double.xmin)
    return(sum(tmp))
  }
} 
if(nbatch==1) {
  temp <- foreach(i=1:nsub, .combine=c) %dopar% {
    tmp <- dnbinom(Ysub[,i],mu=exp(lmu.hat[,i]),size=1/psi,log=TRUE)
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
 alpha.c <- alpha[ctl,,drop=FALSE]
 # estimate W for all samples (initial)
 print('Estimating W for all samples...')
 bef=Sys.time()
 W <- matrix(0,ns,nW)
 W[subsamples.org[subsubsamples.org],] <- Wsub[1:length(subsubsamples.org),]
 # block size variable
 block.size <- min(5000,ncol(Y))
 nb    <- ceiling(ncol(Y)/block.size)
 alpha.c <- alpha[ctl,,drop=FALSE]
 for(block in 1:nb) {
  start.idx <- (block-1)*block.size+1 ; end.idx <- min(ns,block*block.size)
  Ysub  <- as.matrix(Y[ctl,start.idx:end.idx])
  lmu.hat <- gmean[ctl] + Matrix::tcrossprod(alpha.c,W[start.idx:end.idx,,drop=FALSE])
  Z.c  <- lmu.hat - gmean[ctl] + (Ysub+0.01)/(exp(lmu.hat)+0.01) - 1
  W[start.idx:end.idx,1]<-Matrix::crossprod(Z.c,as.matrix(alpha.c[,1]*wt.ctl)) %*% solve(lambda.a + Matrix::crossprod(as.matrix(alpha.c[,1]*wt.ctl),as.matrix(alpha.c[,1]))) 
  if(k>1) {
   for(j in 2:k) {
    Z.c <- Z.c - outer(alpha.c[,j-1],W[start.idx:end.idx,j-1])
    W[start.idx:end.idx,j] <- Matrix::crossprod(Z.c,as.matrix(alpha.c[,j]*wt.ctl)) %*% solve(lambda.a +Matrix::crossprod(as.matrix(alpha.c[,j]*wt.ctl),as.matrix(alpha.c[,j]))) 
   }
  }
 }
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
 iter.W <- 1
 while(!conv.W) {
  W.old <- W
  for(block in 1:nb) {
   start.idx <- (block-1)*block.size+1 ; end.idx <- min(ns,block*block.size)
   Ysub  <- as.matrix(Y[ctl,start.idx:end.idx])
   lmu.hat <- gmean[ctl] + Matrix::tcrossprod(alpha.c,W[start.idx:end.idx,,drop=FALSE])
   Z.c  <- lmu.hat - gmean[ctl] + (Ysub+0.01)/(exp(lmu.hat)+0.01) - 1
   W[start.idx:end.idx,1]<-Matrix::crossprod(Z.c,as.matrix(alpha.c[,1]*wt.ctl)) %*% solve(lambda.a + Matrix::crossprod(as.matrix(alpha.c[,1]*wt.ctl),as.matrix(alpha.c[,1]))) 
   if(k>1) {
    for(j in 2:k) {
     Z.c <- Z.c - outer(alpha.c[,j-1],W[start.idx:end.idx,j-1])
     W[start.idx:end.idx,j] <- Matrix::crossprod(Z.c,as.matrix(alpha.c[,j]*wt.ctl)) %*% solve(lambda.a +Matrix::crossprod(as.matrix(alpha.c[,j]*wt.ctl),as.matrix(alpha.c[,j]))) 
    }
   }
  }

  if(ncol(W)!=k) 
   W <- t(W)

  # normalize W and alpha
  calnorm.W <- sqrt(matrixStats::colSums2(W^2))
  #W     <- sweep(W,2,calnorm.W,'/')
  #alpha <- sweep(alpha,2,calnorm.W,'*')

  # orthogonalize W (optional)
  if(nW>1 & ortho.W)  
    W <- matlib::GramSchmidt(W) 

  # fixed the sign of W
  #for(i in 1:ncol(W)) 
  # W[,i] <- W[,i] * sign(cor(W[subsamples.org[subsubsamples.org],i],Wsub[1:length(subsubsamples.org),i]))

  crit.W <- mean( (abs(W-W.old)/abs(W.old))^2)
  #print(round(crit.W,10))
  iter.W <- iter.W + 1
  conv.W <- crit.W< 1e-5 | iter.W>8
 }
 aft=Sys.time()
 #print(paste0('Time to estimate W for all samples:',difftime(aft,bef,units='secs')))
}
} # end of outer IRLS loop

# now calculate Mb
print('Estimating Mb....')
Mb.all <- matrix(0,nrow(Y),ncol(Y))
# for cells with annotation
idx.annot <- unlist(rep.ind)
Mb.all[,idx.annot] <- Mb[,apply(M[idx.annot,],1,which)]
# for other cells
if(length(idx.annot) < ncol(Y) ) {
 lmu  <- gmean + Matrix::tcrossprod(alpha,W) + Mb.all
 Z    <- Mb.all + ( (as.matrix(Y)+0.01)/(exp(lmu)+0.01) - 1)
 if(nbatch>1)
  wtmat<- 1/(exp(-lmu) + psi[,batch])
 if(nbatch==1)
  wtmat<- 1/(exp(-lmu) + psi)
 Mb.all[,-idx.annot] <- Z[,-idx.annot] * (wtmat[,-idx.annot]/(wtmat[,-idx.annot]+lambda.b))
 Mb.all[ctl,] <- 0
}
# calibrate alpha and W
#for(k in 1:nW) {
# ab.W     <- Rfast::lmfit(y=Wsub[1:length(subsubsamples.org),k],x=cbind(1,W[subsamples.org[subsubsamples.org],k]))$be
# ab.alpha <- Rfast::lmfit(y=best.a[,k],x=cbind(1,alpha[,k]))$be
# W[,k]    <- ab.W[1] + ab.W[2]*W[,k]
# alpha[,k]<- ab.alpha[1] + ab.alpha[2]*alpha[,k]
#}
#alpha <- sweep(alpha,2,scfac.alpha,'/')
cor.check <- diag(as.matrix(cor(Wsub[1:length(subsubsamples.org),],W[subsamples.org[subsubsamples.org],])))

if(use.pseudosample) { 
   psi <- psi[,1:nbatch.org]
}
#output
return( list("counts"=Y,"W"=W, "M"=M, "ctl"=ctl, "logl"=logl.outer, "a"=alpha,"Mb"=Mb.all, "gmean"=gmean,
		"psi"=psi,'L.a'=lambda.a,'L.b'=lambda.b,batch=batch,corW=cor.check,subsamples=subsamples.org,Wsub=Wsub,Mb.sub=Mb) )
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



