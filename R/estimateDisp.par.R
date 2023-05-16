#' Parallelized version of \code{edgeR::estimateDisp} function. 
#'
estimateDisp.par <-
function (y, design = NULL, group = NULL, lib.size = NULL, offset = NULL, 
    prior.df = NULL, trend.method = "locfit",  
    tagwise = TRUE, span = NULL, min.row.sum = 5, grid.length = 21, 
    grid.range = c(-10, 10), robust = FALSE, winsor.tail.p = c(0.05, 
        0.1), tol = 1e-06, weights = NULL, BPPARAM=BPPARAM, ...) 
{
    require(foreach,quietly=TRUE)
    y <- as.matrix(y)
    ntags <- nrow(y)
    if (ntags == 0) 
        return(list(span = span, prior.df = prior.df, prior.n = prior.n))
    nlibs <- ncol(y)
    trend.method <- match.arg(trend.method, c("none", "loess", 
        "locfit", "movingave"))
    if (is.null(group)) 
        group <- rep(1, nlibs)
    if (length(group) != nlibs) 
        stop("Incorrect length of group.")
    group <- edgeR::dropEmptyLevels(group)
    if (is.null(lib.size)) 
        lib.size <- colSums(y)
    if (length(lib.size) != nlibs) 
        stop("Incorrect length of lib.size.")
    offset <- edgeR:::.compressOffsets(y, lib.size = lib.size, offset = offset)
    weights <- edgeR:::.compressWeights(y, weights)
    sel <- rowSums(y) >= min.row.sum
    sely <- edgeR:::.subsetMatrixWithoutCopying(y, i = sel)
    seloffset <- edgeR:::.subsetMatrixWithoutCopying(offset, i = sel)
    selweights <- edgeR:::.subsetMatrixWithoutCopying(weights, i = sel)
    spline.pts <- seq(from = grid.range[1], to = grid.range[2], 
        length = grid.length)
    spline.disp <- 0.1 * 2^spline.pts
    grid.vals <- spline.disp/(1 + spline.disp)
    l0 <- matrix(0, sum(sel), grid.length)
    if (is.null(design)) {
        cat("Design matrix not provided. Switch to the classic mode.\n")
        if (length(levels(group)) == 1) 
            design <- matrix(1, nlibs, 1)
        else design <- model.matrix(~group)
        if (all(tabulate(group) <= 1)) {
            warning("There is no replication, setting dispersion to NA.")
            return(list(common.dispersion = NA, trended.dispersion = NA, 
                tagwise.dispersion = NA))
        }
        eq <- equalizeLibSizes(y, group = group, dispersion = 0.01, 
            lib.size = lib.size)
        y.pseudo <- eq$pseudo.counts[sel, , drop = FALSE]
        y.split <- splitIntoGroups(y.pseudo, group = group)
        delta <- optimize(commonCondLogLikDerDelta, interval = c(1e-04, 
            100/(100 + 1)), tol = tol, maximum = TRUE, y = y.split, 
            der = 0)
        delta <- delta$maximum
        disp <- delta/(1 - delta)
        eq <- equalizeLibSizes(y, group = group, dispersion = disp, 
            lib.size = lib.size)
        y.pseudo <- eq$pseudo.counts[sel, , drop = FALSE]
        y.split <- splitIntoGroups(y.pseudo, group = group)
        foreach(j = 1:grid.length,  .combine=cbind, .packages="edgeR") %:% 
  		 foreach (i = 1:length(y.split), .combine=c, .packages="edgeR") %dopar% {  
		  out <- condLogLikDerDelta(y.split[[i]], grid.vals[j],der = 0) + l0[, j]
		  return(out)
		 }
		
    }
	
    else {
        design <- as.matrix(design)
        if (ncol(design) >= nlibs) {
            warning("No residual df: setting dispersion to NA")
            return(list(common.dispersion = NA, trended.dispersion = NA, 
                tagwise.dispersion = NA))
        }
        glmfit <- edgeR::glmFit(sely, design, offset = seloffset, dispersion = 0.05, 
            prior.count = 0)
        zerofit <- (glmfit$counts < 1e-04 & glmfit$fitted.values < 
            1e-04)
        by.group <- edgeR:::.comboGroups(zerofit)
		#print(dim(by.group))
		#print(length(by.group))
	
        res.out <- foreach(j=1:length(by.group), .combine=rbind, .packages="edgeR") %dopar% {
	    subg      <- by.group[[j]]
            cur.nzero <- !zerofit[subg[1], ]
            if (!any(cur.nzero)) {
                next
            }
            if (all(cur.nzero)) {
                redesign <- design
            }
            else {
                redesign <- design[cur.nzero, , drop = FALSE]
                QR <- qr(redesign)
                redesign <- redesign[, QR$pivot[1:QR$rank], drop = FALSE]
                if (nrow(redesign) == ncol(redesign)) {
                  next
                }
            }
            cury <- edgeR:::.subsetMatrixWithoutCopying(sely, i = subg, 
                j = cur.nzero)
            curo <- edgeR:::.subsetMatrixWithoutCopying(seloffset, i = subg, 
                j = cur.nzero)
            curw <- edgeR:::.subsetMatrixWithoutCopying(selweights, i = subg, 
                j = cur.nzero)
            last.beta <- NULL
            out <- BiocParallel::bplapply(spline.disp,FUN=adjustedProfileLik,y = cury, 
                           design = redesign, offset = curo, weights = curw, 
                           start = last.beta, get.coef = FALSE,BPPARAM=BPPARAM)
   	    return(out)
            
       }
    }   
	
    # reorder properly parallel results
    if(length(by.group)==1) {
     for(j in 1:length(by.group)) 
        l0[by.group[[j]],] <- simplify2array(res.out)    
    }

    if(length(by.group)>1) {
     for(j in 1:length(by.group)) {
      for(i in 1:length(spline.disp))
        l0[by.group[[j]],i] <- c(res.out[j,][[i]])    
     }
    }

	
    overall <- edgeR::maximizeInterpolant(spline.pts, matrix(colSums(l0), 
        nrow = 1))
    common.dispersion <- 0.1 * 2^overall
    if (trend.method != "none") {
        AveLogCPM <- edgeR::aveLogCPM(y, lib.size = lib.size, offset=offset, dispersion = common.dispersion, 
            weights = weights)
        out.1 <- edgeR::WLEB(theta = spline.pts, loglik = l0, covariate = AveLogCPM[sel], 
            trend.method = trend.method,  
            span = span, overall = FALSE, individual = FALSE, 
            m0.out = TRUE)
        span <- out.1$span
        m0 <- out.1$shared.loglik
        disp.trend <- 0.1 * 2^out.1$trend
        trended.dispersion <- rep(disp.trend[which.min(AveLogCPM[sel])], 
            ntags)
        trended.dispersion[sel] <- disp.trend
    }
    else {
        AveLogCPM <- NULL
        m0 <- matrix(colMeans(l0), ntags, length(spline.pts), 
            byrow = TRUE)
        disp.trend <- common.dispersion
        trended.dispersion <- NULL
    }
    if (!tagwise) 
        return(list(common.dispersion = common.dispersion, trended.dispersion = trended.dispersion))
    if (is.null(prior.df)) {
        glmfit <- edgeR::glmFit(sely, offset = seloffset, weight = selweights, 
            design = design, dispersion = disp.trend, prior.count = 0)
        df.residual <- glmfit$df.residual
        zerofit <- (glmfit$counts < 1e-04 & glmfit$fitted.values < 
            1e-04)
        df.residual <- edgeR:::.residDF(zerofit, design)
        s2 <- glmfit$deviance/df.residual
        s2[df.residual == 0] <- 0
        s2 <- pmax(s2, 0)
        s2.fit <- limma::squeezeVar(s2, df = df.residual, covariate = AveLogCPM[sel], 
            robust = robust, winsor.tail.p = winsor.tail.p)
        prior.df <- s2.fit$df.prior
    }
    ncoefs <- ncol(design)
    prior.n <- prior.df/(nlibs - ncoefs)
    if (trend.method != "none") 
        tagwise.dispersion <- trended.dispersion
    else tagwise.dispersion <- rep(common.dispersion, ntags)
    too.large <- prior.n > 1e+06
    if (!all(too.large)) {
        temp.n <- prior.n
        if (any(too.large)) {
            temp.n[too.large] <- 1e+06
        }
        out.2 <- edgeR::WLEB(theta = spline.pts, loglik = l0, prior.n = temp.n, 
            covariate = AveLogCPM[sel], trend.method = trend.method, 
            span = span, overall = FALSE, 
            trend = FALSE, m0 = m0)
        if (!robust) {
            tagwise.dispersion[sel] <- 0.1 * 2^out.2$individual
        }
        else {
            tagwise.dispersion[sel][!too.large] <- 0.1 * 2^out.2$individual[!too.large]
        }
    }
    if (robust) {
        temp.df <- prior.df
        temp.n <- prior.n
        prior.df <- prior.n <- rep(Inf, ntags)
        prior.df[sel] <- temp.df
        prior.n[sel] <- temp.n
    }
    list(common.dispersion = common.dispersion, trended.dispersion = trended.dispersion, 
        tagwise.dispersion = tagwise.dispersion, span = span, 
        prior.df = prior.df, prior.n = prior.n)
}
