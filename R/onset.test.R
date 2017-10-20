## Build a likelihood ratio test of H0 (same onset in both groups) vs
## H1 (different onset allowed)
onset.test <- function(ll, data1, data2, search.range, CI.lvl=0.95, tol=.Machine$double.eps^0.25, warn=FALSE, do.plot=FALSE, plot.file=NULL, grid.len=100){
    ## Maximize the likelihood over all parameters including onset for
    ## group 1
    c1.opt <- optimize(ll, interval=search.range, dataIn = data1, maximum=TRUE, tol=tol)
    
    ## Maximize the likelihood over all parameters including onset for
    ## group 2
    c2.opt <- optimize(ll, interval=search.range, dataIn = data2, maximum=TRUE, tol=tol)

    ## Joint loglikelihood across data1 and data2 with a single onset (ie,
    ## loglikelihood under H0)
    ll.joint <- function(thr, ll, data1, data2){
        ll(thr, data1) + ll(thr, data2)
    }
    ## Maximize the likelihood over all parameters including (single)
    ## onset for both groups    
    joint.opt <- optimize(ll.joint, interval=search.range, ll, data1 = data1, data2 = data2, maximum=TRUE, tol=tol)
    
    ## Loglikelihood under H1 (sum of the loglikelihood obtained
    ## for both groups)    
    lh1 <- as.numeric(c1.opt$objective + c2.opt$objective)
    ## Loglikelihood under H0
    lh0 <- as.numeric(joint.opt$objective)

    ## Loglikelihood ratio (LLR) statistic
    llr <- 2*(lh1 - lh0)

    ## Sanity check: the likelihood under H1 should be larger than
    ## under H0 by construction.    
    cvg.ok <- TRUE
    if(llr < 0){
        ## lh1 <- lh0
        llr <- 0
        if(warn){
            warning('[onset.test] Likelihood under H1 was lower than under H0, probably due to optimization issues. Returning null log likelihood ratio. est.1 and est.2 values are suboptimal and should be both replaced by est.joint.')
        }
        cvg.ok <- FALSE
    }

    ## LLR is asymptotically chi2_1 under H0. Use this to compute an
    ## (asymptotic) p-value.    
    pv <- 1 - pchisq(llr, 1)
    
    ## Optionally get confidence intervals on the ages at onset
    if(!is.na(CI.lvl)){
      CI.1 <- CI.from.opt(c1.opt, ll, lvl=CI.lvl, search.range, data1)
      CI.2 <- CI.from.opt(c2.opt, ll, lvl=CI.lvl, search.range, data2)
      joint.CI <- CI.from.opt(joint.opt, ll.joint, lvl=CI.lvl, search.range, ll, data1, data2)
    }else{
        CI.1 <- CI.2 <- joint.CI <- NA
    }
    
    if(do.plot){

        ll.plot <- function(ll.fun, ..., x.grid, x.opt, ll.opt, CI, data.name, CEX=1.5){
            ll.grid <- sapply(x.grid, ll.fun, ...)
            plot(x.grid, ll.grid, main=data.name, xlab="Age at onset", ylab="Log likelihood", cex.lab=CEX, cex.axis=CEX)
            points(x.opt, ll.opt, pch=8, col='red')
            if(!any(is.na(CI))){
                abline(v=CI[1], lty=2)
                abline(v=CI[2], lty=2)
                abline(h=min(ll.fun(CI[1], ...), ll.fun(CI[2], ...)), lty=2)
            }
        }
        
        if(!is.null(plot.file)){
            pdf(file=plot.file, width=15, height=5)
        }
        mars <- c(c(4.5, 5, 2.5, 1))
        par(mar=mars, mfrow=c(1, 3))
        thr.grid <- seq(search.range[1], search.range[2], length.out=grid.len)
        ll.plot(ll, data1, x.grid=thr.grid, x.opt=c1.opt$maximum, ll.opt=c1.opt$objective, CI=CI.1, data.name='Group 1')
        ll.plot(ll, data2, x.grid=thr.grid, x.opt=c2.opt$maximum, ll.opt=c2.opt$objective, CI=CI.2, data.name='Group 2')
        ll.plot(ll.joint, ll, data1, data2, x.grid=thr.grid, x.opt=joint.opt$maximum, ll.opt=joint.opt$objective, CI=joint.CI, data.name='All data')
        
        if(!is.null(plot.file)){
            dev.off()
        }        
    }
    
    return(list(pv=pv,
                est.1=c1.opt$maximum, est.2=c2.opt$maximum, est.joint=joint.opt$maximum,
                CI.1=CI.1, CI.2=CI.2, joint.CI=joint.CI,
                lh0=lh0, lh1=lh1, llr=llr,
                cvg.ok=cvg.ok))
}


## Build CI by identiying smallest and largest values thr of onset for
## which the null hypothesis true_onset=thr is rejected at level
## 1-lvl. Caveat: assumes the log-likelihood is monotonically
## decreasing, i.e., that the log-likelihood ratio statistic is
## monotonically increasing around the maximum likelihood estimator of
## onset (can be checked visually using do.plot=TRUE in
## onset.test). Otherwise the CI may not be connected and a more
## sophisticated method should be used.

CI.from.opt <- function(opt, ll.fun, lvl=0.95, search.range, ...){
    
    ## Zero of this function is an onset thresholds thr for which
    ## H0:true_onset=thr is rejected at level exactly 1-lvl using a
    ## likelihood ratio test. The ML estimator of onset leads to
    ## rejection at level 1 by construction of the likelihood ratio
    ## statistic (equal to 0).  We assume (monotonicity of the
    ## likelihood) that all values between the ML estimator and this
    ## zero lead to rejection at level >= 1-lvl, making them part of
    ## the lvl confidence interval.    
    dist.from.quantile <- function(thr)
    {
        2*(opt$objective - ll.fun(thr, ...)) - qchisq(lvl, 1)
    }
    
    bp.opt <- opt$maximum # Maximum likelihood onset
    
    if(dist.from.quantile(search.range[1]) < 0){
        ## warning(sprintf('[CI.from.opt] H0: onset=search.range[1] not rejected at level %g. Setting lower bound of the CI to search.range[1].', 1-lvl))
        lb <- search.range[1]
    }else{
        lb <- uniroot(dist.from.quantile, lower=search.range[1], upper=bp.opt)$root
    }
    if(dist.from.quantile(search.range[2]) < 0){
        ## warning(sprintf('[CI.from.opt] H0: onset=search.range[2] not rejected at level %g. Setting upper bound of the CI to search.range[2].', 1-lvl))
        ub <- search.range[2]
    }else{
        ub <- uniroot(dist.from.quantile, lower=bp.opt, upper=search.range[2])$root
    }
    CI <- c(lb, ub)    
    names(CI) <- c('lb', 'ub')
    
    return(CI)
}
