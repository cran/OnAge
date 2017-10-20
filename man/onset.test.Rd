\name{onset.test}
\alias{onset.test}

\title{Tests the null hypothesis that the age at the onset of senescence is
  the same in two groups}

\description{The function takes as input a log-likelihood function
  defining a statistical model of the senescence and two datasets. It
  performs a likelihood ratio test of the null hypothesis that the age
  at the onset of senescence is the same in the two groups.}

\usage{
onset.test(ll, data1, data2, search.range, CI.lvl=0.95, 
    tol=.Machine$double.eps^0.25, warn=FALSE, do.plot=FALSE, 
    plot.file=NULL, grid.len=100)
}

\arguments{
  \item{ll}{A function computing the log-likelihood of data under a
    statistical model of the senescence.}
  \item{data1}{A \link[base]{data.frame} object containing the data for the first
    group.}
  \item{data2}{A \link[base]{data.frame} object containing the data for the second
    group.}
  \item{search.range}{The range over which the log-likelihood should be
    maximized with respect to the age at the onset of senescence.}
  \item{CI.lvl}{An optional confidence level to obtain a confidence
    interval for each onset estimate. No confidence interval is computed
    if \code{NA}. Default \code{0.95}.}
  \item{tol}{A tolerance, to be passed to the optimize function. Optional.}
  \item{warn}{A logical variable indicating whether the function should
    issue a warning when the computed likelihood ratio is negative (which
    may happen for numerical reasons). Optional.}
  \item{do.plot}{A logical variable indicating whether the function should
    produce a plot of the log-likelihood against age at the onset of senescence (with
    vertical dotted lines indicating the confidence interval if \code{CI.lvl} is not
    \code{NA}). Optional (default \code{FALSE}).}
  \item{plot.file}{If not \code{NULL} and \code{do.plot} is \code{TRUE}, the function makes the
    plot in a pdf file with name \code{plot.file}. Optional (default \code{NULL}).}
  \item{grid.len}{Integer value. If \code{do.plot} is \code{TRUE}, specifies the number
    of onset values for which the log-likelihood is evaluated between
    \code{search.range[1]} and \code{search.range[2]}. Optional (default \code{100}).}
}

\value{
  A \code{\link[base]{list}} with elements:
  \item{pv}{A p-value for the tested null hypothesis.}
  \item{est.1}{A maximum likelihood estimate for the age at the onset of
    senescence in the first group.}
  \item{est.2}{A maximum likelihood estimate for the age at the onset of
    senescence in the second group.}
  \item{est.joint}{A maximum likelihood estimate for the age at the
    onset of senescence in the merged groups.}
  \item{CI.1}{A confidence interval for the age at the onset of
    senescence in the first group, with confidence level \code{CI.lvl} 
    (\code{NA} if \code{CI.lvl} is \code{NA}).}
  \item{CI.2}{A confidence interval for the age at the onset of
    senescence in the second group, with confidence level \code{CI.lvl} 
    (\code{NA} if \code{CI.lvl} is \code{NA}).}
  \item{joint.CI}{A confidence interval for the age at the onset of
    senescence in the total group, with confidence level \code{CI.lvl} 
    (\code{NA} if \code{CI.lvl} is \code{NA}).}
  \item{lh0}{The log-likelihood maximized under the null hypothesis.}
  \item{lh1}{The log-likelihood maximized under the alternative hypothesis.}
  \item{llr}{The likelihood ratio statistic.}  
  \item{cvg.ok}{A logical variable indicating whether the computed
    likelihood ratio was negative (the returned value is \code{0} in this case).}
}

\examples{
if(requireNamespace("lme4", quietly=TRUE)) {
    data(RoeDeerMassData)
    RoeDeerMassData$ID <- factor(RoeDeerMassData$ID)
    RoeDeerMassData$cohort <- factor(RoeDeerMassData$cohort)

    dataFCH <- RoeDeerMassData[RoeDeerMassData$sex\%in\%"F"&
                                 RoeDeerMassData$population\%in\%"CH", ]
    dataMCH <- RoeDeerMassData[RoeDeerMassData$sex\%in\%"M"&
                                 RoeDeerMassData$population\%in\%"CH", ]
    dataFTF <- RoeDeerMassData[RoeDeerMassData$sex\%in\%"F"&
                                 RoeDeerMassData$population\%in\%"TF", ]
    dataMTF <- RoeDeerMassData[RoeDeerMassData$sex\%in\%"M"&
                                 RoeDeerMassData$population\%in\%"TF", ]

    ## b1: function for piecewise regression (transforms x into 0 before bp)
    b1 <- function(x, bp) ifelse(x < bp, 0, x - bp)

    ## Use this function to define the model in which the differential
    ## onset hypothesis is tested.
    ll <- function(thr, dataIn){
        logLik(lme4::lmer(body.mass ~ b1(age, thr) + age.at.last.capture + 
            last.year.of.capture + (1|ID) + (1|cohort), data=dataIn, 
            REML="FALSE"))
    }

    search.range <- c(6,12) # data not available before 6 years old
    search.range.TF <- search.range.CH <- search.range

    \dontrun{
      # if you don't run the example in RStudio, create a new graphic 
      # window with suitable dimensions
      if(.Platform$GUI!="RStudio") {dev.new(width=15, height=5)}
      # all graphs are square
      par(pty="s")
      
      # test for Trois Fontaines
      res.tf <- onset.test(ll, dataFTF, dataMTF, search.range.TF, do.plot=TRUE)
      
      if(.Platform$GUI!="RStudio") {dev.new(width=15, height=5)}
      par(pty="s")
      
      # test for Chizé
      res.ch <- onset.test(ll, dataFCH, dataMCH, search.range.CH, do.plot=TRUE)
      
      cat(sprintf("p-value for differential age at onset is \%g in 
          Trois Fontaines, \%g in Chizé", res.tf$pv, res.ch$pv))
    }
    
    # a shorter executable example for Trois Fontaines
    set.seed(40)
    shortdataFTF <- dataFTF[sample(1:NROW(dataFTF), 200), ]
    shortdataMTF <- dataMTF[sample(1:NROW(dataMTF), 100), ]
    res.tf.short <- onset.test(ll, shortdataFTF, shortdataMTF, search.range.TF, 
      do.plot=FALSE, CI.lvl=NA)
  }
}
