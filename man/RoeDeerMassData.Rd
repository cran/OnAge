\name{RoeDeerMassData}
\alias{RoeDeerMassData}
\docType{data}

\title{Data on 454 roe deer}

\description{A dataset containing the age, body weight and other
  attributes of 454 roe deer from Chizé and Trois Fontaines.
}

\usage{data(RoeDeerMassData)}

\format{A data frame with 1428 rows and 16 variables:
  \describe{
    \item{ID}{An individual specific identifier.}
    \item{age}{The age of the individual (in years).}
    \item{sex}{The sex of the individual (F for females, M for males).}
    \item{population}{The population site of the individual: either Chizé (CH)
      or Trois Fontaines (TF).}
    \item{cohort}{The year of birth of the individual.}
    \item{body.mass}{The individual body mass measurement at this
      particular age (in kg).}
    \item{year.of.capture}{The year of capture of the individual for this particular
      body mass measurement.}
    \item{age.at.last.capture}{The age of the individual at its last capture (in years).}
    \item{last.year.of.capture}{A binary variable indicating whether or
      not the body mass measurement was done in the last capture of
      individual's life (no = 0; yes = 1).}
    \item{early.mass}{The predicted body mass during the first winter at
      the median date of capture (27th of January) (in kg).}
    \item{residual.early.mass}{The difference between the individual
      early mass and the average early mass of individuals within each sex
      and population (in kg).}
    \item{adult.body.mass}{Median body mass during the prime age stage
      between four (after body growth has ceased) and six years of age
      (when senescence in body mass can begin)}
    \item{mass.gain}{The body mass gain corresponding to the difference
      between adult body mass and early mass (in kg).}    
    \item{late.growth}{The body mass gain once the effect of early mass
      has been accounted for within each sex and population (i.e. residual
      of the linear regression model between early mas and mass gain) (in
      kg)}
    \item{early.mass.group}{The binary variable indicating whether the
      individual has an early mass less, or equal or greater than the median
      early mass within each sex and population (no = 0; yes =1).}  
    \item{late.growth.group}{The binary variable indicating whether the
      individual has a late growth less, or equal or greater than the
      median early mass within each sex and population (no = 0; yes = 1).}
  }
}

\source{\url{http://onlinelibrary.wiley.com/doi/10.1111/oik.04421/abstract}}

\examples{
data(RoeDeerMassData)
str(RoeDeerMassData)
}

\keyword{datasets}
