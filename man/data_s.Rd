\name{data_s}
\Rdversion{1.0}
\alias{data_s}
\docType{data}
\title{
Synthetic sample dataset for Multivatiate Linear Mixed Model
}
\description{
Synthetic data frame containing a sample of 10.000 individuals along with socio-economic indicators.
}
\usage{
data(data_s)
}
\format{
  A data frame with 10000 observations on 11 variables:
  \describe{
    \item{\code{dom}}{domain of interest codes, corresponding to the munipal codes}
    \item{\code{emp}}{binary variable, 1 for employed 0 otherwise}
    \item{\code{unemp}}{binary variable, 1 for unemployed 0 otherwise}              \item{\code{inact}}{binary variable, 1 for inactive 0 otherwise}
    \item{\code{sexage}}{cross classification of age and sex}
    \item{\code{edu}}{educational level}
    \item{\code{fore}}{binary variable, 2 for foreigner 1 otherwise}
    \item{\code{mun}}{municipal codes}
    \item{\code{pro}}{provincial codes}
    \item{\code{occ_stat}}{occupational status, 1 for employed 2 for unemployed 3 for inactive}

  }
 }
\details{
The informations on the sample unit are the same collected in the syntethic population dataframe \link{univ} appart from the information on the occupational status that are present only for the sample units.
}
\examples{
# Load example data
data(data_s)
summary(data_s)
}
\keyword{datasets}
