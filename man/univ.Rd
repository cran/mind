\name{univ}
\Rdversion{1.0}
\alias{univ}
\docType{data}
\title{
Synthetic population dataset for Multivatiate Linear Mixed Model}
\description{
Synthetic population data frame containing the complete list of the units belonging to the target population along with the corresponding values of the auxiliary variables.
}
\usage{
data(univ)
}
\format{
   A data frame with 514320 observations on 8 variables:
  \describe{
    \item{\code{dom}}{domain of interest codes}
    \item{\code{sexage}}{cross classification of age and sex}
    \item{\code{edu}}{educational level}
    \item{\code{fore}}{bynary variable, 2 for foreigner 1 otherwise}
    \item{\code{mun}}{municipal codes}
    \item{\code{pro}}{provincial codes}
    \item{\code{tot}}{column of 1}
    }
  }
\details{
The informations on the population are the same collected in the syntethic sample \link{data_s} appart from the information on the occupational status that are present only for the sample units.

\code{\link{mind.unit}} allows to use a data frame of known population totals based on the marginal distribution of the profile identified by the auxiliary variables (See 'Examples').
}
\examples{

library(dplyr)

# Load example data
data(data_s);data(univ)
summary(univ)

\donttest{
formula<-as.formula(occ_stat~(1|pro)+factor(sexage)+factor(edu)+factor(fore))

# Drop from the universe data frame variables not referenced in the formula or in the broadarea
univ_1<-univ[,-5]

# 1) Estimation using the complete list of the unit beloging the target population:
example.1<-mind.unit(formula=formula,dom="dom",data=data_s,universe=univ_1)
rm(univ_1)

# Creation of the know population totals object:
univ_ag<-aggregate(tot~-1+factor(dom)+factor(pro)+
factor(sexage)+factor(edu)+factor(fore),univ,sum)
colnames(univ_ag)<-c("dom","pro","sexage","edu","fore","tot")

# Set all variables as numeric.
#Remember that only the domains codes and the random terms must to be numeric variables.
univ_ag <- mutate_all(univ_ag, function(x) as.numeric(as.character(x)))

# 2) Estimation using the know population totals (totals in univ_ag) :
example.2<-mind.unit(formula=formula,dom="dom",data=data_s,universe=univ_ag)
}
}
\keyword{datasets}
