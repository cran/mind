\name{benchSAE}
\alias{benchSAE}
\title{Benchmark for SAE}
\description{
Benchmarked values based on ratio adjustment of the SAE estimates}
\usage{
benchSAE(estim_sae,benchmark_area,area,name_dom,estimator,Nest)
}
\arguments{
\item{estim_sae}{a data frame containing the arguments \code{area}, \code{name_dom}, \code{estimator} and \code{Nest}}
\item{benchmark_area}{a data frame identifing the area to be benchmarked among with the benchmark value}
\item{area}{character, identified the name of the area to be benchmarked}
\item{name_dom}{character, identified the name of the domains}
\item{estimator}{vector, contain the name of the estimates to be benchmarked}
\item{Nest}{character, identified the name of the total population size for domain}
}

\details{
The \code{benchSAE} function allows (i) to benchmark more than one SAE estimates at times and (ii) to  specified the broadarea to be benchmarked, from the national level to a more disaggregated level. 
}

\value{
\code{benchSAE} produces a \code{data.frame} with a column of domain indicator and set of columns (as specified in \code{estimator}) of benchmarked point predictions.
}
\author{
Developed by Michele  D'Alò
}

\references{Datta,G. S., Ghosh, M., Steorts, R. and Maples, J. (2010) Bayesian benchmarking with applications to small area estimation. Test, 20, 574–588.}

\examples{
# Load example data

data(data_s);data(univ)
\donttest{
tot<-aggregate(tot~dom+pro,univ,sum)
  
# One random effect at domain level


  formula<-as.formula(cbind(emp,unemp,inact)~(1|mun)+
                        factor(sexage)+factor(edu)+factor(fore))
  
# Drop from the universe data frame variables not referenced in the formula or in the broadarea

  univ_1<-univ[,-6]

# Estimate mind model benchmarking the unemployment variable

  example.1<-mind.unit(formula=formula,dom="dom",data=data_s,universe=univ_1,MSE=FALSE)
  SAEest<-example.1$EBLUP[,c(1,3)]
  SAEest<-merge(SAEest,tot)
  SAEest$unemp<-SAEest$unemp/SAEest$tot
  bench<-data.frame(pro=c(10,11),bb=c(0.27,0.31))
  
# Benchmark the unemployment point estimates
  
  SAEest_bench<-benchSAE(estim_sae=SAEest,
                          benchmark_area=bench,
                          area="pro",
                          name_dom="dom",
                          estimator="unemp",
                          Nest="tot")
}
}