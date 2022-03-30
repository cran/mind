\name{mind.unit}
\alias{mind.unit}
\title{Fitting Unit level Multivariate Linear Mixed Model}
\description{
\code{mind.unit} is used to fit unit level multivariate linear mixed models [D'Alo', Falorsi 2021, FAO 2021]. It can be used to carry out different estimators (EBLUP, Synthetic and Projection) and the Mean Squared Error (MSE) for unplanned domain, analysis of the random effect and study of the variance components.
}
\usage{
mind.unit(formula,dom,data,universe,weights=NA,broadarea=NA,
max_iter=200,max_diff=1e-05,phi_u0=0.05,REML=TRUE)
}
\arguments{
\item{formula}{an object of class "formula": a symbolic description of the model to be fitted. The details of model specification are given under 'Details'.}
 \item{dom}{numeric, the domain of interest.  See also 'Details'.}
 \item{data}{a data frame containing the variables in the model, e.g. \code{\link{data_s}}.}
 \item{universe}{a data frame containing the complete list of the units belonging to the target population, along with the corresponding values of the auxiliary variables. Also an aggregated version of the universe information is possible, e.g. \code{\link{univ}}. See also 'Details'.}
  \item{weights}{an optional column of weights to be used in the fitting process. Should be NULL or a numeric vector. If non-NULL, weighted least squares is used with weights; otherwise ordinary least squares is used. See also 'Details'.}
\item{broadarea}{an optional character to be used if a broadarea is required in the model. See also 'Details'.}
\item{max_iter}{integer scalar. Numeber of maximum iteration for the optimization of the REML criterion (default=200) .}
\item{max_diff}{double number. Stopping criteria to be satisfied to achieve the REML convergence (defualt=1e-05).}
\item{phi_u0}{double number. Initialization value for the ratio among the variance components effect and the variance of the errors [Saei, Chambers 2003] (defualt=0.05)}
\item{REML}{logical scalar. Should the estimates be chosen to optimize the REML criterion (as opposed to the maximum-likelihood)?}
}
\details{
A typical predictor for a Multivariate Linear Mixed Model has the form \code{responses ~ random.terms+fixed.terms} where \code{responses} is the multivariate response, \code{random.terms} is a series of terms which specifies random intercept and \code{fixed.terms} is a series of terms which specified a linear predictor for \code{responses}.
\cr The \code{responses} can be specified as a column (so the \code{responses} have m different values as the modalities are) or as a m-column with the columns giving the presence and absence of the modalities (using \code{cbind} function).
\cr The \code{random.terms} in the formula will be re-ordered when both domain and marginal effect are presents so that domain effects come first, followed by the marginal. The \code{random.terms} must be \code{numeric} varaibles.
\cr In the actual version of \code{mind} (i) only qualitative \code{fixed.terms} are allowed.

The mandatory argument \code{dom} must be numeric and must not contain any missing value (NA).

The mandatory argument \code{universe} is a \code{data.frame} conteining the auxiliary information referenced in the \code{formula} for each unit of the population of interest.
\cr For computational reason it is possible use an aggregated version of the population information using the profile derived by the \code{random.terms} and \code{fixed.terms}. In this case a column equal to the summation of the population units for each profile is required.
\cr See \link{univ} for more details.

Non\code{-NULL weights} can be used to indicate that different observations have different variances. If no \code{weights} are specified all the units have an unitary weight. If specified must be present in \code{data}.

\code{broadarea} represents the grupping factor specifying the partitioning of the data. If non\code{-NULL broadarea} is includes different \code{mind.unit} fits should be performed according to \code{broadarea}. Must be present in \code{universe}.


}

\value{
\code{mind.unit} returns an object of class "\code{mind}". The object is a list contains 13 objects:
\item{EBLUP}{a data frame containing for the domain of interest the EBLUP estimates [Rao, Molina 2015] for the m-modalities  of the responce variable.}
\item{PROJ}{a data frame containing for the domain of interest the Projection estimates [Kim, Rao 2011] for the m-modalities of the responce variables.}
\item{SYNTH}{a data frame containing for the domain of interest the Synthetic estimates [Rao, Molina 2015] for the m-modalities of the responce variable.}
\item{mse_EBLUP}{a data frame containing for the domain of interest the MSE, along with the single components G1, G2 and G3, for the EBLUP estimator for the m-modalities of the variables of interest.}
\item{cv_EBLUP}{a data frame containing the coefficent of variation for the EBLUP estimator for the m-modalities of the variable of interest.}
\item{Nd}{a data frame with the total population of the domain of interest.}
\item{nd}{a data frame with the sample size of the sampled domain of interest.}
\item{r_effect}{a list containing the random effects for each modes of the \code{responses} and for each \code{broadarea} (if any).}
\item{beta}{a data frame with named columns of coefficents.}
\item{mod_performance}{a list containing fit indices, absolute error metrics, tests of overall model significance (taking into account only the \code{fixed.terms}) for each modes of the \code{responses} and for each \code{broadarea} (if any).}
\item{sigma_e}{a data frame with the residuals standard deviation \eqn{\sigma_e} for each modes of the \code{responses} and for each \code{broadarea} (if any).}
\item{sigma_u}{a data frame with the random effects standard deviation \eqn{\sigma_u} for each modes of the \code{responses}, for each random effect and for each \code{broadarea} (if any).}
\item{ICC}{a data frame with the Intraclass Coefficent Correlation for each modes of the \code{responses}, for each random effect and and for each \code{broadarea} (if any).
The population ICC in this framework is:
\deqn{ICC = \frac{\sigma_u^2}{(\sigma_u^2+\sigma_e^2)}}{%
        ICC = \sigma_u^2/(\sigma_u^2+\sigma_e^2)}
This ICC can be generalized to allow for covariate effects, in which case the ICC is interpreted as capturing the within-class similarity of the covariate-adjusted data values.}
}
\author{
Developed by Michele D'Alo', Stefano Falorsi,  Andrea Fasulo
}
\references{
Battese, G., E., Harter, R., M., Fuller, W., A., (1988). 'An Error-Components Model for Prediction of County Crop Areas Using Survey and
Satellite Data', Journal of the American Statistical Association Vol. 83, No. 401 (Mar., 1988), pp. 28-36.

Datta, G., S., Day, B., Basawa, I., (1999). 'Empirical best linear unbiased and empirical Bayes prediction in multivariate small area estimation', Journal of Statistical Planning and Inference, Volume 75, Issue 2, 1 January 1999, Pages 269-279


D'Alo', M., Falorsi, S., Fasulo, A., (2021). 'MIND, an R package for multivariate small area
estimation with multiple random effects', SAE2021 BIG4small. Book of short papers, Pages 43-48

ESSnet on SAE, (2012). 'Guidelines for the application of the small area estimation methods in NSI sample surveys'

FAO (2021). 'Guidelines on data disaggregation for SDG Indicators using survey data', pp. 105. http://www.fao.org/publications/card/en/c/CB3253EN/

Harmening, S., Kreutzmann, A.K., Pannier, S., Salvati, N., Schmid, T., (2021). 'A Framework for Producing Small Area Estimates Based on Area-Level Models in R, The R package emdi vignette'

Kim, J. K., Rao, J. N., (2011). 'Combining data from two independent surveys: a model-assisted approach', Biometrika 99(1), 85100.

Rao, J.N., Molina, I., (2015). 'Small Area Estimation', John Wiley & Sons

Saei, A., Chambers, R., (2003). 'Small Area Estimation Under Linear and Generalized Linear Mixed Models With Time and Area Effects', S3RI Methodology Working Paper M03/15
}

\seealso{
\link{predict.mind} has examples of fitting multivariate variables.
}
\examples{
# Load example data
data(data_s);data(univ)

# The sample units cover 104 over 333 domains in the population data frame
length(unique(data_s$dom));length(unique(univ$dom))
\donttest{
## Example 1
# One random effect at domain level
# Double possible formulations
formula<-as.formula(occ_stat~(1|mun)+
factor(sexage)+factor(edu)+factor(fore))
#or
formula<-as.formula(cbind(emp,unemp,inact)~(1|mun)+
factor(sexage)+factor(edu)+factor(fore))

# Drop from the universe data frame variables not referenced in the formula or in the broadarea
univ_1<-univ[,-6]

example.1<-mind.unit(formula=formula,dom="dom",data=data_s,universe=univ_1)
summary(example.1$EBLUP)
rm(univ_1)

## Example 2
# One random effect for a marginal domain
formula<-as.formula(occ_stat~(1|pro)+factor(sexage)+factor(edu)+factor(fore))

# Drop from the universe data frame variables not referenced in the formula or in the broadarea
univ_2<-univ[,-5]

example.2<-mind.unit(formula=formula,dom="dom",data=data_s,universe=univ_2)
summary(example.2$EBLUP)
rm(univ_2)

## Example 3
# Two random effects both at domain level and marginal level
formula<-as.formula(occ_stat~(1|mun)+(1|pro)+
factor(sexage)+factor(edu)+factor(fore))

example.3<-mind.unit(formula=formula,dom="dom",data=data_s,universe=univ)
summary(example.3$EBLUP)

## Example 4
# One random effect at domain level and with broadarea
formula<-as.formula(occ_stat~(1|mun)+factor(edu)+factor(fore))

# Drop from the universe data frame variables not referenced in the formula or in the broadarea
univ_4<-univ[,-2]

example.4<-mind.unit(formula=formula,dom="dom",data=data_s,universe=univ_4,broadarea="pro")
summary(example.4$EBLUP)
rm(univ_4)
}
}