\name{predict.mind}
\alias{predict.mind}
\title{Predict method for Multivariate Linear Mixed Model}
\description{
Predicted values based on Multivariate Linear Mixed Model object
}
\usage{
\method{predict}{mind}(object,data,type= "proj",dir_s=NULL,dir_cov=NULL,\dots)
}
\arguments{
\item{object}{an object of class "\code{mind}"}
\item{data}{an object in which to look with which to predict.}
\item{type}{the type of prediction required. ; the predictors type are the \code{"eblup"}, \code{"proj"} and \code{"synth"}.The default to \code{"proj"}.}
\item{dir_s}{optionally, if \code{type} is equal to \code{eblup}  a data frame with the count for each modalities of the responce variable for the covariate patterns must be provided.}
\item{dir_cov}{optionally, if \code{type} is equal to \code{eblup}  a data frame with the sample units count for the covariate patterns must be provided.}
\item{\dots}{arguments based from or to other methods}
}
\details{
\code{predict.mind} produces predicted values, obtained by means the regression parameters on the frame \code{data}.
In the actual version of \code{predict.mind} only unit level predictions are provided. If the type is equal to \code{proj} or  \code{synth} a data.frame with individual predictions for all the modalities of the responce variable will be produced. When the \code{eblup} predictor is chosen two more input data.frame must be provided and the function will produce predictions for all the profiles identified by the cross-classification of the domains and covariate patterns.
}

\value{
When the predictor type is set equal to \code{proj} or \code{synth} the \code{predict.mind} produces a \code{data.frame} of predictions with the columns name equal to the multivariate \code{responses}. If the predictor type is \code{eblup} a data frame with predictions for all the profile (cross-clafficiation of domain of interest and coavariate patterns) is provided.
}
\author{
Developed by Andrea Fasulo
}

\examples{
# Load example data
data(data_s);data(univ)

# The sample units cover 104 over 333 domains in the population data frame
length(unique(data_s$dom));length(unique(univ$dom))
\donttest{

# One random effect at domain level
formula<-as.formula(cbind(emp,unemp,inact)~(1|mun)+
factor(sexage)+factor(edu)+factor(fore))

# Drop from the universe data frame variables not referenced in the formula or in the broadarea
univ_1<-univ[,-6]

example.1<-mind.unit(formula=formula,dom="dom",data=data_s,universe=univ_1)

## Example 1
#Projection predictions
example.1.predict<-predict.mind(object=example.1,data=univ,type= "proj")

# Check if the sum of the unit level predictions at 
# domain level are equal to the mind.unit Projection predictions
ck<-cbind(univ,example.1.predict)
ck<-aggregate(cbind(emp,unemp,inact)~dom,ck,sum)
head(ck);head(example.1$PROJ)

## Example 2
#Synthetic predictions
example.1.synth<-predict.mind(object=example.1,data=univ,type="synth")
  
# Check if the sum of the unit level predictions at 
# domain level are equal to the mind.unit Synthetic predictions
ck<-cbind(univ,example.1.synth)
ck<-aggregate(cbind(emp,unemp,inact)~dom,ck,sum)
head(ck);head(example.1$SYNTH)

## Example 3
#EBLUP predictions
inp_1<-aggregate(cbind(emp,unemp,inact)~dom+mun+sexage + edu + fore,data_s,sum)
  
inp_2<-aggregate(emp+unemp+inact~dom+mun+sexage+edu +fore,data_s,sum)
  
example.1.eblup<-predict.mind(object=example.1,data=univ_1,type="eblup",dir_s=inp_1,dir_cov=inp_2)
  
# Check if the sum of the predictions at 
# profile level are equal to the mind.unit Eblup predictions
ck<-aggregate(cbind(emp,unemp,inact)~dom,example.1.eblup,sum)
head(ck);head(example.1$EBLUP)

}
}