\name{predict.mind}
\alias{predict.mind}
\title{Predict method for Multivariate Linear Mixed Model}
\description{
Predicted values based on Multivariate Linear Mixed Model object
}
\usage{
## S3 method for class 'mind'
\method{predict}{mind}(object,data,\dots)
}
\arguments{
\item{object}{an object of class "\code{mind}"}
\item{data}{an object in which to look with which to predict.}
\item{\dots}{arguments based from or to other methods}
}
\details{
\code{predict.mind} produces predicted values, obtained by means the regression parameters on the frame \code{data}.
In the actual version of \code{predict.mind} (i) only unit level prediction (ii) using the Projection predictor are provided.

}

\value{
\code{predict.mind} produces a \code{data.frame} of predictions setting the columns name equal to the multivariate \code{responses}.
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
## Example 1
# One random effect at domain level
formula<-as.formula(cbind(emp,unemp,inact)~(1|mun)+
factor(sexage)+factor(edu)+factor(fore))

# Drop from the universe data frame variables not referenced in the formula or in the broadarea
univ_1<-univ[,-6]

example.1<-mind.unit(formula=formula,dom="dom",data=data_s,universe=univ_1)

#Predictions
example.1predict<-predict.mind(object=example.1,data=univ)

# Check if the sum of the unit level predictions at 
# domain level are equal to the mind.unit Projection predictions
ck<-cbind(univ,example.1predict)
ck<-aggregate(cbind(emp,unemp,inact)~dom,ck,sum)
head(ck);head(example.1$PROJ)

}
}