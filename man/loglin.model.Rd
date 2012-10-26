\name{loglin.model}
\alias{loglin.model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{define a log-linear model}
\description{
  Function to specify a hierarchical log-linear model. This is a particular case of a hmm model. 
}
\usage{loglin.model(lev, int = NULL, strata = 1, dismarg = 0, type = "b", 
D = TRUE, c.gen = TRUE, printflag = FALSE, names = NULL, formula = NULL)}
%- maybe also 'usage' for other objects documented here.
\arguments{\item{lev}{Vector of number of categories of variables}
  \item{int}{Generating class of the log-linear model (must be a list) or list of all the interactions included}
  \item{strata}{Number of strata}
  \item{dismarg}{List of interactions constrained by inequalities - see `hmmm.model'}
  \item{type}{"b" for baseline logits, "l" for local logits}
  \item{D}{Input argument for inequalities - see `hmmm.model'}
  \item{c.gen}{If FALSE the input \code{int} must be the list of the minimal interaction sets to be excluded}
\item{printflag}{If TRUE information on the included and excluded interactions are given}
\item{names}{A character vector whose elements are the names of the variables}
\item{formula}{A formula describing a log-linear model}}
\details{This function simplifies `hmmm.model' in the case of log-linear models. If \code{formula} is employed, \code{c.gen} and \code{int} 
must not be declared and \code{printflag} is not used.}
\value{An object of the class \code{hmmmmod} defining a log-linear model that can be estimated by `hmmm.mlfit'.}

\references{Agresti A (2012) Categorical data Analysis, (3ed), Wiley, New York.

Bergsma W, Croon M, Hagenaars JA (2009) Marginal Models for Dependent, Clustered, and Longitudinal Categorical Data. Springer.
}
\note{If \code{int} and \code{formula} are not supplied a saturated log-linear model is defined. For log-linear models where the parameters 
depend on covariates first define a saturated log-linear model and then use the function `create.XMAT'.
}
\seealso{
\code{\link{hmmm.model}}, \code{\link{hmmm.mlfit}}, \code{\link{create.XMAT}}
}
\examples{data(madsen)
y<-getnames(madsen)
names<-c("Infl","Sat","Co","Ho")

f<-~Co*Ho+Sat*Co+Infl*Co+Sat*Ho+Infl*Sat
model<-loglin.model(lev=c(3,3,2,4),formula=f,names=names)

# alternatively 
# model<-loglin.model(lev=c(3,3,2,4),
# int=list(c(3,4),c(2,3),c(1,3),c(2,4),c(1,2)),names=names)

mod<-hmmm.mlfit(y,model,maxit=3000)
print(mod,printflag=TRUE)
}
\keyword{models}
