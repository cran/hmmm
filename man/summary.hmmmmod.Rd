\name{summary.hmmmmod}
\alias{summary.hmmmmod}
\alias{print.hmmmmod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{summary and print for the class hmmmmod}
\description{
The generic functions `summary' and `print' are adapted to the objects inheriting from class \code{hmmmod}
(summary.hmmmmod, print.hmmmmod) to display the summary of a model defined by `hmmm.model'.
}
\usage{
\method{summary}{hmmmmod}(object,...)
\method{print}{hmmmmod}(x,...)}
\arguments{
  \item{object,x}{An object of the class \code{hmmmod}, i.e. a result of `hmmm.model'}
\item{\dots}{Further arguments passed to or from other methods}}
\details{The output provides the list of interactions and the marginal distributions where those 
interactions are defined. The names of the involved variables are displayed 
if \code{names} is not NULL. For every interaction, the logit type used for each variable 
in the interaction set and the number of parameters are indicated. 
The last two columns give the position of the parameters in the vector
where all the interactions are arranged.}
\note{Functions `summary' and `print' display the same output.}
\seealso{\code{\link{hmmm.model}}}
\examples{
marginals<-marg.list(c("g-m","m-l","g-l"),mflag="m")
model<-hmmm.model(marg=marginals,lev=c(3,7),names=c("A","B"))
summary(model) # or print(model)}
\keyword{package}
