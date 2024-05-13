\name{print.hmmmfit}
\alias{print.hmmmfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{print for the class hmmmfit}
\description{
The generic function `print' is adapted to the objects inheriting from class \code{hmmmfit} 
(print.hmmmfit) to display the results of the estimation of a hmm model by `hmmm.mlfit'.
}
\usage{\method{print}{hmmmfit}(x, aname = " ", printflag = FALSE, ...)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of the class \code{hmmmfit}, i.e. a result of `hmmm.mlfit'}
  \item{aname}{The name of the fitted object model}
  \item{printflag}{If FALSE only the goodness-of-fit test is displayed,
 if TRUE the estimates of the interaction parameters are also returned}
\item{\dots}{Further arguments passed to or from other methods}}
\details{
The printed output provides the likelihood ratio statistic test to assess the fitting of the model estimated by `hmmm.mlfit'.
Degrees of freedom and pvalues are meaningful only for the hmm models without inequality constraints (see `hmmm.chibar' to 
test hmm models defined under inequality constraints on interactions). 
Moreover, if \code{printflag} is TRUE, the estimated interactions are displayed for every stratum, together with 
the marginal sets where they are defined and the type of logits considered.
}
\value{No return value}
\note{Use `summary' to display a much detailed output.}
\examples{data(relpol)
y<-getnames(relpol,st=12)
# 1 = Religion, 2 = Politics
names<-c("Rel","Pol")
marglist<-c("l-m","m-g","l-g")
marginals<-marg.list(marglist,mflag="m")

# Hypothesis of stochastic independence: all log odds ratios are null 
model<-hmmm.model(marg=marginals,lev=c(3,7),sel=c(9:20),names=names)
fitmodel<-hmmm.mlfit(y,model)
print(fitmodel,aname="independence model",printflag=TRUE)
# summary(fitmodel)
}
\seealso{
\code{\link{hmmm.mlfit}}, \code{\link{summary.hmmmfit}}, \code{\link{anova.hmmmfit}}
}
\keyword{package}
