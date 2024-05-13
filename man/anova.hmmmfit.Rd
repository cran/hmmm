\name{anova.hmmmfit}
\alias{anova.hmmmfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{anova for the class hmmmfit}
\description{The generic function `anova' is adapted to the objects inheriting from class \code{hmmmfit} 
(anova.hmmmfit) to compute the likelihood ratio test for nested hmm models estimated 
by `hmmm.mlfit'.
}
\usage{\method{anova}{hmmmfit}(object,objectlarge,...)}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{object}{Object of the class \code{hmmmfit}, reduced model, i.e. modelA}
\item{objectlarge}{Object of the class \code{hmmmfit}, large model, i.e. modelB}
\item{\dots}{Other models and further arguments passed to or from other methods}}
\details{
Nested models, fitted by `hmmm.mlfit', are compared (e.g. modelA is nested in modelB), 
the likelihood ratio statistic with the degrees of freedom and 
the associated pvalue is printed.
}
\value{ A matrix with information about Likelihood ratio tests degree of freedom and P values for the compared models.}
\examples{data(madsen)
y<-getnames(madsen)
names<-c("Infl","Sat","Co","Ho")

fA<-~Co*Ho+Sat*Co+Sat*Ho
modelA<-loglin.model(lev=c(3,3,2,4),formula=fA,names=names)
fB<-~Co*Ho+Sat*Co+Infl*Co+Sat*Ho+Infl*Sat
modelB<-loglin.model(lev=c(3,3,2,4),formula=fB,names=names)

modA<-hmmm.mlfit(y,modelA)
modB<-hmmm.mlfit(y,modelB)

anova(modA,modB)
}
\seealso{
\code{\link{hmmm.mlfit}}, \code{\link{summary.hmmmfit}}, \code{\link{print.hmmmfit}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{package}