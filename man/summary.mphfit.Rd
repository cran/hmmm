\name{summary.mphfit}
\alias{summary.mphfit}
\alias{print.mphfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
summary and print for the class \code{mphfit}
}
\description{The generic functions `summary' and `print' are adapted to the objects inheriting from class \code{mphfit}
(summary.mphfit, print.mphfit) to display the results of the estimation of a mph model by `mphineq.fit'.}
\usage{\method{summary}{mphfit}(object, ...)
\method{print}{mphfit}(x,...)}
\arguments{
\item{object,x}{An object inheriting from class \code{mphfit}, i.e. a result of `mphineq.fit'}
\item{\dots}{Further arguments passed to or from other methods}}
\details{The output of `summary' provides: 1. the goodness-of-fit of the estimated model tested by the likelihood ratio and
Pearson's Score Statistics, degrees of freedom and pvalues. Note that degrees of freedom and pvalues are meaningful only for 
the mph models without inequality constraints; 
2. cell-specific statistics: 
observed and predicted frequencies of the multi-way table, estimated joint probabilities with standard errors, adjusted residuals.}
\value{No return value}
\note{
Use `print' to display only the goodness-of-fit test.
}
\seealso{
\code{\link{mphineq.fit}}}
\keyword{package}
