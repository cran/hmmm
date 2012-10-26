\name{summary.mphfit}
\alias{summary.mphfit}
\alias{print.mphfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
summary and print for the class mphfit
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
2. the linear predictor model results: estimated betas, standard errors, z-ratios, pvalues; estimated mph parameters,
standard errors, residuals; 3. cell-specific statistics: 
observed and predicted frequencies of the multi-way table, estimated joint probabilities with standard errors, adjusted residuals; 
4. convergence statistics. 5. model information: linear predictor model link function (L.fct), derivative of transpose L.fct,
linear predictor model design matrix  X, orthogonal complement of X, constraint function, derivative of transpose constraint function,
 population matrix Z, sampling constraint matrix ZF.}
\note{
Use `print' to display only the goodness-of-fit test and the convergence statistics.
}
\seealso{
\code{\link{mphineq.fit}}}
\keyword{package}
