\name{summary.hmmmchibar}
\alias{summary.hmmmchibar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{summary for hmmm.chibar}
\description{
 Function to print the results for tests of type A and B (Silvapulle and Sen, 2005) on inequality constraints and to 
tabulate the chi-bar distribution functions of the statistics test.
The generic function `summary'
is adapted to the objects inheriting from class \code{hmmmchibar} (summary.hmmmchibar).
}
\usage{\method{summary}{hmmmchibar}(object, plotflag = 1, step = 0.01, lsup = 0,...)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of the class \code{hmmmchibar}}
  \item{plotflag}{1 to print only pvalues and statistic values,
2 to display the survival functions for type A and type B statistics tests and 3 to provide 
a plot of the survival functions (red: type B, black: type A)}
  \item{step}{ Distance between points at which the distribution functions are evaluated}
  \item{lsup}{ Distribution functions are evaluated in the interval 0 - lsup}
\item{\dots}{Further arguments passed to or from other methods}

}
\value{It provides the output
of the function `hmmm.chibar', that is the results of testing models defined by `hmmm.model' and
 estimated 
by `hmmm.mlfit', with equality and inequality constraints on marginal interactions. The statistics test of type A and B, their simulated pvalues and the chi-bar distribution functions 
are tabulated.}
\references{Silvapulle MJ, Sen PK (2005) Constrained statistical
inference, Wiley, New Jersey.}
\note{Use `print' for a short output.}
\seealso{\code{\link{hmmm.chibar}}, \code{\link{print.hmmmchibar}}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{package}
