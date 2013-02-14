\name{print.hmmmchibar}
\alias{print.hmmmchibar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{print for hmmm.chibar}
\description{
 Function to print the results for tests of type A and B (Silvapulle and Sen, 2005) on 
inequality constraints.
The generic function `print' 
is adapted to the objects inheriting from class \code{hmmmchibar} (print.hmmmchibar).
}
\usage{\method{print}{hmmmchibar}(x,...)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of the class \code{hmmmchibar}}
\item{\dots}{Further arguments passed to or from other methods}
}
\value{It provides the output
of the function `hmmm.chibar', that is the results of testing models defined by `hmmm.model' and
 estimated 
by `hmmm.mlfit', with equality and inequality constraints on marginal interactions. The statistics test of type A and B and 
their simulated pvalues are tabulated.}
\references{Silvapulle MJ, Sen PK (2005) Constrained statistical
inference, Wiley, New Jersey.}
\note{Use `summary' to display a much detailed output.}
\seealso{\code{\link{hmmm.chibar}}, \code{\link{summary.hmmmchibar}}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{package}
