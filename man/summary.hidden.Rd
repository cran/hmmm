\name{summary.hidden}
\alias{summary.hidden}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{summary for the class hidden}
\description{The generic function `summary' is adapted to the objects inheriting from class \code{hidden} 
(summary.hidden) to display the results of the estimation of a hidden model by `hidden.emfit'.}
\usage{\method{summary}{hidden}(object,...)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of the class \code{hidden}, i.e. a result of `hidden.emfit'}
\item{\dots}{Further arguments passed to or from other methods}}
\value{The transition probabilities and the probabilities of observations 
given the latent states are printed.}
\note{Use `print' to display the estimated parameters.}
\seealso{\code{\link{hidden.emfit}}, \code{\link{print.hidden}}}
\keyword{package}
