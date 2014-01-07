\name{print.hidden}
\alias{print.hidden}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{print for the class hidden}
\description{
The generic function `print' is adapted to the objects inheriting from class \code{hidden} 
(print.hidden) to display the results of the estimation of a hidden model by `hidden.emfit'.
}
\usage{\method{print}{hidden}(x,printflag = FALSE,...)}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{An object of the class \code{hidden}, i.e. a result of `hidden.emfit'}
\item{printflag}{If TRUE the estimated parameters are displayed}
\item{\dots}{Further arguments passed to or from other methods}}
\value{The value of the log-likelihood functions and the estimated parameters of 
transition and observation models are printed.}
\seealso{
\code{\link{hidden.emfit}}, \code{\link{summary.hidden}}
}
\keyword{package}
