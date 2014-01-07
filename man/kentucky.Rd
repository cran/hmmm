\name{kentucky}
\alias{kentucky}
\docType{data}
\title{Kentucky traffic accident data}
\description{The traffic accident data collected by the Kentucky State Police from 1995 to 1999. 
The annual numbers of vehicle occupants involved in Kentucky accidents are classified according to 3 
variables: injury, restraint usage and year.
}
\usage{data(kentucky)}
\format{A data frame whose columns contain:
  \describe{\item{\code{Injury}}{A factor with levels: \code{1} = not injured; \code{2} = possible injury;
 \code{3} = nonincapacitating injury;  \code{4} = incapacitating injury;  \code{5} = killed}
\item{\code{Restraint.usage}}{A factor with levels: \code{yes} = restraint used, \code{no} = restraint not used}
\item{\code{Year}}{\code{1995}, \code{1996}, \code{1997}, \code{1998}, \code{1999}}
\item{\code{Freq}}{A numeric vector of frequencies}}}
\source{www.kentuckystatepolice.org/text/data.htm}
\references{Lang JB (2005) Homogeneous linear predictor models for contingency tables. 
Journal of the American Statistical Association, 100, 121-134.}
\examples{data(kentucky)}
\keyword{datasets}
