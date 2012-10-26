\name{madsen}
\alias{madsen}
\docType{data}
\title{Madsen data}
\description{The dataset concerns 1681 rental property residents classified according to their 
satisfaction from the house, perceived influence on the management of the property,  type of rental accommodation, 
and contact with other residents.}
\usage{data(madsen)}
\format{
  A data frame whose columns contain:
  \describe{\item{\code{Influence}}{A factor with levels: \code{low}, \code{medium}, \code{high}}
    \item{\code{Satisfaction}}{A factor with levels: \code{low}, \code{medium}, \code{high}}
    \item{\code{Contact}}{A factor with levels: \code{low}, \code{high}}
    \item{\code{Housing}}{A factor with levels: \code{tower block}, \code{apartment}, \code{atrium house}, \code{terraced house}}
 \item{\code{Freq}}{A numeric vector of frequencies}}}
\references{
Madsen M (1976) Statistical analysis of multiple contingency tables. Two examples. 
Scandinavian Journal of Statistics, 3, 97-106. 
}
\examples{
data(madsen)
}
\keyword{datasets}
