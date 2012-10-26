\name{relpol}
\alias{relpol}
\docType{data}
\title{religion and political orientation data}
\description{Data on religion and political orientation of a sample of 911 U.S. citizens.}
\usage{data(relpol)}
\format{A data frame whose columns contain:
\describe{\item{\code{Religion}}{A factor with levels: \code{Protestant}, \code{Catholic}, \code{None}}
\item{\code{Politics}}{A factor with levels: \code{Extremely liberal}, \code{Liberal}, \code{Slightly liberal}, \code{Moderate},
\code{Slightly conservative}, \code{Conservative}, \code{Extremely conservative}} 
\item{\code{Freq}}{A numeric vector of frequencies}}}
\details{This is a sub-data frame obtained by marginalizing the data frame `relpolbirth'
with respect to the variable opinion on teenage birth control.}
\source{General Social Survey, 1993}
\references{Bergsma W, Croon M, Hagenaars JA (2009) Marginal models for dependent, clustered, and longitudinal categorical data. Springer.}
\examples{data(relpol)}
\keyword{datasets}
