\name{polbirth}
\alias{polbirth}
\docType{data}
\title{political orientation and teenage birth control data}
\description{Data on political orientation and opinion on teenage birth control of a sample of 911 U.S. citizens.}
\usage{data(polbirth)}
\format{A data frame whose columns contain:
\describe{\item{\code{Politics}}{A factor with levels: \code{Extremely liberal}, \code{Liberal}, \code{Slightly liberal}, \code{Moderate},
\code{Slightly conservative}, \code{Conservative}, \code{Extremely conservative}}
\item{\code{Birthcontrol}}{A factor with levels: \code{Strongly agree}, \code{Agree}, \code{Disagree}, \code{Strongly disagree}} 
\item{\code{Freq}}{A numeric vector of frequencies}}}
\details{This is a sub-data frame obtained by marginalizing the data frame `relpolbirth'
with respect to the variable religion.}
\source{General Social Survey, 1993.}
\references{Bergsma W, Croon M, Hagenaars JA (2009) Marginal Models for Dependent, Clustered, and Longitudinal Categorical Data. Springer.}
\examples{data(polbirth)}
\keyword{datasets}
