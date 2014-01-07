\name{accident}
\alias{accident}
\docType{data}
\title{
factory accident data
}
\description{
Data on factory accidents occurred in Bergamo (Italy) in 1998, collected by the Inail (Italian institute for insurance against factory accidents): 1052 workers who suffered an accident and claimed for a compensation are classified according to the type of injury, the time to 
recover (number of working days lost), the age (years), 
and the solar hour (part of the day in which the accident occurred).}
\usage{data(accident)}
\format{
  A data frame whose columns contain:
  \describe{
\item{\code{Type}}{ A factor with levels: \code{uncertain}, \code{avoidable}, \code{not-avoidable}}
 \item{\code{Time}}{A factor with levels: \code{0 |-- 7}, \code{7 |-- 21}, \code{21 |-- 60}, \code{>= 60}}
  \item{\code{Age}}{A factor with levels: \code{<= 25}, \code{26 -- 45}, \code{> 45}}
    \item{\code{Hour}}{A factor with levels: \code{morning}, \code{afternoon}}
     \item{\code{Freq}}{A numeric vector of frequencies}} }
\source{
Inail, Bergamo (Italy) 1998}
\references{
Cazzaro M, Colombi R (2008) Modelling two way contingency tables with recursive logits and odds ratios. Statistical Methods and Applications, 17, 435-453.}
\examples{data(accident)}
\keyword{datasets}
