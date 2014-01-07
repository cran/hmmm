\name{depression}
\alias{depression}
\docType{data}
\title{
longitudinal study of mental depression}
\description{
A longitudinal study comparing a new drug with a standard drug for treatment of 
340 subjects suffering mental depression (Koch et al., 1977, Agresti, 2013). The patients are classified according 
to the severity of the initial diagnosis (mild and severe), the treatment type they received (standard and new drugs)
and the responses on 
the depression assessment (normal and abnormal) at three occasions: after one week (R1), 
two (R2), and four weeks (R3) of treatment.}
\usage{data(depression)}
\format{
  A data frame whose columns contain:
\describe{
\item{\code{R3}}{ A factor with levels: \code{N} as normal, \code{A} as abnormal}
\item{\code{R2}}{ A factor with levels: \code{N} as normal, \code{A} as abnormal}
\item{\code{R1}}{ A factor with levels: \code{N} as normal, \code{A} as abnormal}
 \item{\code{Treatment}}{A factor with levels: \code{standard}, \code{new drug}}
  \item{\code{Diagnosis}}{A factor with levels: \code{mild}, \code{severe}}
 \item{\code{Freq}}{A numeric vector of frequencies}} }
\references{ Agresti A (2013) Categorical Data Analysis (third edition). Wiley.

Koch GG,  Landis JR,  Freeman JL, Freeman DH and Lehnen RG  (1977) A general metodology for the analysis of experiments
with repeated measurement of categorical data. Biometrics, 38, 563-595.}
\examples{
data(depression)
}
\keyword{datasets}
