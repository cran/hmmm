\name{chibar}
\alias{chibar}
\title{simulation of chi-bar pvalues}
\description{Function to simulate weights and pvalues of a chi-bar distribution for tests 
of type A and type B (Silvapulle MJ and Sen PK (2005) Constrained statistical
inference, Wiley, New Jersey) on inequality constraints
on the parameters of hlp, mph models.
}
\usage{
chibar(m, Z, ZF, d.fct = 0, h.fct = 0, 
test0 = 0, test1 = 0, repli = 0, derdt.fct = 0,
derht.fct = 0, formula = NULL, names = NULL, lev)
}
\arguments{
  \item{m}{Estimated expected frequencies}
  \item{Z}{Population matrix - see the help of `mphineq.fit'}
  \item{ZF}{ Sample matrix - see the help of `mphineq.fit'}
  \item{d.fct}{ Inequality constraints function - see the help of `mphineq.fit'}
  \item{h.fct}{ Equality constraints function - see the help of `mphineq.fit'}
  \item{test0}{ Likelihood ratio statistics for testing problems of type A (Silvapulle and Sen, 2005)}
  \item{test1}{ Likelihood ratio statistics for testing problems of type B (Silvapulle and Sen, 2005)}
  \item{repli}{ Number of simulations}
  \item{derdt.fct}{ Derivatives  of inequality constraints - see the help of `mphineq.fit'} 
  \item{derht.fct}{ Derivatives  of equality constraints - see the help of `mphineq.fit'}
  \item{formula}{Formula of the reference log-linear model (see details in `hmmm.model')} 
\item{names}{A character vector whose elements are the names of the variables}
\item{lev}{Number of categories of the variables}
}
\details{
The method "Simulation 2" described in Silvapulle and Sen, 2005, pg. 79 is used.
 }
\value{
 A list with the likelihood ratio statistics for testing hypotheses of type A and B (Silvapulle and Sen, 2005, pg. 61) 
and their simulated pvalues.
}
\references{Silvapulle MJ, Sen PK (2005) Constrained statistical
inference, Wiley, New Jersey.}
\note{Use `summary' to display a detailed output, use `hmmm.chibar' for testing
hmm models.}
\seealso{\code{\link{hmmm.chibar}}, \code{\link{summary.hmmmchibar}}, \code{\link{print.hmmmchibar}}}

\keyword{htest}
