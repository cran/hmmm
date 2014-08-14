\name{hmmm.mlfit}
\alias{hmmm.mlfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{fit a hmm model}
\description{Function to estimate a hierarchical multinomial marginal model.
}
\usage{hmmm.mlfit(y, model, noineq = TRUE, maxit = 1000, 
norm.diff.conv = 1e-05, norm.score.conv = 1e-05,
y.eps = 0, chscore.criterion = 2,
m.initial = y, mup = 1, step = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{A vector of frequencies of the contingency table}
  \item{model}{An object created by `hmmm.model'}
  \item{noineq}{If TRUE inequality constraints specified in the model are ignored}
  \item{maxit}{Maximum number of iterations}
  \item{norm.diff.conv}{Convergence criterium value on the parameters}
  \item{norm.score.conv}{Convergence criterium value on the constraints}
  \item{y.eps}{Non-negative constant to be  added 
                    to the original counts in y}
  \item{chscore.criterion}{If equal to zero, convergence
 information are printed at every iteration}
  \item{m.initial}{Initial estimate of m (expected frequencies)}
  \item{mup}{ Weight for the constraints penalty part of the merit function}
  \item{step}{Interval length for the line search}
}
\details{A sequential quadratic procedure is used to maximize the log-likelihood function
 under inequality and equality constraints. 
This function calls the procedure `mphineq.fit'  which is a generalization of the procedure
`mph.fit' by Lang (2004).}
\value{An object of the class \code{hmmmfit}; an estimate of a marginal model defined by `hmmm.model'.
  The output can be displayed using `summary' or `print'.}
\references{ Bartolucci F, Colombi R, Forcina A (2007) An extended class of marginal link 
functions for modelling contingency tables by equality
 and inequality constraints. Statistica Sinica, 17, 691-711.

Bergsma WP, Rudas T (2002) Marginal models for categorical data. The Annals of Statistics, 30, 140-159. 

Colombi R, Giordano S, Cazzaro M (2014) hmmm: An R Package for hierarchical multinomial marginal models. Journal of Statistical Software, 59(11), 1-25, URL http://www.jstatsoft.org/v59/i11/.

Lang  JB (2004) Multinomial Poisson homogeneous models for contingency tables. The Annals of Statistics, 32, 340-383.}
\seealso{\code{\link{hmmm.model}}, \code{\link{hmmm.model.X}}, \code{\link{summary.hmmmfit}}, \code{\link{print.hmmmfit}}}
\examples{
data(relpol)
y<-getnames(relpol,st=12)
# 1 = Religion, 2 = Politics
names<-c("Rel","Pol")
marglist<-c("l-m","m-g","l-g")
marginals<-marg.list(marglist,mflag="m")

# Hypothesis of stochastic independence: all log odds ratios are null 
model<-hmmm.model(marg=marginals,lev=c(3,7),sel=c(9:20),names=names)
fitmodel<-hmmm.mlfit(y,model)
print(fitmodel, aname="Independence model",printflag=TRUE)
summary(fitmodel)

}
\keyword{ models}% at least one, from doc/KEYWORDS
