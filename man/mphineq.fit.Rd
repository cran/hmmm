\name{mphineq.fit}
\alias{mphineq.fit}
\title{fit mph models under inequality constraints}
\description{Function to maximize the log-likelihood function of multinomial Poisson homogeneous (mph) models
under nonlinear equality and inequality constraints.}
\usage{mphineq.fit(y, Z, ZF = Z, h.fct = 0, derht.fct = 0, d.fct = 0,
derdt.fct = 0, L.fct = 0, derLt.fct = 0, X = NULL, formula = NULL, 
names = NULL, lev = NULL, E = NULL, maxiter = 100, 
step = 1, norm.diff.conv = 1e-05, norm.score.conv = 1e-05, 
y.eps = 0, chscore.criterion = 2, m.initial = y, mup = 1)}
\arguments{\item{y}{Vector of frequencies of the multi-way table}
\item{Z}{Population matrix. The population matrix Z is a c x s zero-one matrix, where c is the number of counts 
and s is the number of strata or populations. Thus, the rows correspond to the number of observations  and 
the columns correspond to the strata.  A 1 in row i and column j means that the ith count comes 
from the jth stratum.  Note that Z has exactly one 1 in each row, and at least one 1 in each 
column.  When Z = matrix(1,length(y),1) is a column vector of 1, all 
the counts come from the same and only stratum.}
\item{ZF}{Sample constraints matrix. For non-zero ZF, the columns are a subset of the columns in the population matrix Z,
then the sample size of the jth stratum is considered fixed, 
otherwise if the jth column of Z is NOT included in ZF, the jth stratum sample size is taken to be a realization 
of a Poisson random variable.  When ZF=0, all of the stratum sample sizes are taken to be realizations of
Poisson random variables.  The default, ZF=Z, means that all the stratum sample sizes are fixed; this 
is the (product-)multinomial setting.  Note that ZF'y = n is the vector of fixed sample sizes}
\item{h.fct}{Function h(m) of equality constraints, m is the vector of expected frequencies. This function of m must return a vector}
\item{derht.fct}{Derivative of h(m), if not supplied numerical derivative are used}
\item{d.fct}{Function for inequality constraints d(m)>0. This function of m must return a vector}
\item{derdt.fct}{Derivative of d(m), if not supplied numerical derivative are used}
\item{L.fct}{Link function for the linear model L(m)=Xbeta}
\item{derLt.fct}{Derivative of L(m), if not supplied numerical derivative are used}
\item{X}{Model matrix for L(m)=Xbeta}
\item{formula}{Formula of the reference log-linear model}
\item{names}{A character vector whose elements are the names of the variables}
\item{lev}{Number of categories of the variables}
\item{E}{If E is a matrix, then X is ignored and E defines the equality contrasts as
EL(m)=0}
\item{maxiter}{Maximum number of iterations}
\item{step}{ Interval length for the linear search}
\item{norm.diff.conv}{Convergence criterium for parameters}
\item{norm.score.conv}{Convergence criterium for constraints}
\item{y.eps}{ Non-negative constant to be temporarily added to the original frequencies in y}
\item{chscore.criterion}{If zero, convergence information are printed at every iteration}
\item{m.initial}{Initial estimate of m}
\item{mup}{Weight for the constraint part of the merit function}}
\details{This function extends `mph.fit' written by JB Lang, Dept of Statistics and Actuarial Science
University of Iowa, in order to include inequality constraints.
In particular, the Aitchison Silvey (AS) algorithm has been replaced by a sequential 
quadratic algorithm which is equivalent to AS when inequalities are not present.
The R functions `quadprog' and `optimize' have been used to implement 
the sequential quadratic algorithm.
More precisely, the AS updating formulas are replaced 
by an equality-inequality constrained quadratic programming problem.
The `mph.fit' step halving linear search is replaced by an optimal
step length search performed by `optimize'. 
}
\references{Lang  JB (2004) Multinomial Poisson homogeneous models for contingency tables. The Annals of Statistics, 32, 340-383.

Lang JB (2005) Homogeneous linear predictor models
for contingency tables. Journal of the American Statistical Association, 100, 121-134.
}
 \seealso{\code{\link{print.mphfit}}, \code{\link{summary.mphfit}}, \code{\link{hmmm.model}}, \code{\link{hmmm.mlfit}}}
\author{Roberto Colombi} 
\examples{
y <- c(104,24,65,76,146,30,50,9,166) # Table 2 (Lang, 2004)
y <- matrix(y,9,1)

# population matrix: 3 strata with 3 observations each
Z <- kronecker(diag(3),matrix(1,3,1))
# the 3rd stratum sample size is fixed
ZF <- kronecker(diag(3),matrix(1,3,1))[,3]

# Let P_ij be the expected number of cross-citations, P_ij=P(A=i,B=j),
# where A = Citing journal and B = Cited journal. 
# The Gini concentrations of citations for each of the journals are: 
# G_i = sum_j=1_3 (P_ij/P_i+)^2  for i=1,2,3.
\donttest{
Gini.fct <- function(m) {
p <- diag(c(1/(Z%*%t(Z)%*%m)))%*%m 
G1 <- sum(diag(matrix(c(p[1]/(p[1]+p[2]+p[3]),0,0,
                  0,p[2]/(p[1]+p[2]+p[3]),0,
                  0,0,p[3]/(p[1]+p[2]+p[3])),3,3,byrow=TRUE)*
         matrix(c(p[1]/(p[1]+p[2]+p[3]),0,0,
                  0,p[2]/(p[1]+p[2]+p[3]),0,
                  0,0,p[3]/(p[1]+p[2]+p[3])),3,3,byrow=TRUE)))
G2 <- sum(diag(matrix(c(p[4]/(p[4]+p[5]+p[6]),0,0,
                  0,p[5]/(p[4]+p[5]+p[6]),0,
                  0,0,p[6]/(p[4]+p[5]+p[6])),3,3,byrow=TRUE)*
         matrix(c(p[4]/(p[4]+p[5]+p[6]),0,0,
                  0,p[5]/(p[4]+p[5]+p[6]),0,
                  0,0,p[6]/(p[4]+p[5]+p[6])),3,3,byrow=TRUE)))
G3 <- sum(diag(matrix(c(p[7]/(p[7]+p[8]+p[9]),0,0,
                  0,p[8]/(p[7]+p[8]+p[9]),0,
                  0,0,p[9]/(p[7]+p[8]+p[9])),3,3,byrow=TRUE)*
         matrix(c(p[7]/(p[7]+p[8]+p[9]),0,0,
                  0,p[8]/(p[7]+p[8]+p[9]),0,
                  0,0,p[9]/(p[7]+p[8]+p[9])),3,3,byrow=TRUE)))
c(G1,G2,G3)}

# h_1 = c(G1,G2,G3)-c(0.410,0.455,0.684) = 0
# HYPOTHESIS: no change in Gini concentrations 
# from the 1987-1989 observed values

h.fct <- function(m) {G<-Gini.fct(m)
c(G[1],G[2],G[3])-c(0.410,0.455,0.684)}

mod_eq <- mphineq.fit(y,Z,ZF,h.fct=h.fct)

print(mod_eq)

# Example of MPH model subject to inequality constraints 

# d_1 = c(G1,G2,G3)-c(0.410,0.455,0.684) >= 0
# HYPOTHESIS: increase in Gini concentrations
# from the 1987-1989 observed values

d.fct <- function(m) {G<-Gini.fct(m)
c(G[1],G[2],G[3])-c(0.410,0.455,0.684)}

mod_ineq <- mphineq.fit(y,Z,ZF,d.fct=d.fct)
print(mod_ineq)

# HYPOTHESES TESTED:
# NB: testA --> H0=(mod_eq) vs H1=(mod_ineq model)
#     testB --> H0=(mod_ineq model) vs H1=(sat_mod model)
# sat_mod --> saturated model (Gsq=0)

m <- mod_ineq$m

chibar(m,Z,ZF,d.fct=d.fct,test0=mod_eq$Gsq-mod_ineq$Gsq,
test1=mod_ineq$Gsq,repli=6000,lev=c(3,3))
}
}
\keyword{models}