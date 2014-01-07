\name{hidden.emfit}
\alias{hidden.emfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ML estimation of a multinomial hidden Markov model}
\description{ Maximum likelihood estimation of a hidden Markov model with several categorical observed and latent variables.
 Observed and latent processes are specified by hmm models.}
\usage{
hidden.emfit(y, model.obs, model.lat, nlat = 1, noineq = TRUE, maxit = 10, 
maxiter = 100, 
norm.diff.conv = 1e-05, norm.score.conv = 1e-05, y.eps = 0, mup = 1, step = 1, 
printflag = 0, old.tran.p = NULL, bb = NULL, q.par=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{The observed multivariate categorical time series}
  \item{model.obs}{The model for the observations specified by `hmmm.model' or `hmmm.model.X'}
  \item{model.lat}{The model for the latent chain specified by `hmmm.model' or `hmmm.model.X'}
\item{nlat}{The number of latent variables}
  \item{noineq}{If TRUE  inequality constraints are not used}
  \item{maxit}{Maximum number of iterations for the M step}
  \item{maxiter}{Maximum number of iterations for the EM algorithm}
  \item{norm.diff.conv}{Convergence criterium for the parameters}
  \item{norm.score.conv}{Convergence criterium for the log-likelihood function}
  \item{y.eps}{Non-negative constant to be added to the original counts in y}
   \item{mup}{Weight for the constraints penalty part of the merit function}
  \item{step}{Interval length for the line search}
  \item{printflag}{If printflag=n the log-likelihood function is displayed any n iterations}
  \item{old.tran.p}{Starting values for the transition matrix}
  \item{bb}{Starting values for the observation probabilities}
\item{q.par}{The percentage of parameters that must satisfy the convergence criterium, q.par has values in [0 1]}}
\value{\item{vecpar}{List of two vectors of parameters of the observed and latent processes}
  \item{model.obs}{Information about the observed process}
\item{model.lat}{Information about the latent process}
\item{initial}{Invariant distribution of the latent process}
   \item{Ptr}{Estimated transition probabilities}
\item{Ptobs}{Estimated observation probabilities}
\item{Ptr.iniz}{Starting values of Ptr}
\item{Ptobs.iniz}{Starting values for Ptobs}
 \item{filter}{Filtered probabilities }
\item{smooth}{Smoothed probabilities}
\item{conv}{List of convergence criteria.}
Use `print' or `summary' to display the output.}
\details{Every column of \code{y} corresponds to an observed variable, every row in \code{y} reports the realization of the observed variable at each time occasion. The r realizations of each observed variable must be coded by the first r integers. The \code{model.lat} and \code{model.obs} are objects inheriting from class \code{hmmmmod}. 
So, the model for the transition matrix of the latent Markov chain and the model for the multinomial process can be marginal models specified by `hmmm.model' or `hmmm.model.X'. Consider a hidden Markov model for p observed variables and q latent variables. 
In defining the hmm models for observed and latent components using `hmmm.model' bear in mind that: for the observations, the first q variables are the latent variables, followed by the p observed variables; in the latent model, the first q variables refer to the latent at time t, while the remaining q indicate the latent variables at time (t-1). Note that in both cases the first marginal set must contain only the latent variables. On the other hand, in defining the hmm models for the observations and latent chain using `hmmm.model.X' consider that: for the observed model, the responses are the observed variables while the latent variables have the role of covariates and the number of their categories is declared in \code{strata}; in the latent model, the responses are the latent variables at time t, while the lagged variables at time (t-1) are considered as covariates. 
Declare the argument \code{nlat} only if \code{model.obs} is specified by `hmmm.model' and \code{model.lat} is not explicitly defined.}
\references{Colombi R, Giordano S (2011) Lumpability for discrete hidden Markov models. Advances in Statistical Analysis, 95(3), 293-311.}
\seealso{\code{\link{print.hidden}}, \code{\link{summary.hidden}}, \code{\link{hmmm.model}}, \code{\link{hmmm.model.X}}}
\examples{data(drinks)
y<-cbind(drinks$lemon.tea,drinks$orange.juice)
fm<-c("l-l-l")
fmargobs<-marg.list(fm)
#initial values of transition matrix and obs distribution given the two latent states
Ptr<-matrix(c(0.941, 0.199,0.059, 0.801),2,2,byrow=TRUE)
Ptobs<-matrix(c(0.053, 0.215, 0.206, 0.001, 0.039, 0.021, 0.020, 0.176, 0.270,
                0.000, 0.000, 0.000, 0.048, 0.263, 0.360, 0.065, 0.053, 0.211),
2,9,byrow=TRUE)
find<-~lat+lat*tea+lat*juice   # lat is the latent variable
model.obsf<-hmmm.model(marg=fmargobs,
lev=c(2,3,3),names=c("lat","tea","juice"),formula=find)

# model of independent observed variables given the latent states
modelind<-hidden.emfit(y,model.obsf,y.eps=0.01,maxit=10,maxiter=2500,
old.tran.p=Ptr,bb=Ptobs)
print(modelind,printflag=TRUE)

#alternative definition based on hmmm.model.X
f<-list(tea=~tea*lat,juice=~juice*lat,tea.juice="zero")
model.obsfX<-hmmm.model.X(marg=marg.list(c("l-l")),names=c("tea","juice"),
fnames=c("lat"),lev=c(3,3),strata=c(2))
modelindX<-hidden.emfit(y,model.obsfX,y.eps=0.01,maxit=10,maxiter=2500,
old.tran.p=Ptr,bb=Ptobs)
modelindX
summary(modelindX)
}
\keyword{ models}
