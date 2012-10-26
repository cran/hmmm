\name{hidden.emfit}
\alias{hidden.emfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ML estimation of a multinomial hidden Markov model}
\description{ Maximum likelihood estimation of a hidden Markov model with several observed categorical variables.
 Observed and latent processes are specified by hmm models.}
\usage{
hidden.emfit(y, model.obs, model.lat, noineq = TRUE, maxit = 10, maxiter = 100, 
norm.diff.conv = 1e-05, norm.score.conv = 1e-05, y.eps = 0, mup = 1, step = 1, 
printflag = 0, old.tran.p = NULL, bb = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{The observed multivariate categorical time series}
  \item{model.obs}{The model for the observations specified by `hmmm.model'}
  \item{model.lat}{The model for the latent chain specified by `hmmm.model'}
  \item{noineq}{If TRUE  inequality constraints are not used}
  \item{maxit}{Maximum number of iterations for the M step}
  \item{maxiter}{Maximum number of iterations for the EM algorithm}
  \item{norm.diff.conv}{Convergence criterium for the parameters}
  \item{norm.score.conv}{Convergence criterium for the log-likelihood function}
  \item{y.eps}{Non-negative constant to be added to the original counts in y}
   \item{mup}{Weight for the constraints penalty part of the merit function}
  \item{step}{Interval length for the line search}
  \item{printflag}{If printflag=n the log-likelihood function is displayed at any n iterations}
  \item{old.tran.p}{Starting values for the transition matrix}
  \item{bb}{Starting values for the observation probabilities}
}
\details{ The model for the transition matrix of the latent Markov chain and the model for the multinomial process
are marginal models specified by `hmmm.model'. In defining the hmm model for the observations, variable 1 indicates the latent 
variable and its logits must not be constrained. In the model for the latent variable,
2 is the variable at time (t-1) and 1 the variable at time t. Logits of variable 2 must not be constrained.}
\value{\item{vecpar}{List of two vectors of parameters of the observed and latent processes}
  \item{initial}{Invariant distribution of the latent process}
   \item{Ptr}{Estimated transition probabilities}
\item{Ptobs}{Estimated observation probabilities}
\item{Ptr.iniz}{Starting values of Ptr}
\item{Ptobs.iniz}{Starting values for Ptobs}
 \item{filter}{Filtered probabilities }
\item{smooth}{Smoothed probabilities}
\item{conv}{List of convergence criteria}
}
\references{Colombi R, Giordano S (2011) Lumpability for discrete hidden Markov models. Advances in Statistical Analysis, 95(3), 293-311.}
\seealso{\code{\link{print.hidden}}, {\code{\link{summary.hidden}}}}
\examples{data(drinks)
y<-cbind(drinks$lemon.tea,drinks$orange.juice)
fm<-c("l-l-l")
fmargobs<-marg.list(fm,mflag="m")
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
}
\keyword{ models}
