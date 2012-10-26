\name{create.XMAT}
\alias{create.XMAT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{design matrix for a hmm model}
\description{
Function to specify the matrix X of the linear predictor Cln(Mm)=Xbeta for a hmm model.}
\usage{
create.XMAT(model, Formula = NULL, 
strata = 1, fnames = NULL, cocacontr = NULL, 
ncocacontr = NULL, replace = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{Object created by  `hmmm.model'}
  \item{Formula}{List of model-formulas; one formula for every marginal interaction}
  \item{strata}{Number of categories of the factors that describe the strata}
  \item{fnames}{Names of the factors that describe the strata}
  \item{cocacontr}{A list of zero-one matrices to build "r" logits created by the function `recursive'}
  \item{ncocacontr}{Number of contrasts for every factor, if NULL the maximum number is used}
  \item{replace}{If TRUE (default) a new model object with design matrix X is produced, if FALSE the list of design matrices associated to each element specified in Formula is returned}
}
\details{When the marginal interactions of a hmm model are defined in terms of
a linear predictor of covariates Cln(Mm)=Xbeta, the list of model-formulas defines additive effects of covariates on the interactions. 
In a case with 
two response variables declared by \code{names<-c("A","B")}
and two covariates, named C and D by \code{fnames=c("C","D")}, the additive effect of the covariates 
on marginal logits of A and B and 
log odds ratios (A.B) of the two responses is specified by the following \code{Formula}:
\code{Formula<-list(A=~A+C+D+A:C+A:D, B=~B+C+D+B:C+B:D,
A.B=~A.B+C+D+A.B:C+A.B:D)}. Use "zero" to constrain to zero all the interactions of a given type.
}
\value{A list of matrices or a hmm model with X as design matrix according to the input argument \code{replace}. The parameters beta in 
the predictor Cln(Mm)=Xbeta are the effects specified in \code{Formula} and correspond to the columns of X.}
\references{
Lang  JB (2004) Multinomial Poisson homogeneous models for contingency tables. The Annals of Statistics, 32, 340-383.

Lang JB (2005) Homogeneous linear predictor models for contingency tables. Journal of the American Statistical 
Association, 100, 121-134.}
\seealso{\code{\link{hmmm.model}}, \code{\link{hmmm.mlfit}}, \code{\link{summary.hmmmfit}}}
\examples{data(accident)
y<-getnames(accident,st=9,sep=";")
# responses: 1 = Type, 2 = Time; covariates: 3 = Age, 4 = Hour

marglist<-c("l-m","m-g","l-g")
marginals<-marg.list(marglist,mflag="m")
names<-c("Type","Time")

modelsat<-hmmm.model(marg=marginals,lev=c(3,4),
strata=6, names=names) 

# Create X to account for additive effect of Age and Hour on the logits of Type and Time
# and constant association between Type and Time
al<-list(Type=~Type+Age+Hour+Type:Age+Type:Hour,
Time=~Time+Age+Hour+Time:Age+Time:Hour,Type.Time=~Type.Time)
# list of matrices (replace=FALSE)
listmat<-create.XMAT(modelsat,Formula=al,strata=c(3,2),fnames=c("Age","Hour"),replace=FALSE)

# the model obtained by the modified X (replace=TRUE)
model<-create.XMAT(modelsat,Formula=al,strata=c(3,2),fnames=c("Age","Hour")) 
fitmodel<-hmmm.mlfit(y,model,y.eps=0.00001,maxit=2000,mup=1,step=1,m.initial=y)
print(fitmodel)
}
\keyword{models}