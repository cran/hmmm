\name{hmmm.chibar}
\alias{hmmm.chibar}
\title{chi-bar statistics test for hmm models}
\description{Function to simulate weights and pvalues of a chi-bar distributed statistic
for testing hypotheses of inequality constraints on parameters of hmm models.
The models in inputs are objects inheriting from class \code{hmmmfit}. 
}
\usage{
hmmm.chibar(nullfit, disfit, satfit, repli = 6000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nullfit}{The estimated model with inequalities turned into equalities, a result of `hmmm.mlfit'}
  \item{disfit}{The estimated model with inequalities, a result of `hmmm.mlfit'}
  \item{satfit}{The estimated model without inequalities, a result of `hmmm.mlfit'}
  \item{repli}{Number of simulations}
}
\details{
  The method "Simulation 2" described in Silvapulle and Sen, 2005, pg. 79 is used.
}
\value{
A list with the statistics test of type A and B (Silvapulle and Sen, 2005, pg. 61) 
and their simulated pvalues.
}
\references{Silvapulle MJ, Sen PK (2005) Constrained statistical
inference, Wiley, New Jersey.}
\examples{data(polbirth)
# 1 = Politics; 2 = Birthcontrol
y<-getnames(polbirth,st=12,sep=";")                     
names<-c("Pol","Birth")
marglist<-c("l-m","m-l","l-l")
marginals<-marg.list(marglist,mflag="m")
ineq<-list(marg=c(1,2),int=list(c(1,2)),types=c("l","l"))

# definition of the model with inequalities on interactions in ineq
model<-hmmm.model(marg=marginals,dismarg=list(ineq),lev=c(7,4),
strata=1,X=diag(1,27),names=names)

# saturated model
msat<-hmmm.mlfit(y,model)

# model with non-negative local log-odds ratios: "Likelihood ratio monotone dependence model"
mlr<-hmmm.mlfit(y,model,noineq=FALSE)

# model with null local log-odds ratios: "Stochastic independence model"
model0<-hmmm.model(marg=marginals,lev=c(7,4),sel=c(10:27),names=names)
mnull<-hmmm.mlfit(y,model0)

# HYPOTHESES TESTED:
#     testA --> H0=(mnull model) vs H1=(mlr model)
#     testB --> H0=(mlr model) vs H1=(msat model)

P<-hmmm.chibar(nullfit=mnull,disfit=mlr,satfit=msat)
summary(P)}
\seealso{\code{\link{chibar}}, \code{\link{summary.hmmmchibar}}, \code{\link{print.hmmmchibar}}}

\keyword{htest}% at least one, from doc/KEYWORDS
