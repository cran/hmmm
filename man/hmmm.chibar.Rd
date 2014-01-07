\name{hmmm.chibar}
\alias{hmmm.chibar}
\title{chi-bar statistic test for hmm models}
\description{Function to calculate weights and pvalues of a chi-bar-square distributed statistic
for testing hypotheses of inequality constraints on parameters of hmm models.
The models in input are objects inheriting from class \code{hmmmfit} or \code{mphfit}. 
}
\usage{
hmmm.chibar(nullfit, disfit, satfit, repli = 6000, 
kudo = FALSE, TESTAB = FALSE, 
alpha = c(0.02,0.03,0), pesi = NULL)
}
\arguments{
  \item{nullfit}{The estimated model with inequalities turned into equalities}
  \item{disfit}{The estimated model with inequalities}
  \item{satfit}{The estimated model without inequalities}
  \item{repli}{Number of simulations}
\item{kudo}{If TRUE, the chi-bar weights are not simulated but  computed by the Kudo's method}
   \item{TESTAB}{If TRUE, the LR tuned testing procedure is performed (see Details)}
   \item{alpha}{Three significance levels c(alpha1, alpha2, alpha12) of the LR tuned testing procedure}
   \item{pesi}{The chi-bar weights if they are known}
}
\details{All the 3 argument models must be obtained by `hmmm.mlfit' or by `mphineq.fit'.
  The method "Simulation 2" described in Silvapulle and Sen, 2005, pg. 79 is used if \code{kudo = FALSE}, otherwise the Kudo's exact method is used
as described by El Barmi and Dykstra (1999). The Kudo's method can be reasonably used with less than 10-15 inequalities.
If TESTA is the LR statistics  for \code{nullfit}  against the \code{disfit} model
 while TESTB  is the LR statistics for \code{disfit} against the \code{satfit} model then 
the LR tuned testing procedure (Colombi and Forcina, 2013) runs as follows:
accept \code{nullfit} if TESTB < y2 and TESTA < y1,
where 
Pr(TESTB > y2) = alpha2-alpha12 and  Pr(TESTA < y1, TESTB < y2) = 1-alpha1-alpha2,
reject \code{nullfit}  in favour of \code{disfit}  if TESTA > y1 and
TESTB < y12, where Pr(TESTA > y1, TESTB < y12) = alpha1,
otherwise reject \code{nullfit} for \code{satfit}.
}
\value{A list with the statistics test of type A and B (Silvapulle and Sen, 2005, pg. 61) 
and their pvalues. If \code{TESTAB = TRUE} details on the LR tuned testing procedure (Colombi and Forcina, 2013) are reported.
}
\references{ Colombi R. Forcina A. (2013) Testing order restrictions in contingency tables. Submitted.

El Barmi H, Dykstra R (1999) Likelihood ratio test against a set of inequality constraints. 
Journal of Nonparametric Statistics, 11, 233-261. 


Silvapulle MJ, Sen PK (2005) Constrained statistical
inference, Wiley, New Jersey.
}
\examples{data(polbirth)
# 1 = Politics; 2 = Birthcontrol
y<-getnames(polbirth,st=12,sep=";")                     
names<-c("Pol","Birth")
marglist<-c("l-m","m-l","l-l")
marginals<-marg.list(marglist,mflag="m")
ineq<-list(marg=c(1,2),int=list(c(1,2)),types=c("l","l"))

# definition of the model with inequalities on interactions in ineq
model<-hmmm.model(marg=marginals,dismarg=list(ineq),lev=c(7,4),names=names)

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
\seealso{\code{\link{summary.hmmmchibar}}, \code{\link{print.hmmmchibar}}}

\keyword{htest}% at least one, from doc/KEYWORDS
