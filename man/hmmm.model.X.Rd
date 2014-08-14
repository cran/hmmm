\name{hmmm.model.X}
\alias{hmmm.model.X}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{hmm model with covariates effect on parameters}
\description{
Function to define a hmm model whose parameters depend on covariates.}
\usage{
hmmm.model.X(marg, lev, names, Formula = NULL, strata = 1, 
fnames = NULL, cocacontr = NULL, ncocacontr = NULL, replace=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{ \item{marg}{A list of the marginal sets and their  marginal interactions as
              described in Bartolucci et al. (2007). See details of \code{hmmm.model}}
  \item{lev}{Number of categories of the response variables}
\item{names}{A character vector whose elements are the names of the response variables}
 \item{Formula}{List of model-formulas; one formula for every marginal interaction}
  \item{strata}{Number of categories of the covariates that describe the strata}
  \item{fnames}{Names of the covariates that describe the strata}
  \item{cocacontr}{A list of zero-one matrices to build "r" logits created by the function `recursive'}
  \item{ncocacontr}{Number of contrasts for every covariate, if NULL the maximum number is used}
  \item{replace}{If TRUE a new model object with design matrix X is produced, if FALSE the list of design matrices associated to each element specified in Formula is returned}
}
\details{The arguments \code{names} and \code{fnames} report the names of responses and covariates according to the order in which the variables are declared, see details of function `hmmm.model'. 

When the marginal interactions of a hmm model are defined in terms of
a linear predictor of covariates Cln(Mm)=Xbeta, the list of model formulas defines additive effects of covariates on the interactions. 
In a case with 
two response variables declared by \code{names<-c("A","B")}
and two covariates, named C and D by \code{fnames=c("C","D")}, the additive effect of the covariates 
on marginal logits of A and B and 
log odds ratios (A.B) of the two responses is specified by the following \code{Formula}:
\code{Formula<-list(A=~A*(C+D), B=~B*(C+D), A.B=~A.B*(C+D))}. Use "zero" to constrain to zero all the interactions of a given type. The saturated model is the default if \code{Formula} is not specified.
}
\value{An object of the class \code{hmmmmod}; it describes a marginal model with effects of covariates on the interactions. This model can be estimated by `hmmm.mlfit'.}
\references{Colombi R, Giordano S, Cazzaro M (2014) hmmm: An R Package for hierarchical multinomial marginal models. Journal of Statistical Software, 59(11), 1-25, URL http://www.jstatsoft.org/v59/i11/.

Glonek GFV, McCullagh P (1995) Multivariate logistic models for
contingency tables. Journal of the Royal Statistical Society, B, 57, 533-546.

Marchetti GM, Lupparelli M (2011) Chain graph models of multivariate regression type for categorical data. Bernoulli, 17, 827-844.}
\seealso{\code{\link{hmmm.model}}, \code{\link{create.XMAT}}, \code{\link{summary.hmmmmod}}, \code{\link{print.hmmmmod}}, 
\code{\link{marg.list}}, \code{\link{recursive}}, \code{\link{hmmm.mlfit}}}
\examples{data(accident)
y<-getnames(accident,st=9,sep=";")
# responses: 1 = Type, 2 = Time; covariates: 3 = Age, 4 = Hour

marginals<-marg.list(c("b-marg","marg-g","b-g"))
al<-list(
Type=~Type*(Age+Hour),
Time=~Time*(Age+Hour),
Type.Time=~Type.Time*(Age+Hour)
)
# model with additive effect of the covariates on logits and log-o.r. of the responses
model<-hmmm.model.X(marg=marginals,lev=c(3,4),names=c("Type","Time"),
Formula=al,strata=c(3,2),fnames=c("Age","Hour"))
mod<-hmmm.mlfit(y,model,y.eps=0.1)
}
\keyword{models}