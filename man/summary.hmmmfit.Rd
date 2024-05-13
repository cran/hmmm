\name{summary.hmmmfit}
\alias{summary.hmmmfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
summary for the class hmmmfit
}
\description{The generic function `summary' is adapted to the objects inheriting from class \code{hmmmfit} 
(summary.hmmmfit) to display the results of the estimation of a hmm model by `hmmm.mlfit'.
}
\usage{\method{summary}{hmmmfit}(object, cell.stats = TRUE, ...)}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{object}{An object of the class \code{hmmmfit}, i.e. a result of `hmmm.mlfit'}
\item{cell.stats}{If TRUE cell-specific statistics are returned}\item{\dots}{Further arguments passed to or from other methods}}
\details{The marginal interactions of a hmm model can be defined in terms of linear predictor of covariates Cln(Mm)=Xbeta, where
 the X matrix is specified by `create.XMAT' and 
the parameters beta indicate the additive effects of covariate on the marginal interactions. 
The function `hmmm.mlfit' estimates either the parameters beta and the interactions; the function `summary' of a fitted model 
(by `hmmm.mlfit') returns the estimated betas and the estimated interactions, while the function `print' provides the 
estimated interactions only. If the model is defined under  equality constraints ECln(Mm)=0, parameters betas are meaningless so they are not printed. 

The printed output of `summary' provides: 1. values of the likelihood ratio and
Pearson's score statistics, degrees of freedom and pvalues. Note that degrees of freedom and pvalues are meaningful only for 
the hmm models without inequality constraints (see `hmmm.chibar' to test hmm models defined under inequality constraints on interactions); 
2. the linear predictor model results: estimated betas, standard errors, z-ratios, pvalues; estimated interactions, 
standard errors, residuals; 3. cell-specific statistics: 
observed and predicted frequencies of the multi-way table, estimated joint probabilities with standard errors, adjusted residuals; 
4. convergence statistics.}
\value{No return value}
\note{Use `print' to display only the goodness-of-fit test and the estimated interactions.}
\examples{
data(relpol)
y<-getnames(relpol,st=12,sep=";")
# 1 = Religion, 2 = Politics
names<-c("Rel","Pol")
marglist<-c("l-m","m-g","l-g")
marginals<-marg.list(marglist,mflag="m")

# Hypothesis of stochastic independence: all log odds ratios are null 
model<-hmmm.model(marg=marginals,lev=c(3,7),sel=c(9:20),names=names)
fitmodel<-hmmm.mlfit(y,model)
# print(fitmodel,aname="Independence model",printflag=TRUE)
summary(fitmodel)
}
\seealso{
\code{\link{hmmm.mlfit}}, \code{\link{print.hmmmfit}}, \code{\link{anova.hmmmfit}}, \code{\link{create.XMAT}}
}
\keyword{package}
