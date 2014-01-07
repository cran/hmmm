\name{anova.hidden}
\alias{anova.hidden}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{anova for the class hidden}
\description{The generic function `anova' is adapted to the objects inheriting from class \code{hidden} 
(anova.hidden) to compute the likelihood ratio test for nested hidden models estimated by `hidden.emfit'.
}
\usage{\method{anova}{hidden}(object,objectlarge,...)}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{object}{Object of the class \code{hidden}, reduced model, i.e. modelA}
\item{objectlarge}{Object of the class \code{hidden}, large model, i.e. modelB}
\item{\dots}{Other models and further arguments passed to or from other methods}}
\details{
Nested models, fitted by `hidden.emfit', are compared (e.g. modelA is nested in modelB), 
the likelihood ratio statistic with the degrees of freedom and 
the associated pvalue is returned.
}
\examples{data(drinks)
y<-cbind(drinks$lemon.tea,drinks$orange.juice)

f<-~lat*tea+lat*juice+tea*juice  # lat indicates the latent variable
fm<-c("l-l-l")
fmargobs<-marg.list(fm,mflag="m")
Ptr<-matrix(c(0.941, 0.199,
               0.059, 0.801),2,2,byrow=TRUE)
Ptobs<-matrix(c(0.053, 0.215, 0.206, 0.001, 0.039, 0.021, 0.020, 0.176, 0.270,
                0.000, 0.000, 0.000, 0.048, 0.263, 0.360, 0.065, 0.053, 0.211)
                ,2,9,byrow=TRUE)

# saturated model (fsat<-~lat*tea*juice is implicit)
model.obsf<-hmmm.model(marg=fmargobs,
lev=c(2,3,3),names=c("lat","tea","juice"))
modelsat<-hidden.emfit(y,model.obsf,y.eps=0.01,maxit=10,
maxiter=2500,norm.diff.conv=0.001,old.tran.p=Ptr,bb=Ptobs)

# model with constant association
model.coass<-hmmm.model(marg=fmargobs,
lev=c(2,3,3),names=c("lat","tea","juice"),formula=f)
modelca<-hidden.emfit(y,model.coass,y.eps=0.01,maxit=10,
maxiter=2500,norm.diff.conv=0.001,old.tran.p=Ptr,bb=Ptobs)

a<-anova(modelca,modelsat)}
\seealso{
\code{\link{hidden.emfit}}, \code{\link{summary.hidden}}, \code{\link{print.hidden}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{package}