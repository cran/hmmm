\name{inv_GMI}
\alias{inv_GMI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute a vector of joint probabilities from a vector of generalizedmarginal interactions (GMI)
}
\description{ Given an hmmm model and the vector of its  generalized interactions eta
 the vector of joint probabilities p is  computed by inverting
eta=C*ln(M*p)
}
\usage{
inv_GMI(etpar, mod, start = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{etpar}{
vector of GMI
}
  \item{mod}{
hmmm model corresponding to etapar; an object of class hmmmod created by hmmm.model
}
  \item{start}{
starting values for log-liner parameters in the non linear equations problem
}
}

\value{
vector of joint probabilities
}

\author{
Colombi R.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
p4<-c(
0.0895, 0.0351 ,0.0004, 0.0003, 0.0352, 0.2775, 0.0619, 0.0004, 0.0004, 0.0620, 0.2775,
0.0351 ,0.0001,
0.0004, 0.0352 ,0.089)

marg<-marg.list(c("l-m","m-l","l-l"), mflag="m") 
labelrisp<-c("R1","R2")
modello<-hmmm.model(marg=marg,lev=c(4,4),names=labelrisp)
etpar<-GMI(c(p4),c("l-m","m-l","l-l"),c(4,4),labelrisp,mflag="m")
etpar$gmi
p4rec<-inv_GMI(etpar$gmi,modello)
P<-cbind(p4rec,c(p4),c(p4)-p4rec)
colnames(P)<-c("prob","prob from eta","check")
P
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ generalized marginal interactions }
\keyword{marginal models }% __ONLY ONE__ keyword per line