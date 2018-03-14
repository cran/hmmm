\name{GMI}
\alias{GMI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function  to compute the generalized marginal interactions associated to a hierarchical family of marginal sets
}
\description{
Given a vector of joint probabilities, the generalized marginal interactions (gmi) associated to a hierarchical family of marginal sets are computed.
If the input is a matrix, gmi are computed for every column.
}

\usage{
GMI(freq, marg, lev, names, mflag = "M")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{freq}{
Matrix of joint probabilities. Every column describes a joint pdf.
}
  \item{marg}{
A character vector  decribing  the marginal sets and the logits used to build the interactions. See \code{marg.list}
}
  \item{lev}{Number of categories of the categorical variables. See the help of \code{hmmm.model} 
%%     ~~Describe \code{lev} here~~
}
  \item{names}{Names of the categorical variables
%%     ~~Describe \code{names} here~~
}
  \item{mflag}{
The symbol used to denote variables that are marginalized, default "M". See \code{marg.list}
}
}
\value{
A list with two components: marginals and gmi; marginals is a legend that explains the interactions, gmi is a vector or a matrix that contains the interactions.
}
\references{
Colombi R, Giordano S, Cazzaro M (2014) hmmm: An R Package for hierarchical multinomial marginal models. Journal of Statistical Software, 59(11), 1-25, URL http://www.jstatsoft.org/v59/i11/.
}


\seealso{\code{\link{inv_GMI}}, \code{\link{hmmm.model}}, \code{\link{marg.list}}}

\examples{
# joint frequencies for two ordinal variables 
# H: level of happiness on a scale from 1 to 5
# S: level of satisfaction on a scale from 1 to 5

y<-c(50,36,15,15,13,15,84,60,42,
     35,6,26,105,113,57,5,26,62,
     465,334,4,10,34,186,1404)
  
lev<-c(5,5)
marg<-c("g-m","m-g","g-g")
names<-c("H","S")

o<-GMI(cbind(c(y),c(y/sum(y))),marg,lev,names,mflag="m")
o
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{generalized marginal interactions}
\keyword{marginal models}% __ONLY ONE__ keyword per line
