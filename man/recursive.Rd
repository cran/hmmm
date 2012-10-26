\name{recursive}
\alias{recursive}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
recursive marginal interactions
}
\description{
A function to define logits of recursive (or nested) type.
}
\usage{
recursive(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{\dots}{As many inputs as there are variables in the multi-way table. Each input is a matrix of values -1,0,1 to define recursive logits or 0 for logits of different type}
}
\details{This function is used when logits of type "r" are used for at least one variable.
 An input argument for each categorical variable is necessary. Inputs are ordered according
to the order of the variables.

For a categorical variable with k categories, k-1 recursive logits can be defined using a matrix 
with k-1 rows and k columns. 
The rows of this matrix specify the categories whose probabilities constitute numerator and 
denominator of every recursive logit. Specifically,  
in every row, a value among -1,0,1 is associated to every category: value 1 (-1) corresponds to the category whose probability 
is cumulated at the numerator (denominator), 0 if the 
category is not involved.}
\value{
A zero-one matrix to be assigned to the \code{cocacontr} input argument in defining a model by `hmmm.model' when logits "r" 
are used for at least one variable in the multi-way table.
}
\references{Cazzaro M, Colombi R (2008) Modelling two way contingency tables with recursive logits and odds ratios. Statistical Methods and Applications, 17, 435-453.

Cazzaro M, Colombi R (2013)
Marginal nested interactions for contingency tables. Communications in Statistics - Theory and Methods, to appear.
}
\seealso{\code{\link{hmmm.model}}, \code{\link{create.XMAT}}}
\examples{
data(kentucky)
# 1 = injury 2 = restraint 3 = year
y<-getnames(kentucky,st=4)

marglist<-marg.list(c("m-m-l","m-l-l","r-l-l"),mflag="m")
R1<-matrix(c(1,1,1,-1,-1,
             0,0,0,1,-1,
             1,1,-1,0,0,
             1,-1,0,0,0),4,5,byrow=TRUE)
# logits of recursive (or nested) type for variable 1:
# log p(injury<=3)/p(injury>3); log p(injury=4)/p(injury=5);
# log p(injury<=2)/p(injury=3); log p(injury=1)/p(injury=2);
rec<-recursive(R1,0,0) # only variable 1 has recursive logits 

# additive effect of variables 2,3 on the recursive logits of variable 1
model<-hmmm.model(marg=marglist,lev=c(5,2,5),sel=c(34:49),cocacontr=rec)
mod<-hmmm.mlfit(y,model)
print(mod,printflag=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{logit}