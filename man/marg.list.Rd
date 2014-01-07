\name{marg.list}
\alias{marg.list}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ lists of marginal sets}
\description{
  An easy option to define the first input argument 
\code{marg} of the function `hmmm.model' which specifies the list of marginal sets of a hmm model.
}
\usage{
marg.list(all.m, sep = "-", mflag = "marg")
}
\arguments{
  \item{all.m}{A character vector with one element for every marginal set, see below}
  \item{sep}{ The separator used between logits type, default "-" }
  \item{mflag}{The symbol used to denote variables that are marginalized, default "marg"}
}
\details{\code{all.m} is a string indicating the logit types used to build the interactions 
in each marginal set. For each variable in the marginal set the corresponding logit symbol is inserted 
("b" baseline, "g" global, "c" continuation, "rc" 
reverse continuation, "r" recursive, "l" local). Symbols are separated 
by \code{sep} and the variables not included in the marginal set are denoted by \code{mflag}.
So, for example, "marg-g-c" indicates a marginal set involving variables 2, 3 with global and continuation logits respectively.}
\value{The list \code{marg} used as first input argument in `hmmm.model' -- see the function `hmmm.model'.}
\note{This function creates the complete list of the interactions that can be
defined in a marginal set. Therefore, it cannot be used to specify only the
interactions subject to inequality constraints. When inequalities are involved in the model, marginal sets and types of logits are declared as illustrated in the details of function `hmmm.model'.}
\seealso{\code{\link{hmmm.model}}}
\examples{
data(madsen)
marginals<-c("m-m-b-b","g-m-b-b","m-g-b-b","g-g-b-b")
margi<-marg.list(marginals,mflag="m")
names<-c("Inf","Sat","Co","Ho")
model<-hmmm.model(marg=margi,lev=c(3,3,2,4),names=names)
print(model)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
