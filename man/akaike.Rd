\name{akaike}
\alias{akaike}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ akaike criterium }
\description{
Compute AIC value for a list of hmm and hidden Markov models
}
\usage{
akaike(..., LRTEST = FALSE, ORDERED = FALSE, NAMES = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{objects  created by `hmmm.mlfit' or 
`hidden.emfit'
}
  \item{LRTEST}{
If TRUE, the first model must include all the others models as special cases.
For every model, the  likelihood ratio statistic test with respect to the first model is computed
}
  \item{ORDERED}{
If TRUE, in the output the models are ordered according to the  Akaike criterium starting from the lowest AIC value
}
  \item{NAMES}{
Optional character vector with the names of the models. If it is NULL (the default) model names are created as model1, model2..... 
}
}
\details{
The models in input must be at least two objects of the classes \code{hmmmfit} or \code{hidden}.
}
\value{A matrix with row names given by \code{NAMES} and column names describing the 
output for every model (position of the model in the input list \code{#model}, 
the loglikelihood function \code{loglik}, the number of parameters \code{npar}, 
the number of constraints of the model \code{dfmodel}, likelihood ratio test \code{LRTEST}, degrees of freedom \code{dftest}, \code{PVALUE}, \code{AIC}, \code{DELTAAIC}). The \code{DELTAAIC} is the difference between the AIC value of every model and the lowest AIC value.
}
\references{Konishi S, Kitagawa G (2008)
Information criteria and statistical modeling. Springer.
}

\examples{data(madsen)
# 1 = Influence; 2 = Satisfaction;  3 = Contact; 4 = Housing 
names<-c("Inf","Sat","Co","Ho")
y<-getnames(madsen,st=6)

margin <- marg.list(c("marg-marg-l-l", "g-marg-l-l", "marg-g-l-l", "g-g-l-l"))

# additive effect of 3 and 4 on logits of 1 in marginal  
# distribution {1, 3, 4}, conditional independence 2_||_3|4 
modelA <- hmmm.model(marg = margin, lev = c(3, 3, 2, 4), names = names)
modA <- hmmm.mlfit(y, modelA)
modA

# additive effect of 3 and 4 on logits of 1 in marginal
# distributions {1, 3, 4} and {2, 3, 4}
modelB <- hmmm.model(marg = margin, lev = c(3, 3, 2, 4),
names = names, sel = c(18:23, 34:39))
modB <- hmmm.mlfit(y, modelB)
modB

# 1 and 2 do not depend on the levels of 3 and 4
modelC <- hmmm.model(marg = margin, lev = c(3, 3, 2, 4),
names = names, sel = c(18:23, 34:39, 44:71))
modC <- hmmm.mlfit(y, modelC)
modC

akaike(modB, modA, modC, ORDERED = TRUE, NAMES = c("modB", "modA", "modC"))
akaike(modB, modA, modC, LRTEST = TRUE, NAMES = c("modB", "modA", "modC"))}

\keyword{ models}% 

