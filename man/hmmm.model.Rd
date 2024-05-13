\name{hmmm.model}
\alias{hmmm.model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{define a hmm model}
\description{Function to define a hierarchical multinomial marginal model.}
\usage{
hmmm.model(marg = NULL, dismarg = 0, lev, cocacontr = NULL, strata = 1,
Z = NULL, ZF = Z, X = NULL, D = NULL, E = NULL, 
names = NULL, formula = NULL, sel = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{marg}{
 A list of the marginal sets and their  marginal interactions as
              described in Bartolucci et al. (2007). See below}
  \item{dismarg}{Similar to marg but used to define inequalities Kln(Am)>0. Default 0 if there are no
 inequalities}
\item{lev}{Number of categories of the variables}
  \item{cocacontr}{A list of zero-one matrices to build "r" logits
 created by the function `recursive'}
  \item{strata}{
Number of strata defined by the combination of the categories of the covariates}
  \item{Z}{
Zero-one matrix describing the strata}
  \item{ZF}{
Zero-one matrix for strata with fixed number of observations}
  \item{X}{
Design matrix for Cln(Mm)=Xbeta. Identity matrix if not declared.  
It can be defined later or changed only  by using the function `create.XMAT'}
  \item{D}{
If the matrix D is declared, the inequalities are expressed as  DKln(Am)>0.
           Useful for changing the sign of inequalities 
           or for selecting a subset of inequalities}
  \item{E}{If E is a matrix, then E defines the equality contrasts as
ECln(Mm)=0}
  \item{names}{A character vector whose elements are the names of the variables
}
  \item{formula}{
Formula of the reference log-linear model}
\item{sel}{
Vector reporting the positions of the interactions constrained to be zero}
}
\details{Variables are denoted by integers, the lower the number identifying the variable 
the faster its category subscript changes in the vectorized contingency table. Suppose that the variables are 1 and 2 
with categories k_1, k_2, the joint frequencies yij, where i=1,...,k_1, j=1,...,k_2, 
are arranged in a vector so that the subscript i changes faster than j. If \code{strata} is greater than one, the
vectorized contingency tables must be entered strata by strata. So that, for example, if the variables are distinguished in 
responses and covariates, the categories of the covariates determine the strata and the data are arranged 
 in such a way that the categories of the response variable changes faster than the categories of covariate. The names of the variables in \code{names} must be declared according to the order of the variables.  

The list \code{marg} of the marginal sets of a complete hierarchical marginal parameterization, together with the types of logits for the variables,
must be  created by the function `marg.list'. See the help of this function for more details.
If \code{marg} is not specified the multivariate logit model by Glonek and McCullagh (1995) 
with interactions of type local is used. The list  \code{marg} is used to create the link function 
Cln(Mm) and its derivative (m is the vector of expected frequencies).

If the model is defined in the form Cln(Mm) = Xbeta, the matrix X 
has to be declared (see the function `create.XMAT'). If there are only nullity constraints on parameters, the model is in the form ECln(Mm)=0 and X is ignored. In such a case, E can be declared as matrix or it is automatically constructed if \code{sel} is declared. If \code{sel} is not NULL, then the model is defined under equality constraints, i.e. ECln(Mm)=0.
When \code{X}, \code{E} and \code{sel} are left at default level, a saturated model is defined.  
        

For models with inequality constraints on marginal parameters, the input argument \code{dismarg} is declared as a list whose components are of type: \code{list(marg=c(1,2),int=list(c(1),c(1,2)),} \code{types=c("g","l" ))}, with elements
  \code{marg}: the marginal set, \code{int}: the list of the interaction subject to inequality constraint, and \code{types}: the logit 
used for every variable ("g"=global, "l"=local, "c"=continuation, "rc"=reverse continuation, 
"r"=recursive, "b"=baseline, "marg" is assigned  to each variable not belonging 
to the marginal set).  This list is used to create the link function 
Cln(Mm) and its derivative for the inequality constraints. 

 The matrix Z is of dimension c x s, where c is the number of counts 
               and s is the number of strata or populations. 
                   Thus, the rows correspond to the number of observations  and 
                   the columns correspond to the strata.  A 1 in 
                   row i and column j means that the ith count comes 
                   from the jth stratum.  Note that Z has exactly 
                   one 1 in each row, and at least one 1 in each 
                   column.  When the population matrix Z 
                   is a column vector of 1 indicates that all 
                   the counts come from the same and only stratum.
                    For hmm models, it is assumed that all 
                     the strata have  the same number of response levels. 
                     If Z is not given, a population Z matrix corresponding 
                     to data entered by \code{strata} is defined and ZF=Z. 
                     For non-zero ZF, the columns 
                   are a subset of the columns in Z.
                   If the jth column of Z is included in ZF, then the 
                   sample size of the jth stratum is considered fixed, otherwise
                    if the jth column of Z is NOT included in ZF, the 
                   jth stratum sample size is taken to be a realization 
                   of a Poisson random variable.  As ZF=Z the sample size in every stratum 
                 is fixed; this 
                   is the (product-)multinomial setting.  
                   
The formula  of the reference log-linear model must be defined using the names of the variables declared in \code{names}, for example
\code{names<-c("A","B","C","D")}, \code{formula=~A*C*D+B*C*D+A:B}. The interactions not involved in \code{formula} cannot be further constrained in 
the marginal model. The default \code{formula = NULL} indicates the saturated log-linear model as reference model.
The likelihood function of the reference model is maximized by `hmmm.mlfit' under the constraints ECln(Mm)=0 on the marginal parameters.

The arguments \code{dismarg} and \code{formula} can be used only if \code{strata=1}.}
\value{
  An object of the class \code{hmmmmod}; it describes a marginal model that can be estimated by `hmmm.mlfit'.
}
\references{Bartolucci F, Colombi R, Forcina A (2007) An extended class of marginal link functions for modelling
 contingency tables by equality and inequality constraints. Statistica Sinica, 17, 691-711.

Bergsma WP, Rudas T (2002) Marginal models for categorical data. 
The Annals of Statistics, 30, 140-159.

Cazzaro M, Colombi R (2009) Multinomial-Poisson models subject to inequality 
constraints. Statistical Modelling, 9(3), 215-233.

Colombi R, Giordano S, Cazzaro M (2014) hmmm: An R Package for hierarchical multinomial marginal models. Journal of Statistical Software, 59(11), 1-25, URL http://www.jstatsoft.org/v59/i11/.

Glonek GFV, McCullagh P (1995) Multivariate logistic models for
contingency tables. Journal of the Royal Statistical Society, B, 57, 533-546.
}
\seealso{\code{\link{hmmm.model.X}}, \code{\link{create.XMAT}}, \code{\link{summary.hmmmmod}}, \code{\link{print.hmmmmod}}, 
\code{\link{marg.list}}, \code{\link{recursive}}, \code{\link{hmmm.mlfit}}}
\examples{data(madsen)
# 1 = Influence; 2 = Satisfaction;  3 = Contact; 4 = Housing 
names<-c("Inf","Sat","Co","Ho")
y<-getnames(madsen,st=6)

# hmm model -- marginal sets: {3,4} {1,3,4} {2,3,4} {1,2,3,4}
margi<-c("m-m-l-l","l-m-l-l","m-l-l-l","l-l-l-l")
marginals<-marg.list(margi,mflag="m")
model<-hmmm.model(marg=marginals,lev=c(3,3,2,4),names=names)
summary(model)

# hmm model with equality constraints
# independencies 1_||_4|3 and 2_||_3|4 impose equality constraints 
sel<-c(12:23,26:27,34:39) # positions of the zero-constrained interactions 
model_eq<-hmmm.model(marg=marginals,lev=c(3,3,2,4),sel=sel,names=names)
summary(model_eq)

# hmm model with inequality constraints
# the distribution of 1 given 4 is stochastically decreasing wrt the categories of 3;
# the distribution of 2 given 3 is stochastically decreasing wrt the categories of 4:
marg134ineq<-list(marg=c(1,3,4),int=list(c(1,3)),types=c("l","marg","l","l"))
marg234ineq<-list(marg=c(2,3,4),int=list(c(2,4)),types=c("marg","l","l","l"))
ineq<-list(marg134ineq,marg234ineq)
model_ineq<-hmmm.model(marg=marginals,lev=c(3,3,2,4),dismarg=ineq,D=diag(-1,8),names=names)
summary(model_ineq)
# The argument D is used to turn the 8 inequalities from 
# non-negative (default) into non-positive constraints}
\keyword{ models}
