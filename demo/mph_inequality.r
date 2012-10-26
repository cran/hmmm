library(hmmm)

# 
# Example of MPH model subject to equality constraints: 
# 
# see Section 4.2, pg. 350,
# "Multinomial-Poisson homogeneous models for contingency tables",
# Lang, J.B.
# The Annals of Statistics, (2004)
# 
# Table 2 - 1999 statistics journals citation pattern counts (n_3=225)
# Citing x Cited statistics journals: JASA, BMCS, ANNS

y <- c(
104 , 24 , 65,
 76 ,146,  30,
 50 ,  9 ,166 )
y <- matrix(y,9,1)

# population matrix: 3 strata with 3 observations each
Z <- kronecker(diag(3),matrix(1,3,1))

# the 3rd stratum sample size is fixed
ZF <- kronecker(diag(3),matrix(1,3,1))[,3]

####################################################################
# Let P_ij the expected number of cross-citations, P_ij=P(A=i,B=j),
# where A = Citing journal and B = Cited journal. 
# The Gini concentrations of citations for each of the journals are: 
# G_i = sum_j=1_3 (P_ij/P_i+)^2  for i=1,2,3.
 
Gini.fct <- function(m) {
p <- diag(c(1/(Z%*%t(Z)%*%m)))%*%m 
G1 <- sum(diag(matrix(c(p[1]/(p[1]+p[2]+p[3]),0,0,
                  0,p[2]/(p[1]+p[2]+p[3]),0,
                  0,0,p[3]/(p[1]+p[2]+p[3])),3,3,byrow=TRUE)*
         matrix(c(p[1]/(p[1]+p[2]+p[3]),0,0,
                  0,p[2]/(p[1]+p[2]+p[3]),0,
                  0,0,p[3]/(p[1]+p[2]+p[3])),3,3,byrow=TRUE)))
G2 <- sum(diag(matrix(c(p[4]/(p[4]+p[5]+p[6]),0,0,
                  0,p[5]/(p[4]+p[5]+p[6]),0,
                  0,0,p[6]/(p[4]+p[5]+p[6])),3,3,byrow=TRUE)*
         matrix(c(p[4]/(p[4]+p[5]+p[6]),0,0,
                  0,p[5]/(p[4]+p[5]+p[6]),0,
                  0,0,p[6]/(p[4]+p[5]+p[6])),3,3,byrow=TRUE)))
G3 <- sum(diag(matrix(c(p[7]/(p[7]+p[8]+p[9]),0,0,
                  0,p[8]/(p[7]+p[8]+p[9]),0,
                  0,0,p[9]/(p[7]+p[8]+p[9])),3,3,byrow=TRUE)*
         matrix(c(p[7]/(p[7]+p[8]+p[9]),0,0,
                  0,p[8]/(p[7]+p[8]+p[9]),0,
                  0,0,p[9]/(p[7]+p[8]+p[9])),3,3,byrow=TRUE)))
c(G1,G2,G3)
}
####################################################################

## 
## --> h_2 = c(G1,G1)-c(G2,G3) = 0
## HYPOTHESIS: G1 = G2 = G3 
## 

h.fct<-function(m) {G<-Gini.fct(m)
c(G[1],G[1])-c(G[2],G[3])
}

mod_eq <- mphineq.fit(y,Z,ZF,h.fct=h.fct)

#summary(mod_eq)

print(mod_eq)


# OVERALL GOODNESS OF FIT: TEST of   Ho: h(m)=0 vs. Ha: not #Ho...


# Example of MPH model subject to inequality constraints: 


## --> d_2 = c(G1,G1)-c(G2,G3) >= 0
## HYPOTHESIS: G1 > G2, G1 > G3 

d.fct <-function(m) {G<-Gini.fct(m)
c(G[1],G[1])-c(G[2],G[3])
}

mod_ineq <- mphineq.fit(y,Z,ZF,d.fct=d.fct)

#summary(mod_ineq)

print(mod_ineq)

#OVERALL GOODNESS OF FIT: TEST of   Ho: h(m)=0 vs. Ha: not Ho...
#    Likelihood Ratio Stat (df= 0 ):  Gsq =  21.78688
#    Pearson's Score Stat  (df= 0 ):  Xsq =  20.81808



# HYPOTHESES TESTED:
# NB: testA --> H0=(mod_eq) vs H1=(mod_ineq model)
#     testB --> H0=(mod_ineq model) vs H1=(sat_mod model)
# sat_mod --> saturated model (Gsq=0)

m <- mod_ineq$m

chibar(m,Z,ZF,d.fct=d.fct,test0=mod_eq$Gsq-mod_ineq$Gsq,test1=mod_ineq$Gsq,repli=6000,lev=c(3,3))


# Example of MPH model subject to equality and inequality #constraints: 


## --> h = c(G1)-c(G2) = 0
## --> d = c(G3)-c(G2) >= 0
## HYPOTHESES: G1 = G2, G3 > G2 


h.fct2 <-function(m) {G<-Gini.fct(m)
G[1]-G[2]
}
 
d.fct2  <-function(m) {G<-Gini.fct(m)
G[3]-G[2]
}

mod_eq_ineq <- mphineq.fit(y,Z,ZF,h.fct=h.fct2,d.fct=d.fct2)

#summary(mod_eq_ineq)

print(mod_eq_ineq)




# HYPOTHESES TESTED:
# NB: testA --> H0=(mod_eq) vs H1=(mod_eq_ineq model)
#     testB --> H0=(mod_eq_ineq model) vs H1=(sat_mod model)
# sat_mod --> saturated model (Gsq=0)

m <- mod_eq_ineq$m

chibar(m,Z,ZF,h.fct=h.fct2,d.fct=d.fct2,test0=mod_eq$Gsq-mod_eq_ineq$Gsq,test1=mod_eq_ineq$Gsq,repli=2000,lev=c(3,3))



