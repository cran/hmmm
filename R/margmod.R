
block.fct <-  function(...)
{
        #  Author: Joseph B. Lang 
        #          Dept of Stat and Act Sci 
        #          Univ of Iowa  (6/16/99, last update: 3/30/04)
        #   
        #  Direct sum function.  Creates a block diagonal matrix.
        #   
        n <- nargs()
        all.m <- list(...)
        rows <- 0
        cols <- 0
        for(i in 1:n) {
                rows <- rows + nrow(all.m[[i]])
                cols <- cols + ncol(all.m[[i]])
        }
        m <- matrix(0, rows, cols)
        r1 <- 1
        r2 <- 0
        c1 <- 1
        c2 <- 0
        for(i in 1:n) {
                r2 <- r2 + nrow(all.m[[i]])
                c2 <- c2 + ncol(all.m[[i]])
                m[r1:r2, c1:c2] <- all.m[[i]]
                r1 <- r2 + 1
                c1 <- c2 + 1
        }
        m

}


pop <- function(npop,nlev) {
 # 
 # Creates a population (Z) matrix corresponding 
 # to data entered by strata. It is assumed that all 
 # the strata have  the same number of response levels.
 #
 # Author: Joseph B. Lang
 # Date  : June 5, 2002  (last updated: March 30, 2004)
 # 
 # Input:
 #    npop = number of populations (aka strata)
 #    nlev = number of response levels per stratum
 # Output:
 #    Z = population matrix of dimension (npop*nlev X npop).
 #
   kronecker(diag(npop),matrix(1,nlev,1))
  }


rmult <- function(N = 1, n, p) 
{
#  Author: Joseph B. Lang 
#          Dept of Stat and Act Sci 
#          Univ of Iowa  (6/16/99, last update: 3/30/04)
#      
#  Create a matrix of multinomial realizations.
#  There will be N columns, each column contains a realization of 
#      a multinomial(n,p) random vector.
#  Note:  Components in p must add up to 1.0. (e.g. p = c(0.2,0.8))
#
        K <- length(p)
        x <- matrix(0, K, N)
        su <- x[1,  ]
        sup <- 0
        for(i in 1:(K - 1)) {
                x[i, n > su] <- rbinom(length(su[n > su]), 
                n - su[n >su], p[i]/(1 - sup))
                x[i, n <= su] <- 0
                su <- su + x[i,  ]
                sup <- sup + p[i]
        }
        x[K,  ] <- n - su
        x
}

 #  Author: Roberto Colombi
        #          Dept IGIF
        #          Univ of Bergamo ( last update: 25/07/06)

   #definition of HMM models


cocamod<-function(model,Z){
nmarg<-length(model)-3
strata<-c(model$strata)
levs<-c(model$livelli)
if(is.null(model$cocacontr)){model$cocacontr<-as.list(rep(0,length(levs)))}
	
C<-list()
M<-list()

for (mi in 1:nmarg){

marginal<-model[[mi]]
margindex<-c(marginal$marg)
margset<-marginal$int
types<-c(marginal$types)


nint<-length(margset)
CM<-list()
MM<-list()


for (ii in 1:nint){
margint<-c(margset[[ii]])
matrici<-cocamat(margint,margindex,types,levs,model$cocacontr)
CM[[ii]]<-matrici$CMAT
MM[[ii]]<-matrici$MMAT


}

M[[mi]]<-MM[[1]]
C[[mi]]<-CM[[1]]

if(nint>1){
for(is in 2:nint){
M[[mi]]<-rbind(M[[mi]],MM[[is]])
C[[mi]]<-block.fct(C[[mi]],CM[[is]])
}
}
remove(list=c("MM","CM"))
}
MG<-M[[1]]
CG<-C[[1]]

if(nmarg>1){
for(i in 2:nmarg){
MG<-rbind(MG,M[[i]])
CG<-block.fct(CG,C[[i]])
}
}
I<-which(Z==1,arr.ind=TRUE)
IM<-diag(1,dim(Z)[1])
PZZ<-IM[I[,1],]

mmat<-kronecker(diag(1,strata),MG)%*%PZZ
list(CMAT=kronecker(diag(1,strata),CG),MMAT=mmat)

}
cocamat<-function(margint,margindex,types,levs,rmat)
{

###################


for (i in 1:length(levs)){
if(types[i]=="b"){
a<-gl(levs[i],1)
a<-t(contr.treatment(a))
rownames(a)<-NULL
colnames(a)<-NULL
m<-matrix(0,levs[i]-1,levs[i])
m[1:(levs[i]-1),1]<-1


m<-rbind(m,a)
rmat[[i]]<-m

}
}


########################



m<-list()
c<-list()
for(i in 1:length(levs)){
c[[i]]<-1
if (types[i]=="marg" ){ m[[i]]<-matrix(1,1,levs[i])}
else {
if ((types[i]=="r")|(types[i]=="b") ){m[[i]]<-matrix(c(rmat[[i]][1,]),1,levs[i],byrow=TRUE)}
else{m[[i]]<-matrix(c(1,rep(0,levs[i]-1)),1,levs[i])}
}
}
for(i in margint){
c[[i]]<-cbind(-diag(1,levs[i]-1),diag(1,levs[i]-1))
mid<-diag(1,levs[i])
mupper<-mid
mlower<-mid
mupper[upper.tri(mupper,diag=TRUE)]<-1
mlower[lower.tri(mlower,diag=TRUE)]<-1
if  (types[i]=="l") {m[[i]]<- t(cbind(mid[,-levs[i]],mid[,-1]))}    
if  (types[i]=="c") {m[[i]] <-rbind(t(diag(1,levs[i])[,-levs[i]]),t(mlower[,-1]))}     
if  (types[i]=="rc") {m[[i]]<- rbind(t(mupper[,-levs[i]]),t(diag(1,levs[i])[,-1]))} 
 if  (types[i]=="g")									  
 {m[[i]]<- rbind(t(mupper[,-levs[i]]),t(mlower[,-1])) }
 if  ((types[i]=="r")|(types[i]=="b")) {m[[i]]<-rmat[[i]]}

}
M<-m[[1]]
C<-c[[1]]
if(length(levs)>1){

for(i in 2:length(levs)){
M<-kronecker(m[[i]],M)
C<-kronecker(c[[i]],C)
}
}
matrici<-list(CMAT=C,MMAT=M)
}


#design matrix of the  working log-linear model used by the main function and by the chibar functions

LDMatrix<-function(level,formula,names=NULL){
if (is.null(names)) {names<-LETTERS[1:length(level)]
}
clev<-cumprod(level)
totlev<-prod(level)
C<-gl(level[1],1,totlev)
for (i in 2:length(level)) {
c<-gl(level[i],clev[i-1],totlev)
C<-cbind(C,c)
}

colnames(C)<-names
C<-as.data.frame(C)
C<-data.frame(lapply(C,as.factor))
C<-model.matrix(as.formula(formula),C)
C<-C[,-1]
matrici<-list(IMAT=solve(t(C)%*%C)%*%t(C),DMAT=C)
}


#########
cocadise<-function(Z=NULL,names=NULL,formula=NULL,lev)
{
if (is.null(formula))
{
ncz<-dim(Z)[2]
zi<-Z[,1]
zi<-zi[zi>0]
DD<-diag(c(zi))[,2:length(zi)]
zi<-zi[-1]
DI<-cbind(-zi,diag(c(zi)))
if( ncz>1){
for(i in 2:ncz){
zi<-Z[,i]
zi<-zi[zi>0]
DMATI<-diag(c(zi))[,2:length(zi)]
zi<-zi[-1]
CMATI<-cbind(-zi,diag(c(zi)))
DD<-block.fct(DD,DMATI)
DI<-block.fct(DI,CMATI)
}
}
matrici<-list(IMAT=DI,DMAT=DD)
}
else
LDMatrix(lev,formula,names)

}
#  Author: Roberto Colombi

#  Author: Roberto Colombi
        #          Dept IGIF
        #          Univ of Bergamo ( last update: 25/07/06)

#function of constraints and their derivatives for HM models

make.h.fct<-function(models,E=TRUE)
{
if(all(E)){
function(m){models$matrici$CMAT%*%log(models$matrici$MMAT%*%m)}
}
else{

function(m){models$matrici$E%*%models$matrici$CMAT%*%log(models$matrici$MMAT%*%m)}
}
}

make.d.fct<-function(dismod,D=TRUE)
{
if(is.null(D)){
function(m){dismod$matrici$CMAT%*%log(dismod$matrici$MMAT%*%m)}
}
else{

function(m){dismod$matrici$D%*%
dismod$matrici$CMAT%*%log(dismod$matrici$MMAT%*%m)}
}
}



make.derht.fct<-function(models,E=TRUE)
{
if(all(E)){

function(m){
t(models$matrici$CMAT%*%diag(1/c(models$matrici$MMAT%*%m))%*%models$matrici$MMAT)}
}
else{

function(m){t(models$matrici$E%*%models$matrici$CMAT%*%diag(1/c(models$matrici$MMAT%*%m))%*%models$matrici$MMAT)}
}
}

make.derdt.fct<-function(models,D=TRUE)
{
if(is.null(D)){

function(m){
t(models$matrici$CMAT%*%diag(1/c(models$matrici$MMAT%*%m))
%*%models$matrici$MMAT)}
}
else{

function(m){t(models$matrici$D%*%
models$matrici$CMAT%*%diag(1/c(models$matrici$MMAT%*%m))%*%models$matrici$MMAT)}
}
}

make.L.fct<-function(models,E=TRUE)
{
if(E==TRUE){
function(m){models$matrici$CMAT%*%log(models$matrici$MMAT%*%m)}
}
else{

function(m){t(create.U(models$matrici$X))%*%models$matrici$CMAT%*%log(models$matrici$MMAT%*%m)}
}
}

make.derLt.fct<-function(models,E=TRUE)
{

if(E==TRUE){
function(m){
t(models$matrici$CMAT%*%diag(1/c(models$matrici$MMAT%*%m))%*%models$matrici$MMAT)}
}
else{
function(m){t(t(create.U(models$matrici$X))%*%models$matrici$CMAT%*%
diag(1/c(models$matrici$MMAT%*%m))%*%models$matrici$MMAT)}
}
}

#  Author: Roberto Colombi
        #          Dept IGIF
        #          Univ of Bergamo ( last update: 25/07/06)


# chibar pvalue

chibar<-function(m,Z,ZF,d.fct=0,h.fct=0,test0=0,test1=0,repli=0,
derdt.fct=0,derht.fct=0,formula=NULL,names=NULL,lev){
Zlist<-cocadise(Z,formula=formula,lev=lev,names=names)
  p<-m*c(1/Z%*%t(Z)%*%m)
   Dm <- diag(c(m))
#-((ZF*c(m))%*%t(ZF*c(p)))
   Hmat<- t(Zlist$DMAT)%*%( diag(c(m))-((ZF*c(m))%*%t(ZF*c(p))))%*%Zlist$DMAT
 if (is.function(derdt.fct)==F)
   {
     DH <- num.deriv.fct(d.fct,m) 
   }
   else {
     DH <- derdt.fct(m)
 

   }
   DH<-	t(Zlist$DMAT)%*%Dm%*%DH
   	  
	 D<-t(DH) 
	 if (is.function(h.fct)==TRUE)
  {
   
 if (is.function(derht.fct)==F) {
     H <- num.deriv.fct(h.fct,m) 
   }
   else {
     H <- derht.fct(m)
   }
   	H<-t(Zlist$DMAT)%*%Dm%*%H
	E<-t(H)
	#X<-Null(t(E))
       X<-create.U(t(E))
     D<-D%*%X
   }

else{X<-diag(1,dim(D)[2])}

Hmat<-t(X)%*%Hmat%*%X
l<-dim(Hmat)[1]
 #prov
mi<-solve(chol(Hmat))

DD<-D%*%mi
#ZZ<-matrix(rnorm(repli*l,0,1),l,repli)

ur<-qr(D)$rank
qq<-qq0<-w<-matrix(0,ur+1,1)


#for(i in 1:repli){
#z<-matrix(ZZ[,i],l,1)
Z<-matrix(rnorm(l*repli,0,1),l,repli)
qq<-qq0<-matrix(0,ur+1,1)
pesiw<-function(z){
#solve.QP(Dmat, dvec, Amat, bvec, meq=0, factorized=FALSE)
pw<-matrix(0,ur+1,1)
ddl<-NULL

ddl<-
solve.QP(Dmat=diag(1,l),dvec= z, Amat=t(DD), meq=0, factorized=FALSE)

if (is.null(ddl)){next}

#dd<-ddl$solution
#dd0<-ddl$unconstrainted.solution
d<-NULL
d<-qr(matrix(D[!(DD%*%ddl$solution>1e-16),],ncol=l,byrow=TRUE))$rank


pw[d+1]<-pw[d+1]+1
if (all(!is.na(pw))){
qq[d+1]<<-qq[d+1]+t(ddl$solution-z)%*%(ddl$solution-z)
qq0[d+1]<<-qq0[d+1]+ur-t(ddl$solution-z)%*%(ddl$solution-z)}
rm(ddl)
if (any(is.na(pw))){
pw<-matrix(0,ur+1,1)
print( "failed")}

pw
}

w<-apply(Z,MARGIN=2,FUN=pesiw)
w<-apply(w,1,sum)
qq0<-qq0/w
qq<-qq/w
w<-w/sum(w)

gdl0<-matrix(0,ur+1,1)
if( (test0 >0)){
q1<-matrix(0,ur+1,1)}
else{
q1<-matrix(1,ur+1,1)}
if( (test0 >0)){
for(i in 1:(ur+1)){
gdl0[i]<-ur+1-i
if ((w[i] >0)&(gdl0[i]==0)){
q1[i]=0}
if ((w[i] >0)&(gdl0[i]>0)){
q1[i]<-1-pchisq(test0,gdl0[i])}
}
}

p0<-t(w)%*%q1
gdl1<-matrix(0,ur+1,1)
if( (test1 >0)){
q1<-matrix(0,ur+1,1)}
else{
q1<-matrix(1,ur+1,1)}
if( (test1 >0)){
for(i in 1:(ur+1)){
gdl1[i]<-i-1
if ((w[i] >0)&(gdl1[i]==0)){
q1[i]=0}
if ((w[i] >0)&(gdl1[i]>0)){
q1[i]<-1-pchisq(test1,gdl1[i])}
}
}
p1<-t(w)%*%q1

lista<-list(testA=test0,pvalA=p0,testB=test1,pvalB=p1,pesi=cbind(w,gdl0,gdl1,qq0,qq))
class(lista)<-"hmmmchibar"
lista

}
#  Author: Roberto Colombi
        #          Dept IGIF
        #          Univ of Bergamo ( last update: 25/07/06)
#------------------------------



print.hmmmchibar<-function(x,...){chibar.summary(x)}
summary.hmmmchibar<-function(object,plotflag=1,step=0.01,lsup=0,...){chibar.summary(object,plotflag=plotflag,step=step,lsup=lsup)}

chibar.summary<-function(P,plotflag=0,step=0.01,lsup=0){
results<-
matrix(c(P$testA,P$testB,P$pvalA,P$pvalB),2,2)
rownames(results)<-c("testA","testB")
colnames(results)<-c("test","pvalue")
colnames(P$pesi)<-c("weights","df A","df B","sim df A","sim df B")
rownames(P$pesi)<-as.character(P$pesi[,3])
cat("\n chibar simulated pvalues\n")
cat("\n")
print(results)
cat("\n")
if(plotflag>0){
cat("\n simulated weights of the chibar-distribution\n")
cat("\n")
print(P$pesi)}
if(plotflag>1){

ur<-length(P$pesi[,1])-1
if(lsup==0){lsup=2*ur}
Z<-seq(from=0,to=lsup,by=step)
Z<-matrix(Z,length(Z),1)
#print(Z)

Fchibar<-function(x){



test0<-x
test1<-x


gdl0<-matrix(0,ur+1,1)
if( (test0 >0)){
q1<-matrix(0,ur+1,1)}
else{
q1<-matrix(1,ur+1,1)}
if( (test0 >0)){
for(i in 1:(ur+1)){
gdl0[i]<-ur+1-i
if ((P$pesi[,1][i] >0)&(gdl0[i]==0)){
q1[i]=0}
if ((P$pesi[,1][i] >0)&(gdl0[i]>0)){
q1[i]<-1-pchisq(test0,gdl0[i])}
}
}
p0<-t(P$pesi[,1])%*%q1
gdl1<-matrix(0,ur+1,1)
if( (test1 >0)){
q1<-matrix(0,ur+1,1)}
else{
q1<-matrix(1,ur+1,1)}
if( (test1 >0)){
for(i in 1:(ur+1)){
gdl1[i]<-i-1
if ((P$pesi[,1][i] >0)&(gdl1[i]==0)){
q1[i]=0}
if ((P$pesi[,1][i] >0)&(gdl1[i]>0)){
q1[i]<-1-pchisq(test1,gdl1[i])}
}
}
p1<-t(P$pesi[,1])%*%q1
c(p0,p1)
}
F<-apply(Z,1,FUN=Fchibar)
F<-t(F)
if(plotflag>2){
matplot(Z,F,type="l",ylim=c(0,1))
abline(P$pvalA,0,col="black")
abline(P$pvalB,0,col="red")
if(P$testA<lsup){
abline(v=P$testA,col="black")}
if(P$testB<lsup){
abline(v=P$testB,col="red")}

}
F<-cbind(Z,F)
colnames(F)<-c("x","FA","FB")
F
}
}




#  Author: Roberto Colombi
        #          Dept IGIF
        #          Univ of Bergamo ( last update: 25/07/06)

bcf.interactions<-function(marglist)
{
nmarg<-length(marglist)
cumintset<-NULL
for(i in 1:nmarg){
if (length(marglist[[i]]$marg)==1){
allintset<-marglist[[i]]$marg}

else{
allintset<-bar.col.for(marglist[[i]]$marg)}

marglist[[i]]$int<-setdiff(allintset,cumintset)
cumintset<-union(cumintset,marglist[[i]]$int)
}
marglist
}

#rm(list=ls())

bar.col.for<-function(marg)
{
p<-length(marg)
cifre<-p
bincol<-function(n){
if (n<=1) v=n
else if (n%%2==0) v=c(Recall(n/2),0)
else v=c(Recall((n-1)/2),1)

#n=length(v);
#print(cifre)

#if (n<cifre)  c(rep(0,cifre-n),v)
#else v
}

b<-matrix(0,1,p)
s<-2^p-1
for (i in 1:s){

bb<-bincol(i)
n=length(bb)
if (n<p)  bb<-c(rep(0,cifre-n),bb)
b<-rbind(b,matrix(bb,1,p))
}


b<-b[-1,]

b<-b[order(rowSums(b),b%*%matrix(1:p,p,1)),]

int<-list()
for(ii in c(1:dim(b)[1])){
int[[ii]]<-marg[which(b[ii,]==1)]}
int

}
#marg1<-list(marg=c(2),types=c("marg","c"))
#marg12<-list(marg=c(1,2),types=c("g","l"))
#marg123<-list(marg=c(1,2,3),types=c("g","l","l"))
#marginali<-list(marg1,marg12,marg123)
#newmmmm<-bcf.interactions(marginali)
#a<-bar.col.for(c(2,4,))
#b<-bar.col.for(c(6,10))

#  Author: Roberto Colombi
        #          Dept IGIF
        #          Univ of Bergamo ( last update: 25/07/06)
marg.list<-function(all.m,sep="-",mflag="marg") {
n <- length(all.m)


marglist<-list()
for(i in 1:n) {
ca<-all.m[i]

ca<-unlist(strsplit(ca,split=sep))
ca[ca==mflag]<-"marg"
ci<-which(ca!="marg")

marglist[[i]]<-list(marg=ci,types=ca)
}
marglist
}

#  Author: Roberto Colombi
        #          Dept IGIF
        #          Univ of Bergamo ( last update: 25/07/06)


loglin.model<-function(lev,int=NULL,strata=1,dismarg=0,type="b",D=TRUE,
c.gen=TRUE,printflag=FALSE,names=NULL,formula=NULL){
if(is.null(formula)){

cocacontr=NULL
if(type=="b"){
cocacontr<-list()
for (i in 1:length(lev)){
a<-gl(lev[i],1)
a<-t(contr.treatment(a))
rownames(a)<-NULL
colnames(a)<-NULL
m<-matrix(0,lev[i]-1,lev[i])
m[1:(lev[i]-1),1]<-1


m<-rbind(m,a)
cocacontr[[i]]<-m

}
}


MARG<-c(1:length(lev))
all.int<-bar.col.for(1:length(lev))
s<-length(all.int)
type.temporary<-rep(type,length(lev))
marg<-list()
marg<-list(marg=MARG,int=all.int,types=type.temporary)
model<-hmmm.model(marg=list(marg),lev=lev,cocacontr=cocacontr)
dscr<-hmmm.model.summary(model,printflag=printflag)
XX<-diag(1,prod(lev)-1)
if(!c.gen==TRUE&!is.null(int)){
keep<-list()

for (i in 1:s){

for (ii in 1:length(int)){
a<-intersect(all.int[[i]],int[[ii]])
if(setequal(a,int[[ii]])){keep[[length(keep)+1]]<-all.int[[i]]



}
}


}
int<-setdiff(all.int,keep)
}

if(!is.null(int)){




sel<-0
for(i in 1:length(int)){
if(length(int[[i]]) >=2){
included.index<-bar.col.for(int[[i]])}
else{included.index<-int[[i]]}
for(ii in 1:length(included.index)){
inint<-paste(c(included.index[[ii]]),collapse="")
inint<-as.numeric(dscr[which(dscr[,1]==inint,arr.ind=TRUE),][5:6])
sel<-c(sel,c(inint[1]:inint[2]))}
}
XX<-unique(XX[,sel],MARGIN=2)

}


if(printflag==TRUE){

i<-rowSums(XX)
i<-i[as.numeric(dscr[,6])]
print("included interactions")
print(dscr[which(i==1),],quote=FALSE)
print("exluded interactions")
print(dscr[which(i==0),],quote=FALSE)
}
XX=kronecker(diag(1,strata),XX)
mod.loglin<-hmmm.model(marg=list(marg),dismarg=dismarg,lev=lev,strata=strata,X=XX,cocacontr=cocacontr,names=names)
}
else{
marginali<-paste(c(rep(type,length(lev))),collapse="-")
marginali<-marg.list(marginali)

hmmm.model(marg=marginali
,lev=lev,
names=names,formula=formula)
}

}
 








#  Author: Roberto Colombi
        #          Dept IGIF
        #          Univ of Bergamo ( last update: 25/07/06)
hmmm.model<-
function(marg=NULL,dismarg=0,lev,cocacontr=NULL,strata=1,Z=NULL,ZF=Z,X=NULL,D=NULL,
E=NULL,names=NULL,formula=NULL,sel=NULL)
{if(strata>1){formula<-NULL
dismarg<-0}

if(!is.null(sel)){ E<-t(diag(1,(prod(lev)-1))[,sel])}
if(is.matrix(E)){X<-0}
if((is.null(X))&!is.matrix(E)){ X<-diag(1,(prod(lev)-1))  }   
if (is.null(Z)){
 Z<-pop(strata,prod(lev))}
if (is.null(ZF)){
 ZF<-Z}
if(is.null(marg)){
MARG<-bar.col.for(1:length(lev))
s<-length(MARG)
marg<-list()
for (i in 1:s){
type.temporary<-rep("marg",length(lev))
type.temporary[MARG[[i]]]<-"l"
marg[[i]]<-list(marg=MARG[[i]],types=type.temporary)
}
}


s<-length(marg)

for (i in 1:s){
if (is.null(marg[[i]]$int)){
marg<-bcf.interactions(marg)
break
}
else next
}
marginali<-c(marg,list(livelli=lev,cocacontr=cocacontr,strata=strata))
#print(str(marginali))

model<-list(modello=marginali,matrici=cocamod(marginali,Z),formula=formula)
model$matrici$Z<-Z
model$matrici$ZF<-ZF

model$matrici$X<-X
model$matrici$E<-0
if (is.matrix(X)&is.null(E)) {
E<-FALSE
#t(create.U(models$matrici$X))
model$matrici$E<-t(create.U(model$matrici$X))
}
if(is.matrix(E)){
model$matrici$E<-E
model$matrici$X<-create.U(t(E))
E<-FALSE }

if ( sum(abs(model$matrici$E)) == 0){
model$functions$h.fct<-0
model$functions$derht.fct<-0}
else{
model$functions$h.fct<-make.h.fct(model,E)
 model$functions$derht.fct<-make.derht.fct(model,E)
 }
model$functions$L.fct<-make.L.fct(model,TRUE)
model$functions$derLt.fct<-make.derLt.fct(model,TRUE)




if(is.list(dismarg)){
if(!is.list(dismarg[[1]])){dismarg<-list(dismarg)}

model$dismod<-c(dismarg,list(livelli=lev,cocacontr=cocacontr,strata=strata))

model$dismod<-list(modello=model$dismod,matrici=cocamod(model$dismod,Z))
model$dismod$matrici$D<-D
model$functions$d.fct=make.d.fct(model$dismod,D)
model$functions$derdt.fct=make.derdt.fct(model$dismod,D)
}
model$names<-names

class(model)<-"hmmmmod"
model
}
#  Author: Roberto Colombi
        #          Dept IGIF
        #          Univ of Bergamo ( last update: 25/07/06)

hmmm.mlfit<-
function(y,model,noineq=TRUE,maxit=1000,norm.diff.conv=1e-5,norm.score.conv=1e-5,
                   y.eps=0,chscore.criterion=2,
                   m.initial=y,mup=1,step=1){





if(noineq){
a <- mphineq.fit(y,Z=model$matrici$Z,ZF=model$matrici$ZF,E=model$matrici$E,
L.fct=model$functions$L.fct,
derLt.fct=model$functions$derLt.fct,
h.fct=model$functions$h.fct,derht.fct=model$functions$derht.fct,
X=model$matrici$X,formula=model$formula,names=model$names,lev=model$modello$livelli,
maxiter=maxit,norm.diff.conv=norm.diff.conv,norm.score.conv=norm.score.conv,
y.eps=y.eps,
chscore.criterion=chscore.criterion,m.initial=m.initial,mup=mup,
step=step)


}
else{
a<-mphineq.fit(y,Z=model$matrici$Z,ZF=model$matrici$ZF,E=model$matrici$E,
L.fct=model$functions$L.fct,
derLt.fct=model$functions$derLt.fct,d.fct=model$functions$d.fct,
h.fct=model$functions$h.fct,derht.fct=model$functions$derht.fct,
derdt.fct=model$functions$derdt.fct,
X=model$matrici$X,formula=model$formula,names=model$names,lev=model$modello$livelli,
maxiter=maxit,norm.diff.conv=norm.diff.conv,norm.score.conv=norm.score.conv,
y.eps=y.eps,
chscore.criterion=chscore.criterion,m.initial=m.initial,mup=mup,step=step)

}
a$model<-model
class(a)<-"hmmmfit"
a
}
#  Author: Roberto Colombi
        #          Dept IGIF
        #          Univ of Bergamo ( last update: 25/07/06)

hmmm.chibar<-function(nullfit,disfit,satfit,repli=6000){
model<-disfit$model
P<-chibar(m=nullfit$m,Z=model$matrici$Z,ZF=model$matrici$ZF,
d.fct=model$functions$d.fct,derdt.fct=model$functions$derdt.fct,
h.fct<-model$functions$h.fct,derht.fct=model$functions$derht.fct,
test0=c(nullfit$Gsq)-c(disfit$Gsq),test1=c(disfit$Gsq)-c(satfit$Gsq),repli=repli,
formula=model$formula,names=model$names,lev=model$modello$livelli)
class(P)<-"hmmmchibar"
P
}



#  Author: Roberto Colombi
        #          Dept IGIF
        #          Univ of Bergamo ( last update: 25/07/06)


#definition of HMM models
#
print.hmmmmod<-function(x,...){
hmmm.model.summary(x,printflag=TRUE)}
summary.hmmmmod<-function(object,...){
hmmm.model.summary(object,printflag=TRUE)}
print.hmmmfit<-function(x,aname=" ",printflag=FALSE,...){
hmmm.model.summary(x$model,x,aname=aname,printflag=printflag,printhidden=0)}

hmmm.model.summary<-function(modelfull,fitmod=NULL,printflag=TRUE,aname="modfit",printhidden=0){
names<-modelfull$names
if (!is.null(fitmod)) {
a<-fitmod
a$df=a$df+dim(a$Zlist$DMAT)[1]-dim(a$Zlist$DMAT)[2]-a$model$modello$strata

if(printhidden==0){
cat("\n")
cat("SUMMARY of MODEL:", aname )


cat("\nOVERALL GOODNESS OF FIT:")
cat("\n")
cat("    Likelihood Ratio Stat (df=",a$df,"):  Gsq = ",
    round(a$Gsq,5))
if (a$df > 0) cat(" (p = ",signif(1-pchisq(a$Gsq,a$df),5),")")

cat("\n")
sm <- 100*length(a$m[a$m < 5])/length(a$m)  
if ((sm > 75)&(a$df > 0)) {
    cat("\n    WARNING:", paste(sm,"%",sep=""),
    "of expected counts are less than 5. \n")
    cat("             Chi-square approximation may be questionable.")




}
cat("\n")

}


if(printhidden==1){
cat("\n")
cat("SUMMARY of MODEL:", aname )



cat("   
 (df=",a$df,"):  Loglik = ",
    round(a$Gsq,5))





}
cat("\n")


if(printhidden==2){
cat("\n")
cat("SUMMARY of MODEL:", aname )



cat("
df=",a$df)
    






}
cat("\n")








}
printflagT<-TRUE
if(printflagT){
model<-modelfull$modello
nmarg<-length(model)-3
#strata<-c(model$strata)
levs<-c(model$livelli)
#print(nmarg)
#print(levs)
np<-prod(levs)-1	
C<-matrix("",np,1)
M<-matrix("",np,1)
T<-matrix("",np,1)

npar<-matrix(0,np,1)
iiii<-0
for (mi in 1:nmarg){

marginal<-model[[mi]]
margindex<-c(marginal$marg)
margset<-marginal$int
#print(margset)
types<-c(marginal$types)


nint<-length(margset)
#print(nint)
for (ii in 1:nint){
margint<-c(margset[[ii]])
npar[ii+iiii,1]<-prod(levs[margint]-1)
#print(c(marginal$marg),collapse="")
M[ii+iiii,1]<-paste(c(marginal$marg),collapse="") 

C[ii+iiii,1]<-paste(c(margset[[ii]]),collapse="")
T[ii+iiii,1]<-paste(c(types[margset[[ii]]]),collapse="")
}
iiii<-iiii+nint
}

npar2<-cumsum(npar)
npar1<-npar2-npar+1
if(is.null(fitmod$L)){
MCTL<-cbind(C,M,T,npar,npar1,npar2 )
colnames(MCTL)<-c("inter.","marg.","type","npar","start","end")
MCTL<-MCTL[1:iiii,]
MCTL<-MCTL

if(!is.null(names)){
C<-matrix(MCTL[,1])
M<-matrix(MCTL[,2])

C1<-matrix("",dim(MCTL)[1],1)
M1<-matrix("",dim(MCTL)[1],1)


      for(j in 1:dim(MCTL)[1]){


      NC<-as.numeric(unlist(strsplit(C[j],split="")))

      NM<-as.numeric(unlist(strsplit(M[j],split="")))

      C1[j]<-paste(names[NC],collapse=".")
      M1[j]<-paste(names[NM],collapse=",")}

MCTL<-cbind(C,C1,M,M1,MCTL[,3:6] )

colnames(MCTL)<-c("inter.","inter.names","marg.","marg.names","type","npar","start","end")


}

if (printflag) {print(MCTL,quote=FALSE)}
MCTL<-MCTL
}
else{
if (model$strata==1){

fitmod<-fitmod$L}
else{

a<-fitmod
fitmod<-matrix(fitmod$L,max(npar2),model$strata)}

if(!is.null(names)){


C1<-matrix("",dim(C)[1],1)
M1<-matrix("",dim(M)[1],1)


      for(j in 1:dim(C)[1]){


      NC<-as.numeric(unlist(strsplit(C[j],split="")))

      NM<-as.numeric(unlist(strsplit(M[j],split="")))

      C1[j]<-paste(names[NC],collapse=".")
      M1[j]<-paste(names[NM],collapse=",")}
C<-C1
M<-M1
}

C<-rep(C,as.numeric(npar))
M<-rep(M,as.numeric(npar))
T<-rep(T,as.numeric(npar))
MCTL<-cbind(C,M,T,round(fitmod,6))
colnames(MCTL)<-c("inter.","marg.","type",paste("STRATA",(1:model$strata),sep="_"))
if (printflag) {print(MCTL,quote=FALSE)}


}
MCTL<-MCTL
}

}


#  Author: Roberto Colombi
        #          Dept IGIF
        #          Univ of Bergamo ( last update: 12/04/2012)

#used to define recursive logits 
recursive<-function(...){
n<-nargs()
all.logit<-list(...)
cocacontr<-list()
 for(i in 1:n) {
if(!is.matrix(all.logit[[i]])){cocacontr[[i]]=0}
else{
X<-all.logit[[i]]
XX<-X
for(j in 1:nrow(X)){
X[j,]<-as.numeric(all.logit[[i]][j,]==1)
XX[j,]<-as.numeric(all.logit[[i]][j,]==-1)

}

cocacontr[[i]]<-rbind(X,XX)
}

}
cocacontr
}


hmmm.model.X<-function(marg,lev,names,Formula=NULL,strata=1,fnames=NULL,cocacontr=NULL,ncocacontr=NULL,replace=TRUE){
str<-prod(strata)
model<-hmmm.model(marg=marg,lev=lev,names=names,strata=str,cocacontr=cocacontr)
model<-create.XMAT(model,Formula=Formula,strata=strata,fnames=fnames,
cocacontr=cocacontr,ncocacontr=ncocacontr,replace=replace)
}
#  Author: Roberto Colombi
        #          Dept IGIF
        #          Univ of Bergamo ( last update: 25/07/06)


create.XMAT<-function(model,Formula=NULL,strata=1,fnames=NULL,cocacontr=NULL,ncocacontr=NULL,replace=TRUE){
if(is.null(cocacontr)){
cocacontr<-as.list(rep("contr.treatment",length(strata) ))
ncocacontr<-strata-1
}
if(is.null(ncocacontr)){ncocacontr<-strata-1}
descr<-hmmm.model.summary(model,printflag=FALSE)
if(is.null(model$names)){
intnames<-paste("int",descr[,1],sep="_")
descr[,1]<-intnames}
else{intnames<-descr[,2]
descr[,1]<-intnames}

if(!is.null(names(Formula))){
reo<-match(descr[,1],names(Formula))
Formula<-Formula[reo]}
npar<-as.numeric(descr[,"npar"])
if (is.null(fnames)){
fnames=paste("f",1:length(strata),sep="_")}
#else {fnames<-c("f_0",fnames)}
factlist<-list()
for (i in 1:dim(descr)[1]){
    np<-npar[i]*prod(strata)
    #fact<-matrix(1,np)
     #rep1<-1
    intfact<-data.frame(gl(npar[i],1,np))
    rep2<-npar[i]
       for(ii in 1:length(strata)){
           factii<-gl(strata[ii],rep2,np)
          contrasts(factii)<-eval(cocacontr[[ii]])
           contrasts(factii,ncocacontr[ii])<-contrasts(factii)[,1:ncocacontr[ii]]
           rep2<-rep2*strata[ii]
          intfact<-cbind(intfact,factii)
       }
      names(intfact)<-c(intnames[i],fnames)
     factlist[[i]]<-intfact
}
XL<-as.list(rep("zero",dim(descr)[1]))
px<-0
Xnames<-"zero"
 for (i in 1:dim(descr)[1]){
     
 if (Formula[[i]]=="zero"){
      px<-c(px,0)}

      else{
      XL[[i]]<-model.matrix(as.formula(Formula[[i]]),data=factlist[[i]])
      px<-c(px,dim(XL[[i]])[2])
Xnames<-c(Xnames,colnames(XL[[i]]))
}
}
px<-px[-1]
Xnames<-Xnames[-1]
XX<-matrix(0,sum(npar)*prod(strata),sum(px))
or<-rep(rep((1:dim(descr)[1]),times=npar),prod(strata))
pstart<-0
  for(i in 1:dim(descr)[1]){
  if (Formula[[i]]!="zero"){
  XX[or==i,(pstart+1):(pstart+px[i])]<-XL[[i]]
  pstart<-pstart+px[i]}
  }
if (replace) {
colnames(XX)<-Xnames
model$matrici$X<-XX
model$matrici$E<-t(create.U(model$matrici$X))
if ( sum(abs(model$matrici$E)) == 0){
model$functions$h.fct<-0
model$functions$derht.fct<-0}

else{
model$functions$h.fct<-make.h.fct(model,E=FALSE)
 model$functions$derht.fct<-make.derht.fct(model,E=FALSE)
 }
model$functions$L.fct<-make.L.fct(model,TRUE)
model$functions$derLt.fct<-make.derLt.fct(model,TRUE)
class(model)<-"hmmmmod"

model}
else{XL}
}
     
            

  


anova.hmmmfit<-function(object,objectlarge,...){t<-
hmmm.hmmm.anova(object,objectlarge)
t}

hmmm.hmmm.anova<-function(modelA,modelB){
a<-modelA
if(class(a)=="hmmmfit"){
modelA$df=a$df+dim(a$Zlist$DMAT)[1]-dim(a$Zlist$DMAT)[2]-a$model$modello$strata}
pA<-signif(1-pchisq(a$Gsq,modelA$df),5)
a<-modelB
if(class(a)=="hmmmfit"){
modelB$df=a$df+dim(a$Zlist$DMAT)[1]-dim(a$Zlist$DMAT)[2]-a$model$modello$strata}
pB<-signif(1-pchisq(a$Gsq,modelB$df),5)
Gsq<-abs(modelA$Gsq-modelB$Gsq)
dof<-abs(modelA$df-modelB$df)
pvalue<-(1-pchisq(abs(Gsq),dof))
anova.table<-matrix(c(modelA$Gsq,modelB$Gsq,Gsq,modelA$df,modelB$df,dof,pA,pB,pvalue),3,3,
dimnames = list(c("model A", "model B","LR test"), c("statistics value",
"df", "pvalue")))
#class(anova.table)<-"hmmmanova"
#anova.table<-anova.table

}

