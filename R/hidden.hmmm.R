

#Filtering of an arbitrary Hidden Markov Model

f.filter <- function(b,pr,p0)
{
  ## Purpose: Computing the filter probabilities fx and the predictions fy
  ##          of the data in a HMM
  ## ----------------------------------------------------------------------
  ## Arguments: b=conditional observation densities (b(k,t)=p(y_t|x_t=k))
  ##            pr=transition probabilites, p0=initial probability
  ## ----------------------------------------------------------------------
  ## Author: Hans-Ruedi Kuensch, Date:  2.4.03/12.9.03
  Tt <- ncol(b)
  m <- nrow(b)
  fx <- matrix(p0,nrow=m,ncol=Tt)
  fy <- rep(1,Tt)
  fx[,1] <- b[,1]*fx[,1]
  fy[1] <- sum(fx[,1])
  fx[,1] <- fx[,1]/fy[1]
  for (ti in (2:Tt)) {
    fx[,ti] <- b[,ti]*(fx[,ti-1]%*%pr)
    fy[ti] <- sum(fx[,ti])
    fx[,ti] <- fx[,ti]/fy[ti]
  }
  list(fx=fx,fy=fy)
}



#Smoothing

f.smooth <- function(b,fx,fy,pr)
{
  ## Purpose: Smoothing in a HMM
  ## ----------------------------------------------------------------------
  ## Arguments: b=observation densities, fx=filter distributions,
  ##            fy=predictions of observations,
  ##            pr=transition probabilities
  ## Output: margin=marginal smoother probabilities
  ##         pairs=smoother probabilities for successive pairs
  ## ----------------------------------------------------------------------
  ## Author: Hans-Ruedi Kuensch, Date:  3 Apr 2003, 14:30
  Tt <- ncol(b)
  m <- nrow(b)
  r <- matrix(1,nrow=m,ncol=Tt)
  sm2 <- array(1,dim=c(m,m,Tt))
  for (ti in (Tt:2)) {
    mat <- t(matrix(r[,ti]*b[,ti],m,m))
    mat <- pr*mat/fy[ti]
    r[,ti-1] <- apply(mat,1,sum)
    sm2[,,ti-1] <- t( matrix(fx[,ti-1],m,m)*mat)
  }
  sm2[,,1]<-NA
  list(margin=r*fx,pairs=sm2)
}





hidden.emfit<-function(y,model.obs,model.lat=NULL,
      noineq = TRUE, maxit = 10, maxiter=100, norm.diff.conv = 1e-05, 
    norm.score.conv = 1e-05, y.eps = 0,  mup = 1, step = 1,printflag=0,old.tran.p=NULL,
bb=NULL)
{
#-------------------------------------------------------------------
#CORPO FUNZIONE
#INIZIALIZZAZIONE
#modaltà osservate più latenti
if(!class(model.lat)=="hmmmmod"){
nstrata<-model.obs$modello$livelli[1]
marglat<-c("l-l")
marglat<-marg.list(marglat,mflag="m")
model.lat<-hmmm.model(marg=marglat,lev=c(nstrata,nstrata),names=c("lat","laggedlat"),
X=diag(1,nstrata*nstrata-1)
)

}

lev<-model.obs$modello$livelli
#p-1 sono osservate attuali
p<-length(lev)

#modalità osservate
levobs<-lev[-1]
#numero stati latenti
strata<-model.obs$modello$livelli[1]

#inizializzo prob transizione con distrib equiprobabilità
if(is.null(old.tran.p)){ 
#old.tran.p<-matrix(runif(strata^2,0,1),strata,strata)
old.tran.p<-matrix(1,strata,strata)
old.tran.p<-prop.table(old.tran.p,2)
}
A<-diag(1,strata)-old.tran.p
A<-rbind(A,matrix(1,1,strata))
old.iniz.p<-solve(t(A)%*%A)%*%t(A)
old.iniz.p<-old.iniz.p[,strata+1]


#inizializzo working data set  osservazioni complete: 
#latenti attuale e osservate
r<-apply(y,1,function(x) matrix(rep(x,strata),ncol=length(x)))
r<-matrix(r,ncol=dim(y)[2],byrow=TRUE)
oggi<-gl(strata,1,dim(r)[1])
ieri<-gl(strata,strata,dim(r)[1])
r<-cbind(oggi,r)
r<-as.data.frame(r)


#inizializzo working data per letente attuale e ritardata
#la attuale varia più velocemente
oggi<-gl(strata,1,strata^2*dim(y)[1])
ieri<-gl(strata,strata,strata^2*dim(y)[1])
rrr<-cbind(oggi,ieri)
rrr<-as.data.frame(rrr)
#inizializzo prob a posteriori degli stati
#ESTEP
#bb<-c(table(as.data.frame(y)))
if(is.null(bb)){
bb<-matrix(runif(prod(lev),0,1),strata,prod(levobs))
bb<-prop.table(bb,1)
}
#bb<-c(table(as.data.frame(y)))
#bb<-matrix(bb,strata,length(bb),byrow=TRUE)
#bb<- prop.table(bb,1)

ii<-cumprod(c(1,levobs[-1]))
#ii<-t(apply(iiy,1,cumprod))
#/(max(y[,1]-1,1))
#ii<-cbind(ii,rep(1,dim(ii)[1])
#ii<-(cumprod(rev(ii))


yii<-(y-1)%*%ii+1
bb<-bb[,yii]


#bb<-matrix(1/dim(y)[1],nrow=2,ncol=dim(y)[1],byrow=TRUE)

ff<-f.filter(bb,t(old.tran.p),old.iniz.p)

ffs<-f.smooth(bb,ff$fx,ff$fy,t(old.tran.p))
#probabilità congiunte stati tempo t e tempo t-1
#smoothed
smooth.p<-array(
apply(X=ffs$pairs,FUN=function(x) t(x),
MARGIN=c(3))
,dim(ffs$pairs))

#inizilizzazioni per  condizioni convergenza

oldpar<-rep(999,prod(lev)-1+strata^2-1-2*(strata-1))
Lold<-999999

par.conv<-tran.conv<-L.conv<-99999999
#EM iterations----------------------------------------------
iter<-1
while ( ((par.conv> norm.diff.conv)||(L.conv > norm.score.conv))
     &(iter< maxiter))
   {


#MSTEP

#tabella frequenze   congiunte delle osservate e latente attuale  per M step
r$count<-c(ffs$margin)
oldm<-xtabs(count~. ,data =r)
#edit(as.data.frame(oldm))

#tabella frequenze congiunte latente attuale e ritardata per M step
rrr$count<-c(smooth.p)
oldmmm<-xtabs(count~. ,data =rrr)





#stima modello hmmm definito da hmm.hmm.model usando come vettore frequenze congiunte oldm 
#sink("o")
fit.obs<-hmmm.mlfit(c(oldm), model.obs, maxit =maxit, norm.diff.conv =  norm.diff.conv, 
    norm.score.conv = norm.score.conv, y.eps =y.eps, 
    mup = mup, step = step) 
fit.lat<-hmmm.mlfit(c(oldmmm), model.lat, maxit =maxit, norm.diff.conv =  norm.diff.conv, 
    norm.score.conv = norm.score.conv, y.eps =y.eps, 
    mup = mup, step = step) 
#sink()

p1<-fit.lat$L[-c(strata:(2*strata-2))]
p2<-fit.obs$L[-c(1:strata-1)]
newpar<-c(p1,p2)
#forma tabellare frequenze teoriche congiunte obs + latente attuale
newm<-array(fit.obs$m,lev)
newmmm<-array(fit.lat$m,c(strata,strata))
#da forma tabellare frequenze teoriche congiunte
#ricavo probabilità di transizione
new.tran.p<-prop.table(newmmm,2)
#calcolo distr.  invariante della latente
A<-diag(1,strata)-new.tran.p
A<-rbind(A,matrix(1,1,strata))
new.iniz.p<-solve(t(A)%*%A)%*%t(A)
new.iniz.p<-new.iniz.p[,strata+1]


#ESTEP
#in forma tabellare frequenze teoriche congiunte
#marginalizzo rispetto latente ritardata
#newm<-apply(newm,c(1:p,(p+2):length(lev)),sum)
#calcolo probabilità osservate dato stato latente attuale
b<-prop.table(newm,1)

#per ogni istante temporale calcolo:
#probabilità osservazioni tempo t 
#dato ogni stato latente
b<-matrix(b,nrow=strata)
Ptobs<-b
ii<-cumprod(c(1,levobs[-1]))
#ii<-t(apply(iiy,1,cumprod))
#/(max(y[,1]-1,1))
#ii<-cbind(ii,rep(1,dim(ii)[1])
#ii<-(cumprod(rev(ii))


yii<-(y-1)%*%ii+1
bb<-b[,yii]
#FILTERING + SMOOTHING
# bb[,1]<-1
ff<-f.filter(bb,t(new.tran.p),new.iniz.p)

ffs<-f.smooth(bb,ff$fx,ff$fy,t(new.tran.p))
#probabilità congiunte stati tempo t e tempo t-1
#smoothed
smooth.p<-array(
apply(X=ffs$pairs,FUN=function(x) t(x),
MARGIN=c(3))
,dim(ffs$pairs))
#verosimiglianza osservazioni
Lnew<-sum(log(ff$fy))
#CRITERI CONVERGENZA
par.conv<-max(abs(newpar-oldpar))
L.conv<-quantile(abs(c(Lnew-Lold)),0.75)
tran.conv<-max(abs(c(new.tran.p-old.tran.p)))
oldpar<-newpar
old.tran.p<-new.tran.p
Lold<-Lnew
if(printflag > 0){
if(iter%%printflag==0){ cat("loglik= ",Lnew," Lconv= ",L.conv, " parconv= ",par.conv, " niter= ",iter)
cat("\n")
}
}

iter<-iter+1
}  
#END EM iterations----------------------------------------------------------------------------------------------------
#restituisco working models, prob transizione, distr iniziale, 
#prob stati smmothed e filtrate, probabilità osservazioni tempo t dato passato , max verosim marginale, condizioni convergenza
#hmlist<-list(model=list(lat=fit.lat,obs=fit.obs),initial=new.iniz.p,
#Ptr<-new.tran.p,filter=ff,smooth=ffs, 
#Gsq=Lnew,conv=c(iter,par.conv,tran.conv,L.conv))
#dl<-hmmm.model.summary(model.lat,fit.lat,printflag=FALSE)
#do<-hmmm.model.summary(model.obs,fit.obs,printflag=FALSE)
#do<-do[-c(1:strata-1),]
#dl<-dl[-c(strata:(2*strata-2)),]
#dof<-(dim(model.obs$matrici$X)[1]+dim(model.lat$matrici$X)[1]-
#dim(model.obs$matrici$X)[2]-dim(model.lat$matrici$X)[2])
dof<-fit.obs$df+fit.lat$df
fit.obs$Gsq<-Lnew
fit.obs$Xsq<-NULL
fit.obs$Wsq<-NULL
fit.obs$df<-dof
fit.lat$Gsq<-NaN
fit.lat$Xsq<-NaN
fit.lat$Wsq<-NaN
#fit.lat$df<-


hmlist<-list(model.lat=fit.lat,model.obs=fit.obs,vecpar=list(lat=p1,obs=p2),initial=new.iniz.p,
Ptr=new.tran.p,Ptobs=Ptobs  ,filter=ff,smooth=ffs,conv=list(niter=iter,par.conv=par.conv,transprob.conv=tran.conv,Lconv=L.conv))
class(hmlist)<-"hidden"
hmlist
}
print.hidden<-function(x,printflag=FALSE,...){fitted<-x
#print(fitted$model.lat,printflag=printflag,aname="transition model",printhidden=TRUE)
hmmm.model.summary(fitted$model.lat$model,fitted$model.lat,aname="transition model",printflag=printflag,printhidden=2)
if(printflag==TRUE){
cat(" Principal effects relating  to the lagged latent variable, (#2)","\n", "must not be considered" , "\n")}
#print(fitted$model.obs,printflag=printflag,aname="observation model",printhidden=1)
hmmm.model.summary(fitted$model.obs$model,fitted$model.obs,aname="observation model",printflag=printflag,printhidden=TRUE)
if(printflag==TRUE){
cat(" Principal effects relating  to the latent variable, (#1)","\n", "must not be considered" , "\n")}
}
summary.hidden<-function(object,...){
fitted<-object
#print(fitted$model.lat,printflag=TRUE,aname="transition model")
cat("\n","Transition probabilities","\n")
print(fitted$Ptr)
#print(fitted$model.obs,printflag=TRUE,aname="observation model")
 cat("\n","Probabilities of observations (by rows) given the latent states (by columns)","\n")
print(t(fitted$Ptobs))}

hmm.hmm.anova<-function(modelloA,modelloB){
Gsq<-2*abs(modelloA$model.obs$Gsq-modelloB$model.obs$Gsq)
a<-modelloA$model.obs
if(class(a)=="hmmmfit"){
dfA=a$df+dim(a$Zlist$DMAT)[1]-dim(a$Zlist$DMAT)[2]-a$model$modello$strata}
pA<-signif(1-pchisq(a$Gsq,dfA),5)

b<-modelloB$model.obs
if(class(b)=="hmmmfit"){
dfB=b$df+dim(b$Zlist$DMAT)[1]-dim(b$Zlist$DMAT)[2]-b$model$modello$strata}
pB<-signif(1-pchisq(b$Gsq,dfB),5)


dof<-abs(dfA-dfB)
P<-signif(1-pchisq(Gsq,dof),5)
print(matrix(c(a$Gsq,b$Gsq,Gsq,dfA,dfB,dof,"","",P),3,3,
dimnames = list( c("model A", "model B","LR test"), c("statistics value", "dof","pvalue")  )    ),quote=FALSE)

}
anova.hidden<-function(object,objectlarge,...){t<-
hmm.hmm.anova(object,objectlarge)
t}
