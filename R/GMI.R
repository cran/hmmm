GMI<-function(freq,marg,lev,names,mflag="M"){
  marg<-marg.list(marg, mflag="m") 
interf<-hmmm.model(marg=marg,lev=lev,names=names)
interfun<-interf$functions$L.fct 
freq<-as.matrix(freq)
gmi<-apply(freq,2,interfun)
x<-hmmm.model.summary(interf,printflag=FALSE)
  #print(interf,quote=FALSE)
na<-rep(x[,2],times=x[,6])
for(i in 1:dim(x)[1]){
na[x[i,7]:x[i,8]]<-paste(x[i,2],1:x[i,6],sep="")}
rownames(gmi)<-na
colnames(gmi)<-paste("F",1:dim(gmi)[2],sep="")
gmi<-list(marginals=interf,gmi=gmi)
}

inv_GMI<-function(etpar,mod,start=rep(0,prod(mod$modello$livelli))){
  #given a vector of complete hyerarchical generalized marginal interactions  computes the joint probabilities
  #etpar=vettorog generalized marginal interactions
  #modello marginale oggetto classe hmmmmod
  #vector of starting values for the loglinear parameters zero vector is the default
  #if(is.null(start)) start<-rep(0,prod(mod$modello$livelli))
  interf<-mod$functions$L.fct 
  interfder<-mod$functions$derLt.fct 
  
  Zlist<-cocadise(mod$matrici$Z)  
  myfun<-function(x){x<-as.matrix(exp(Zlist$DMAT%*%x))
                     x<-x/sum(x)          
                     interf(x)-etpar}
  myder<-function(x){
    x<-as.matrix(exp(Zlist$DMAT%*%x))
    x<-x/sum(x)  
    (
      t( t(Zlist$DMAT)%*%
           diag(c(x+1e-08)) %*%interfder(x)
      )) }
  r<-nleqslv(rep(0,15),myfun,myder)$x
  r<-as.matrix(exp(Zlist$DMAT%*%r))
  r<-r/sum(r) 
}
