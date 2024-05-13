
akaike<-function(...,LRTEST=FALSE,ORDERED=FALSE,NAMES=NULL){
MOD<-list(...)
n<-length(MOD)
VAK<-matrix(0,1,8)
if(LRTEST){
#if(class(MOD[[1]])=="hidden"){

if(is(MOD[[1]],"hidden")){
	
a<-MOD[[1]]$model.obs
b<-MOD[[1]]$model.lat
adf1<-a$df+dim(a$Zlist$DMAT)[1]-dim(a$Zlist$DMAT)[2]-a$model$modello$strata
bdf1<-b$df+dim(b$Zlist$DMAT)[1]-dim(b$Zlist$DMAT)[2]-b$model$modello$strata
df1<-adf1+bdf1
Gsq1<-2*MOD[[1]]$model.obs$Gsq
}

else{
df1=MOD[[1]]$df+dim(MOD[[1]]$Zlist$DMAT)[1]-dim(MOD[[1]]$Zlist$DMAT)[2]-MOD[[1]]$model$modello$strata
y<-MOD[[1]]$y 
Gsq1 <--MOD[[1]]$Gsq+ 2*sum(y[y>0]*log(y[y>0]))
}
}
for(i in 1:n) {
GSQ<-df<-P<-0
#if(class(MOD[[i]])=="hidden"){
if(is(MOD[[i]],"hidden")){
	
if(is.null(MOD[[i]])){ VAK<-rbind(VAK,c(i,rep(NA,7))) }
if(!is.null(MOD[[i]]))
{
a<-MOD[[i]]$model.obs
b<-MOD[[i]]$model.lat
adf<-a$df+dim(a$Zlist$DMAT)[1]-dim(a$Zlist$DMAT)[2]-a$model$modello$strata
bdf<-b$df+dim(b$Zlist$DMAT)[1]-dim(b$Zlist$DMAT)[2]-b$model$modello$strata
npar<-length(MOD[[i]]$vecpar$obs)+length(MOD[[i]]$vecpar$lat)-adf-bdf
AK<-2*npar-2*a$Gsq

#################
if(LRTEST){
GSQ<- -2*a$Gsq+Gsq1
df<-adf+bdf-df1
P<-1-pchisq(GSQ,df)
}
#########################
#AK<-c(i,a$Gsq,npar,AK)
VAK<-rbind(VAK,c(i,a$Gsq,adf+bdf,npar,GSQ,df,P,AK))
}
}
else{
dfi=MOD[[i]]$df+dim(MOD[[i]]$Zlist$DMAT)[1]-dim(MOD[[i]]$Zlist$DMAT)[2]-MOD[[i]]$model$modello$strata
npar<-length(MOD[[i]]$y)-MOD[[i]]$model$modello$strata-dfi


y<-MOD[[i]]$y   
Gsqi <--MOD[[i]]$Gsq/2+ sum(y[y>0]*log(y[y>0]))
AK<--2*Gsqi+2*npar
###########################################################
if(LRTEST){
df<-dfi-df1
GSQ<--2*Gsqi+Gsq1
P<-1-pchisq(GSQ,df)
}
########################################################
VAK<-rbind(VAK,c(i,Gsqi,MOD[[i]]$df,npar,GSQ,df,P,AK))

}
}
VAK<-VAK[-1,]

Delta<-VAK[,8]-min(VAK[,8])
VAK<-cbind(VAK,matrix(Delta,length(Delta),1))
VAK[,-1]
if(is.null(NAMES)){
rownames(VAK)<-paste("model",1:n,sep="")}
else{rownames(VAK)<-NAMES}
if(LRTEST){
colnames(VAK)<-c("#model","loglik","dfmodel","npar","LRTEST","dftest","PVALUE","AIC","DeltaAIC")
}
else{

VAK<-VAK[,c(1,2,3,4,8,9)]
colnames(VAK)<-c("#model","loglik","dfmodel","npar","AIC","DeltaAIC")
}
if(ORDERED){r<-order(VAK[,dim(VAK)[2]])
VAK[r,]}
else{VAK}


}

