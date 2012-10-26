#
# 
#
######################   Begin mphineq.fit ########################################## 
mphineq.fit  <- function(y, Z, ZF=Z, h.fct=0,derht.fct=0,d.fct=0,derdt.fct=0,
                   L.fct=0,derLt.fct=0,X=NULL,formula=NULL,names=NULL,lev=NULL,E=NULL,
maxiter=100,step=1,
                   norm.diff.conv=1e-5,norm.score.conv=1e-5,
                   y.eps=0,chscore.criterion=2,
                   m.initial=y,mup=1)
{

start.time <- proc.time()[3]  
version <- "mphineq.fit, version 1.0.1, 5/10/06"
Zlist<-cocadise(Z,formula=formula,lev=lev,names=names)
if(is.null(X)){X<-0}
inv <- solve        
y<-as.matrix(y) 
 
lenh<-0
if ((missing(h.fct))&(sum(abs(X)) != 0))
# 
# if the user inputs the linear predictor model design matrix,
# but does not input the constraint function h.fct then do... 

{
    if(is.null(E)){U <- create.U(X)}
else{U<-t(E)}
	
    if (sum(abs(U)) == 0) {h.fct <- 0}
    else {
       h.fct <- function(m) {
          t(U)%*%L.fct(m)
           }
   }

} 
# end "if missing h.fct and X != 0"
else 

{
   U <- "Not created within the program."
}
#
	
if ((is.function(derht.fct)==FALSE)&(sum(abs(X)) != 0)&(is.function(derLt.fct)==TRUE))
# 
# if the user inputs the linear predictor model design matrix,
# but does not input the constraint function h.fct then do... 
{
    U <- create.U(X)
    if (sum(abs(U)) == 0) {derht.fct <- 0}
    else {
       derht.fct <- function(m) {
          derLt.fct(m)%*%U
          

       }
   }
} 




 
	 
#########################################################################
 
if ((is.function(h.fct)==TRUE)||(is.function(d.fct)==TRUE)||(class(formula)
=="formula"))  {

#
   lenm <- length(y);
   m <- as.matrix(c(m.initial)) + y.eps 
   m[m==0] <- 0.01
  
   p<-m*c(1/Z%*%t(Z)%*%m)

		  

    
  # Zlist<-cocadise(Z,formula=formula,lev=lev)
   xi <- Zlist$IMAT%*%log(m)
if ((is.function(d.fct)==FALSE)&(is.function(h.fct)==FALSE)) {
p  <- as.matrix(exp(Zlist$DMAT%*%xi))
         p<-(p/sum(p))
        m<-p*c(Z%*%t(Z)%*%y)
}

   if (is.function(h.fct)==TRUE){
   	h <- hobs <- h.fct(m)
   lenh <- length(h) 

if (is.function(derht.fct)==FALSE) {
     H <- num.deriv.fct(h.fct,m) 
   }
   else {
     H <- derht.fct(m)
   }  
   #HtDHinv <-
 HtDHinvobs <- inv(t(H)%*%(H*c(m))) #HtDHinv <- inv(t(H)%*%Dm%*%H)
  # HHtDHinv <- H%*%HtDHinv
    }
 	if (is.function(d.fct)==TRUE)
	{
#
# if d(.) is not the zero function, i.e. there is at least one constraint, then do...
# 
   
   d <- dhobs <- d.fct(m)
 lend <- length(d)

 
 ##}
  
  if (is.function(derdt.fct)==FALSE)
   {
     DH <- num.deriv.fct(d.fct,m) 
   }
   else {
     DH <- derdt.fct(m)
 

   }  
}	
#p  <- as.matrix(exp(Zlist$DMAT%*%xi))
        # p<-(p/sum(p))
        # m<-p*c(Z%*%t(Z)%*%y)
		  
  
		
  	
	
  	 
   #lagrange multiplers not needed and computed by this program
   lam <- matrix("NA",lenh,1)
   Dm <- diag(c(m+1e-08))-((ZF*c(m))%*%t(ZF*c(p)))
   
   if (!is.matrix(Dm)) {return ("unable to reach convergence")}






   
   norm.score <- 999999   
   
   theta<-xi
  
   iter <- 0
   step.iter <- 0
   norm.diff <-  10
  
  # cat("\n",version,", running...\n")
   while ( ((norm.diff > norm.diff.conv)||(norm.score > norm.score.conv))
     &(iter< maxiter))
   {
     
       
       

	 

	   			   #programazione quadratica programazione quadratica  programazione quadratica programazione quadratica
	   
	qpmatr<-t(Zlist$DMAT)%*%Dm%*%Zlist$DMAT
  
	   if ((is.function(d.fct)==TRUE)&(is.function(h.fct)==TRUE)) {
          Amat<-cbind(t(Zlist$DMAT)%*%Dm%*%H,t(Zlist$DMAT)%*%Dm%*%DH)
          bvec<-rbind(-h,-d ) 
          }
########################################################################
		  else {
               if (is.function(h.fct)==TRUE){ 			 
              Amat<- t(Zlist$DMAT)%*%Dm%*%H
               bvec<- -h         }
               if (is.function(d.fct)==TRUE) {
             Amat<-t(Zlist$DMAT)%*%Dm%*%DH
              bvec<- -d }
if ((is.function(d.fct)==FALSE)&(is.function(h.fct)==FALSE)) {
           Amat<-matrix(0,nrow(qpmatr),1)
           bvec<-0
}
              }
     #           qpmatr<-t(Zlist$DMAT)%*%Dm%*%Zlist$DMAT
		
if (any(is.null(qpmatr))||any(is.na(qpmatr))) {
print("matrix in quadratic programming not positive def.")
return("matrix in quadratic programming not positive def.")}	  
	   	
As<- solve.QP(qpmatr,t(Zlist$DMAT)%*%(y-m), Amat, bvec, meq=lenh, factorized=FALSE)
	
#function used by optimize

 ff.fct<-function(steptemp){

 theta.temp <- theta + steptemp*matrix(As$solution)
        
         p  <- as.matrix(exp(Zlist$DMAT%*%theta.temp))
         #p<-p/sum(p) MODIFICA IN PROVA PER STRATI
         p<-p*c(1/Z%*%t(Z)%*%p)
         m<-p*c(Z%*%t(Z)%*%y) 
        
      if (is.function(h.fct)==FALSE) {
          h<-0}
      else{
         h  <- h.fct(m)
           }
        if (is.function(d.fct)==FALSE) {
         dd<-100 
          }
         else {
        dd <- d.fct(m)
         }  




        
        

        
         norm.score.temp <-
		 as.matrix(2/sum(y)*sum(y[y>0]*(log(y[y>0])-log(m[y>0]))))+
		 mup*sum(abs(h))
            -mup*sum(pmin(dd,dd*0))
         norm.score.temp
		 
}




stepco<-optimize(ff.fct, c(step*0.5^5, step), tol = 0.0001)
step.temp<-stepco$minimum
step.iter<-step.temp



     
         theta.temp <- theta + step.temp*matrix(As$solution)
		 
        

         norm.diff  <- sqrt(sum((theta-theta.temp)*(theta-theta.temp)))
        p  <- as.matrix(exp(Zlist$DMAT%*%theta.temp))
         #p<-(p/sum(p))  MODINPROVAPERSTRATI
          p<-p*c(1/Z%*%t(Z)%*%p)
         m<-p*c(Z%*%t(Z)%*%y)
		  
        if (is.function(h.fct)==FALSE) {
          h<-0}
        else{
         h  <- h.fct(m)
if (is.function(derht.fct)==FALSE) {
            H <- num.deriv.fct(h.fct,m)
         }
         else {
           H <- derht.fct(m)
         }

}

  if (is.function(d.fct)==FALSE) {
         d<-100 
          }
         else {
        d <- d.fct(m)
         }  


         Dm <- diag(c(m+1e-08))-((ZF*c(m))%*%t(ZF*c(p)))
       
         if (!is.matrix(Dm)) {return ("unable to reach convergence")}

        
		
        norm.score <- sum(abs(h))-sum(pmin(d,d*0))
         	 

       

       theta <- theta.temp 


       iter <- iter + 1

       if(chscore.criterion==0){
      # cat("  iter=",iter, "[",step.iter,"]", 
       #    " norm.diff=",norm.diff," norm.score=", norm.score,"\n")
                                 }
   }
}

satflag<-dim(Zlist$DMAT)[1]-dim(Zlist$DMAT)[2]

if ((is.function(h.fct)==TRUE)||(  (class(formula)=="formula")&(satflag >1  ) )){
##################

if ((is.function(h.fct)==FALSE)&(class(formula)=="formula")){

M<-cbind(Zlist$DMAT,matrix(1,nrow(Zlist$DMAT)  ))
H<-create.U(M)
H<-diag(1/c(m))%*%H
hobs<-t(H)%*%log(m)
lenh <- length(hobs)
#print(lenh)
 lam <- matrix("NA",lenh,1)

HtDHinvobs <- inv(t(H)%*%(H*c(m))) 
#H<-t(H)
}


if ((is.function(h.fct)==TRUE)&(class(formula)=="formula")){
#lenh <- length(t(H)%*%log(m))
M<-cbind(Zlist$DMAT,matrix(1,nrow(Zlist$DMAT)  ))
H2<-create.U(M)
H2<-diag(1/c(m))%*%H2
H<-cbind(H,H2)
hobs<-t(H)%*%log(m)
lenh <- length(hobs)

 lam <- matrix("NA",lenh,1)
lenh<-lenh-dim(H2)[2]
#print(lenh)
HtDHinvobs <- inv(t(H)%*%(H*c(m))) 
#H<-t(H)
}


 HtDHinv <- inv(t(H)%*%(H*c(m)))     #HtDHinv <- inv(t(H)%*%Dm%*%H)
         HHtDHinv <- H%*%HtDHinv


                                          #Ninv <- diag(c(1/Z%*%t(Z)%*%y))
   p <- m*c(1/Z%*%t(Z)%*%y)                #p <- Ninv%*%m
   resid <- y-m
   covresid <- (H*c(m))%*%HtDHinv%*%t(H*c(m)) 
                                    #covresid <- Dm%*%H%*%HtDHinv%*%t(H)%*%Dm
   covm.unadj <- covm <- Dm -  covresid
   if (sum(ZF) != 0) {
       covm <- covm.unadj -  ((ZF*c(m))%*%t(ZF*c(m)))*c(1/Z%*%t(Z)%*%y)
      #covm <- covm.unadj - Ninv%*%Dm%*%ZF%*%t(ZF)%*%Dm
   }
    covp <- t(t((covm.unadj-((Z*c(m))%*%t(Z*c(m)))*c(1/Z%*%t(Z)%*%y))*
               c(1/Z%*%t(Z)%*%y))* c(1/Z%*%t(Z)%*%y)) 
   #covp <- Ninv%*%(covm.unadj-((Z*c(m))%*%t(Z*c(m)))*c(1/Z%*%t(Z)%*%y))%*%Ninv
   #covp <- Ninv%*%(covm.unadj-Ninv%*%Dm%*%Z%*%t(Z)%*%Dm)%*%Ninv 
   
 # Compute adjusted residuals...
     dcovresid <- diag(covresid)
     dcovresid[abs(dcovresid)<1e-8] <- 0
     adjresid <- resid
     adjresid[dcovresid > 0] <- resid[dcovresid>0]/sqrt(dcovresid[dcovresid>0])
 # end compute adjusted residuals.

   presid <- resid/sqrt(m) 
#----------------------------
   covlam <- HtDHinv
   Gsq <- as.matrix(2*sum(y[y>0]*(log(y[y>0])-log(m[y>0]))))
   Xsq <- as.matrix(t(y-m)%*%((y-m)*c(1/m)))
   #Xsq <- as.matrix(t(y-m)%*%Dminv%*%(y-m))
###do not compute wsd if inequalities are present
if(is.function(d.fct)==FALSE){
   Wsq <- as.matrix(t(hobs)%*%HtDHinvobs%*%hobs)}
else  {Wsq<-as.matrix("NA")}
   beta <- "NA"
   covbeta <-  "NA"
   covL <-  "NA"
   L <- "NA"
   Lobs <- "NA"
   Lresid <- "NA"
   
   if (sum(abs(X)) != 0) {
       L <- L.fct(m)
       Lobs <- L.fct(y+y.eps)
       if (is.function(derLt.fct)==FALSE) {
          derLt <- num.deriv.fct(L.fct,m)
       }
       else {
          derLt <- derLt.fct(m)
       }
       PX <- inv(t(X)%*%X)%*%t(X) 
       beta <- PX%*%L 
       covL <- t(derLt)%*%covm%*%derLt
       covbeta <- PX%*%covL%*%t(PX)
       Lres <- Lobs - L
       covLres <- t(derLt)%*%covresid%*%derLt
       dcovLres <-   diag(covLres)
       dcovLres[abs(dcovLres)<1e-8] <- 0
       Lresid <- Lres
       Lresid[dcovLres > 0] <- Lres[dcovLres>0]/sqrt(dcovLres[dcovLres>0])  
       # end compute adjusted Link residuals.
       

       lbeta <- ll <- c()
       for (i in 1:length(beta)) { 
          lbeta <- c(lbeta,paste("beta",i,sep=""))
       }  
       for (i in 1:length(L)) {
          ll <- c(ll,paste("link",i,sep=""))
       } 
       dimnames(beta) <- list(lbeta,"BETA")  
     if(!is.null(colnames(X))){dimnames(beta) <- list(colnames(X),"BETA")}
       dimnames(covbeta) <- list(lbeta,lbeta)   
       dimnames(L) <- list(ll,"ML LINK") 
       dimnames(Lobs) <- list(ll,"OBS LINK") 
       dimnames(covL) <- list(ll,ll)
       dimnames(Lresid) <- list(ll,"LINK RESID")

   }
} # end "if h.fct is not the zero function"
else {
#


# else, if h.fct = 0, i.e. there are no constraints, then do...
# 
   lenh <- 0 
   lenm <- length(y)
   if(is.function( d.fct)==FALSE){
   m <- as.matrix(c(m.initial))+y.eps 
   m[m==0] <- 0.01 
   xi <- log(m) 
   Dm <- diag(c(m))
   Dminv <- diag(c(1/m))
   s <- y-m
   norm.score <- sqrt(sum(s*s))
   theta <- xi 
   lentheta <- length(theta)
   iter <- 0
   norm.diff <- 10
  # cat("\n",version,", running...\n") 
   while ( ((norm.diff > norm.diff.conv)||(norm.score > norm.score.conv))
     &(iter< maxiter))
   {     
       
       A <- Dminv      
       thetanew <- theta + step*(s*c(1/m))   #thetanew <- theta + step*A%*%s
       norm.diff <- sqrt(sum((theta-thetanew)*(theta-thetanew)))
       theta <- thetanew
       m <- exp(theta)          
       Dm <- diag(c(m))
       Dminv <- diag(c(1/m))
       s <- y-m 
       norm.score <- sqrt(sum(s*s))   
       iter <- iter + 1
      # cat("  iter=",iter, " norm.diff=",norm.diff," norm.score=", norm.score,"\n")
   } 
}
                                            #Ninv <- diag(c(1/Z%*%t(Z)%*%y)) 
   p <- m*c(1/Z%*%t(Z)%*%y)                  #p <- Ninv%*%m
   resid <- 0*y 
   covm.unadj <- covm <- Dm 
   covresid <- 0*covm 
   if (sum(ZF) != 0) {
       covm <- covm.unadj -  ((ZF*c(m))%*%t(ZF*c(m)))*c(1/Z%*%t(Z)%*%y)
      #covm <- covm.unadj - Ninv%*%Dm%*%ZF%*%t(ZF)%*%Dm
   }
   covp <- t(t((covm.unadj-((Z*c(m))%*%t(Z*c(m)))*c(1/Z%*%t(Z)%*%y))*
               c(1/Z%*%t(Z)%*%y))* c(1/Z%*%t(Z)%*%y)) 
   #covp <- Ninv%*%(covm.unadj-((Z*c(m))%*%t(Z*c(m)))*c(1/Z%*%t(Z)%*%y))%*%Ninv
   #covp <- Ninv%*%(covm.unadj-Ninv%*%Dm%*%Z%*%t(Z)%*%Dm)%*%Ninv 
                                            
   adjresid <-  0*y 
   presid <- 0*y
   covlam <- as.matrix(0);
#-----------------------------------
lam <- as.matrix(0)
   Gsq <- as.matrix(2*sum(y[y>0]*(log(y[y>0])-log(m[y>0]))))
   Xsq <- as.matrix(t(y-m)%*%((y-m)*c(1/m)))   
   #Xsq <- as.matrix(t(y-m)%*%Dminv%*%(y-m))
  #do not compute wsd if inequalities are present
if(is.function(d.fct)==FALSE){
 Wsq <- as.matrix(0) }
else  {Wsq<-as.matrix("NA")}
   beta <- "NA"
   covbeta <-  "NA"
   covL <-   "NA"
   L <-  "NA"
   Lresid <- "NA"
   Lobs <- "NA"
   if (sum(abs(X)) != 0) {
       L <- L.fct(m)
       Lobs <- L.fct(y)
       if (is.function(derLt.fct)==FALSE)  {
         derLt <- num.deriv.fct(L.fct,m)
       }
       else {
         derLt <- derLt.fct(m)
       }
       PX <- inv(t(X)%*%X)%*%t(X) 
       beta <- PX%*%L 
       covL <- t(derLt)%*%covm%*%derLt
       Lresid <- 0*L
       covbeta <- PX%*%covL%*%t(PX) 
       lbeta <- ll <- c() 
       for (i in 1:length(beta)) { 
          lbeta <- c(lbeta,paste("beta",i,sep=""))
       }  
       for (i in 1:length(L)) {
          ll <- c(ll,paste("link",i,sep=""))
       } 
       dimnames(beta) <- list(lbeta,"BETA")  
       dimnames(covbeta) <- list(lbeta,lbeta)
       dimnames(Lobs) <- list(ll,"OBS LINK")   
       dimnames(L) <- list(ll,"ML LINK")  
       dimnames(covL) <- list(ll,ll) 
       dimnames(Lresid) <- list(ll,"LINK RESID") 
   }
}# end "else if h.fct = 0"
#
# ASSIGN LABELS...
#
  lm <- ly <- lp <- lbeta <- lr <- lar <- lpr <- ll <- llam <- c()
 


if(is.null(rownames(y))){
  for (i in 1:lenm) {
     lm <- c(lm,paste("m",i,sep=""))
     ly <- c(ly,paste("y",i,sep=""))
     lp <- c(lp,paste("p",i,sep=""))
     lr <- c(lr,paste("r",i,sep=""))
     lar <- c(lar,paste("adj.r",i,sep=""))
     lpr <- c(lpr,paste("pearson.r",i,sep=""))
  }
 }

else{ly<-c( paste("y(",rownames(y),")"))
     lm<- c( paste("m(",rownames(y),")"))
     lp <-c( paste("p(",rownames(y),")"))
     lr<- c(paste("r(",rownames(y),")"))
     lar<-c(paste("a.r(",rownames(y),")"))
     lpr<-c(paste("p.r(",rownames(y),")"))}





  for (i in 1:length(lam)) {
     llam <- c(llam,paste("lambda",i,sep=""))
  } 
  dimnames(y) <- list(ly,"OBS")
  dimnames(m) <- list(lm,"FV")
  dimnames(p) <- list(lp,"PROB")
  dimnames(resid) <- list(lr,"RAW RESIDS")
  dimnames(presid) <- list(lpr,"PEARSON RESIDS")
  dimnames(adjresid) <- list(lar, "ADJUSTED RESIDS")  
  dimnames(lam) <- list(llam,"LAGRANGE MULT") 
  dimnames(covm) <- list(lm,lm)
  dimnames(covp) <- list(lp,lp)
  dimnames(covresid) <- list(lr,lr)
  dimnames(covlam) <- list(llam,llam)
  dimnames(Xsq) <- list("","PEARSON SCORE STATISTIC")
  dimnames(Gsq) <- list("","LIKELIHOOD RATIO STATISTIC")
  dimnames(Wsq) <- list("","GENERALIZED wALD STATISTIC")

if (is.function(derht.fct)==FALSE) {derht.fct <- "Numerical derivatives used."}
if (is.function(derLt.fct)==FALSE) {derLt.fct <- "Numerical derivatives used."}
#cat("\n")
#cat(" Time Elapsed:", proc.time()[3]-start.time,"seconds")
#cat("\n")
lenh<-lenh*(is.function(h.fct)==TRUE)
modlist<-list(y=y,m=m,covm=covm,p=p,covp=covp, 
lambda=lam,covlambda=covlam,
resid=resid,presid=presid,adjresid=adjresid,covresid=covresid,
Gsq=Gsq,Xsq=Xsq,Wsq=Wsq,df=lenh,
beta=beta,covbeta=covbeta, Lobs=Lobs, L=L,covL=covL,Lresid=Lresid,
iter=iter, norm.diff=norm.diff,norm.score=norm.score,
h.fct=h.fct,derht.fct=derht.fct,L.fct=L.fct,derLt.fct=derLt.fct,X=X,
U=U,Z=Z,ZF=ZF,Zlist=Zlist,version=version)
class(modlist)="mphfit"
modlist
}
######################   end mph.fit   ########################################## 
summary.hmmmfit<-function(object,cell.stats=TRUE,...){mph.summary(object,cell.stats=cell.stats,model.info=FALSE)}
summary.mphfit<-function(object,...){mph.summary(object,cell.stats=TRUE,model.info=TRUE)}
print.mphfit<-function(x,...){mph.summary(x,cell.stats=FALSE,model.info=FALSE)}

######################   Begin mph.summary ######################################  
mph.summary <- function(mph.out,cell.stats=FALSE,model.info=FALSE) {
#
#  This function is used in conjunction with the ML fitting 
#  function `mph.fit'.  It computes and prints a collection 
#  of summary statistics of the fitted MPH model given in `mph.out',
#  which is the result of `mph.fit'.  That is, 
#  mph.out <- mph.fit(y,Z,ZF,...)
#  
#  Author:  Joseph B. Lang, 
#           Dept of Statistics 
#              and Actuarial Science
#           Univ of Iowa, Iowa City, IA 52242
#           8/16/01
#  Last Updated:  3/30/04
#
#  INPUT 
#     Required:
#          mph.out = result of `mph.fit'
#
#     Optional:
#          cell.stats = logical variable indicating whether cell
#                       specific statistics are to be output
#                       (default: cell.stats=FALSE)?
#          model.info = logical variable indicating whether model
#                       information is to be output
#                       (default: model.info=FALSE)
#                    
# 
a <- mph.out 
if(class(a)=="hmmmfit"){
a$df=a$df+dim(a$Zlist$DMAT)[1]-dim(a$Zlist$DMAT)[2]-a$model$modello$strata}
cat("\nOVERALL GOODNESS OF FIT:")
cat("\n")
cat("    Likelihood Ratio Stat (df=",a$df,"):  Gsq = ",
    round(a$Gsq,5))
if (a$df > 0) cat(" (p = ",signif(1-pchisq(a$Gsq,a$df),5),")")
cat("\n")
cat("    Pearson's Score Stat  (df=",a$df,"):  Xsq = ",
    round(a$Xsq,5))
if (a$df > 0) cat(" (p = ",signif(1-pchisq(a$Xsq,a$df),5),")")
cat("\n")
#if(a$Wsq != "NA"){
#cat("    Generalized Wald Stat (df=",a$df,"):  Wsq = ",
#    round(a$Wsq,5))
#if (a$df > 0) cat(" (p = ",signif(1-pchisq(a$Wsq,a$df),5),")")}
cat("\n")
sm <- 100*length(a$m[a$m < 5])/length(a$m)  
if ((sm > 75)&(a$df > 0)) {
    cat("\n    WARNING:", paste(sm,"%",sep=""),
    "of expected counts are less than 5. \n")
    cat("             Chi-square approximation may be questionable.")
}
cat("\n") 
if (cell.stats==TRUE) {

if (a$L[1] != "NA") {
  cat("\nLINEAR PREDICTOR MODEL RESULTS...")
  cat("\n")
  sbeta <- as.matrix(sqrt(abs(diag(a$covbeta))))
  z <- a$beta/sbeta
  pval <- 2*(1-pnorm(abs(z)))
  dimnames(sbeta)[2] <- "StdErr(BETA)"
  dimnames(z)[2] <- "Z-ratio"
  dimnames(pval)[2] <- "p-value"
 # betanames<-colnames(a$matrici$X)

if(class(a)=="mphfit"||a$model$modello$strata >1){
  print(cbind(a$beta,sbeta,z,pval))

  cat("\n")}
  stdL <- as.matrix(sqrt(diag(a$covL)))
  dimnames(stdL)[2] <- "StdErr(L)"
  LLL<-round(cbind(a$Lobs,a$L,stdL,a$Lresid),4)
if(class(a)=="hmmmfit"){
if(is.null(a$model$names)){
descr<-hmmm.model.summary(a$model,printflag=FALSE)
descr<-rep(descr[,1],descr[,4])}
else{
descr<-hmmm.model.summary(a$model,printflag=FALSE)

descr<-hmmm.model.summary(a$model,printflag=FALSE)
descr<-rep(descr[,2],descr[,6])}

 # descr<-hmmm.model.summary(modello,print=FALSE)
  #intnames<-paste("int",descr[,1],sep="_")
  intnames<-descr
  rownames(LLL)<-rep(intnames,dim(a$model$matrici$Z)[2])}
  print(LLL)
  cat("\n")
}
############################if (cell.stats==TRUE) {
   stdm <- as.matrix(sqrt(diag(a$covm)))
   stdp <- as.matrix(sqrt(diag(a$covp)))
   dimnames(stdm)[2] <- "StdErr(FV)"
   dimnames(stdp)[2] <- "StdErr(PROB)"
   cat("\nCELL-SPECIFIC STATISTICS...")
   cat("\n")
   print(round(cbind(a$y,a$m,stdm,a$p,stdp,a$adjresid),5))
   cat("\n")
}
cat("\nCONVERGENCE STATISTICS...")
cat("\n")
cat("    iterations =",a$iter)
cat("\n")
cat("    norm.diff  =",signif(a$norm.diff,6))
cat("\n")
cat("    norm.score =",signif(a$norm.score,6))
cat("\n")


if (model.info==TRUE) {
cat("\nMODEL INFORMATION...")
cat("\n")
  if (a$L[1] != "NA") {
    cat("Linear Predictor Model Link Function  L.fct:\n")
    print(a$L.fct)
    cat("\n")
    cat("Derivative of Transpose Link Function derLt.fct: \n")
    print(a$derLt.fct)
    cat("\n")
    cat("Linear Predictor Model Design Matrix  X: \n")
    print(a$X)
    cat("\n")
    cat("U = Orthogonal Complement of X: \n")
    if (is.matrix(a$U)) print(a$U)
    else cat("\n  ",a$U,"\n")
  }
  cat("\n")  
  cat("Constraint Function  h.fct:\n")
  print(a$h.fct)
  cat("\n")
  cat("Derivative of Transpose Constraint Function derht.fct:\n")
  print(a$derht.fct)
  cat("\n")
  cat("Population Matrix Z: \n")
  print(a$Z)
  cat("\n")
  cat("Sampling Constraint Matrix ZF:\n")
  print(a$ZF)
  cat("\n")
}
cat("\nFITTING PROGRAM USED: ", a$version,"\n\n")
}
######################   end mph.summary     ######################################

######################   begin num.deriv.fct ######################################       
num.deriv.fct <- function(f.fct,m) {
# 
#   Author: Joseph B. Lang, Univ of Iowa
#   Created:  c. 8/25/00 (last update: 3/30/04)
#
#   The numerical derivative of the transpose of function f.fct is 
#   computed at the value m.  If f.fct is a mapping from 
#   Rp to Rq then the result is a pxq matrix.
#     I.e. Result is approximation to 
#         d f.fct(m)^T/d m.
# 
  eps <- (.Machine$double.eps)^(1/3)
  d <- eps * m + eps  
  lenm <- length(m)
  E <- diag(c(d)) 
  f1 <- f.fct(m+E[,1])
  lenf <- length(f1)
  Ft <- (f1-f.fct(m-E[,1]))/(2*d[1])
  for (j in 2:lenm) {
     Ft <- cbind(Ft,((f.fct(m+E[,j])-f.fct(m-E[,j]))/(2*d[j])))
  }
  dimnames(Ft) <- NULL
  t(Ft)
}
######################   end num.deriv.fct ######################################
   

######################   begin create.U    ###################################### 
create.U <- function(X) {
#
#   Author: Joseph B. Lang, Univ of Iowa
#   Created:  8/19/01  (last update: 3/30/04)
#
# This program creates a full-column rank matrix, U, with column space 
# equal to the orthogonal complement of the column space of X.  That is,
# U has column space equal to the null space of X^T.
#
#  Input:  X must be of full column rank
#  
  nrowX <- nrow(X)
  u <- nrowX - ncol(X)
  if (u == 0) {U <- 0}
  else {w.mat <- matrix(runif(nrowX*u,1,10),nrowX,u)
    U <- w.mat - X%*%solve(t(X)%*%X)%*%t(X)%*%w.mat
  }
  U
}
######################   end create.U     ######################################
 getnames<-function(dat,st=3,sep=" "){
#dat=aggregated data frame with frequencies in the last column
#st length of the string for every category name
#sep separetor of category names of a cell
fa<-dim(dat)[2]-1
a<-substr(dat[[1]],1,st)
for(i in 2:fa){
a<-paste(a,substr(dat[[i]],1,st),sep=sep)
}
y<-as.matrix(dat[[fa+1]])
rownames(y)<-a
y}

