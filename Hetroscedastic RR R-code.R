#R-Code for Heteroscedastic Ridge estimators (RR)
#Proposed vs existing Hetroscedastic RR
#Authors: 1. Mashal, 2. Dr Syed Muhammad Asim and 3. Dr. Muhammad Suhail

rm(list=ls())
set.seed(1988)
library(MASS)

ssize=c(20,50,100)
pred=c(10)
corr=c(0.90,0.95,0.99,0.999)

MSE.r=matrix(0,15,length(corr))
MSE.n=array(0,dim = c(15,length(corr),length(ssize)))
MSE.p=array(0,dim = c(15,length(corr),length(ssize),length(pred)))


for(b in 1:length(pred)){
  p=pred[b]
  I=diag(1,p,p)
  N=5000                                   #Simulation runs
  
  for(a in 1:length(ssize)){
    n=ssize[a]
    
    z.sn=matrix(0,n,p)
    for (i in 1:p) {
      z.sn[,i]=rnorm(n,0,1)                 #Random no generation   
    }
    
    for(f in 1:length(corr)){
      rho=corr[f]
      x=matrix(0,n,p)
      
      #incorporating varying levels of correlation
      #Note, this code is for only P=10, researchers can choose the following terms accordinlgy
      
      x[,1]=sqrt(1-corr[1]^2)*z.sn[,1]+corr[1]*z.sn[,p]
      x[,2]=sqrt(1-corr[2]^2)*z.sn[,2]+corr[2]*z.sn[,p]
      x[,3]=sqrt(1-corr[3]^2)*z.sn[,3]+corr[3]*z.sn[,p]
      x[,4]=sqrt(1-corr[4]^2)*z.sn[,4]+corr[4]*z.sn[,p]
      x[,5]=sqrt(1-corr[1]^2)*z.sn[,5]+corr[1]*z.sn[,p]
      x[,6]=sqrt(1-corr[2]^2)*z.sn[,6]+corr[2]*z.sn[,p]
      x[,7]=sqrt(1-corr[3]^2)*z.sn[,7]+corr[3]*z.sn[,p]
      x[,8]=sqrt(1-corr[4]^2)*z.sn[,8]+corr[4]*z.sn[,p]
      x[,9]=sqrt(1-corr[3]^2)*z.sn[,9]+corr[3]*z.sn[,p]
      x[,10]=sqrt(1-corr[4]^2)*z.sn[,10]+corr[4]*z.sn[,p]
      
      x=scale(x,center = TRUE,scale = TRUE)    #standardizing x-matrix
      ev=eigen(cor(x))
      e.val=ev$values
      e.vec=ev$vectors
      beta=e.vec[,which.max(e.val)]
      c=t(x)%*%x
      D=eigen(c)$vectors
      Z=x%*%D
      lam=t(Z)%*%Z
      lam=round(lam,4)
      lamda=diag(lam)
      alpha=t(D)%*%beta
      
      #Defining Vectors
      K1=rep(0,N)
      K2=rep(0,N)
      K3=rep(0,N)
      K4=rep(0,N)
      K5=rep(0,N)
      K6=rep(0,N)
      K7=rep(0,N)
      K8=rep(0,N)
      K9=rep(0,N)
      K10=rep(0,N)
      K11=rep(0,N)
      K12=rep(0,N)
      K13=rep(0,N)
      K14=rep(0,N)
      
      #Defining matrices
      alphahat.ols=matrix(0,nrow = p,ncol = N)
      alphahat.AGH=matrix(0,nrow = p,ncol = N)
      alphahat.AMGH=matrix(0,nrow = p,ncol = N)
      alphahat.D=matrix(0,nrow = p,ncol = N) 
      alphahat.KBG.HC4m.IR=matrix(0,nrow = p,ncol = N)
      alphahat.KBG.HC4m.IR2=matrix(0,nrow = p,ncol = N)
      alphahat.HSL.IR_scale=matrix(0,nrow = p,ncol = N)
      alphahat.LW.HC2=matrix(0,nrow = p,ncol = N)
      alphahat.LW.HC3=matrix(0,nrow = p,ncol = N)
      alphahat.LW.HC4=matrix(0,nrow = p,ncol = N)
      alphahat.LW.HC4m=matrix(0,nrow = p,ncol = N)
      alphahat.KBG=matrix(0,nrow = p,ncol = N)
      alphahat.KBG.HC2=matrix(0,nrow = p,ncol = N)
      alphahat.KBG.HC3=matrix(0,nrow = p,ncol = N)
      alphahat.KBG.HC4=matrix(0,nrow = p,ncol = N)
      alphahat.KBG.HC4m=matrix(0,nrow = p,ncol = N)
      
      #For Generating Heteroscedasticity
      a.hetro=0.0692      
      
      #For low Low levels of hetro=0.0260,Med=0.0692 and High=0.08432.
      #Note that the researchers needs to choose the a.hetro value so that delta.h close to 4 (low), 36 (medium) and 100 (high)
      
      x.row=apply(x, 1, sum)
      sigma.i=exp(a.hetro*x.row)
      delta.h=max(sigma.i)/min(sigma.i)
      
      #Simulation loop starts here
      for (i in 1:N) {
        
        e.hetro=rep(0,n)
        for(j in 1:n){
          e.hetro[j]=rnorm(1,0,sigma.i[j])  #researcher can consider other error distributions here
        }
        e=e.hetro
        y=Z%*%alpha+e
        
        beta.hat=solve(t(x)%*%x)%*%t(x)%*%y
        y.hat=x%*%beta.hat
        e.res=y-y.hat
        sigma.hat=sum((e.res)^2)/(n-p)
        alpha.hat=solve(lam)%*%t(Z)%*%y
        
        #OLS estimator
        alphahat.ols[,i]=alpha.hat
        
        #Alkhamisi (2012) generalized harmonic hetroscedastic Ridge Estimators 
        AGH=(p^2)/(sum(abs(alpha.hat)))^2        
        K1[i]=AGH
        T=solve(lam+I*K1[i])%*%lam
        alphahat.AGH[,i]=T%*%alphahat.ols[,i]
        
        #Alkhamisi (2012) Modified generalized Harmonic hetroscedastic Ridge Estimators 
        fact=((n-p)*(p+2))/(n+2)
        AMGH=(fact*(1/sum(abs(alpha.hat))))^2       
        K2[i]=AMGH
        T=solve(lam+I*K2[i])%*%lam
        alphahat.AMGH[,i]=T%*%alphahat.ols[,i]
        
        #Dorugade (2016) hetroscedastic Ridge Estimators 
        D=sqrt(sigma.hat)
        K3[i]=D
        T=solve(lam+I*K3[i])%*%lam
        alphahat.D[,i]=T%*%alphahat.ols[,i]
        
        #Dar and Chand (2023) Hetroscedastic RR Estimators
        P=solve(c)%*%t(x)
        e.res=c(e.res)
        omga=diag((e.res)^2)
        H=x%*%P
        Hjj=diag(H)
        E2=diag(1/(1-Hjj))
        E3=diag(1/(1-Hjj)^2)
        u1=min(4,(n*Hjj/p))
        E4=diag(1/(1-Hjj)^u1)
        u2=min(1.0,(n*Hjj/p))+min(1.5,(n*Hjj/p))
        E4m=diag(1/(1-Hjj)^u2)
        
        W=t(x)%*%x
        
        # Maximum of Diagonal Values
        HC4m.IR=P%*%E4m%*%omga%*%t(P)
        Sci_1=HC4m.IR%*%W
        sigma.hatro=max(diag(Sci_1))
        
        
        #Hetroscedastic HSL 1976, Ridge estimator 
        KBG.HC4m.IR=sigma.hatro/(prod(alpha.hat^2))^(1/p)
        K4[i]=KBG.HC4m.IR
        T.hatro=solve(lam+I*K4[i])%*%lam
        alphahat.KBG.HC4m.IR[,i]=T.hatro%*%alphahat.ols[,i]
        
        
        # Mean of Diagonal Values
        sigma.hatros=mean(diag(Sci_1))
        
        
        #Hetroscedastic HSL 1976, Ridge estimator 
        KBG.HC4m.IR2=sigma.hatros/(prod(alpha.hat^2))^(1/p)
        K5[i]=KBG.HC4m.IR2
        T.hatro=solve(lam+I*K5[i])%*%lam
        alphahat.KBG.HC4m.IR2[,i]=T.hatro%*%alphahat.ols[,i]
        
        #DAR (2023), Scaling Factor, Ridge estimator
        tau=sqrt(p*(abs(max(e.res))/abs(min(e.res))))
        sigma.hatro_sc=tau*sigma.hat
        
        #Dar (2023) Hetroscedastic HSL, Ridge estimator
        HSL.IR_scale= (sigma.hatro_sc)* sum((e.val*alpha.hat)^2)/(sum(e.val*alpha.hat^2))^2
        K6[i]=HSL.IR_scale
        T.hatro=solve(lam+I*K6[i])%*%lam
        alphahat.HSL.IR_scale[,i]=T.hatro%*%alphahat.ols[,i]
        
        ######Proposed Hetroscedastic RR estimators########
        
        # Proposed Estimators based on HC2
        HC2=P%*%E2%*%omga%*%t(P)
        Sci_2=HC2%*%W
        sigma.hatro1=sum(diag(Sci_2))
        
        # Proposed based on Hetroscedastic LW 1976, Ridge estimator
        LW.HC2= (p*sigma.hatro1)/(sum(e.val*(alpha.hat)^2))
        K7[i]=LW.HC2
        T.hatro=solve(lam+I*K7[i])%*%lam
        alphahat.LW.HC2[,i]=T.hatro%*%alphahat.ols[,i]
        
        # Proposed Hetroscedastic based on Kibria 2003 GM, Ridge estimator
        KBG.HC2= sigma.hatro1/(prod(alpha.hat^2))^(1/p)
        K8[i]=KBG.HC2
        T.hatro=solve(lam+I*K8[i])%*%lam
        alphahat.KBG.HC2[,i]=T.hatro%*%alphahat.ols[,i]
        
        
        # Proposed Estimators based on HC3
        HC3=P%*%E3%*%omga%*%t(P)
        Sci_3=HC3%*%W
        sigma.hatro2=sum(diag(Sci_3))
        
        # Proposed Hetroscedastic based on LW 1976, Ridge estimator
        LW.HC3= (p*sigma.hatro2)/(sum(e.val*(alpha.hat)^2))
        K9[i]=LW.HC3
        T.hatro=solve(lam+I*K9[i])%*%lam
        alphahat.LW.HC3[,i]=T.hatro%*%alphahat.ols[,i]
        
        # Proposed Hetroscedastic based on Kibria 2003 GM, Ridge estimator
        KBG.HC3= sigma.hatro2/(prod(alpha.hat^2))^(1/p)
        K10[i]=KBG.HC3
        T.hatro=solve(lam+I*K10[i])%*%lam
        alphahat.KBG.HC3[,i]=T.hatro%*%alphahat.ols[,i]
        
        
        # Proposed Estimators based on HC4
        HC4=P%*%E4%*%omga%*%t(P)
        Sci_4=HC4%*%W
        sigma.hatro3=sum(diag(Sci_4))
        
        # Proposed Hetroscedastic based on LW 1976, Ridge estimator
        LW.HC4= (p*sigma.hatro3)/(sum(e.val*(alpha.hat)^2))
        K11[i]=LW.HC4
        T.hatro=solve(lam+I*K11[i])%*%lam
        alphahat.LW.HC4[,i]=T.hatro%*%alphahat.ols[,i]
        
        
        # Proposed Hetroscedastic based on Kibria 2003 GM, Ridge estimator
        KBG.HC4= sigma.hatro3/(prod(alpha.hat^2))^(1/p)
        K12[i]=KBG.HC4
        T.hatro=solve(lam+I*K12[i])%*%lam
        alphahat.KBG.HC4[,i]=T.hatro%*%alphahat.ols[,i]
        
        
        # Proposed Estimators based on HC4m
        HC4m=P%*%E4m%*%omga%*%t(P)
        Sci_4m=HC4m%*%W
        sigma.hatro4=sum(diag(Sci_4m))
        
        # Proposed Hetroscedastic based on LW 1976, Ridge estimator
        LW.HC4m= (p*sigma.hatro4)/(sum(e.val*(alpha.hat)^2))
        K13[i]=LW.HC4m
        T.hatro=solve(lam+I*K13[i])%*%lam
        alphahat.LW.HC4m[,i]=T.hatro%*%alphahat.ols[,i]
        
        # Proposed Hetroscedastic based on Kibria 2003 GM, Ridge estimator
        KBG.HC4m= sigma.hatro4/(prod(alpha.hat^2))^(1/p)
        K14[i]=KBG.HC4m
        T.hatro=solve(lam+I*K14[i])%*%lam
        alphahat.KBG.HC4m[,i]=T.hatro%*%alphahat.ols[,i]
        
      } #End of simulation loop for N=5000 runs.
      
      MSE.ols=sum((alphahat.ols-c(alpha))^2)/N
      MSE.AGH=sum((alphahat.AGH-c(alpha))^2)/N
      MSE.AMGH=sum((alphahat.AMGH-c(alpha))^2)/N
      MSE.D=sum((alphahat.D-c(alpha))^2)/N
      MSE.KBG.HC4m.IR2=sum((alphahat.KBG.HC4m.IR2-c(alpha))^2)/N
      MSE.KBG.HC4m.IR=sum((alphahat.KBG.HC4m.IR-c(alpha))^2)/N
      MSE.HSL.IR_scale=sum((alphahat.HSL.IR_scale-c(alpha))^2)/N
      MSE.LW.HC2=sum((alphahat.LW.HC2-c(alpha))^2)/N
      MSE.LW.HC3=sum((alphahat.LW.HC3-c(alpha))^2)/N
      MSE.LW.HC4=sum((alphahat.LW.HC4-c(alpha))^2)/N
      MSE.LW.HC4m=sum((alphahat.LW.HC4m-c(alpha))^2)/N
      MSE.KBG.HC2=sum((alphahat.KBG.HC2-c(alpha))^2)/N
      MSE.KBG.HC3=sum((alphahat.KBG.HC3-c(alpha))^2)/N
      MSE.KBG.HC4=sum((alphahat.KBG.HC4-c(alpha))^2)/N
      MSE.KBG.HC4m=sum((alphahat.KBG.HC4m-c(alpha))^2)/N
      
      MSE=c(MSE.ols,MSE.AGH,MSE.AMGH,MSE.D,
            MSE.KBG.HC4m.IR2,MSE.KBG.HC4m.IR,MSE.HSL.IR_scale,
            MSE.LW.HC2,MSE.LW.HC3,MSE.LW.HC4,MSE.LW.HC4m,
            MSE.KBG.HC2,MSE.KBG.HC3,MSE.KBG.HC4,MSE.KBG.HC4m)
      
      MSE=round(MSE,4)
      as.matrix(MSE)
      MSE.r[,f]=MSE
      
    } #End of correlation loop
    
    MSE.n[,,a]=MSE.r
  } #End of sample size
  
  MSE.p[,,,b]=MSE.n
} #End of predictors loop

col.names=c("Repl1","Repl2","Repl3","Repl4")
row.names=c("OLS","GH","MGH","DRR","DAR_HC4m","DAR_MHC4m","DAR_HA",
            "LW_NHC2","LW_NHC3","LW_NHC4","LW_NHC4m",
            "KBG_NHC2","KBG_NHC3","KBG_NHC4","KBG_NHC4m")

matrix.names2=c("20","50","100")
matrix.names3=c("10")
dimnames(MSE.p)=list(row.names,col.names,matrix.names2,matrix.names3)
write.csv(MSE.p,file = "Enter the path to save the output")

delta.h           #this term needs to be 4, 36 and 100 for incorporating required levels of hetroscedasticty

