rm(list = ls()) 
source("R/stab_knock_functions.R")
source("R/ds_mds_functions.R")
library(knockoff)
library(glmnet)
library(Rfast)
library(hdi)
outdir <- "results/simulation"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
####################################################
set.seed(17482)
q = 0.1
####################################################
##############simulation_FDR&Power##################
####################################################
p = 150; n = 400; k = 30; A=0.6; rho=0.4; m=10###A=0.6,0.8; rho=0.4,0.5,0.6,0.7,0.8
mu = rep(0,p); Sigma= toeplitz(rho^(1:p-1))
nIterations = 100
r1=matrix(rep(0,nIterations*22),nrow=nIterations)
colnames(r1) = c(
  "Knockoff_fdp", "Knockoff_power",
  "Lasso_fdp", "Lasso_power",
  "SIS_fdp", "SIS_power",
  "StabKnock_fdp", "StabKnock_power",
  "StabEN_fdp", "StabEN_power",
  "StabAL_fdp", "StabAL_power",
  "StabLasso_fdp", "StabLasso_power",
  "DS_fdp", "DS_power",
  "MDS_fdp", "MDS_power",
  "BH_fdp", "BH_power",
  "MBH_fdp", "MBH_power"
)
for(j in 1:nIterations){
  # Generate the Data
  X = matrix(rnorm(n*p),n) %*% chol(Sigma)
  # X = rmvt(n, mu=mu, sigma = Sigma, v = 3)
  nonzero = sample(p, k)
  beta.true = rep(0, p)
  beta.true[nonzero] = sample(c(A,-A),k,replace=T)
  y = X %*% beta.true+ rnorm(n)   
  cv.fit=cv.glmnet(x=X,y=y, type.measure = "mse", nfolds = 10)
  lambda=cv.fit$lambda.min
  ### Knockoff
  knockoffs = function(X) create.fixed(X, method='equi')
  result_X_1 = knockoff.filter(X = X, y = y, knockoffs=knockoffs, statistic=stat.lasso_lambdasmax, fdr= q)
  fdp_power1 = fdp_power(selected_index=result_X_1$selected,signal_index=nonzero)
  ### Lasso
  lasso_mod=glmnet(x=X,y=y,alpha=1,lambda=lambda)
  nonzero_lasso=which(lasso_mod$beta!=0)
  fdp_power2 = fdp_power(selected_index=nonzero_lasso,signal_index=nonzero)
  ### SIS
  c=0.1
  correlation=cor(X,y)
  nonzero_SIS=which(abs(correlation)>c)
  fdp_power3 <- fdp_power(selected_index=nonzero_SIS, signal_index = nonzero)
  ### Stab-Knock
  SPD_result=Knock_SPD(X=X,y=y,L=10,q=q,cv=1)
  nonzero_StabKnock=SPD_result$S_plus
  fdp_power4 <- fdp_power(selected_index=nonzero_StabKnock,signal_index=nonzero)
  ## Stab-EN
  L_EN=50
  SP_EN=rep(0,p)
  for(i in 1:L_EN){
    index=sample(n,size =floor(n/2),replace = F )
    lasso_mod=glmnet(x=X[index,],y=y[index],alpha=0.5,lambda=lambda)
    SP_EN[which(lasso_mod$beta!=0)]=SP_EN[which(lasso_mod$beta!=0)]+1
  }
  SP_EN=SP_EN/L_EN
  fdp_power5 <- fdp_power(selected_index=which(SP_EN>0.6),signal_index=nonzero)
  ## Stab-AL
  L_AL=50
  SP_AL=rep(0,p)
  ada_weight=1/abs(as.vector(lasso_mod$beta))
  for(i in 1:L_AL){
    index=sample(n,size =floor(n/2),replace = F )
    fit.AL=glmnet(x=X[index,],y=y[index],penalty.factor = ada_weight,lambda=lambda)
    SP_AL[which(fit.AL$beta!=0)]=SP_AL[which(fit.AL$beta!=0)]+1
  }
  SP_AL=SP_AL/L_AL
  fdp_power6 <- fdp_power(selected_index=which(SP_AL>0.6),signal_index=nonzero)
  ##Stab-Lasso
  L_Lasso=50
  SP_Lasso=rep(0,p)
  for(i in 1:L_Lasso){
    index=sample(n,size =floor(n/2),replace = F )
    fit.Lasso=glmnet(x=X[index,],y=y[index],lambda=lambda)
    SP_Lasso[which(fit.Lasso$beta!=0)]=SP_Lasso[which(fit.Lasso$beta!=0)]+1
  }
  SP_Lasso=SP_Lasso/L_Lasso
  fdp_power7 <- fdp_power(selected_index=which(SP_Lasso>0.6),signal_index=nonzero)
  ### DS
  DS_result = DS(X = X, y = y, num_split = 1, q = q)
  fdp_power8 <- fdp_power(selected_index=DS_result$DS_feature,signal_index=nonzero)
  ### MDS
  MDS_result = DS(X = X, y = y, num_split = 20, q = q)
  fdp_power9 <- fdp_power(selected_index=MDS_result$MDS_feature,signal_index=nonzero)
  ### BH
  adh = (p.adjust(univglms(y=y,x=X, oiko='normal')[,2],method="BH") < q)
  nonzero_BH = which(as.vector(adh)==TRUE)
  fdp_power10 <- fdp_power(selected_index=nonzero_BH,signal_index=nonzero)
  ### MBH
  fit2 <- multi.split(X, y, B=20)
  nonzero_MBH <- which(fit2$pval.corr < q)
  fdp_power11 <- fdp_power(selected_index=nonzero_MBH,signal_index=nonzero)
  r1[j,] = c(
    fdp_power1$fdp, fdp_power1$power,
    fdp_power2$fdp, fdp_power2$power,
    fdp_power3$fdp, fdp_power3$power,
    fdp_power4$fdp, fdp_power4$power,
    fdp_power5$fdp, fdp_power5$power,
    fdp_power6$fdp, fdp_power6$power,
    fdp_power7$fdp, fdp_power7$power,
    fdp_power8$fdp, fdp_power8$power,
    fdp_power9$fdp, fdp_power9$power,
    fdp_power10$fdp, fdp_power10$power,
    fdp_power11$fdp, fdp_power11$power
  )
}
colMeans(r1)
save.image(file.path(outdir, 'FDR&Power_g_rho04_A06.RData'))
####################################################
###############simulation_AIC&MS####################
####################################################
set.seed(17482)
nIterations = 100
p = 100; p1 = 20; A=0.8; rho=0.2; m=10
mu = rep(0,p); Sigma= toeplitz(rho^(1:p-1))
r1=matrix(rep(0,nIterations*14),nrow=nIterations)
k=nIterations
q=0.1
N=c(200,300,400,500)
l=length(N)
AIC_S_knock=matrix(rep(0,l*k),nrow=k)
AIC_S_Lasso=matrix(rep(0,l*k),nrow=k)
AIC_S_SIS=matrix(rep(0,l*k),nrow=k)
AIC_S_StabKnock=matrix(rep(0,l*k),nrow=k)
AIC_S_EN=matrix(rep(0,l*k),nrow=k)
AIC_S_AL=matrix(rep(0,l*k),nrow=k)
AIC_S_StabLasso=matrix(rep(0,l*k),nrow=k)
MS_S_knock=matrix(rep(0,l*k),nrow=k)
MS_S_Lasso=matrix(rep(0,l*k),nrow=k)
MS_S_SIS=matrix(rep(0,l*k),nrow=k)
MS_S_StabKnock=matrix(rep(0,l*k),nrow=k)
MS_S_EN=matrix(rep(0,l*k),nrow=k)
MS_S_AL=matrix(rep(0,l*k),nrow=k)
MS_S_StabLasso=matrix(rep(0,l*k),nrow=k)
AIC_S_DS=matrix(rep(0,l*k),nrow=k)
AIC_S_MDS=matrix(rep(0,l*k),nrow=k)
AIC_S_BH=matrix(rep(0,l*k),nrow=k)
AIC_S_MBH=matrix(rep(0,l*k),nrow=k)
MS_S_DS=matrix(rep(0,l*k),nrow=k)
MS_S_MDS=matrix(rep(0,l*k),nrow=k)
MS_S_BH=matrix(rep(0,l*k),nrow=k)
MS_S_MBH=matrix(rep(0,l*k),nrow=k)
for(j in 1:k){
  for(i in 1:l){
    n=N[i]
    # Generate the Data
    X = matrix(rnorm(n*p),n) %*% chol(Sigma)
    # X = rmvt(n, mu=mu, sigma = Sigma, v = 4)
    nonzero = sample(p, p1)
    beta.true = rep(0, p)
    beta.true[nonzero] = sample(c(A,-A),p1,replace=T)
    y = X %*% beta.true+ rnorm(n)  
    cv.fit=cv.glmnet(x=X,y=y, type.measure = "mse", nfolds = 10)
    lambda=cv.fit$lambda.min
    ### Knockoff
     if(n<2*p){
      index_knock=sort(order(cor(X,y),decreasing=T)[1:(n/2-1)])
      knockoffs = function(X) create.fixed(X, method='equi')
      result_X_1 = knockoff.filter(X = X[,index_knock], y = y, knockoffs=knockoffs, statistic=stat.lasso_lambdasmax, fdr= q)
      MS_S_knock[j,i]=length(result_X_1$selected)
      if(length(result_X_1$selected)>0){
        X_Knock=X[,index_knock][,result_X_1$selected]
        Knock_mod=glmnet(x=X_Knock,y=y,alpha=1,lambda=lambda)
        pred_Knock=predict(Knock_mod,newx = X_Knock)
        AIC_S_knock[j,i]=AIC(y,pred_Knock,n=N[i],d=length(result_X_1$selected))
      }else{
        AIC_S_knock[j,i]=NA
      }
    }else{
      knockoffs = function(X) create.fixed(X, method='equi')
      result_X_1 = knockoff.filter(X = X, y = y, knockoffs=knockoffs, statistic=stat.lasso_lambdasmax, fdr= q)
      MS_S_knock[j,i]=length(result_X_1$selected)
      if(length(result_X_1$selected)>0){
        X_Knock=X[,result_X_1$selected]
        cv.fit=cv.glmnet(x=X_Knock,y=y, type.measure = "mse", nfolds = 10)
        Knock_mod=glmnet(x=X_Knock,y=y,alpha=1,lambda=cv.fit$lambda.min)
        pred_Knock=predict(Knock_mod,newx = X[,result_X_1$selected])
        AIC_S_knock[j,i]=AIC(y,pred_Knock,n=N[i],d=length(result_X_1$selected))
      }else{
        AIC_S_knock[j,i]=NA
      }
    }
    ### Lasso
    lasso_mod=glmnet(x=X,y=y,alpha=1,lambda=lambda)
    nonzero_lasso=which(lasso_mod$beta!=0)
    MS_S_Lasso[j,i]=length(nonzero_lasso)
    if(length(nonzero_lasso)>0){
      X_Lasso=X[,nonzero_lasso]
      Lasso_mod=glmnet(x=X_Lasso,y=y,alpha=1,lambda=lambda)
      pred_Lasso=predict(Lasso_mod,newx = X_Lasso)
      AIC_S_Lasso[j,i]=AIC(y,pred_Lasso,n=N[i],d=length(nonzero_lasso))
    }else{
      AIC_S_Lasso[j,i]=NA
    }
    ### SIS
    c=0.1
    correlation=cor(X,y)
    nonzero_SIS=which(abs(correlation)>c)
    MS_S_SIS[j,i]=length(nonzero_SIS)
    X_SIS=X[,nonzero_SIS]
    SIS_mod=glmnet(x=X_SIS,y=y,alpha=1,lambda=lambda)
    pred_SIS=predict(SIS_mod,newx = X_SIS)
    AIC_S_SIS[j,i]=AIC(y,pred_SIS,n=N[i],d=length(nonzero_SIS))
    ### Stab-Knock
    if(n<2*p){
    SPD_result=Knock_SPD(X=X[,index_knock],y=y,L=10,q=q,cv=1)
    nonzero_StabKnock=SPD_result$S_plus
    MS_S_StabKnock[j,i]=length(nonzero_StabKnock)
    if(length(nonzero_StabKnock)>0){
      X_StabKnock=X[,index_knock][,nonzero_StabKnock]
      StabKnock_mod=glmnet(x=X_StabKnock,y=y,alpha=1,lambda=lambda)
      pred_StabKnock=predict(StabKnock_mod,newx = X_StabKnock)
      AIC_S_StabKnock[j,i]=AIC(y,pred_StabKnock,n=N[i],d=length(nonzero_StabKnock))
    }else{AIC_S_StabKnock[j,i]=NA}
    }else{
      SPD_result=Knock_SPD(X=X,y=y,L=10,q=q,cv=1)
      nonzero_StabKnock=SPD_result$S_plus
      MS_S_StabKnock[j,i]=length(nonzero_StabKnock)
      if(length(nonzero_StabKnock)>0){
        X_StabKnock=X[,nonzero_StabKnock]
        StabKnock_mod=glmnet(x=X_StabKnock,y=y,alpha=1,lambda=lambda)
        pred_StabKnock=predict(StabKnock_mod,newx = X_StabKnock)
        AIC_S_StabKnock[j,i]=AIC(y,pred_StabKnock,n=N[i],d=length(nonzero_StabKnock))
      }else{AIC_S_StabKnock[j,i]=NA}
    }
    ## Stab-EN
    L_EN=20
    SP_EN=rep(0,p)
    for(kk in 1:L_EN){
      index=sample(n,size =floor(n/2),replace = F )
      lasso_mod=glmnet(x=X[index,],y=y[index],alpha=0.5,lambda=lambda)
      SP_EN[which(lasso_mod$beta!=0)]=SP_EN[which(lasso_mod$beta!=0)]+1
    }
    SP_EN=SP_EN/L_EN
    MS_S_EN[j,i]=length(which(SP_EN>0.8))
    if(length(which(SP_EN>0.8))>0){
      X_EN=X[,which(SP_EN>0.8)]
      EN_mod=glmnet(x=X_EN,y=y,alpha=1,lambda=lambda)
      pred_EN=predict(EN_mod,newx = X_EN)
      AIC_S_EN[j,i]=AIC(y,pred_EN,n=N[i],d=length(which(SP_EN>0.8)))
    }else{AIC_S_EN[j,i]=NA}
    ## Stab-AL
    L_AL=20
    SP_AL=rep(0,p)
    ada_weight=1/abs(as.vector(lasso_mod$beta))
    for(kk in 1:L_AL){
      index=sample(n,size =floor(n/2),replace = F )
      fit.AL=glmnet(x=X[index,],y=y[index],penalty.factor = ada_weight,lambda=lambda)
      SP_AL[which(fit.AL$beta!=0)]=SP_AL[which(fit.AL$beta!=0)]+1
    }
    SP_AL=SP_AL/L_AL
    MS_S_AL[j,i]=length(which(SP_AL>0.8))
    if(length(which(SP_AL>0.8))>0){
      X_AL=X[,which(SP_AL>0.8)]
      AL_mod=glmnet(x=X_AL,y=y,alpha=1,lambda=lambda)
      pred_AL=predict(AL_mod,newx = X_AL)
      AIC_S_AL[j,i]=AIC(y,pred_AL,n=N[i],d=length(which(SP_AL>0.8)))
    }else{AIC_S_AL[j,i]=NA}
    ##Stab-Lasso
    L_Lasso=20
    SP_Lasso=rep(0,p)
    for(kk in 1:L_Lasso){
      index=sample(n,size =floor(n/2),replace = F )
      fit.Lasso=glmnet(x=X[index,],y=y[index],lambda=lambda)
      SP_Lasso[which(fit.Lasso$beta!=0)]=SP_Lasso[which(fit.Lasso$beta!=0)]+1
    }
    SP_Lasso=SP_Lasso/L_Lasso
    MS_S_StabLasso[j,i]=length(which(SP_Lasso>0.8))
    if(length(which(SP_Lasso>0.8))>0){
      X_StabLasso=X[,which(SP_Lasso>0.8)]
      StabLasso_mod=glmnet(x=X_StabLasso,y=y,alpha=1,lambda=lambda)
      pred_StabLasso=predict(StabLasso_mod,newx = X_StabLasso)
      AIC_S_StabLasso[j,i]=AIC(y,pred_StabLasso,n=N[i],d=length(which(SP_Lasso>0.8)))
    }else{AIC_S_StabLasso[j,i]=NA}
    ### DS
    DS_result=DS(X = X, y = y, num_split = 1, q = q)
    nonzero_DS=DS_result$DS_feature
    MS_S_DS[j,i]=length(nonzero_DS)
    if(length(nonzero_DS)>0){
      X_DS=X[,nonzero_DS]
      DS_mod=glmnet(x=X_DS,y=y,alpha=1,lambda=lambda)
      pred_DS=predict(DS_mod,newx = X_DS)
      AIC_S_DS[j,i]=AIC(y,pred_DS,n=N[i],d=length(nonzero_DS))
    }else{AIC_S_DS[j,i]=NA}
    ### MDS
    MDS_result=DS(X = X, y = y, num_split = 20, q = q)
    nonzero_MDS=MDS_result$MDS_feature
    MS_S_MDS[j,i]=length(nonzero_MDS)
    if(length(nonzero_MDS)>0){
      X_MDS=X[,nonzero_MDS]
      MDS_mod=glmnet(x=X_MDS,y=y,alpha=1,lambda=lambda)
      pred_MDS=predict(MDS_mod,newx = X_MDS)
      AIC_S_MDS[j,i]=AIC(y,pred_MDS,n=N[i],d=length(nonzero_MDS))
    }else{AIC_S_MDS[j,i]=NA}
    ### BH
    adh = (p.adjust(univglms(y=y,x=X, oiko='normal')[,2],method="BH") < q)
    nonzero_BH=which(as.vector(adh)==TRUE)
    MS_S_BH[j,i]=length(nonzero_BH)
    if(length(nonzero_BH)>0){
      X_BH=X[,nonzero_BH]
      BH_mod=glmnet(x=X_BH,y=y,alpha=1,lambda=lambda)
      pred_BH=predict(BH_mod,newx = X_BH)
      AIC_S_BH[j,i]=AIC(y,pred_BH,n=N[i],d=length(nonzero_BH))
    }else{AIC_S_BH[j,i]=NA}
    ### MBH
    fit2 <- multi.split(X, y, B=20)
    nonzero_MBH <- which(fit2$pval.corr < q)
    MS_S_MBH[j,i]=length(nonzero_MBH)
    if(length(nonzero_MBH)>0){
      X_MBH=X[,nonzero_MBH]
      MBH_mod=glmnet(x=X_MBH,y=y,alpha=1,lambda=lambda)
      pred_MBH=predict(MBH_mod,newx = X_MBH)
      AIC_S_MBH[j,i]=AIC(y,pred_MBH,n=N[i],d=length(nonzero_MBH))
    }else{AIC_S_MBH[j,i]=NA}
  }}
save.image(file.path(outdir, "AIC_simu_p100_k20_A08_rho02_g.RData"))
save.image(file.path(outdir, "MS_simu_p100_k20_A08_rho02_g.RData"))
####################################################
########simulation_different_threshold_SS###########
####################################################
p = 150; n = 400; k = 30; A=0.3; rho=0.8; m=10;thre=0.1##thre=0.1,0.3,0.5,0.7,0.9
mu = rep(0,p); Sigma= toeplitz(rho^(1:p-1))
nIterations = 100
r1=matrix(rep(0,nIterations*6),nrow=nIterations)
for(j in 1:nIterations){
  # Generate the Data
  X = matrix(rnorm(n*p),n) %*% chol(Sigma)
  # X = rmvt(n, mu=mu, sigma = Sigma, v = 3)
  nonzero = sample(p, k)
  beta.true = rep(0, p)
  beta.true[nonzero] = sample(c(A,-A),k,replace=T)
  y = X %*% beta.true+ rnorm(n)   
  cv.fit=cv.glmnet(x=X,y=y, type.measure = "mse", nfolds = 10)
  lambda=cv.fit$lambda.min
  ## Stab-EN
  L_EN=50
  SP_EN=rep(0,p)
  for(i in 1:L_EN){
    index=sample(n,size =floor(n/2),replace = F )
    lasso_mod=glmnet(x=X[index,],y=y[index],alpha=0.5,lambda=lambda)
    SP_EN[which(lasso_mod$beta!=0)]=SP_EN[which(lasso_mod$beta!=0)]+1
  }
  SP_EN=SP_EN/L_EN
  fdp_power1 <- fdp_power(selected_index=which(SP_EN>thre),signal_index=nonzero)
  ## Stab-AL
  L_AL=50
  SP_AL=rep(0,p)
  ada_weight=1/abs(as.vector(lasso_mod$beta))
  for(i in 1:L_AL){
    index=sample(n,size =floor(n/2),replace = F )
    fit.AL=glmnet(x=X[index,],y=y[index],penalty.factor = ada_weight,lambda=lambda)
    SP_AL[which(fit.AL$beta!=0)]=SP_AL[which(fit.AL$beta!=0)]+1
  }
  SP_AL=SP_AL/L_AL
  fdp_power2 <- fdp_power(selected_index=which(SP_AL>thre),signal_index=nonzero)
  ##Stab-Lasso
  L_Lasso=50
  SP_Lasso=rep(0,p)
  for(i in 1:L_Lasso){
    index=sample(n,size =floor(n/2),replace = F )
    fit.Lasso=glmnet(x=X[index,],y=y[index],lambda=lambda)
    SP_Lasso[which(fit.Lasso$beta!=0)]=SP_Lasso[which(fit.Lasso$beta!=0)]+1
  }
  SP_Lasso=SP_Lasso/L_Lasso
  fdp_power3 <- fdp_power(selected_index=which(SP_Lasso>thre),signal_index=nonzero)
  r1[j,] = c(
    fdp_power1$fdp, fdp_power1$power,
    fdp_power2$fdp, fdp_power2$power,
    fdp_power3$fdp, fdp_power3$power
  )
  
}
colMeans(r1)
save.image(file.path(outdir, 'FDR&Power_diffthreshold_A03_rho08_thre01.RData'))







