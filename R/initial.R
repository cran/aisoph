initial.ft=function(df,Q,z1.obs,zk1,z2.obs,zk2,m1,m2,shape1,shape2){
  oldwarn <- options(warn=2)
  on.exit(options(oldwarn))
  
  fit1<-isoSurv::isoph(Surv(df$time,df$status)~iso(df$z1,shape=shape1))
  psi1=rep(NA,m1)
  for(j in 1:m1)
    psi1[j]=(fit1$iso.cov[which(z1.obs[j]==fit1$iso.cov[,1]),2])[1] #all same
  psi1[is.infinite(psi1) & psi1<0]=min(psi1[is.finite(psi1)])
  psi1[is.infinite(psi1) & psi1>0]=max(psi1[is.finite(psi1)])

  fit2<-isoSurv::isoph(Surv(df$time,df$status)~iso(df$z2,shape=shape2))
  psi2=rep(NA,m2)
  for(j in 1:m2)
    psi2[j]=(fit2$iso.cov[which(z2.obs[j]==fit2$iso.cov[,1]),2])[1] #all same
  psi2[is.infinite(psi2) & psi2<0]=min(psi2[is.finite(psi2)])
  psi2[is.infinite(psi2) & psi2>0]=max(psi2[is.finite(psi2)])

  beta=0
  if(Q>=1){
    if(Q==1){
      w=matrix(df[,-c(1:4)],ncol=1)
    }else{
      w=as.matrix(df[,-c(1:4)])
    }
    beta=coxph(Surv(df$time,df$status)~w)$coefficient
  }
    
  return(list(psi1=psi1,psi2=psi2,beta=beta))
}
