initial.ft=function(df,Q,z1.obs,zk1,z2.obs,zk2,m1,m2,shape1,shape2){
  oldwarn <- options(warn=2)
  on.exit(options(oldwarn))

  beta=0
  if(Q==0){
    beta.hat=coxph(Surv(df$time,df$status)~df$z1+df$z2)$coefficient
  }else{
    x=matrix(df[,!(colnames(df)%in%c("time","status","z1","z2","id","id1","id2"))])
    beta.hat=coxph(Surv(df$time,df$status)~df$z1+df$z2+x)$coefficient
    beta=beta.hat[-c(1:2)]
  }

  if(shape1=='increasing'){
    psi1= abs(beta.hat[1])*(z1.obs-zk1)
  }else if(shape1=='decreasing'){
    psi1=-abs(beta.hat[1])*(z1.obs-zk1)
  }
    
  if(shape2=='increasing'){
    psi2= abs(beta.hat[2])*(z2.obs-zk2)
  }else if(shape2=='decreasing'){
    psi2=-abs(beta.hat[2])*(z2.obs-zk2)
  }

  return(list(psi1=psi1,psi2=psi2,beta=beta))
}
