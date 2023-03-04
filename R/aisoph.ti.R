aisoph.ti=function(df.full, Q, shape1, shape2, K1, K2, maxiter, eps){
  #1. sorted by z
  n.full=nrow(df.full)
  z1.full=sort(df.full$z1)
  z2.full=sort(df.full$z2)
  df=df.full
  
  #2. remove subjects whose cov is less than z^*_(1) #result in -Inf
  delta.z1=aggregate(df$status,by=list(z1=df$z1),sum)
  if(shape1=='increasing'){
    z1.star=delta.z1[which(delta.z1[,2]>0)[1],1]
  }else if(shape1=='decreasing'){
    z1.star=delta.z1[max(which(delta.z1[,2]>0)),1]
  }
  
  delta.z2=aggregate(df$status,by=list(z2=df$z2),sum)
  if(shape2=='increasing'){
    z2.star=delta.z2[which(delta.z2[,2]>0)[1],1]
  }else if(shape2=='decreasing'){
    z2.star=delta.z2[max(which(delta.z2[,2]>0)),1]
  }
  
  if(shape1=="increasing"){;       df=df[which(df$z1>=z1.star),]
  }else if(shape1=="decreasing"){; df=df[which(df$z1<=z1.star),]
  }  

  if(shape2=="increasing"){;       df=df[which(df$z2>=z2.star),]
  }else if(shape2=="decreasing"){; df=df[which(df$z2<=z2.star),]
  }

  #3. redefine variables
  n=nrow(df);
  df1=df[order(df$z1),];     df1$id1=1:n #id is used for sorting
  df2=df1[order(df1$z2),];   df2$id2=1:n
  df1=df2[order(df2$z1),];   
  
  id1=df1$id2;                id2=df2$id1
  status1=df1$status;         status2=df2$status
  
  t1=df1$time;                t2=df2$time;                    
  t.obs=sort(unique(t1));
  nt=length(t.obs);       
  
  z1=df1$z1;                      z2=df2$z2
  z1.obs=unique(z1[status1==1]);  z2.obs=unique(z2[status2==1])  
  m1=length(z1.obs);              m2=length(z2.obs)    

  if(Q==0){
    x1=x2=matrix(rep(0,n),ncol=1)
  }else{
    x1=as.matrix(df1[,!(colnames(df1)%in%c("time","status","z1","z2","id1","id2"))])
    x2=as.matrix(df2[,!(colnames(df2)%in%c("time","status","z1","z2","id1","id2"))])
  }
  
  #anchor
  k1=sum(z1.obs<K1);              k2=sum(z2.obs<K2)
  if(k1==0) k1=1;                 if(k2==0) k2=1
  zk1=z1.obs[k1];                 zk2=z2.obs[k2]   
  
  #counting process
  Y1=Y2=dN1=dN2=matrix(0,n,nt)       #row is the subj, col is the time corresponding z_(i);
  for(i in 1:n){
    rank.t1=which(t1[i]==t.obs)
    Y1[i,][1:rank.t1]=1
    if(status1[i]==1) dN1[i,][rank.t1]=1
  
    rank.t2=which(t2[i]==t.obs)
    Y2[i,][1:rank.t2]=1
    if(status2[i]==1) dN2[i,][rank.t2]=1    
  }
  
  #initial value
  int.val=initial.ft(df,Q,z1.obs,zk1,z2.obs,zk2,m1,m2,shape1,shape2)
  psi1=int.val$psi1
  psi2=int.val$psi2
  beta=int.val$beta
  
  #upper / lower bounds
  #1/a<HR<1; log(1/a)<HR<log(a)
  #lb=lw2=log(1/10)
  #ub=up2=log(10)
  
  #cycling picm & Newton Raphson algorithm
  iter=0;  dist=1;
  psi2.full=BTF.ft(m=m2,n=n,z=z2,z.obs=z2.obs,psi.new=psi2,shape=shape2)
  psi2.full1=psi2.full[id2] #matched to df1
  #psi2.full[which(psi2.full>up2)]=up2
  #psi2.full[which(psi2.full<lw2)]=lw2
  
  if(Q==0){
    while(dist>=eps){
      iter=iter+1
      if(iter>maxiter) break    
      
      #estimate psi1
      psi1.new=cpicm.ft(psi1,psi2.full1,z1,z1.obs,m1,k1,n,Y1,dN1,nt, x1,beta ,eps,maxiter, shape1)
      #  psi1.new[which(psi1.new>up1)]=up1
      #  psi1.new[which(psi1.new<lw1)]=lw1
      psi1.full=BTF.ft(m=m1,n=n,z=z1,z.obs=z1.obs,psi.new=psi1.new,shape=shape1)
      psi1.full2=psi1.full[id1] #matched to df2
      
      #estimate psi2;
      psi2.new=cpicm.ft(psi2,psi1.full2,z2,z2.obs,m2,k2,n,Y2,dN2,nt, x2,beta ,eps,maxiter, shape2)
      #  psi2.new[which(psi2.new>up2)]=up2
      #  psi2.new[which(psi2.new<lw2)]=lw2
      psi2.full=BTF.ft(m=m2,n=n,z=z2,z.obs=z2.obs,psi.new=psi2.new,shape=shape2)
      psi2.full1=psi2.full[id2] #matched to df1
      
      #update;
      dist=sqrt(sum((psi1.new-psi1)^2))+sqrt(sum((psi2.new-psi2)^2))
      if(is.infinite(dist) | is.na(dist))
        break
      psi1=psi1.new
      psi2=psi2.new
    }
  }else{
    while(dist>=eps){
      iter=iter+1
      if(iter>maxiter) break    
      
      #estimate psi1
      psi1.new=cpicm.ft(psi1,psi2.full1,z1,z1.obs,m1,k1,n,Y1,dN1,nt, x1,beta ,eps,maxiter, shape1)
      #  psi1.new[which(psi1.new>up1)]=up1
      #  psi1.new[which(psi1.new<lw1)]=lw1
      psi1.full=BTF.ft(m=m1,n=n,z=z1,z.obs=z1.obs,psi.new=psi1.new,shape=shape1)
      psi1.full2=psi1.full[id1] #switch psi1 according to z2
      
      #estimate psi2;
      psi2.new=cpicm.ft(psi2,psi1.full2,z2,z2.obs,m2,k2,n,Y2,dN2,nt, x2,beta ,eps,maxiter, shape2)
      #  psi2.new[which(psi2.new>up2)]=up2
      #  psi2.new[which(psi2.new<lw2)]=lw2
      psi2.full=BTF.ft(m=m2,n=n,z=z2,z.obs=z2.obs,psi.new=psi2.new,shape=shape2)
      psi2.full1=psi2.full[id2] #switch psi2 according to z1
      
      #estimate beta;
      beta.new=NR.ft(x1,beta,Q,psi1.full,psi2.full1,n,Y1,dN1, maxiter,eps) #does not matter based on df1 or df2
      
      #update;
      dist=sqrt(sum((psi1.new-psi1)^2))+sqrt(sum((psi2.new-psi2)^2))+sqrt(sum(((beta.new-beta)^2)))
      if(is.infinite(dist) | is.na(dist))
        break
      psi1=psi1.new
      psi2=psi2.new
      beta=beta.new
    }
  }
      
  #cyclic picm result
  conv="converged"
  if(is.na(dist)){; conv="not converged"
  }else if(is.nan(dist)){; conv="not converged"
  }else if(is.infinite(dist)){; conv="not converged"
  #}else if(max(diff(psi1))<eps | max(diff(psi2))<eps){; conv="not converged" #no jump
  }
  
  if(conv=="not converged")
    psi1=psi2=beta=NA
  
  psi1.full=BTF.ft(m=m1,n=n.full,z=z1.full,z.obs=z1.obs,psi.new=psi1,shape=shape1) #for full data
  psi2.full=BTF.ft(m=m2,n=n.full,z=z2.full,z.obs=z2.obs,psi.new=psi2,shape=shape2)
  
  #add removed points
  if(shape1=="increasing"){
    idx=which(df.full$z1<z1.star)
    if(length(idx)>0){
      z1.full=c(df.full$z1[idx],z1.full)
      psi1.full=c(rep(-Inf,length(idx)),psi1.full)
    }
  }else if(shape1=="decreasing"){
    idx=which(df.full$z1>z1.star)
    if(length(idx)>0){
      z1.full=c(z1.full,df.full$z1[idx])
      psi1.full=c(psi1.full,rep(-Inf,length(idx)))
    }
  }
  
  if(shape2=="increasing"){
    idx=which(df.full$z2<z2.star)
    if(length(idx)>0){
      z2.full=c(df.full$z2[idx],z2.full)
      psi2.full=c(rep(-Inf,length(idx)),psi2.full)
    }
  }else if(shape2=="decreasing"){
    idx=which(df.full$z2>z2.star)
    if(length(idx)>0){
      z2.full=c(z2.full,df.full$z2[idx])
      psi2.full=c(psi2.full,rep(-Inf,length(idx)))
    }
  }
  
  iso1=data.frame(z=z1.full,psi=psi1.full,HR=exp(psi1.full))
  iso2=data.frame(z=z2.full,psi=psi2.full,HR=exp(psi2.full))
  
  iso1$cens="yes" #potential jump points
  iso2$cens="yes"
  
  for(i in 1:nrow(iso1))
    if(iso1$z[i]%in%z1.obs) iso1$cens[i]="no"
  
  for(i in 1:nrow(iso2))
    if(iso2$z[i]%in%z2.obs) iso2$cens[i]="no"
  
  est=data.frame(beta=NA,HR=NA)
  if(Q>0){
    est=data.frame(beta=beta,HR=exp(beta.new))
    rownames(est)=paste0("x",1:Q)
  }
  
  if(conv=="not converged"){
    iso1$psi=iso1$HR=NA
    iso2$psi=iso2$HR=NA  }
  

  return(list(iso1=iso1,iso2=iso2,est=est,
              conv=conv,shape1=shape1,shape2=shape2,K1=zk1,K2=zk2))
}
