aisoph.ti=function(df.full, Q, shape1, shape2, K1, K2, maxiter, eps){
  #1. sorted by z
  n.full=nrow(df.full)
  z1.full=sort(df.full$z1)
  z2.full=sort(df.full$z2)
  df=df.full
  
  #2. remove subjects whose cov is less than z^*_(1)
  df=df[order(df$z1),]
  delta.z1=aggregate(df$status,by=list(z1=df$z1),sum)
  if(shape1=='increasing'){
    if(delta.z1[1,2]==0){ #first obs is censored
      z1.star=delta.z1[which(delta.z1[,2]>0)[1],1]
      df=df[df$z1>=z1.star,]
    }
  }else if(shape1=='decreasing'){
    if(delta.z1[nrow(delta.z1),2]==0){ #last obs is censored
      z1.star=delta.z1[max(which(delta.z1[,2]>0)),1]
      df=df[df$z1<=z1.star,]
    }
  }
  
  df=df[order(df$z2),]
  delta.z2=aggregate(df$status,by=list(z2=df$z2),sum)
  if(shape2=='increasing'){
    if(delta.z2[1,2]==0){ #first obs is censored
      z2.star=delta.z2[which(delta.z2[,2]>0)[1],1]
      df=df[df$z2>=z2.star,]
    }
  }else if(shape2=='decreasing'){
    if(delta.z2[nrow(delta.z2),2]==0){
      z2.star=delta.z2[max(which(delta.z2[,2]>0)),1]
      df=df[df$z2<=z2.star,]
    }
  }
  
  #3. redefine variables
  df1=df[order(df$z1),];          df2=df[order(df$z2),]
  n1=nrow(df1);                   n2=nrow(df2) #n1 & n2 are the same
  status1=df1$status;             status2=df2$status
  
  t1=df1$time;                    t2=df2$time
  t1.obs=sort(unique(t1));        t2.obs=sort(unique(t2)) #same
  nt1=length(t1.obs);             nt2=length(t2.obs)      #same
  
  z1=df1$z1;                      z2=df2$z2
  z1.obs=unique(z1[status1==1]);  z2.obs=unique(z2[status2==1])  
  m1=length(z1.obs);              m2=length(z2.obs)    

  if(Q==0){
    w1=matrix(rep(0,n1),ncol=1);      w2=matrix(rep(0,n2),ncol=1)
  }else if(Q==1){
    w1=matrix(df1[,-c(1:4)],ncol=1);  w2=matrix(df2[,-c(1:4)],ncol=1)
  }else{
    w1=as.matrix(df1[,-c(1:4)]);      w2=as.matrix(df2[,-c(1:4)])
  }
  
  #anchor
  k1=sum(z1.obs<K1);              k2=sum(z2.obs<K2)
  if(k1==0) k1=1;                 if(k2==0) k2=1
  zk1=z1.obs[k1];                 zk2=z2.obs[k2]   
  
  #counting process
  Y1=dN1=matrix(0,n1,nt1)       #row is the subj, col is the time corresponding z_(i);
  Y2=dN2=matrix(0,n2,nt2)
  for(i in 1:n1){
    rank.t1=which(t1[i]==t1.obs)
    Y1[i,][1:rank.t1]=1
    if(status1[i]==1) dN1[i,][rank.t1]=1
  }
  for(i in 1:n2){
    rank.t2=which(t2[i]==t2.obs)
    Y2[i,][1:rank.t2]=1
    if(status2[i]==1) dN2[i,][rank.t2]=1    
  }
  
  #initial value
  int.val=initial.ft(df,Q,z1.obs,zk1,z2.obs,zk2,m1,m2,shape1,shape2)
  psi1=int.val$psi1
  psi2=int.val$psi2
  beta=int.val$beta
  
  #cycling picm & Newton Raphson algorithm
  iter=0;  dist=1;
  psi2.full=BTF.ft(m=m2,n=n2,z=z2,z.obs=z2.obs,psi.new=psi2,shape=shape2)

  if(Q==0){
    while(dist>=eps){
      iter=iter+1
      if(iter>maxiter) break    
      
      #estimate psi1
      psi1.new=cpicm.ft(psi1,psi2.full,z1,z1.obs,m1,k1,n1,Y1,dN1,nt1, w1,beta ,eps,maxiter, shape1)
      psi1.full=BTF.ft(m=m1,n=n1,z=z1,z.obs=z1.obs,psi.new=psi1.new,shape=shape1)

      #estimate psi2;
      psi2.new=cpicm.ft(psi2,psi1.full,z2,z2.obs,m2,k2,n2,Y2,dN2,nt2, w2,beta ,eps,maxiter, shape2)
      psi2.full=BTF.ft(m=m2,n=n2,z=z2,z.obs=z2.obs,psi.new=psi2.new,shape=shape2)

      #update;
      dist=sqrt(sum(psi1.new-psi1)^2)+sqrt(sum(psi2.new-psi2)^2)
      psi1=psi1.new
      psi2=psi2.new
      
      #Note: speed can be improved, by ignoring wbeta, so BTF can be removed
    }
  }else{
    while(dist>=eps){
      iter=iter+1
      if(iter>maxiter) break    
      
      #estimate psi1
      psi1.new=cpicm.ft(psi1,psi2.full,z1,z1.obs,m1,k1,n1,Y1,dN1,nt1, w1,beta ,eps,maxiter, shape1)
      psi1.full=BTF.ft(m=m1,n=n1,z=z1,z.obs=z1.obs,psi.new=psi1.new,shape=shape1)
      
      #estimate psi2;
      psi2.new=cpicm.ft(psi2,psi1.full,z2,z2.obs,m2,k2,n2,Y2,dN2,nt2, w2,beta ,eps,maxiter, shape2)
      psi2.full=BTF.ft(m=m2,n=n2,z=z2,z.obs=z2.obs,psi.new=psi2.new,shape=shape2)
      
      #estimate beta;    
      beta.new=NR.ft(w1,beta,Q,psi1.full,psi2.full,n1,Y1,dN1, maxiter,eps) #does not matter based on df1 or df2
      
      #update;
      dist=sqrt(sum(psi1.new-psi1)^2)+sqrt(sum(psi2.new-psi2)^2)+sqrt(sum((beta.new-beta)^2))
      psi1=psi1.new
      psi2=psi2.new
      beta=beta.new
    }
  }
      
  #cyclic picm result
  conv="converged"
  if(dist>=eps) conv="not converged"
  
  psi1.full=BTF.ft(m=m1,n=n.full,z=z1.full,z.obs=z1.obs,psi.new=psi1,shape=shape1) #for full data
  psi2.full=BTF.ft(m=m2,n=n.full,z=z2.full,z.obs=z2.obs,psi.new=psi2,shape=shape2)
  
  psi1.full[is.infinite(psi1.full)]=NA #change -Inf or Inf to NA 
  psi2.full[is.infinite(psi2.full)]=NA
  
  iso1=data.frame(z=z1.full,psi=psi1.full,HR=exp(psi1.full))
  iso2=data.frame(z=z2.full,psi=psi2.full,HR=exp(psi2.full))
  
  iso1$pj="no" #potential jump points
  iso2$pj="no"
  
  for(i in 1:nrow(iso1))
    if(iso1$z[i]%in%z1.obs) iso1$pj[i]="yes"
  
  for(i in 1:nrow(iso2))
    if(iso2$z[i]%in%z2.obs) iso2$pj[i]="yes"
  
  est=NA
  if(Q>0){
    est=data.frame(beta=beta,HR=exp(beta.new))
    rownames(est)=paste0("w",1:Q)
  }

  return(list(iso1=iso1,iso2=iso2,est=est,
              conv=conv,shape1=shape1,shape2=shape2))
}
