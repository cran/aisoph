aisoph.ti=function(TIME, STATUS, Z1, Z2, W, shape, K1, K2, maxdec, maxiter, eps){
  #sorted by z
  n=length(STATUS)
  order.z1=order(Z1);             order.z2=order(Z2)  
  t1=TIME[order.z1];              t2=TIME[order.z2]  
  status1=STATUS[order.z1];       status2=STATUS[order.z2]

  z1=sort(Z1);                    z2=sort(Z2)  
  z1.obs=unique(z1[status1==1]);  z2.obs=unique(z2[status2==1])  
  m1=length(z1.obs);              m2=length(z2.obs)    
  
  #t1.obs=sort(unique(t1));       t2.obs=sort(unique(t2))
  #nt1=length(t1.obs);            nt2=length(t2.obs)  

  #anchor
  k1=sum(z1.obs<K1);              k2=sum(z2.obs<K2)  
  if(k1==0) k1=1;                 if(k2==0) k2=1
  zk1=z1.obs[k1];                 zk2=z2.obs[k2]    
  
  #counting process
  rank.t1=rank(t1,ties.method='min')
  rank.t2=rank(t2,ties.method='min')
  Y1=dN1=matrix(,n,n)       #row is the subj, col is the time corresponding z_(i);
  Y2=dN2=matrix(,n,n)
  for(i in 1:n){
    Y1[i,]=rep(0,n)
    Y1[i,][1:rank.t1[i]]=1
    dN1[i,]=rep(0,n)
    if(status1[i]==1) dN1[i,][rank.t1[i]]=1
    
    Y2[i,]=rep(0,n)
    Y2[i,][1:rank.t2[i]]=1
    dN2[i,]=rep(0,n)
    if(status2[i]==1) dN2[i,][rank.t2[i]]=1    
  }  
  
  #interval (RPA)
  int1=list()
  int1[[1]]=c(-Inf,z1.obs[2])
  int1[[m1]]=c(z1.obs[m1],Inf)
  for(i in 2:(m1-1))
    int1[[i]]=c(z1.obs[i],z1.obs[i+1])
  
  int2=list()
  int2[[1]]=c(-Inf,z2.obs[2])
  int2[[m2]]=c(z2.obs[m2],Inf)
  for(i in 2:(m2-1))
    int2[[i]]=c(z2.obs[i],z2.obs[i+1])

  #initial value
  if(is.null(W)){ #no trt group
    w1=rep(0,n); w2=rep(0,n)
    beta=0
    
    try(beta.hat<-coxph(Surv(TIME,STATUS)~Z1+Z2)$coefficient, silent=T)
    if(!is.numeric(beta.hat)) beta.hat=c(0.01,0.01)
    
    if(shape=='increasing'){
      psi1= abs(beta.hat[1])*(z1.obs-zk1);    psi2= abs(beta.hat[2])*(z2.obs-zk2)
    }else if(shape=='decreasing'){
      psi1=-abs(beta.hat[1])*(z1.obs-zk1);    psi2=-abs(beta.hat[2])*(z2.obs-zk2)
    }
  }else{
    w1=W[order.z1];                 w2=W[order.z2]    
    
    try(beta.hat<-coxph(Surv(TIME,STATUS)~Z1+Z2+W)$coefficient, silent=T)
    if(!is.numeric(beta.hat)){; beta.hat=c(0.01,0.01); beta=0;
    }else{; beta=beta.hat[3]
    }
    if(shape=='increasing'){
      psi1= abs(beta.hat[1])*(z1.obs-zk1);    psi2= abs(beta.hat[2])*(z2.obs-zk2)
    }else if(shape=='decreasing'){
      psi1=-abs(beta.hat[1])*(z1.obs-zk1);    psi2=-abs(beta.hat[2])*(z2.obs-zk2)
    }
  }

  #cycling picm & Newton Raphson algorithm
  psi2.full=BTFft(m2, n, int2, z2, psi2)
  iter=0;  dist=1;
  
  while(dist>=eps){
    iter=iter+1
    if(iter>maxiter) break    
    
    #estimate psi1
    picm1=cpicm.ft(psi1,psi2.full,z1,m1,k1,n,int1,Y1,dN1, w1,beta ,eps,maxiter, shape)
    if(picm1$conv==0) stop
    psi1.new=picm1$psi.new
    psi1.full=BTFft(m1, n, int1, z1, psi1.new)

    #estimate psi2;
    picm2=cpicm.ft(psi2,psi1.full,z2,m2,k2,n,int2,Y2,dN2, w2,beta ,eps,maxiter, shape)
    if(picm2$conv==0) stop
    psi2.new=picm2$psi.new
    psi2.full=BTFft(m2, n, int2, z2, psi2.new)

    #estimate beta (Y1&w1 or Y2&w2 should be the same);    
    if(is.null(W)){;  beta.new=0
    }else{;           beta.new=NR.ft(w1,beta,psi1.full,psi2.full,n,Y1,dN1, maxiter,eps)
    }
      
    #update;
    dist=sqrt(sum(psi1.new-psi1)^2)+
         sqrt(sum(psi2.new-psi2)^2)+
         (beta.new-beta)^2
    psi1=psi1.new;    psi2=psi2.new;         beta=beta.new;
  }
    
  #cyclic picm result
  conv="converged"
  if(dist>=eps) conv="not converged"
  
  #back to full rank (later)
  psi1.obs=round(psi1.new, maxdec)
  psi2.obs=round(psi2.new, maxdec)
  
  #level sets for psi1
  z=z1;  z.obs=z1.obs;  psi.obs=psi1.obs
  psi.uniq=unique(psi.obs)
  n.lv=length(psi.uniq)
  
  lv.sets=c()
  zmin=formatC( round(min(z),maxdec), format='f', digits=maxdec)
  zmax=formatC( round(max(z),maxdec), format='f', digits=maxdec)
  
  if(n.lv==1){ #only one level sets
    lv.sets[1]=paste('[',zmin,',',zmax,']', sep='')
  }else{
    lv=c()
    for(i in 1:n.lv)
      lv[[i]]=formatC( round(z.obs[which(psi.obs==psi.uniq[i])],maxdec)[1], format='f', digits=maxdec)
    for(i in 1:n.lv){
      if(i==1){
        lv.sets[1]=paste('[',zmin,', ',lv[2],')', sep='')
      }else if(i==n.lv){
        lv.sets[i]=paste('[',lv[i][1],', ',zmax,']', sep='')
      }else{
        lv.sets[i]=paste('[',lv[i],', ',lv[i+1],')', sep='')
      }
    }
  }
  psi.hat=formatC( unique(psi.obs), format='f', digits=maxdec)
  HR.hat=formatC( unique(exp(psi.obs)), format='f', digits=maxdec)
  
  est1=data.frame(psi.hat=psi.hat, HR.hat=HR.hat, lv.set=lv.sets)
  names(est1)=c("psi.hat","exp(psi.hat)","level set of psi.hat")
  
  #level sets for psi2
  z=z2;  z.obs=z2.obs;   psi.obs=psi2.obs
  psi.uniq=unique(psi.obs)
  n.lv=length(psi.uniq)
  
  lv.sets=c()
  zmin=formatC( round(min(z),maxdec), format='f', digits=maxdec)
  zmax=formatC( round(max(z),maxdec), format='f', digits=maxdec)
  
  if(n.lv==1){ #only one level sets
    lv.sets[1]=paste('[',zmin,',',zmax,']', sep='')
  }else{
    lv=c()
    for(i in 1:n.lv)
      lv[[i]]=formatC( round(z.obs[which(psi.obs==psi.uniq[i])],maxdec)[1], format='f', digits=maxdec)
    for(i in 1:n.lv){
      if(i==1){
        lv.sets[1]=paste('[',zmin,', ',lv[2],')', sep='')
      }else if(i==n.lv){
        lv.sets[i]=paste('[',lv[i][1],', ',zmax,']', sep='')
      }else{
        lv.sets[i]=paste('[',lv[i],', ',lv[i+1],')', sep='')
      }
    }
  }
  psi.hat=formatC( unique(psi.obs), format='f', digits=maxdec)
  HR.hat=formatC( unique(exp(psi.obs)), format='f', digits=maxdec)
  
  est2=data.frame(psi.hat=psi.hat, HR.hat=HR.hat, lv.set=lv.sets)
  names(est2)=c("psi.hat","exp(psi.hat)","level set of psi.hat")  
  
  #range  
  z1.range=range(z1)
  z2.range=range(z2)
  
  if(is.null(W)){; exp.beta=NA
  }else{; exp.beta=round(exp(beta.new), maxdec)
          exp.beta=formatC( exp.beta, format='f', digits=maxdec)
  }
  
  return(list(est1=est1, est2=est2,
              psi1=psi1.obs, psi2=psi2.obs, exp.beta=exp.beta, z1=z1.obs, z2=z2.obs,
              z1.range=z1.range, z2.range=z2.range, conv=conv,
              K1=K1, K2=K2, shape=shape,
              n=n, nevent=sum(STATUS), njump1=m1, njump2=m2))
}
