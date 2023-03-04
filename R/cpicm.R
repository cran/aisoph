cpicm.ft=function(psi1,psi2.full1,z1,z1.obs,m1,k1,n,Y1,dN1,nt, x1,beta, eps,maxiter, shape1){
  #Yest1
  x1b=x1%*%beta
  exp.psi2.x1b=exp(psi2.full1+x1b) 
  Yest1=matrix(NA,n,nt)
  for(j in 1:nt)
    Yest1[,j]=Y1[,j]*exp.psi2.x1b
  
  Y1.new=matrix(0,m1,nt)
  dN1.new=matrix(0,m1,nt)
  
  if(shape1=="increasing"){
    intv1=c(z1.obs,Inf) #right continuous
    for(h in 1:m1){
      idx=which(intv1[h]<=z1 & z1<intv1[h+1]) #right continuous
      if(length(idx)==1){
        Y1.new[h,]=Yest1[idx,]
        dN1.new[h,]=dN1[idx,]
      }else{
        Y1.new[h,]=colSums(Yest1[idx,])
        dN1.new[h,]=colSums(dN1[idx,])
      }
    }
  }else if(shape1=="decreasing"){
    intv1=c(-Inf,z1.obs) #left continuous
    for(h in 1:m1){
      idx=which(intv1[h]<z1 & z1<=intv1[h+1]) #left continuous
      if(length(idx)==1){
        Y1.new[h,]=Yest1[idx,]
        dN1.new[h,]=dN1[idx,]
      }else{
        Y1.new[h,]=colSums(Yest1[idx,])
        dN1.new[h,]=colSums(dN1[idx,])
      }
    }
  }
  
  #picm
  dNsum=colSums(dN1.new)
  Delta=rowSums(dN1.new)     
  
  iter=0
  d.e=1
  
  while(d.e>=eps){  
    iter=iter+1
    if(iter>maxiter) break    
    
    #picm
    den=colSums(Y1.new*exp(psi1))
    index.zero=which(den>0) #0/0=0
    
    weight=rep(NA,m1)
    for(s in 1:m1)
      weight[s]=sum( (Y1.new[s,]*dNsum/den)[index.zero] )
    
    if((sum(is.na(weight))+sum(is.infinite(weight))+sum(is.infinite(1/weight)))>=1)
      return(rep(NA,m1))
    
    
    if(shape1=='increasing'){
      exp.psi1.new=pava((Delta/weight), weight)
    }else if(shape1=='decreasing'){
      exp.psi1.new=-pava(-(Delta/weight), weight)
    }
    psi1.new=log(exp.psi1.new)    
    
    #distance
    d.e=sum(abs(exp(psi1.new)-exp(psi1)))
    if(is.infinite(d.e) | is.na(d.e))
      return(rep(NA,m1))
    psi1=psi1.new
  }
  
  #impose the anchor
  psi1=psi1-psi1[k1] #psi is the same as psi.new;
  
  if(d.e>eps)
    return(rep(NA,m1))

  return(psi1)
}
