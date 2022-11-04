cpicm.ft=function(psi1,psi2.full, z1,m1,k1,n,int1,Y1,dN1, w1,beta, eps,maxiter, shape){

  Yest1=matrix(NA,n,n)
  for(j in 1:n)
    Yest1[,j]=Y1[,j]*exp(psi2.full+w1*beta)
  rpa.Y1=Y2ft(m1, n, int1, z1, Yest1, dN1)

  Y1.new=rpa.Y1$Y2
  dN1.new=rpa.Y1$dN2
  picm=picm.ft(Y1.new, dN1.new, psi1, k1, m1, n ,eps,maxiter, shape)
  
  psi.new=picm$psi.new
  conv=picm$conv     
  
  return(list(psi.new=psi.new, conv=conv))
}
  
picm.ft=function(Y2, dN2, psi, k, m, n ,eps,maxiter, shape){
  dNsum=colSums(dN2)
  Delta=rowSums(dN2)     
  
  #picm
  iter=0
  d.e=1
  while(d.e>=eps){  
    iter=iter+1
    if(iter>maxiter) break    
    
    #picm
    den=colSums(Y2*exp(psi))
    index.zero=which(den>0) #0/0=0
    
    weight=c()
    for(s in 1:m)
      weight[s]=sum( (Y2[s,]*dNsum/den)[index.zero] )
    
    if(sum(is.na(weight))+sum(is.infinite(weight))>=1)
      break;
    
    if(shape=='increasing'){
      exp.psi.new=pava(Delta/weight, weight)
    }else if(shape=='decreasing'){
      exp.psi.new=-pava(-Delta/weight, weight)
    }
    psi.new=log(exp.psi.new)    
    
    #distance
    d.e=sum(abs(exp(psi.new)-exp(psi)))
    psi=psi.new
  }
  
  #impose the anchor
  psi.new=psi.new-psi.new[k] #psi is the same as psi.new;
  
  conv=0
  if(d.e<eps) conv=1
  
  return(list(psi.new=psi.new, conv=conv))
}