#newton-raphson algo with p=1
NR.ft=function(w1,beta,psi1.full,psi2.full,n,Y1,dN1, maxiter,eps){
  iter=0
  dist=1
  while(dist>=eps){   
    iter=iter+1
    if(iter>maxiter) break  
    
    w1b=w1*beta  
    Y1est=matrix(,n,n)
    S1=rep(0,n)    
    for(j in 1:n){
      Y1est[,j]=Y1[,j]*exp(psi1.full+psi2.full+w1b)
      S1[j]=sum(Y1est[,j]*w1)
    }
    S2=S1^2
    S0=colSums(Y1est)
    
    #0/0=0
    idx=which(S0>0)
    S2=S2[idx]
    S1=S1[idx]
    S0=S0[idx]
    dN1=dN1[,idx]
    
    #update beta;
    U=0;    I=0
    for(i in 1:n){
      U=U+sum( ((w1[i]-S1/S0)*dN1[i,]) ) #sum over time
      I=I+sum( ((-S2/S0+1)*dN1[i,]) )
    }
    beta.new = beta - U/I
    
    #distance
    dist=(beta.new-beta)^2
    beta=beta.new    
  }
  
  return(beta.new)
}