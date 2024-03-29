#newton-raphson algo with p=1
NR.ft=function(x1,beta,Q,psi1.full,psi2.full1,n,Y1,dN1, maxiter,eps){
  iter=0
  dist=1
  
  exp.psi=exp(psi1.full+psi2.full1)
  while(dist>=eps){
    iter=iter+1
    if(iter>maxiter) break
    
    x1b=x1%*%beta
    Yest=matrix(0,n,1); S1=matrix(0,Q,1)
    U=matrix(0,Q,1);    H=matrix(0,Q,Q)
    
    exp.psi.x1b=exp.psi*exp(x1b) #Y*estimated exp(psi.hat+wb)

    for(i in 1:n){
      for(j in which(dN1[i,]>=1)){
        Yest=Y1[,j]*exp.psi.x1b
        Yest=Yest/sum(Yest) #normalization for stablizting computations
        
        S0=sum(Yest)
        S1=matrix(t(Yest) %*% x1)
        E1=S1/S0
        U=U+(x1[i,]-E1)*dN1[i,j]
              
        S2=matrix(0,Q,Q) 
        #for(m in 1:n) #slow
        #  S2=S2+Yest[m]*(x1[m,]%*%t(x1[m,]))
        
        S2=matrix(0,Q,Q) 
        if(Q==1){
          S2=sum(Yest*x1^2) #Q=1
        }else{
          for(s in 1:Q)
              S2[s,s]=sum(x1[,s]*Yest)
          for(s in 1:(Q-1))
            for(t in (s+1):Q)
              S2[s,t]=S2[t,s]=sum(x1[,s]*x1[,t]*Yest)
        }
        H=H+(-S2/S0+E1%*%t(E1))*dN1[i,j]
      }
    }
    beta.new = beta - solve(H)%*%U
    
    if(is.infinite(beta.new) | is.na(beta.new)) #exp.psi.x1b could be infinite
      return(rep(NA,Q))
    
    #distance
    dist=sqrt(sum((beta.new-beta)^2))
    if(is.infinite(dist) | is.na(dist))
      return(rep(NA,Q))
    
    beta=beta.new
  }
  return(beta.new)
}
