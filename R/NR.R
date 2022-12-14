#newton-raphson algo with p=1
NR.ft=function(w1,beta,Q,psi1.full,psi2.full,n1,Y1,dN1, maxiter,eps){
  iter=0
  dist=1
  
  exp.psi=exp(psi1.full+psi2.full)
  while(dist>=eps){
    iter=iter+1
    if(iter>maxiter) break
    
    w1b=w1%*%beta
    Yest=matrix(0,n1,1); S1=matrix(0,Q,1)
    U=matrix(0,Q,1);    H=matrix(0,Q,Q)
    
    exp.psi.w1b=exp.psi*exp(w1b) #Y*estimated exp(psi.hat+wb)

    for(i in 1:n1){
      for(j in which(dN1[i,]>=1)){
        Yest=Y1[,j]*exp.psi.w1b
        
        S0=sum(Yest)
        S1=matrix(t(Yest) %*% w1)
        E1=S1/S0
        U=U+(w1[i,]-E1)*dN1[i,j]
              
        S2=matrix(0,Q,Q) 
        #for(m in 1:n1) #slow
        #  S2=S2+Yest[m]*(w1[m,]%*%t(w1[m,]))
        
        S2=matrix(0,Q,Q) 
        if(Q==1){
          S2=sum(Yest*w1^2) #Q=1
        }else{
          for(s in 1:Q)
              S2[s,s]=sum(w1[,s]*Yest)
          for(s in 1:(Q-1))
            for(t in (s+1):Q)
              S2[s,t]=S2[t,s]=sum(w1[,s]*w1[,t]*Yest)
        }
        H=H+(-S2/S0+E1%*%t(E1))*dN1[i,j]
      }
    }
    
    beta.new = beta - solve(H)%*%U
    
    #distance
    dist=sqrt(sum((beta.new-beta)^2))
    beta=beta.new
  }
  return(beta.new)
}
