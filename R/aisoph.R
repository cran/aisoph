aisoph=function(time, status, z1, z2, x=NULL, shape1="increasing", shape2="increasing", K1=NULL, K2=NULL, maxiter=10^5, eps=10^-3){
  
  if(any(grep("inc",tolower(shape1)))==TRUE){
    shape1="increasing"
  }else if(any(grep("dec",tolower(shape1)))==TRUE){
    shape1="decreasing"
  }
  
  if(any(grep("inc",tolower(shape2)))==TRUE){
    shape2="increasing"
  }else if(any(grep("dec",tolower(shape2)))==TRUE){
    shape2="decreasing"
  }
  
  if(is.null(x)){
    Q=0
    df.full=data.frame(time=time,status=status,z1=z1,z2=z2)
  }else{
    df.full=data.frame(time=time,status=status,z1=z1,z2=z2,x=x)
    if(is.vector(x)){
      Q=1
    }else{
      Q=ncol(x)
    }
  }
  df.full=na.omit(df.full) #complete case analysis
  
  if(is.null(K1))
    K1=median(df.full$z1)
  
  if(is.null(K2))
    K2=median(df.full$z2)
  
  #4. aisoph
  res=aisoph.ti(df.full=df.full, Q=Q, shape1=shape1, shape2=shape2, K1=K1, K2=K2, maxiter=maxiter, eps=eps)
  class(res)="aisoph"
  return(res)
}