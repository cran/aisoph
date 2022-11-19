aisoph=function(time, status, z1, z2, trt=NULL, shape1="increasing", shape2="increasing", K1=0, K2=0, maxdec=2, maxiter=10^5, eps=10^-3){
  

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
  
  #4. isoph
  est=aisoph.ti(TIME=time, STATUS=status, Z1=z1, Z2=z2, W=trt, shape1=shape1, shape2=shape2, K1=K1, K2=K2, maxdec=maxdec, maxiter=maxiter, eps=eps)
    
  est$call=match.call()
  class(est)="aisoph"
  
  est
}