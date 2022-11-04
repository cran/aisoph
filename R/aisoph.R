aisoph=function(time, status, z1, z2, trt=NULL, shape="increasing", K1=0, K2=0, maxdec=2, maxiter=10^5, eps=10^-3){
  
  if ( any(is.na(time)+is.na(status)+is.na(z1)+is.na(z2)) ) stop("Data included NA")
  if ( any(is.infinite(time)+is.infinite(status)+is.infinite(z1)+is.infinite(z2)) )  stop("Data included infinite values")
  
  if(min(time)<0) stop("'time must be greater than zero")
  if ( length(unique(status))>=3 ) stop("status must be either 0 or 1")
  if ( !all(status %in% c(0,1)) )  stop("status must be either 0 or 1")
  
  if(any(grep("inc",tolower(shape)))==TRUE){
    shape="increasing"
  }else if(any(grep("dec",tolower(shape)))==TRUE){
    shape="decreasing"
  }
  if(shape!='increasing' && shape!='decreasing')
    stop("shape must be either increasing or decreasing")  
  
  if(!is.null(trt)){
    if( any(is.na(trt)) ) stop("trt included NA")
    if( any(is.infinite(trt)) ) stop("trt included NA")
    
    uniq.trt=unique(trt)
    if(length(uniq.trt)!=2)  stop("trt must be coded by 0 and 1")
    if(min(uniq.trt)!=0)       stop("trt must be coded by 0 and 1")
    if(max(uniq.trt)!=1)       stop("trt must be coded by 0 and 1")
  }  
  
  
  #4. isoph
  est=aisoph.ti(TIME=time, STATUS=status, Z1=z1, Z2=z2, W=trt, shape=shape, K1=K1, K2=K2, maxdec=maxdec, maxiter=maxiter, eps=eps)
    
  est$call=match.call()
  class(est)="aisoph"
  
  est
}