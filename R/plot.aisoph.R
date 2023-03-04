plot.aisoph=function(x,...){
  iso1=x$iso1;  iso2=x$iso2
  hr1=iso1$HR;  z1=iso1$z
  hr2=iso2$HR;  z2=iso2$z
  
  iso1.cens=iso1[which(iso1$cens=="yes"),]
  iso2.cens=iso2[which(iso2$cens=="yes"),]
  hr1.cens=iso1.cens$HR;  z1.cens=iso1.cens$z
  hr2.cens=iso2.cens$HR;  z2.cens=iso2.cens$z
  
  xlab1="Cov1"
  xlab2="Cov2"
  opar=par(mfrow=c(1,2))
  on.exit(par(opar))
  
  ylab1="Estimated hazard ratio"
  ylab2="Estimated hazard ratio"
  
  type1=type2='s'
  if(x$shape1=="decreasing")
    type1='S'
  if(x$shape2=="decreasing")  
    type2='S'
  
  plot(hr1~z1,type=type1, xlab=xlab1, ylab=ylab1)
  points(hr1.cens~z1.cens,pch=3)
  
  plot(hr2~z2,type=type2, xlab=xlab2, ylab=ylab2)
  points(hr2.cens~z2.cens,pch=3)
}
