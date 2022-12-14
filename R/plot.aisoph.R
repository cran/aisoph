plot.aisoph=function(x,...){
  iso1=x$iso1;  iso2=x$iso2
  hr1=iso1$HR;  z1=iso1$z
  hr2=iso2$HR;  z2=iso2$z
  
  iso1.pj=iso1[which(iso1$pj=="yes"),]
  iso2.pj=iso2[which(iso2$pj=="yes"),]
  hr1.pj=iso1.pj$HR;  z1.pj=iso1.pj$z
  hr2.pj=iso2.pj$HR;  z2.pj=iso2.pj$z
  
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
  points(hr1.pj~z1.pj,pch=19)
  
  plot(hr2~z2,type=type2, xlab=xlab2, ylab=ylab2)
  points(hr2.pj~z2.pj,pch=19)
}
