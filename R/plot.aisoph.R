plot.aisoph=function(x, lty=1, lcol=1, lwd=1, pch=19, pcol=1, pcex=0.7,
                   ylab1=NULL, xlab1=NULL, ylim1=NULL, xlim1=NULL,
                  ylab2=NULL, xlab2=NULL, ylim2=NULL, xlim2=NULL,
                   ...){
  y1.obs=x$psi1;        y2.obs=x$psi2
  z1.obs=x$z1;          z2.obs=x$z2
  z1.range=x$z1.range;  z2.range=x$z2.range
  hr1.obs=exp(y1.obs);  hr2.obs=exp(y2.obs)
  
  n1=length(y1.obs);    n2=length(y2.obs)
  
  y1=c(y1.obs[1],y1.obs,y1.obs[n1])
  hr1=c(hr1.obs[1],hr1.obs,hr1.obs[n1])
  z1=c(z1.range[1],z1.obs,z1.range[2])
  
  y2=c(y2.obs[1],y2.obs,y2.obs[n2])
  hr2=c(hr2.obs[1],hr2.obs,hr2.obs[n2])
  z2=c(z2.range[1],z2.obs,z2.range[2])
  
  if(is.null(xlab1)) xlab1="cov1"
  if(is.null(xlab2)) xlab2="cov2"
  
  if(is.null(xlim1)) xlim1=z1.range
  if(is.null(xlim2)) xlim2=z2.range
  
  opar=par(mfrow=c(1,2))
  on.exit(par(opar))

  if(is.null(ylab1)) ylab1=expression(exp(hat(psi)))
  if(is.null(ylab2)) ylab2=expression(exp(hat(psi)))
  if(is.null(ylim1)) ylim1=range(exp(y1))
  if(is.null(ylim2)) ylim2=range(exp(y2))
    
  plot(hr1~z1,type='s', xlab=xlab1, ylab=ylab1, xlim=xlim1, ylim=ylim1, col=lcol, lty=lty, lwd=lwd)
  plot(hr2~z2,type='s', xlab=xlab2, ylab=ylab2, xlim=xlim2, ylim=ylim2, col=lcol, lty=lty, lwd=lwd)
}
