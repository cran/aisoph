print.aisoph=function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\nEstimated monotone function, hazard ratio and level set for psi1:\n")
  print(x$est1)
  cat("Estimated monotone function, hazard ratio and level set for psi2:\n")
  print(x$est2)
  if(!is.na(x$exp.beta)){
    a1.0=paste("\n exp(beta): ",x$exp.beta,"\n",sep='')
    cat(a1.0)
  }

  a1.1=paste("\nNumber of events/subjects was ",x$nevent,"/",x$n,".",sep='')
  a1.21=paste("\nNumber of distinct covariates associated with observed events for psi1 was ",x$njump1,".",sep='')
  a1.22=paste("\nNumber of distinct covariates associated with observed events for psi2 was ",x$njump2,".",sep='')
  a1.3=paste("\nShape restriction was monotone ",x$shape,".",sep='')
  cat(a1.1)
  cat(a1.21)
  cat(a1.22)
  cat(a1.3)
  
  if(x$conv=="converged"){
    a2=paste("\n\nAlgorithm was ",x$conv,".\n",sep='')
  }else if(x$conv=="not converged"){
    a2=paste("\n\nAlgorithm was ",x$conv,". Results were questinable.\n",sep='')
  }
  cat(a2)
}


