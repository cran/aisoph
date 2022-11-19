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

  if(x$conv=="converged"){
    a2=paste("\n\nAlgorithm was ",x$conv,".\n",sep='')
  }else if(x$conv=="not converged"){
    a2=paste("\n\nAlgorithm was ",x$conv,". Results were questinable.\n",sep='')
  }
  cat(a2)
}


