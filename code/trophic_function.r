f_troph <- function(Ni,is,P1,alpha,dfrac,sigma){
  P2=P3=P4=P5=P6 <- numeric(Ni) 
  
  for(i in 1:Ni){
    size <- exp(is[i])                                          
    g    <- dnorm(exp(is),mean=dfrac*size,sd=sigma*dfrac*size)  
    g    <- if(sum(g)>0){g <- g/sum(g)}else{g <- rep(0,Ni)}    
    P2[i] <- alpha*sum(g*P1)                                   
  }
  for(i in 1:Ni){
    size <- exp(is[i])
    g    <- dnorm(exp(is),mean=dfrac*size,sd=sigma*dfrac*size)
    g    <- if(sum(g)>0){g <- g/sum(g)}else{g <- rep(0,Ni)}
    P3[i] <- alpha*sum(g*P2)
  }
  for(i in 1:Ni){
    size <- exp(is[i])
    g    <- dnorm(exp(is),mean=dfrac*size,sd=sigma*dfrac*size)
    g    <- if(sum(g)>0){g <- g/sum(g)}else{g <- rep(0,Ni)}
    P4[i] <- alpha*sum(g*P3)
  }
  for(i in 1:Ni){
    size <- exp(is[i])
    g    <- dnorm(exp(is),mean=dfrac*size,sd=sigma*dfrac*size)
    g    <- if(sum(g)>0){g <- g/sum(g)}else{g <- rep(0,Ni)}
    P5[i] <- alpha*sum(g*P4)
  }
  for(i in 1:Ni){
    size <- exp(is[i])
    g    <- dnorm(exp(is),mean=dfrac*size,sd=sigma*dfrac*size)
    g    <- if(sum(g)>0){g <- g/sum(g)}else{g <- rep(0,Ni)}
    P6[i] <- alpha*sum(g*P5)
  }
  return(P4+P5+P6)
}



