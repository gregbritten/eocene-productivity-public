###################################################################
## SUM OF SQUARES FUNCTIONS #######################################
###################################################################
ss_alpha <- function(par){
  y_pred   <- f_troph(Ni=input$Ni,is=input$is,P1=P0*input$P1,alpha=par,dfrac=dfrac0,sigma=sigma0)
  y_pred_i <- y_pred[b*is+loga >size_min & b*is+loga <= size_max]
  y_pred_i <- approx(x=i2,    y=y_pred_i,xout=xout)$y
  e <- sum((y_pred_i-y_obs_i)^2,na.rm=TRUE)
  return(e)    
}

ss_P <- function(par){
  y_pred   <- f_troph(Ni=input$Ni,is=input$is,P1=par*input$P1,alpha=alpha0,dfrac=dfrac0,sigma=sigma0)
  y_pred_i <- y_pred[b*is+loga >size_min & b*is+loga <= size_max]
  y_pred_i <- approx(x=i2,    y=y_pred_i,xout=xout)$y
  e <- sum((y_pred_i-y_obs_i)^2,na.rm=TRUE)
  return(e)    
}

ss_P_size <- function(par){
  mu       <- -75 + 0.003*par
  P_size   <- c(dnorm(c(1:100),mean=mu,sd=5),rep(0,Ni-100))
  y_pred   <- f_troph(Ni=input$Ni,is=input$is,P1=par*P_size,alpha=alpha0,dfrac=dfrac0,sigma=sigma0)
  y_pred_i <- y_pred[b*is+loga >size_min & b*is+loga <= size_max]
  y_pred_i <- approx(x=i2,    y=y_pred_i,xout=xout)$y
  e <- sum((y_pred_i-y_obs_i)^2,na.rm=TRUE)
  return(e)    
}

ss_dfrac <- function(par){
  y_pred   <- f_troph(Ni=input$Ni,is=input$is,P1=P0*input$P1,alpha=alpha0,dfrac=par,sigma=sigma0)
  y_pred_i <- y_pred[b*is+loga >size_min & b*is+loga <= size_max]
  y_pred_i <- approx(x=i2,    y=y_pred_i,xout=xout)$y
  e <- sum((y_pred_i-y_obs_i)^2,na.rm=TRUE)
  return(e)    
}

ss_sigma <- function(par){
  y_pred   <- f_troph(Ni=input$Ni,is=input$is,P1=P0*input$P1,alpha=alpha0,dfrac=dfrac0,sigma=par)
  y_pred_i <- y_pred[b*is+loga >size_min & b*is+loga <= size_max]
  y_pred_i <- approx(x=i2,    y=y_pred_i,xout=xout)$y
  e <- sum((y_pred_i-y_obs_i)^2,na.rm=TRUE)
  return(e)    
}

ss_P_alpha <- function(par){
  y_pred   <- f_troph(Ni=input$Ni,is=input$is,P1=par[1]*input$P1,alpha=par[2],dfrac=dfrac0,sigma=sigma0)
  y_pred_i <- y_pred[b*is+loga >size_min & b*is+loga <= size_max]
  y_pred_i <- approx(x=i2,    y=y_pred_i,xout=xout)$y
  e <- sum((y_pred_i-y_obs_i)^2,na.rm=TRUE)
  return(e)    
}

ss_alpha_factor <- function(par){
  y_pred   <- f_troph_alpha(Ni=input$Ni,is=input$is,P1=P0*input$P1,alpha=par[1],dfrac=dfrac0,sigma=sigma0,alpha_factor=par[2])
  y_pred_i <- y_pred[b*is+loga >size_min & b*is+loga <= size_max]
  y_pred_i <- approx(x=i2,    y=y_pred_i,xout=xout)$y
  e <- sum((y_pred_i-y_obs_i)^2,na.rm=TRUE)
  return(e)    
}


ss_all <- function(par){
  P1    <- c(dnorm(c(1:100),mean=par[4],sd=par[5]),rep(0,Ni-100)) 
  y_pred   <- f_troph(Ni=input$Ni,is=input$is,P1=P0*P1,alpha=par[1],dfrac=par[2],sigma=par[3])
  y_pred_i <- y_pred[b*is+loga >size_min & b*is+loga <= size_max]
  y_pred_i <- approx(x=i2,    y=y_pred_i,xout=xout)$y
  e <- sum((y_pred_i-y_obs_i)^2,na.rm=TRUE)
  return(e)    
}

