rm(list=ls())
setwd('d:/dropbox/working/eocene/eocene-productivity-public/')
source('code/process_data.r')
source('code/trophic_function.r')
source('code/parameters.r')
source('code/ss_functions.r')

ages   <- unique(size$age)
cramer <- read.csv('data/cramer_2009_avg.csv')

t_min <- 62
t_max <- 46

c_age <- cramer$Age_ma 
c_o18 <- cramer$Pacific_d18O_trend
c_ind <- which(c_age<t_min & c_age>t_max)

ft <- function(do18){return(16.5 - 5*(do18+1.2) + 0.14*(do18+1.2)^2)}

par(mfrow=c(1,1),cex.axis=0.9,cex.lab=1.0,mar=c(4,4,4,6))
  plot(ages[ages<=t_min&ages>=t_max],accum$ich_accum[ages<=t_min&ages>=t_max],xlim=c(62,46),bty='n',pch=19,xlab='',ylab='',ylim=c(0,390),xaxt='n')
    axis(side=1,at=seq(62,46,-2))
    mtext(side=2,line=2.5,expression('Ichthyolith Accumulation Rate (teeth/cm'^2*'/Myr)'))
    mtext(side=1,line=2.5,'Myr')
par(new=TRUE)
  plot(c_age[c_ind],ft(c_o18[c_ind]) ,xlim=c(62,46),xaxt='n',yaxt='n',xlab='',ylab='',bty='n',
	cex=1.1,type='l',col='red',ylim=c(9,16))
      axis(side=4,col='red')
  legend('topleft',pch=c(19,NA),lty=c(NA,1),col=c('black','red'),legend=c('Ichthyolith Accumulation',expression('Temperature ('*degree*'C)')),bty='n',cex=1.2)
  mtext(side=4,expression('Deep Water Temperature ('*degree*'C)'),line=2.5)

################################################
## ACCUMS ######################################
################################################
temp           <- read.csv('data/temp-cramer.csv',header=TRUE)
colnames(temp) <- c('age','temp')

int <- 62
fin <- 46
del <- 0.5   
accumss=accumssd=temps2=temps2sd <- c()

i <- 1
while(int>fin){
	tmp <- size[size$age < int & size$age >= int - del,]
	accumss <- c(accumss,mean(tmp$accum))
	accumssd <- c(accumssd,sd(tmp$accum))
	tmp_temp <- temp[temp$age <int & temp$age >= int -del,]
	temps2 <- c(temps2,mean(tmp_temp$temp))
	temps2sd <- c(temps2sd,sd(tmp_temp$temp))
	int <- int - del
	i <- i + 1
}
accumssd[accumssd==0] <- sort(unique(accumssd))[2]

library(mgcv)
mod  <- gam(accumss ~ temps2 + I(temps2^2))
xin  <- seq(9,14.5,0.01)
pred <- predict.gam(mod,newdata=list(temps2=xin),se.fit=TRUE)
    
	par(mfrow=c(1,1),cex.axis=0.9,cex.lab=1.0,mar=c(4,4,4,4))
	plot(temps2,accumss,bty='n',pch=19,ylim=c(0,350),xlim=c(9,15),xlab='',ylab='')
	segments(x0=temps2-temps2sd,x1=temps2+temps2sd,
			 y0=accumss,y1=accumss,lty=1,lwd=0.5)
	segments(x0=temps2,x1=temps2,
			 y0=accumss-accumssd,y1=accumss+accumssd,lty=1,lwd=0.5)
	lines(xin,pred$fit)
	lines(xin,pred$fit+2*pred$se.fit,lty=2)
	lines(xin,pred$fit-2*pred$se.fit,lty=2)
	mtext(side=1,expression('Deep Water Temperature ('*degree*'C)'),line=2.5)
    mtext(side=2,line=2.5,expression('Ichthyolith Accumulation Rate (teeth/cm'^2*'/Myr)'))

####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
xout <- seq(size_min,size_max,0.01)
old <- 62
new <- 46
i2 <- (b*is+loga)[b*is+loga >size_min & b*is+loga <= size_max]

ages   <- unique(size$age) 
cramer <- read.csv('data/cramer_2009_avg.csv')
#########################################################################
## FIT MODELS ###########################################################
#########################################################################
int <- 62
fin <- 46
del <- 1

es_alpha=es_P=es_dfrac=es_sigma=es_const=es_P_alpha=es_alpha_factor=es_P_size <- c()

lens=sds=means_size=sds_size <- numeric()

	while(int>fin){
	  print(int)
	  dat     <- size[size$age < int & size$age >= int-del,]
	  accumm   <- mean(dat$accum)
	  den     <- density(log(dat$Feret/1000))
	  y_obs   <- (den$y/sum(den$y))*accumm
	  x_obs   <- den$x
	  y_obs_i <- y_obs[x_obs >size_min & x_obs <= size_max]
	  x_obs_i <- x_obs[x_obs >size_min & x_obs <= size_max]
	  y_obs_i <- approx(x=x_obs_i,y=y_obs_i, xout=xout)$y
	 
	  lens <- c(lens,length(unique(dat$age)))
	  sds  <- c(sds,sd(dat$accum))
	 
	  means_size <- c(means_size,mean(dat$Feret/1000))	
	  sds_size <- c(sds_size,sd(dat$Feret/1000)) 
	 
	  int <- int - del
	  
	  alpha_hat <- optimize(f=ss_alpha, c(0,1))
	  P_hat     <- optimize(f=ss_P,     c(0,1000000))
	  dfrac_hat <- optimize(f=ss_dfrac, c(0,1))
	  sigma_hat <- optimize(f=ss_sigma, c(0,2))
	  P_size_hat<- optimize(f=ss_P_size,c(0,100000))
	  
	  alphas <- c(alphas,alpha_hat$minimum)
	  Ps     <- c(Ps,    P_hat$minimum)
	  dfracs <- c(dfracs,dfrac_hat$minimum)
	  sigmas <- c(sigmas,sigma_hat$minimum)
	  Ps_size<- c(Ps_size,P_size_hat$minimum)
	  
	  es_const <- c(es_const,ss_sigma(par=sigma0))
	  es_alpha <- c(es_alpha,alpha_hat$objective)
	  es_P     <- c(es_P,    P_hat$objective)
	  es_dfrac <- c(es_dfrac,dfrac_hat$objective)
	  es_sigma <- c(es_sigma,sigma_hat$objective)
	  es_P_size<- c(es_P_size,P_size_hat$objective)
	  
	}

t_min <- 62
t_max <- 46

xins <- seq(62,46,length.out=16)
par(mfrow=c(3,2),cex.axis=0.9,cex.lab=1.0,mar=c(2,4,2,2),oma=c(2,2,2,2))
  plot(ages[ages<=t_min&ages>=t_max],accum$ich_accum[ages<=t_min&ages>=t_max],xlim=c(62,46),bty='n',pch=19,xlab='',ylab='',ylim=c(0,390),xaxt='n')
    axis(side=1,at=seq(62,46,-2))
    mtext(side=2,line=2.5,expression('Teeth/cm'^2*'/Myr'),cex=0.7)
par(new=TRUE)
  plot(c_age[c_ind],ft(c_o18[c_ind]) ,xlim=c(62,46),xaxt='n',yaxt='n',xlab='',ylab='',bty='n',
	cex=1.1,type='l',col='red',ylim=c(9,16))
      axis(side=4,col='red')
  legend('topleft',pch=c(19,NA),lty=c(NA,1),col=c('black','red'),legend=c('Ichthyolith Accumulation',expression('Temperature ('*degree*'C)')),bty='n',cex=1.2)

plot(xins,Ps[-1]/P0,type='b',pch=19,col='blue',bty='n',xlim=c(62,46),xaxt='n',ylim=c(0,1.5),ylab='');  axis(side=1,at=seq(62,46,-2))
	mtext(side=2,line=2.5,'Primary Production Scale Factor',cex=0.7)
plot(xins,alphas[-1],pch=19,type='b',col='red',bty='n',xlim=c(62,46),xaxt='n',ylab='');  axis(side=1,at=seq(62,46,-2))
	mtext(side=2,line=2.5,'Trophic Transfer Efficiency',cex=0.7)
plot(xins,sigmas[-1],pch=19,type='b',col='orange',bty='n',xlim=c(62,46),xaxt='n',ylab=''); axis(side=1,at=seq(62,46,-2))
	mtext(side=2,line=2.5,'Prey Size Range',cex=0.7)
plot(xins,dfracs[-1],pch=19,type='b',col='dark green',bty='n',xlim=c(62,46),xaxt='n',ylab=''); axis(side=1,at=seq(62,46,-2))
	mtext(side=2,line=2.5,'Mean Prey Range',cex=0.7)
plot(xins,Ps_size[-1]/P0,type='b',col='purple',bty='n',xlim=c(62,46),xaxt='n',ylab=''); axis(side=1,at=seq(62,46,-2))
	mtext(side=2,line=2.5,'Size-Dependent Prim. Prod. Factor',cex=0.7)
mtext(side=1,outer=TRUE,'Myr',line=1)


####################################################################################
## PLOTS ###########################################################################
####################################################################################
ylims <- c(0,2.5)
xlims <- c(-10,-5)

int <- 62; fin <- 46; del <- 1
i <- 1
j <- 0

par(mfrow=c(4,4),mar=c(0,0,0,0),oma=c(5,9,5,5),cex.axis=0.9,xpd=TRUE)
while(int>fin){
  print(int)
  dat     <- size[size$age < int & size$age >= int-del,]
  accum   <- mean(dat$accum)
  den     <- density(log(dat$Feret/1000),bw=bw)
  y_obs   <- (den$y/sum(den$y))*accum
  x_obs   <- den$x
  y_obs_i <- y_obs[x_obs >size_min & x_obs <= size_max]
  x_obs_i <- x_obs[x_obs >size_min & x_obs <= size_max]
  y_obs_i <- approx(x=x_obs_i,y=y_obs_i, xout=xout)$y
  
  plot(xout,y_obs_i,xlim=xlims,type='l',ylim=ylims,bty='n',lty=2,lwd=1.5,xaxt='n',yaxt='n')
  
  pred_alpha <- f_troph(Ni=input$Ni,is=input$is,P1=P0*input$P1,alpha=alphas[i+1],dfrac=dfrac0,sigma=sigma0)
  pred_alpha <- pred_alpha[b*is+loga >size_min & b*is+loga <= size_max]
  pred_alpha <- approx(x=i2,  y=pred_alpha,xout=xout)$y
  lines(xout,pred_alpha,col='red',lwd=1.5)

  pred_P <- f_troph(Ni=input$Ni,is=input$is,P1=Ps[i+1]*input$P1,alpha=alpha0,dfrac=dfrac0,sigma=sigma0)
  pred_P <- pred_P[b*is+loga >size_min & b*is+loga <= size_max]
  pred_P <- approx(x=i2,  y=pred_P,xout=xout)$y
  lines(xout,pred_P,col='blue',lty=2,lwd=1.5)

  pred_sigma <- f_troph(Ni=input$Ni,is=input$is,P1=P0*input$P1,alpha=alpha0,dfrac=dfrac0,sigma=sigmas[i+1])
  pred_sigma <- pred_sigma[b*is+loga >size_min & b*is+loga <= size_max]
  pred_sigma <- approx(x=i2,  y=pred_sigma,xout=xout)$y
  lines(xout,pred_sigma,col='orange',lwd=1.5)

  pred_dfrac <- f_troph(Ni=input$Ni,is=input$is,P1=P0*input$P1,alpha=alpha0,dfrac=dfracs[i+1],sigma=sigma0)
  pred_dfrac <- pred_dfrac[b*is+loga >size_min & b*is+loga <= size_max]
  pred_dfrac <- approx(x=i2,  y=pred_dfrac,xout=xout)$y
  lines(xout,pred_dfrac,col='dark green',lwd=1.5)

  P_size      <- c(dnorm(c(1:100),mean=-75 + 0.003*Ps_size[i+1],sd=5),rep(0,Ni-100))
  pred_P_size <- f_troph(Ni=input$Ni,is=input$is,P1=Ps_size[i+1]*P_size,alpha=alpha0,dfrac=dfrac0,sigma=sigma0)
  pred_P_size <- pred_P_size[b*is+loga >size_min & b*is+loga <= size_max]
  pred_P_size <- approx(x=i2,  y=pred_P_size,xout=xout)$y
  lines(xout,pred_P_size,col='purple',lwd=1.5)

  
  lines(xout,y_obs_i,lty=2,lwd=1.5)
  
   if(i%in%c(1,5,9,13)){axis(side=2)}else{axis(side=2,labels=FALSE,tcl=-0.25)}
   if(i%in%c(13:16)){axis(side=1)}else{axis(side=1,labels=FALSE,tcl=-0.25)}
   if(i%in%c(4,8,12,16)){axis(side=4,labels=FALSE,tcl=0.25)}
   
   legend(-10.5,2.5,legend=c(paste(round(es_P[i],digits=1)),
							 paste(round(es_alpha[i],digits=1)),
							 paste(round(es_sigma[i],digits=1)),
							 paste(round(es_dfrac[i],digits=1)),
							 paste(round(es_P_size[i],digits=1))),text.col=c('blue','red','orange','dark green','purple'),bty='n')
 

 
   int <- int - del
   i <- i + 1
   box()
   
   mtext(paste(int+del,' - ',int+del-del, 'Myr'),adj=0.8,line=-2,cex=0.9)
  }
  mtext(side=2,'Productivity Density',outer=TRUE,line=3,cex=1.25)
  mtext(side=1,'log(Size [m])',outer=TRUE,line=3.5,cex=1.25)
  
	par(fig = c(0, 1, 0, 1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = TRUE) 
	plot(0,0, type = 'n', axes = F, bty = 'n')

	legend(horiz=TRUE,'top',cex=1.2,lwd=1.3,legend=c(''),lty=2,bty='n',y=-1)

	legend(horiz=TRUE,'top',cex=1.2,lwd=1.3,legend=c(expression(italic('Primary Prod. (62.7)')),
                                      expression(italic('Trophic Transfer (61.9)')),
                                      expression(italic('Prey Size Range (226.1)')),
                                      expression(italic('Mean Prey Size (355.2)')),
									  expression(italic('Size-Dep. Prod. (288.6)'))),
                   lty=c(1,1,1,1,1),bty='n',col=c('blue','red','orange','dark green','purple'))

  
  
 
  