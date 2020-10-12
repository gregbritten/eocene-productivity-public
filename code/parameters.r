
Ni      <- 512   #number of size classes
imin    <- 10^-7  #minimum size
imax    <- 10     #maximum size
Np      <- 100    #number of phytoplankton size classes
L       <- 5      #number of trophic levels
is      <- seq(log(imin),log(imax),length.out=Ni) #uniform spacing in log space

P1      <- c(dnorm(c(1:100),mean=50,sd=5),rep(0,Ni-100)) #primary production distribution across size classes
alpha   <- 0.2
dfrac   <- 0.1
sigma   <- 0.1

input   <- list(Ni=Ni,is=is,P1=P1)
par     <- c(alpha,dfrac,sigma)

alphas  <- 0.07  #initial alpha
Ps      <- 25591 #initial P
dfracs  <- 0.2   #initial dfrac
sigmas  <- 1.8   #initial sigma
Ps_size <- Ps    #same initial for Ps_size as for Ps

P0      <- 65000 #reference P
alpha0  <- 0.094 #reference alpha
dfrac0  <- 0.096 #reference dfrac
sigma0  <- 0.236 #reference sigma0

size_min <- -9 #min size for sum of squares
size_max <- 1  #max size for sum of squares

loga     <- -2.35 #Power law pre-factor
b        <- 0.80   #Power law power

bw       <- 0.05	