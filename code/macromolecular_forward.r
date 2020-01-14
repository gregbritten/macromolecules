
library(deSolve)
library(rstan)
library(viridis)

theta <- list(KN   = 0.02,    #half-saturation constant 
              mu   = 0.234,   #max protein synthesis rate
              CHsyn= 2.19,    #max photosynthetic rate
              m_ex = 1.43,    #maximum excretion rate
              R_ex = 13.49,   #CN above which excretion happens
              tau  = 9.59,    #photoacclimation time
			  b    = 0.05705) #change in Chl:N relative to N:C

x <- c(CH    = 0.112,
       PR    = 0.0169,
       theta = (4.1/1E6)/0.0169,
	   N     = 0.0911)

dt <- 0.1
n_0  <- 5
n_t  <- 25
t  <- seq(5,25,dt)

dXdt <- function(t,x,theta){
	with(as.list(c(x,theta)),{

	PRsynth = mu*N/(KN + N);                     #protein synthesis
    theta0  = b*(PR/CH);                         #target Chl/PR ratio
    Chl     = theta*PR;                          #cellular Chl
    Rcell   = CH/PR;                             #C:N ratio
    excr    = 0.5*m_ex*(1+tanh(Rcell - R_ex));   #excretion of excess C
    
    dCH     = PR*(CHsyn - excr);                 #Euler CH increment
    dtheta  = (1/tau)*(theta0-theta);            #theta increment
    dPR     = PR*PRsynth;                        #PR increment            
    dN      = -dPR/(1+exp(-10000*N));            #N increment 
	
	list(c(dCH, dPR, dtheta, dN))
	})
}

for (i in 1

X <- as.data.frame(ode(y=x, times=t, func=dXdt, parms=theta))

par(mfrow=c(2,2),mar=c(2,2,2,2))
for(i in 2:5) plot(t,X[,i])




