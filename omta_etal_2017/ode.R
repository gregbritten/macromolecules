library(deSolve)


parameters <- c(CNpro=6.6, #CN ratio of protein
                KN=0.33,   #half-saturation constant  
                mu=0.3,    #max protein synthesis rate
                CHsyn=10,  #max photosynthetic rate
                m_ex=2,    #maximum excretion rate
                R_ex=10,   #CN above which excretion happens
                tau=1,     #photoacclimation time
                b=1,
                Nin=1)       #change in Chl:N relative to N:C

state <- c(CH    = 0.1242,
           PR    = 0.0186,
           theta = 15/0.0186,
           N     = 1.411)

dxdt <- function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    
    PRsynth = mu*N/(KN+N)
    theta0  = b*(PR/CH)
    Chl     = theta*PR
    Rcell   = CH/(PR + Chl)
    excr    = (1/2)*m_ex*(1-tanh(Rcell - R_ex))
    
    dCH    = PR*(CHsyn - CNpro*PRsynth - excr)
    dtheta = (1/tau)*(theta0-theta)
    dPR    = PR*PRsynth
    dN     = (1/(1+exp(-1000*N)))*(-dCH + -dPR) 
    
    list(c(dCH,dPR,dtheta,dN))})}

times <- seq(0, 25, by = 0.01)

out <- as.data.frame(ode(y=state, times=times, func=dxdt, parms=parameters))

par(mfrow=c(2,2),mar=c(2,2,2,2))
plot(out[,1],out[,2],type='l')
plot(out[,1],out[,3],type='l')
plot(out[,1],out[,4],type='l')
plot(out[,1],out[,5],type='l')


#################################################################
##
################################################################

parameters <- c(CNpro=6.6, #CN ratio of protein
                KN=0.02,   #half-saturation constant  
                mu=0.234,    #max protein synthesis rate
                CHsyn=5,  #max photosynthetic rate
                m_ex=3,    #maximum excretion rate
                R_ex=13.49,   #CN above which excretion happens
                tau=9.59,     #photoacclimation time
                b=57/(1E6))   #change in Chl:N relative to N:C

state <- c(CH    = 0.1242,
           PR    = 0.0169,
           theta = (41/1E6)/0.0169,
           N     = 0.09)

dxdt <- function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    
    PRsynth = mu*N/(KN+N)
    theta0  = b*(PR/CH)
    Chl     = theta*PR
    Rcell   = CH/(PR + Chl)
    excr    = (1/2)*m_ex*(1-tanh(Rcell - R_ex))
    
    dCH    = PR*(CHsyn - CNpro*PRsynth - excr)
    dtheta = (1/tau)*(theta0-theta)
    dPR    = PR*PRsynth
    dN     = (1/(1+exp(-1000*N)))*(-dPR - (dtheta*dPR)) 
    
    list(c(dCH,dPR,dtheta,dN))})}

times <- seq(0, 25, by = 0.01)

out <- as.data.frame(ode(y=state, times=times, func=dxdt, parms=parameters))

par(mfrow=c(2,2),mar=c(2,2,2,2))
plot(out[,1],out[,2],type='l')
plot(out[,1],out[,3],type='l')
plot(out[,1],out[,4]*out[,3]*1E6,type='l')
plot(out[,1],out[,5],type='l')

###########################################################################
##
###########################################################################
