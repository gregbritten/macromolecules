library(rstan)
library(lattice)
options(mc.cores = parallel::detectCores())

mod_code <- "functions {
  real[] macro(real   t,           // time
               real[] x,           // state x[1]:CH  x[2]:PR, x[3]:Chl , x[4]:N
               real[] theta,
               real[] x_r,
               int[]  x_i) {       // parameters

    real KN    = theta[1];    
    real mu    = theta[2]; 
    real CHsyn = theta[3]; 
    real m_ex  = theta[4];  
    real R_ex  = theta[5];  
    real tau   = theta[6];
    real b     = theta[7];
    
    real CH = fmax(x[1],0.0);
    real PR = fmax(x[2],0.0);
    real Chl = fmax(x[3],0.0);
    real N = fmax(x[4],0.0);

    real PRsynth = theta[2]*N/(theta[1]+ pow(theta[1]*N,0.5) + N);
    real r0      = theta[7]*(PR/CH);
    //real Chl     = Chl*P;
    real Rcell   = x[1]/x[2];
    real excr    = 0.5*theta[4]*(1+tanh(Rcell - theta[5]));
    
    real dCH    = x[2]*(theta[3] - excr);
    real dr     = (1/theta[6])*(r0-x[3]);
    real dPR    = x[2]*PRsynth;
    real dN     = -dPR/(1+exp(-10000*x[4])); 

    return {dCH,dPR,dr,dN};
  }
}
data {
  int<lower = 0> n;           // num obs
  real t_obs[n];              // obs times
  real<lower = 0> y[n,4];     // observed variable at measurement times
}
parameters {
  real<lower = 0> theta[7];   // parameters
  real<lower = 0> x0[4];      // initial population
  real<lower = 1E-15> sigma[4]; // obs error
}
transformed parameters {
  real<lower=1E-15> x[n,4] = integrate_ode_rk45(macro, x0, 0, t_obs, theta, rep_array(0.0,0), rep_array(0,0), 1e-6, 1e-5, 1e3) ;
  for(i in 1:n){
    x[i,3] = x[i,3]*x[i,2]*1E6;
  }
}
model {
  x0[1]    ~ normal(0.1,1);
  x0[2]    ~ normal(0.1,1);
  x0[3]    ~ normal(10,10);
  x0[4]    ~ normal(0.1,1);  
  theta[1] ~ normal(0.002,3);
  theta[2] ~ normal(0.3,1);
  theta[3] ~ normal(5,10);
  theta[4] ~ normal(10,10);
  theta[5] ~ normal(13,10);
  theta[6] ~ normal(10,10); 
  theta[7] ~ normal(0.05,5); 
  
  for(i in 1:4){
    y[1:n,i] ~ normal(x[1:n,i], sigma[i]);
  }
}"

dat       <- read.csv('data/macromolecules.csv')

dat_flynn <- read.csv('data/flynn_macromolecules.csv')

dat[,2] <- dat[,2]*0.001                 #PR;
dat[,3] <- dat[,3]*1000                  #Chl;
dat[,4] <- dat[,4]*0.001                 #C
dat     <- merge(dat,dat_flynn[,c(1,4)]) #add environmental nitrogen to the data
dat[,5] <- dat[,5]*0.08325909            

head(dat)

options(repr.plot.width=6, repr.plot.height=4)
par(mfrow=c(2,2),mar=c(2,3,2,2),oma=c(2,2,2,2))
for(i in 2:5){
    plot(dat[,1],dat[,i],ylab='',xlab='')
    mtext(side=2,colnames(dat)[i],line=2.5)
}

data <- list(n    =nrow(dat),
             t_obs=dat[,1],
             y    =dat[,c(4,2,3,5)])

mod <- stan_model(model_code=mod_code)

mcmc <- sampling(mod,data=data,iter=2000,chains=4,open_progress=TRUE)

post <- extract(mcmc)

mcmc

post <- extract(mcmc)

plot_theta <- function(post,theta_names){   #theta_names has to be in the the same order 
	for(i in 1:dim(post$theta)[2]){
		hist(post$theta[,i],main='')
		abline(v=mean(post$theta[,i]),lwd=2)
		mtext(theta_names[i])
	}
}
plot_CV <- function(post,state_names){
	for(i in 1:dim(post$x)[3]){
		hist(post$sigma[,i]/mean(post$x[,,i]),main='')
		mtext(state_names[i])
	}; mtext(outer=TRUE,expression(sigma['x']/mu['x']),line=-1,cex=2)
}
plot_state <- function(post,state_names){
	for(i in 1:dim(post$x)[3]){
		plot(data$t_obs,colMeans(post$x[,,i]),type='l',bty='n')
		lines(data$t_obs,apply(post$x[,,i],2,function(x) quantile(x,probs=0.025)),lty=2)
		lines(data$t_obs,apply(post$x[,,i],2,function(x) quantile(x,probs=0.975)),lty=2)
		points(data$t_obs,data$y[,i])
		mtext(state_names[i])
	}
}	
plot_trace <- function(post,theta_names){
	for(i in 1:dim(post$theta)[2]){
		plot(post$theta[,i],pch=19,cex=0.2)
		abline(h=mean(post$theta[,i]),lwd=2)
		mtext(theta_names[i])
	}
}                           

theta_names <- c('KN','mu','CHsyn','m_ex','R_ex','tau','b')
state_names <- c('CH','PR','Chl','N')

options(repr.plot.width=6, repr.plot.height=6)
par(mfrow=c(4,2),mar=c(2,2,2,2),oma=c(2,2,2,2))
	plot_theta(post=post,theta_names=theta_names)

options(repr.plot.width=6, repr.plot.height=6)
par(mfrow=c(4,2),mar=c(2,2,2,2),oma=c(2,2,2,2))
	plot_trace(post=post,theta_names=theta_names)

options(repr.plot.width=6, repr.plot.height=4)
par(mfrow=c(2,2),oma=c(2,2,2,2),mar=c(2,2,2,2))
	plot_CV(post=post,state_names=state_names)

par(mfrow=c(2,2),mar=c(2,2,2,2),mar=c(2,2,2,2))
	plot_state(post=post,state_names=state_names)

options(repr.plot.width=6, repr.plot.height=6)
pairs(post$theta[sample(1:nrow(post$theta),size=500),],pch=16,cex=0.5,labels=theta_names)

mod_code2 <- "functions {
  real[] macro(real   t,           // time
               real[] x,           // state x[1]:CH  x[2]:PR, x[3]:Chl , x[4]:N
               real[] theta,
               real[] x_r,
               int[]  x_i) {       // parameters

    // real KN    = theta[1];    
    real mu    = theta[1]; 
    real CHsyn = theta[2]; 
    real m_ex  = theta[3];  
    real R_ex  = theta[4];  
    real tau   = theta[5];
    real b     = theta[6];

    real PRsynth = theta[1]*x[4]/(x_r[1]+x[4]);
    real r0      = theta[6]*(x[2]/x[1]);
    real Chl     = x[3]*x[2];
    real Rcell   = x[1]/x[2];
    real excr    = (1.0/2.0)*theta[3]*(1+tanh(Rcell - theta[4]));  //avoid integer division
    
    real dCH    = x[2]*(theta[2] - excr);
    real dr     = (1/theta[5])*(r0-x[3]);
    real dPR    = x[2]*PRsynth;
    real dN     = -dPR/(1+exp(-10000*x[4])); 

    return {dCH,dPR,dr,dN};
  }
}
data {
  int<lower = 0> n;           // num obs
  real t_obs[n];              // obs times
  real<lower = 0> y[n,4];     // observed variable at measurement times
}
parameters {
  real<lower = 0> theta[6];   // parameters
  real<lower = 0> x0[4];      // initial population
  real<lower = 1E-15> sigma[4]; // obs error
}
transformed parameters {
  real x[n,4] = integrate_ode_rk45(macro, x0, -1, t_obs, theta, {0.02}, rep_array(0,0), 
    1e-6, 1e-5, 1e3) ;
  for(i in 1:n){
    x[i,3] = x[i,3]*x[i,2]*1E6;
  }
}
model {
  x0[1]    ~ normal(0.1,1);
  x0[2]    ~ normal(0.1,1);
  x0[3]    ~ normal(10,10);
  x0[4]    ~ normal(0.1,1);  
  // theta[1] ~ normal(0.002,3);
  theta[1] ~ normal(0.3,1);
  theta[2] ~ normal(5,10);
  theta[3] ~ normal(10,10);
  theta[4] ~ normal(13,10);
  theta[5] ~ normal(10,10); 
  theta[6] ~ normal(0.05,5); 
  
  for(i in 1:4){
    y[1:n,i] ~ normal(x[1:n,i], sigma[i]);
  }
}"

mod2 <- stan_model(model_code=mod_code2)

mcmc2 <- sampling(mod2,data=data,open_progres=TRUE)

post2 <- extract(mcmc2)

mcmc2

pairs(post$theta[sample(1:4000,200),])
pairs(post2$theta[sample(1:4000,200),])

options(repr.plot.width=8, repr.plot.height=2)
par(mfrow=c(1,3),mar=c(2,2,2,2))
for(i in 1:3){
    hist(post$theta[,i+1],col=adjustcolor('blue',alpha.f=0.3),main='')
    hist(post2$theta[,i],add=TRUE,col=adjustcolor('red',alpha.f=0.3))
}

