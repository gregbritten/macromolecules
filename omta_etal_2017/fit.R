library(rstan)
options(mc.cores = parallel::detectCores())
library(zoo)

dat_flynn <- read.csv('d:/dropbox/data/Flynn_etal_1993/table.csv')

dat <- read.csv('d:/dropbox/teaching/stan_examples/macromolecular/omta_etal_2017/data_omta.csv')
dat[,2] <- dat[,2]*0.001 #PR 
dat[,3] <- dat[,3]*1000  #Chl
dat[,4] <- dat[,4]*0.001 #C
dat <- merge(dat,dat_flynn[,c(1,4)])
dat[,5] <- dat[,5]*0.08325909

par(mfrow=c(2,2))
plot(dat[,1],dat[,2])
plot(dat[,1],dat[,3])
plot(dat[,1],dat[,4])
plot(dat[,1],dat[,5])


data <- list(n    =nrow(dat),
             t_obs=dat[,1],
             y    =dat[,c(4,2,3,5)])


#########################################################################
## FIT MODEL ############################################################
#########################################################################
mod <- stan_model('d:/dropbox/teaching/stan_examples/macromolecular/omta_etal_2017/mod.stan')

mcmc <- sampling(mod,
	data=data,
	iter=2000,chains=4,
	open_progress=TRUE)
	
	
post <- extract(mcmc)

theta_names <- c('CNpro','KN','mu','CHsyn','m_ex','R_ex','tau','b')
par(mfrow=c(4,2),mar=c(2,2,2,2))
for(i in 1:8){
	hist(post$theta[,i],main='')
	mtext(theta_names[i])
}

par(mfrow=c(2,2))
for(i in 1:4){
	plot(data$t_obs,colMeans(post$x[,,i]),type='l')
	lines(data$t_obs,colMeans(post$x[,,i]) + 2*apply(post$x[,,i],2,sd))
	lines(data$t_obs,colMeans(post$x[,,i]) - 2*apply(post$x[,,i],2,sd))
	points(data$t_obs,data$y[,i])
}


