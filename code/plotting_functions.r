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
