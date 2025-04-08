# main file

library(discreteRV)
library(abind)

# load data (matrix of baseline covariates X; matrices of time dependent covariates Z1,Z2,...; vector of observed states and times)

# create the array of time dependnet covariates
m=2 # number of time varying covariates
TVC=array(NA,dim=c(n,K+1,m))
TVC[,,1]=as.matrix(Z1)
TVC[,,2]=as.matrix(Z2)

approx_order=3  #order of polynomial approximation of time varying intercept (if 0 it means constant intercept) for the transition probabilities 1->2, 1->3
d1=dim(X)[2]
d2=(m+approx_order)
initial_par= rep(0,2*d2+3*d1) # initial parameter for the optimization procedure

# perform maximum likelihood estimation
DMS_fit=DMS(baseline_cov=X[,-1],TVC=TVC,approx_order=approx_order,obs_State,obs_time,optim_method ="BFGS",maxit=100,initial_par)
est=c(DMS_fit$alpha_0,DMS_fit$beta_12,DMS_fit$beta_13,DMS_fit$gamma_12,DMS_fit$gamma_13)

# repeat the previous procedure with different random starting points to ensure identification of the global maximum

# bootstrap confidence intervals

conf_lev=0.05
B=200 # number of bootstrap samples

est_boot=matrix(NA,nrow = B,ncol=(3*d1+2*d2))
j=1
while(j<=B){
    index_sample=sample(1:n,size=n,replace = TRUE)
    tryCatch({
      DMS_fit=DMS(X[index_sample,-1],TVC[index_sample,,],obs_State[index_sample],obs_time[index_sample],optim_method="BFGS",maxit=100,initial_par=est)
      est_boot[j,]=c(DMS_fit$alpha_0,DMS_fit$beta_12,DMS_fit$beta_13,DMS_fit$gamma_12,DMS_fit$gamma_13)
    }, error = function(e){j<<- j-1}, finally ={j<<- j+1})
  }

boot_sd=apply(est_boot,2,sd)
L=est-boot_sd*qnorm(1-conf_lev/2)
U=est+boot_sd*qnorm(1-conf_lev/2)

# results dataframe
results=data.frame(cbind(est,L,U,boot_sd))
colnames(results)=c("par_est","Lower bound CI", "Upper bound CI","standard error")


