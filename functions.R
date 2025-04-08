# required functions


# multinomial logistic functions
ml12=function(g12,g13,x){
  r1=g12%*%t(x)
  r2=g13%*%t(x)
  r=exp(r1)/(1+exp(r1)+exp(r2))
  return(r)
}
ml13=function(g12,g13,x){
  r1=g12%*%t(x)
  r2=g13%*%t(x)
  r=exp(r2)/(1+exp(r1)+exp(r2))
  return(r)
}
# logistic function
p=function(gamma,x){
  r=gamma%*%t(x)
  r=1/(1+exp(-r))
  return(r)}


# main function DMS requiring in input: 
# baseline_cov = matrix (n\times d1) of baseline covariates 
# TVC = array of time dependent covariates (n\times K+1 \times m) where m is the number of time varying covariates
# approx_order = order of polynomial approximation of time varying intercept (if 0 it means constant intercept) for the transition probabilities 1->2, 1->3
# obs_State = vector (length n) of observed states for each subject, takes value 3 if event is observed and NA otherwise (meaning subject could be in state 1 or 2)
# obs_time = vector (length n) of observed event times for each subject, takes value in {0,...,K-1} if event (transition 1->3) happens during (t_i,t_{i+1}] and NA otherwise
# optim_method = optimization method used for likelihood maximization in the optim function (can be "BFGS","Nelder-Mead","SANN"...)
# maxit = maximum number of iterations for the optimization algorithm
# initial_par = initial value of the parameters for the optimization procedure
# The output of the function contains the parameter estimates, log-likelihood value and whether the optimization procedure converged (conv=1 mean successful)

DMS=function(baseline_cov,TVC,approx_order,obs_State,obs_time,optim_method="BFGS",maxit=100,initial_par){
  
  n=dim(baseline_cov)[1]
  X=cbind(rep(1,n),baseline_cov)
  Z=TVC
  m=dim(TVC)[3] # number of time varying covariates
  K=dim(TVC)[2]-1
  
  d1=dim(X)[2] # dimension of alpha and beta parameters for baseline covariates
  d2=m+approx_order #dimension of the gamma parameters for time-dependent covariates (including the deterministic polynomial terms of time)
  for(r in 1:approx_order){
    U=matrix((0:K)^r,nrow = n,ncol=K+1,byrow = TRUE)
    Z=abind(Z,array(NA,dim = c(n,K+1,1)))
    Z[,,m+r]=U
  }
  
  
  log_Lik=function(par){
    
    
    alpha_0=par[1:d1]
    
    beta_12=par[(d1+1):(2*d1)]
    beta_13=par[(2*d1+1):(3*d1)]
    gamma_12=par[(3*d1+1):(3*d1+d2)]
    gamma_13=par[(3*d1+d2+1):(3*d1+2*d2)]
    
    pi_1=c(p(alpha_0,X))
    
    
    l=0*(1:n)
    
    for(i in 1:n){
      if(is.na(obs_State[i])){ # censored case
        k=obs_time[i]
        l1=(1-pi_1[i]) #if starts in 2 stays in 2
        l2=pi_1[i] #starts in 1 and stays in 1
        {
          for(h in 0:(k)){  
            p_11=1/c(1+exp(c(beta_12,gamma_12)%*%(c(X[i,],Z[i,h+1,])))+exp(c(beta_13,gamma_13)%*%(c(X[i,],Z[i,h+1,]))))
            l2=l2*p_11
          }
        }
        l3=pi_1[i] #starts in 1 and moves to 2, then stays in 2
        temp=0
        {
          for(r in 0:(k)){
            p_12=c(exp(c(beta_12,gamma_12)%*%(c(X[i,],Z[i,r+1,]))))/c(1+exp(c(beta_12,gamma_12)%*%(c(X[i,],Z[i,r+1,])))+exp(c(beta_13,gamma_13)%*%(c(X[i,],Z[i,r+1,]))))
            product=p_12
            if(r>0){
              for(h in 0:(r-1)){
                p_11=1/c(1+exp(c(beta_12,gamma_12)%*%(c(X[i,],Z[i,h+1,])))+exp(c(beta_13,gamma_13)%*%(c(X[i,],Z[i,h+1,]))))
                product=product*p_11
              }
            }
            temp=temp+product
          }
        }
        l3=l3*temp
        { l[i]=log(l1+l2+l3)}
      }
      if(is.na(obs_State[i])==FALSE){ # so observed mover
        k=obs_time[i]
        p_13=c(exp(c(beta_13,gamma_13)%*%(c(X[i,],Z[i,k+1,]))))/c(1+exp(c(beta_12,gamma_12)%*%(c(X[i,],Z[i,k+1,])))+exp(c(beta_13,gamma_13)%*%(c(X[i,],Z[i,k+1,]))))
        l[i]=l[i]+log(pi_1[i])+log(p_13)
        if(k>0){
          h=0
          while(h<=k-1){
            p_11=1/c(1+exp(c(beta_12,gamma_12)%*%(c(X[i,],Z[i,h+1,])))+exp(c(beta_13,gamma_13)%*%(c(X[i,],Z[i,h+1,]))))
            l[i]=l[i]+log(p_11)
            h=h+1
          }
        }
      }
    }
    return(sum(l))
  }
  
  MLE=optim(initial_par,log_Lik,method=optim_method,control=list(fnscale=-1,maxit=maxit))
  
  alpha_0=MLE$par[1:d1]
  beta_12=MLE$par[(d1+1):(2*d1)]
  beta_13=MLE$par[(2*d1+1):(3*d1)]
  gamma_12=MLE$par[(3*d1+1):(3*d1+d2)]
  gamma_13=MLE$par[(3*d1+d2+1):(3*d1+2*d2)]
  
  
  conv=MLE$convergence
  
  return(list(alpha_0=alpha_0,beta_12=beta_12,beta_13=beta_13,gamma_12=gamma_12,gamma_13=gamma_13,conv=conv,log_lik=MLE$value))
  
}

