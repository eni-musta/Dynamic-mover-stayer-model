# Data generation


n=1000 # sample size

X=matrix(0,n,3) # baseline covariate matrix 
X[,1]=X[,1]+1   # intercept
X[,2]=rnorm(n,0,1)
X[,3]=rbinom(n,1,0.4)

T_alpha_0=c(0.8,0.5,-1) # coefficients for initial probabilities of being susceptible at baseline (starting in state 1)
pi_1=c(p(T_alpha_0,X)) # initial probabilities for starting in state 1 (potential mover)

# initial mover/stayer status
B = c(runif(n) <= 1-pi_1)+1 # B=1 starts in state 1 (potential mover), B=2 starts in state 2 (stayer)

initial_stayer_prop=length(which(B==2))/n  # proportion of stayers at time 0

T_gamma_12=c(0.11,-0.2) # coefficients for transitions 1 -> 2 corresponding to time varying covariates
T_gamma_13=c(-0.5,0.3) # coefficients for transitions 1 -> 3 corresponding to time varying covariates

T_beta_12=c(-1,0.6,-0.1) # coefficients for transitions 1 -> 2 corresponding to the baseline covariates
T_beta_13=c(-2,-0.4,0.1) # coefficients for transitions 1 -> 3 corresponding to the baseline covariates

# intercept term is the first component of the beta coefficients


K=5 #years of follow-up

S=matrix(NA,nrow=n,ncol=K+1) # matrix of states for each observation and each time, column i corresponds to the start of year i (so 1st columns is for baseline t0=0=start of year 1)
S[,1]=B


Z1=matrix(0,nrow = n,ncol=K+1) # first time-dependent covariate
Z2=matrix(0,nrow = n,ncol=K+1) # second time-dependent covariate
Z2[,1]=sample(1:5,n,replace = TRUE)

for (k in 2:(K+1)){
  
  # transition probabilities during (t_k-2,t_k-1]
  P_11=1/c((1+exp(c(T_beta_12,T_gamma_12)%*%t(cbind(X,Z1[,k-1],Z2[,k-1])))+exp(c(T_beta_13,T_gamma_13)%*%t(cbind(X,Z1[,k-1],Z2[,k-1])))))
  P_12=c(exp(c(T_beta_12,T_gamma_12)%*%t(cbind(X,Z1[,k-1],Z2[,k-1]))))/c(1+exp(c(T_beta_12,T_gamma_12)%*%t(cbind(X,Z1[,k-1],Z2[,k-1])))+exp(c(T_beta_13,T_gamma_13)%*%t(cbind(X,Z1[,k-1],Z2[,k-1]))))
  P_13=c(exp(c(T_beta_13,T_gamma_13)%*%t(cbind(X,Z1[,k-1],Z2[,k-1]))))/c(1+exp(c(T_beta_12,T_gamma_12)%*%t(cbind(X,Z1[,k-1],Z2[,k-1])))+exp(c(T_beta_13,T_gamma_13)%*%t(cbind(X,Z1[,k-1],Z2[,k-1]))))
  
  transit_1=rep(1,n)
  u = runif(n)
  transit_1[which(u >P_11 & u<=(P_11+P_12))]=2
  transit_1[which(u >(P_11+P_12))]=3
  
  S[which(S[,k-1]==2),k]=2
  S[which(S[,k-1]==1),k]=transit_1[which(S[,k-1]==1)]
  S[which(S[,k-1]==3),k]=3  # if the state was 3 we leave as it is, we are not interested in what happens after
  
  
  Z1[,k]=Z1[,k-1]+rnorm(n,0.5,1)
  Z2[,k]=Z2[,k-1]+rbinom(n,2,0.5)-1
  
  
}
transition_time=rep(NA,n) # time transition to 3 (if it happens) and NA otherwise; if equal to i, it means that the transition happened during (t_i,t_i+1]

for(i in 1:n){
  if(S[i,K+1]==3){
    transition_time[i]=which(S[i,]==3)[1]-2
  }
}

# censoring

C=1+floor(pmin(rexp(n,0.03),rep(K-2,n))) #censoring <=t_{K}; if equal to i means that censoring happened in (t_i,t_{i+1}] (at the end of the interval)

# observations

obs_time=rep(NA,n) # observed times of transition 3, or censoring: If equal to i it means that the transition/censoring happened during (t_i,t_i+1]
obs_State=rep(NA,n) # state 1 and 2 are not observed
for(i in 1:n){
  if(is.na(transition_time[i])==FALSE){
    obs_time[i]=min(C[i],transition_time[i])
    if(transition_time[i]<=C[i]){
      obs_State[i]=3
    }
  }
  else {
    obs_time[i]=C[i]
  }
}


# true cumulative probability of being a stayer or moving

baseline_true_stayer_prob=(1-c(p(T_alpha_0,X)))

prob_moving_true=matrix(NA,nrow=n,ncol=K)
prob_stayer_true=matrix(NA,nrow=n,ncol=K)

prob_stayer_true[,1]=(1-c(p(T_alpha_0,X)))+c(p(T_alpha_0,X))*c(ml12(c(T_beta_12,T_gamma_12),c(T_beta_13,T_gamma_13),cbind(X,Z1[,1],Z2[,1])))
prob_moving_true[,1]=c(p(T_alpha_0,X))*c(ml13(c(T_beta_12,T_gamma_12),c(T_beta_13,T_gamma_13),cbind(X,Z1[,1],Z2[,1])))

for(i in 1:n){
  
  for(j in 2:K){
    prob_moving_true[i,j]=prob_moving_true[i,(j-1)]+(1-prob_moving_true[i,(j-1)]-prob_stayer_true[i,j-1])*c(ml13(c(T_beta_12,T_gamma_12),c(T_beta_13,T_gamma_13),t(c(X[i,],Z1[i,j],Z2[i,j]))))
    prob_stayer_true[i,j]=prob_stayer_true[i,(j-1)]+(1-prob_moving_true[i,(j-1)]-prob_stayer_true[i,j-1])*c(ml12(c(T_beta_12,T_gamma_12),c(T_beta_13,T_gamma_13),t(c(X[i,],Z1[i,j],Z2[i,j]))))
  }
}

true_par=c(T_alpha_0,T_beta_12,T_beta_13,T_gamma_12,T_gamma_13) 

df=data.frame(cbind(X[,-1]),Z1,Z2,obs_State,obs_time)
colnames(df)=c("X1","X2","Z1_0","Z1_1","Z1_2","Z1_3","Z1_4","Z1_5","Z2_0","Z2_1","Z2_2","Z2_3","Z2_4","Z2_5","Obs_state","Obs_time")
