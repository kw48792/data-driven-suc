###############################################
## Generate posterior samples of the IMS-VAR model
## data: input data
## warm_up: # of warm up iterations
## num_iter: # of posterior samples generated
## alpha,eta,kappa: hyper parameters for (sticky) HDP
## p: VAR lag of order
## Output: A list of posterior samples, each sample
##    is specified by:
##    count: # of obs in each state (the length varies, indicating the numeber of states)
##    Psi,sigma: VAR parameters
##    S: state vectors (S_1,...,S_TT)
##    N: transition matrix (N_ij is the # of tranition from i to j)
##    resid: prediction residues
###############################################
IMS_AR_P = function(data,warm_up,num_iter,alpha,eta,kappa,p){
  results = list()
  x = data
  TT = length(data)
  #N=matrix(0,round(TT/3),round(TT/3))
  #count = c(rep(3,floor(TT/3)),TT%%3)
  #K = dim(N)[1]
  #pi = rep(1/(K+1),K+1)
  count = rep(1,TT)
  S = seq(1,TT,1)
  N=matrix(0,TT,TT)
  K = TT
  pi = rep(1/(TT+1),TT+1)
  assign = list(c(1))
  mu_psi = 0
  sigma_psi = 1
  aa = 1
  bb = 1
  Psi = matrix(c(rnorm(1,mu_psi,sigma_psi),rtnorm(p,0,1,-1,1)),1,p+1)
  sigma = c(rigamma(1,aa/2,bb/2))
  resid = rep(0,TT)
  for(t in 2:TT){
    N[t-1,t] = 1
    if(t<TT){
      N[t,t+1] = 1
    }
    assign[[length(assign)+1]] = c(t)
    Psi = rbind(Psi,c(rnorm(1,mu_psi,sigma_psi),rtnorm(p,mu_psi,sigma_psi,-1,1)))
    sigma = c(sigma,rigamma(1,aa/2,bb/2))
  }
  for(rep in 1:(warm_up+num_iter)){
    #print(c(rep,K))
    for(t in (p+1):TT){
      N[S[t-1],S[t]] = N[S[t-1],S[t]]-1
      count[S[t]] = count[S[t]]-1
      if(t<TT){
        N[S[t],S[t+1]] = max(N[S[t],S[t+1]]-1,0)
      }
      assign[[S[t]]] = setdiff(assign[[S[t]]],t)
      post_prob = rep(0,K+1)
      for(k in 1:K){
        if(t<TT){
          pp = (eta*pi[k]+N[S[t-1],k]+kappa*(k==S[t-1]))*(eta*pi[S[t+1]]+N[k,S[t+1]]+kappa*(k==S[t+1]))/
            (eta+count[k]+kappa)*dnorm(x[t],sum(c(1,x[(t-p):(t-1)])%*%Psi[k,]),sigma[k])
        }
        else{
          pp = (eta*pi[k]+N[S[t-1],k]+kappa*(k==S[t-1]))*dnorm(x[t],sum(c(1,x[(t-p):(t-1)])%*%Psi[k,]),sigma[k])
          #print(dnorm(x[t],sum(c(1,x[t-1])*Psi[k])))
        }
        post_prob[k] = pp 
      }
      new_Psi = c(rnorm(1,mu_psi,sigma_psi),rtnorm(p,mu_psi,sigma_psi,-1,1))
      #new_Psi = c(rnorm(1,mu_psi[1],sigma_psi[1]),0)
      new_sigma = sqrt(rigamma(1,aa/2,bb/2))
      pp = (eta*pi[K+1])*dnorm(x[t],sum(c(1,x[(t-p):(t-1)])%*%new_Psi),new_sigma)
      if(t<TT){
        pp = pp*eta*pi[S[t+1]]/(eta+kappa)
      }
      post_prob[K+1] = pp
      if(sum(post_prob)==0){
        post_prob = rep(1,K+1)
      }
      post_prob = post_prob/sum(post_prob)
      S[t] = sample(1:(K+1),1,replace = F, post_prob)
      if(S[t] == K+1){
        N = cbind(N,matrix(rep(0,K),K,1))
        N = rbind(N,rep(0,K+1))
        count = c(count,0)
        assign[[K+1]] = c(t)
        Psi = rbind(Psi,new_Psi)
        sigma = c(sigma,new_sigma)
        b = rbeta(1,alpha,1)
        pi=c(pi,pi[K+1]*(1-b))
        pi[K+1] = pi[K+1]*b
        K = K+1
      }
      else{
        assign[[S[t]]] = c(assign[[S[t]]],t)
      }
      count[S[t]] = count[S[t]]+1
      N[S[t-1],S[t]] = N[S[t-1],S[t]]+1
      if(t<TT){
        N[S[t],S[t+1]] = N[S[t],S[t+1]]+1
      }
      ind = assign[[S[t]]]
      ind = ind[ind>p]
      data_in = x[ind]
      data_prev = matrix(0,length(data_in),p)
      for(i in 1:length(data_in)){
        for(j in 1:p){
          data_prev[i,j] = x[ind[i]-j]
        }
      }
      size_in = length(data_in)
      curSigma = sigma[S[t]]
      mean_psi_0 = (sigma_psi^2*sum(data_in-data_prev %*% Psi[S[t],2:(p+1)])+curSigma^2*mu_psi)/(curSigma^2+sigma_psi^2*size_in)
      var_psi_0 = (curSigma^2*sigma_psi^2)/(curSigma^2+sigma_psi^2*size_in)
      Psi[S[t],1] = rnorm(1,mean_psi_0,sqrt(var_psi_0))
      err = data_in-Psi[S[t],1]-data_prev %*% Psi[S[t],2:(p+1)]
      for(j in 1:p){
        lag_data = x[ind-j]
        mean_psi = (sigma_psi^2*sum(lag_data*(err+lag_data*Psi[S[t],j+1]))+curSigma^2*mu_psi)/(curSigma^2+sigma_psi^2*sum(lag_data^2))
        var_psi = (curSigma^2*sigma_psi^2)/(curSigma^2+sigma_psi^2*sum(lag_data^2))
        Psi[S[t],(j+1)] = rtnorm(1,mean_psi,sqrt(var_psi),-1,1)
      }
      sigma[S[t]] = sqrt(rigamma(1,(aa+size_in)/2,(bb+sum(err^2))/2))
      resid[t] = data[t]-data[(t-p):(t-1)] %*% Psi[S[t],2:(p+1)]-Psi[S[t],1]
    }
    
    removed = c()
    remove = 0
    K_removed = c()
    for(k in 1:K){
      if(count[k]==0){
        K = K-1
        K_removed = c(K_removed,k)
        remove = remove+1
      }
      removed = c(removed,remove)
    }
    if(length(K_removed) > 0){
      N = N[-K_removed,-K_removed]
      if(!is.matrix(N)){
        N = as.matrix(N)
      }
      assign = assign[-K_removed]
      Psi = matrix(Psi[-K_removed, ],K,(p+1))
      sigma = sigma[-K_removed]
      count = count[-K_removed]
    }
    for(tt in 1:TT){
      S[tt] = S[tt]-removed[S[tt]]
    }
    paraDiri =c(alpha*apply(N,1,sum),1) 
    paraDiri[1] = paraDiri[1]+1*alpha
    weights = rgamma(K+1,paraDiri,1)
    pi = weights/sum(weights)
    if(rep>warm_up){
      results[[length(results)+1]] = list(count,Psi,sigma,S,N,resid) 
    }
  }
  return(results)
}

getPrediction_IMSAR = function(prev, p, prev_S, P, Psi, sigma, step){
  cur = prev
  res = rep(0,step)
  for(i in 1:step){
    prob = P[prev_S,]
    prob[is.na(prob)] = 1e-5
    prob = prob/sum(prob)
    cur_S = sample(1:dim(P)[1],1,replace = F,prob)
    prev_S = cur_S
    pred = sum(c(1,cur)*Psi[cur_S,])+rnorm(1,0,sigma[cur_S])
    if(p == 1){
      cur = c(pred)
    }
    else{
      cur = c(pred,cur[1:(p-1)])
    }
    res[i] = pred
  }
  return(res)
}

#######################################
## Determine the wind power data contains ramp event or not
#######################################
isRamp = function(data){
  TT = length(data)
  for(i in 1:(TT-1)){
    if(abs(data[i+1]-data[i])>0.25){
      return (TRUE)
    }
  }
  return (FALSE)
}