#### Coalescent likelihood functions ####
### Adapted to sample Tajima genealogies 
# Compute log likelihood of coalescent model, energy function in HMC algorithms, and metric tensor needed in aMALA.

coal_lik_init = function(samp_times, n_sampled, coal_times, grid)
{
  ns = length(samp_times)
  nc = length(coal_times)
  ng = length(grid)-1
  
  
  if (length(samp_times) != length(n_sampled))
    stop("samp_times vector of differing length than n_sampled vector.")
  
  if (length(coal_times) != sum(n_sampled) - 1)
    stop("Incorrect length of coal_times: should be sum(n_sampled) - 1.")
  
  if (max(samp_times, coal_times) > max(grid))
    stop("Grid does not envelop all sampling and/or coalescent times.")
  
  t = sort(unique(c(samp_times, coal_times, grid)))
  l = rep(0, length(t))
  
  for (i in 1:ns)
    l[t >= samp_times[i]] = l[t >= samp_times[i]] + n_sampled[i]
  
  for (i in 1:nc)
    l[t >= coal_times[i]] = l[t >= coal_times[i]] - 1
  
  #print(l)
  
  if (sum((l < 1) & (t >= min(samp_times))) > 0)
    stop("Number of active lineages falls below 1 after the first sampling point.")
  
  mask = l > 0
  t = t[mask]
  l = utils::head(l[mask], -1)
  
  gridrep = rep(0, ng)
  for (i in 1:ng)
    gridrep[i] = sum(t > grid[i] & t <= grid[i+1])
  
  C = 0.5 * l * (l-1)
  D = diff(t)
  
  y = rep(0, length(D))
  y[t[-1] %in% coal_times] = 1
  
  bins = cut(x = samp_times, breaks = t,
                include.lowest = TRUE)
  tab <- stats::aggregate(n_sampled ~ bins, FUN = sum, labels = FALSE)
  count <- rep(0, length(D))
  count[as.numeric(tab$bins)] <- tab$n_sampled
  count[utils::head(t, -1) >= max(samp_times)] <- NA
  
  rep_idx = cumsum(gridrep)
  rep_idx = cbind(rep_idx-gridrep+1,rep_idx)
  
  return(list(t=t, l=l, C=C, D=D, y=y, count=count, gridrep=gridrep, ns=sum(n_sampled), nc=nc, ng=ng, rep_idx=rep_idx, args=list(samp_times=samp_times, n_sampled=n_sampled, coal_times=coal_times, grid=grid)))
}

coal_loglik = function(init, f, grad=FALSE)
{
  if (init$ng != length(f))
    stop(paste("Incorrect length for f; should be", init$ng))
  
  f = rep(f, init$gridrep)
  
  llnocoal = init$D * init$C * exp(-f)
  
  if (!grad)
  {  
    lls = -init$y * f - llnocoal
    #print(lls)
    
    ll = sum(lls[!is.nan(lls)])
    
    return(ll)
  }
  else
  {  
    dll = apply(init$rep_idx,1,function(idx)sum(-init$y[idx[1]:idx[2]]+llnocoal[idx[1]:idx[2]])) # gradient of log-likelihood wrt f_midpts
    
    return(dll)
  }
}

coal_loglik2 = function(init, f, grad=FALSE)
{ #computes the gradient wrt coalescent times
  if (init$ng != length(f))
    stop(paste("Incorrect length for f; should be", init$ng))
  
  f = rep(f, init$gridrep)
  
  llnocoal = init$D * init$C * exp(-f)
  
  if (!grad)
  {  
    lls = -init$y * f - llnocoal
    #print(lls)
    
    ll = sum(lls[!is.nan(lls)])
    
    return(ll)
  }
  else
  {  
    movey = c(0,init$y[-length(init$y)])
    if (init$y[length(init$y)]==1){
      dll = c((movey*init$C*exp(-f))[movey==1],0) - (init$y *init$C * exp(-f))[init$y==1]
     # dll = - (init$y *init$C * exp(-f))[init$y==1]
      
    }else{
    #needed to correct because my interval doesn't end with 1 in the last position. Need to generalize this
    dll = (movey*init$C*exp(-f))[movey==1] - (init$y *init$C * exp(-f))[init$y==1]
    #dll = - (init$y *init$C * exp(-f))[init$y==1]
    }#dll = apply(init$rep_idx,1,function(idx)sum(-init$y[idx[1]:idx[2]]+llnocoal[idx[1]:idx[2]])) # gradient of log-likelihood wrt f_midpts
    
    return(dll)
  }
}

samp_loglik = function(init, fs, betas)
{
  fs = as.matrix(fs)
  
  if (init$ng != dim(fs)[1])
    stop(paste("Incorrect number of rows for fs; should be", init$ng))
  
  if (length(betas) != dim(fs)[2] + 1)
    stop(paste("Incompatible number of betas and columns for fs"))
  
  beta0 = betas[1]
  betas = utils::tail(betas, -1)
  
  #f = rep(f, init$gridrep)
  fs = as.matrix(fs[rep(1:init$ng, init$gridrep),])
  
  #llsampevents = beta1 * init$count * f
  #llsampnoevents = init$D * exp(beta0) * exp(f)^beta1
  #llsamp = init$ns * beta0 + sum(llsampevents[!is.na(init$count)]) - sum(llsampnoevents[!is.na(init$count)])
  
  fs_betas = fs %*% betas
  
  llsampevents = init$count * fs_betas
  llsampnoevents = init$D * exp(beta0 + fs_betas)
  llsamp = init$ns * beta0 + sum(llsampevents[!is.na(init$count)]) - sum(llsampnoevents[!is.na(init$count)])
  #print(init$ns * beta0)
  #print(sum(llsampevents[!is.na(init$count)]))
  #print(-1 * sum(llsampnoevents[!is.na(init$count)]))
  
  return(llsamp)
}

coal_samp_loglik = function(init, f, beta0, beta1)
{
  if (init$ng != length(f))
    stop(paste("Incorrect length for f; should be", init$ng))
  
  f = rep(f, init$gridrep)
  
  llnocoal = init$D * init$C * exp(-f)
  
  lls = -init$y * f - llnocoal
  
  llsampevents = beta1 * init$count * f
  llsampnoevents = init$D * exp(beta0) * exp(f)^beta1
  llsamp = init$ns * beta0 + sum(llsampevents[!is.na(init$count)]) - sum(llsampnoevents[!is.na(init$count)])
  
  llcoal = sum(lls[!is.nan(lls)])
  
  return(llcoal + llsamp)
}

coal_samp_fns_loglik = function(init, f, fs, beta0, beta1, betas)
{
  if (init$ng != length(f))
    stop(paste("Incorrect length for fs; should be", init$ng))
  
  if (init$ng != dim(fs)[1])
    stop(paste("Incorrect number of rows for fs; should be", init$ng))
  
  f = rep(f, init$gridrep)
  fs = fs[rep(1:init$ng, init$gridrep),]
  
  llnocoal = init$D * init$C * exp(-f)
  
  lls = -init$y * f - llnocoal
  
  llsampevents = beta1 * init$count * f + diag(init$count) %*% fs %*% betas
  llsampnoevents = init$D * exp(beta0 + f * beta1 + fs %*% betas)
  llsamp = init$ns * beta0 + sum(llsampevents[!is.na(init$count)]) - sum(llsampnoevents[!is.na(init$count)])
  
  llcoal = sum(lls[!is.nan(lls)])
  
  return(llcoal + llsamp)
}

U = function(theta, init, invC, alpha, beta, grad=FALSE)
{
  D = length(theta)
  f = theta[-D]
  tau = theta[D]
  invCf = invC %*% f
  if(!grad)
  {
    loglik = coal_loglik(init, f)
    logpri = ((D-1)/2+alpha)*tau - (t(f)%*%invCf/2+beta)*exp(tau)
    return(list(loglik = -loglik, logpri = -logpri, logpos = -(loglik+logpri)))
  }
  else
  {
    dloglik = c(coal_loglik(init, f, grad),0)
    dlogpri = c(-invCf*exp(tau),((D-1)/2+alpha)-(t(f)%*%invCf/2+beta)*exp(tau))
    return(list(dloglik = -dloglik, dlogpri = -dlogpri, dlogpos = -(dloglik+dlogpri)))
  }
}

#python.call("calcPF", result$change_F, result$F_nodes, result$family_size, oldsuff,u,"True")

U_times = function(times,nsites, result, oldsuff, theta, grid,logliktot,indicator,const=1,meantheta)
  #print(paste0("times:",sum(times)))
  #print(paste0("grid:",sum(grid)))

  #For HMC of coalescent times, times are u
  ##Dynamic grid construction
{
  newtheta<-theta
  if (max(grid)<sum(times)) {
    #I need to update theta and the grid
    oldsize<-length(grid)
    intl = grid[2]-grid[1]
    ngrid<-ceiling(sum(times)/intl)
    grid_bds = range(ngrid*intl, 0); 
    #print(sum(times))
    #print(intl)
    #print(ngrid)
    #print(grid_bds)
    grid = seq(grid_bds[1], grid_bds[2], by=intl)
    midpts = grid[-1]-intl/2
    if (length(grid)>length(meantheta)){
      meantheta<-c(meantheta,rep(meantheta[length(theta)-1],length(grid)-length(meantheta)))
    }
    if ((length(grid) <= length(meantheta)) & length(grid)>length(theta) ){
      newtheta<-c(theta[-length(theta)],meantheta[length(theta):(length(grid)-1)],theta[length(theta)])
    }
   # newtheta<-c(theta[-length(theta)],rep(theta[length(theta)-1],length(grid)-oldsize),theta[length(theta)])
    }
  init = coal_lik_init(samp_times =0, n_sampled = length(times)+1, 
                           coal_times = cumsum(times), grid = grid)
  if (indicator){
  logliktot<-python.call("calcPF", result$change_F, result$F_nodes, result$family_size, oldsuff,times*nsites,"True")
  }
  #loglik = logliktot[[2]]-sum(seq(length(times)+1,2)*times*nsites)
    loglik = logliktot[[2]]/const
    logpri = coal_loglik(init, newtheta[-length(newtheta)])
    dloglik=(logliktot[[3]])*nsites/const #gives the gradient wrt times
   # dloglik=(logliktot[[3]]-seq(length(times)+1,2))*nsites #gives the gradient wrt times
    dlogpri = coal_loglik2(init, newtheta[-length(newtheta)],TRUE) 
    return(list(logliktot=logliktot,dloglik = -dloglik, dlogpri = -dlogpri, dlogpos = -(dloglik+dlogpri),loglik = -loglik, logpri = -logpri, logpos = -(loglik+logpri),newtheta=newtheta,grid=grid,meantheta=meantheta))
 }
# 
# result<-resultL
# oldsuff<-initial$oldsuffL
# times<-initial$timesL
# 
# U_times2L = function(times,tmrca,nsites, result, oldsuff, theta, grid,logliktot,indicator,const=1,sig=.01)
#   #For HMC of coalescent times, times are u and fixing the rjMCMC 
#   ##Dynamic grid construction
# {
#   newtheta<-theta
#   oldsize<-length(grid)
#   #print(oldsize)
#   intl = grid[2]-grid[1]
#   ngrid<-ceiling(tmrca/intl)
#   grid_bds = range(ngrid*intl, 0); 
#   grid = seq(grid_bds[1], grid_bds[2], by=intl)
#   rjmcmc<-0
#   if (length(grid)>oldsize){
#     #print("")
#     #newval<-rnorm(length(grid)-oldsize,0,sig)
#     newval<-rep(0,length(grid)-oldsize)
#     #  -0.5*t(newval)%*%diag(1/sig^2,length(newval))%*%newval-3*log(sqrt(2*pi)*sig)
#     
#     # rjmcmc<--dmvnorm(newval,sigma=diag(sig^2,length(newval)),log=TRUE)
#     #rjmcmc<--dnorm(newval,mean=0,sd=sig,log=TRUE)
#     newtheta<-c(theta[-length(theta)],newval+theta[length(theta)-1],theta[length(theta)])
#   }
#   
#   if (length(grid)<oldsize){
#     newval<-theta[length(grid):(oldsize-1)]-theta[length(grid)-1]
#     # rjmcmc<-dmvnorm(newval,sigma=diag(sig^2,length(newval)),log=TRUE)
#     #rjmcmc<-dnorm(newval,mean=0,sd=sig,log=TRUE)
#     newtheta<-c(theta[1:(length(grid)-1)],theta[length(theta)])
#   }
#   
#   #print(length(grid))
#   lik_initL<-list()
#   loglik<-0
#   dloglik=rep(0,length(times[[1]]))
#   dlogpri=rep(0,length(times[[1]]))
#   for (j in 1:length(times)){
#   lik_initL[[j]] = coal_lik_init(samp_times =0, n_sampled = length(times[[j]])+1, 
#                        coal_times = cumsum(times[[j]]), grid = grid)
#   if (indicator){
#     logliktot[[j]]<-python.call("calcPF", result[[j]]$change_F, result[[j]]$F_nodes, result[[j]]$family_size, oldsuff[[j]],timesL[[j]]*nsites,"True")
#   }
#   
#   loglik=loglik+lik_call[[j]][[2]]/const
#   logpri=logpri+coal_loglik(lik_initL[[j]], newtheta[-length(newtheta)])
#   dloglik=dloglik+(lik_call[[j]][[3]])*nsites 
#   dlogpri=coal_loglik2(lik_initL[[j]], newtheta[-length(newtheta)],TRUE) 
#   }
#   # print(length(grid))
#   #if (indicator){
#   #  logliktot<-python.call("calcPF", result[[j]]$change_F, result$F_nodes, result$family_size, oldsuff[[j]],timesL[[j]]*nsites,"True")
#   #}
#   #loglik = logliktot[[2]]-sum(seq(length(times)+1,2)*times*nsites)
#   #loglik = logliktot[[2]]/const
#   #dloglik=(logliktot[[3]])*nsites #gives the gradient wrt times
#   # dloglik=(logliktot[[3]]-seq(length(times)+1,2))*nsites #gives the gradient wrt times
#   #dlogpri = coal_loglik2(init, newtheta[-length(newtheta)],TRUE) 
#   return(list(rjmcmc=rjmcmc,init=lik_initL,logliktot=logliktot,dloglik = -dloglik, dlogpri = -dlogpri, dlogpos = -(dloglik+dlogpri),loglik = -loglik, logpri = -logpri, logpos = -(loglik+logpri),newtheta=newtheta,grid=grid))
# }

U_times2 = function(times,nsites, result, oldsuff, theta, grid,logliktot,indicator,const=1,sig=.01)
  #For HMC of coalescent times, times are u and fixing the rjMCMC 
  ##Dynamic grid construction
{
  newtheta<-theta
  oldsize<-length(grid)
  #print(oldsize)
  intl = grid[2]-grid[1]
  ngrid<-ceiling(sum(times)/intl)
  grid_bds = range(ngrid*intl, 0); 
  grid = seq(grid_bds[1], grid_bds[2], by=intl)
  rjmcmc<-0
    if (length(grid)>oldsize){
      #print("")
      #newval<-rnorm(length(grid)-oldsize,0,sig)
      newval<-rep(0,length(grid)-oldsize)
      #  -0.5*t(newval)%*%diag(1/sig^2,length(newval))%*%newval-3*log(sqrt(2*pi)*sig)
      
     # rjmcmc<--dmvnorm(newval,sigma=diag(sig^2,length(newval)),log=TRUE)
      #rjmcmc<--dnorm(newval,mean=0,sd=sig,log=TRUE)
      newtheta<-c(theta[-length(theta)],newval+theta[length(theta)-1],theta[length(theta)])
    }
  
  if (length(grid)<oldsize){
    newval<-theta[length(grid):(oldsize-1)]-theta[length(grid)-1]
   # rjmcmc<-dmvnorm(newval,sigma=diag(sig^2,length(newval)),log=TRUE)
    #rjmcmc<-dnorm(newval,mean=0,sd=sig,log=TRUE)
    newtheta<-c(theta[1:(length(grid)-1)],theta[length(theta)])
  }
  
  #print(length(grid))
  init = coal_lik_init(samp_times =0, n_sampled = length(times)+1, 
                       coal_times = cumsum(times), grid = grid)
 # print(length(grid))
  if (indicator){
    logliktot<-python.call("calcPF", result$change_F, result$F_nodes, result$family_size, oldsuff,times*nsites,"True")
    }
  
  #loglik = logliktot[[2]]-sum(seq(length(times)+1,2)*times*nsites)
  loglik = logliktot[[2]]/const
  logpri = coal_loglik(init, newtheta[-length(newtheta)])
  dloglik=(logliktot[[3]])*nsites #gives the gradient wrt times
  # dloglik=(logliktot[[3]]-seq(length(times)+1,2))*nsites #gives the gradient wrt times
  dlogpri = coal_loglik2(init, newtheta[-length(newtheta)],TRUE) 
  return(list(rjmcmc=rjmcmc,init=init,logliktot=logliktot,dloglik = -dloglik, dlogpri = -dlogpri, dlogpos = -(dloglik+dlogpri),loglik = -loglik, logpri = -logpri, logpos = -(loglik+logpri),newtheta=newtheta,grid=grid))
}


log_f_prior_kappa = function(f, kappa, invC, alpha, beta)
{
  D = length(f) + 1
  invCf = invC %*% f
  
  return(((D-1)/2+alpha-1)*log(kappa) - (t(f)%*%invCf/2+beta)*kappa)
}



log_f_prior_kappa = function(f, kappa, invC, alpha, beta)
{
  D = length(f) + 1
  invCf = invC %*% f
  
  return(((D-1)/2+alpha-1)*log(kappa) - (t(f)%*%invCf/2+beta)*kappa)
}
U_kappa = function(theta, init, invC, alpha, beta, grad=FALSE)
{
  D = length(theta)
  f = theta[-D]
  kappa=theta[D]
  invCf = invC %*% f
  if(!grad)
  {
    loglik = coal_loglik(init, f)
    logpri = ((D-1)/2+alpha-1)*log(kappa) - (t(f)%*%invCf/2+beta)*kappa
    return(list(loglik = -loglik, logpri = -logpri, logpos = -(loglik+logpri)))
  }
  else
  {
    dloglik = c(coal_loglik(init, f, grad),0)
    dlogpri = c(-invCf*kappa,((D-1)/2+alpha-1)/kappa-(t(f)%*%invCf/2+beta))
    return(list(dloglik = -dloglik, dlogpri = -dlogpri, dlogpos = -(dloglik+dlogpri)))
  }
}

U_splitL = function(theta, init, invC, alpha, beta, grad=FALSE)
{
  D=length(theta)
  f=theta[-D]
  tau=theta[D]
  invCf=invC%*%f
  if(!grad)
  {
    loglik_list<-list()
    loglik = coal_loglik(init[[1]], f)
    loglik_list[[1]]<-loglik
    logpri = ((D-1)/2+alpha)*tau - (t(f)%*%invCf/2+beta)*exp(tau)
    for (j in 2:length(init)){
      loglik=loglik+coal_loglik(init[[j]], f)
      loglik_list[[j]]<-loglik
    }
    return(list(loglik = -loglik, logpri = -logpri, logpos = -(loglik+logpri), loglik_list=loglik_list))
  }
  else
    {
    #This gradient is wrt log effective pop size
    #loglik_grad_list<-list()
    totgrad<-coal_loglik(init[[1]], f, grad)
    #loglik_grad_list[[1]]<-totgrad
    for (j in 2:length(init)){
      totgrad<-totgrad+coal_loglik(init[[j]], f, grad)
      #loglik_grad_list[[j]]<-coal_loglik(init[[j]], f, grad)
      #print(loglik_grad_list[[j]])
    }
    dU_res = -c(totgrad,((D-1)/2+alpha)-beta*exp(tau))
    return(dU_res)	
  }
}

U_split = function(theta, init, invC, alpha, beta, grad=FALSE)
{
  D=length(theta)
  f=theta[-D]
  tau=theta[D]
  invCf=invC%*%f
  if(!grad)
  {
    loglik = coal_loglik(init, f)
    logpri = ((D-1)/2+alpha)*tau - (t(f)%*%invCf/2+beta)*exp(tau)
    return(list(loglik = -loglik, logpri = -logpri, logpos = -(loglik+logpri)))
  }
  else
  {
    dU_res = -c(coal_loglik(init, f, grad),((D-1)/2+alpha)-beta*exp(tau))
    return(dU_res)	
  }
}


precBM = function(times, delta=1e-6)
{
  D = length(times)
  
  diff1 <- diff(times)
  diff1[diff1==0] <- delta
  
  diff <- 1/diff1
  Q <- spam::spam(0,D,D)
  if (D>2)
    Q[cbind(1:D,1:D)] <- c(diff[1]+ifelse(times[1]==0,1/delta,1/times[1]),diff[1:(D-2)]+diff[2:(D-1)],diff[D-1])
  else
    Q[cbind(1:D,1:D)] <- c(diff[1]+ifelse(times[1]==0,1/delta,1/times[1]),diff[D-1])
  
  Q[cbind(1:(D-1),2:D)] = -diff[1:(D-1)]; Q[cbind(2:D,1:(D-1))]=-diff[1:(D-1)]
  return(Q)
}

Met = function(theta, init, invC)
{
  D = length(theta)
  f = theta[-D]
  kappa=theta[D]
  
  f = rep(f, init$gridrep)
  llnocoal = init$D * init$C * exp(-f)
  diagF = apply(init$rep_idx,1,function(idx)sum(llnocoal[idx[1]:idx[2]]))
  
  G = invC*kappa + diag(diagF)
  return(G)
}
