###Functions for Bounded Coalescent 
###Piecewise constant prior and  MLE
###Authors: Mackenzie Simper and Julia Palacios
###Oct 2024


#' Compute the MLE of effective population size
#' under the bounded coalescent
#' 
#' @param data A list of coal_times, samp_times and n_sampled
#'        bound Value of upper bound on TMRCA
#' @return grid and Ne estimates
#' @export
MLE_BoundCoal<-function(data,lengthout=5,bound,eps=.02,eta=.01){
  if (class(data) == "phylo") {
    phy <- summarize_phylo(data)
  }
  else if (all(c("coal_times", "samp_times", "n_sampled") %in% 
               names(data))) {
    phy <- with(data, list(samp_times = samp_times, coal_times = coal_times, 
                           n_sampled = n_sampled))
  }  
  grid_bds<-range(0, bound) #range(c(coal_times,samp_times))
  grid = seq(grid_bds[1],grid_bds[2], length.out = lengthout)
  lik_init = phylodyn:::coal_lik_init(samp_times = c(0), n_sampled = n, 
                                      coal_times = coal_times1, grid=grid)
  f_init = rep(0, lik_init$ng)
  eta = 0.01
  eps = 0.02
  gradResult <- Ne_gradient_ascent(f_init, lik_init, bound, eps, eta)
  return(list(grid=grid,x=grid,effpop=exp(c(gradResult[[1]][1],gradResult[[1]]))))
}


#computes all values of probability2
bound_prob_base2<-function(fgrid,deltagrid,n){
  lineages<-seq(n,2)
  ngrid<-length(fgrid)
  precal2<-list()
  for (j in 1:(n-2)){
    #j=1, indicates first coalecent (t_{n})
    precal2[[j]]<-list()
    for (i in 1:(n-j-1)){#indicates the length of equal indices
      precal2[[j]][[i]]<-list()
      for (k in (j+i):(n-1)){ #indicates, extra coal. integrating out
        
        temp<-bound_prob2_2(fgrid,deltagrid,lineages[j:(j+i-1)],lineages[k])
        precal2[[j]][[i]][[k]]<-matrix(0,nrow=ngrid,ncol=2)
        precal2[[j]][[i]][[k]][,1]<-temp$prob
        precal2[[j]][[i]][[k]][,2]<-temp$gradient
        
      }
    }
  }
  return(precal2)
}

#vector version
bound_prob2_2<-function(fgrid,deltagrid,lineages,l2){
  B<-length(lineages)
  if (B==1){
    fact1<-choose(lineages[1],2)/(choose(lineages[1],2)-choose(l2,2))
    prob<-fact1*(exp(-choose(l2,2)*deltagrid*exp(-fgrid))-exp(-choose(lineages[1],2)*deltagrid*exp(-fgrid)))
    gradient<-fact1*(exp(-choose(l2,2)*deltagrid*exp(-fgrid))*choose(l2,2)*deltagrid*exp(-fgrid)-exp(-choose(lineages[1],2)*deltagrid*exp(-fgrid))*choose(lineages[1],2)*deltagrid*exp(-fgrid))
    return(list(prob=prob,gradient=gradient))
    
  }else{
    fact<-choose(lineages[B],2)/(choose(lineages[B],2)-choose(l2,2))
    temp1<-bound_prob2_2(fgrid,deltagrid,lineages[1:(B-1)],l2)
    temp2<-bound_prob2_2(fgrid,deltagrid,lineages[1:(B-1)],lineages[B])
    prob<-fact*(temp1$prob-temp2$prob)
    gradient<-fact*(temp1$gradient-temp2$gradient)
    return(list(prob=prob,gradient=gradient))
  }
}


bound_prob1_2<-function(fgrid,deltagrid,lineages){
  B<-length(lineages)
  if (B==1){prob<-(1-exp(-choose(lineages,2)*deltagrid*exp(-fgrid)))
  gradient<--exp(-choose(lineages,2)*deltagrid*exp(-fgrid))*choose(lineages,2)*deltagrid*exp(-fgrid)
  }else{
    temp1<-bound_prob1_2(fgrid,deltagrid,lineages[1:(B-1)])
    temp2<-bound_prob2_2(fgrid,deltagrid,lineages[1:(B-1)],lineages[B])
    ##Calculation of P_{1,1,1,1..,1}
    prob<-temp1$prob-temp2$prob
    gradient<-temp1$gradient-temp2$gradient
  }
  return(list(prob=prob,gradient=gradient))
}

bound_prob_base<-function(fgrid,deltagrid,n){
  lineages<-seq(n,2)
  ngrid<-length(fgrid)
  precal<-list()
  for (j in 1:(n-1)){
    #j=1, indicates first coalecent (t_{n})
    precal[[j]]<-list()
    for (i in 1:(n-j)){
      temp<-bound_prob1_2(fgrid,deltagrid,lineages[j:(j+i-1)])
      precal[[j]][[i]]<-matrix(0,nrow=ngrid,ncol=2)
      precal[[j]][[i]][,1]<-temp$prob
      precal[[j]][[i]][,2]<-temp$gradient
    }
  }
  return(precal)
}


bound_prob0<-function(fgrid,deltagrid,n){
  ngrid<-length(deltagrid)
  lineages<-seq(n,2)
  precal<-bound_prob_base(fgrid,deltagrid,n)
  precal2<-bound_prob_base2(fgrid,deltagrid,n)
  totprob<-0
  
  totgradient<-rep(0,length(fgrid))
  state<-rep(1,n-1)
  #print(state)
  totprob<-totprob+precal[[1]][[n-1]][1,1]
  totgradient[1]<-totgradient[1]+precal[[1]][[n-1]][1,2]
  
  #state denotes the subindex of P
  state_1<-as.numeric(paste(state,collapse=""))
  statemax<-rep(ngrid,n-1)
  # statemaxnum<-as.numeric(paste(statemax,collapse=""))+0
  while (sum(state==statemax)<(n-1)){
    #state_1<-state_1+1
    #state<-as.numeric(strsplit(as.character(state_1), "")[[1]])
    state[n-1]<-state[n-1]+1
    maxstate<-max(state)
    while  (maxstate>ngrid){
      n_ind<-min(seq(1,length(state))[state==maxstate])-1
      state[n_ind:(n-1)]<-state[n_ind]+1
      #state_1<-as.numeric(paste(state,collapse=""))+0
      
      maxstate<-max(state)
    }
    #print(state)
    
    if (sum(diff(state)>=0)>0){
      #print(state)
      maxstate<-max(state)
      changepoints<-seq(2,n-1)[diff(state)>=1] #position of changepoint
      if (sum(diff(state)>=1)==0){
        #print("This happens")
        #no changepoints, for example 2222
        factor1<-1
        factor1grad<-0
        if (state[1]>1){ #do not coalesce all the way
          factor1<-exp(-choose(n,2)*sum(deltagrid[1:(state[1]-1)]*exp(-fgrid[1:(state[1]-1)])))
          factor1vec<-exp(-choose(n,2)*deltagrid[1:(state[1]-1)]*exp(-fgrid[1:(state[1]-1)]))
          factor1grad<-factor1vec*choose(n,2)*deltagrid[1:(state[1]-1)]*exp(-fgrid[1:(state[1]-1)])
        }
        #f<-rep(fgrid[state[1]],n-1)
        #delta<-rep(deltagrid[state[1]],n-1)
        ###Here
        #temp<-bound_prob1(f,delta,lineages)
        factor<-factor1*precal[[1]][[n-1]][state[1],1]
        if (state[1]>1 & max(state)<=ngrid){
          for (k in 1:(state[1]-1)){
            totgradient[k]<-totgradient[k]+precal[[1]][[n-1]][state[1],1]*prod(factor1vec[-k])*factor1grad[k]
          }
          # print(totgradient)
        }
        
        #totgradient[1]<-totgradient[1]+factor1*precal[[1]][[n-1]][state[1],2]
        #print(totgradient)
        
      }else{
        #there are  changepoints
        factor1<-1
        factor1grad<-0
        factors<-c(1)
        factorsgrad<-c(0)
        indices<-c(0)##indicates states
        if (state[1]>1){
          factor1<-exp(-choose(n,2)*sum(deltagrid[1:(state[1]-1)]*exp(-fgrid[1:(state[1]-1)])))
          factor1vec<-exp(-choose(n,2)*deltagrid[1:(state[1]-1)]*exp(-fgrid[1:(state[1]-1)]))
          factor1grad<-factor1vec*choose(n,2)*deltagrid[1:(state[1]-1)]*exp(-fgrid[1:(state[1]-1)])
          factors<-c(factors,factor1vec)
          factorsgrad<-c(factorsgrad,factor1grad)
          indices<-c(indices,seq(1,(state[1]-1)))
        } 
        #f<-rep(fgrid[state[1]],changepoints[1]-1)
        #delta<-rep(deltagrid[state[1]],(changepoints[1]-1))
        #lineages1<-lineages[1:(changepoints[1]-1)]
        #temp<-bound_prob2(f,delta,lineages1,lineages[changepoints[1]])
        
        factor<-factor1*precal2[[1]][[(changepoints[1]-1)]][[changepoints[1]]][state[1],1]
        factors<-c(factors,precal2[[1]][[(changepoints[1]-1)]][[changepoints[1]]][state[1],1])
        factorsgrad<-c(factorsgrad,precal2[[1]][[(changepoints[1]-1)]][[changepoints[1]]][state[1],2])
        indices<-c(indices,state[1])
        for (i in 1:length(changepoints)){
          #now to the right of changepoints and left if jumps
          index1<-changepoints[i]
          #if there is a jump greater than one between changepoints
          factor2<-1
          factor2vec<-1
          factor2grad<-0
          if ((state[index1]-state[(index1-1)])>1){ #need to multiply by no coalescing there
            factor2<-exp(-choose(lineages[index1],2)*sum(deltagrid[(state[(index1-1)]+1):(state[index1]-1)]*exp(-fgrid[(state[(index1-1)]+1):(state[index1]-1)])))
            #a vector
            factor2vec<-exp(-choose(lineages[index1],2)*deltagrid[(state[(index1-1)]+1):(state[index1]-1)]*exp(-fgrid[(state[(index1-1)]+1):(state[index1]-1)]))
            factor2grad<-factor2vec*choose(lineages[index1],2)*deltagrid[(state[(index1-1)]+1):(state[index1]-1)]*exp(-fgrid[(state[(index1-1)]+1):(state[index1]-1)])
            factors<-c(factors,factor2vec)
            factorsgrad<-c(factorsgrad,factor2grad)
            indices<-c(indices,(state[(index1-1)]+1):(state[index1]-1))##need to figure out the index
          }
          
          if (i==length(changepoints)){index2<-n-1}else{index2<-changepoints[(i+1)]}
          
          if (i==length(changepoints)){# last change point
            #temp<-bound_prob1(f,delta,lineages1)
            factor<-factor*factor2*precal[[index1]][[index2-index1+1]][state[index1],1]
            indices<-c(indices,state[index1])
            factors<-c(factors,precal[[index1]][[index2-index1+1]][state[index1],1])
            factorsgrad<-c(factorsgrad,precal[[index1]][[index2-index1+1]][state[index1],2])
          }else{
            #f<-rep(fgrid[state[index1]],index2-index1)
            #delta<-rep(deltagrid[state[index1]],index2-index1)
            #lineages1<-lineages[index1:(index2-1)]
            #need fix here
            #temp<-bound_prob2(f,delta,lineages1,lineages[index2])
            #bound_prob2_2(fgrid,deltagrid,4,lineages[index2])$prob[state[index1]]
            
            factor<-factor*factor2*precal2[[index1]][[(index2-index1)]][[index2]][state[index1],1]
            indices<-c(indices,state[index1])
            factors<-c(factors,precal2[[index1]][[(index2-index1)]][[index2]][state[index1],1])
            factorsgrad<-c(factorsgrad,precal2[[index1]][[(index2-index1)]][[index2]][state[index1],2])}
        }
        indices<-indices[-1]
        factorsgrad<-factorsgrad[-1]
        factors<-factors[-1]
        for (k in 1:length(indices)){
          totgradient[indices[k]]<-totgradient[indices[k]]+prod(factors[-k])*factorsgrad[k]
          #  print(totgradient)
        }
      }
      
      totprob=totprob+factor
    }
  }
  return(list(totprob=totprob,totgradient=totgradient/totprob))
}


coal_loglik_bounded = function(init, f, bound)
{
  if (init$ng != length(f))
    stop(paste("Incorrect length for f; should be", init$ng))
  fext =f
  f = rep(f, init$gridrep)
  #print(fext)
  Dext = diff(init$args$grid)
  Dext[length(Dext)]<-Dext[length(Dext)]+(bound-max(init$args$grid))
  #Dext = c(init$D,(bound-max(init$args$grid)))
  #print(Dext)
  #gridext=c(init$args$grid,bound)
  llnocoal = init$D * init$C * exp(-f)
 # print(Dext)
#  print(fext)
#  print(init$ns)
  bound_probability<-bound_prob0(fext,Dext,init$ns)
  
  lls = -init$y * f - llnocoal -log(bound_probability$totprob)
  
  ll = sum(lls[!is.nan(lls)])
  dll = apply(init$rep_idx,1,function(idx)sum(-init$y[idx[1]:idx[2]]+llnocoal[idx[1]:idx[2]]))-bound_probability$totgradient # gradient of log-likelihood wrt f_midpts
  return(list(ll=ll,dll=dll))
}

Ne_gradient_ascent <- function(f_init, lik_init, bound, eps, eta) {
  diff = 1
  numSteps = 1
  
  currF = f_init
  ll = c()
  while (diff > eps) {
    
    result = coal_loglik_bounded(lik_init, currF, bound)
    ll = c(ll, result$ll)
    newF = currF + eta*result$dll
    
    diff = sum(abs(newF - currF))
    if (numSteps %% 10 == 0) { #Print out results every 10 steps
      print(newF)
      print(diff)
    }
    numSteps = numSteps + 1
    print(numSteps)
    currF = newF
    print(newF)
    print(diff)
    print(result$ll)
  }
  
  return(list(newF, ll))
}