####Functions needed for simulating data for lineage tracing inference methods
#author: Mackenzie Simper, Julia Palacios
#library(hash)
#library(phylodyn)
#library(stringr)

#' Simulation of alleles at S target sites evolving along the tree with M types
#' @description This function simulates alleles for all nodes according to a CTMC
#' 
#' @param numSim number of simulated datasets
#' @param tree genealogical tree in ape format 
#' @param S number of target sites
#' @param M number of different alleles
#' @param theta a matrix of mutation rates
#' It returns a list where the first
#' entry is the matrix of mutation states for the nodes, second entry is matrix of allele states. 
simulate_SLData<-function(numSim=1,tree,S=2,M=4,theta=diag(1,S)){
  M = matrix(M, nrow = S, ncol = S) #allele types
  states = state_space(S)
  states_matrix = state_space_matrix(S) #just because I don't want to change everything over to hash yet
  Q = Q_matrixSL(states_matrix, theta)
  eigenQ<-eigen(Q)
  eigenQ$inv<-solve(eigenQ$vectors)
  
  ##vectors to store output
  theta_vec <- rep(1,(S*(S+1))*.5)
  theta_vec_list_k <- list()
  naive_theta_list <- list()
  datalist <- list()
  Dlist <- list()
  theta_vec_list <- matrix(0, nrow = numSim, ncol=length(theta_vec))  
  for (i in 1:numSim){
    
    datalist[[i]]= simulate_data_finite_alleles(tree, Q, states_matrix, M)
  }
  return(datalist)
}

simulate_data_finite_alleles <- function(tree, Q, states, M) {
  n = tree$Nnode + 1
  S = ncol(states)
  alleles_lumped = matrix(0, 2*n-1, S) #To store alleles of all nodes
  alleles_full = matrix(0, 2*n-1, S)
  
  #The leaves are labeled 1, ..., n and the interior nodes are
  #labeled n+1, ..2n - 1, with n+1 as the root
  edges <- tree$edge #Returns a list of (parent, child)
  edge_lengths <- tree$edge.length #2n - 2 edges
  #loop through all edges
  for (i in 1:(2*n - 2)) {
    parent = edges[i, 1]
    child = edges[i, 2]
    t = edge_lengths[i]
    
    #Assuming that we're going down in the tree,
    #the parent allele should already be set
    parent_lumped = alleles_lumped[parent, ]
    child_lumped = transition2(parent_lumped, t, Q, states)
    alleles_lumped[child, ] = child_lumped
    
    #Now fill in the specific mutations in a consistent way
    parent_full = alleles_full[parent, ]
    child_full = rep(0, S)
    j = 1
    while (j < S + 1) {
      if (child_lumped[j] > 0) {
        if (child_lumped[j] == parent_lumped[j]) {
          child_full[j] = parent_full[j]
          j = j + 1
        }
        else { #A new  mutation
          if (child_lumped[j] == 1) { #a new mutation at a single site
            m = sample(1:M[j, j], 1) #sample uniformly from the number of possible mutations
            
            #Just want a way to indicate this is a single mutation at site j
            #One way of doing that is to store as a string j:m
            child_full[j] = paste(paste(j, ":", sep = ""), m, sep = "")
            j = j + 1
          }
          else { #an overlapping mutation
            all = which(child_lumped == child_lumped[j]) #To assign all sites in the overlap the same mutation
            i = max(all) #assume the mutation arose fro a simultaneous cut at sites i and j
            m = sample(1:M[i, j], 1)
            
            #The overlap mutation will be stored ji:m (where j is start site, i is end site)
            child_full[all] = paste(paste(j*10 + i, ":", sep = ""), m, sep = "")
            j = j + length(all)
          }
        }
      }
      else {j = j +1}
    }
    alleles_full[child, ] = child_full
  }
  
  return(list(alleles_lumped, alleles_full))
}

#Generates the rate matrix for the list of states
#theta is an SxS (S is numer of states) matrix, with theta[i, j] the rate of simultaneous
#cut at site i and j (could be theta[i, i] + theta[j, j] or more general), and theta[i, i] the rate of single cut at site i
Q_matrixSL <- function(states, theta) {
  L = nrow(states)
  S = ncol(states)
  M = matrix(0, nrow = L, ncol = L)
  for (i in 1:L) {
    s1 = states[i, ]
    active = which(s1 == 0)
    if (length(active) > 0) {
      #First, all possible single mutations
      for (k in 1:length(active)) {
        site = active[k]
        s2 = s1
        s2[site] = 1
        #This tells me which index the new state corresponds to,
        #so I can update the appropriate entry in the rate matrix
        #(Like a cheat version of a dictionary)
        j = which(apply(states, 1, function(x) all.equal(x[1:S], s2)) == "TRUE")
        M[i, j] = theta[site, site]  
        
        if (length(active) - k > 0) {
          #now all pairs of mutations
          for (l in (k+1):length(active)) {
            site2 = active[l]
            s2 = s1
            #This mutation assignment is consistent with what we did before
            s2[site:site2] = site + 1 #we should always have site2 > site
            j = which(apply(states, 1, function(x) all.equal(x[1:S], s2)) == "TRUE")
            M[i, j] = theta[site, site2]
          }
        }
      }
    }
  }
  
  #Add negatives along diagonal so rows sum to 1
  rows = rowSums(M)
  diagonal = diag(rows, nrow = nrow(M), ncol = ncol(M))
  M = M - diagonal
  return(M)
}


#turns the list into a matrix, and relables
#returns a matrix
state_space_matrix <- function(S) {
  states = state_space_list(S)
  return(relabel(states))
}

#Returns a hash table/dictionary
#states are stored as strings, and returns the index
state_space <- function(S) {
  states = state_space_list(S)
  states = relabel(states) #matrix
  h = hash()
  for (i in 1:nrow(states)) {
    state = paste(states[i, ], collapse = "")
    h[[state]] = i
  }
  
  return(h)
}

#Returns a list of the lumped state-space: the possible
#vectors of length S, where 0 indicates no mutation,
#1 indicates mutation just at that site, and larger numbers
#are mutations that are shared between multiple adjacent sites
#S is the number of sites
state_space_list <- function(S) {
  #recursive function
  if (S == 1) {
    return(list(c(1), c(0)))
  }
  else {
    SS_1 = state_space_list(S-1)
    
    new = list()
    for (state in SS_1) {
      new[[length(new) + 1]] = c(1, state)
      new[[length(new) + 1]] = c(0, state)
      
      if (state[1] == 1) { #can merge
        m = max(state) + 1
        new[[length(new) + 1]] = c(m, m, state[-1])
      }
      if (state[1] > 1) { #1 is already part of an overlap, merge with that
        new[[length(new) + 1]] = c(state[1], state)
      }
      
    }
    return(new)
  }
}

#To keep consistent, rename the mutations in an overlap
#I will name the mutation (start + 1), where start is the
#leftmost index.
#This takes as output the list from the previous function
#and returns a matrix
relabel <- function(states) {
  L = length(states) #states is a list
  new = matrix(0, nrow = L, ncol = length(states[[1]]))
  
  for (i in 1:L) {
    s = states[[i]]
    newS = s
    if (max(s) > 1) {
      for (k in 2:max(s)) {
        overlap = which(s == k)
        first = overlap[1]
        newS[overlap] <- first + 1
      }
    }
    new[i, ]  = newS
  }
  return(new)
}

#Returns transition in time t for given parent and rate matrix Q
#Also need matrix of states for indexing
library(expm)
transition2 <- function(parent, t, Q, states) {
  S = ncol(states)
  #print(parent)
  i = which(apply(states, 1, function(x) all.equal(x[1:S], parent)) == "TRUE")
  
  ##check that this is faster:eigenQ$vectors%*%diag(exp(eigenQ$values*t))%*%eigenQ$inv
  T_prob = expm(t*Q) ##This doesn't need to be called every single time
  j = sample(1:nrow(states), 1, replace = FALSE, prob = T_prob[i, ])
  
  return(states[j, ])
}


