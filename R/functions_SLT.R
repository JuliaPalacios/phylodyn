##All functions for reconstructing cell-lineage trees
#author: Mackenzie Simper, Julia Palacios

###################### Functions for estimation of theta
#' Estimation of mutation rates via Method-of-Moments
#' @description This function estimates the rates of the CTMC
#' 
#' @param data rows correspond to cells and columns to S target sites
#' @param t time to the most recent common ancestor
#' It returns a vector where the first S values are individual target mutations, followed by 
#' S(S-1)/2 off-diagonal(by row)
estimate_mu<-function(data,t){
  n<-nrow(data)
  S<-ncol(data)
  states = state_space(S)
  lumpedD<-t(apply(data,1,lumped_state))
  A<- observed_mutations2(data,S, lumpedD)$obsM
  A<-A
  B <- A/n
  #Need to store theta as a vector: The first S entries are the single-cut rates,
  #The next entries are the entries of the upper-triangular entries in matrix theta
  #theta0 =  c(1, 1, 1, 2, 2, 2) 
  ##Initial value of theta for optimization
  theta0<-rep(1,(S*(S+1)/2))*.01
  states_matrix<-state_space_matrix(S)
  #probs_function computes the sum of (observed - expected)^2 for each mutation type
  est_theta1 = optim(theta0, fn = probs_function,  observedMuts = B, t = t, S = S, 
                     states = states, states_matrix = states_matrix, lower = rep(0, length(theta0)), method = "L-BFGS-B")
  
  theta<-matrix(0,nrow=S,ncol=S)
  diag(theta)<-est_theta1$par[1:S]
  theta[lower.tri(theta,diag=FALSE)]<-est_theta1$par[-c(1:S)]
  theta=t(theta)
  
  return(theta)
  
}

#' Reformats data from 10x to be used with our method
#' @param data a table with the following columns: cellBC, intBC, lineageGrp, Allele,r1,r2,..
#' @param S number of target sites 
#' Although the column Allele is not used, the function assumes target information starts in column 5. 
#' 
#' @returns Dlist a list of cells x targets of alleles. One list per intBC 
reformat_slt<-function(data,S){
  
  intBCS = unique(data$intBC) #For datsets with multiple intBCs
  #Make a dictionary for cellBCs, to keep track across intBCs
  h = hash()
  cellBCS = unique(data$cellBC)
  for (i in 1:length(cellBCS)) {
    h[[cellBCS[i]]] = i
  }
  n = length(cellBCS)
  
  
  Dlist = list()
  
  for (j in 1:length(intBCS)) {
    #filter to one intBC
    datanew = data[which(data$intBC == intBCS[j]), ]
    #datanew = datanew[, 3:6] #14:17]
    rownames(datanew) <- NULL 
    # n = nrow(datanew)
    
    #----GET D for this intBC------------
    
    #the rows correspond to the cellBC indices stored in the hash table
    #if a cellBC does not have this intBC, just -1 in that row
    D = matrix(-1, nrow = n, ncol = S)
    #Keep track of how many unique alleles we see of each type
    numAlleles = matrix(0, nrow = S, ncol = S)
    
    uniqueA = list()
    for (s in 1:S) {
      edits = c()
      for (k in 1:nrow(datanew)) { #urg, must be a better way
        if (!grepl("None", datanew[k, 4 + s], fixed = TRUE)) { #WARNING: change index of datanew
          edits = c(edits,  unlist(datanew[k, 4 + s], use.names = FALSE))
        }
      }
      # temp = datanew[which( !grepl( "None", c(datanew[, 3 + s]), fixed = TRUE)), 3 + s] 
      uniqueA[[s]] <- unique(edits)
    }
    
    #loop through each cell
    for (i in 1:nrow(datanew)) {
      c = h[[datanew$cellBC[i]]] #the index of the cell
      allele = unlist(datanew[i, 5:(5+S-1) ], use.names = FALSE) #WARNING: CHANGE THIS depending on ALLELE COLUMNS; 4:7 for daisy data, 4:6 for Cassio data
      #Use 4:7 for filtered Nick data, 4:6 for cassio data
      if (any(is.na(allele))) { break } #skip if the cell happens to have na?
      s = 1
      while (s <= S) {
        if (!grepl( "None", allele[s], fixed = TRUE))  { #there is an edit
          #check if this is an overlap, or a single site
          temp = which(allele == allele[s])
          if (length(temp) > 1 ) { #it was an overlap
            r = max(temp)
            
            #Note: this formatting only works if we assume number of sites < 10
            D[c, s:r] = rep(paste(c(s*10 + r, which(uniqueA[[s]] == allele[s])), collapse = ":"), length(temp))
            s = r + 1
          }
          else {
            D[c, s] = paste(c(s, which(uniqueA[[s]] == allele[s])), collapse = ":")
            s = s + 1
          }
        }
        else { 
          D[c, s] = "0" #indicating "None" mutation
          s = s + 1 }
      }
    }
    
    Dlist[[j]] <- D
  }
  return(Dlist)
}




##counts the number of observed mutations from the matrix D, which stores alleles of data on tips
#(Note: we don't care about the actual alleles here, just the mutation state)
#Returns an SxS matrix where S[i,i] is number of cells with a single mutation in site i,
#S[i, j], i < j, is number of cells with an overlapping mutation going from i to j
observed_mutations2 <- function(D, S, lumpedD = NA) {
  obsM <- matrix(0, nrow=S, ncol=S)
  obsA <- matrix(0, nrow=S, ncol=S)
  n = nrow(D)
  
  if (max(is.na(lumpedD)) > 0) { #check if lumped states provided, otherwise have to get that now
    lumpedD = matrix(0, nrow = n, ncol = S)
    for (r in 1:n) {
      lumpedD[r, ] = lumped_state(D[r, ])
    }
  }
  
  for (i in 1:S) {
    haveMut = which(lumpedD[, i] == 1)
    obsM[i, i] = length(haveMut)
    obsA[i,i]<-length(unique(D[haveMut,i]))
  }
  for (i in 1:(S-1)) {
    for (j in (i+1):S) {
      haveMut = which(lumpedD[, i] == (i+1) & lumpedD[, j] == (i+1))
      if (j < S) {
        haveMut = intersect(haveMut, which(lumpedD[, j+1] != i+1))
        #this identifies the rows which have an overlap exactly between i, j (not extending any further)
      }
      #haveMut = which(lumpedD[, i:j] == rep(i+1, (j - i + 1)))
      obsM[i, j] = length(haveMut)
      obsA[i,j]= length(unique(D[haveMut,i]))
    }
  }
  vecA<-0
  for (j in 1:S){
    
    vecA<-c(vecA,obsA[cbind(1:(S-j+1),j:S)])
  }
  
  return(list(obsM=obsM,obsA=vecA[-1]))
}

#given rate vector theta, and time t, returns the probability of observing a mutation
#of a given type in time t
#Returns an SxS matrix, with (i, i) entry probability of single cut at site i,
#(i, j), i < j, probability of overlapping cut at sites i, j
#theta is a vector of length S + S(S-1)/2 -- the first S entries are the diagonal
#entries of the theta_metrix, the next S(S-1)/2 are the off-diagonal(by row)
probs_function <- function(theta, observedMuts, t, S, states, states_matrix) {
  
  theta_matrix = matrix(0, nrow = S, ncol = S)
  pos<-c(0,0)
  for (j in 1:S){
    pos<-rbind(pos,cbind(1:(S-j+1),j:S))
  }
  
  theta_matrix[pos[-1,]]<-theta
  #theta_matrix[lower.tri(theta_matrix, diag=FALSE)] <-  theta[(S+1):length(theta)]
  #tail(theta, as.integer((S-1)*(S)/2)) #theta[-(1:(S))] #get the last entries, starting from S+1
  #print(theta_matrix[lower.tri(theta_matrix, diag=FALSE)])
  #print(theta[-(1:(S+1))])
  #theta_matrix = t(theta_matrix)
  #theta_matrix = theta_matrix + diag(theta[1:S])
  
  
  Q  <- Q_matrixSL(states_matrix, theta_matrix)
  # print(Q)
  #eigenQ <- eigen(Q)
  #eigenQ$inv <- solve(eigenQ$vectors)
  #P=eigenQ$vectors%*%diag(exp(eigenQ$values*t))%*%eigenQ$inv
  P = expm(t*Q) #could try eigenvectors to speed up? But still have to compute eigenvectors for each new theta..
  
  probs = matrix(0, nrow = S, ncol = S)
  
  #Should always be the last entry in the states list, but just to check
  zeroStateIndex = states[[paste(rep(0, S), collapse = "")]]
  #which(apply(states, 1, function(x) all.equal(x[1:S], rep(0, S))) == "TRUE")
  for (i in 1:S) {
    #Find the indices of the states which have 1 in ith position
    indices = which(states_matrix[, i] == 1)
    probs[i, i] = sum(P[zeroStateIndex, indices])
  }
  
  for (i in 1:(S-1)) {
    for (j in (i+1):S) {
      #Based on how we assign mutation states, these are the states with an overlap
      #going from i to j (exactly)
      indices = which(states_matrix[, i] == (i+1) & states_matrix[, j] == (i+1))
      probs[i, j] = sum(P[zeroStateIndex, indices])
    }
  }
  #adding the probability of staying at (0,...0)
  #probs[S,1]<-P[zeroStateIndex,zeroStateIndex]
  #print(probs)
  diff = probs - observedMuts
  return(sum(diff^2))
}


#From the full states, removes mutation labels to give lumped state
lumped_state <- function(state) {
  S = length(state)                               
  lumped = rep(0, S)
  
  state[state == "None"| state == "0"] <- 0       
  
  if (S == 0) {
    print("lumped state error")
    print(state)
  }
  
  i = 1
  while (i < S + 1) {      
    if (state[i] != 0) {                          
      if (i < S && state[i] == state[i+1]) {                                
        all = which(state == state[i])
        lumped[all] <- i + 1
        i = i + length(all)
      } else{
        lumped[i] = 1                           
        i = i + 1
      }
    } else{
      i = i + 1
    }
  }
  return(lumped)
}




#######################
##Functions for estimating coalescent times

##counts the number of observed alleles from the matrix D, which stores alleles of data on tips
#So here the number of unique mutations observed
#Returns an SxS matrix where S[i,i] is number of cells with a single mutation in site i,
#S[i, j], i < j, is number of cells with an overlapping mutation going from i to j
observed_alleles <- function(D, lumpedD = NA) {
  S = ncol(D)
  obsM <- matrix(0, nrow=S, ncol=S)
  n = nrow(D)
  
  if (max(is.na(lumpedD)) > 0) { #check if lumped states provided, otherwise have to get that now
    lumpedD = matrix(0, nrow = n, ncol = S)
    for (r in 1:n) {
      lumpedD[r, ] = lumped_state(D[r, ])
    }
  }
  
  for (i in 1:S) {
    haveMut = which(lumpedD[, i] == 1)
    obsM[i, i] = length(unique(D[haveMut, i]))
  }
  for (i in 1:(S-1)) {
    for (j in (i+1):S) {
      haveMut = which(lumpedD[, i] == (i+1) & lumpedD[, j] == (i+1))
      if (j < S) {
        haveMut = intersect(haveMut, which(lumpedD[, j+1] != i+1))
        #this identifies the rows which have an overlap exactly between i, j (not extending any further)
      }
      #haveMut = which(lumpedD[, i:j] == rep(i+1, (j - i + 1)))
      obsM[i, j] = length(unique(D[haveMut, i]))
      obsM[j, i] = obsM[i, j]
    }
  }
  
  return(obsM)
}

#returns vector of gradients of log-likelihood,
#with respect to internal tree lengths
#p, q are matrices where each column corresponds
#to interior node; leaves are removed, so first
#column corresponds to n+1, n+2, ... 2n - 1
gradient <- function(p, q, Qnew) {
  like = t(q[, 1])%*%p[, 1]
  C = ncol(p) - 1 #exclude the last column which is the root, and doesn't have a branch 
  G = rep(0, C)
  for (c in 1:C) {
    p1 = p[, c]
    q1 = q[, c]
    
    
    G[c] = t(q1)%*%Qnew%*%p1
  }
  
  #G[i] is the gradient with respect to the branch which ends at node i  
  return(G/c(like))
}


#From a matrix of full allele states, builds the
#rate matrix specifically for those alleles
#So note that this will not be a complete rate matrix (the rows won't sum to 0),
#it just stores the rates specifically for the alleles
#which have non-zero probaiblity in the likelihood
Q_alleles <- function(alleles, Q, eigenQ, states, M, T_prob) {
  m = nrow(alleles)
  Qnew = matrix(0, nrow = m, ncol = m)
  
  for (i in 1:m) {
    allele1 = alleles[i, ]
    for (j in 1:m) {
      allele2 = alleles[j, ]
      Pt1= T_prob
      # Pt1 = eigenQ$vectors%*%diag(exp(eigenQ$values*1))%*%eigenQ$inv
      if (transition_prob_finite_alleles(allele1, allele2, Pt1, states, M) > 0) {
        lumped1 = lumped_state(allele1)
        lumped2 = lumped_state(allele2)
        I1 = states[[paste(lumped1, collapse = "")]]
        I2 = states[[paste(lumped2, collapse = "")]]
        # I1 = which(apply(states, 1, function(x) all.equal(x[1:S], lumped1)) == "TRUE")
        # I2 = which(apply(states, 1, function(x) all.equal(x[1:S], lumped2)) == "TRUE")
        Qnew[i, j] = Q[I1, I2]
      }
    }
  }
  return(Qnew)
}


#Returns the vectors which are used for the gradient:
#p_i(a) = P(descendants of i in data | Y_i = a)
#The p-vectors are stored as a matrix A x 2n - 1,
#columns representing nodes and rows alleles
#the number of alleles A is added on as we discover 
#new possible alleles (from possible overlaps?)
#So we also return a matrix alleles which indices the A alleles
#M is a matrix where M[i, i] is the number of alleles allowed in site i
#and M[i, j] is number of alleles for an overlap going from site i --> j
#(Even for approximate likelihood, we need this for computing transition probabilities)
#if mu = NULL, then uniform dist on M is used
#Otherwise, mu is a list of length S where mu[[i]][[j - i + 1]] for i <= j is a hash storing the probability distribution
#on alleles that go (i, j)
#approx is 0 if exact likelihood is calculated (in case of overlap, every possibility considtered),
#1 if approximate likelihood (in case of overlap, only consider alleles seen in the leaves, or one "wild card")
pq_vectors <- function(D, hap, pairs, edge_lengths, eigenQ, states, states_matrix, M, approx, T_prob, mu = NULL) {
  n = nrow(D)
  S = ncol(D)
  
  #turn ape tree into matched pairs for easier traversal
  p <- hap$p
  alleles <- hap$alleles
  
  #Now go up the tree using the pairs
  if (approx == 0) {
    for (j in 1:(n-1)) {
      
      node = n + j
      
      #The pairs should be stored in decreasing order, but just to be safe...
      #row = which(pairs[, 3] == node) 
      row = n-j
      c1 = pairs[row, 1]
      c2 = pairs[row, 2]
      Pt1 = eigenQ$vectors%*%diag(exp(eigenQ$values*edge_lengths[c1]))%*%eigenQ$inv
      Pt2 = eigenQ$vectors%*%diag(exp(eigenQ$values*edge_lengths[c2]))%*%eigenQ$inv
      
      
      #There must be a better way to do this??
      pos_child1 = which(p[, c1] > 0)
      pos_child2 = which(p[, c2] > 0)
      #if (sum(pos_child1) == 0) { print(c1)}
      #if (sum(pos_child2) == 0) { print(c2)}
      
      
      for (i in 1:length(pos_child1)) {
        child1 = alleles[pos_child1[i], ]
        c1_lumped = lumped_state(child1)
        pc1 = p[pos_child1[i], c1]
        parents1 = possible_parents_FA(child1, states, states_matrix, T_prob, M) 
        for (l in 1:length(pos_child2)) {
          child2 = alleles[pos_child2[l], ] 
          c2_lumped = lumped_state(child2)
          pc2 = p[pos_child2[l], c2]
          parents2 = possible_parents_FA(child2, states, states_matrix, T_prob, M)
          
          possible = Reduce(intersect, list(parents1, parents2))
          if (length(possible) == 0) {
            print(error )
          }
          for (k in 1:length(possible)) {
            allele = possible[[k]]
            allele_lumped = lumped_state(allele)
            
            p1 = transition_prob_finite_alleles(allele, child1, Pt1, states, M, allele_lumped, c1_lumped, mu = mu)
            p2 = transition_prob_finite_alleles(allele, child2, Pt2, states, M, allele_lumped, c2_lumped, mu = mu)
            
            if (p1*p2 == 0) {
              print ("PROBLEM")
              #print(edge_lengths[c1])
              #print(edge_lengths[c2])
            }
            
            allele_index = which(apply(alleles, 1, function(x) all.equal(x[1:S], allele)) == "TRUE")
            if (length(allele_index) > 0) { #The allele has already been added to list
              p[allele_index, node] = p[allele_index, node] + p1*p2*pc1*pc2
            }
            else {
              alleles = rbind(alleles, allele)
              p = rbind(p, rep(0, ncol(p)))
              p[nrow(p), node] = p1*p2*pc1*pc2
            }
            
          }
          
        }
      }
      if (sum(p[, node]) == 0) {
        print("problem")
        print(node)
      }
    }
  }
  else{
    
    for (j in 1:(n-1)) {
      
      node = n + j
      
      #The pairs should be stored in decreasing order, but just to be safe...
      row = which(pairs[, 3] == node) 
      c1 = pairs[row, 1]
      c2 = pairs[row, 2]
      Pt1 = eigenQ$vectors%*%diag(exp(eigenQ$values*edge_lengths[c1]))%*%eigenQ$inv
      Pt2 = eigenQ$vectors%*%diag(exp(eigenQ$values*edge_lengths[c2]))%*%eigenQ$inv
      
      #There must be a better way to do this?? #alleles
      pos_child1 = which(p[, c1] > 0)
      pos_child2 = which(p[, c2] > 0)
      ###I commented these two lines out (julia)
      #if (sum(pos_child1) == 0) { print(c1)}
      #if (sum(pos_child2) == 0) { print(c2)}
      
      
      for (i in 1:length(pos_child1)) {
        child1 = alleles[pos_child1[i], ]
        c1_lumped = lumped_state(child1)
        c1_index = states[[stringi::stri_join(c1_lumped, collapse = "")]]
        
        ##commented out, not used (julia), it is a scalar, don't see the point
        pc1 = p[pos_child1[i], c1]
        if(any(is.na(child1))) {
          print("Na in pq")
          print(child1)
          print(alleles)
          print(pos_child1)
          print(p)
          print(c1)
          print(c2)
          print(pairs)
        }
        
        parents1 = possible_parents_approximate(child1, states, states_matrix, T_prob, D) 
        for (l in 1:length(pos_child2)) {
          child2 = alleles[pos_child2[l], ] 
          c2_lumped = lumped_state(child2)
          c2_index =  states[[stringi::stri_join(c2_lumped, collapse = "")]]
          
          pc2 = p[pos_child2[l], c2]
          
          parents2 = possible_parents_approximate(child2, states, states_matrix, T_prob, D) 
          possible = Reduce(intersect, list(parents1, parents2))
          if (length(possible) == 0) {
            print(error )
          }
          for (k in 1:length(possible)) {
            allele = possible[[k]]
            allele_lumped = lumped_state(allele)
            parent_index = states[[stringi::stri_join(allele_lumped, collapse = "")]]
            
            p1 = transition_prob_finite_alleles(allele, child1, Pt1, states, M, allele_lumped, c1_lumped, parent_index, c1_index, mu = mu)
            p2 = transition_prob_finite_alleles(allele, child2, Pt2, states, M, allele_lumped, c2_lumped, parent_index, c2_index, mu = mu)
            
            if (p1*p2 == 0) {
              print ("PROBLEM")
              print(edge_lengths[c1])
              print(edge_lengths[c2])
            }
            
            allele_index = which(apply(alleles, 1, function(x) all.equal(x[1:S], allele)) == "TRUE")
            if (length(allele_index) > 0) { #The allele has already been added to list
              p[allele_index, node] = p[allele_index, node] + p1*p2*pc1*pc2
            }
            else {
              alleles = rbind(alleles, allele)
              p = rbind(p, rep(0, ncol(p)))
              p[nrow(p), node] = p1*p2*pc1*pc2
            }
            
          }
          
        }
      }
      if (sum(p[, node]) == 0) {
        print("problem")
        print(node)
      }
    }
  }
  
  
  ##========================Upgrade========================##
  # Initialize q matrix
  q = matrix(0, nrow = nrow(alleles), ncol = 2 * n - 1)
  q[1, 2 * n - 1] = 1 # Root as ancestral condition
  
  # Precompute eigenvalue-based matrices for each node to avoid redundant calculations
  Pt_matrices <- lapply(1:(2*n - 2), function(node) {
    eigenQ$vectors %*% diag(exp(eigenQ$values * edge_lengths[node])) %*% eigenQ$inv
  })
  
  for (k in 1:(2*n - 2)) {
    node <- (2*n - 1) - k
    Pt1 <- Pt_matrices[[node]]
    
    row <- c(which(pairs[, 1] == node), which(pairs[, 2] == node))
    parent <- pairs[row, 3]
    sibling <- pairs[row, c(which(pairs[row, 1:2] != node))]
    Pt2 <- Pt_matrices[[sibling]]
    
    # Compute p1_vec for current allele to parent once
    p1_list <- lapply(1:nrow(alleles), function(b) {
      sapply(1:nrow(alleles), function(a) {
        transition_prob_finite_alleles(alleles[b, ], alleles[a, ], Pt1, states, M, mu = mu)
      })
    })
    
    # Compute p2_vec between sibling alleles once
    p2_list <- lapply(1:nrow(alleles), function(b) {
      sapply(1:nrow(alleles), function(c) {
        transition_prob_finite_alleles(alleles[b, ], alleles[c, ], Pt2, states, M, mu = mu)
      })
    })
    
    for (a in 1:nrow(alleles)) {
      #allele <- alleles[a, ]
      q_total <- 0
      
      # Vectorize calculations over alleles to reduce nested loops
      qb_vec <- q[, parent]
      
      for (b in 1:nrow(alleles)) {
        p1_val <- p1_list[[b]][a]
        p2_vals <- p2_list[[b]]
        pc_vec <- p[, sibling]
        
        q_total <- q_total + sum(qb_vec[b] * p1_val * p2_vals * pc_vec)
      }
      q[a, node] <- q_total
    }
  }
  
  return(list(alleles, p, q))
  
}




# pq_vectors <- function(D, hap, pairs, edge_lengths, eigenQ, states, states_matrix, M, approx, T_prob, mu = NULL) {
#   n = nrow(D)
#   S = ncol(D)
#   
#   #turn ape tree into matched pairs for easier traversal
#   #pairs = matching(tree)
#   p <- hap$p
#   alleles <- hap$alleles
#   #T_prob = eigenQ$vectors%*%diag(exp(eigenQ$values))%*%eigenQ$inv
#   #Initialize the p-vectors with the data
#   # alleles = matrix("0", nrow = 1, ncol = S)
#   # p = matrix(0, ncol = 2*n - 1, nrow = 1)
#   # for (r in 1:n) {
#   #   node = D[r, ]
#   #   allele_index = which(apply(alleles, 1, function(x) all.equal(x[1:S], node)) == "TRUE")
#   #   if (length(allele_index) > 0) { #The allele has already been added to list
#   #     p[allele_index, r] = 1
#   #   }
#   #   else {
#   #     alleles = rbind(alleles, node)
#   #     p = rbind(p, rep(0, ncol(p)))
#   #     p[nrow(p), r] = 1
#   #   }
#   # }
#   # 
#   
#   #Now go up the tree using the pairs
#   if (approx == 0) {
#     for (j in 1:(n-1)) {
#       
#       node = n + j
#       
#       #The pairs should be stored in decreasing order, but just to be safe...
#       #row = which(pairs[, 3] == node) 
#       row = n-j
#       c1 = pairs[row, 1]
#       c2 = pairs[row, 2]
#       Pt1 = eigenQ$vectors%*%diag(exp(eigenQ$values*edge_lengths[c1]))%*%eigenQ$inv
#       Pt2 = eigenQ$vectors%*%diag(exp(eigenQ$values*edge_lengths[c2]))%*%eigenQ$inv
#       
#       
#       #There must be a better way to do this??
#       pos_child1 = which(p[, c1] > 0)
#       pos_child2 = which(p[, c2] > 0)
#       #if (sum(pos_child1) == 0) { print(c1)}
#       #if (sum(pos_child2) == 0) { print(c2)}
#       
#       
#       for (i in 1:length(pos_child1)) {
#         child1 = alleles[pos_child1[i], ]
#         c1_lumped = lumped_state(child1)
#         pc1 = p[pos_child1[i], c1]
#         parents1 = possible_parents_FA(child1, states, states_matrix, T_prob, M) 
#         for (l in 1:length(pos_child2)) {
#           child2 = alleles[pos_child2[l], ] 
#           c2_lumped = lumped_state(child2)
#           pc2 = p[pos_child2[l], c2]
#           parents2 = possible_parents_FA(child2, states, states_matrix, T_prob, M)
#           
#           possible = Reduce(intersect, list(parents1, parents2))
#           if (length(possible) == 0) {
#             print(error )
#           }
#           for (k in 1:length(possible)) {
#             allele = possible[[k]]
#             allele_lumped = lumped_state(allele)
#             
#             p1 = transition_prob_finite_alleles(allele, child1, Pt1, states, M, allele_lumped, c1_lumped, mu = mu)
#             p2 = transition_prob_finite_alleles(allele, child2, Pt2, states, M, allele_lumped, c2_lumped, mu = mu)
#             
#             if (p1*p2 == 0) {
#               print ("PROBLEM")
#               #print(edge_lengths[c1])
#               #print(edge_lengths[c2])
#             }
#             
#             allele_index = which(apply(alleles, 1, function(x) all.equal(x[1:S], allele)) == "TRUE")
#             if (length(allele_index) > 0) { #The allele has already been added to list
#               p[allele_index, node] = p[allele_index, node] + p1*p2*pc1*pc2
#             }
#             else {
#               alleles = rbind(alleles, allele)
#               p = rbind(p, rep(0, ncol(p)))
#               p[nrow(p), node] = p1*p2*pc1*pc2
#             }
#             
#           }
#           
#         }
#       }
#       if (sum(p[, node]) == 0) {
#         print("problem")
#         print(node)
#       }
#     }
#   }
#   else{
#     
#     for (j in 1:(n-1)) {
#       
#       node = n + j
#       
#       #The pairs should be stored in decreasing order, but just to be safe...
#       row = which(pairs[, 3] == node) 
#       c1 = pairs[row, 1]
#       c2 = pairs[row, 2]
#       Pt1 = eigenQ$vectors%*%diag(exp(eigenQ$values*edge_lengths[c1]))%*%eigenQ$inv
#       Pt2 = eigenQ$vectors%*%diag(exp(eigenQ$values*edge_lengths[c2]))%*%eigenQ$inv
#       
#       #There must be a better way to do this?? #alleles
#       pos_child1 = which(p[, c1] > 0)
#       pos_child2 = which(p[, c2] > 0)
#       ###I commented these two lines out (julia)
#       #if (sum(pos_child1) == 0) { print(c1)}
#       #if (sum(pos_child2) == 0) { print(c2)}
#       
#       
#       for (i in 1:length(pos_child1)) {
#         child1 = alleles[pos_child1[i], ]
#         c1_lumped = lumped_state(child1)
#         c1_index = states[[stringi::stri_join(c1_lumped, collapse = "")]]
#         
#         ##commented out, not used (julia), it is a scalar, don't see the point
#         pc1 = p[pos_child1[i], c1]
#         if(any(is.na(child1))) {
#           print("Na in pq")
#           print(child1)
#           print(alleles)
#           print(pos_child1)
#           print(p)
#           print(c1)
#           print(c2)
#           print(pairs)
#         }
#         
#         parents1 = possible_parents_approximate(child1, states, states_matrix, T_prob, D) 
#         for (l in 1:length(pos_child2)) {
#           child2 = alleles[pos_child2[l], ] 
#           c2_lumped = lumped_state(child2)
#           c2_index =  states[[stringi::stri_join(c2_lumped, collapse = "")]]
#           
#           pc2 = p[pos_child2[l], c2]
#           
#           parents2 = possible_parents_approximate(child2, states, states_matrix, T_prob, D) 
#           possible = Reduce(intersect, list(parents1, parents2))
#           if (length(possible) == 0) {
#             print(error )
#           }
#           for (k in 1:length(possible)) {
#             allele = possible[[k]]
#             allele_lumped = lumped_state(allele)
#             parent_index = states[[stringi::stri_join(allele_lumped, collapse = "")]]
#             
#             p1 = transition_prob_finite_alleles(allele, child1, Pt1, states, M, allele_lumped, c1_lumped, parent_index, c1_index, mu = mu)
#             p2 = transition_prob_finite_alleles(allele, child2, Pt2, states, M, allele_lumped, c2_lumped, parent_index, c2_index, mu = mu)
#             
#             if (p1*p2 == 0) {
#               print ("PROBLEM")
#               print(edge_lengths[c1])
#               print(edge_lengths[c2])
#             }
#             
#             allele_index = which(apply(alleles, 1, function(x) all.equal(x[1:S], allele)) == "TRUE")
#             if (length(allele_index) > 0) { #The allele has already been added to list
#               p[allele_index, node] = p[allele_index, node] + p1*p2*pc1*pc2
#             }
#             else {
#               alleles = rbind(alleles, allele)
#               p = rbind(p, rep(0, ncol(p)))
#               p[nrow(p), node] = p1*p2*pc1*pc2
#             }
#             
#           }
#           
#         }
#       }
#       if (sum(p[, node]) == 0) {
#         print("problem")
#         print(node)
#       }
#     }
#   }
#   
#   #The q vectors are discovered from the root down
#   q = matrix(0, nrow = nrow(alleles), ncol = 2*n - 1)
#   #q
#   q[1, 2*n-1] = 1 #This condition indicates the root equals ancestral, which
#   #should be the first entry in the alleles matrix anyway
#   
#   for (k in 1:(2*n-2)) {
#     node = (2*n - 1) - k
#     Pt1 = eigenQ$vectors%*%diag(exp(eigenQ$values*edge_lengths[node]))%*%eigenQ$inv
#     #Hacky way since it could be in first or second column
#     row = c(which(pairs[, 1] == node), which(pairs[, 2] == node))
#     parent = pairs[row, 3]
#     sibling = pairs[row, c(which(pairs[row, 1:2] != node))]
#     Pt2 = eigenQ$vectors%*%diag(exp(eigenQ$values*edge_lengths[sibling]))%*%eigenQ$inv
#     
#     #Super inefficient ...
#     for (a in 1:nrow(alleles)) {
#       allele = alleles[a, ]
#       q_total = 0
#       for (b in 1:nrow(alleles)) {
#         allele2 = alleles[b, ]
#         qb = q[b, parent]
#         p1 = transition_prob_finite_alleles(allele2, allele, Pt1, states, M, mu = mu)
#         for (c in 1:nrow(alleles)) {
#           allele3 = alleles[c, ]
#           p2 = transition_prob_finite_alleles(allele2, allele3, Pt2, states, M, mu = mu)
#           pc = p[c, sibling]
#           q_total = q_total + qb*p1*p2*pc
#         }
#       }
#       q[a, node] = q_total
#     }
#     
#   }
#   
#   
#   return(list(alleles, p, q))
#   
# }

##Generates summary table of the data, input for pq_vectors
haplotypes<-function(D,n,S){
  #Because it is function of the data only, we only need to run it once
  #Initialize the p-vectors with the data
  alleles = matrix("0", nrow = 1, ncol = S)
  p = matrix(0, ncol = 2*n - 1, nrow = 1)
  for (r in 1:n) {
    node = D[r, ]
    allele_index = which(apply(alleles, 1, function(x) all.equal(x[1:S], node)) == "TRUE")
    if (length(allele_index) > 0) { #The allele has already been added to list
      p[allele_index, r] = 1
    }
    else {
      alleles = rbind(alleles, node)
      p = rbind(p, rep(0, ncol(p)))
      p[nrow(p), r] = 1
    }
  }
  return(list(alleles=alleles,p=p))
}


#The probability of parent --> child (allele states) along an edge of length t
#NOTE: this uses M as input for the total number of alleles, assuming that each individual
#allele is equally likely 
#transition_prob_finite_alleles <- function(parent, child, Pt, eigenQ, states, M, parent_m = NULL, child_m = NULL) {
transition_prob_finite_alleles <- function(parent, child, Pt, states, M, parent_m = NULL, child_m = NULL, k = NULL, j = NULL, mu = NULL) {
  S = nchar(names(states)[1]) #ncol(states)
  #First compute probability of mutation state transition
  #comment this out (julia) and call it explicitly to avoid this check?
  if (is.null(parent_m)) {
    parent_m = lumped_state(parent)
  }
  if (is.null(child_m)) {
    child_m = lumped_state(child)
  }
  
  if (is.null(k)) {
    #stringi is supposed to be faster than paste?? doesn't seem to make much difference
    k = states[[stringi::stri_join(parent_m, collapse = "")]] #paste(parent_m, collapse = "")]]
  }
  if (is.null(j)) {
    j = states[[stringi::stri_join(child_m, collapse = "")]]
  }
  
  
  #k = which(apply(states, 1, function(x) all.equal(x[1:S], parent_m)) == "TRUE")
  #j = which(apply(states, 1, function(x) all.equal(x[1:S], child_m)) == "TRUE")
  
  #Pt = eigenQ$vectors%*%diag(exp(eigenQ$values*t))%*%eigenQ$inv
  
  prob = Pt[k, j]
  
  
  #If the mutation states are not compatible, zero probability
  if (prob == 0) {
    # print(k)
    #print(j)
    # print("reutnring 0")
    return(0)
  }
  #still have to check the allele states are consistent and multiply by allele mutation probabilities
  s = 1
  while (s < S + 1) {
    
    #If mutation states are the same, allele states should be the same
    if (parent_m[s] == child_m[s]) {
      if (parent[s] != child[s]) {
        return(0)
      }
      s = s + 1
    }
    else {
      #If parent has inactive site, but child site is active, not possible
      #Actually, we don't need to include this because it should already be checked by 
      #the transition probabiliy of the mutation state
      #if (parent_m[s] > 0 & child_m[s] == 0) {return(0)}
      
      if (child_m[s] == 1) { #a new single mutation at site s
        if (is.null(mu)) { #Assume all mutations equally likely
          prob = prob*(1/M[s, s]) 
        }
        else {
          mu_site = mu[[s]][[1]]
          prob = prob*mu_site[[child[s]]]
        }
        s = s + 1
      }
      else { #in this case child_m[s] > 1 means a simultaneous cut
        end_site = max(which(child_m == child_m[s])) #the end site of the mutation
        
        if (is.null(mu)) {
          prob = prob*(1/M[s, end_site])
        }
        else {
          mu_site = mu[[s]][[end_site - s + 1]]
          prob = prob*mu_site[[child[s]]]
        }
        s = end_site + 1
      }
      
    }
  }
  
  
  return(prob)
}

pairwise_prob <- function(a1, a2, a1_lumped, a2_lumped,pa1,pa2, Pt, PT, states, states_matrix, T_prob, D, M, approx = 1, mu = NULL) {
  S = length(a1)

  #First, find all possible parents of each allele state
  #....maybe there's a way to do this more quickly for both at the same time?
  # if (approx == 0) {
  #   parents1 = possible_parents_FA(a1, states, states_matrix, T_prob, M)
  #   parents2 = possible_parents_FA(a2, states, states_matrix, T_prob, M)
  # }
  # else {
  #   parents1 = possible_parents_approximate(a1, states, states_matrix, T_prob, D)
  #   parents2 = possible_parents_approximate(a2, states, states_matrix, T_prob, D)
  # }
  #
  #These are the possible allele states that work for both
  #possible = Reduce(intersect,list(parents1, parents2))
  possible = Reduce(intersect,list(pa1, pa2))

  total_p = 0
  for (k in 1:length(possible)) {
    allele = possible[[k]]
    allele_lumped = lumped_state(allele)

    #Probability of going from the ancestral state (0, 0, 0) to this allele in time T - t
    p = transition_prob_finite_alleles(rep("0", S), allele, PT, states, M, rep(0, S), allele_lumped, mu = mu)

    #Probabilities of going from allele--> a1, a2 in time t
    p1 = transition_prob_finite_alleles(allele, a1, Pt, states, M, allele_lumped, a1_lumped, mu = mu)
    p2 = transition_prob_finite_alleles(allele, a2, Pt, states, M, allele_lumped, a2_lumped, mu = mu)

    if (p1*p2 == 0) {
      print ("PROBLEM")
      print(allele)
      print(p)
      print(PT)
    }

    total_p = total_p + p*p1*p2
  }

  return(total_p)
}


#For this function a1, a2 could be lists, D could be a list, each entry results from independent intBC
#If D is a list, then eigenQ and T_prob are also lists to allow for different rates per intBC
#Returns sum of log probabilites for each intBC
pairwise_prob_wrapper <- function(r, T, a1, a2, a1_lumped, a2_lumped,pa1,pa2, states, states_matrix, eigenQlist, T_problist, D, M, approx = 1, mu_list = NULL) {

  t = T*1/(1 + exp(-r)) #to avoid problems with optimization
  #This should ensure 0 < t < T

  if (is.list(a1)) {
    total_log_prob = 0
    for (k in 1:length(a1)) {

      Dcurr = D[[k]]
      ##slow because it re-computes all per pair
      Pt = eigenQlist[[k]]$vectors%*%diag(exp(eigenQlist[[k]]$values*t))%*%eigenQlist[[k]]$inv
      PT = eigenQlist[[k]]$vectors%*%diag(exp(eigenQlist[[k]]$values*(T - t)))%*%eigenQlist[[k]]$inv

      p = pairwise_prob(a1[[k]], a2[[k]], a1_lumped[[k]], a2_lumped[[k]], pa1[[k]],pa2[[k]],
                        Pt, PT, states, states_matrix, T_problist[[k]], Dcurr, M, approx, mu = mu_list[[k]])

      total_log_prob = total_log_prob + log(p)
    }
  }
  else {

    Pt = eigenQlist$vectors%*%diag(exp(eigenQlist$values*t))%*%eigenQlist$inv
    PT = eigenQlist$vectors%*%diag(exp(eigenQlist$values*(T - t)))%*%eigenQlist$inv

    p = pairwise_prob(a1, a2, a1_lumped, a2_lumped,pa1,pa2, Pt, PT, states, states_matrix, T_problist, D, M, approx, mu = mu_list)
    total_log_prob = log(p)
  }

  return(-total_log_prob) #returns negative b/c optim searches for minimum

}


optim_for_par<-function(index,cell_list,parent_list,t,r0,states,states_matrix, eigenQlist,T_problist, Dlist, M, mu_list){

  a1_lumped = cell_list[[index[1]]][[1]]
  a1 = cell_list[[index[1]]][[2]]
  pa1 = parent_list[[index[1]]]
  a2_lumped = cell_list[[index[2]]][[1]]
  a2 = cell_list[[index[2]]][[2]]
  pa2 = parent_list[[index[2]]]

  mle_t = optim(r0, fn = pairwise_prob_wrapper, T = t, a1 = a1, a2 = a2, a1_lumped = a1_lumped, a2_lumped = a2_lumped, pa1=pa1,pa2=pa2,
                states = states, states_matrix = states_matrix,
                eigenQlist = eigenQlist, T_problist = T_problist, D = Dlist, M = M, approx = 1, mu_list = mu_list, method = "CG") #BFGS

  #Go back from r to t
  r = mle_t$par
  mles = t*1/(1 + exp(-r)) + runif(1, 0, 0.000001)
  return(mles)
}

#For a given mutation state, returns all possible
#allele states, using only alleles which occur in the leaves
#Or a `wildcard allele`

all_allele_states_leaves<-function(m_state, leaves) {
  S = length(m_state)
  
  alleles = c()
  s = 1
  while (s < S + 1) {                                             
    if (m_state[s] == 0) {                                                     
      site_possible = c(0)                                                     
      news = s + 1
      end_site = s                                                             
    }
    
    ##========================Upgrade========================##
    else if (m_state[s] == 1) {                                                 
      leaves_site = leaves[, s]
      possible = c()
      non_zero_leaves <- leaves_site[leaves_site != "0"]
      allele_numbers <- sapply(strsplit(non_zero_leaves, ":"), `[`, 1)
      valid_indices <- which(as.numeric(allele_numbers) <= S)                 
      possible <- non_zero_leaves[valid_indices]
      
      site_possible = c(unique(possible), paste(s, "WC", sep = ""))             
      news = s + 1
      end_site = s                                                             
    }
    
    ##========================Upgrade========================##
    else if (m_state[s] > 1) {                                                 
      end_site = max(which(m_state == m_state[s]))    
      
      # Use Boolean indexing to replace the loop for extracting possible allele states that meet the criteria
      # We select rows in leaves where the values in columns s and end_site are equal and not equal to "0"
      valid_rows <- leaves[, s] == leaves[, end_site] & leaves[, s] != "0"
      possible <- unique(leaves[valid_rows, s])  
      
      m = 10*s + end_site
      site_possible = c(unique(possible), paste(m, "WC", sep = ""))             
      news = end_site + 1                                                       
    }
    
    ##========================Upgrade========================##
    P = length(site_possible)                                                   
    
    if (s == 1) {                                                               
      add_on_matrix <- matrix(rep(site_possible, each = end_site - s + 1), nrow = P, ncol = end_site - s + 1, byrow = TRUE)
      new_matrix <- add_on_matrix
    } else {
      existing_rows <- nrow(alleles)
      
      # Expand alleles by a factor of P in terms of row count
      row_expanded <- alleles[rep(1:existing_rows, each = P), , drop = FALSE]
      
      # Create add_on_matrix to contain all possible state combinations
      add_on_matrix <- matrix(rep(site_possible, each = end_site - s + 1), nrow = P * existing_rows, ncol = end_site - s + 1, byrow = TRUE)
      
      # Use cbind() to concatenate row_expanded and add_on_matrix horizontally
      new_matrix <- cbind(row_expanded, add_on_matrix)
    }
    alleles = new_matrix
    s = news
  }
  return(alleles)                                                               
}

#For a given mutation state, returns all possible
#allele states, using only alleles which occur in the leaves
#Or a `wildcard allele`
# all_allele_states_leaves <- function(m_state, leaves) {
#   S = length(m_state)
#   
#   alleles = c()
#   s = 1
#   while (s < S + 1) {
#     
#     if (m_state[s] == 0) {
#       site_possible = c(0)
#       news = s + 1
#       end_site = s
#     }
#     else if (m_state[s] == 1) {
#       leaves_site = leaves[, s]
#       possible = c()
#       for (j in 1:length(leaves_site)) { #surely a better way to do this
#         if (leaves_site[j] != "0") {
#           temp = strsplit(leaves_site[j], ":")[[1]]
#           A = temp[1]
#           if (as.numeric(A) <= S) {#make sure it's not an overlap
#             possible = c(possible, leaves_site[j])
#           }
#         }
#       }
#       site_possible = c(unique(possible), paste(s, "WC", sep = ""))
#       news = s + 1
#       end_site = s
#     }
#     else if (m_state[s] > 1) {
#       end_site = max(which(m_state == m_state[s]))
#       site_possible = site_alleles(s*10 + end_site, M[s, end_site])
#       
#       possible = c()
#       for (j in 1:nrow(leaves)) { #surely a better way to do this
#         row = leaves[j, ]
#         if (row[s] == row[end_site] & row[s] != "0") {
#           possible = c(possible, row[s])
#         }
#       }
#       m = 10*s + end_site
#       site_possible = c(unique(possible), paste( m, "WC", sep = ""))
#       news = end_site + 1
#     }
#     
#     P = length(site_possible)
#     if (s == 1) { #first site
#       new_matrix = matrix(0, nrow = P, ncol = (news -1))
#       for (j in 1:P) {
#         add_on = rep(site_possible[j], end_site - s + 1)
#         new_matrix[j, ] = add_on
#       }
#     }
#     else {
#       new_matrix = matrix(0, nrow = nrow(alleles)*P, ncol = news - 1)
#       for (r in 1:nrow(alleles)) {
#         row = alleles[r, ]
#         for (j in 1:P) {
#           add_on = rep(site_possible[j], end_site - s + 1)
#           new_state = c(row, add_on)
#           new_matrix[P*(r - 1) + j, ] = new_state
#         }
#       }
#     }
#     alleles = new_matrix
#     s = news
#     #print(alleles)
#   }
#   
#   return(alleles)
# }

possible_parents_approximate <- function(child, states, states_matrix, T_prob, leaves) {
  #print("Leaves in ppa")
  #print(leaves)
  
  lumped_parents = possible_parents_lumped(lumped_state(child), states, states_matrix, T_prob) #This a matrix, rows are states
  lumped_child = lumped_state(child)
  
  S = length(child)
  
  parents = list()
  
  #If the child is ancestral type, only the ancestral type is possible
  if (sum(lumped_child) == 0) {
    return(list(as.character(lumped_parents))) #Note: important because alleles are characters, mutations are integers
  }
  else if (nrow(lumped_parents) > 0) {
    for (r in 1:(nrow(lumped_parents))) {
      state = lumped_parents[r, ]
      
      #This could be a matrix, if there is more than one 
      #possible parent due to masking mutations
      parent = state
      s = 1
      while (s < S + 1) {
        if (lumped_child[s] == state[s]) {
          if (!is.null(nrow(parent))) { #already more possibilities from a previous overlap
            parent[, s] = rep(child[s], nrow(parent))
          }
          else {
            parent[s] = child[s]
          }
          news = s + 1
        }
        else {
          if (lumped_child[s] == 1) {
            if (state[s] > 0) { print("Error in parent state")}
            news = s + 1
          }
          else { #child has a new overlap mutation 
            #This case is if the parent has a mutation which vanished in the child due to an overlapping mutation
            #We will generate all possible configurations, but only
            #from alleles which occurred in the leaves, and so show up elsewhere in the tree
            #We also add on one `wildcard' indicating something not seen in the tree
            #!!!!NOTE!!! There could be more than one WC introduced at the same site
            #Currently the code treats all WC as the same allele, ignoring the possibility (which seems more likely)
            #That the WC are different alleles ... I could modify this later to treat all WC as unique
            endstate = max(which(lumped_child == lumped_child[s]))
            news = endstate + 1
            if (sum(state[s:endstate]) > 0) { #indicating an overlap
              
              temp_parent = rep(0, S)
              temp_parent[s:endstate] = state[s:endstate]
              
              possible_masked = all_allele_states_leaves(temp_parent, leaves) #THIS NEEDS TO CHANGE!
              
              ##========================Upgrade========================##
              if (!is.null(nrow(parent))) {
                new_parents <- vector("list", nrow(parent) * nrow(possible_masked))
                counter <- 1
                for (r in seq_len(nrow(parent))) {
                  for (k in seq_len(nrow(possible_masked))) {
                    new_parent <- parent[r, ]
                    new_parent[s:endstate] <- possible_masked[k, s:endstate]
                    new_parents[[counter]] <- new_parent
                    counter <- counter + 1
                  }
                }
                parent <- do.call(rbind, new_parents)
              } else {
                new_parents <- vector("list", nrow(possible_masked))
                for (k in seq_len(nrow(possible_masked))) {
                  new_parent <- parent
                  new_parent[s:endstate] <- possible_masked[k, s:endstate]
                  new_parents[[k]] <- new_parent
                }
                parent <- do.call(rbind, new_parents)
              }
            }
            
          }
        }
        s = news
      }
      
      if (!is.null(nrow(parent))) { #more than one option
        #print("masked alleles")
        for (r in 1:nrow(parent)) {
          parents[[length(parents) + 1]] <- as.character(parent[r, ])
        }
      }
      else{
        parents[[length(parents) + 1]] <- as.character(parent)
      }
      
    }
    
  }
  return(parents)
}
# possible_parents_approximate <- function(child, states, states_matrix, T_prob, leaves) {
#   #print("Leaves in ppa")
#   #print(leaves)
#   lumped_parents = possible_parents_lumped(lumped_state(child), states, states_matrix, T_prob) #This a matrix, rows are states
#   lumped_child = lumped_state(child)
#   S = length(child)
#   
#   parents = list()
#   
#   #If the child is ancestral type, only the ancestral type is possible
#   if (sum(lumped_child) == 0) {
#     return(list(as.character(lumped_parents))) #Note: important because alleles are characters, mutations are integers
#   }
#   else if (nrow(lumped_parents) > 0) {
#     for (r in 1:(nrow(lumped_parents))) {
#       state = lumped_parents[r, ]
#       
#       #This could be a matrix, if there is more than one 
#       #possible parent due to masking mutations
#       parent = state
#       s = 1
#       while (s < S + 1) {
#         if (lumped_child[s] == state[s]) {
#           if (!is.null(nrow(parent))) { #already more possibilities from a previous overlap
#             parent[, s] = rep(child[s], nrow(parent))
#           }
#           else {
#             parent[s] = child[s]
#           }
#           news = s + 1
#         }
#         else {
#           if (lumped_child[s] == 1) {
#             if (state[s] > 0) { print("Error in parent state")}
#             news = s + 1
#           }
#           else { #child has a new overlap mutation 
#             #This case is if the parent has a mutation which vanished in the child due to an overlapping mutation
#             #We will generate all possible configurations, but only
#             #from alleles which occurred in the leaves, and so show up elsewhere in the tree
#             #We also add on one `wildcard' indicating something not seen in the tree
#             #!!!!NOTE!!! There could be more than one WC introduced at the same site
#             #Currently the code treats all WC as the same allele, ignoring the possibility (which seems more likely)
#             #That the WC are different alleles ... I could modify this later to treat all WC as unique
#             endstate = max(which(lumped_child == lumped_child[s]))
#             news = endstate + 1
#             if (sum(state[s:endstate]) > 0) { #indicating an overlap
#               
#               temp_parent = rep(0, S)
#               temp_parent[s:endstate] = state[s:endstate]
#               possible_masked = all_allele_states_leaves(temp_parent, leaves) #THIS NEEDS TO CHANGE!
#               
#               if (!is.null(nrow(parent))) { #already more possibilities from a previous overlap
#                 results = c()
#                 for (r in 1:nrow(parent)) {
#                   for (k in 1:nrow(possible_masked)) {
#                     new_parent = parent[r, ]
#                     new_parent[s:endstate] = possible_masked[k, s:endstate]
#                     results = c(results, new_parent)
#                   }
#                 }
#                 parent = matrix(results, nrow = nrow(possible_masked)*nrow(parent), ncol = S, byrow = TRUE)
#               }
#               else {
#                 results = c()
#                 for (r in 1:nrow(possible_masked)) {
#                   new_parent = parent
#                   new_parent[s:endstate] = possible_masked[r, s:endstate]
#                   results = c(results, new_parent)
#                 }
#                 parent = matrix(results, nrow = nrow(possible_masked), ncol = S, byrow = TRUE)
#               }
#             }
#             
#           }
#         }
#         s = news
#       }
#       
#       if (!is.null(nrow(parent))) { #more than one option
#         #print("masked alleles")
#         for (r in 1:nrow(parent)) {
#           parents[[length(parents) + 1]] <- as.character(parent[r, ])
#         }
#       }
#       else{
#         parents[[length(parents) + 1]] <- as.character(parent)
#       }
#       
#     }
#     
#   }
#   return(parents)
# }



#Returns the tree (ape format) estimated with the phylotime algorithm,
#from input Dlist
#total time T
phylotime_tree <- function(Dlist, haplist, Ti, states, states_matrix, 
                           eigenQlist, T_problist, M, mu_list = NULL) {
  
  #From the Dlist, get the lists of states for each individual cell
  cell_list = list()
  parent_list = list()
  K = length(Dlist)
  n = nrow(Dlist[[1]])
  for (i in 1:n) {
    a1 = list()
    a1_lumped = list()
    a1_parent = list()
    
    for (k in 1:K) {
      a1[[k]] = Dlist[[k]][i, ]
      a1_lumped[[k]] = lumped_state(a1[[k]]) #the first in the list is lumped state
      a1_parent[[k]] = possible_parents_approximate(a1[[k]], states, states_matrix, T_problist[[k]], Dlist[[k]])
      
    }
    
    a1_total = list(a1_lumped, a1)
    cell_list[[i]] = a1_total
    parent_list[[i]]<-a1_parent
  }
  
  
  #Loop through each pair of cells
  probs = matrix(0, nrow = n, ncol = n)
  mlem = matrix(0, nrow = n, ncol = n)
  converge =  matrix(0, nrow = n, ncol = n)
  vals = matrix(0, nrow = n, ncol = n)
  
  t0 = Ti/3 #Starting for optimization
  #Huh, for some reason this doesn't work for -log prob (gives NA),
  #so I transform variables to optimize for r --> ensures 0 < t < T
  r0 = -log(Ti/t0 - 1)
  if (r0 == 0) {
    print("Problem: r0 = 0")
    print(t0)
    prin(T)
  }
  
  ##Obtain the list of parents once, this is used in pairwise_prob, there is one per K
  
  index<-t(combn(1:n, 2))
  indexlist<-as.list(data.frame(t(index)))
  mles<-mclapply(indexlist,optim_for_par,cell_list=cell_list,parent_list=parent_list,t=Ti,r0=r0,states = states,states_matrix = states_matrix, eigenQlist = eigenQlist, T_problist=T_problist,Dlist=Dlist,M=M,mu_list=mu_list) 
  for (i in 1:nrow(index)){
    mlem[index[i,1],index[i,2]]<-mlem[index[i,2],index[i,1]]<-mles[i][[1]]
  }
  
   #also observe that often MLE is max time T
  
  upgma_tree <- upgma(2*mlem) #or is is 2*(T - mles)?
  currloglik = 0
  return(list(tree = upgma_tree, loglik = currloglik))
}

#From the lumped state of child, returns list of parents that have non-zero probability
#Q is the rate matrix and states is the list of all states
possible_parents_lumped <- function(child, states, states_matrix, T_prob) {
  S = length(child) #number of sites
  
  j = states[[paste(child, collapse = "")]]
  #which(apply(states, 1, function(x) all.equal(x[1:S], child)) == "TRUE")
  
  # T_prob = expm(Q)
  
  parents_indices = which(T_prob[, j] > 0 )
  #Make sure that (0, 0, ...0) is there, though not sure why I need to do this???
  parents_indices = unique(c(parents_indices, nrow(states_matrix))) 
  
  return(states_matrix[parents_indices, ])
}


matching <- function(tree) {
  n = tree$Nnode + 1
  
  #if any edges are 0, then interior nodes might not be ordered which would give an invalid tree
  #a quick fix, just add tiny number to any 0s
  tree$edge.length[which(tree$edge.length == 0)] = 0.0000000001
  #min(tree$edge.length[which(tree$edge.length) != 0])/2
  
  t.tot <- max(ape::node.depth.edgelength(tree)) #the total tree length
  #node.depth.edgelength returns distance from root to each node
  n.t <- t.tot - ape::node.depth.edgelength(tree) #this is distance from leaves to each node
  
  #Sorts according to distance from leaves ... You would hope that the nodes are already labeled
  #in logical order according to their coal time, e.g. interior nodes are ordered
  #(2n-1), (2n-2), ..., (n+1)
  #but apparently this doesn't always happen so we need to check
  xx <- sort(n.t, index.return=TRUE)
  newlabel <- seq(n+1, 2*n-1)
  oldlabel <- xx$ix[xx$ix>n]
  
  for (j in 1:nrow(tree$edge)) {
    #the first entry in each edge is always an interior node
    tree$edge[j, 1] <- newlabel[oldlabel==tree$edge[j,1]]
    if (tree$edge[j, 2] > n) { #if entry is a leaf, don't need to change
      tree$edge[j, 2] <- newlabel[oldlabel==tree$edge[j,2]] 
    }
  }
  
  #the inner nodes are already relabeled, now form the matches
  match = matrix(0, nrow = n-1, ncol = 3)
  for (r in 1:(n-1)) {
    child <- tree$edge[tree$edge[,1] == (2*n-r), 2] #find which children form 2*n - r merge,
    #this will be a vector with two entries
    match[r, ] = c(child, 2*n - r)
  }
  
  return(match)
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


#turns the list into a matrix, and relables
#returns a matrix
state_space_matrix <- function(S) {
  states = state_space_list(S)
  return(relabel(states))
}
library(expm)

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


#branchlengths is a vector of length 2*n-2, where the ith entry is the length of the edge with child i
#Also relabel vertices from what ape does, so that 2n-1 is root
tree_from_pairs <- function(pairs, coal_times){
  n = nrow(pairs) + 1
  edges = matrix(0, nrow = 2*n - 2, ncol = 2)
  edge_lengths = rep(0, 2*n - 2 )
  
  #each pair defines 2 edges
  for (r in 1:nrow(pairs)) {
    edges[2*r - 1, 1] = pairs[r, 3] #parent
    edges[2*r - 1, 2] = pairs[r, 1] #child
    if (pairs[r, 1] <= n) {
      edge_lengths[2*r - 1] = coal_times[pairs[r, 3] - n]
    }
    else {
      edge_lengths[2*r - 1] = coal_times[pairs[r, 3] - n] - coal_times[pairs[r, 1] - n] 
    }
    
    edges[2*r, 1] = pairs[r, 3]
    edges[2*r, 2] = pairs[r, 2]
    if (pairs[r, 2] <= n) {
      edge_lengths[2*r] = coal_times[pairs[r, 3] - n]
    }
    else {
      edge_lengths[2*r] = coal_times[pairs[r, 3] - n] - coal_times[pairs[r, 2] - n] 
    }
    
  }
  
  #relabel the inner nodes: n+1 is now the root
  for (r in 1:nrow(edges)) {
    if (edges[r, 1] > n) {
      edges[r, 1] = 2*n - (edges[r, 1] - n)
    }
    if (edges[r, 2] > n) {
      edges[r, 2] = 2*n - (edges[r, 2] - n)
    }
  }
  
  
  tree <- list(edge = edges, edge.length = edge_lengths)
  tree$tip.label <- paste("t", 1:n, sep = "")
  tree$Nnode <- n - 1
  class(tree) <- "phylo"
  return(tree)
}


#The ape tree interior nodes are labeled with n+1 the root
#I want to rearrange so that the root is 2*n - 1
relabel_interior <- function(alleles) {
  new = alleles
  n = (nrow(alleles)+1)/2
  
  for (j in (n+1):(2*n - 1)) {
    newIndex = (2*n - 1) - j + n + 1
    new[j, ] = alleles[newIndex, ]
  }
  return(new)
}

#From here, used for projection to l1 ball https://rdrr.io/github/vguillemot/sparseMCA/src/R/projl1.R
#' Compute the projection of a vector onto the l1 ball.
#'
#' @param v A vector of numerics
#' @param a The radius (>0) of the l1 ball
#' @return The projection of \eqn{v} onto the l1 ball of radius \eqn{a}.
#' @examples
#' projl1(1:10, 21)
#' @export
projl1 <- function(v, a) {
  #v <- v
  #if (sum(abs(v)) <= a)
  #  return(v)
  u <- sort(abs(v), decreasing = TRUE)
  n <- length(v)
  ukmaok <- (cumsum(u) - a) / (1:n)
  K <- max(which(ukmaok < u))
  tau <- ukmaok[K]
  
  new = sign(v) * pmax(abs(v) - tau, 0)
  
  #Add a small number which might shift it back out of the 
  #ball, but at least we don't get 0 branch lengths
  # for (i in 1:length(new)) {
  #   if (new[i] == 0) {new[i] = 0.0001}
  # }
  return(new)
}

#eta is learning rate
#eps is stopping condition
#hap is the combination of alleles from the data D, which is input in the p, q vectors function to save some time
#Performs gradient ascent for the inter-coal times -- this is done by computing the gradient for branch lengths, stepping in that direction for branch lengths,
#then generating inter-coal times from branch lengths
#ict[i] = intercoalescent-time: difference between (ith) merge and (i-1) merge (so first entry is first merge time) -- length (n-1)
#T is total intercoaltime constraint (if not provided, set to be large = 10)
gradient_ascent_intercoaltimes_adagrad <- function(pairs, D, hap, ict, eta, eps, Q, eigenQ, states, states_matrix, M, T = 10, approx = 1, maxSteps = 1000, T_prob, mu_list = NULL) {
  n = nrow(pairs) + 1
  
  #get the inter-coalescent times -- u_k = t_k - t_{k-1}
  #intercoal = c(coal_times, 0) - c(0, coal_times)
  #intercoal = intercoal[1:(length(intercoal)-1)]
  
  #matrix which stores transformation from intercoal times to branch lenghts -- BL = t(ICT)*A -- A is (n-1)x(2n - 2)
  A = decompose_branches(pairs) 
  currICT = matrix(ict, nrow = n-1, ncol = 1) #store as column vector
  currB = A%*%currICT
  #print(currB)
  
  #branch lengths from coal_times
  #temp_tree = update_time(tree, coal_times)
  #b = branch_lengths(temp_tree)
  
  currC = log(currICT) #We will reparameterize so that b = e^c is positive
  
  ps = c()
  logps=c()
  if (is.list(D)) {
    Qalleles = list()
    plist = list()
    qlist = list()
    currloglik = 0
    for (i in 1:length(D)) {
      currD = D[[i]]
      mu = NULL
      if (!is.null(mu_list)) { mu = mu_list[[i]]}
      currhap = hap[[i]]
      
      starting = pq_vectors(currD, currhap, pairs, currB, eigenQ[[i]], states, states_matrix, M, approx, T_prob[[i]], mu = mu)
      print(i)
      alleles = starting[[1]]
      Qnew = Q_alleles(alleles, Q[[i]], eigenQ[[i]], states, M,T_prob[[i]])
      
      Qalleles[[i]] = Qnew
      p = starting[[2]]
      q = starting[[3]]
      #currP = p[1, ncol(p)]
      plist[[i]] = p
      qlist[[i]] = q
      currloglik = currloglik + log(p[1, ncol(p)])
      print(currloglik)
    }
    currP = exp(currloglik)
    ps = c(ps, currP)
    logps=c(logps,currloglik)
    
  }
  else {
    starting = pq_vectors(D, hap, pairs, currB, eigenQ, states, states_matrix, M, approx, T_prob, mu = mu_list)
    alleles = starting[[1]]
    Qnew = Q_alleles(alleles, Q, eigenQ, states, M,T_prob)
    p = starting[[2]]
    q = starting[[3]]
    
    
    currP = p[1, ncol(p)]
    ps = c(ps, currP)
    logps=c(logps,log(currP))
  }
  
  diff = 1
  numSteps = 1
  
  Gt = rep(0, length(currC))
  while (diff > eps) {
    
    if (is.list(D)) {
      gradB = 0
      for (i in 1:length(D)) {
        gradB = gradB + gradient(plist[[i]], qlist[[i]], Qalleles[[i]])
      }
    }
    else {
      gradB = gradient( p, q, Qnew) #gradient with respect to b
    }
    #print(gradB)
    
    gradICT = t(A)%*%matrix(gradB, nrow = 2*n - 2, ncol = 1) #Since (dL/dICT) = (dL/dB)*(dB/dI), and dB/dI is determined by A
    #update in normal scale
    Gt = Gt + gradICT*gradICT
    if (numSteps==1){rhot=eta}else{rhot = eta*1/sqrt(Gt)}
    # rhot = eta*1/sqrt(Gt)
    oldICT=currICT
    currICT =currICT+rhot*gradICT
    
    
    ##original code
    #print(gradICT)
    
    # #Now for gradient with respect to c, just chain rule
    # gradC = gradICT*currICT #Since b(c) = e^c, b'(c) = e^c = b
    # #Update the learning rate
    # Gt = Gt + gradC*gradC #entry-wise product
    # rhot = eta*1/sqrt(Gt)
    # newC = currC + rhot*gradC #entrywise product
    # #newICT = exp(newC) 
    # oldICT =currICT
    # #newICT = currICT + eta*gradICT
    # currICT = exp(newC) 
    # 
    #Check if the new ICT satisfies constraint (that sum(ICT) < T); 
    #if not, project
    #WARNING: This might set some of the entries to 0, which will not work
    #!!!!!!!
    if ((sum(currICT) > T) & (sum(currICT>0)==(n-1))) {
      #currICT=currICT-(sum(currICT)-T)/(n-1)
      # print("projects") 
      #}
      
      #if (sum(currICT)>T & sum(currICT=<0)>0){
      currICT = projl1(currICT, T)
      #currICT<-currICT*T/sum(currICT)
    }
    
    while (sum(currICT<=0)>0){
      eta<-eta/5
      # print("eta")
      #  print(eta)
      currICT<-oldICT
      if (numSteps==1){rhot=eta}else{rhot = eta*1/sqrt(Gt)}
      currICT =currICT+rhot*gradICT
      if ((sum(currICT) > T) & (sum(currICT>0)==(n-1))) {
        #print("projects") 
        currICT = projl1(currICT, T)
      }
      # return(gradient_ascent_intercoaltimes_adagrad(pairs, D , hap , start, eta = eta, eps = eps, Q = Q, eigenQ = eigenQ, 
      #                                         states_matrix = states_matrix, states = states, M = M,T=t,approx=1,maxSteps = 1000,T_prob))
    }
    
    # if (sum(currICT) > T) {
    #    print(newICT)
    #    currICT = newICT*(T/sum(newICT)) #rescale so sums to T
    #    print(currICT)
    #  }
    
    
    currB = A%*%currICT
    currC = log(currICT)
    
    
    if (is.list(D)) {
      newloglik = 0
      for (i in 1:length(D)) {
        currD = D[[i]]
        mu = NULL
        if (!is.null(mu_list)) { mu = mu_list[[i]]}
        currhap = hap[[i]]
        starting = pq_vectors(currD, currhap, pairs, currB, eigenQ[[i]], states, states_matrix, M, approx,T_prob[[i]], mu = mu)
        
        p = starting[[2]]
        q = starting[[3]]
        #currP = p[1, ncol(p)]
        plist[[i]] = p
        qlist[[i]] = q
        newloglik = newloglik + log(p[1, ncol(p)])
      }
      newP = exp(newloglik)
      diff = sum(abs(logps[length(logps)] - newloglik))
      realdiff<-newloglik-logps[length(logps)]
      #currP = newP
      #ps = c(ps, currP)
      #logps=c(logps,newloglik)
      
      if (realdiff>=0){
        currP = newP
        ps = c(ps, currP)
        logps=c(logps,newloglik)
      }
      while(realdiff<0){
        # print("decrease")
        eta<-eta/5
        # print("eta")
        # print(eta)
        currICT<-oldICT
        if (numSteps==1){rhot=eta}else{rhot = eta*1/sqrt(Gt)}
        currICT =currICT+rhot*gradICT
        if ((sum(currICT) > T) & (sum(currICT>0)==(n-1))) {
          # print("projects") 
          currICT = projl1(currICT, T)
        }
        currB = A%*%currICT
        currC = log(currICT)
        newloglik = 0
        for (i in 1:length(D)) {
          currD = D[[i]]
          currhap = hap[[i]]
          mu = NULL
          if (!is.null(mu_list)) { mu = mu_list[[i]]}
          
          starting = pq_vectors(currD, currhap, pairs, currB, eigenQ[[i]], states, states_matrix, M, approx,T_prob[[i]], mu = mu)
          
          p = starting[[2]]
          q = starting[[3]]
          #currP = p[1, ncol(p)]
          plist[[i]] = p
          qlist[[i]] = q
          newloglik = newloglik + log(p[1, ncol(p)])
        }
        newP = exp(newloglik)
        diff = sum(abs(logps[length(logps)] - newloglik))
        realdiff<-newloglik-logps[length(logps)]
        #return(gradient_ascent_intercoaltimes_adagrad(pairs, D , hap , start, eta = eta, eps = eps, Q = Q, eigenQ = eigenQ, 
        #                                             states_matrix = states_matrix, states = states, M = M,T=t,approx=1,maxSteps = 1000,T_prob))
        
        #return(list(ps, oldICT))
      }
    }
    else {
      starting = pq_vectors(D, hap, pairs, currB,eigenQ, states, states_matrix, M, approx,T_prob, mu = mu_list)
      #alleles = starting[[1]]
      #Qnew = Q_alleles(alleles, eigenQ, states, M)
      p = starting[[2]]
      q = starting[[3]]
      diff = sum(abs(logps[length(logps)] - log(p[1, ncol(p)])))
      
      currP = p[1, ncol(p)]
      ps = c(ps, currP)
      logps=c(logps,log(p[1,ncol(p)]))
    }
    #print(logps)
    
    #print(newICT)
    #print(currP)
    
    
    
    
    if (numSteps %% 10 == 0) { #Print out results every 10 steps
      print(currICT)
      #print(logps[length(logps)])
      print("Likelihood")
      print(currP)
      #print(p[1,ncol(p)])
      print(diff)
      print(numSteps)
      
    }
    numSteps = numSteps + 1
    if (numSteps > maxSteps) {
      break
    }
  }
  
  print("Steps")
  print(numSteps)  
  return(list(ps, currICT))
}

#This new version only computes the parents once, not every time pq_vectors is called
#Julia Sep 3, warning: it only uses approx
#eta is learning rate
#eps is stopping condition
#hap is the combination of alleles from the data D, which is input in the p, q vectors function to save some time
#Performs gradient ascent for the inter-coal times -- this is done by computing the gradient for branch lengths, stepping in that direction for branch lengths,
#then generating inter-coal times from branch lengths
#ict[i] = intercoalescent-time: difference between (ith) merge and (i-1) merge (so first entry is first merge time) -- length (n-1)
#T is total intercoaltime constraint (if not provided, set to be large = 10)
gradient_ascent_intercoaltimes_adagrad2 <- function(pairs, Dlist, haplist, cell_list,parent_list, ict, eta, eps, Qlist, eigenQlist, states, states_matrix, M, T = 10, approx = 1, maxSteps = 1000, T_problist, mu_list = NULL) {
  n = nrow(pairs) + 1
  A = decompose_branches(pairs) 
  currICT = matrix(ict, nrow = n-1, ncol = 1) #store as column vector
  currB = A%*%currICT
  #branch lengths from coal_times
  currC = log(currICT) #We will reparameterize so that b = e^c is positive
  ps = c()
  logps=c()
  if (is.list(Dlist)) {
    Qalleles = list()
    plist = list()
    qlist = list()
    currloglik = 0
    for (i in 1:length(Dlist)) {
      currD = Dlist[[i]]
      mu = NULL
      if (!is.null(mu_list)) { mu = mu_list[[i]]}
      currhap = haplist[[i]]
      currindex = i #index of dataset
      
      #starting = pq_vectors(currD, hap=currhap,  pairs=pairs, edge_lengths=currB, eigenQ=eigenQlist[[i]], states=states, states_matrix=states_matrix, M=M, approx=approx, T_prob=T_problist[[i]], mu = mu)
      
      starting = pq_vectors2(currD, hap=currhap, currindex=currindex, cell_list=cell_list, parent_list=parent_list, pairs=pairs, edge_lengths=currB, eigenQ=eigenQlist[[i]], states=states, states_matrix=states_matrix, M=M, approx=approx, T_prob=T_problist[[i]], mu = mu)
      print(i)
      alleles = starting[[1]]
      Qnew = Q_alleles(alleles, Qlist[[i]],eigenQlist[[i]], states, M,T_problist[[i]])
      
      Qalleles[[i]] = Qnew
      p = starting[[2]]
      q = starting[[3]]
      #currP = p[1, ncol(p)]
      plist[[i]] = p
      qlist[[i]] = q
      currloglik = currloglik + log(p[1, ncol(p)])
      print(currloglik)
    }
    currP = exp(currloglik)
    ps = c(ps, currP)
    logps=c(logps,currloglik)
    
  }
  else {
    starting = pq_vectors(Dlist, haplist, pairs, currB, eigenQlist, states, states_matrix, M, approx, T_problist, mu = mu_list)
    alleles = starting[[1]]
    Qnew = Q_alleles(alleles, Qlist, eigenQlist, states, M,T_problist)
    p = starting[[2]]
    q = starting[[3]]
    
    
    currP = p[1, ncol(p)]
    ps = c(ps, currP)
    logps=c(logps,log(currP))
  }
  
  diff = 1
  numSteps = 1
  
  Gt = rep(0, length(currC))
  while (diff > eps) {
    
    if (is.list(Dlist)) {
      gradB = 0
      for (i in 1:length(Dlist)) {
        gradB = gradB + gradient(plist[[i]], qlist[[i]], Qalleles[[i]])
      }
    }
    else {
      gradB = gradient( p, q, Qnew) #gradient with respect to b
    }
    #print(gradB)
    
    gradICT = t(A)%*%matrix(gradB, nrow = 2*n - 2, ncol = 1) #Since (dL/dICT) = (dL/dB)*(dB/dI), and dB/dI is determined by A
    #update in normal scale
    Gt = Gt + gradICT*gradICT
    if (numSteps==1){rhot=eta}else{rhot = eta*1/sqrt(Gt)}
    # rhot = eta*1/sqrt(Gt)
    oldICT=currICT
    currICT =currICT+rhot*gradICT
    
    
    ##original code
    #print(gradICT)
    
    # #Now for gradient with respect to c, just chain rule
    # gradC = gradICT*currICT #Since b(c) = e^c, b'(c) = e^c = b
    # #Update the learning rate
    # Gt = Gt + gradC*gradC #entry-wise product
    # rhot = eta*1/sqrt(Gt)
    # newC = currC + rhot*gradC #entrywise product
    # #newICT = exp(newC) 
    # oldICT =currICT
    # #newICT = currICT + eta*gradICT
    # currICT = exp(newC) 
    # 
    #Check if the new ICT satisfies constraint (that sum(ICT) < T); 
    #if not, project
    #WARNING: This might set some of the entries to 0, which will not work
    #!!!!!!!
    if ((sum(currICT) > T) & (sum(currICT>0)==(n-1))) {
      #currICT=currICT-(sum(currICT)-T)/(n-1)
      # print("projects") 
      #}
      
      #if (sum(currICT)>T & sum(currICT=<0)>0){
      currICT = projl1(currICT, T)
      #currICT<-currICT*T/sum(currICT)
    }
    
    while (sum(currICT<=0)>0){
      eta<-eta/5
      # print("eta")
      #  print(eta)
      currICT<-oldICT
      if (numSteps==1){rhot=eta}else{rhot = eta*1/sqrt(Gt)}
      currICT =currICT+rhot*gradICT
      if ((sum(currICT) > T) & (sum(currICT>0)==(n-1))) {
        #print("projects") 
        currICT = projl1(currICT, T)
      }
      # return(gradient_ascent_intercoaltimes_adagrad(pairs, D , hap , start, eta = eta, eps = eps, Q = Q, eigenQ = eigenQ, 
      #                                         states_matrix = states_matrix, states = states, M = M,T=t,approx=1,maxSteps = 1000,T_prob))
    }
    
    
    
    currB = A%*%currICT
    currC = log(currICT)
    
    
    if (is.list(Dlist)) {
      newloglik = 0
      for (i in 1:length(Dlist)) {
        currD = Dlist[[i]]
        mu = NULL
        if (!is.null(mu_list)) { mu = mu_list[[i]]}
        currhap = haplist[[i]]
        starting= pq_vectors2(currD, hap=currhap, currindex=i, cell_list=cell_list, parent_list=parent_list, pairs=pairs, edge_lengths=currB, eigenQ=eigenQlist[[i]], states=states, states_matrix=states_matrix, M=M, approx=approx, T_prob=T_problist[[i]], mu = mu)
        
        #starting = pq_vectors(currD, currhap, pairs, currB, eigenQ[[i]], states, states_matrix, M, approx,T_prob[[i]], mu = mu)
        
        p = starting[[2]]
        q = starting[[3]]
        #currP = p[1, ncol(p)]
        plist[[i]] = p
        qlist[[i]] = q
        newloglik = newloglik + log(p[1, ncol(p)])
      }
      newP = exp(newloglik)
      diff = sum(abs(logps[length(logps)] - newloglik))
      realdiff<-newloglik-logps[length(logps)]
      #currP = newP
      #ps = c(ps, currP)
      #logps=c(logps,newloglik)
      
      if (realdiff>=0){
        currP = newP
        ps = c(ps, currP)
        logps=c(logps,newloglik)
      }
      while(realdiff<0){
        # print("decrease")
        eta<-eta/5
        # print("eta")
        # print(eta)
        currICT<-oldICT
        if (numSteps==1){rhot=eta}else{rhot = eta*1/sqrt(Gt)}
        currICT =currICT+rhot*gradICT
        if ((sum(currICT) > T) & (sum(currICT>0)==(n-1))) {
          # print("projects") 
          currICT = projl1(currICT, T)
        }
        currB = A%*%currICT
        currC = log(currICT)
        newloglik = 0
        for (i in 1:length(Dlist)) {
          currD = Dlist[[i]]
          currhap = haplist[[i]]
          mu = NULL
          if (!is.null(mu_list)) { mu = mu_list[[i]]}
          
          # starting = pq_vectors(currD, currhap, pairs, currB, eigenQ[[i]], states, states_matrix, M, approx,T_prob[[i]], mu = mu)
          starting =pq_vectors2(currD, hap=currhap, currindex=i, cell_list=cell_list, parent_list=parent_list, pairs=pairs, edge_lengths=currB, eigenQ=eigenQlist[[i]], states=states, states_matrix=states_matrix, M=M, approx=approx, T_prob=T_problist[[i]], mu = mu)
          
          p = starting[[2]]
          q = starting[[3]]
          #currP = p[1, ncol(p)]
          plist[[i]] = p
          qlist[[i]] = q
          newloglik = newloglik + log(p[1, ncol(p)])
        }
        newP = exp(newloglik)
        diff = sum(abs(logps[length(logps)] - newloglik))
        realdiff<-newloglik-logps[length(logps)]
        #return(gradient_ascent_intercoaltimes_adagrad(pairs, D , hap , start, eta = eta, eps = eps, Q = Q, eigenQ = eigenQ, 
        #                                             states_matrix = states_matrix, states = states, M = M,T=t,approx=1,maxSteps = 1000,T_prob))
        
        #return(list(ps, oldICT))
      }
    }
    else {
      starting = pq_vectors(Dlist, haplist, pairs, currB,eigenQlist, states, states_matrix, M, approx,T_problist, mu = mu_list)
      #alleles = starting[[1]]
      #Qnew = Q_alleles(alleles, eigenQ, states, M)
      p = starting[[2]]
      q = starting[[3]]
      diff = sum(abs(logps[length(logps)] - log(p[1, ncol(p)])))
      
      currP = p[1, ncol(p)]
      ps = c(ps, currP)
      logps=c(logps,log(p[1,ncol(p)]))
    }
    #print(logps)
    
    #print(newICT)
    #print(currP)
    
    
    
    
    if (numSteps %% 10 == 0) { #Print out results every 10 steps
      print(currICT)
      #print(logps[length(logps)])
      print("Likelihood")
      print(currP)
      #print(p[1,ncol(p)])
      print(diff)
      print(numSteps)
      
    }
    numSteps = numSteps + 1
    if (numSteps > maxSteps) {
      break
    }
  }
  
  print("Steps")
  print(numSteps)  
  return(list(ps, currICT))
}


#Returns a matrix of size (2n -2)x(n-1)
#Cols represent intercoalescent times (starting from leaves), ros branch lengths (ith column is branch that ends with node i)
#Each branch length b_j can be written as sum of intercoalescent times
#A[j, i] = 1 if branch j uses intercoal time i
decompose_branches <- function(pairs) {
  #pairs = matching(tree)
  n = nrow(pairs) + 1
  
  A = matrix(0, nrow = n - 1, ncol = 2*n - 2)
  
  #A[1, 1:n] = 1 #The first intercoal-time u_k contributes to the branch lengths for all the leaves
  for (r in 1:(n-1)) { #loop through pairs
    parent = pairs[r, 3]
    
    child1 = pairs[r, 1]
    child2 = pairs[r, 2]
    
    if (child1 <= n) {
      A[1:(parent - n), child1] = 1
    }
    else {
      A[(child1 - n+1):(parent - n), child1] = 1
    }
    
    if (child2 <= n) {
      A[1:(parent - n), child2] = 1
    }
    else {
      A[(child2 - n+1):(parent - n), child2] = 1
    }
  }
  
  return(t(A)) #I decided it would be better to transpose
  
}

#' MLE estimation of coalescent times given tree and data
#' @param Dlist a list  per barcode, of alleles. Each column is a target site
#' @param thetas are the mutation rate parameters per barcode
#' @param tree is the known ape genealogy
#' 
#' @return a vector of inter-coalescent times 
#' 
#' 
inference_times <- function(Dlist, thetas, M = 20,tree, mu_list = NULL) {
  
  numIntBC = length(Dlist)
  n = nrow(Dlist[[1]])
  S = ncol(Dlist[[1]])
  states = state_space(S)
  states_matrix<-state_space_matrix(S)
  
  t0<-coalescent.intervals(tree)$total.depth
  
  
  #If M is not given, estimate an M based on the number of unique alleles observed in dataset
  if (is.null(M)) {
    M = matrix(0, nrow = S, ncol = S)
    for (j in 1:length(Dlist)) {
      M0 = observed_alleles(Dlist[[j]])
      M = pmax(M0, M) #pairwise maximum
    }
    print(M) 
    #add on one to each entry to avoid errors, assume at least one mutation possible
    M = M + matrix(1, nrow = S, ncol = S)
  }else{
    M = matrix(M, nrow = S, ncol = S) #allele types
  }
  
  #Step 1 -- Estimate theta for each intBC
  if (is.null(thetas)) {
    thetas <- list()
    for (j in 1:numIntBC) {
      D = Dlist[[j]]
      
      #This is counting the number of cells which have mutations of each type
      A = observed_mutations2(D, S)$obsM
      B = A/nrow(D)
      
      #Need to store theta as a vector: The first S entries are the single-cut rates,
      #The next entries are the entries of the upper-triangular entries in matrix theta
      theta0 =  rep(1, S*(S+1)/2) 
      theta0 = theta0*0.02
      
      #probs_function computes the sum of (observed - expected)^2 for each mutation type
      est_theta1 = optim(theta0, fn = probs_function,  observedMuts = B, t = t0, S = S, 
                         states = states, states_matrix = states_matrix, lower = rep(0, length(theta0)), method = "L-BFGS-B")
      
      result = est_theta1$par
      #This is stored as a vector, change it back to a matrix
      theta_matrix = matrix(0, nrow = S, ncol = S)
      theta_matrix[lower.tri(theta_matrix, diag=FALSE)] <-result[(S+1):length(result)]
      theta_matrix = t(theta_matrix)
      theta_matrix = theta_matrix + diag(result[1:S])
      
      thetas[[j]] <- theta_matrix
    }
  }
  
  #Create the Q_matrix list and T_prob for each theta
  Q_list <- list()
  eigenQ_list <- list()
  T_prob_list <- list()
  for (j in 1:numIntBC) {
    Q = Q_matrixSL(states_matrix, thetas[[j]])
    eigenQ <- eigen(Q)
    eigenQ$inv <- solve(eigenQ$vectors)
    T_prob = eigenQ$vectors%*%diag(exp(eigenQ$values))%*%eigenQ$inv
    
    Q_list[[j]] = Q
    eigenQ_list[[j]] = eigenQ
    T_prob_list[[j]] = T_prob
  }
  
  #Create haplotypes
  haplist = list()
  for (i in 1:numIntBC) {
    haplist[[i]] = haplotypes(Dlist[[i]], n, S)
  }
  
  # pt_result <- phylotime_tree(Dlist, haplist, t0, states, states_matrix, 
  #                             eigenQ_list, T_prob_list, M, mu_list)
  # 
  # tree <- pt_result$tree
  # 
  #Step 3 -- Run gradient ascent to find best branch lengths
  pairs = matching(tree)
  eta = 0.0001 #learning parameter
  eps = 0.001 #stopping condition 
  #(algorithm will stop when L2 difference of ict from one step to the next is less than eps, or 1000 steps)
  
  ##Initialize GA with the times from upgma
  #start = coalescent.intervals(tree)$interval.length
  start = rep(1, n-1) + runif(n-1, min = -0.2, max = 0.2)
  #start = 1/choose(seq(n,2),2)
  start = start*(t0/sum(start))
  #constrained gradient ascent 
  ga_result = gradient_ascent_intercoaltimes_adagrad(pairs, D = Dlist, hap = haplist, ict = start, eta = eta, eps = eps, Q = Q_list, eigenQ = eigenQ_list, 
                                                     states_matrix = states_matrix, states = states, M = M, T = t0, approx=1, 
                                                     maxSteps = 1000, T_prob = T_prob_list, mu_list = mu_list) # T = t + 0.5)
  print("Gradient ascent")
  print(ga_result[[1]]) #how the likelihood increased
  
  #Package all results
  mle_times = ga_result[[2]]
  coal_times = cumsum(mle_times)
  final_tree = phylodyn:::tree_from_pairs(pairs, coal_times)
  
  result = list(tree = final_tree, 
                curr_ll = log(ga_result[[1]][length(ga_result[[1]])]), mle_times = mle_times,
                thetas = thetas)
  
  return(result)
}

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

#' MLE estimation of topology and coalescent times 
#' @param Dlist a list  per barcode, of alleles. Each column is a target site
#' @param t0 time to the most recent common ancestor
#' @param thetas are the mutation rate parameters per barcode
#' @param M is the number of alleles
#' 
#' @return a genealogy
#' 
#' #Dlist of data (fulle alleles) matrices, corresponding to intBCs -- assume all lists the same
#' Uses the input M parameter (assumed to be the same for every intBC) if provided; otherwise
#' uses the counts of number of unique alleles in the dataset as M (max across intBCs)
#'This function does the full inference procedure:
#'.....1) Infers rates for each intBC (NOT assumed to be the same), if thetas is not input
#'.....2) Gets MLE upgma tree + initial branchlengths
#'.....3) Gradient ascent for branch lengths
inference_all <- function(Dlist, t0, thetas, M = 20,mu_list=NULL) {
  
  numIntBC = length(Dlist)
  n = nrow(Dlist[[1]])
  S = ncol(Dlist[[1]])
  states = state_space(S)
  states_matrix<-state_space_matrix(S)
  
  
  
  #If M is not given, estimate an M based on the number of unique alleles observed in dataset
  if (is.null(M)) {
    M = matrix(0, nrow = S, ncol = S)
    for (j in 1:length(Dlist)) {
      M0 = observed_alleles(Dlist[[j]])
      M = pmax(M0, M) #pairwise maximum
    }
    print(M) 
    #add on one to each entry to avoid errors, assume at least one mutation possible
    M = M + matrix(1, nrow = S, ncol = S)
  }else{
    M = matrix(M, nrow = S, ncol = S) #allele types
  }
  
  #Step 1 -- Estimate theta for each intBC
  if (is.null(thetas)) {
    thetas <- list()
    for (j in 1:numIntBC) {
      D = Dlist[[j]]
      
      #This is counting the number of cells which have mutations of each type
      A = observed_mutations2(D, S)$obsM
      B = A/nrow(D)
      
      #Need to store theta as a vector: The first S entries are the single-cut rates,
      #The next entries are the entries of the upper-triangular entries in matrix theta
      theta0 =  rep(1, S*(S+1)/2) 
      theta0 = theta0*0.02
      
      #probs_function computes the sum of (observed - expected)^2 for each mutation type
      est_theta1 = optim(theta0, fn = probs_function,  observedMuts = B, t = t0, S = S, 
                         states = states, states_matrix = states_matrix, lower = rep(0, length(theta0)), method = "L-BFGS-B")
      
      result = est_theta1$par
      #This is stored as a vector, change it back to a matrix
      theta_matrix = matrix(0, nrow = S, ncol = S)
      theta_matrix[lower.tri(theta_matrix, diag=FALSE)] <-result[(S+1):length(result)]
      theta_matrix = t(theta_matrix)
      theta_matrix = theta_matrix + diag(result[1:S])
      
      thetas[[j]] <- theta_matrix
    }
  }
  
  #Create the Q_matrix list and T_prob for each theta
  Q_list <- list()
  eigenQ_list <- list()
  T_prob_list <- list()
  for (j in 1:numIntBC) {
    Q = Q_matrixSL(states_matrix, thetas[[j]])
    eigenQ <- eigen(Q)
    eigenQ$inv <- solve(eigenQ$vectors)
    T_prob = eigenQ$vectors%*%diag(exp(eigenQ$values))%*%eigenQ$inv
    
    Q_list[[j]] = Q
    eigenQ_list[[j]] = eigenQ
    T_prob_list[[j]] = T_prob
  }
  
  #Create haplotypes
  haplist = list()
  for (i in 1:numIntBC) {
    haplist[[i]] = haplotypes(Dlist[[i]], n, S)
  }
  
  #Step 2 -- Estimate topology + initial times by calculating MLE TMRCA for each pair
  #and in putting into upgma
  pt_result <- phylotime_tree(Dlist, haplist, t0, states, states_matrix, 
                              eigenQ_list, T_prob_list, M, mu_list)
  
  tree <- pt_result$tree
  
  #Step 3 -- Run gradient ascent to find best branch lengths
  pairs = matching(tree)
  eta = 0.0001 #learning parameter
  eps = 0.001 #stopping condition 
  #(algorithm will stop when L2 difference of ict from one step to the next is less than eps, or 1000 steps)
  
  ##Initialize GA with the times from upgma
  start = coalescent.intervals(tree)$interval.length
  
  #constrained gradient ascent 
  ga_result = gradient_ascent_intercoaltimes_adagrad(pairs, D = Dlist, hap = haplist, ict = start, eta = eta, eps = eps, Q = Q_list, eigenQ = eigenQ_list, 
                                                     states_matrix = states_matrix, states = states, M = M, T = t0, approx=1, 
                                                     maxSteps = 1000, T_prob = T_prob_list, mu_list = mu_list) # T = t + 0.5)
  print("Gradient ascent")
  print(ga_result[[1]]) #how the likelihood increased
  
  #Package all results
  mle_times = ga_result[[2]]
  coal_times = cumsum(mle_times)
  final_tree = tree_from_pairs(pairs, coal_times)
  
  result = list(tree = final_tree, inital_ll = pt_result$loglik, initial_times = start,
                curr_ll = log(ga_result[[1]][length(ga_result[[1]])]), mle_times = mle_times,
                thetas = thetas)
  
  return(result)
}




  
  
  
  
  
  