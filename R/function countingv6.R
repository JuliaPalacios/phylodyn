###################################
#### FUNCTION COUNTING ############
###################################
generate_data<-function(n,mu){
  tree1<-rcoal(n,tip.label = NULL, br = "coalescent")
  data1<-simulate_data(mu,tree1)
  oldsuff<-sufficient_stats(data1)
  return(oldsuff)
}


siscount_tree<-function(n,mu,N,oldsuff){
  
  numbnodes<-length(unique(as.vector(oldsuff$nodes[,1:2])))
  numbleafnodes<-length(oldsuff$nodes[oldsuff$nodes[,5]==1,1])
  
  ##Kingman initialization 
  nodesTaj<-oldsuff$nodes 
  nodes<-kingman_perphylo(nodesTaj) #Kingman perfect phylogeny
  sampnodes<-king_initial_sampnodes(nodes) #sampling nodes
  labels<-label_king(nodes,sampnodes)#sampling nodes
  
  
  #SIS initialization
  t<-rep(1,n-1)
  saveTaj<-c()
  saveKing<-c()
  saveLT<-c()
  saveTS<-c()
  
  
  for (i in 1:N){
    #tajima
    result <- python.call("F_sample", oldsuff) 
    where_save2<-python.call("calcPF", result$change_F, result$F_nodes, result$family_size, oldsuff,t,"True")
    saveTaj[i]<-1/where_save2[[1]]
    TS<-tree_shape_discount(result)
    saveTS[i]<-saveTaj[i]/TS
    
    #Kingman
    kinglist<-kingtree_sample_v3(nodes,n,sampnodes,labels)
    saveKing[i]<-1/prod(kinglist$prob)
    LT<-labelled_trees_discount(kinglist)
    saveLT[i]<-c(saveKing[i]/LT)
    
  }
  
  countTaj<-sum(saveTaj)/N
  countKing<-sum(saveKing)/N
  countLT<-sum(saveLT)/N
  countTS<-sum(saveTS)/N
  
  diagTaj<-diagnotics(saveTaj,N)
  diagKing<-diagnotics(saveKing,N)
  diagLT<-diagnotics(saveLT,N)
  diagTS<-diagnotics(saveTS,N)
  
  
  
  return(list(n=n,mu=mu,numbnodes=numbnodes,numbleafnodes=numbleafnodes,countTaj=countTaj,countKing=countKing,
              countTS=countTS,countLT=countLT,cv2Taj=diagTaj$cv2,
              essTaj=diagTaj$ess,seTaj=diagTaj$se,cv2King=diagKing$cv2,
              essKing=diagKing$ess,seKing=diagKing$se,cv2LT=diagLT$cv2,
              essLT=diagLT$ess,seLT=diagLT$se,cv2TS=diagTS$cv2,
              essTS=diagTS$ess,seTS=diagTS$se,
              qTaj=diagTaj$qmax,qKing=diagKing$qmax,qTS=diagTS$qmax,qLT=diagLT$qmax))
}

realtree<-function(n,mu){
  R<-100000
  
  ##sample a data set
  tree1<-rcoal(n,tip.label = NULL, br = "coalescent")
  data1<-simulate_data(mu,tree1)
  oldsuff<-sufficient_stats(data1)
  
  
  ##Kingman initialization 
  nodesTaj<-oldsuff$nodes 
  nodes<-kingman_perphylo(nodesTaj) #Kingman perfect phylogeny
  sampnodes<-king_initial_sampnodes(nodes) #sampling nodes
  labels<-label_king(nodes,sampnodes)#sampling nodes
  

  #cycle initialization
  t<-rep(1,n-1)
  save_Tajshapes<-list()
  save_Tajshapes[[1]]<-matrix(0,nrow=(n-1),ncol=n-1)
  save_Kingshapes<-list()
  save_Kingshapes[[1]]<-matrix(0,nrow=(n-1),ncol=n)
  j=1 #index the number of shapes(the 1st actually doesn't count)
  k=1
  
  
  
  for (i in 1:R){
    #tajima
    result <- python.call("F_sample", oldsuff) 
    result$F<-matrix(unlist(result$F_mat),nrow=length(result$F_mat[[1]]),byrow=TRUE)

    
    #count Tajima 
    a=FALSE
    z=0
    while ((a==TRUE | z==j)==FALSE) {#sweeps the matrix checking two identical objects. 
      z=z+1
      a=identical(save_Tajshapes[[z]],result$F)
    }
    if (z==j & a==FALSE) {#save the new shapes
      j=j+1
      save_Tajshapes[[j]]=result$F
    }
    
    #Kingman
    kinglist<-kingtree_sample_v3(nodes,n,sampnodes,labels)
    
    
    
    #count Kingman 
    a=FALSE
    z=0
    while ((a==TRUE | z==k)==FALSE) {#sweeps the matrix checking two identical objects. 
      z=z+1
      a=identical(save_Kingshapes[[z]],kinglist$kingtree)
    }
    if (z==k & a==FALSE) {#save the new shapes
      k=k+1
      save_Kingshapes[[k]]=kinglist$kingtree
    }
    
    
  }
  
  realTaj=j-1
  
  realKing=k-1
  
  return(list(realTaj=realTaj,realKing=realKing))
}

###################################
### CONVERGENCE CRITERIAS #########
####### + UNCERTAINTY #############

#qmax (Chattarjee Diaconis Section 3)


#effective sample size  cv2 and se (based on section 2Chen et al. 2005) Chatterjee Diaconis metric (1 run)
diagnotics<-function(invprob,N){
cv2_num<-sum((invprob-(sum(invprob)/N))^2)/(N-1)
cv2_den<-(sum(invprob)/N)^2
cv2<-cv2_num/cv2_den
ess<-N/(1+cv2)
se<-sqrt(cv2_num)/sqrt(N)
qmax<-max(invprob)/sum(invprob)
return(list(cv2=cv2,ess=ess,se=se,qmax=qmax))
}

#estimate standard deviation. 


#this sample function is to avoid the behaviour of sample in R odd 
sampleWithoutSurprises <- function(x) {
  if (length(x) <= 1) {
    return(x)
  } else {
    return(sample(x, replace=FALSE))
  }
}


#this sample function is to avoid the behaviour of sample in R odd and consider input probabilities
sampleWithoutSurprisesProb <- function(x,prob) {
  if (length(x) <= 1) {
    return(x)
  } else {
    return(sample(x,replace=FALSE,prob=prob))
  }
}



#################################
## FUNCTION TO SAMPLE KINGMAN'S COALESCENT ##
##################################

#remove the collapsed nodes and label new nodes 
kingman_perphylo<-function(nodes){
listnode<-unique(c(nodes[,1],nodes[,2]))
leafnodes<-nodes[nodes[,5]==1,1]
collap<-nodes[nodes[,4]>1,1]
if (length(collap)>=1){
for (i in 1:length(collap)){
  newline<-c()
  idcollap<-which(nodes[,1]==collap[i])
  for (j in 1:(nodes[idcollap,4]-1)){
    newnode<-max(listnode)+1
    newline<-c(newline, c(newnode,nodes[idcollap,2:3],1,nodes[idcollap,5]))
    listnode<-c(listnode,newnode)
  }
  nodes[idcollap,4]=1
  newmat<-matrix(newline,ncol=5,byrow=TRUE)
  nodes<-rbind(newmat,nodes)
}
}
return(nodes)
}

#create list of nodes to sample from 
king_initial_sampnodes<-function(nodes) {
#09/28/2018: Modified because not all nodes where in sampling nodes
sampnodes<-list()
oldsuff1<-matrix(nodes[which(nodes[,5]==1),],ncol=5, byrow=FALSE)#only leaf nodes
oldsuff2<-matrix(oldsuff1[which(oldsuff1[,3]>1),],ncol=5, byrow=FALSE)#non singleton nodes
oldsuff3<-matrix(oldsuff2[which(oldsuff2[,1] %in% oldsuff2[,2]==FALSE),],ncol=5,byrow=FALSE) #remove internal connected to sampling nodes
col<-nodes[which(nodes[,5]==1 & nodes[,3]==1),2]#parent nodes of the leaf nodes with 1 element
qq<-count(col)
samp<-c(oldsuff3[,1],qq[which(qq[,2]>1),1])
reser<-qq[which(qq[,2]==1),1]
sampnodes[[1]]<-list(samp, reser) #sampling nodes that are available at step 1
names(sampnodes[[1]])<-c("samp","reser")
return(sampnodes)
}


#label the node in the Kingman phylogeny 
label_king<-function(nodes,sampnodes){
#09/28/2018: due to the changes in initial_sampnodes there was a bug here that has been fixed. 
allnodes<-unique(c(nodes[,1],nodes[,2]))
lab<-list()
lab$node<-0
lab$single<-0
lab$tree<-0
labels<-rep(list(lab),length(allnodes))

label<-0
for (i in 1:length(sampnodes[[1]]$samp)) {
  id=sampnodes[[1]]$samp[i]+1
  lab$node<-sampnodes[[1]]$samp[i]
  if (length(nodes[nodes[,1]==sampnodes[[1]]$samp[i] & nodes[,5]==1,3])>=1){#this distinguish between collapsed and/or leaf nodes
  lab$single<-seq(nodes[nodes[,1]==sampnodes[[1]]$samp[i],3])+label
  } else {
  lab$single<-seq(length(which(nodes[,5]==1 & nodes[,3]==1 & nodes[,2]==lab$node)))+label 
  }
  lab$tree<-0
  labels[[id]]<-lab 
  label=max(lab$single)
}
if (length(sampnodes[[1]]$reser)>=1){
for (i in 1:length(sampnodes[[1]]$reser)) {
  id=sampnodes[[1]]$reser[i]+1
  lab$node<-sampnodes[[1]]$reser[i]
  lab$single<-seq(sum(nodes[nodes[,2]==lab$node & nodes[,3]==1,4]))+label
  lab$tree<-0
  labels[[id]]<-lab 
  label<-max(lab$single)
}
}
return(labels)
}


#sample Kingman
kingtree_sample_v3<-function(nodes,n,sampnodes,labels){
  #09/27/18: When choosing the node it accounts for node size
  #09/28/18: Add also that there isn't this half move between joining two trees etc. 
  #09/29/18: I add in the collapsing the possibility that the collapsed node is already there. 
  #09/29/18: I add in the collapsing the possibility that I need to collapse twice.  
  
  #sampling for cycle
  total_vintage<-c()
  indicator<-c()
  change_F<-c()
  prob<-c()
  F_nodes<-list()
  kingtree<-matrix(0,ncol=n,nrow=(n-1))
  
  
  
  for (stp in 1:(n-1)){
    
    samp<-sampnodes[[stp]]$samp
    #compute sampling probabilities.
    probsamp<-c()
    for (i in 1:length(samp)){
      id<-samp[i]+1
      probsamp<-c(probsamp,(length(labels[[id]]$single[labels[[id]]$single>0])+length(labels[[id]]$tree[labels[[id]]$tree>0])))
    }
    normprob<-probsamp/sum(probsamp)
    #sample nodes proportional to their size
    x<-sampleWithoutSurprisesProb(samp,normprob)[1]
    prob<-c(prob,normprob[which(samp==x)])
    singlelist<-labels[[(x+1)]]$single
    singlelist<-singlelist[singlelist!=0]
    sizex<-length(singlelist)
    treelist<-labels[[(x+1)]]$tree #it includes 0s now though
    treelist<-treelist[treelist!=0]
    sizetx<-length(treelist)
    #sample the two elements to coalesce. Note: they are indexed by the entry number in single-/tree-list
    move<-sample((sizex+sizetx),2,replace = FALSE)
    #Look at first element
    sample_sing<-c()
    sample_tree<-c()
    if (move[1]<=sizex){
      sample_sing<-c(sample_sing,singlelist[move[1]])
    } else {
      sample_tree<-c(sample_tree,treelist[(move[1]-sizex)])
    }
    if (move[2]<=sizex){
      sample_sing<-c(sample_sing,singlelist[move[2]])
    } else {
      sample_tree<-c(sample_tree,treelist[(move[2]-sizex)])
    }
    
    prob<-c(prob,2/(sizex+sizetx),1/(sizex+sizetx-1)) #it is 2 because you could have picked anything first. 
    
    ##start iteration to update whether it is a cherry, tree+ sing, 2 trees. 
    if (length(sample_sing)==2){#cherry
      change_F<-c(change_F,2)
      F_nodes[[stp]]<-c(-1,stp)
      kingtree[stp,1:2]<-sort(sample_sing,decreasing=TRUE) #I am labelling obs in this odd way. It has a sense though.
      #remove from the singlelist
      singlelist<-singlelist[-which(singlelist %in% sample_sing)]
    } else if (length(sample_sing)==1) {#singleton+ tree
      change_F<-c(change_F,1)
      F_nodes[[stp]]<-c(sample_tree,stp)
      newtree<-c(sample_sing,kingtree[sample_tree,1:(n-1)])
      kingtree[stp,]<-sort(newtree,decreasing = TRUE)
      #remove from singlelist and treelist
      singlelist<-singlelist[-which(singlelist %in% sample_sing)]
      treelist<-treelist[-which(treelist %in% sample_tree)]
    } else {
      change_F<-c(change_F,0)
      F_nodes[[stp]]<-list(c(sort(sample_tree)[1],stp),c(sort(sample_tree)[2],stp))
      newtree<-c(kingtree[sample_tree[1],],kingtree[sample_tree[2],])
      newtree<-sort(newtree,decreasing = TRUE)[1:n]
      kingtree[stp,]<-sort(newtree,decreasing = TRUE)
      #remove from tree list
      treelist<-treelist[-which(treelist %in% sample_tree)]
    }
    
    
    
    #attach new singlelist and treelist
    labels[[(x+1)]]$single<-singlelist
    labels[[(x+1)]]$tree<-treelist
    #record I have sampled from it.
    total_vintage<-c(total_vintage,x)
    
    #define indicator
    if (x!=0) {
      if (length(kingtree[stp,kingtree[stp,]!=0])==nodes[nodes[,1]==x,3]) {
        indicator<-c(indicator,1)
      } else {
        indicator<-c(indicator,0)
      }
    } else {
      indicator<-c(indicator,0)
    }
    
    #Update sampling nodes now and add the new tree to labels
    if (indicator[length(indicator)]==1){#If I am collapsing WRONG!!
      newsampnodes<-sampnodes[[stp]]$samp[-which(sampnodes[[stp]]$samp==x)]#remove x
      newsampnode<-nodes[which(nodes[,1]==x),2]#find candidate new samp node
      #\\Edit : newsampnode may require to be immediately collapsed. One need to check for that. 
      #Note: if it already in reser, it clearly does not require a 2nd collapsing. 
      if (newsampnode!=0){
        element_newsampnode<-nodes[nodes[,1]==newsampnode,3]*nodes[nodes[,1]==newsampnode,4]
        if (length(kingtree[stp,kingtree[stp,]!=0])==element_newsampnode){
          newsampnode<-nodes[which(nodes[,1]==newsampnode),2]#find candidate 2nd new samp node
        } 
      } #\\
      #Actually collapsing
      if (newsampnode %in% sampnodes[[stp]]$reser) { #if I  have already collapsed once there
        newsampnodes<-c(newsampnodes,newsampnode)
        newsampnodes<-unique(newsampnodes) #doesn't add it if it is there already
        reser<-sampnodes[[stp]]$reser[-which(sampnodes[[stp]]$reser==newsampnode)]
        sampnodes[[(stp+1)]]<-list(newsampnodes, reser)
      }  else if (newsampnode %in% sampnodes[[stp]]$samp) {
        sampnodes[[(stp+1)]]<-list(newsampnodes, sampnodes[[stp]]$reser)
        }else {#never collapsed and it is empty
        sampnodes[[(stp+1)]]<-list(newsampnodes, c(sampnodes[[stp]]$reser,newsampnode))
      }
      #update the labels list
      labels[[(newsampnode+1)]]$tree<-c(labels[[(newsampnode+1)]]$tree,stp)#the label of the new tree is stp
    } else { #If I am not collapsing
      labels[[(x+1)]]$tree<-c(labels[[(x+1)]]$tree,stp) #I add the tree to the list I am sampling from.
      if  ((length(labels[[(x+1)]]$tree)+length(labels[[(x+1)]]$single))>1){#if I should keep sampling from it
        sampnodes[[(stp+1)]]<-sampnodes[[stp]]
      } else {
        newsampnodes<-sampnodes[[stp]]$samp[(sampnodes[[stp]]$samp==x)==FALSE]
        reser<-c(sampnodes[[stp]]$reser,x)
        sampnodes[[(stp+1)]]<-list(newsampnodes, reser)
      }
    }
    
    names(sampnodes[[(stp+1)]])<-c("samp","reser")
    
  }#end of for cycle stp
  
  kingtree<-list(kingtree,prob,change_F,total_vintage,indicator,F_nodes)
  names(kingtree)<-c("kingtree","prob","change_F","total_vintage","indicator","F_nodes")
  return(kingtree)
}

####### END OF KINGMAN'S FUNCTIONS #####

##count the unranked options 

labelled_trees_discount<-function(result){
  
  permt<-c()
  treelist<-c()
  idx_0<-which(result$change_F==0)
  
  #
  if(length(idx_0)!=0){
    for (i in 1:length(idx_0)){
      idt<-idx_0[i]
      tree1<-result$F_nodes[[idt]][[1]][1]
      treelist<-count_vintages(tree1,treelist,result)
      tree2<-result$F_nodes[[idt]][[2]][1]
      treelist<-count_vintages(tree2,treelist,result)
      size1_2<-c(treelist[treelist[,2]==tree1,1],treelist[treelist[,2]==tree2,1])
      permt<-c(permt,factorial(sum(size1_2))/(factorial(size1_2[1])*factorial(size1_2[2])))
      treelist<-rbind(treelist,c((sum(size1_2)+1),idt))
    }
    permu<-prod(permt)
  } else {
    permu<-1
  }
  return(permu)
}

tree_shape_discount<-function(result){
  
  permt<-c()
  treelist<-c()
  idx_0<-which(result$change_F==0)
  
  #
  if(length(idx_0)!=0){
    for (i in 1:length(idx_0)){
      idt<-idx_0[i]
      tree1<-result$F_nodes[[idt]][[1]][1]
      treelist<-count_vintages(tree1,treelist,result)
      tree2<-result$F_nodes[[idt]][[2]][1]
      treelist<-count_vintages(tree2,treelist,result)
      size1_2<-c(treelist[treelist[,2]==tree1,1],treelist[treelist[,2]==tree2,1])
      if (length(unique(size1_2))==2){#the tree have different sizes
        permt<-c(permt,factorial(sum(size1_2))/(factorial(size1_2[1])*factorial(size1_2[2])))
      } else {
        permt<-c(permt,factorial(sum(size1_2))/(2*factorial(size1_2[1])*factorial(size1_2[2])))
      }
      treelist<-rbind(treelist,c((sum(size1_2)+1),idt))
    }
    permu<-prod(permt)
  } else {
    permu<-1
  }
  return(permu)
}

#function to count vintages 
count_vintages<-function(tree,treelist,result){
  if (tree %in% treelist[,2]){
  } else {
    vtree<-min(which(unlist(result$F_nodes)==tree))-1
    j=1
    while (unlist(result$F_nodes)[vtree]!=-1) {
      vtree<-min(which(unlist(result$F_nodes)==unlist(result$F_nodes)[vtree]))-1
      if (unlist(result$F_nodes)[vtree] %in% treelist[,2]){
        j=j+treelist[treelist[,2]==unlist(result$F_nodes)[vtree],1]
      } else {   
        j=j+1
      }
    }
    treelist<-rbind(treelist,c(j,tree))
  }
  
  return(treelist)
}  

