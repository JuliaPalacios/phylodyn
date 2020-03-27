####################################################################################################
######################### TAJIMA HETEROCHRONOUS COALESCENT INFERENCE ################################
####################################################################################################



########################################################
####### Data Simulation & Utilities for plotting #######
########################################################

#' Sample heterochronous topology.
#' 
#' @param SimTimeHet is the output of coalsim
#'   
#' @return heterochronous genealogy
#' @export

sample_genealogy<-function(SimTimeHet){
  n=sum(SimTimeHet$n_sampled)
  NNode=n-1
  addtime<-c() #index of the coalescent events in which I am adding the other sampling groups. 
  for (i in 1:length(SimTimeHet$samp_times)){
    addtime<-c(addtime,max((SimTimeHet$coal_times<SimTimeHet$samp_times[i])*seq(1,NNode)))
  }
  ##### Inputs ####
  #1st column: node name, 2nd col:vintage where it is created (1 if leaf, regardless of sampling time)
  # 3rd col: sampling group
  pool<-cbind(seq(1,sum(SimTimeHet$n_sampled[which(addtime==0)])),rep(1,sum(SimTimeHet$n_sampled[which(addtime==0)])),rep(seq(1,length(addtime[addtime==0])),SimTimeHet$n_sampled[which(addtime==0)]))
  NewickLab<-as.character(seq(1,sum(SimTimeHet$n_sampled[which(addtime==0)])))
  label<-sum(SimTimeHet$n_sampled[which(addtime==0)])
  edge<-c()
  edge.length<-c()
  lab<-c()
  Treenode<-NNode+n
  for (i in 1:NNode){
    s<-sample(pool[,1],2,replace = FALSE)
    lab<-c(lab,s)
    edge<-c(edge,c(Treenode,s[1]),c(Treenode,s[2]))
    edge.length<-c(edge.length,(sum(SimTimeHet$intercoal_times[pool[pool[,1]==s[1],2]:i])-SimTimeHet$samp_times[pool[pool[,1]==s[1],3]]),
                   (sum(SimTimeHet$intercoal_times[pool[pool[,1]==s[2],2]:i])-SimTimeHet$samp_times[pool[pool[,1]==s[2],3]]))
    NewickLab<-c(NewickLab,paste("(",NewickLab[pool[,1]==s[1]],":",edge.length[length(edge.length)-1],",",NewickLab[pool[,1]==s[2]],":",edge.length[length(edge.length)],")",sep=""))
    pool<-rbind(pool,c(Treenode,(i+1),1)) #I put 1 in column 3 because it doesn't matter for node added later.
    #update for next cycle.
    Treenode<-Treenode-1
    ids<-pool[,1]!=s[1] & pool[,1]!=s[2]
    pool<-pool[ids,]
    NewickLab<-NewickLab[ids]
    while (i %in% addtime){
      idx=min(which(addtime==i))
      nadd=SimTimeHet$n_sampled[idx]
      addPool=cbind(seq((label+1),(label+nadd)),rep(1,nadd),rep(idx,nadd))
      addChar<-as.character(seq((label+1),(label+nadd)))
      pool=rbind(pool,addPool)
      NewickLab<-c(NewickLab,addChar)
      addtime[idx]=0
      label=label+nadd
    }
  }
  
  edge<-matrix(edge,ncol=2,byrow=TRUE)
  lab<-lab[lab<=n]
  NewickLab<-paste(NewickLab,";",sep="")
  
  #prepare output
  tree<-list()
  tree$edge<-as.matrix(edge)
  tree$edge.length<-edge.length
  tree$Nnode<-as.integer(NNode)
  tree$tip.label<-as.character(lab)
  tree$Newick<-NewickLab
  
  
  return(tree)
}

## Plot the genealogy
plot_genealogy<-function(tree){
  edges<-tree$edge
  g1<-graph_from_edgelist(edges)
  E(g1)$weight<-tree$edge.length
  g2<-layout_as_tree(g1, root = 11, circular = FALSE, mode = "out", flip.y = TRUE)
  plot(g1,layout=g2,edge.arrow.size=.2)
}


graph_gene_tree<-function(suff){
  #suff<-oldsuff$nodes
  labels<-sort(unique(c(suff[,1],suff[,2])))
  new1<-suff[,1:2]
  for (i in 1:nrow(suff)){
    new1[i,1]<-seq(1,length(labels))[labels==suff[i,1]]
    new1[i,2]<-seq(1,length(labels))[labels==suff[i,2]]
  }
  descendants<-rep(0,length(labels))
  for (j in 1:length(labels)){
    if (sum(new1[,2]==j)>0){
      descendants[j]<-sum(suff[new1[,2]==j,3]*suff[new1[,2]==j,4])
    }else{
      descendants[j]<-sum(suff[new1[,1]==j,3]*suff[new1[,1]==j,4])
    }
  }
  
  library("igraph")
  edges<-cbind(new1[,2],new1[,1])
  g1<-graph_from_edgelist(edges)
  g2<-layout_as_tree(g1, root = 1, circular = FALSE, mode = "out", flip.y = TRUE)
  plot(g1,layout=g2,edge.arrow.size=.2,vertex.color="lightblue", vertex.label.cex=.7, vertex.label=descendants,
       edge.label=suff[,4])
}

graph_gene_tree2<-function(suff,mutations){
  #suff<-oldsuff$nodes
  labels<-sort(unique(c(suff[,1],suff[,2])))
  new1<-suff[,1:2]
  for (i in 1:nrow(suff)){
    new1[i,1]<-seq(1,length(labels))[labels==suff[i,1]]
    new1[i,2]<-seq(1,length(labels))[labels==suff[i,2]]
  }
  descendants<-rep(0,length(labels))
  for (j in 1:length(labels)){
    if (sum(new1[,2]==j)>0){
      descendants[j]<-sum(suff[new1[,2]==j,3]*suff[new1[,2]==j,4])
    }else{
      descendants[j]<-sum(suff[new1[,1]==j,3]*suff[new1[,1]==j,4])
    }
  }
  
  library("igraph")
  edges<-cbind(new1[,2],new1[,1])
  g1<-graph_from_edgelist(edges)
  g2<-layout_as_tree(g1, root = 1, circular = FALSE, mode = "out", flip.y = TRUE)
  plot(g1,layout=g2,edge.arrow.size=.2,vertex.color="lightblue", vertex.label.cex=.7, vertex.label=descendants,
       edge.label=mutations,edge.label.cex=.6,edge.label.color="red",vertex.size=6)
}



###############################
#### Superimpose mutations ####
###############################

#' Sample and superimpose mutations on a heterochronous genealogy
#' 
#' @param mu mutation rate
#' @param tree output of sample_genealogy
#' @param SimTimeHet output of coalsim
#'   
#' @return sequence data in an incidence matrix format(rows are segregating sites, columns are samples)
#' @export

simulate_data_het<-function(mu,tree,SimTimeHet){
  n<-tree$Nnode+1
  t<-sort(c(SimTimeHet$coal_times,SimTimeHet$samp_times))
  interT<-t[2:length(t)]-t[1:(length(t)-1)]
  cT<-SimTimeHet$intercoal_times
  tol<-.0000001
  cT[cT==0L]<-tol
  addtime<-0 #index of the coalescent events in which I am adding the other sampling groups. 
  nLin<-SimTimeHet$lineages
  indicator<-seq(n,2)
  add<-c()
  if (length(SimTimeHet$samp_times)==1){#isochronous
    #do nothing
  } else {#heterochronous
    for (i in 2:length(SimTimeHet$samp_times)){
      addtime<-c(addtime,max((SimTimeHet$coal_times<SimTimeHet$samp_times[i])*seq(1,tree$Nnode)))
      if (addtime[i]==0){#sampling time is before the first coalescent event
        nLin<-append(nLin,nLin[1],(addtime[i]+i-2))
        nLin[1]<-nLin[1]-SimTimeHet$n_sampled[i]
        indicator<-append(indicator,indicator[1],(addtime[i]+i-2))
        add<-c(add,addtime[i]+i)
      } else if (addtime[i]==addtime[i-1]){#two consecutive sampling events
        nLin<-append(nLin,(nLin[addtime[i]+i-2]+SimTimeHet$n_sampled[i-1]),(addtime[i]+i-2))
        indicator<-append(indicator,(indicator[addtime[i]+i-2]),(addtime[i]+i-2))
        add<-c(add,addtime[i]+i)
      } else {
        nLin<-append(nLin,(nLin[addtime[i]+i-2]-1),(addtime[i]+i-2))
        indicator<-append(indicator,(indicator[addtime[i]+i-2]-1),(addtime[i]+i-2))
        add<-c(add,addtime[i]+i) #check
      }
    }
  }
  iL<-sum(cT*seq(n,2)) #total isochrnous length.
  hL<-iL-sum(SimTimeHet$samp_times*SimTimeHet$n_sampled)
  # mu<-20 ##2*N*l, l is the number of loci
  mut<-rpois(1,hL*mu)
  where<-runif(mut,0,hL)
  wheresort<-sort(where)
  stretchedtree<-cumsum(rep(interT,nLin))  #length of all the branches (divi)
  indEpoch<-rep(indicator,nLin)
  indBranch<-rep(seq(length(interT)+1,2),nLin)
  ind2<-matrix(0,nrow=n,ncol=n-1)
  ind3<-matrix(0,nrow=n,ncol=n-1)
  for (i in 1:(n-1)){
    node<-2*n-i  #Claim: 1--12 are the leaf branches
    #(parent,offspring)
    children<-tree$edge[tree$edge[,1]==node,2]
    for (k in children){
      if (k<=n){ #It means that it is a leaf node so it is assigned 1
        ind2[k,i]<-1
      }else{
        ind2[,i]<-ind2[,i]+ind2[,2*n-k]
      }
    }
    # #not 100% sure what it is but we should be able to define in the same way
    if (i==1){
      vv<-ind2[,i]
      id<-which(vv==0)
      vv[id]<-seq(2,n-1)
      ind3[,]<-rep(vv,n-1)
    }else{
      val<-min(ind3[ind2[,i]==1,i-1])
      ind3[ind2[,i]==1,i]<-val
      ind3[ind3[,i-1]>min(ind3[ind3[,i]-ind3[,i-1]<0,i-1]),i]<-ind3[ind3[,i-1]>min(ind3[ind3[,i]-ind3[,i-1]<0,i-1]),i]-1
      ind3[,(i+1):n-1]<-rep(ind3[,i],n-i)
    }
  }
  
  #create leaf branches indicator
  nbranches<-length(indEpoch)
  leaf<-rep(0,nbranches)
  leaf[1:SimTimeHet$n_sampled[1]]<-1
  leafLab<-seq(1,SimTimeHet$n_sampled[1])
  if (length(add)>0){#heterochrnous case
    for (i in 1:length(add)){
      if(addtime[i+1]==0){
        l=sum(nLin[1:(add[i]-1)])+1+nLin[add[i]-1]
        u=sum(nLin[1:add[i]])
        leaf[1:u]<-1
        leafLab<-c(leafLab,seq(1,max(leafLab)),seq(max(leafLab)+1,max(leafLab)+SimTimeHet$n_sampled[i+1]))
      } else{
        l=sum(nLin[1:(add[i]-1)])+1+nLin[add[i]-1]
        u=sum(nLin[1:add[i]])
        leaf[l:u]<-1
        leafLab<-c(leafLab,seq(max(leafLab)+1,max(leafLab)+SimTimeHet$n_sampled[i+1]))
      }
    }
  }
  
  data<-matrix(0,nrow=n,ncol=mut)
  for (j in 1:mut){
    if (wheresort[j]<min(stretchedtree)) {who<-1}
    else{who<-max(seq(1,nbranches)[stretchedtree<wheresort[j]])+1}
    if (leaf[who]==1){data[leafLab[sum(leaf[1:who])],j]<-1} # Note sum(leaf[1:who]) is practically a way to recover the label
    else {
      epoch<-n-indEpoch[who]
      branch_order<-sum(indBranch[1:who]==indBranch[who])
      data[seq(1,n)[ind3[,epoch]==branch_order],j]<-1
    }
    if (j==1){data1<-as.character(paste(data[,j],collapse = '', sep = ''))}
    else{data1<-rbind(data1,as.character(paste(data[,j],collapse = '', sep = '')))}
  }
  
  return(data1)
}



##################################################
##### Perfect phylogeny heterochronous ##### ####
##################################################

sufficient_stats_het<-function(data){
  #Update Dec 2017
  #This is the old function needed for 
  sorted<-sort(data[,1],index.return=T)
  n<-as.integer(nchar(data[1,]))
  groups<-group_data(sorted)
  
  ##Dan Gusfield Algorithm 1.1 and 1.2 
  l<-length(groups$mut.groups)
  O<-matrix(0,nrow=n,ncol=l)
  frequency<-rep(1,n)
  haplotypes<-paste0(rep(0,l),collapse="")
  Index<-matrix(0,nrow=n,ncol=l)
  Lmat<-matrix(0,nrow=n,ncol=l)
  for (j in 1:l){
    O[,j]<-as.numeric(strsplit(groups$mut.groups[l-j+1],NULL)[[1]])
  }
  #L:I believe O is a more traditional view of the incidence matrix 
  Index<- O%*%diag(seq(1,l))  #L:numbers/labels the mutations
  #(3)For each cell
  for (j in 1:n){
    haplotypes<-c(haplotypes,paste0(O[j,],collapse="")) #L:all the haplotypes from the incidence matrix + a vector of 0s (root? which we then remove)
    for (i in 2:l){
      if (O[j,i]==1){Lmat[j,i]<-max(Index[j,1:(i-1)])} #L: Lmat is the closest mutations to the left of the one you are activating
    }
  }
  #correct for multiple haplotypes
  sort.hap<-sort(haplotypes[-1],index.return=T) #L:note: sort.hap
  base<-sort.hap$x[1]
  base.ind<-sort.hap$ix[1]
  remove<-0
  for (j in 2:n){
    if (base==sort.hap$x[j]){
      remove<-c(remove,sort.hap$ix[j])
      frequency[base.ind]<-frequency[base.ind]+1  #L:refers to the vector haplotypes I believe. 
      frequency[sort.hap$ix[j]]<-0
      
    }else{
      base<-sort.hap$x[j]
      base.ind<-sort.hap$ix[j]
    }
  }
  leaves<-apply(Index,1,max) #L:looks like leaves mutations labels
  card<-rev(groups$cardinality)# L: 1 s are not only 1 mutations but the it is a position
  carriers<-rev(groups$carriers) #I don't use this
  y<-0
  mylist<-list(list(x=0,y=0))
  orderlist<-0
  parentlist<-0
  famsize<-0
  numsize<-0
  L<-apply(Lmat,2,max) #this vector has the nesting information, it has parents nodes L: the closest mutations that they gotta share 
  #L: I believe that N divides by the loci because if you have different loci you must be on separate branches
  parents<-sort(unique(L),d=TRUE)
  i<-2
  labSAVE<-c()
  for (j in parents){
    offspring<-seq(1,l)[L==j] #L:this reference to the column, so the loci
    offspringsize<-0
    indic_offspring<-rep(0,length(offspring))
    s<-1
    for (no in offspring){
      offspringsize<-c(offspringsize,sum(leaves==no))
      if (sum(parents==no)>0) {indic_offspring[offspring==no]<-1} # L:this means that there is something attached to it. 
      s<-s+1
    }
    offspringsize<-offspringsize[-1]
    #offspringsize<-frequency[offspring]
    offspringsizelist<-unique(offspringsize) ### I merge nodes with the same size
    #this should only be allowed for leaf nodes.
    ##Per-size, see if there are parents there
    for (k in offspringsizelist){
      s<-1
      #  if  (sum(indic_offspring[offspringsize==k])==0){
      ind2<-indic_offspring[offspringsize==k]
      #those where ind2 is 0
      if (sum(ind2==0)>0){
        lab1<-which(leaves %in% offspring[offspringsize==k & indic_offspring==0])
        mylist[[i]]<-list(x=sum(offspringsize==k & indic_offspring==0),y=card[offspring[offspringsize==k & indic_offspring==0]],z=lab1) 
        labSAVE<-c(labSAVE,lab1)
        famsize<-c(famsize,k)
        y<-c(y,1) #just indicating if it is a leaf
        numsize<-c(numsize,sum(offspringsize==k & indic_offspring==0))
        parentlist<-c(parentlist,j)
        orderlist<-c(orderlist,min(offspring[offspringsize==k & indic_offspring==0]))
        for (le in offspring[offspringsize==k & indic_offspring==0]){
          leaves[leaves==le]<-j
        }
        i<-i+1
      }
      if (sum(ind2==1)>0){
        who<-seq(1:length(offspring))[offspringsize==k & indic_offspring==1]
        for (w in who){
          lab1<-which(leaves %in% offspring[w])
          labSAVE<-c(labSAVE,lab1)
          mylist[[i]]<-list(x=1,y=card[offspring[w]],z=lab1) 
          famsize<-c(famsize,k)
          numsize<-c(numsize,1)
          y<-c(y,0)
          parentlist<-c(parentlist,j)
          orderlist<-c(orderlist,offspring[w])
          leaves[leaves==offspring[w]]<-j
          i<-i+1
          s<-s+1
        }
      }
      
    }
    
  }
  
  mylist[[1]]<-NULL
  suff<-list(mylist=mylist,nodes=cbind(orderlist[-1],parentlist[-1],famsize[-1],numsize[-1],y[-1]))
  ##I need to correct for when the root doesn't add up the total of nodes
  tot2<-cbind(suff$nodes[,2],suff$nodes[,3]*suff$nodes[,4])
  tot<-sum(suff$nodes[suff$nodes[,2]==0,3]*suff$nodes[suff$nodes[,2]==0,4])
  labADD<-sort(seq(1,n)[seq(1,n)%in% unique(labSAVE)==FALSE])
  ######$#### !!!!!!!!!!  ####3
  ### Check for the labels ###
  #print(tot)
  if (tot<n){
    #if there are singletons descending from 0, add it there if not, add a singleton
    if (sum(suff$nodes[suff$nodes[,2]==0 & suff$nodes[,3]==1,1])>0){
      where<-seq(1,nrow(suff$nodes))[suff$nodes[,2]==0 & suff$nodes[,3]==1]
      suff$nodes[where,4]<-suff$nodes[where,4]+n-tot
      suff$mylist[[where]]$x<-suff$mylist[[where]]$x+n-tot
      suff$mylist[[where]]$y<-c(suff$mylist[[where]]$y,rep(0,n-tot))
      suff$mylist[[where]]$z<-c(suff$mylist[[where]]$z,labADD)
    }else{
      #I am not sure I actually need this, actually yes! It is needed
      suff$nodes<-rbind(suff$nodes,c(max(suff$nodes[,1])+1,0,1,n-tot,1))
      where<-nrow(suff$nodes)
      suff$mylist[[where]]<-list(x=n-tot,y=rep(0,n-tot))
    }
  }
  ##I need to make sure everything adds up. If something is off, I need
  #to add 0 mutations
  leaves<-apply(Index,1,max) #L:need it again
  parents<-parents[parents>0]
  for (j in parents){
    toadd<-suff$nodes[suff$nodes[,1]==j,3]*suff$nodes[suff$nodes[,1]==j,4]-sum(tot2[tot2[,1]==j,2])
    ExtraZ=which(leaves==j)
    if (toadd>0){
      #just adds a duplicate
      if (sum(suff$nodes[suff$nodes[,2]==j & suff$nodes[,3]==1,1])>0){#L: the above considers only group with 
        where<-seq(1,nrow(suff$nodes))[suff$nodes[,2]==j & suff$nodes[,3]==1]
        suff$nodes[where,4]<-suff$nodes[where,4]+toadd
        suff$mylist[[where]]$x<-suff$mylist[[where]]$x+toadd
        suff$mylist[[where]]$y<-c(suff$mylist[[where]]$y,rep(0,toadd))
        suff$mylist[[where]]$z<-c(suff$mylist[[where]]$z,ExtraZ)
      }else{#L:it is missing the line
        where<-nrow(suff$nodes)+1
        who<-max(suff$nodes[,1:2])+1
        suff$mylist[[where]]<-list(x=toadd,y=rep(0,toadd),z=ExtraZ)
        suff$nodes<-rbind(suff$nodes,c(who,j,1,toadd,1))#with no mutations are singletons repeated
      }
      #I need to add an else for when I need to complement of different size
    }else{
      #I am not sure I need this
      # suff$nodes<-rbind(suff$nodes,c(max(suff$nodes[,1])+1,j,1,toadd))
      # where<-nrow(suff$nodes)
      # suff$mylist[[where]]<-list(x=toadd,y=rep(0,toadd)) 
    }
    
  }
  return(suff)
}




pp_heterochronous<-function(n_sampled,oldsuff,AddTime){
  #Perfect phylogeny augmented and conditional on the sampling times 
  #create oldsuff function that separates nodes from different sampling groups. 
  oldsuff<-undo_twin_nodes(oldsuff,n_sampled)
  MaxNodeLab<-max(oldsuff$nodes[,c(1,2)])
  samp_group<-rep(seq(1,length(n_sampled)),n_sampled)
  
  #Loop to split nodes that are multiples
  # 1) !!!!!!! there should be a difference between leaf and non leaf nodes.
  nN<-length(oldsuff$mylist)
  ExtraNewCol6<-c()
  NewCol6=rep(0,nN)
  ExtraNewCol7<-c()
  NewCol7=rep(0,nN)
  for (i in 1:nN){
    group.list<-unique(samp_group[unique(oldsuff$mylist[[i]]$z)])
    if (length(group.list)>1 & oldsuff$nodes[i,5]==1){#more groups within the same node
      if (oldsuff$nodes[i,4]>1) {#If it is
        for (z in 1:(length(group.list)-1)) {
          groupLab<-which(samp_group==group.list[z])
          newnodeLab<-oldsuff$mylist[[i]]$z[oldsuff$mylist[[i]]$z %in% groupLab]
          #create new nodes
          newline<-oldsuff$nodes[i,]
          newline[1]<-MaxNodeLab+1
          MaxNodeLab<-MaxNodeLab+1
          if (oldsuff$nodes[i,3]==1){newline[4]<-length(newnodeLab) #you have singletons twin nodes
          } else {newline[4]<-1}
          oldsuff$nodes<-rbind(oldsuff$nodes,newline)
          #update the old line nodes
          if (oldsuff$nodes[i,3]==1){oldsuff$nodes[i,4]<-oldsuff$nodes[i,4]-length(newnodeLab) #you have singletons twin nodes
          } else {oldsuff$nodes[i,4]<-oldsuff$nodes[i,4]-1}
          #create new Mylist entry
          oldsuff$mylist[[(nN+1)]]<-list()
          if (oldsuff$nodes[i,3]==1){ #singletons twin node
            oldsuff$mylist[[(nN+1)]]$x<-length(newnodeLab)
            oldsuff$mylist[[(nN+1)]]$y<-oldsuff$mylist[[(i)]]$y[oldsuff$mylist[[i]]$z %in% newnodeLab]
          } else{
            oldsuff$mylist[[(nN+1)]]$x<-1
            oldsuff$mylist[[(nN+1)]]$y<-oldsuff$mylist[[(i)]]$y[1] #it is one because I am going to remove one at a time
          }
          oldsuff$mylist[[(nN+1)]]$z<-newnodeLab
          nN<-nN+1
          #update the previousMylist entry
          if (oldsuff$nodes[i,3]==1){
            oldsuff$mylist[[(i)]]$x<-oldsuff$mylist[[(i)]]$x-length(newnodeLab)
            oldsuff$mylist[[(i)]]$y<-oldsuff$mylist[[(i)]]$y[(oldsuff$mylist[[i]]$z %in% newnodeLab)==FALSE]
          }else{
            oldsuff$mylist[[(i)]]$x<-oldsuff$mylist[[(i)]]$x-1
            oldsuff$mylist[[(i)]]$y<-oldsuff$mylist[[(i)]]$y[-1] #it is one because I am going to remove one at a time
          }  
          oldsuff$mylist[[(i)]]$z<-oldsuff$mylist[[(i)]]$z[(oldsuff$mylist[[i]]$z %in% newnodeLab)==FALSE]
          #New Column indicator
          ExtraNewCol6<-c(ExtraNewCol6,AddTime[group.list[z]])
          ExtraNewCol7<-c(ExtraNewCol7,group.list[z])
        }
        #create here the vector for the column I am missing(last in the group list)
        NewCol6[i]<-AddTime[group.list[z+1]]
        NewCol7[i]<-group.list[z+1]
      } else if (oldsuff$nodes[i,4]==1){
        oldsuff$nodes[i,5]<-0 #It is not anymore a leaf node
        for (z in 1:length(group.list)) {#no mutations in the group but people from different labels
          groupLab<-which(samp_group==group.list[z])
          newnodeLab<-oldsuff$mylist[[i]]$z[oldsuff$mylist[[i]]$z %in% groupLab]
          #create new nodes
          newline<-oldsuff$nodes[i,]
          newline[1]<-MaxNodeLab+1
          newline[2]<-oldsuff$nodes[i,1]
          MaxNodeLab<-MaxNodeLab+1
          newline[3]<-1
          newline[4]<-length(newnodeLab)
          newline[5]<-1
          oldsuff$nodes<-rbind(oldsuff$nodes,newline)
          #create new Mylist entry
          oldsuff$mylist[[(nN+1)]]<-list()
          oldsuff$mylist[[(nN+1)]]$x<-length(newnodeLab)
          oldsuff$mylist[[(nN+1)]]$y<-rep(0,length(newnodeLab))
          oldsuff$mylist[[(nN+1)]]$z<-newnodeLab
          nN<-nN+1
          #New Column indicator
          ExtraNewCol6<-c(ExtraNewCol6,AddTime[group.list[z]])
          ExtraNewCol7<-c(ExtraNewCol7,group.list[z])
        }
        #create here the vector for the column I am missing(last in the group list)
        NewCol6[i]<-min(AddTime[group.list]) #Think!! Now I am not splitting internal node. I may need to do that. 
        NewCol7[i]<-min(group.list)
      }
    } else {#Only one sampling group or internal node, I create an additional column. 
      NewCol6[i]<-min(AddTime[group.list]) #Think!! Now I am not splitting internal node. I may need to do that. 
      NewCol7[i]<-min(group.list)
    }
  }
  column6<-c(NewCol6,ExtraNewCol6)
  column7<-c(NewCol7,ExtraNewCol7)-1
  oldsuff$nodes<-cbind(oldsuff$nodes,column6,column7)
  colnames(oldsuff$nodes)<-NULL
  rownames(oldsuff$nodes)<-NULL
  
  return(oldsuff)
}

undo_twin_nodes<-function(oldsuff,n_sampled){
  
  #Function that remove all the twin nodes before starting to crete the final PP
  
  MaxNodeLab<-max(oldsuff$nodes[,c(1,2)])
  samp_group<-rep(seq(1,length(n_sampled)),n_sampled)
  
  id_twins<-which(oldsuff$nodes[,3]>1 & oldsuff$nodes[,4]>1)
  
  for ( i in  id_twins){
    qt<-oldsuff$nodes[i,4]
    node_size<-oldsuff$nodes[i,3]
    for (j in 1:(qt-1)){
      #create a new node
      newline<-oldsuff$nodes[i,]
      newline[1]<-MaxNodeLab+1
      MaxNodeLab<-MaxNodeLab+1
      newline[4]<-1
      oldsuff$nodes<-rbind(oldsuff$nodes,newline)
      #update the old line nodes
      oldsuff$nodes[i,4]<-oldsuff$nodes[i,4]-1
      #create new Mylist entry
      entry<-dim(oldsuff$nodes)[1]
      oldsuff$mylist[[entry]]<-list()
      oldsuff$mylist[[entry]]$x<-1
      oldsuff$mylist[[entry]]$y<-oldsuff$mylist[[(i)]]$y[1] #it is one because I am going to remove one at a time
      oldsuff$mylist[[entry]]$z<-oldsuff$mylist[[(i)]]$z[1:node_size]
      #update the old Mylist entry
      oldsuff$mylist[[i]]$x<-oldsuff$mylist[[i]]$x-1
      oldsuff$mylist[[i]]$y<-oldsuff$mylist[[(i)]]$y[-1] #it is one because I am going to remove one at a time
      oldsuff$mylist[[i]]$z<-oldsuff$mylist[[(i)]]$z[-c(seq(1,node_size))]
    }
  }
  colnames(oldsuff$nodes)<-NULL
  rownames(oldsuff$nodes)<-NULL
  return(oldsuff)
}




group_data<-function(sort.main){
  #Update Dec 2017
  #This function summarizes sufficient statistics and also provides site frequency spectra
  #The site frequency spectra can be used to generate a first estimate of Ne and the TMRCA 
  
  mut.groups<-unique(sort.main$x)
  new.label<-seq(1,length(mut.groups))
  cardinality<-rep(0,length(mut.groups))
  carriers<-rep(0,length(mut.groups))
  for (j in 1:length(mut.groups)){
    cardinality[j]<-sum(sort.main$x==mut.groups[j])
    carriers[j]<-sum(as.numeric(strsplit(mut.groups[j],NULL)[[1]]))
  }
  if (sum(cardinality)!=max(sort.main$ix)){print("Error"); break}
  #if (max(carriers)>n){print("Error");break}
  return(list(carriers=carriers,cardinality=cardinality,mut.groups=mut.groups))
}



#######################################################
##### Initialization and/or general functions #########
#######################################################






# function to trim branches from UPGMA tree
trimBranches <- function(tree, n, sign,name_samp) {
  #Define the corr vector
  #Define samp_times of 
  idx<-c()
  for (i in 1:length(tree$tip.label)){
    idx<-c(idx,which(name_samp[,2]==tree$tip.label[i]))
  }
  corr<-as.numeric(name_samp[idx,1])
  if (sign=="minus"){
    corr<-corr*-1
  }
  idLeafEDGE<-which(tree$edge[,2]<=sum(n))
  tree$edge.length[idLeafEDGE]<-tree$edge.length[idLeafEDGE]-corr[tree$edge[idLeafEDGE,2]]
  return(tree)
}  


#correct distance matrix if we are in the heterochrnous case. 
correct_distance_het<-function(samp_times,n,dm1){
  n1<-sum(n)
  samplingdates=rep(samp_times,n)
  c1=matrix(rep(samplingdates,n1),n1,n1,byrow=FALSE)
  c2=matrix(rep(samplingdates,n1),n1,n1,byrow=TRUE)
  correction=c1+c2-2*diag(samplingdates)
  CorrDistance=dm1+correction
  return(CorrDistance)
}  


#construct UPGMA tree, possibly heterochrnous 
upgma_tree<-function(CorrDistance,n,samp_times,name_samp){
  #Create UPGMA tree
  treeUPGMA1<-upgma(CorrDistance)
  if(length(samp_times)>1){#heterochronous samples
    treeUPGMA1_het<-trimBranches(treeUPGMA1,n,"plus",name_samp)
    treeUPGMA1_het$edge.length[treeUPGMA1_het$edge.length<.00001]<-max(min(treeUPGMA1_het$edge.length/10),0.001)
    #correct for possibly negative edges
    while(sum(treeUPGMA1_het$edge.length<=0*1)>=1){
      MinEdge=min(abs(treeUPGMA1_het$edge.length))/10 #MinEdge length in the UPGMA tree: we assign it to the negative edge
      idn=min(which(treeUPGMA1_het$edge.length<=0))
      parNode=treeUPGMA1_het$edge[idn,1]
      idPar=which(treeUPGMA1_het$edge[,2]==parNode)
      idChild=which(treeUPGMA1_het$edge[,1]==parNode) #Look for the child node complementary to the negative one
      chilNode=treeUPGMA1_het$edge[idChild[idChild!=idn],2]
      fill=abs(treeUPGMA1_het$edge.length[idn])+MinEdge #Edge length to compensate
      treeUPGMA1_het$edge.length[idn]=MinEdge #correct the negative edge
      treeUPGMA1_het$edge.length[idChild[idChild!=idn]]=treeUPGMA1_het$edge.length[idChild[idChild!=idn]]+fill #Compensate 
      treeUPGMA1_het$edge.length[idPar]=treeUPGMA1_het$edge.length[idPar]-fill #Compensate 
    }
  } else {treeUPGMA1_het=treeUPGMA1}
  treeUPGMA1_het$edge.length[treeUPGMA1_het$edge.length<.00001]<-max(min(treeUPGMA1_het$edge.length)/10,0.001)
  return(tree=treeUPGMA1_het)
}

#compute coalescent times from the UPGMA tree
coaltimes_upgma<-function(tree,n,samp_times,max.nCoal,name_samp){
  if (length(samp_times)>1){#Heterochronous
    tree_hetRev<-trimBranches(tree,n,"minus",name_samp)
    coal_timesI<-cumsum(coalescent.intervals(tree_hetRev)$interval.length) #It seems like the changes in the tree lengths are not recorded. 
    times=times_lik(coal_timesI,samp_times)
    times$t[times$t[,1]<0.001,1]=0.001
    coal_times<-cumsum(times$t[,1])[times$t[,2]==0]
    times=times_lik(coal_times,samp_times)
    #the cycle belows correct for uncompatible times with max.nCoal
    while (sum((times$addtimes[-1]>max.nCoal[-length(max.nCoal)])*1)>0) {
      coal_times<-cumsum(times$t[,1])[times$t[,2]==0]
      idx<-min(which(times$addtimes[-1]>max.nCoal[-length(max.nCoal)]))
      gap<-min(times$t[,1])/2#arbitrary sampling time interval
      add<-(samp_times[idx+1]-gap-coal_times[max.nCoal[idx]])/max.nCoal[idx]
      coal_times[1:max.nCoal[idx]]<-coal_times[1:max.nCoal[idx]]+seq(1,max.nCoal[idx])*add
      coal_times[(max.nCoal[idx]+1):length(coal_times)]<-coal_times[(max.nCoal[idx]+1):length(coal_times)]+max.nCoal[idx]*add
      times=times_lik(coal_times,samp_times)
      #I am doing for now no correction below.
    }
  } else{#Isochronous
    coal_timesI<-cumsum(coalescent.intervals(tree)$interval.length)
    times=times_lik(coal_timesI,samp_times)
    times$t[times$t[,1]<0.001,1]=0.001
  }
  return(times)
}



#' Sample heterochronous topology.
#' 
#' @param data1 sequence data as sampled by simulate_data_het
#' @param name name of the file that converts sequence data into nucleotide basis
#' @param mu mutation rate
#' @param npoints #grid points to approximate N_e(t)
#' @param fact scaling factor (if time or N_e(t) is rescaled)
#' @param alpha precision parameter cov matrix Gaussian process
#' @param samp_times vector sampling times
#' @param n vector sample size
#' @param input_tree input tree
#' @param name_samp two columns in which one is the sampling times paired with an arbitrary sequence label.
#'   
#' @return initialization list
#' @export
#' 
initial_tajima_het<-function(data1,name="newSim10Bottle",mu,npoints=49,fact=1,alpha=0.2,samp_times,n,input_tree,name_samp){
  #Initialization effective population size, topology, and coalescent times via UPGMA
  
  oldsuff<-sufficient_stats_het(data1)
  max.nCoal=intraCoalCheck(oldsuff,n)
  beastfile_het(data1,name,samp_times,n) 
  if (length(input_tree)==0){#no Initialization tree is given
    tip.labels<-paste(seq(1,sum(n)),rep("_",sum(n)),rep(samp_times,n),sep="")
    name_samp<-cbind(rep(samp_times,n),tip.labels)
    fastaformat<-read.FASTA(name)
    fastafile<-as.phyDat(fastaformat)
    dm <- dist.ml(fastafile)
    treeUPGMA <- upgma(dm)
    treeUPGMA$edge.length[treeUPGMA$edge.length==0]<-max(min(treeUPGMA$edge.length/10,0.001))
    #Correct the UPGMA distance matrix
    rescale=nrow(data1)/(mu*sum(treeUPGMA$edge.length))
    dm1=as.matrix(dm*rescale)
    CorrDistance<-correct_distance_het(samp_times,n,dm1)
    tree<-upgma_tree(CorrDistance,n,samp_times,name_samp)
  } else{
    tree<-input_tree
  }
  #Effective population size estimate
  res2b<-BNPR(data=tree,lengthout=npoints,prec_alpha = alpha)
  #plot_BNPR(res2b)
  #Coalescent times estimate
  times<-coaltimes_upgma(tree,n,samp_times,max.nCoal,name_samp)
  return(list(res2b=res2b,oldsuff=oldsuff,times=times,max.nCoal=max.nCoal[-length(max.nCoal)],theta=c(log(res2b$effpopmean*fact),1)))
}

#' Initial MCMC fuction
#' 
#' @param initial output of initial_tajima_het
#' @param ngrid grid size
#' @param mu mutation rate
#' @param n vector sample size
#' @param alpha precision parameter cov matrix Gaussian process
#' @param fact scaling factor (if time or N_e(t) is rescaled)
#' @param Nate number of initialization steps
#' @param times input initial$times
#' @param samp_times vector sampling times
#'   
#' @return initialization MCMC
#' @export

initial_MCMC_het<-function(initial,ngrid=50,mu,n,alpha=.1,fact,Nate,times,samp_times){
  
  #Initialization MCMC for Tajima inference. 
  
  t<-initial$times$t[,1]
  indt<-initial$times$t[,2]
  nsites<-mu/fact
  oldsuff_het<-pp_heterochronous(n_sampled=n,initial$oldsuff,initial$times$addtimes)
  result <- python.call("F_sample_het", oldsuff_het) 
  likH<-python.call("calcPF_het",result$F_nodes, result$family_size,oldsuff_het,result$node_group,t*nsites,indt,"True")
  currentlik<-likH[[2]]
  for (i in 1:Nate){
    prop <- python.call("F_sample_het", oldsuff_het)
    lik1<-python.call("calcPF_het",prop$F_nodes, prop$family_size,oldsuff_het,prop$node_group,t*nsites,indt,"True")
    if (!is.null(lik1[[2]])){
      if (lik1[[2]]>currentlik){
        result<-prop
        currentlik<-lik1[[2]]
        print(lik1[[2]])
      }
    }  
  }
  #print(currentlik[[2]])
  tmrca<-sum(t)
  grid_bds = range(tmrca, 0); grid = seq(grid_bds[1], grid_bds[2], length.out = ngrid); intl = grid[2]-grid[1]; midpts = grid[-1]-intl/2
  prior_F<-tajimaPrior_het(result$F_nodes,n,initial$times$addtimes) #TO CHANGE!
  #proposal_F<-where_save2[1][[1]]
  probs_list<-c(currentlik,prior_F,0,0)
  times_list<-cumsum(times$t[,1])[times$t[,2]==0] #keep only the coalescent times 
  times_list<-times_list-c(0,times_list[1:(length(times_list)-1)])
  theta_list<-initial$theta[-length(initial$theta)]
  #theta_list<-cbind(theta_list,theta_list)
  prec_list<- 1 ##this is to record tau
  coal_times=cumsum(times$t[,1])[times$t[,2]==0]
  lik_init = coal_lik_init(samp_times, n,coal_times, grid = grid)
  invC <- Q_matrix(as.matrix(midpts), 0, 1); library("spam"); # WHat 's this??
  invC[1,1]<-invC[1,1]+.0001  # WHat 's this??
  eig  = spam::eigen.spam(invC, TRUE) # WHat 's this??
  rtEV = sqrt(eig$values)
  EVC  = eig$vectors
  us1  = U_split(initial$theta,lik_init,invC,.01,.01)$logpos
  dus1 = U_split(initial$theta,lik_init,invC,.01,.01, TRUE)
  Ufun1 = function(theta, lik_init, invC, grad=FALSE) U_split(theta = theta, init = lik_init, invC = invC,
                                                              alpha = alpha, beta = .01, grad = grad)
  tLik=U_times_het(currentlik,n,samp_times,initial$times,nsites, result, initial$theta, grid,const=1,sig=.01)
  currentval<-list(coalprior=0,gaussprior=0,grid=grid,result=result,times=initial$times,theta=initial$theta,us1=us1,dus1=dus1,Ufun1=Ufun1,lik_init=lik_init,invC=invC,rtEV=rtEV,EVC=EVC,currentlik=currentlik,prior_F=prior_F,tLik=tLik)
  return(list(currentval=currentval,probs_list=probs_list,times_list=times_list,theta_list=theta_list,prec_list=prec_list))
}





##############################################
##### Likelihood and Prior functions #########
##############################################


times_lik<-function(coal_times,samp_times){
  #Prepare the times matrix needed for the likelihood. 
  # $t 1st column is intracoalescent/intrasampling times
  # $t 2nd column is 0 if coalescent even. A number to tell which sampling group enters
  
  # $addtime: after how many coalescent events the ith group enters into the sampling
  co<-c()
  pos<-c()
  t<-sort(c(coal_times,samp_times))
  interT<-t[2:length(t)]-t[1:(length(t)-1)]
  time<-rep(0,length(interT))
  if (length(samp_times)>1){
    for (i in 2:length(samp_times)){
      co<-c(co,sum((samp_times[i]>coal_times)*1))
      pos<-c(pos,co[i-1]+i-1)
    }
  }
  time[pos]<-seq(1,length(co))
  times<-list()
  times$t<-cbind(interT,time)
  times$addtimes<-c(0,co)
  return(times)
}



intraCoalCheck<-function(oldsuff,n){
  #
  #This function check by greedy search how many coalescent vents there can be in between sampling points.
  #
  sn=sum(oldsuff$nodes[oldsuff$nodes[,2]==0,3]*oldsuff$nodes[oldsuff$nodes[,2]==0,4]) #total sample size
  max.nCoal=cumsum(n)-1
  if (length(n)>1){
    for (i in 1:(length(n)-1)){
      test=FALSE
      while (test==FALSE){
        addtimes=c(rep(0,i),rep(max.nCoal[i],length(n)-i))
        oldsuff_het<-pp_heterochronous(n_sampled=n,oldsuff,addtimes)
        result <- python.call("F_sample_het", oldsuff_het)
        if (length(result$family_size)>0){
          test=result$family_size[length(result$family_size)]==sn
        } 
        if (test==FALSE){
          max.nCoal[i]=max.nCoal[i]-1
        }
      }
    }
  }
  return(max.nCoal)
}


# tajimaPrior_het<-function(F_nodes,n,addtimes){
#   #Compute Prior Tajima
# 
#   #@first coalescent event
#   ns=n[addtimes==0]
#   nt=sum(ns)
#   vint<-c()
#   prob<-1
#   for (i in 1:(sum(n)-1)){
#     #Transition probability
#     new.ns=ns
#     if (F_nodes[[i]][1]<=0){#if a singleton
#       new.ns[-F_nodes[[i]][1]+1]=ns[-F_nodes[[i]][1]+1]-1
#       if (is.na(ns[-F_nodes[[i]][1]+1])){
#         prob=0
#         break
#       }
#     } else {
#       vint=vint[vint!=F_nodes[[i]][1]]
#     }
#     if (F_nodes[[i]][2]<=0){#if a singleton
#       new.ns[-F_nodes[[i]][2]+1]=ns[-F_nodes[[i]][2]+1]-1
#       if (is.na(ns[-F_nodes[[i]][2]+1])){
#         prob=0
#         break
#       }
#     } else {
#       vint=vint[vint!=F_nodes[[i]][2]]
#     }
#     num=prod(choose(ns,ns-new.ns))
#     denom=choose(nt,2)
#     prob=prob*num/denom
#     #Update
#     ns=new.ns
#     vint<-c(vint,i)
#     if (i %in% addtimes){
#       ns<-c(ns,n[addtimes==i])
#     }
#     nt=sum(ns)+length(vint)
#   }
#   return(prior=prob)
# }


tajimaPrior_het<-function(F_nodes,n,addtimes){
  #Compute Prior Tajima
  
  #@first coalescent event
  ns=n[addtimes==0]
  nt=sum(ns)
  vint<-c()
  prob<-1
  for (i in 1:(sum(n)-1)){
    #Transition probability
    new.ns=ns
    if (F_nodes[[i]][1]<=0){#if a singleton
      new.ns[-F_nodes[[i]][1]+1]=new.ns[-F_nodes[[i]][1]+1]-1
      if (is.na(ns[-F_nodes[[i]][1]+1])){
        prob=0
        break
      }
    } else {
      vint=vint[vint!=F_nodes[[i]][1]]
    }
    if (F_nodes[[i]][2]<=0){#if a singleton
      new.ns[-F_nodes[[i]][2]+1]=new.ns[-F_nodes[[i]][2]+1]-1
      if (is.na(ns[-F_nodes[[i]][2]+1])){
        prob=0
        break
      }
    } else {
      vint=vint[vint!=F_nodes[[i]][2]]
    }
    num=prod(choose(ns,ns-new.ns))
    denom=choose(nt,2)
    prob=prob*num/denom
    #Update
    ns=new.ns
    vint<-c(vint,i)
    if (i %in% addtimes){
      ns<-c(ns,n[addtimes==i])
    }
    nt=sum(ns)+length(vint)
  }
  return(prior=prob)
}



coal_lik_init = function(samp_times, n_sampled, coal_times, grid)
{
  
  # Prepare quantities needed to ompute log likelihood of coalescent model, energy function in HMC algorithms, and metric tensor needed in aMALA.
  
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



U_split = function(theta, init, invC, alpha, beta, grad=FALSE)
{
  #This function has to do with the grid of the Heterochronous Tajima. 
  
  
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



U_times_het = function(newlik,n,samp_times,times,mu, result, theta, grid,const=1,sig=.01)
  
  
  ##Dynamic grid construction for the times + some coal lik computations. 
{
  coal_times=cumsum(times$t[,1])[times$t[,2]==0]
  newtheta<-theta
  oldsize<-length(grid)
  #print(oldsize)
  intl = grid[2]-grid[1]
  ngrid<-ceiling(coal_times[length(coal_times)]/intl)
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
  init = coal_lik_init(samp_times, n, coal_times, grid = grid)
  # print(length(grid))
  #loglik = logliktot[[2]]-sum(seq(length(times)+1,2)*times*nsites)
  loglik = newlik/const
  logpri = coal_loglik(init, newtheta[-length(newtheta)])
  #dloglik=(newlik[[3]])*nsites #gives the gradient wrt times
  # dloglik=(logliktot[[3]]-seq(length(times)+1,2))*nsites #gives the gradient wrt times
  #dlogpri = coal_loglik2(init, newtheta[-length(newtheta)],TRUE) 
  # return(list(rjmcmc=rjmcmc,init=init,dloglik = -dloglik, dlogpri = -dlogpri, dlogpos = -(dloglik+dlogpri),loglik = -loglik, logpri = -logpri, logpos = -(loglik+logpri),newtheta=newtheta,grid=grid))
  return(list(rjmcmc=rjmcmc,init=init,loglik = -loglik, logpri = -logpri, logpos = -(loglik+logpri),newtheta=newtheta,grid=grid))
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



##############################################
###########   MCMC functions    ##############
##############################################

#' HMC step 
#' 
#' @param currentval current val of the chain
#' @param theta_list array in which each column a realization of the chain.
#' @param prec_list vector recording precision parameters 
#' @param probs_list array recording probability 
#' @param const scaling factor
#' @param j 
#' @param nsites mutation parameter
#' @param eps step size
#'   
#' @return MC steps for effective population size
#' @export


updateTheta_het<-function(currentval,theta_list,prec_list,probs_list,const=1,j=1,nsites=nsites,eps=.3){  
  res1 = splitHMC2(currentval$theta, currentval$us1, currentval$dus1, currentval$Ufun1, currentval$lik_init, currentval$invC, currentval$rtEV, currentval$EVC, eps=eps, L=5, rand_leap=TRUE)
  res1$Ind
  if (res1$Ind==1){
    m<-length(res1$q)
    tl<-nrow(as.matrix(theta_list))
    addna<-tl-m+1
    if (addna>=0){theta2<-c(res1$q[-m],rep(NA,addna));theta_list<-cbind(as.matrix(theta_list),theta2)}
    if (addna<0){theta2<-res1$q[-m]; theta_list<-rbind(as.matrix(theta_list),matrix(NA,nrow=abs(addna),ncol=ncol(as.matrix(theta_list))));theta_list<-cbind(theta_list,theta2)}
    
    
    #theta2<-c(res1$q[-m],rep(NA,addna))
    prec_list<-c(prec_list,res1$q[m]) ##precision tau
    currentval$us1<-res1$u;currentval$dus1=res1$du;currentval$theta=res1$q;acp1<-acp1+1;
    
    
    #theta_list<-cbind(theta_list,theta2);
    # if(j==1){
    #   tus2=currentval$Ufun2(currentval$times,nsites,currentval$result,currentval$grid,currentval$theta,logliktot=currentval$tus2$logliktot,TRUE)
    # }else{
    #   tus2=currentval$Ufun2(currentval$times,nsites,currentval$result,currentval$grid,currentval$theta,logliktot=currentval$tus2$logliktot,FALSE)
    # }
    # currentval$us2=tus2$logpos
    # currentval$dus2=tus2$dlogpos
    # currentval$tus2=tus2
    probs_list<-rbind(probs_list,c(currentval$currentlik,currentval$prior_F,-res1$pos_summ$loglik,-res1$pos_summ$logpri))
  }else{
    m<-length(currentval$theta)
    tl<-nrow(as.matrix(theta_list))
    addna<-tl-m+1
    if (addna>=0){theta2<-c(res1$q[-m],rep(NA,addna));theta_list<-cbind(theta_list,theta2)}
    if (addna<0){theta2<-res1$q[-m]; theta_list<-rbind(theta_list,matrix(NA,nrow=abs(addna),ncol=ncol(theta_list)));theta_list<-cbind(theta_list,theta2)}
    
    #theta2<-c(currentval$theta[-m],rep(NA,addna))
    #theta_list<-cbind(theta_list,theta2);
    prec_list<-c(prec_list,currentval$theta[m])
    probs_list<-rbind(probs_list,c(currentval$currentlik,currentval$prior_F,-res1$pos_summ$loglik,-res1$pos_summ$logpri))
  }
  currentval$coalprior<--res1$pos_summ$loglik
  currentval$gaussprior<--res1$pos_summ$logpri
  return(list(currentval=currentval,probs_list=probs_list,prec_list=prec_list,theta_list=theta_list, acp=res1$Ind))
}



#HMC for the effective pop size
splitHMC2 = function (q_cur, u_cur, du_cur, U, lik_init,invC, rtEV, EVC, eps=.1, L=5, rand_leap=TRUE)
{
  # initialization
  q = q_cur
  D = length(q)
  u = u_cur
  du = du_cur
  
  # sample momentum
  p = stats::rnorm(D)
  
  # calculate current energy
  E_cur = u + sum(p^2)/2
  
  
  if (rand_leap)
    randL = ceiling(stats::runif(1)*L)
  else
    randL = ceiling(L)
  
  p = p - eps/2*du
  qT = rtEV*(t(EVC)%*%q[-D])
  pT = t(EVC)%*%p[-D]
  A = t(qT)%*%qT
  # Alternate full steps for position and momentum
  for (l in 1:randL)
  {
    p[D] <- p[D] - eps/2*A/2*exp(q[D])
    q[D] <- q[D] + eps/2*p[D]
    
    # Make a full step for the middle dynamics
    Cpx = complex(modulus = 1, argument = -rtEV*exp(q[D]/2)*eps)*complex(real = qT*exp(q[D]/2), imaginary = pT)
    qT = Re(Cpx)*exp(-q[D]/2)
    pT = Im(Cpx)
    q[-D] = EVC%*%(qT/rtEV)
    
    # Make a half step for the last half dynamics
    A=t(qT)%*%qT
    
    q[D] <- q[D] + eps/2*p[D]
    p[D] <- p[D] - eps/2*A/2*exp(q[D])
    
    du = U(q,lik_init,invC, grad = TRUE)
    if(l!=randL)
    {
      pT = pT - eps*(t(EVC)%*%du[-D])
      p[D] = p[D] - eps*du[D]
    }
  }
  p[-D] = EVC%*%pT - eps/2*du[-D]
  p[D] = p[D] - eps/2*du[D]
  
  # Evaluate potential and kinetic energies at start and end of trajectory
  pos_summ = U(q,lik_init,invC)
  u = pos_summ$logpos
  E_prp = u + sum(p^2)/2
  
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  logAP = -E_prp + E_cur
  
  if( is.finite(logAP) && (log(stats::runif(1))<min(0,logAP)) )
    return (list(q = q, u = u, du = du, Ind = 1, pos_summ = pos_summ))
  else
    return (list(q = q_cur, u = u_cur, du = du_cur, Ind = 0, pos_summ = U(q_cur,lik_init,invC)))
}


#' Global Local step for times 
#' 
#' @param n vector of samples size
#' @param oldsuff perfect phylogeny
#' @param currentval current val of the chain
#' @param times_list array in which each column a realization of the chain.
#' @param theta_list array in which each column a realization of the chain.
#' @param prec_list vector recording precision parameters 
#' @param probs_list array recording probability 
#' @param const scaling factor
#' @param mu mutation parameter
#' @param samp_times vector of sampling times
#' @param a Z1
#' @param sigma  sdev of the truncated normal.
#'   
#' @return MC steps for coalescent times
#' @export

updateTimes_hetLocal<-function(n,oldsuff,currentval,times_list,theta_list,prec_list,probs_list,const=1,mu,samp_times,a=2,sigma=0.2){  
  
  #New proposal
  proposal=timesLocalProposal(currentval,initial,a,sigma)  #TO REDO
  propTimes=proposal$propTimes
  #compute quantities for MH ratios
  #1)Likelihood
  newoldsuff=pp_heterochronous(n_sampled=n,oldsuff,proposal$propTimes$addtimes)
  newlik=python.call("calcPF_het",currentval$result$F_nodes, currentval$result$family_size,newoldsuff,currentval$result$node_group,proposal$propTimes$t[,1]*mu,proposal$propTimes$t[,2],"True")
  #2)Prior & posterior
  if (!is.null(newlik[[2]])){
    updatePos=U_times_het(newlik[[2]],n,samp_times,propTimes,mu, currentval$result, currentval$theta, currentval$grid,const=1,sig=.01)
    logMHR=min(0,-updatePos$logpos+log(proposal$revDen)-currentval$currentlik-currentval$coalprior-log(proposal$forDen))
    if (log(runif(1))<=logMHR){Ind=1} else {Ind=0}
  } else {Ind=0}
  #Updates depending on the acceptance or not.
  if (Ind==1){ #Note: I have not been updating at all tLik
    currentval$times=propTimes
    add_times_list<-cumsum(propTimes$t[,1])[propTimes$t[,2]==0] #keep only the coalescent times 
    add_times_list<-add_times_list-c(0,add_times_list[1:(length(add_times_list)-1)])
    times_list<-cbind(times_list,add_times_list);
    currentval$theta<-updatePos$newtheta  #why does this one changes??
    tmrca<-sum(propTimes$t[,1])
    if (max(currentval$grid)<tmrca) {
      #break
      currentval$grid<-updatePos$grid
      newl<-length(currentval$grid)-1
      oldl<-nrow(as.matrix(theta_list))
      intl = currentval$grid[2]-currentval$grid[1]; 
      currentval$grid[length(currentval$grid)]<-tmrca
      midpts<-((c(0,currentval$grid[-length(currentval$grid)])+currentval$grid)/2)[-1]
      #midpts = grid[-1]-intl/2
      currentval$invC <- Q_matrix(as.matrix(midpts), 0, 1)+1e-04;
     # diag(currentval$invC) <- diag(currentval$invC) + 1e-04;
      eig  = spam::eigen.spam(as.spam(currentval$invC), TRUE); currentval$rtEV = sqrt(eig$values);currentval$EVC  = eig$vectors
      if (newl>oldl){
        theta_list<-rbind(as.matrix(theta_list),matrix(NA,nrow=newl-oldl,ncol=ncol(as.matrix(theta_list))))
        theta_list[,ncol(theta_list)]<-currentval$theta[1:newl]
      }else{
        
        theta_list[1:newl,ncol(theta_list)]<-currentval$theta[1:newl]
        if (newl<oldl){theta_list[(newl+1):oldl,ncol(as.matrix(theta_list))]<-NA}
      }
    }
    
    if (max(currentval$grid)>tmrca) {
      #break; #I need to shrink the grid, reduce theta and update calculations
      currentval$grid<-currentval$grid[1:(max(seq(1,length(currentval$grid))[currentval$grid<tmrca])+1)]
      intl = currentval$grid[2]-currentval$grid[1]; midpts = currentval$grid[-1]-intl/2
      #midpts<-((c(0,currentval$grid[-length(currentval$grid)])+currentval$grid)/2)[-1]
      currentval$invC <- Q_matrix(as.matrix(midpts), 0, 1)+1e-04;
      #diag(currentval$invC) <- diag(currentval$invC) + 1e-04;
      eig  = spam::eigen.spam(currentval$invC, TRUE); currentval$rtEV = sqrt(eig$values);currentval$EVC  = eig$vectors
      where<-ncol(as.matrix(theta_list))
      where2<-nrow(as.matrix(theta_list))
      if (where==1){
        if (where2>length(currentval$grid)){theta_list[(length(currentval$grid)):where2]<-NA}
      }else{
        if (where2>length(currentval$grid)){theta_list[(length(currentval$grid)):where2,where]<-NA}  
      }
    }
    
    #why do I need to always call this part? what is it doing? 
    coal_times=cumsum(currentval$times$t[,1])[currentval$times$t[,2]==0]
    currentval$lik_init = coal_lik_init(samp_times, n_sampled = n,coal_times, grid = currentval$grid)
    tus1<-U_split(currentval$theta,currentval$lik_init,currentval$invC,.01,.01)
    currentval$us1  = tus1$logpos
    currentval$dus1 = U_split(currentval$theta,currentval$lik_init,currentval$invC,.01,.01, TRUE)
    currentval$currentlik<--updatePos$loglik
    currentval$gaussprior<--tus1$logpri
    probs_list<-rbind(probs_list,c(-updatePos$loglik,currentval$prior_F,-updatePos$logpri,currentval$gaussprior))
  }else{#Ind=0, i.e. rejecting the update. 
    times_list<-cbind(times_list,cumsum(currentval$times$t[,1])[currentval$times$t[,2]==0])
    ##WHY do I update this if I am rejecting? I don't think it is necessary. 
    #currentval$currentlik<--res2$pos_summ$loglik
    #currentval$coalprior<--res2$pos_summ$logpri
    #CHECK!!!!
    probs_list<-rbind(probs_list,c(currentval$currentlik,currentval$prior_F,currentval$coalprior,currentval$gaussprior))
    
  }
  return(list(currentval=currentval,probs_list=probs_list,prec_list=prec_list,times_list=times_list,theta_list=theta_list, acp=Ind))
  
}


timesLocalProposal<-function(currentval,initial,a,sigma){
  #
  #Times update proposal according to the new scheme: Local choice of number of increments to change and truncated normal upates
  #
  TotConstraint=initial$max.nCoal+1
  n1=length(currentval$lik_init$args$coal_times)
  tchange=sample(a,1)
  change=sample(n1,tchange)
  coal_times_incre=cumsum(currentval$times$t[,1])[currentval$times$t[,2]==0]-c(0,cumsum(currentval$times$t[,1])[currentval$times$t[,2]==0][1:(n1-1)])
  samp_times1=samp_times[-1]
  muGauss=coal_times_incre[change]
  forDen<-c()
  revDen<-c()
  for (i in 1:length(muGauss)){
    constraints=TotConstraint[TotConstraint>=change[i]]
    idCon=which(TotConstraint>=change[i])
    lBound=0
    #if(length(constraints)>0){for (j in 1:length(constraints)){lBound=c(lBound,samp_times1-cumsum(coal_times_incre[-change[i]])[constraints[j]-1])}}
    #TO REDO
    lBound=max(0,-(cumsum(coal_times_incre)-muGauss[i])[constraints]+samp_times1[idCon])
    #lBound=abs(min(lBound))
    newtime=0
    while(newtime<=0.1*10^(-5)){newtime=rtruncnorm(1,a=lBound,b=Inf,mean=muGauss[i],sd=sigma*muGauss[i])}#ensuress I am not sampling sthg to small
    coal_times_incre[change[i]]=newtime
    forDen<-c(forDen,dtruncnorm(coal_times_incre[change[i]],a=lBound,b=Inf,mean=muGauss[i],sd=sigma*muGauss[i]))
    revDen<-c(revDen,dtruncnorm(muGauss[i],a=lBound,b=Inf,mean=coal_times_incre[change[i]],sd=sigma*coal_times_incre[change[i]]))
  }
  forDen=prod(forDen)
  revDen=prod(revDen)
  coal_time_prop=cumsum(coal_times_incre)
  propTimes=times_lik(coal_time_prop,samp_times)
  proposal=list(propTimes=propTimes,forDen=forDen,revDen=revDen)
  return(proposal)
}



localUpdate_RT_FREE_het<-function(result,n){
  #Implement local update algorithm suitably adapapted to ranked tree shape.
  #Note: It is uniform on the tree level, which may not be the best thing to do.   
  #Note: lots of constraints!!! to check irreducibility. 
  
  newresult=list()
  newresult$change_F=result$change_F
  newresult$family_size=result$family_size
  newresult$F_nodes=result$F_nodes
  newresult$node_group=result$node_group #For a given data it never changes. 
  
  
  #sample level (now uniform). Note: 1 means 1st vintage. 
  
  l<-sample((sum(n)-2),1)
  
  
  #transition probabilities first part
  newresult$ForwProb<-1/(sum(n)-2)
  newresult$BackProb<-1/(sum(n)-2)
  #order F_nodes
  result$F_nodes[[l]]<-sort(result$F_nodes[[l]])
  result$F_nodes[[l+1]]<-sort(result$F_nodes[[l+1]])
  
  #########################################
  ########## MARKOV  MOVES ############### 
  ########################################
  #read type A or B (FALSE MEANS B)
  
  
  #TO THINK ABOUT WHAT TO DO WITH TOTAL_VINTAGE AND INDICATOR IN CASE A!!!
  if (l %in% unlist(result$F_nodes[[(l+1)]])){ #Case A   
    a=1
    #sample move
    #1) do not change 2) 1 with 3 3)2 with 3  [note I consider 1 2 the initial cherry]
    m=sample(c(2,3),size=1) #Now stop it from doing sthg. 
    #transition probabilities.
    newresult$ForwProb<-newresult$ForwProb/2
    newresult$BackProb<-newresult$BackProb/2
    if (m==1){#do nothing
    } else{#move 2 and 3
      # Assumption: the vintaged are always ordered in incresing order: vintage l is always gonna be in [[l+1]][2]
      lineages=c(result$F_nodes[[l]][1:2],result$F_nodes[[l+1]][1]) 
      if (length(unique(lineages))==1){#singleton+cherry=nothing happens
      } else if (lineages[3]<=0){#joining to a singleton of some kind
        if (lineages[1]==lineages[2]) {#cherry joining a vintage or a singleton=it doesn't mate m=2 or 3
          newresult$ForwProb<-newresult$ForwProb*2
          newresult$F_nodes[[l]]<-c(lineages[1],lineages[3],l)
          newresult$F_nodes[[l+1]]<-c(lineages[2],l,l+1)
        } else if (lineages[1]==lineages[3]){#Note: they must be singletons of the same group
          if (m==2){#lineages[1] is involved=nothing happens
          } else {#m=3=it is not symmetric
            newresult$BackProb<-newresult$BackProb*2
            newresult$F_nodes[[l]]<-c(lineages[1],lineages[3],l)
            newresult$F_nodes[[l+1]]<-c(lineages[2],l,l+1)
            newresult$family_size[l]<-2 #l+1 position does not change
          }
        } else if (lineages[2]==lineages[3]){
          if (m==2){#not symmetric
            newresult$BackProb<-newresult$BackProb*2
            newresult$F_nodes[[l]]<-c(lineages[2],lineages[3],l)
            newresult$F_nodes[[l+1]]<-c(lineages[1],l,l+1)
            newresult$family_size[l]<-2 #l+1 position does not change
          } else {#m=3=they are equal=nothing happens
          }
        } else {#they are all different: lineages[3] is a singleton, the others can be singletons or vintages
          if (m==2){
            newresult$F_nodes[[l]]<-c(lineages[2],lineages[3],l)
            newresult$F_nodes[[l+1]]<-c(lineages[1],l,l+1)
            if (lineages[1]>0){
              newresult$family_size[l]<-result$family_size[l+1]-result$family_size[lineages[1]]+1 #from first tree
            }#if lineages[1] is a singleton, I am swapping with a singleton, so np.
          } else {#m=3
            newresult$F_nodes[[l]]<-c(lineages[1],lineages[3],l)
            newresult$F_nodes[[l+1]]<-c(lineages[2],l,l+1)
            if (lineages[2]>0){
              newresult$family_size[l]<-result$family_size[l+1]-result$family_size[lineages[2]]+1 #from first tree
            }#if lineages[1] is a singleton, I am swapping with a singleton, so np.
          }
        }
      } else if (lineages[3]>0){#joining a vintage in l+1
        if (lineages[1]==lineages[2]) {#cherry joining a vintage or a singleton=it doesn't mate m=2 or 3
          #It is not symmetric
          newresult$ForwProb<-newresult$ForwProb*2
          newresult$F_nodes[[l]]<-c(lineages[1],lineages[3],l)
          newresult$F_nodes[[l+1]]<-c(lineages[2],l,l+1)
          newresult$family_size[l]<-result$family_size[l+1]-1
        } else {#they must be all distinct (lineages[3] is unique)
          if (m==2){
            newresult$F_nodes[[l]]<-c(lineages[2],lineages[3],l)
            newresult$F_nodes[[l+1]]<-c(lineages[1],l,l+1)
            if (lineages[1]>0){
              newresult$family_size[l]<-result$family_size[l+1]-result$family_size[lineages[1]]+result$family_size[lineages[3]] 
            } else{
              newresult$family_size[l]<-result$family_size[l+1]-1
            }
          } else {#m=3
            newresult$F_nodes[[l]]<-c(lineages[1],lineages[3],l)
            newresult$F_nodes[[l+1]]<-c(lineages[2],l,l+1)
            if (lineages[2]>0){
              newresult$family_size[l]<-result$family_size[l+1]-result$family_size[lineages[2]]+result$family_size[lineages[3]] 
            } else{
              newresult$family_size[l]<-result$family_size[l+1]-1
            }
          }
        }
      }
    }    
  }  else { #Case B
    a=0
    newresult$family_size[c(l,l+1)]=result$family_size[c(l+1,l)]#invert the position
    newresult$F_nodes[[l]][1:2]=result$F_nodes[[l+1]][1:2]
    newresult$F_nodes[[l+1]][1:2]=result$F_nodes[[l]][1:2]
    #fix the result of F_node: what's l has to become l+1 and viceversa (not necessarily for Case A?? check)
    l_epoch=max(which(unlist(result$F_nodes)==l))%/%3+1
    l_epoch_position=max(which(unlist(result$F_nodes)==l))%%3
    newresult$F_nodes[[l_epoch]][l_epoch_position]=l+1
    newresult$F_nodes[[l_epoch]]=sort(newresult$F_nodes[[l_epoch]])
    l1_epoch=max(which(unlist(result$F_nodes)==l+1))%/%3+1
    l1_epoch_position=max(which(unlist(result$F_nodes)==l+1))%%3
    newresult$F_nodes[[l1_epoch]][l1_epoch_position]=l
    newresult$F_nodes[[l1_epoch]]=sort(newresult$F_nodes[[l1_epoch]])
  }
  newresult$l<-l
  newresult$a<-a
  
  
  return(newresult)
}

#' Moves topology Cappello et al. (2020)
#' 
#' @param currentval current val of the chain
#' @param probs_list array recording probability 
#' @param const scaling factor
#' @param mu mutation parameter
#' @param oldsuff perfect phylogeny
#' @param n vector of samples
#'   
#' @return MC steps for topology
#' @export


updateFmat_Markov_FREE_het<-function(currentval,probs_list,const=1,mu,oldsuff,n){
  result<-currentval$result
  result_new <- localUpdate_RT_FREE_het(result,n)
  newprior<-tajimaPrior_het(result_new$F_nodes,n,currentval$times$addtimes)
  t<-currentval$times$t[,1]
  indt<-currentval$times$t[,2]
  oldsuff_het<-pp_heterochronous(n_sampled=n,oldsuff,currentval$times$addtimes)
  likH<-python.call("calcPF_het",result_new$F_nodes, result_new$family_size,oldsuff_het,result_new$node_group,t*mu,indt,"True")
  newlike<-likH[[2]]/const
  if (!is.null(likH[[2]])){
    AR<--currentval$currentlik+newlike+log(newprior)-log(currentval$prior_F)-log(result_new$ForwProb)+log(result_new$BackProb)
    if (runif(1)<exp(AR)){
      acp<-1
      currentval$result<-result_new
      #currentval$tus2=currentval$Ufun2(currentval$times,nsites,currentval$result,currentval$grid,currentval$theta,logliktot=0,TRUE)
      #currentval$us2=currentval$tus2$logpos
      #currentval$dus2=currentval$tus2$dlogpos
      currentval$currentlik<-newlike
      #currentval$lik_call<-where_save
      currentval$prior_F<-newprior
      #currentval$proposal_F<-newproposal
    } else{
      acp<-0
    }
  } else{
    acp<-0
  }
  #F_list[[length(F_list)+1]]<- currentval$result$F
  probs_list<-rbind(probs_list,c(currentval$currentlik,currentval$prior_F,currentval$coalprior,currentval$gaussprior))
  return(list(currentval=currentval,probs_list=probs_list, acp=acp))
}


#######################################################
################  Support functions ###################
#######################################################

##Code for Tajima-based inference from a single locus
#All the working functions are updated here
#Update December 2017
beastfile_het<-function(data1,fname,samp_times,n){
  ##For the Wald MtDNA dataset
  nn<-n
  s<-nrow(data1)
  n<-nchar(data1[1,1])
  #  set.seed(1234)
  ancestral<-rmultinom(s,1,p=c(1,1,1,1))
  
  anc_allele<-rep("A",s)
  der_allele<-rep("A",s)
  for (j in 1:s){
    anc_allele[j]<-c("A","T","C","G")[ancestral[,j]==1]
    der_allele[j]<-c("A","T","C","G")[rmultinom(1,1,p=1-ancestral[,j])==1]
  }
  
  tabledata<-matrix("0",nrow=n,ncol=s)
  for (j in 1:s){
    x<-strsplit(data1[j,],NULL)[[1]]
    for (i in 1:length(x)){
      x[i]<-ifelse(x[i]=="0",anc_allele[j],der_allele[j])
    }
    tabledata[,j]<-x
  }
  la<-seq(1,n)
  lb<-rep(samp_times,nn)
  labels<-paste(la,"_",lb,sep="")
  rownames(tabledata)<-labels
  #generate fasta file
  write.dna(tabledata,file=fname,format="fasta")
}


# Intrinsic precision matrix
Q_matrix <- function(input, s_noise = 0, signal = 1)
{
  n2 <- length(input)
  diff1 <- diff(input)
  diff1[diff1==0] <- s_noise #correction for dividing over 0
  diff <- (1/(signal*diff1))
  
  Q<-spam::spam(0,n2,n2)  
  if (n2>2)
  {
    Q[cbind(seq(1,n2),seq(1,n2))] <- c(diff[1], diff[1:(n2-2)] + diff[2:(n2-1)],
                                       diff[n2-1]) + (1/signal)*rep(s_noise, n2)
  }
  else
  {
    Q[cbind(seq(1,n2),seq(1,n2))] <- c(diff[1],diff[n2-1])+(1/signal)*rep(s_noise,n2)
  }
  Q[cbind(seq(1,n2-1),seq(2,n2))] <- -diff[1:(n2-1)]
  Q[cbind(seq(2,n2),seq(1,n2-1))] <- -diff[1:(n2-1)]
  
  return(Q)
}

# backwards compatibility (deprecate soon)
Q.matrix = function(...)
{
  return(Q_matrix(...))
}



#' Subsample
#' 
#' @param res realization of the chain
#' @param burnin burn of the chain 
#' @param subsample subsampling of the MC chain
#'   
#' @return approximate sample from the posterior
#' @export

subsample<-function(res,burnin=0,subsample=1){
  nsamp=dim(res)[2]
  indices=seq(burnin+1,nsamp,by=subsample)
  pos_res=res[,indices]
  return(pos_res)
}



envelope<-function(INLA_out, traj, hilim=0.07, lolim=0, yhilim=Inf, ylolim=0)
{
  #mod = INLA_out$result$summary.random$time
  Infun<-stepfun(INLA_out[,1],c(INLA_out[1,2],INLA_out[,2]))
  InfunLo<-stepfun(INLA_out[,1],c(INLA_out[1,4],INLA_out[,4]))
  InfunHi<-stepfun(INLA_out[,1],c(INLA_out[1,3],INLA_out[,3]))
  grid_pts = seq(0,hilim,length.out = 100)
  n = length(grid_pts)
  lo = InfunLo(grid_pts)
  hi = InfunHi(grid_pts)
  truth = traj(grid_pts)
  result = sum(truth < hi & truth > lo)
  return(list(tot = result, avg = result / n))
}

dev = function(INLA_out, traj, hilim=Inf, lolim=0, yhilim=Inf, ylolim=0)
{
  Infun<-stepfun(INLA_out[,1],c(INLA_out[1,2],INLA_out[,2]))
  InfunLo<-stepfun(INLA_out[,1],c(INLA_out[1,4],INLA_out[,4]))
  InfunHi<-stepfun(INLA_out[,1],c(INLA_out[1,3],INLA_out[,3]))
  grid_pts = seq(0,hilim,length.out = 100)
  n = length(grid_pts)
  lo = InfunLo(grid_pts)
  hi = InfunHi(grid_pts)
  truth = traj(grid_pts)
  med = Infun(grid_pts)
  result = sum( abs(med - truth)/truth )
  return(list(tot = result, avg = result / n ))
}

relwid = function(INLA_out, traj, hilim=Inf, lolim=0, yhilim=Inf, ylolim=0)
{
  Infun<-stepfun(INLA_out[,1],c(INLA_out[1,2],INLA_out[,2]))
  InfunLo<-stepfun(INLA_out[,1],c(INLA_out[1,4],INLA_out[,4]))
  InfunHi<-stepfun(INLA_out[,1],c(INLA_out[1,3],INLA_out[,3]))
  grid_pts = seq(0,hilim,length.out = 100)
  n = length(grid_pts)
  lo = InfunLo(grid_pts)
  hi = InfunHi(grid_pts)
  truth = traj(grid_pts)
  result = sum( (hi - lo)/truth )
  return(list(tot = result, avg = result / n ))
}


