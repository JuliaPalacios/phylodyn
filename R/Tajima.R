##Code for Tajima-based inference from a single locus
#All the working functions are updated here
#Update December 2017
beastfile<-function(data1,fname){
  ##For the Wald MtDNA dataset
  
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
  
  #generate fasta file
  write.dna(tabledata,file=fname,format="fasta")
  
}

differentF<-function(res1,res2){
  #update May 2018
  #Function to summarize posterior distribution of F alone
  sample1<-res1[[1]]
  Flist<-list(sample1)
  tList<-list(res2[,1])
  tlist2<-list(res2[,1])
  freqF<-1
  for (j in 2:(length(res1))){
    sample1<-res1[[j]]
    exists<-0
    for (i in 1:length(Flist)){
      if ((sum(abs(Flist[[i]]-sample1)>0)!=0)) { #This F is new
        exists<-exists+1
      }else{
        freqF[i]<-freqF[i]+1;
        tList[[i]]<-tList[[i]]+res2[,j]
        tlist2[[i]]<-cbind(tlist2[[i]],res2[,j])
      }
      if (exists==length(Flist)){
        Flist[[length(Flist)+1]]<-sample1
        freqF<-c(freqF,1)
        tList[[length(tList)+1]]<-res2[,j]
        tlist2[[length(tlist2)+1]]<-res2[,j]
      }
    }
  }
  return(list(freqF=freqF,tList=tList,Flist=Flist,tlist2=tlist2))
}

sufficient_stats<-function(data){
  #Update Dec 2017
  #This is the old function needed for 
  sorted<-sort(data[,1],index.return=TRUE)
  n<-nchar(data[1,])
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
  Index<- O%*%diag(seq(1,l))
  #(3)For each cell
  for (j in 1:n){
    haplotypes<-c(haplotypes,paste0(O[j,],collapse=""))
    for (i in 2:l){
      if (O[j,i]==1){Lmat[j,i]<-max(Index[j,1:(i-1)])}
    }
  }
  #correct for multiple haplotypes
  sort.hap<-sort(haplotypes[-1],index.return=T)
  base<-sort.hap$x[1]
  base.ind<-sort.hap$ix[1]
  remove<-0
  for (j in 2:n){
    if (base==sort.hap$x[j]){
      remove<-c(remove,sort.hap$ix[j])
      frequency[base.ind]<-frequency[base.ind]+1
      frequency[sort.hap$ix[j]]<-0
      
    }else{
      base<-sort.hap$x[j]
      base.ind<-sort.hap$ix[j]
    }
  }
  leaves<-apply(Index,1,max)
  card<-rev(groups$cardinality)
  carriers<-rev(groups$carriers) #I don't use this
  y<-0
  mylist<-list(list(x=0,y=0))
  orderlist<-0
  parentlist<-0
  famsize<-0
  numsize<-0
  L<-apply(Lmat,2,max) #this vector has the nesting information, it has parents nodes
  parents<-sort(unique(L),d=TRUE)
  i<-2
  for (j in parents){
    offspring<-seq(1,l)[L==j]
    offspringsize<-0
    indic_offspring<-rep(0,length(offspring))
    for (no in offspring){
      offspringsize<-c(offspringsize,sum(leaves==no))
      if (sum(parents==no)>0) {indic_offspring[offspring==no]<-1}
    }
    offspringsize<-offspringsize[-1]
    #offspringsize<-frequency[offspring]
    offspringsizelist<-unique(offspringsize) ### I merge nodes with the same size
    #this should only be allowed for leaf nodes.
    ##Per-size, see if there are parents there
    for (k in offspringsizelist){
      #  if  (sum(indic_offspring[offspringsize==k])==0){
      ind2<-indic_offspring[offspringsize==k]
      #those where ind2 is 0
      if (sum(ind2==0)>0){
        mylist[[i]]<-list(x=sum(offspringsize==k & indic_offspring==0),y=card[offspring[offspringsize==k & indic_offspring==0]]) 
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
          mylist[[i]]<-list(x=1,y=card[offspring[w]]) 
          famsize<-c(famsize,k)
          numsize<-c(numsize,1)
          y<-c(y,0)
          parentlist<-c(parentlist,j)
          orderlist<-c(orderlist,offspring[w])
          leaves[leaves==offspring[w]]<-j
          i<-i+1
        }
      }
      
    }
    
  }
  
  mylist[[1]]<-NULL
  suff<-list(mylist=mylist,nodes=cbind(orderlist[-1],parentlist[-1],famsize[-1],numsize[-1],y[-1]))
  ##I need to correct for when the root doesn't add up the total of nodes
  tot2<-cbind(suff$nodes[,2],suff$nodes[,3]*suff$nodes[,4])
  tot<-sum(suff$nodes[suff$nodes[,2]==0,3]*suff$nodes[suff$nodes[,2]==0,4])
  if (tot<n){
    #if there are singletons descending from 0, add it there if not, add a singleton
    if (sum(suff$nodes[suff$nodes[,2]==0 & suff$nodes[,3]==1,1])>0){
      where<-seq(1,nrow(suff$nodes))[suff$nodes[,2]==0 & suff$nodes[,3]==1]
      suff$nodes[where,4]<-suff$nodes[where,4]+n-tot
      suff$mylist[[where]]$x<-suff$mylist[[where]]$x+n-tot
      suff$mylist[[where]]$y<-c(suff$mylist[[where]]$y,rep(0,n-tot))
    }else{
      #I am not sure I actually need this, actually yes! It is needed
       suff$nodes<-rbind(suff$nodes,c(max(suff$nodes[,1])+1,0,1,n-tot,1))
       where<-nrow(suff$nodes)
       suff$mylist[[where]]<-list(x=n-tot,y=rep(0,n-tot))
    }
  }
  ##I need to make sure everything adds up. If something is off, I need
  #to add 0 mutations
  parents<-parents[parents>0]
  for (j in parents){
    toadd<-suff$nodes[suff$nodes[,1]==j,3]*suff$nodes[suff$nodes[,1]==j,4]-sum(tot2[tot2[,1]==j,2])
    if (toadd>0){
      #just adds a duplicate
      if (sum(suff$nodes[suff$nodes[,2]==j & suff$nodes[,3]==1,1])>0){
        where<-seq(1,nrow(suff$nodes))[suff$nodes[,2]==j & suff$nodes[,3]==1]
        suff$nodes[where,4]<-suff$nodes[where,4]+toadd
        suff$mylist[[where]]$x<-suff$mylist[[where]]$x+toadd
        suff$mylist[[where]]$y<-c(suff$mylist[[where]]$y,rep(0,toadd))
      }else{
        where<-nrow(suff$nodes)+1
        who<-max(suff$nodes[,1:2])+1
        suff$mylist[[where]]<-list(x=toadd,y=rep(0,toadd))
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

sufficient_stats2<-function(data){
  #Update March 15, 2017
  ##This is a modified version that does not group nodes with the same size
  #data is a matrix of m rows (snps) and 1 column of concateted 0s and 1s. The number of characters in each row is n, the number of individuals
  sorted<-sort(data[,1],index.return=T)
  n<-nchar(data[1,])
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
  Index<- O%*%diag(seq(1,l))
  #(3)For each cell
  for (j in 1:n){
    haplotypes<-c(haplotypes,paste0(O[j,],collapse=""))
    for (i in 2:l){
      if (O[j,i]==1){Lmat[j,i]<-max(Index[j,1:(i-1)])}
    }
  }
  #correct for multiple haplotypes
  sort.hap<-sort(haplotypes[-1],index.return=T)
  base<-sort.hap$x[1]
  base.ind<-sort.hap$ix[1]
  remove<-0
  for (j in 2:n){
    if (base==sort.hap$x[j]){
      remove<-c(remove,sort.hap$ix[j])
      frequency[base.ind]<-frequency[base.ind]+1
      frequency[sort.hap$ix[j]]<-0
      
    }else{
      base<-sort.hap$x[j]
      base.ind<-sort.hap$ix[j]
    }
  }
  leaves<-apply(Index,1,max)
  card<-rev(groups$cardinality)
  carriers<-rev(groups$carriers) #I don't use this
  y<-0
  orderlist<-0
  parentlist<-0
  famsize<-0
  L<-apply(Lmat,2,max) #this vector has the nesting information, it has parents nodes
  parents <- sort(unique(L), decreasing = TRUE)
  i<-2
  for (j in parents){
    offspring<-seq(1,l)[L==j]
    offspringsize<-0
    for (no in offspring){
      # offspringsize<-c(offspringsize,sum(leaves==no))
      y<-c(y,card[no])
      famsize<-c(famsize,sum(leaves==no))
      parentlist<-c(parentlist,j)
      orderlist<-c(orderlist,no)
      leaves[leaves==no]<-j
      i<-i+1
    }
  }
  nodes=cbind(orderlist[-1],parentlist[-1],famsize[-1],y[-1])
  temp<-sufficient_remove(nodes)
  temp2<-correct_sufficient(temp)

  #  mylist[[1]]<-NULL
  return(list(nodes=temp2,n=n))
}

sufficient_remove<-function(nodes){
  #Update March 14, 2017
  #This function groups leaf nodes with the same size
  
  parents<-unique(nodes[,2])
  total<-nrow(nodes)
  #leaves<-seq(1,total)[-parents[parents>0]]
  leaf<-rep(1,total)
  for (j in parents){
    leaf[nodes[,1]==j]<-0
  }
  nodes<-cbind(nodes,leaf)
  subtable<-nodes[leaf==1,]
  subtable<-cbind(subtable,paste(subtable[,2],"-",subtable[,3],sep=""))
  parents<-as.data.frame(table(subtable[,6]))
  parents<-parents[parents[,2]>1,]
  ##This function removes repeats
  for (j in parents[,1]){
    pa<-as.numeric(unique(subtable[subtable[,6]==j,2]))
    si<-as.numeric(unique(subtable[subtable[,6]==j,3]))
    nodes[nodes[,2]==pa & nodes[,3]==si,1]<-min(nodes[nodes[,2]==pa & nodes[,3]==si,1])
  }
  return(nodes)
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



sufficient_extend<-function(oldsuff){
  rep_idx<-rep(seq(1,nrow(oldsuff$nodes)),oldsuff$nodes[,4])
  nodes_extended<-oldsuff$nodes[rep_idx,]
  return(nodes_extended)
}




bring_branch_lengths<-function(u,F){
  #u is the vector of intercoalescent times
  #given u and F, returns branch lengths ln+1,ln,ln-1,...,l3 and family sizes
  dimf<-nrow(F)
  diffM<-F[2:(dimf),2:dimf]-F[2:dimf,1:(dimf-1)]
  singletons<-seq(1,dimf)[diff(F[,1])<0]
  d<-rep(u[1],dimf) #for corresponding ln+1,ln,ln-1,...,l3
  firstzero<-rep(0,dimf-2) 
  familysize<-rep(2,dimf)#for corresponding Sn,Sn-1,..S2
  clades_offspring<-0
  clades_parents<-0
  coal_times<-cumsum(u)
  for (j in 1:(dimf-2)){
    condition<-diffM[j:(dimf-1),j]==0
    if (sum(condition)>0){ firstzero[j]<-min(seq(j,(dimf-1))[condition])}else{firstzero[j]<-dimf}
    
    d[j+1]<-coal_times[firstzero[j]]-coal_times[j]
  }
  
  d[dimf]<-u[dimf]
  firstzerotrans<-dimf+2-firstzero
  for (i in 1:(dimf-2)){
    count<-sum(firstzerotrans==(dimf-i+1))
    if (count==0){
      familysize[i+1]<-2
    }
    if (count==1){
      familysize[i+1]<-familysize[seq(1,dimf-2)[firstzerotrans==(dimf-i+1)]]+1
      clades_offspring<-c(clades_offspring,seq(1,dimf-2)[firstzerotrans==(dimf-i+1)])
      clades_parents<-c(clades_parents,i+1)
      
    }
    if (count==2){
      familysize[i+1]<-familysize[min(seq(1,(dimf-2))[firstzerotrans==(dimf-i+1)])]+familysize[max(seq(1,(dimf-2))[firstzerotrans==(dimf-i+1)])]
      clades_offspring<-c(clades_offspring,min(seq(1,(dimf-2))[firstzerotrans==(dimf-i+1)]),max(seq(1,(dimf-2))[firstzerotrans==(dimf-i+1)]))
      clades_parents<-c(clades_parents,i+1,i+1)
    }
    
  }
  familysize<-familysize[-dimf]
  familysize<-c(1,familysize)
  return(list(d=d,familysize=familysize,nodes=cbind(clades_offspring[-1],clades_parents[-1]),singletons=singletons))
}

##Generate F matrices with constrained number of cherries
generateFconst<-function(n,cherries){
  ##According to Genetics paper, Oct 2015
  F<-matrix(0,nrow=n-1,ncol=n-1)
  diag(F)<-seq(n,2)
  F[2,1]<-n-2
  current.block.size<-c(2)
  vintage.size<-rep(1,n-1)
  vintage<-n
  c<-1
  for (j in 3:(n-1)){
    if (c>=cherries){
      ran<-sum(rmultinom(1,1,p=c(0,F[j-1,1]*(n+2-j-F[j-1,1]),choose(n+1-j-F[j-1,1],2)))*c(2,1,0))
    }else{
      ran<-sum(rmultinom(1,1,p=c(choose(F[j-1,1],2),F[j-1,1]*(n+2-j-F[j-1,1]),choose(n+1-j-F[j-1,1],2)))*c(2,1,0))
    }
    
    F[j,1]<-F[j-1,1]-ran
    if (ran==2){c<-c+1;vintage<-c(vintage,n+2-j);vintage.size[j-1]<-length(vintage);F[j,2:(j-1)]<-F[j-1,2:(j-1)]-2; current.block.size<-c(current.block.size,2)}
    if (ran==1){
      who<-sample(1:vintage.size[j-2],1)
      vintage.who<-vintage[who]
      F[j,2:(n+1-vintage.who)]<-F[j-1,2:(n+1-vintage.who)]-1
      if ((n+2-vintage.who)<=(j-1)){
        F[j,(n+2-vintage.who):(j-1)]<-F[j-1,(n+2-vintage.who):(j-1)]-2}
      vintage<-c(vintage[-who],n+2-j)
      vintage.size[j-1]<-vintage.size[j-2]
    }
    if (ran==0){
      who<-sample(1:vintage.size[j-2],2)
      vintagewho<-vintage[who]
      whomin<-min(vintagewho)
      whomax<-max(vintagewho)
      F[j,2:(n+1-whomax)]<-F[j-1,2:(n+1-whomax)]
      if ((n+2-whomax)<=(n+1-whomin)){
        F[j,(n+2-whomax):(n+1-whomin)]<-F[j-1,(n+2-whomax):(n+1-whomin)]-1}
      if ((n+2-whomin)<=(j-1)){
        F[j,(n+2-whomin):(j-1)]<-F[j-1,(n+2-whomin):(j-1)]-2}
      vintage<-c(vintage[-who],n+2-j)
      vintage.size[j-1]<-vintage.size[j-2]-1
    }
  }
  return(F)
}

correct_sufficient<-function(sufficient){
  parents<-unique(sufficient[,2])
  parents<-parents[parents>0]
  for (j in parents){
    expected<-sufficient[sufficient[,1]==j,3]
    observed<-sum(sufficient[sufficient[,2]==j,3])
    difference<-expected-observed
    if (difference>0){
      label<-max(sufficient[,1])+1
      sufficient<-rbind(sufficient,c(label,j,1,0,1))
      if (difference>1){
        for (l in 2:difference){
        sufficient<-rbind(sufficient,c(label,j,1,0,1))
      }}
      #I either add one line per node or not do anything --need to check this
    }
  }
return(sufficient)
}

##Generate F proposals

F_calc<-function(PP,ftotal){
  ##Correction: If size==1, node label should be the parent 
  PP[PP[,3]==1,1]<-PP[PP[,3]==1,2]
  free<-PP[,5]*PP[,3] ##Takes leaves and number of offspring
  potential<-unique(PP[free>=1,1])
  pot<-0
  tot<-0
  fre<-0
  for (j in potential){
    total<-sum(free[PP[,1]==j])
    tot2<-sum(ftotal[PP[,1]==j])
    if (total>=2){pot<-c(pot,j); tot<-c(tot,total);fre<-c(fre,tot2)}
  }
  return(list(pot=pot[-1],tot=tot[-1],fre=fre[-1],PP=PP))
}


F_iter<-function(PP){
  ##Update: Nov 29, 2017
  ##Correction: If size==1, node label should be the parent 
  ##New correction, if parent=0,leaf=1, node label should the parent
  PP[PP[,2]==0 & PP[,5]==1,1]<-0
  PP[PP[,3]==1,1]<-PP[PP[,3]==1,2]
  free<-PP[,5]*PP[,3] ##Takes leaves and number of offspring
  potential<-unique(PP[free>=1,1])
  pot<-0
  tot<-0
  for (j in potential){
    total<-sum(free[PP[,1]==j])
    if (total>=2){pot<-c(pot,j); tot<-c(tot,total)}
  }
  return(list(pot=pot[-1],tot=tot[-1],PP=PP))
}


# F_iter<-function(PP){
#   ##Update: Nov 29, 2017
#   ##Correction: If size==1, node label should be the parent 
#   ##New correction, if parent=0,leaf=1, node label should the parent
#   PP[PP[,2]==0 & PP[,5]==1,1]<-0
#   PP[PP[,3]==1,1]<-PP[PP[,3]==1,2]
#   free<-PP[,5]*PP[,3]*PP[,4] ##Takes leaves and number of offspring
#   potential<-unique(PP[free>=1,1])
#   pot<-0
#   tot<-0
#   for (j in potential){
#     total<-sum(free[PP[,1]==j])
#     if (total>=2){pot<-c(pot,j); tot<-c(tot,total)}
#   }
#   return(list(pot=pot[-1],tot=tot[-1],PP=PP))
# }
  

familysize<-function(PP,vintage){
  tot<-PP[,3]*PP[,5]
  res<-sum(tot[PP[,1]==vintage])
  res2<-sum(PP[PP[,1]>vintage & PP[,2]==vintage,3])
  return(res+res2)
}


check<-function(PP,sufficient,vintage){
  ##If there is one, I should not eliminate it
  where<-seq(1,nrow(PP))[PP[,1]==vintage & PP[,3]==1 & PP[,5]==1]
  famsize<-familysize(PP,vintage)
  
  if (famsize==1){
    parent<-min(PP[PP[,1]==vintage,2])
    #the number of active children from parent is also 1
  #  original_remove<-sum(sufficient[PP[,1]==vintage & PP[,2]>=vintage,3])
    where2<-seq(1,nrow(PP))[PP[,1]==vintage & PP[,2]==parent]
    
    if (where!=where2){
 #   PP[PP[,1]==vintage & PP[,2]==parent,3]<-sufficient[PP[,1]==vintage & PP[,2]==parent,3]-original_remove+1
      if (PP[where2,5]==1){
        PP[where2,3]<- PP[where2,3]+1
      }else{
        PP[where2,3]<-1
        PP[where2,5]<-1
      }
      PP[where,5]<-0
    }
  }
  return(PP)
}



check2<-function(PP,sufficient,vintage,vintageRecord,labelvintage){
  ##If there is one, I should not eliminate it
  where<-seq(1,nrow(PP))[PP[,1]==vintage & PP[,3]==1 & PP[,5]==1]
  famsize<-familysize(PP,vintage)
  
  if (famsize==1){
    parent<-min(PP[PP[,1]==vintage,2])
    #the number of active children from parent is also 1
    #  original_remove<-sum(sufficient[PP[,1]==vintage & PP[,2]>=vintage,3])
    where2<-seq(1,nrow(PP))[PP[,1]==vintage & PP[,2]==parent]
    
    if (where!=where2){
      #   PP[PP[,1]==vintage & PP[,2]==parent,3]<-sufficient[PP[,1]==vintage & PP[,2]==parent,3]-original_remove+1
      if (PP[where2,5]==1){
        PP[where2,3]<- PP[where2,3]+1
      }else{
        PP[where2,3]<-1
        PP[where2,5]<-1
      }
      PP[where,5]<-0
    }
    current_record<-sort(unique(c(vintageRecord[labelvintage==parent,],vintageRecord[labelvintage==vintage,])))
    diffcero<-sum(current_record>0)
    vintageRecord[labelvintage==parent,1:diffcero]<-current_record[current_record>0]
  }
  return(list(PP=PP,vintageRecord=vintageRecord))
}
F_iter2<-function(sufficient){ 
  #Update April 27, 2017
 PP<-sufficient
 total_free<-PP[,3]*PP[,5]
 ##In sufficient_back I will retain the original data to identify original free
  total_vintage<-0
  prob<-0
  change_F<-0
  n<-sum(sufficient[sufficient[,2]==0,3]) #descending from root node 0
  
  while(length(total_vintage)<n){
     potential<-F_iter(PP)
  
     if (length(potential$tot)>1){
      potential_size2<-choose(potential$tot,2)
      choose<-rmultinom(1,1,potential_size2)
      vintage<-potential$pot[choose==1]
      prob<-c(prob,potential_size2[choose==1]/sum(potential_size2))
    }else{
      vintage<-potential$pot
      prob<-c(prob,1)
    }
    ##lineages to choose from for the coalescence
    available<-potential$PP[potential$PP[,1]==vintage,3]*potential$PP[potential$PP[,1]==vintage,5]
     ##this are lineages available, some vintaged and some free
    ava_free<-total_free[potential$PP[,1]==vintage]
    tochange<-sum(available)
 
    if (tochange>1){
      free<-sum(ava_free)
      counts_Free<-c(choose(free,2),free*(tochange-free),choose(tochange-free,2))
      therandom<-rmultinom(1,1,counts_Free)
      prob<-c(prob,counts_Free[therandom==1]/sum(counts_Free))
      change_F<-c(change_F,c(2,1,0)[therandom==1])
    }else{
      #do something to define therandom
    }
      if (therandom[1]==1){
        #I will coalesce two free, it doesn't matter which two
        if (sum(ava_free>1)>=1){
          ##Just pick the one with more than one and substract two and also just substract two from the family
          total_free[potential$PP[,1]==vintage & total_free>1][1]<- total_free[potential$PP[,1]==vintage & total_free>1][1]-2
          PP<-potential$PP
          PP[PP[,1]==vintage & PP[,5]==1,3][1]<-PP[PP[,1]==vintage & PP[,5]==1 & PP[,3]>0,3][1]-1
          if (PP[PP[,1]==vintage & PP[,5]==1,3][1]==0){
            PP[PP[,1]==vintage & PP[,5]==1,5][1]<-0
          }
        }else{
          #eliminate two free
          #just pick the first two buckets and substract 1
          where<-min(seq(1,length(total_free))[potential$PP[,1]==vintage & total_free>0])
          total_free[potential$PP[,1]==vintage & total_free>0]<-total_free[potential$PP[,1]==vintage & total_free>0]-c(1,1,rep(0,sum(ava_free)-2))
          #eliminate just one lineage. I need to take one of the leaves and make it 0 and also size equals 0
          PP<-potential$PP
          PP[where,3]<-PP[where,3]-1
          #I make it not leaf if indead this disapears
          if (PP[where,3]==0){
              PP[where,5]<-0
          }
        
        }
      }
      if (therandom[2]==1){
        #there is one free and one vintaged
        #Eliminate one free
        where<-min(seq(1,length(total_free))[potential$PP[,1]==vintage & total_free>0])
        
        total_free[where]<-total_free[where]-1
        ##From PP, I eliminate one "free" lineage
        PP<-potential$PP
        PP[where,3]<-PP[where,3]-1
        if (PP[where,3]==0){
          PP[where,5]<-0
        }
    }
    if (therandom[3]==1){
      #just eliminate one that is not free from PP
      where<-min(seq(1,length(total_free))[potential$PP[,1]==vintage & potential$PP[,3]>0 & total_free==0])
      PP<-potential$PP
      PP[where,3]<-PP[where,3]-1
      if (PP[where,3]==0){
        PP[where,5]<-0
      }
    }
  
    if (sum(PP[,5])>1){PP<-check(PP,sufficient,vintage)}
    total_vintage<-c(total_vintage,vintage)
  }
  return(list(total_vintage=total_vintage[-1],change_F=change_F[-1],prob=prob[-1]))
}

consolidateF<-function(Fju){
  ##generates an F matrix consistent with paper notation
  newF<-matrix(0,nrow=nrow(Fju),ncol=nrow(Fju))
  for (j in 1:nrow(newF)){
    newF[nrow(Fju)-j+1,]<-rev(Fju[,j])
  }
  newF2<-matrix(0,nrow=Fju[1,1],ncol=Fju[1,1])
  newF2[2:Fju[1,1],2:Fju[1,1]]<-newF
  return(newF2)
}



F_sample<-function(suf_ext){ 
  #Update November 29, 2017. This function samples F and allocation given the PP like F_iter3 but uniformly
  n<-sum(suf_ext[suf_ext[,2]==0,3])
  #n<-sufficient_tot$n
  # n<-sum(sufficient[sufficient[,2]==0,3]) #descending from root node 0
  F<-matrix(0,nrow=n-1,ncol=n-1)
  diag(F)<-seq(n,2)
  labelvintage<-sort(unique(c(suf_ext[,2],suf_ext[,1])))
  vintageRecord<-matrix(0,nrow=length(labelvintage),ncol=n-1)
  PP<-suf_ext
  #1 if it is a leaf
  total_free<-PP[,3]*PP[,5]
  ##In sufficient_back I will retain the original data to identify original free
  total_vintage<-0
  prob<-0
  change_F<-0
  countIndex<-2
  
  while(length(total_vintage)<(n-1)){
    potential<-F_iter(PP)
    
    if (length(potential$tot)>1){
      potential_size2<-choose(potential$tot,2)
      choose<-rmultinom(1,1,rep(1,length(potential_size2)))
      #choose<-rmultinom(1,1,potential_size2)
      vintage<-potential$pot[choose==1]
      prob<-c(prob,1/length(potential_size2))
      # prob<-c(prob,sum(potential_size2*choose)/sum(potential_size2))
    }else{
      vintage<-potential$pot
      prob<-c(prob,1)
    }
    ##lineages to choose from for the coalescence
    available<-potential$PP[potential$PP[,1]==vintage,3]*potential$PP[potential$PP[,1]==vintage,5]
    ##this are lineages available, some vintaged and some free
    ava_free<-total_free[potential$PP[,1]==vintage]
    tochange<-sum(available)
    
    free<-sum(ava_free)
    counts_Free<-c(choose(free,2),free*(tochange-free),choose(tochange-free,2))
    therandom<-rmultinom(1,1,counts_Free)
    prob<-c(prob,counts_Free[therandom==1]/sum(counts_Free))
    change_F<-c(change_F,c(2,1,0)[therandom==1])
    # print(therandom)
    if (therandom[1]==1){
      F[countIndex,1:(countIndex-1)]<-F[countIndex-1,1:(countIndex-1)]-2
      countIndex<-countIndex+1
      where_add<-min(seq(1,ncol(vintageRecord),1)[vintageRecord[labelvintage==vintage,]==0])
      vintageRecord[labelvintage==vintage,where_add]<-max(vintageRecord)+1
      #I will coalesce two free, it doesn't matter which two
      if (sum(ava_free>1)>=1){
        ##Just pick the one with more than one and substract two and also just substract two from the family
        total_free[potential$PP[,1]==vintage & total_free>1][1]<- total_free[potential$PP[,1]==vintage & total_free>1][1]-2
        PP<-potential$PP
        PP[PP[,1]==vintage & PP[,5]==1,3][1]<-PP[PP[,1]==vintage & PP[,5]==1 & PP[,3]>0,3][1]-1
        if (PP[PP[,1]==vintage & PP[,5]==1,3][1]==0){
          PP[PP[,1]==vintage & PP[,5]==1,5][1]<-0
        }
      }else{
        #eliminate two free
        #just pick the first two buckets and substract 1
        where<-min(seq(1,length(total_free))[potential$PP[,1]==vintage & total_free>0])
        total_free[potential$PP[,1]==vintage & total_free>0]<-total_free[potential$PP[,1]==vintage & total_free>0]-c(1,1,rep(0,sum(ava_free)-2))
        #eliminate just one lineage. I need to take one of the leaves and make it 0 and also size equals 0
        PP<-potential$PP
        PP[where,3]<-PP[where,3]-1
        #I make it not leaf if indead this disapears
        if (PP[where,3]==0){
          PP[where,5]<-0
        }
        
      }
    }
    if (therandom[2]==1){
      whoV<-sample(seq(1,sum(vintageRecord[labelvintage==vintage,]>0)),1)
      whovalV<-vintageRecord[labelvintage==vintage,][whoV]
      
      where_add<-min(seq(1,ncol(vintageRecord),1)[vintageRecord[labelvintage==vintage,]==0])
      prob<-c(prob,1/sum(vintageRecord[labelvintage==vintage,]>0))
      
      vintageRecord[labelvintage==vintage,where_add]<-max(vintageRecord)+1
      vintageRecord[labelvintage==vintage,]<-c(vintageRecord[labelvintage==vintage,-whoV],0)
      # if (whovalV>1) {
      F[countIndex,1:(whovalV)]<-F[(countIndex-1),1:(whovalV)]-1
      F[countIndex,(whovalV+1):(countIndex-1)]<-F[(countIndex-1),(whovalV+1):(countIndex-1)]-2
      # }else{
      #   F[countIndex,1]<-F[countIndex-1,1]-1
      #   F[countIndex,(1+whovalV):(countIndex-1)]<-F[(countIndex-1),(1+whovalV):(countIndex-1)]-2
      # }
      countIndex<-countIndex+1
      #there is one free and one vintaged
      #Eliminate one free
      where<-min(seq(1,length(total_free))[potential$PP[,1]==vintage & total_free>0])
      total_free[where]<-total_free[where]-1
      ##From PP, I eliminate one "free" lineage
      PP<-potential$PP
      PP[where,3]<-PP[where,3]-1
      if (PP[where,3]==0){
        PP[where,5]<-0
      }
    }
    if (therandom[3]==1){
      whoV<-sample(seq(1,sum(vintageRecord[labelvintage==vintage,]>0)),2)
      whoV1<-min(whoV)
      whoV2<-max(whoV)
      whovalV1<-vintageRecord[labelvintage==vintage,][whoV1]
      whovalV2<-vintageRecord[labelvintage==vintage,][whoV2]
      denom<-sum(vintageRecord[labelvintage==vintage,]>0)
      prob<-c(prob,2/(denom*(denom-1)))
      where_add<-min(seq(1,ncol(vintageRecord),1)[vintageRecord[labelvintage==vintage,]==0])
      vintageRecord[labelvintage==vintage,where_add]<-max(vintageRecord)+1
      vintageRecord[labelvintage==vintage,]<-c(vintageRecord[labelvintage==vintage,-whoV],0,0)
      F[countIndex,1:(whovalV1)]<-F[(countIndex-1),1:(whovalV1)]
      F[countIndex,(whovalV1+1):(whovalV2)]<-F[(countIndex-1),(whovalV1+1):(whovalV2)]-1
      F[countIndex,(whovalV2+1):(countIndex-1)]<-F[(countIndex-1),(whovalV2+1):(countIndex-1)]-2
      
      countIndex<-countIndex+1
      
      #just eliminate one that is not free from PP
      #where<-min(seq(1,length(total_free))[potential$PP[,1]==vintage & potential$PP[,3]>0 & total_free==0])
      where<-min(seq(1,length(total_free))[potential$PP[,1]==vintage & potential$PP[,3]>0])
      PP<-potential$PP
      #print(paste("total free:", total_free,sep=""))
      #print(paste(where," o where ",PP,sep=""))
      
       PP[where,3]<-PP[where,3]-1
      if (PP[where,3]==0){
        PP[where,5]<-0
      }
    }
    
    if (sum(PP[,5])>1){
      check2_2<-check2(PP,sufficient,vintage,vintageRecord,labelvintage)
      PP<-check2_2$PP
      vintageRecord<-check2_2$vintageRecord
    }
    total_vintage<-c(total_vintage,vintage)
  }
  return(list(total_vintage=total_vintage[-1],change_F=change_F[-1],prob=prob[-1],F=F))
}

tree_from_F <- function(matF, coal_times){
  #generate an ape tree (Jaime's code)
  #F is the actual form used in the code that differs from the paper's notation
  
  n= dim(matF)[1]
  edge=matrix(rep(0,6*n-6),ncol=3)
  edge[,2]= (1:(2*n-2))
  vintages = c()
  times=c(rep(0,n),coal_times)
  for (j in n:2){
    new_node = 2*n-j+1
    next_leaf = intersect(which(edge[,1]==0),1:n)[1]
    F_difference = rev(matF[,j]-matF[,j-1])
    
    if (F_difference[1]==2){
      edge[next_leaf:(next_leaf+1),1]=new_node
      vintages=c(vintages, new_node)
    }
    else if (F_difference[1]==1){
      selected_vintage = which(F_difference == 2)[1]+n-1
      edge[selected_vintage,1]=new_node
      edge[next_leaf,1]=new_node
      vintages = c(vintages[vintages != selected_vintage],new_node)
    }
    else {
      selected_vintage1 =which(F_difference == 1)[1]+n-1
      selected_vintage2 =which(F_difference == 2)[1]+n-1
      edge[selected_vintage1,1]=new_node
      edge[selected_vintage2,1]=new_node
      vintages = vintages[vintages!=selected_vintage1]
      vintages = vintages[vintages!=selected_vintage2]
      vintages<-c(vintages,new_node)
    }
  }
  #edge[5:8,2]=c(6,7,8,5)
  edge[1:n,]=edge[order(edge[1:n,2]),]
  #edge=edge[order(edge[,1]),]
  
  for (j in 1:(2*n-2)) {
    #I don't understand this
    edge[j,3]=times[edge[j,1]]-times[edge[j,2]]
  }
  edge[,1]=3*n-edge[,1]
  edge[-(1:n),2]=3*n-edge[-(1:n),2]
  
  final_tree=rcoal(n,br=coal_times)
  final_tree$edge=edge[,-3]
  final_tree$edge.length=edge[,3]
  final_tree$Nnode=n-1
  class(final_tree) <- "phylo"
  final_tree <- reorder(final_tree,"postorder")
  final_tree$edge[final_tree$edge[, 2] <= n, 2] <- 1:n
  return(final_tree)
}

#in the same proposal keep track vintages per potential family and construct F directly

F_iter3<-function(sufficient){ 
  #Update April 27, 2017. This function samples F and allocation given the PP 
  n<-sum(sufficient[sufficient[,2]==0,3]) #descending from root node 0
  F<-matrix(0,nrow=n-1,ncol=n-1)
  diag(F)<-seq(n,2)
  labelvintage<-sort(unique(c(sufficient[,2],sufficient[,1])))
  vintageRecord<-matrix(0,nrow=length(labelvintage),ncol=n-1)
  PP<-sufficient
  total_free<-PP[,3]*PP[,5]
  ##In sufficient_back I will retain the original data to identify original free
  total_vintage<-0
  prob<-0
  change_F<-0
  countIndex<-2
  
  while(length(total_vintage)<(n-1)){
    potential<-F_iter(PP)
    
    if (length(potential$tot)>1){
      potential_size2<-choose(potential$tot,2)
      choose<-rmultinom(1,1,potential_size2)
      vintage<-potential$pot[choose==1]
      prob<-c(prob,potential_size2[choose==1]/sum(potential_size2))
    }else{
      vintage<-potential$pot
      prob<-c(prob,1)
    }
    ##lineages to choose from for the coalescence
    available<-potential$PP[potential$PP[,1]==vintage,3]*potential$PP[potential$PP[,1]==vintage,5]
    ##this are lineages available, some vintaged and some free
    ava_free<-total_free[potential$PP[,1]==vintage]
    tochange<-sum(available)
    
    if (tochange>1){
      free<-sum(ava_free)
      counts_Free<-c(choose(free,2),free*(tochange-free),choose(tochange-free,2))
      therandom<-rmultinom(1,1,counts_Free)
      prob<-c(prob,counts_Free[therandom==1]/sum(counts_Free))
      change_F<-c(change_F,c(2,1,0)[therandom==1])
    }else{
      #do something to define therandom
    }
    if (therandom[1]==1){
      F[countIndex,1:(countIndex-1)]<-F[countIndex-1,1:(countIndex-1)]-2
      countIndex<-countIndex+1
      where_add<-min(seq(1,ncol(vintageRecord),1)[vintageRecord[labelvintage==vintage,]==0])
      vintageRecord[labelvintage==vintage,where_add]<-max(vintageRecord)+1
      #I will coalesce two free, it doesn't matter which two
      if (sum(ava_free>1)>=1){
        ##Just pick the one with more than one and substract two and also just substract two from the family
        total_free[potential$PP[,1]==vintage & total_free>1][1]<- total_free[potential$PP[,1]==vintage & total_free>1][1]-2
        PP<-potential$PP
        PP[PP[,1]==vintage & PP[,5]==1,3][1]<-PP[PP[,1]==vintage & PP[,5]==1 & PP[,3]>0,3][1]-1
        if (PP[PP[,1]==vintage & PP[,5]==1,3][1]==0){
          PP[PP[,1]==vintage & PP[,5]==1,5][1]<-0
        }
      }else{
        #eliminate two free
        #just pick the first two buckets and substract 1
        where<-min(seq(1,length(total_free))[potential$PP[,1]==vintage & total_free>0])
        total_free[potential$PP[,1]==vintage & total_free>0]<-total_free[potential$PP[,1]==vintage & total_free>0]-c(1,1,rep(0,sum(ava_free)-2))
        #eliminate just one lineage. I need to take one of the leaves and make it 0 and also size equals 0
        PP<-potential$PP
        PP[where,3]<-PP[where,3]-1
        #I make it not leaf if indead this disapears
        if (PP[where,3]==0){
          PP[where,5]<-0
        }
        
      }
    }
    if (therandom[2]==1){
      whoV<-sample(seq(1,sum(vintageRecord[labelvintage==vintage,]>0)),1)
      whovalV<-vintageRecord[labelvintage==vintage,][whoV]
      
      where_add<-min(seq(1,ncol(vintageRecord),1)[vintageRecord[labelvintage==vintage,]==0])
      prob<-c(prob,1/sum(vintageRecord[labelvintage==vintage,]>0))
      
      vintageRecord[labelvintage==vintage,where_add]<-max(vintageRecord)+1
      vintageRecord[labelvintage==vintage,]<-c(vintageRecord[labelvintage==vintage,-whoV],0)
     # if (whovalV>1) {
        F[countIndex,1:(whovalV)]<-F[(countIndex-1),1:(whovalV)]-1
        F[countIndex,(whovalV+1):(countIndex-1)]<-F[(countIndex-1),(whovalV+1):(countIndex-1)]-2
      # }else{
      #   F[countIndex,1]<-F[countIndex-1,1]-1
      #   F[countIndex,(1+whovalV):(countIndex-1)]<-F[(countIndex-1),(1+whovalV):(countIndex-1)]-2
      # }
      countIndex<-countIndex+1
      #there is one free and one vintaged
      #Eliminate one free
      where<-min(seq(1,length(total_free))[potential$PP[,1]==vintage & total_free>0])
      total_free[where]<-total_free[where]-1
      ##From PP, I eliminate one "free" lineage
      PP<-potential$PP
      PP[where,3]<-PP[where,3]-1
      if (PP[where,3]==0){
        PP[where,5]<-0
      }
    }
    if (therandom[3]==1){
      whoV<-sample(seq(1,sum(vintageRecord[labelvintage==vintage,]>0)),2)
      whoV1<-min(whoV)
      whoV2<-max(whoV)
      whovalV1<-vintageRecord[labelvintage==vintage,][whoV1]
      whovalV2<-vintageRecord[labelvintage==vintage,][whoV2]
      denom<-sum(vintageRecord[labelvintage==vintage,]>0)
      prob<-c(prob,2/(denom*(denom-1)))
      where_add<-min(seq(1,ncol(vintageRecord),1)[vintageRecord[labelvintage==vintage,]==0])
      vintageRecord[labelvintage==vintage,where_add]<-max(vintageRecord)+1
      vintageRecord[labelvintage==vintage,]<-c(vintageRecord[labelvintage==vintage,-whoV],0,0)
      F[countIndex,1:(whovalV1)]<-F[(countIndex-1),1:(whovalV1)]
      F[countIndex,(whovalV1+1):(whovalV2)]<-F[(countIndex-1),(whovalV1+1):(whovalV2)]-1
      F[countIndex,(whovalV2+1):(countIndex-1)]<-F[(countIndex-1),(whovalV2+1):(countIndex-1)]-2
      
      countIndex<-countIndex+1
      
      #just eliminate one that is not free from PP
      where<-min(seq(1,length(total_free))[potential$PP[,1]==vintage & potential$PP[,3]>0 & total_free==0])
      PP<-potential$PP
      PP[where,3]<-PP[where,3]-1
      if (PP[where,3]==0){
        PP[where,5]<-0
      }
    }
    
    if (sum(PP[,5])>1){
      check2_2<-check2(PP,sufficient,vintage,vintageRecord,labelvintage)
      PP<-check2_2$PP
      vintageRecord<-check2_2$vintageRecord
      }
     total_vintage<-c(total_vintage,vintage)
  }
  return(list(total_vintage=total_vintage[-1],change_F=change_F[-1],prob=prob[-1],F=F))
}

F_iter5<-function(sufficient_tot){ 
  #Update April 27, 2017. This function samples F and allocation given the PP like F_iter3 but uniformly
  #Allocation is sampled uniformly
  n<-sufficient_tot$n
  sufficient<-sufficient_tot$nodes
 # n<-sum(sufficient[sufficient[,2]==0,3]) #descending from root node 0
  F<-matrix(0,nrow=n-1,ncol=n-1)
  diag(F)<-seq(n,2)
  tot<-sum(sufficient[sufficient[,2]==0,3])
  if (tot<n){
    label<-max(sufficient[,1])+1
    for (i in 1:(n-tot)){
    sufficient<-rbind(sufficient,c(label,0,1,0,1))
    }}
  labelvintage<-sort(unique(c(sufficient[,2],sufficient[,1])))
  vintageRecord<-matrix(0,nrow=length(labelvintage),ncol=n-1)
  PP<-sufficient
  total_free<-PP[,3]*PP[,5]
  ##In sufficient_back I will retain the original data to identify original free
  total_vintage<-0
  prob<-0
  change_F<-0
  countIndex<-2
  
  while(length(total_vintage)<(n-1)){
    potential<-F_iter(PP)
    
    if (length(potential$tot)>1){
      potential_size2<-choose(potential$tot,2)
      choose<-rmultinom(1,1,rep(1,length(potential_size2)))
      #choose<-rmultinom(1,1,potential_size2)
      vintage<-potential$pot[choose==1]
      prob<-c(prob,1/length(potential_size2))
     # prob<-c(prob,sum(potential_size2*choose)/sum(potential_size2))
    }else{
      vintage<-potential$pot
      prob<-c(prob,1)
    }
    ##lineages to choose from for the coalescence
    available<-potential$PP[potential$PP[,1]==vintage,3]*potential$PP[potential$PP[,1]==vintage,5]
    ##this are lineages available, some vintaged and some free
    ava_free<-total_free[potential$PP[,1]==vintage]
    tochange<-sum(available)
    
      free<-sum(ava_free)
      counts_Free<-c(choose(free,2),free*(tochange-free),choose(tochange-free,2))
      therandom<-rmultinom(1,1,counts_Free)
      prob<-c(prob,counts_Free[therandom==1]/sum(counts_Free))
      change_F<-c(change_F,c(2,1,0)[therandom==1])
   # print(therandom)
    if (therandom[1]==1){
      F[countIndex,1:(countIndex-1)]<-F[countIndex-1,1:(countIndex-1)]-2
      countIndex<-countIndex+1
      where_add<-min(seq(1,ncol(vintageRecord),1)[vintageRecord[labelvintage==vintage,]==0])
      vintageRecord[labelvintage==vintage,where_add]<-max(vintageRecord)+1
      #I will coalesce two free, it doesn't matter which two
      if (sum(ava_free>1)>=1){
        ##Just pick the one with more than one and substract two and also just substract two from the family
        total_free[potential$PP[,1]==vintage & total_free>1][1]<- total_free[potential$PP[,1]==vintage & total_free>1][1]-2
        PP<-potential$PP
        PP[PP[,1]==vintage & PP[,5]==1,3][1]<-PP[PP[,1]==vintage & PP[,5]==1 & PP[,3]>0,3][1]-1
        if (PP[PP[,1]==vintage & PP[,5]==1,3][1]==0){
          PP[PP[,1]==vintage & PP[,5]==1,5][1]<-0
        }
      }else{
        #eliminate two free
        #just pick the first two buckets and substract 1
        where<-min(seq(1,length(total_free))[potential$PP[,1]==vintage & total_free>0])
        total_free[potential$PP[,1]==vintage & total_free>0]<-total_free[potential$PP[,1]==vintage & total_free>0]-c(1,1,rep(0,sum(ava_free)-2))
        #eliminate just one lineage. I need to take one of the leaves and make it 0 and also size equals 0
        PP<-potential$PP
        PP[where,3]<-PP[where,3]-1
        #I make it not leaf if indead this disapears
        if (PP[where,3]==0){
          PP[where,5]<-0
        }
        
      }
    }
    if (therandom[2]==1){
      whoV<-sample(seq(1,sum(vintageRecord[labelvintage==vintage,]>0)),1)
      whovalV<-vintageRecord[labelvintage==vintage,][whoV]
      
      where_add<-min(seq(1,ncol(vintageRecord),1)[vintageRecord[labelvintage==vintage,]==0])
      prob<-c(prob,1/sum(vintageRecord[labelvintage==vintage,]>0))
      
      vintageRecord[labelvintage==vintage,where_add]<-max(vintageRecord)+1
      vintageRecord[labelvintage==vintage,]<-c(vintageRecord[labelvintage==vintage,-whoV],0)
      # if (whovalV>1) {
      F[countIndex,1:(whovalV)]<-F[(countIndex-1),1:(whovalV)]-1
      F[countIndex,(whovalV+1):(countIndex-1)]<-F[(countIndex-1),(whovalV+1):(countIndex-1)]-2
      # }else{
      #   F[countIndex,1]<-F[countIndex-1,1]-1
      #   F[countIndex,(1+whovalV):(countIndex-1)]<-F[(countIndex-1),(1+whovalV):(countIndex-1)]-2
      # }
      countIndex<-countIndex+1
      #there is one free and one vintaged
      #Eliminate one free
      where<-min(seq(1,length(total_free))[potential$PP[,1]==vintage & total_free>0])
      total_free[where]<-total_free[where]-1
      ##From PP, I eliminate one "free" lineage
      PP<-potential$PP
      PP[where,3]<-PP[where,3]-1
      if (PP[where,3]==0){
        PP[where,5]<-0
      }
    }
    if (therandom[3]==1){
      whoV<-sample(seq(1,sum(vintageRecord[labelvintage==vintage,]>0)),2)
      whoV1<-min(whoV)
      whoV2<-max(whoV)
      whovalV1<-vintageRecord[labelvintage==vintage,][whoV1]
      whovalV2<-vintageRecord[labelvintage==vintage,][whoV2]
      denom<-sum(vintageRecord[labelvintage==vintage,]>0)
      prob<-c(prob,2/(denom*(denom-1)))
      where_add<-min(seq(1,ncol(vintageRecord),1)[vintageRecord[labelvintage==vintage,]==0])
      vintageRecord[labelvintage==vintage,where_add]<-max(vintageRecord)+1
      vintageRecord[labelvintage==vintage,]<-c(vintageRecord[labelvintage==vintage,-whoV],0,0)
      F[countIndex,1:(whovalV1)]<-F[(countIndex-1),1:(whovalV1)]
      F[countIndex,(whovalV1+1):(whovalV2)]<-F[(countIndex-1),(whovalV1+1):(whovalV2)]-1
      F[countIndex,(whovalV2+1):(countIndex-1)]<-F[(countIndex-1),(whovalV2+1):(countIndex-1)]-2
      
      countIndex<-countIndex+1
      
      #just eliminate one that is not free from PP
      where<-min(seq(1,length(total_free))[potential$PP[,1]==vintage & potential$PP[,3]>0 & total_free==0])
      PP<-potential$PP
      PP[where,3]<-PP[where,3]-1
      if (PP[where,3]==0){
        PP[where,5]<-0
      }
    }
    
    if (sum(PP[,5])>1){
      check2_2<-check2(PP,sufficient,vintage,vintageRecord,labelvintage)
      PP<-check2_2$PP
      vintageRecord<-check2_2$vintageRecord
    }
    total_vintage<-c(total_vintage,vintage)
  }
  return(list(total_vintage=total_vintage[-1],change_F=change_F[-1],prob=prob[-1],F=F))
}
# 
# cal_q<-function(Finfo,sufficient,oldsuff){
#   #Update May 10, 2017.
#   #This function calculates the probability Q(F|Data)
#   # 
#   # allocs<-qF(oldsuff,Finfo$F)
#   # for (j in 1:ncol(allocs$alloc)){
#   #  index<-sort(allocs$alloc[,j],index.return=T)
#   #  if (j==1){vintage_list<-allocs$parent_list[index$ix]}else{
#   #    vintage_list<-cbind(vintage_list,allocs$parent_list[index$ix])}
#   #  }
# 
# #n<-sum(sufficient[sufficient[,2]==0,3]) #descending from root node 0
# #F<-matrix(0,nrow=n-1,ncol=n-1)
# #diag(F)<-seq(n,2)
# labelvintage<-sort(unique(c(sufficient[,2],sufficient[,1])))
# vintageRecord<-matrix(0,nrow=length(labelvintage),ncol=n-1)
# PP<-sufficient
# total_free<-PP[,3]*PP[,5]
# ##In sufficient_back I will retain the original data to identify original free
# total_vintage<-0
# prob<-0
# #change_F<-0
# countIndex<-2
# 
# while(length(total_vintage)<(n-1)){
#   potential<-F_calc(PP,total_free)
#   
#   # ind<-rep(1,length(potential$pot))
#   # if (countIndex==2){
#   # for (j in 1:length(potential$pot)){
#   #   ind[j]<-ifelse(sum(vintage_list[countIndex-1,]==potential$pot[j])==0,0,1)
#   # }}else{
#   #   if (sum(vintage_list[countIndex-1,]==0)==0){
#   #   for (j in 1:length(potential$pot)){
#   #     ref<-total_vintage[countIndex-1]
#   #     ind[j]<-ifelse(sum(vintage_list[countIndex-1,vintage_list[countIndex-2,]==ref]==potential$pot[j])==0,0,1)
#   #   } 
#   # }
#   # }
#   # potential$tot<-potential$tot[ind==1]
#   # potential$pot<-potential$pot[ind==1]
#   if (length(potential$tot)>1){
#     #potential_size2<-choose(potential$tot,2)
#     #choose<-rmultinom(1,1,rep(1,length(potential_size2)))
#     vintage<-Finfo$total_vintage[countIndex-1]
#    # prob<-c(prob,1/length(potential_size2))
#   #  prob<-c(prob,1/length(potential$tot))
#   }else{
#     vintage<-potential$pot
#    # prob<-c(prob,1)
#   }
#   available<-potential$PP[potential$PP[,1]==vintage,3]*potential$PP[potential$PP[,1]==vintage,5]
#   ava_free<-total_free[potential$PP[,1]==vintage]
#   free<-potential$fre
#   tochange<-potential$tot
#    # free<-sum(ava_free)
#     counts_Free<-c(sum(choose(free,2)),sum(free*(tochange-free)),sum(choose(tochange-free,2)))
#     therandom<-c(0,0,1)
#     if (Finfo$change_F[countIndex-1]==2){ therandom<-c(1,0,0)}
#     if (Finfo$change_F[countIndex-1]==1){ therandom<-c(0,1,0)}
#     
#     prob<-c(prob,counts_Free[therandom==1]/sum(counts_Free))
#     if (therandom[1]==1){
#     countIndex<-countIndex+1
#     where_add<-min(seq(1,ncol(vintageRecord),1)[vintageRecord[labelvintage==vintage,]==0])
#     vintageRecord[labelvintage==vintage,where_add]<-max(vintageRecord)+1
#     if (sum(ava_free>1)>=1){
#       total_free[potential$PP[,1]==vintage & total_free>1][1]<- total_free[potential$PP[,1]==vintage & total_free>1][1]-2
#       PP<-potential$PP
#       PP[PP[,1]==vintage & PP[,5]==1,3][1]<-PP[PP[,1]==vintage & PP[,5]==1 & PP[,3]>0,3][1]-1
#       if (PP[PP[,1]==vintage & PP[,5]==1,3][1]==0){
#         PP[PP[,1]==vintage & PP[,5]==1,5][1]<-0
#       }
#     }else{
#       where<-min(seq(1,length(total_free))[potential$PP[,1]==vintage & total_free>0])
#       total_free[potential$PP[,1]==vintage & total_free>0]<-total_free[potential$PP[,1]==vintage & total_free>0]-c(1,1,rep(0,sum(ava_free)-2))
#       PP<-potential$PP
#       PP[where,3]<-PP[where,3]-1
#       if (PP[where,3]==0){
#         PP[where,5]<-0
#       }
#       
#     }
#   }
#   if (therandom[2]==1){
#     whoV<-sample(seq(1,sum(vintageRecord[labelvintage==vintage,]>0)),1)
#     whovalV<-vintageRecord[labelvintage==vintage,][whoV]
#     
#     where_add<-min(seq(1,ncol(vintageRecord),1)[vintageRecord[labelvintage==vintage,]==0])
#     prob<-c(prob,1/sum(vintageRecord[labelvintage==vintage,]>0))
#     
#     vintageRecord[labelvintage==vintage,where_add]<-max(vintageRecord)+1
#     vintageRecord[labelvintage==vintage,]<-c(vintageRecord[labelvintage==vintage,-whoV],0)
#     countIndex<-countIndex+1
#     where<-min(seq(1,length(total_free))[potential$PP[,1]==vintage & total_free>0])
#     total_free[where]<-total_free[where]-1
#     PP<-potential$PP
#     PP[where,3]<-PP[where,3]-1
#     if (PP[where,3]==0){
#       PP[where,5]<-0
#     }
#   }
#   if (therandom[3]==1){
#     whoV<-sample(seq(1,sum(vintageRecord[labelvintage==vintage,]>0)),2)
#     whoV1<-min(whoV)
#     whoV2<-max(whoV)
#     whovalV1<-vintageRecord[labelvintage==vintage,][whoV1]
#     whovalV2<-vintageRecord[labelvintage==vintage,][whoV2]
#     denom<-sum(vintageRecord[labelvintage==vintage,]>0)
#     prob<-c(prob,2/(denom*(denom-1)))
#     where_add<-min(seq(1,ncol(vintageRecord),1)[vintageRecord[labelvintage==vintage,]==0])
#     vintageRecord[labelvintage==vintage,where_add]<-max(vintageRecord)+1
#     vintageRecord[labelvintage==vintage,]<-c(vintageRecord[labelvintage==vintage,-whoV],0,0)
#     countIndex<-countIndex+1
#     where<-min(seq(1,length(total_free))[potential$PP[,1]==vintage & potential$PP[,3]>0 & total_free==0])
#     PP<-potential$PP
#     PP[where,3]<-PP[where,3]-1
#     if (PP[where,3]==0){
#       PP[where,5]<-0
#     }
#   }
#   
#   if (sum(PP[,5])>1){
#     check2_2<-check2(PP,sufficient,vintage,vintageRecord,labelvintage)
#     PP<-check2_2$PP
#     vintageRecord<-check2_2$vintageRecord
#   }
#   total_vintage<-c(total_vintage,vintage)
# }
# return(list(total_vintage=total_vintage[-1],prob=prob[-1]))
# }


cal_q<-function(Finfo,sufficient,oldsuff){
  #Update May 21, 2017.
  #This function calculates the probability Q(F|Data)=Sum_{A}Q(F,A|Data)
   ##I need to check this is right
   bb_use<-tajima_info2(Finfo$F)
   allocs<-qF(oldsuff,Finfo$F,bb_use)
   vintage_list_temp<-allocs$parent_list
   for (i in 1:length(allocs$node_list)){
     where<-seq(1,nrow(sufficient))[sufficient[,1]==allocs$node_list[i]][1]
     if (sufficient[where,3]>1 && sufficient[where,5]==1){
       vintage_list_temp[i]<-allocs$node_list[i]
     }
   }
   if (ncol(allocs$alloc)>1){
     #here I am having a mistake when I have a missing node
     for (j in 1:ncol(allocs$alloc)){

       tmp<-seq(1,nrow(allocs$alloc)-1)[diff(allocs$alloc[,j])==0]+1
       if (length(tmp)==0){
         index<-sort(allocs$alloc[,j],index.return=T)
         parent_list<-vintage_list_temp
       }else{
       index<-sort(allocs$alloc[-tmp,j],index.return=T)
       parent_list<-vintage_list_temp[-tmp]
       }
       #index<-sort(allocs$alloc[,j],index.return=T)
       if (j==1){vintage_list<-as.matrix(parent_list[index$ix])}else{
         vintage_list<-cbind(vintage_list,parent_list[index$ix])}
     }
     if (nrow(vintage_list)<nrow(Finfo$F)){
       all<-seq(1,nrow(Finfo$F))
       there<-index$x
       for (j in 1:nrow(Finfo$F)){
         if (there[j]!=j){
           val<-sufficient[sufficient[,1]==vintage_list[max(allocs$ancestries[,j]*seq(1,nrow(allocs$ancestries))),1],2]
           if (val==0){
             val<-sufficient[sufficient[,1]==vintage_list[max(allocs$ancestries[,j]*seq(1,nrow(allocs$ancestries))),1],1]
           }
           # extended_vintage_list[j,]<-rep(0,ncol(vintage_list))
           if (j>1){there<-c(there[1:(j-1)],j,there[j:length(there)])
           vintage_list<-rbind(vintage_list[1:(j-1),],rep(val,ncol(vintage_list)),vintage_list[j:nrow(vintage_list),])}else{
             there<-c(val,there)
             vintage_list<-rbind(rep(0,ncol(vintage_list)),vintage_list)}
         }
       }
     }
     
    
  tot_prob<-rep(0,ncol(vintage_list))
  for (s in 1:ncol(vintage_list)){
    
   total_vintage2<-vintage_list[,s]
  # for (l in unique(total_vintage2)){
  #   if (unique(sufficient[sufficient[,1]==l,5])==1){
  #   total_vintage2[total_vintage2==l]<-unique(sufficient[sufficient[,1]==l,2])
  # }}
  labelvintage<-sort(unique(c(sufficient[,2],sufficient[,1])))
  vintageRecord<-matrix(0,nrow=length(labelvintage),ncol=n-1)
  PP<-F_iter(sufficient)$PP
  total_free<-PP[,3]*PP[,5]
  ##In sufficient_back I will retain the original data to identify original free
  total_vintage<-0
  prob<-0
  #change_F<-0
  countIndex<-2
  
  while(length(total_vintage)<(n-1)){
    potential<-F_iter(PP)
    
    if (length(potential$tot)>1){
      potential_size2<-choose(potential$tot,2)
      #choose<-rmultinom(1,1,rep(1,length(potential_size2)))
      vintage<-total_vintage2[countIndex-1]
      #prob<-c(prob,1/length(potential_size2))
      
      #choose<-rmultinom(1,1,potential_size2)
      #vintage<-potential$pot[choose==1]
      choose<-rep(0,length(potential$tot))
      choose[potential$pot==vintage]<-1
      prob<-c(prob,1/length(potential_size2))
      #prob<-c(prob,sum(potential_size2*choose)/sum(potential_size2))
      
      
      #  prob<-c(prob,1/length(potential$tot))
    }else{
      vintage<-potential$pot
      # prob<-c(prob,1)
    }
    available<-potential$PP[potential$PP[,1]==vintage,3]*potential$PP[potential$PP[,1]==vintage,5]
    ava_free<-total_free[potential$PP[,1]==vintage]
    #free<-potential$fre
   # tochange<-potential$tot
    tochange<-sum(available)
    free<-sum(ava_free)
    counts_Free<-c(sum(choose(free,2)),sum(free*(tochange-free)),sum(choose(tochange-free,2)))
    therandom<-c(0,0,1)
    if (Finfo$change_F[countIndex-1]==2){ therandom<-c(1,0,0)}
    if (Finfo$change_F[countIndex-1]==1){ therandom<-c(0,1,0)}
    
    prob<-c(prob,counts_Free[therandom==1]/sum(counts_Free))
    if (therandom[1]==1){
      countIndex<-countIndex+1
      where_add<-min(seq(1,ncol(vintageRecord),1)[vintageRecord[labelvintage==vintage,]==0])
      vintageRecord[labelvintage==vintage,where_add]<-max(vintageRecord)+1
      PP<-potential$PP
      if (sum(ava_free>1)>=1){
        total_free[potential$PP[,1]==vintage & total_free>1][1]<- total_free[potential$PP[,1]==vintage & total_free>1][1]-2
        PP<-potential$PP
        PP[PP[,1]==vintage & PP[,5]==1,3][1]<-PP[PP[,1]==vintage & PP[,5]==1 & PP[,3]>0,3][1]-1
        if (PP[PP[,1]==vintage & PP[,5]==1,3][1]==0){
          PP[PP[,1]==vintage & PP[,5]==1,5][1]<-0
        }
      }else{
        #I will coalesce two singletons
        where<-min(seq(1,length(total_free))[potential$PP[,1]==vintage & total_free>0])
        total_free[potential$PP[,1]==vintage & total_free>0]<-total_free[potential$PP[,1]==vintage & total_free>0]-c(1,1,rep(0,sum(ava_free)-2))
        PP<-potential$PP
        PP[where,3]<-PP[where,3]-1
        if (PP[where,3]==0){
          PP[where,5]<-0
        }
        
      }
    }
    
    if (therandom[2]==1){
      whoV<-sample(seq(1,sum(vintageRecord[labelvintage==vintage,]>0)),1)
      whovalV<-vintageRecord[labelvintage==vintage,][whoV]
      
      where_add<-min(seq(1,ncol(vintageRecord),1)[vintageRecord[labelvintage==vintage,]==0])
      prob<-c(prob,1/sum(vintageRecord[labelvintage==vintage,]>0))
      
      vintageRecord[labelvintage==vintage,where_add]<-max(vintageRecord)+1
      vintageRecord[labelvintage==vintage,]<-c(vintageRecord[labelvintage==vintage,-whoV],0)
      countIndex<-countIndex+1
      where<-min(seq(1,length(total_free))[potential$PP[,1]==vintage & total_free>0])
      total_free[where]<-total_free[where]-1
      PP<-potential$PP
      PP[where,3]<-PP[where,3]-1
      if (PP[where,3]==0){
        PP[where,5]<-0
      }
    }
    if (therandom[3]==1){
      whoV<-sample(seq(1,sum(vintageRecord[labelvintage==vintage,]>0)),2)
      whoV1<-min(whoV)
      whoV2<-max(whoV)
      whovalV1<-vintageRecord[labelvintage==vintage,][whoV1]
      whovalV2<-vintageRecord[labelvintage==vintage,][whoV2]
      denom<-sum(vintageRecord[labelvintage==vintage,]>0)
      prob<-c(prob,2/(denom*(denom-1)))
      where_add<-min(seq(1,ncol(vintageRecord),1)[vintageRecord[labelvintage==vintage,]==0])
      vintageRecord[labelvintage==vintage,where_add]<-max(vintageRecord)+1
      vintageRecord[labelvintage==vintage,]<-c(vintageRecord[labelvintage==vintage,-whoV],0,0)
      countIndex<-countIndex+1
      where<-min(seq(1,length(total_free))[potential$PP[,1]==vintage & potential$PP[,3]>0 & total_free==0])
      PP<-potential$PP
      PP[where,3]<-PP[where,3]-1
      if (PP[where,3]==0){
        PP[where,5]<-0
      }
    }
    
    if (sum(PP[,5])>1){
      check2_2<-check2(PP,sufficient,vintage,vintageRecord,labelvintage)
      PP<-check2_2$PP
      vintageRecord<-check2_2$vintageRecord
    }
    total_vintage<-c(total_vintage,vintage)
  }
  tot_prob[s]<-prod(prob[-1])
  }
  return(list(vintage_list=vintage_list,tot_prob=tot_prob))}
else{return(list(vintage_list=Finfo$total_vintage,tot_prob=prod(Finfo$prob)))}
}

F_iter4<-function(sufficient,change_F){ 
  #Update May 8, 2017. This function samples F and allocation given the PP and the first column of F
  #The harder constrains are the cherries so I will coalesce the cherries first
  n<-sum(sufficient[sufficient[,2]==0,3]) #descending from root node 0
  F<-matrix(0,nrow=n-1,ncol=n-1)
  diag(F)<-seq(n,2)
  F[2:(n-1),1]<-n-cumsum(change_F)
  labelvintage<-sort(unique(c(sufficient[,2],sufficient[,1])))
  vintageRecord<-matrix(0,nrow=length(labelvintage),ncol=n-1)
  PP<-sufficient
  total_free<-PP[,3]*PP[,5]
  ##In sufficient_back I will retain the original data to identify original free
  total_vintage<-0
  prob<-0
  #change_F<-0
  countIndex<-2
  ##Work with the cherries in the data
  cherries<-sufficient[sufficient[,3]==2,1]
  changeIndex<-seq(1,length(change_F))[change_F==2]
  change<-change_F[1]
  while(length(total_vintage)<(n-1)){
    potential<-F_iter(PP)
    potential$pot<-potential$pot[potential$tot>=change]
    potential$tot<-potential$tot[potential$tot>=change]
    if (change==2 & length(cherries)>0){
      for (j in potential$pot){
        if (sum(cherries==j)==0){
          where<-seq(1,length(potential$pot))[potential$pot==j]
          potential$pot<-potential$pot[-where]
          potential$tot<-potential$tot[-where]
        }
      }
    }
    
    if (length(potential$tot)>1){
      potential_size2<-choose(potential$tot,2)
      choose<-rmultinom(1,1,potential_size2)
      vintage<-potential$pot[choose==1]
      prob<-c(prob,potential_size2[choose==1]/sum(potential_size2))
    }else{
      vintage<-potential$pot
      prob<-c(prob,1)
    }
    
    ##lineages to choose from for the coalescence
    available<-potential$PP[potential$PP[,1]==vintage,3]*potential$PP[potential$PP[,1]==vintage,5]
    ##this are lineages available, some vintaged and some free
    ava_free<-total_free[potential$PP[,1]==vintage]
    tochange<-sum(available)
    
   
      #free<-sum(ava_free)
     # counts_Free<-c(choose(free,2),free*(tochange-free),choose(tochange-free,2))
      #therandom<-rmultinom(1,1,counts_Free)
      #prob<-c(prob,counts_Free[therandom==1]/sum(counts_Free))
      #change_F<-c(change_F,c(2,1,0)[therandom==1])
   
    if (change==2){
      if (length(cherries)>0){cherries<-cherries[cherries!=vintage]}
      if (countIndex>2){F[countIndex,2:(countIndex-1)]<- F[(countIndex-1),2:(countIndex-1)]-2}
      countIndex<-countIndex+1
      where_add<-min(seq(1,ncol(vintageRecord),1)[vintageRecord[labelvintage==vintage,]==0])
      vintageRecord[labelvintage==vintage,where_add]<-max(vintageRecord)+1
      #I will coalesce two free, it doesn't matter which two
      if (sum(ava_free>1)>=1){
        ##Just pick the one with more than one and substract two and also just substract two from the family
        total_free[potential$PP[,1]==vintage & total_free>1][1]<- total_free[potential$PP[,1]==vintage & total_free>1][1]-2
        PP<-potential$PP
        PP[PP[,1]==vintage & PP[,5]==1,3][1]<-PP[PP[,1]==vintage & PP[,5]==1 & PP[,3]>0,3][1]-1
        if (PP[PP[,1]==vintage & PP[,5]==1,3][1]==0){
          PP[PP[,1]==vintage & PP[,5]==1,5][1]<-0
        }
      }else{
        #eliminate two free
        #just pick the first two buckets and substract 1
        where<-min(seq(1,length(total_free))[potential$PP[,1]==vintage & total_free>0])
        total_free[potential$PP[,1]==vintage & total_free>0]<-total_free[potential$PP[,1]==vintage & total_free>0]-c(1,1,rep(0,sum(ava_free)-2))
        #eliminate just one lineage. I need to take one of the leaves and make it 0 and also size equals 0
        PP<-potential$PP
        PP[where,3]<-PP[where,3]-1
        #I make it not leaf if indead this disapears
        if (PP[where,3]==0){
          PP[where,5]<-0
        }
        
      }
    }
    if (change==1){
      
      whoV<-sample(seq(1,sum(vintageRecord[labelvintage==vintage,]>0)),1)
      whovalV<-vintageRecord[labelvintage==vintage,][whoV]
      
      where_add<-min(seq(1,ncol(vintageRecord),1)[vintageRecord[labelvintage==vintage,]==0])
      prob<-c(prob,1/sum(vintageRecord[labelvintage==vintage,]>0))
      
      vintageRecord[labelvintage==vintage,where_add]<-max(vintageRecord)+1
      vintageRecord[labelvintage==vintage,]<-c(vintageRecord[labelvintage==vintage,-whoV],0)
      # if (whovalV>1) {
      F[countIndex,1:(whovalV)]<-F[(countIndex-1),1:(whovalV)]-1
      F[countIndex,(whovalV+1):(countIndex-1)]<-F[(countIndex-1),(whovalV+1):(countIndex-1)]-2
      # }else{
      #   F[countIndex,1]<-F[countIndex-1,1]-1
      #   F[countIndex,(1+whovalV):(countIndex-1)]<-F[(countIndex-1),(1+whovalV):(countIndex-1)]-2
      # }
      countIndex<-countIndex+1
      #there is one free and one vintaged
      #Eliminate one free
      where<-min(seq(1,length(total_free))[potential$PP[,1]==vintage & total_free>0])
      total_free[where]<-total_free[where]-1
      ##From PP, I eliminate one "free" lineage
      PP<-potential$PP
      PP[where,3]<-PP[where,3]-1
      if (PP[where,3]==0){
        PP[where,5]<-0
      }
    }
    if (change==0){
      whoV<-sample(seq(1,sum(vintageRecord[labelvintage==vintage,]>0)),2)
      whoV1<-min(whoV)
      whoV2<-max(whoV)
      whovalV1<-vintageRecord[labelvintage==vintage,][whoV1]
      whovalV2<-vintageRecord[labelvintage==vintage,][whoV2]
      denom<-sum(vintageRecord[labelvintage==vintage,]>0)
      prob<-c(prob,2/(denom*(denom-1)))
      where_add<-min(seq(1,ncol(vintageRecord),1)[vintageRecord[labelvintage==vintage,]==0])
      vintageRecord[labelvintage==vintage,where_add]<-max(vintageRecord)+1
      vintageRecord[labelvintage==vintage,]<-c(vintageRecord[labelvintage==vintage,-whoV],0,0)
      F[countIndex,1:(whovalV1)]<-F[(countIndex-1),1:(whovalV1)]
      F[countIndex,(whovalV1+1):(whovalV2)]<-F[(countIndex-1),(whovalV1+1):(whovalV2)]-1
      F[countIndex,(whovalV2+1):(countIndex-1)]<-F[(countIndex-1),(whovalV2+1):(countIndex-1)]-2
      
      countIndex<-countIndex+1
      
      #just eliminate one that is not free from PP
      where<-min(seq(1,length(total_free))[potential$PP[,1]==vintage & potential$PP[,3]>0 & total_free==0])
      PP<-potential$PP
      PP[where,3]<-PP[where,3]-1
      if (PP[where,3]==0){
        PP[where,5]<-0
      }
    }
    
    if (sum(PP[,5])>1){
      check2_2<-check2(PP,sufficient,vintage,vintageRecord,labelvintage)
      PP<-check2_2$PP
      vintageRecord<-check2_2$vintageRecord
    }
    total_vintage<-c(total_vintage,vintage)
    change<-change_F[(countIndex-1)]
  }
  return(list(total_vintage=total_vintage[-1],change_F=change_F,prob=prob[-1],F=F))
}

F_iter6<-function(sufficient,change_F){ 
  #allocation is uniform
  #Update May 8, 2017. This function samples F and allocation given the PP and the first column of F
  #The harder constrains are the cherries so I will coalesce the cherries first
  n<-sum(sufficient[sufficient[,2]==0,3]) #descending from root node 0
  F<-matrix(0,nrow=n-1,ncol=n-1)
  diag(F)<-seq(n,2)
  F[2:(n-1),1]<-n-cumsum(change_F)
  labelvintage<-sort(unique(c(sufficient[,2],sufficient[,1])))
  vintageRecord<-matrix(0,nrow=length(labelvintage),ncol=n-1)
  PP<-sufficient
  total_free<-PP[,3]*PP[,5]
  ##In sufficient_back I will retain the original data to identify original free
  total_vintage<-0
  prob<-0
  #change_F<-0
  countIndex<-2
  ##Work with the cherries in the data
  cherries<-sufficient[sufficient[,3]==2,1]
  cherries2<-sum(change_F==2)
  cor<-0
  if (length(cherries)==cherries2) {cor=1}
  changeIndex<-seq(1,length(change_F))[change_F==2]
  change<-change_F[1]
  while(length(total_vintage)<(n-1)){
    potential<-F_iter(PP)
    potential$pot<-potential$pot[potential$tot>=change]
    potential$tot<-potential$tot[potential$tot>=change]
    #if (length(cherries)==length(cherries2)) {cor=1}
     if (change==2 & length(cherries)>0 & cor==1){
      for (j in potential$pot){
        if (sum(cherries==j)==0){
          where<-seq(1,length(potential$pot))[potential$pot==j]
          potential$pot<-potential$pot[-where]
          potential$tot<-potential$tot[-where]
        }
      }
    }
    
    if (length(potential$tot)>1){
      potential_size2<-choose(potential$tot,2)
      choose<-rmultinom(1,1,rep(1,length(potential_size2)))
 #     choose<-rmultinom(1,1,potential_size2)
      vintage<-potential$pot[choose==1]
      prob<-c(prob,1/length(potential_size2))
    }else{
      vintage<-potential$pot
      prob<-c(prob,1)
    }
    
    ##lineages to choose from for the coalescence
    available<-potential$PP[potential$PP[,1]==vintage,3]*potential$PP[potential$PP[,1]==vintage,5]
    ##this are lineages available, some vintaged and some free
    ava_free<-total_free[potential$PP[,1]==vintage]
    tochange<-sum(available)
    
    
    #free<-sum(ava_free)
    # counts_Free<-c(choose(free,2),free*(tochange-free),choose(tochange-free,2))
    #therandom<-rmultinom(1,1,counts_Free)
    #prob<-c(prob,counts_Free[therandom==1]/sum(counts_Free))
    #change_F<-c(change_F,c(2,1,0)[therandom==1])
    
    if (change==2){
      if (sum(ava_free)<2){return(0)}
      if (length(cherries)>0){cherries<-cherries[cherries!=vintage]}
      if (countIndex>2){F[countIndex,2:(countIndex-1)]<- F[(countIndex-1),2:(countIndex-1)]-2}
      countIndex<-countIndex+1
      where_add<-min(seq(1,ncol(vintageRecord),1)[vintageRecord[labelvintage==vintage,]==0])
      vintageRecord[labelvintage==vintage,where_add]<-max(vintageRecord)+1
      #I will coalesce two free, it doesn't matter which two
      if (sum(ava_free>1)>=1){
        ##Just pick the one with more than one and substract two and also just substract two from the family
        total_free[potential$PP[,1]==vintage & total_free>1][1]<- total_free[potential$PP[,1]==vintage & total_free>1][1]-2
        PP<-potential$PP
        PP[PP[,1]==vintage & PP[,5]==1,3][1]<-PP[PP[,1]==vintage & PP[,5]==1 & PP[,3]>0,3][1]-1
        if (PP[PP[,1]==vintage & PP[,5]==1,3][1]==0){
          PP[PP[,1]==vintage & PP[,5]==1,5][1]<-0
        }
      }else{
        #eliminate two free
        #just pick the first two buckets and substract 1
        where<-min(seq(1,length(total_free))[potential$PP[,1]==vintage & total_free>0])
        total_free[potential$PP[,1]==vintage & total_free>0]<-total_free[potential$PP[,1]==vintage & total_free>0]-c(1,1,rep(0,sum(ava_free)-2))
        #eliminate just one lineage. I need to take one of the leaves and make it 0 and also size equals 0
        PP<-potential$PP
        PP[where,3]<-PP[where,3]-1
        #I make it not leaf if indead this disapears
        if (PP[where,3]==0){
          PP[where,5]<-0
        }
        
      }
    }
    if (change==1){
      if (sum(ava_free>=1)<1){return(0)}
      whoV<-sample(seq(1,sum(vintageRecord[labelvintage==vintage,]>0)),1)
      whovalV<-vintageRecord[labelvintage==vintage,][whoV]
      
      where_add<-min(seq(1,ncol(vintageRecord),1)[vintageRecord[labelvintage==vintage,]==0])
      prob<-c(prob,1/sum(vintageRecord[labelvintage==vintage,]>0))
      
      vintageRecord[labelvintage==vintage,where_add]<-max(vintageRecord)+1
      vintageRecord[labelvintage==vintage,]<-c(vintageRecord[labelvintage==vintage,-whoV],0)
      # if (whovalV>1) {
      F[countIndex,1:(whovalV)]<-F[(countIndex-1),1:(whovalV)]-1
      F[countIndex,(whovalV+1):(countIndex-1)]<-F[(countIndex-1),(whovalV+1):(countIndex-1)]-2
      # }else{
      #   F[countIndex,1]<-F[countIndex-1,1]-1
      #   F[countIndex,(1+whovalV):(countIndex-1)]<-F[(countIndex-1),(1+whovalV):(countIndex-1)]-2
      # }
      countIndex<-countIndex+1
      #there is one free and one vintaged
      #Eliminate one free
      where<-min(seq(1,length(total_free))[potential$PP[,1]==vintage & total_free>0])
      total_free[where]<-total_free[where]-1
      ##From PP, I eliminate one "free" lineage
      PP<-potential$PP
      PP[where,3]<-PP[where,3]-1
      if (PP[where,3]==0){
        PP[where,5]<-0
      }
    }
    if (change==0){
      whoV<-sample(seq(1,sum(vintageRecord[labelvintage==vintage,]>0)),2)
      whoV1<-min(whoV)
      whoV2<-max(whoV)
      whovalV1<-vintageRecord[labelvintage==vintage,][whoV1]
      whovalV2<-vintageRecord[labelvintage==vintage,][whoV2]
      denom<-sum(vintageRecord[labelvintage==vintage,]>0)
      prob<-c(prob,2/(denom*(denom-1)))
      where_add<-min(seq(1,ncol(vintageRecord),1)[vintageRecord[labelvintage==vintage,]==0])
      vintageRecord[labelvintage==vintage,where_add]<-max(vintageRecord)+1
      vintageRecord[labelvintage==vintage,]<-c(vintageRecord[labelvintage==vintage,-whoV],0,0)
      F[countIndex,1:(whovalV1)]<-F[(countIndex-1),1:(whovalV1)]
      F[countIndex,(whovalV1+1):(whovalV2)]<-F[(countIndex-1),(whovalV1+1):(whovalV2)]-1
      F[countIndex,(whovalV2+1):(countIndex-1)]<-F[(countIndex-1),(whovalV2+1):(countIndex-1)]-2
      
      countIndex<-countIndex+1
      
      #just eliminate one that is not free from PP
      where<-min(seq(1,length(total_free))[potential$PP[,1]==vintage & potential$PP[,3]>0 & total_free==0])
      PP<-potential$PP
      PP[where,3]<-PP[where,3]-1
      if (PP[where,3]==0){
        PP[where,5]<-0
      }
    }
    
    if (sum(PP[,5])>1){
      check2_2<-check2(PP,sufficient,vintage,vintageRecord,labelvintage)
      PP<-check2_2$PP
      vintageRecord<-check2_2$vintageRecord
    }
    total_vintage<-c(total_vintage,vintage)
    change<-change_F[(countIndex-1)]
  }
  return(list(total_vintage=total_vintage[-1],change_F=change_F,prob=prob[-1],F=F))
}
tajima_info<-function(u,F){
  #returns basic info from the Tajima's data
  #u is the vector of intercoalescent times
  #given u and F, returns branch lengths ln+1,ln,ln-1,...,l3 and family sizes
  dimf<-nrow(F)
  diffM<-F[2:(dimf),2:dimf]-F[2:dimf,1:(dimf-1)]
  singletons<-seq(1,dimf-1)[diff(F[,1])<0]
  if (F[dimf,1]>0){singletons<-c(singletons,dimf)}
  dindicator<-matrix(0,nrow=dimf,ncol=dimf)
  d<-rep(u[1],dimf) #for corresponding ln+1,ln,ln-1,...,l3
  dindicator[1,1]<-1
  firstzero<-rep(0,dimf-1) 
  familysize<-rep(2,dimf)#for corresponding Sn,Sn-1,..S2
  clades_offspring<-0
  clades_parents<-0
  coal_times<-cumsum(u)
  for (j in 1:(dimf-1)){
    condition<-diffM[j:(dimf-1),j]==0
    if (sum(condition)>0){ firstzero[j]<-min(seq(j,(dimf-1))[condition])}else{firstzero[j]<-dimf}
    
    d[j+1]<-coal_times[firstzero[j]]-coal_times[j]
    dindicator[(j+1):firstzero[j],(j+1)]<-1
  }
  if (sum(dindicator[dimf,])==0){dindicator[dimf,dimf]<-1}
  d[dimf]<-u[dimf]
  firstzerotrans<-dimf+2-firstzero
  for (i in 1:(dimf-1)){
    count<-sum(firstzerotrans==(dimf-i+1))
    # if (count==0){
    #   familysize[i+1]<-2
    # }
    if (count==1){
      familysize[i+1]<-familysize[seq(1,dimf-1)[firstzerotrans==(dimf-i+1)]]+1
      clades_offspring<-c(clades_offspring,seq(1,dimf-1)[firstzerotrans==(dimf-i+1)])
      clades_parents<-c(clades_parents,i+1)
      
    }
    if (count==2){
      familysize[i+1]<-familysize[min(seq(1,(dimf-1))[firstzerotrans==(dimf-i+1)])]+familysize[max(seq(1,(dimf-1))[firstzerotrans==(dimf-i+1)])]
      clades_offspring<-c(clades_offspring,min(seq(1,(dimf-1))[firstzerotrans==(dimf-i+1)]),max(seq(1,(dimf-1))[firstzerotrans==(dimf-i+1)]))
      clades_parents<-c(clades_parents,i+1,i+1)
    }
    
  }
  #familysize[dimf]<-dimf+1
  #clades_offspring<-c(clades_offspring,max(clades_parents[clades_parents<dimf]))
  #clades_parents<-c(clades_parents,dimf)
  #familysize<-familysize[-dimf]
  #familysize<-c(1,familysize)
  return(list(d=d,dindicator=dindicator,familysize=familysize,nodes=cbind(clades_offspring[-1],clades_parents[-1]),singletons=singletons))
}

tajima_info2<-function(F){
  #returns basic info from the Tajima's data
  #u is the vector of intercoalescent times
  #given u and F, returns branch lengths ln+1,ln,ln-1,...,l3 and family sizes
  dimf<-nrow(F)
  diffM<-F[2:(dimf),2:dimf]-F[2:dimf,1:(dimf-1)]
  singletons<-seq(1,dimf-1)[diff(F[,1])<0]
  if (F[dimf,1]>0) {singletons<-c(singletons,dimf)}
  #d<-rep(u[1],dimf) #for corresponding ln+1,ln,ln-1,...,l3
  firstzero<-rep(0,dimf-1) 
  familysize<-rep(2,dimf)#for corresponding Sn,Sn-1,..S2
  clades_offspring<-0
  clades_parents<-0
  #coal_times<-cumsum(u)
  for (j in 1:(dimf-1)){
    condition<-diffM[j:(dimf-1),j]==0
    if (sum(condition)>0){ firstzero[j]<-min(seq(j,(dimf-1))[condition])}else{firstzero[j]<-dimf}
    
   # d[j+1]<-coal_times[firstzero[j]]-coal_times[j]
  }
  
#  d[dimf]<-u[dimf]
  firstzerotrans<-dimf+2-firstzero
  for (i in 1:(dimf-1)){
    count<-sum(firstzerotrans==(dimf-i+1))
    if (count==0){
      familysize[i+1]<-2
    }
    if (count==1){
      familysize[i+1]<-familysize[seq(1,dimf-1)[firstzerotrans==(dimf-i+1)]]+1
      clades_offspring<-c(clades_offspring,seq(1,dimf-1)[firstzerotrans==(dimf-i+1)])
      clades_parents<-c(clades_parents,i+1)
      
    }
    if (count==2){
      familysize[i+1]<-familysize[min(seq(1,(dimf-1))[firstzerotrans==(dimf-i+1)])]+familysize[max(seq(1,(dimf-1))[firstzerotrans==(dimf-i+1)])]
      clades_offspring<-c(clades_offspring,min(seq(1,(dimf-1))[firstzerotrans==(dimf-i+1)]),max(seq(1,(dimf-1))[firstzerotrans==(dimf-i+1)]))
      clades_parents<-c(clades_parents,i+1,i+1)
    }
    
  }
  nodes<-cbind(clades_offspring[-1],clades_parents[-1])
  ancestries<-matrix(0,nrow=max(nodes[,1]),ncol=length(familysize))
  j<-max(nodes[,1])
  #j<-nrow(tajima$nodes)
  while (j>0){
    where<-seq(1,nrow(nodes))[nodes[,1]==j]
    if (sum(nodes[,1]==j)>0){
      node<-nodes[where,2]
      ancestries[j,node]<-1
      if (node<ncol(ancestries)){
        ancestries[j,(node+1):ncol(ancestries)]<-ancestries[node,(node+1):ncol(ancestries)]
      }
    }
    j<-j-1
  }
 # familysize[dimf]<-dimf+1
#  clades_offspring<-c(clades_offspring,max(clades_parents[clades_parents<dimf]))
 # clades_parents<-c(clades_parents,dimf)
  #familysize<-familysize[-dimf]
  #familysize<-c(1,familysize)
  return(list(familysize=familysize,nodes=nodes,singletons=singletons,ancestries=ancestries))
}

Fstats<-function(F,t){
  nr<-nrow(F)
  Fdif<-rbind(cbind(F[2:nr,2:nr],rep(0,nr-1)),rep(0,nr))-rbind(F[2:nr,1:nr],rep(0,nr))
  
  wherezero<-rep(0,nr-1)
  #  wherezero2<-rep(0,nr-1)
  for (j in 1:(nr-1)){
    wherezero[j]<-nr+2-min(seq(j,nr)[Fdif[j:nr,j]==0])
    #   wherezero2[j]<-nr+2-max(seq(1,nr-j)[Fdif[nr-j+1,1:(nr-j+1)]==0])
  }
  #edge n,n-1,n-2,...,3
  branch.length<-t[11-wherezero]-t[-(nr)]
  branch.length.complement<-t[11-wherezero]
  #  wheresort<-table(sort(wherezero))
  family.size<-rep(2,nr)
  subtree.length<-2*t
  descendants<-matrix(0,nr,2)
  for (j in 2:nr){
    if (sum((wherezero==nr+2-j))==1){
      where<-seq(1:(nr-1))[wherezero==nr+2-j]
      family.size[j]<-1+family.size[where]
      subtree.length[j]<-subtree.length[where]+t[j]+branch.length[where]
      descendants[j,]<-c(11-where,0)
    }else{
      if (sum((wherezero==nr+2-j))==2){
        where<-seq(1:(nr-1))[wherezero==nr+2-j]
        family.size[j]<-family.size[where[1]]+family.size[where[2]]
        subtree.length[j]<-subtree.length[where[1]]+subtree.length[where[2]]+branch.length[where[1]]+branch.length[where[2]]
        descendants[j,]<-c(11-where[1],11-where[2])
      }
    }
  }
  return(list(family.size=family.size,subtree.length=subtree.length,branch.length=branch.length,branch.length.complement=branch.length.complement,descendants=descendants))
}

f_ancestries<-function(tajima){
  #Auxiliary function for Likelihood
 # ancestries<-matrix(0,nrow=nrow(tajima$nodes),ncol=length(tajima$familysize))
  ancestries<-matrix(0,nrow=max(tajima$nodes[,1]),ncol=length(tajima$familysize))
  j<-max(tajima$nodes[,1])
  #j<-nrow(tajima$nodes)
  while (j>0){
    where<-seq(1,nrow(tajima$nodes))[tajima$nodes[,1]==j]
    if (sum(tajima$nodes[,1]==j)>0){
    node<-tajima$nodes[where,2]
    ancestries[j,node]<-1
    if (node<ncol(ancestries)){
      ancestries[j,(node+1):ncol(ancestries)]<-ancestries[node,(node+1):ncol(ancestries)]
    }
    }
    j<-j-1
  }
  return(ancestries)
}


cond<-function(x,exclude){
  y<-rep(0,length(exclude))
  y[as.matrix(x)]<-1
  if (sum(y[exclude==1])>=1){return(1)}else{return(0)}
}

expand_grid_tajima<-function(table1,table2,ancestries){
  ne1<-nrow(table1) 
  nc1<-ncol(table1)
  ne2<-nrow(table2)
  nc2<-ncol(table2)
  result<-as.matrix(rep(0,ne1+ne2))
  for (j in 1:nc1){
    exclude<-rep(0,nrow(ancestries)+1)
    who<-table1[,j]
    exclude[who]<-1
    for (i in 1:ne1){
      exclude[c(ancestries[,who[i]]==1,FALSE)]<-1
    }
    rev<-apply(table2,2,cond,exclude)
    result<-cbind(result,rbind(matrix(rep(table1[,j],sum(rev==0)),nrow=ne1,ncol=sum(rev==0)),as.matrix(table2[,rev==0])))
    # if (j==1){
    #   result<-as.matrix(rbind(rep(table1[,j],sum(rev)),table2[,rev==0]))
    # }else{
    # result<-cbind(result,c(table1[,j],table2[,rev==0]))
    # }
  }
  return(result[,-1])
}


exponentiate_permutations<-function(x,exponents){
  X<-matrix(rep(x,nrow(exponents)),ncol=ncol(exponents),nrow=nrow(exponents),byrow=T)
  # if (ncol(x)!=ncol(exponents)){print("There is a problem with dimensions")}
  return(sum(exp(apply(log(X^exponents),1,sum))))
}

mini_expon<-function(y,x,exponents){
  exponents%*%(y/x)
}

exponentiate_permutations2<-function(x,exponents,indi){
  X<-matrix(rep(x,nrow(exponents)),ncol=ncol(exponents),nrow=nrow(exponents),byrow=T)
  weight<-apply(X^exponents,1,prod)
  values<-apply(indi,2,mini_expon,x=x,exponents=exponents)
  res<-weight%*%values/sum(weight)
  # if (ncol(x)!=ncol(exponents)){print("There is a problem with dimensions")}
  return(res)
}
#Ev<-E
#ind<-0
#suffdata<-suff_use
#offspring<-offspring_where


update<-function(E,parent,ancestors){
  
  values<-rep(0,ncol(ancestors))
  values[E$candidates]<-1
  ancestors[parent,]<-values
  return(list(ancestors=ancestors))  
}



expand.grid2<-function(assignment,permut){
  if (is.null(nrow(assignment))){
    if (is.null(nrow(permut))){
      return(c(assignment,permut))
    }else{
    return(cbind(matrix(rep(assignment,nrow(permut)),nrow=nrow(permut),ncol=length(assignment),byrow = T),permut))}
  }else{
    if (is.null(nrow(permut))){
      return(cbind(assignment,matrix(rep(permut,nrow(assignment)),nrow=nrow(assignment),ncol=length(permut),byrow=T)))
      
    }else{
    nrw<-nrow(assignment)
    assignment2<-cbind(matrix(rep(assignment[1,],nrow=nrow(permut)),nrow=nrow(permut),ncol=ncol(assignment),byrow=T),permut)
  for (j in 2:nrow(assignment)){
    assignment2<-rbind(assignment2,
                       cbind(matrix(rep(assignment[j,],nrow=nrow(permut)),nrow=nrow(permut),ncol=ncol(assignment),byrow=T),permut))
  }
    return(assignment2)
    }}
  
}


TajimaLik2<-function(u,suff,Finfo){
  bb_use<-tajima_info(u,Finfo$F)
  ##for given F with partial allocation
  #the only uncertainty is in leaves
  tab1<-table(suff[,1])
  #tab2<-table(suff[,2]*suff[,5]) #for leaves
  listnodes<-sort(unique(suff[,1]))
  #tabset<-listnodes[tab1>1]
  assignment<-1
  #labels<-1
  ##This works when there are not extra cherries
  #totrows<-1
  for (j in listnodes){
    nrep<-sum(suff[,1]==j)
    if (nrep==1){
      permut<-suff[suff[,1]==j,4]
    #  labels<-c(labels,j)
    }else{
    permut<-allPerm(initMC(suff[suff[,1]==j,4]))
    #labels<-c(labels,rep(j,nrep))
    }
    assignment<-expand.grid2(assignment,permut)
   
    #totrows<-totrows*nrow()
  }
  assignment<-assignment[,-1]
  labels<-sort(suff[,1])
  #if there is only one assingment then it is only one row
  ##assigments will have the different assigment of mutations
  #I only need to bring one row for the length of branch, the order corresponds to the label of the node
  branch_length<-rep(0,nrow(suff))
  ancestries<-f_ancestries(bb_use) #tells me the length of branches by finding the first one
  #lik<-0
  #singletons
  leaves<-unique(suff[suff[,5]==1 & suff[,3]==1,1]) ##only singletons
  for (j in 1:length(leaves)){
    parent<-unique(suff[suff[,1]==leaves[j],2])
    if (parent==0){
      time<-sum(u)
    }else{
      if (sum(suff[suff[,1]==leaves[j],5]==1)>1){
      ord1<-seq(1,length(Finfo$total_vintage))[Finfo$total_vintage==parent & Finfo$change_F==1]
      ord2<-seq(1,length(Finfo$total_vintage))[Finfo$total_vintage==parent & Finfo$change_F==2]
      time<-cumsum(u)[c(ord1,rep(ord2,2))]
      }else{
        ord<-seq(1,length(Finfo$total_vintage))[Finfo$total_vintage==parent & Finfo$change_F>0]
        time<-cumsum(u)[ord]
      }
      }
   branch_length[labels==leaves[j]]<-time
  }
  
  #cherries in data
  cherry<-unique(suff[suff[,5]==1 & suff[,3]==2,1])
  #cherry is an internal node
  # for (j in 1:length(cherry)){
  #     ord<-seq(1,length(Finfo$total_vintage))[Finfo$total_vintage==cherry[j] & Finfo$change_F==2]
  #     time<-cumsum(u)[ord]
  #     branch_length[labels==cherry[j]]<-time
  # }
  #internal nodes
  #internal<-unique(suff[suff[,5]==0,1]) ##only singletons
  internal<-suff[suff[,5]==0,1]
  desc<-suff[suff[,5]==0,3]
  internal<-c(cherry,internal)
  desc<-c(rep(2,length(cherry)),desc)
  famsize<-bb_use$familysize[-length(bb_use$familysize)]
  for (j in 1:length(internal)){
    ord<-seq(1,length(Finfo$total_vintage))[Finfo$total_vintage==internal[j]]
    upp<-min(seq(1,ncol(ancestries))[ancestries[ord,]==1])
    if (length(ord)>1){
    ord<-seq(1,length(Finfo$total_vintage))[Finfo$total_vintage==internal[j] & Finfo$change_F==0 & famsize==desc[j]]
    upp<-min(seq(1,ncol(ancestries))[ancestries[ord,]==1])
     }
    
    time<-cumsum(u)[upp]-cumsum(u)[ord]
    branch_length[labels==internal[j]]<-time
  }
  #I need to correct for the singletons when branch lenghts are the same
  ##correction
  repts<-table(labels)
  problem<-seq(1,length(repts))[repts>1]
  for (j in problem){
    temp_branch<-rep(0,length(branch_length))
    temp_branch[labels==j]<-branch_length[labels==j]
    assignment<-check_combs(assignment,temp_branch)
  }  
   
  #I need to repeat branch_lengths vectors with the permutations of 
  L<-sum(u*seq(length(u)+1,2))
  if (!is.null(nrow(assignment))){
  optional<-rep(0,nrow(assignment))
  for (j in 1:nrow(assignment)){
    optional[j]<-exp(sum(assignment[j,]*log(branch_length)))
    }
  }else{
    optional<-exp(sum(assignment*log(branch_length)))
  }
  likelihood<-log(sum(optional))-L #returns log likelihood
  return(likelihood)  
}  

check_combs<-function(exponents,Llengths){
  #I am sure there is a better way to do this but it will be manual for now
  #This is meant for singletons and cherries. The only way two lineages have the 
  #same length is for when 
  for (j in unique(Llengths[Llengths>0])){
    if (sum(Llengths==j)==2){ ##only cherries
      where<-seq(1,length(Llengths))[Llengths==j]
      exp<-apply(exponents[,where],1,sort)
      paris<-exponents
      paris[,where]<-t(exp)
      paris<-apply(paris,1,paste,collapse = "-")
      exponential2<-strsplit(unique(paris),"-")
      exponential3<-matrix(0,nrow=length(exponential2),ncol=length(exponential2[[1]]))
      for (i in 1:length(exponential2)){
        exponential3[i,]<-as.numeric(exponential2[[i]])
      }
    }
    exponents<-exponential3
  }
  return(exponents)
}


##tajima<-bb  
#Ev<-E
#suffdata<-suff_use
#offspring<-offspring_where


CalculateConditional2<-function(tajima,parent_where,offspring,Ev,suffdata,many,m,gd,path,ind){
 # for this purposes the cherrie counts as two singletons so I need to somehow bring back
#  old many 
 # times<-rep(0,length(Ev$times)+sum(Ev$newsizes==2 & Ev$flag==1))
  #times<-Ev$times
  if (ncol(Ev$times)==1){
    indicator<-rep(1,length(Ev$times))
    indicator[Ev$flag==1 & Ev$newsizes==2]<-2
    times<-as.matrix(Ev$times[rep(seq(1,nrow(Ev$times)),indicator)])
    allocations<-as.matrix(Ev$allocations[rep(seq(1,nrow(Ev$allocations)),indicator)])
  }else{
    
  indicator<-rep(1,nrow(Ev$times))
  indicator[Ev$flag==1 & Ev$newsizes==2]<-2
  times<-matrix(0,nrow=sum(indicator),ncol=ncol(Ev$times))
  allocations<-times
  flag<-times
  for (j in 1:ncol(Ev$times)){
    times[,j]<-Ev$times[rep(seq(1,nrow(Ev$times)),indicator),j]
    allocations[,j]<-Ev$allocations[rep(seq(1,nrow(Ev$allocations)),indicator),j]
    }
  }
  flag<-Ev$flag[rep(seq(1,length(Ev$flag)),indicator)]
  
  #times<-matrix(0,nrow=nrow(Ev$times)+sum(Ev$flag==1 & Ev$newsizes==2),ncol=ncol(Ev$times))
  #offspring<-Ev$offspring_where2
  #newmany<-rep(Ev$manyo,Ev$manyo)
  #newmany[Ev$newsizes==2 & Ev$flag==1]<-newmany[Ev$newsizes==1]
  start_c<-1
  end_c<-Ev$many[1]
  optional<-rep(0,length(Ev$candidates))
  for (i in 1:length(Ev$candidates)){
    start<-1
    end<-many[1] 
    factorsM<-matrix(0,nrow=length(offspring),ncol=Ev$many[i])
    probM<-matrix(0,nrow=length(offspring),ncol=Ev$many[i])
    factorsG<-matrix(0,nrow=length(offspring)*Ev$many,ncol=ncol(gd[[1]]))
    start2<-cumsum(c(0,rep(length(offspring),Ev$many[i]-1)))
    for (j in 1:length(offspring)){
      indicators<-matrix(0,nrow=Ev$many[i]*many[j],ncol=ncol(gd[[1]]))
      prev_factors<-matrix(0,nrow=many[j],ncol=Ev$many[i])
      #position<-rep(0,length=Ev$many[i]*many[j])
      for (k in 1:many[j]){
        for (l in 1:Ev$many[i]){
        #  position[(l-1)*many[j]+k]<-start2[l]+start+k-1
          if (flag[start+k-1]==1){indicators[(l-1)*many[j]+k,1:allocations[(start+k-1),start_c+l-1]]<-1}else{
            indicators[(l-1)*many[j]+k,]<-tajima$dindicator[,(1+allocations[(start+k-1),start_c+l-1])]
          }  
        }
        
        prev_factors[k,]<-m[offspring[j],allocations[(start+k-1),start_c:end_c]]
      }
      
      if (many[j]>1){
        m1 = initMC(suffdata$mylist[[offspring[j]]]$y)
        exponents<-allPerm(m1)
        #Note here I am only checking one of the columns, I am assuming the pattern is the same
        if (suffdata$nodes[offspring[j],3]==1 & nrow(as.matrix(times[start:end,start_c:start_c]))>length(unique(times[start:end,start_c:start_c]))){
          exponents<-check_combs(exponents,as.matrix(times[start:end,start_c:start_c]))
        }
        probM[j,]<-log(apply(as.matrix(times[start:end,start_c:end_c]),2,exponentiate_permutations,exponents))
        factorsM[j,]<-probM[j,]+apply(log(prev_factors),2,sum)
        for (l in 1:length(start2)){
          factorsG[start2[l]+j,]<-exponentiate_permutations2(times[start:end,l],exponents,indicators[(start2[l]+1):(start2[l]+many[j]),])
        }
        #factorsG[start2+j,]<-factorsG[start2+j,]+indicators/times[start:end,start_c:end_c]
      }else{
        mutations<-suffdata$mylist[[offspring[j]]]$y
        probM[j,]<-log(times[start:end,start_c:end_c]^mutations)
        factorsM[j,]<-probM[j,]+log(prev_factors)
        factorsG[start2+j,]<-factorsG[start2+j,]+indicators*(mutations/times[start:end,start_c:end_c])
        }
      start<-end+1
      end<-start+many[j+1]-1
    }
    where<-seq(1,length(path))[path==Ev$candidates[i]]
    weight<-exp(apply(probM,2,sum))
    weight<-weight/sum(weight)
    
    #if (Ev$many[i]>1){
     if (nrow(factorsG)>1 ){
       for (l in 1:Ev$many[i]){
        if (length(offspring)>1){
        gd[[where]][parent_where,]<-gd[[where]][parent_where,]+apply(factorsG[(start2[l]+1):(start2[l]+length(offspring)),],2,sum)*weight[l] 
        }else{
        gd[[where]][parent_where,]<-gd[[where]][parent_where,]+factorsG[(start2[l]+1):(start2[l]+1),]*weight[l] 
        }
       }    
    }else{
      gd[[where]][parent_where,]<-gd[[where]][parent_where,]+factorsG
    }
    
    resp<-sum(exp(apply(factorsM,2,sum)))
    if (ind==0){ m[parent_where,Ev$candidates[i]]<-resp}else{
      optional[i]<-resp
    }
    
    start_c<-end_c+1
    end_c<-start_c+Ev$many[i+1]-1
    
  }
  
  ##Here is where I will calculate all possible allocations when suffdata$mylist[[offspring[j]]]$x>1
  #Here I need to permute branches
  # result<-Ev$times^mutations[-1]
  #  gradient<-0
  return(list(m=m,optional=optional,gd=gd))
  
}
#offspring<-offspring_where
#Ev<-E
#suffdata<-suff_use
#ind<-0
# itional2<-function(parent_where,offspring,Ev,suffdata,many,m,ind){
#   start_c<-1
#   end_c<-Ev$many[1]
#   optional<-rep(0,length(Ev$candidates))
#   for (i in 1:length(Ev$candidates)){
#     start<-1
#     end<-many[1]
#     factorsM<-matrix(0,nrow=length(offspring),ncol=Ev$many[i])
#     
#     for (j in 1:length(offspring)){
#       prev_factors<-matrix(0,nrow=many[j],ncol=Ev$many[i])
#       for (k in 1:many[j]){
#         prev_factors[k,]<-m[offspring[j],Ev$allocations[(start+k-1),start_c:end_c]]
#       }
#       if (many[j]>1){
#         m1 = initMC(suffdata$mylist[[offspring[j]]]$y)
#         exponents<-allPerm(m1)
#         factorsM[j,]<-log(apply(as.matrix(Ev$times[start:end,start_c:end_c]),2,exponentiate_permutations,exponents))+apply(log(prev_factors),2,sum)
#       }else{
#         mutations<-suffdata$mylist[[offspring[j]]]$y
#         factorsM[j,]<-log(Ev$times[start:end,start_c:end_c]^mutations)+log(prev_factors)
#       }
#       start<-end+1
#       end<-start+many[j+1]-1
#     }
#     if (ind==0){ m[parent_where,Ev$candidates[i]]<-sum(exp(apply(factorsM,2,sum)))}else{
#       optional[i]<-sum(exp(apply(factorsM,2,sum)))
#     }
#     
#     start_c<-end_c+1
#     end_c<-start_c+Ev$many[i+1]-1
#     
#   }
#   
#   ##Here is where I will calculate all possible allocations when suffdata$mylist[[offspring[j]]]$x>1
#   #Here I need to permute branches
#   # result<-Ev$times^mutations[-1]
#   #  gradient<-0
#   return(list(m=m,optional=optional))
#   
# }
#offspring<-offspring_where
#Ev<-E
#suffdata<-suff_use
#ind<-0
#tajima<-bb

CalculateConditional<-function(tajima,parent_where,offspring,Ev,suffdata,many,m,gd,path,ind){
  start_c<-1
  end_c<-Ev$many[1]
  optional<-rep(0,length(Ev$candidates))
  for (i in 1:length(Ev$candidates)){
    start<-1
    end<-many[1]
    factorsM<-matrix(0,nrow=length(offspring),ncol=Ev$many[i])
    probM<-matrix(0,nrow=length(offspring),ncol=Ev$many[i])
    factorsG<-matrix(0,nrow=length(offspring)*Ev$many,ncol=ncol(gd[[1]]))
    start2<-cumsum(c(0,rep(length(offspring),Ev$many[i]-1)))
    
    for (j in 1:length(offspring)){
      indicators<-matrix(0,nrow=Ev$many[i]*many[j],ncol=ncol(gd[[1]]))
      prev_factors<-matrix(0,nrow=many[j],ncol=Ev$many[i])
      for (k in 1:many[j]){
        for (l in 1:Ev$many[i]){
        if (Ev$flag[start+k-1]==1){indicators[(l-1)*many[j]+k,1:Ev$allocations[(start+k-1),start_c+l-1]]<-1}else{
          indicators[(l-1)*many[j]+k,]<-tajima$dindicator[,(1+Ev$allocations[(start+k-1),start_c+l-1])]
        }
        }
        
        prev_factors[k,]<-m[offspring[j],Ev$allocations[(start+k-1),start_c:end_c]]
      }
      if (many[j]>1){
        m1 = initMC(suffdata$mylist[[offspring[j]]]$y)
        exponents<-allPerm(m1)
        if (suffdata$nodes[offspring[j],3]==1 & nrow(as.matrix(Ev$times[start:end,start_c:end_c]))>length(unique(Ev$times[start:end,start_c:end_c]))){
          exponents<-check_combs(exponents,as.matrix(Ev$times[start:end,start_c:end_c]))
        }
        probM[j,]<-log(apply(as.matrix(Ev$times[start:end,start_c:end_c]),2,exponentiate_permutations,exponents))
        factorsM[j,]<-probM[j,]+apply(log(prev_factors),2,sum)
        for (l in 1:length(start2)){
          factorsG[start2[l]+j,]<-exponentiate_permutations2(Ev$times[start:end,l],exponents,indicators[(start2[l]+1):(start2[l]+many[j]),])
        }
        #gr[Ev$candidates[i],]<-gr[Ev$candidates[i],]+mutations*(Ev$times[start:end,start_c:end_c]^(mutations-1))*indicators
        
        }else{
        mutations<-suffdata$mylist[[offspring[j]]]$y
        probM[j,]<-log(Ev$times[start:end,start_c:end_c]^mutations)
        factorsM[j,]<-probM[j,]+log(prev_factors)
        factorsG[start2+j,]<-factorsG[start2+j,]+indicators*(mutations/Ev$times[start:end,start_c:end_c])
        
       # factorsM[j,]<-log(Ev$times[start:end,start_c:end_c]^mutations)+log(prev_factors)
        #where<-seq(1,length(path))[path==Ev$candidates[i]]
        #gd[[where]][parent_where,]<-gd[[where]][parent_where,]+(mutations/Ev$times[start:end,start_c:end_c])*indicators
      }
      start<-end+1
      end<-start+many[j+1]-1
    }
    where<-seq(1,length(path))[path==Ev$candidates[i]]
    weight<-exp(apply(probM,2,sum))
    weight<-weight/sum(weight)
    if (nrow(factorsG)>1){
    for (l in 1:Ev$many[i]){
      if (length(offspring)>1){
      gd[[where]][parent_where,]<-gd[[where]][parent_where,]+apply(factorsG[(start2[l]+1):(start2[l]+length(offspring)),],2,sum)*weight[l] 
      }else{
        gd[[where]][parent_where,]<-gd[[where]][parent_where,]+factorsG[(start2[l]+1):(start2[l]+length(offspring)),]*weight[l] 
        
      }
      }
    }else{
      gd[[where]][parent_where,]<-gd[[where]][parent_where,]+factorsG
    }
    if (ind==0){ m[parent_where,Ev$candidates[i]]<-sum(exp(apply(factorsM,2,sum)))}else{
      optional[i]<-sum(exp(apply(factorsM,2,sum)))
    }
    
    start_c<-end_c+1
    end_c<-start_c+Ev$many[i+1]-1
    
  }
  
  ##Here is where I will calculate all possible allocations when suffdata$mylist[[offspring[j]]]$x>1
  #Here I need to permute branches
  # result<-Ev$times^mutations[-1]
  #  gradient<-0
  return(list(m=m,optional=optional,gd=gd))
  
}
#parent<-parent_where
update<-function(E,parent,ancestors){
  
  values<-rep(0,ncol(ancestors))
  values[E$candidates]<-1
  ancestors[parent,]<-values
  return(list(ancestors=ancestors))  
}


generateFprop<-function(Finfo){
  #Update: March 17, 2017
  ##Creates an actual F matrix from the output of Finitial$change_F
  ##It doesn't take into account Finitial$total_vintage
  #that already includes partial raking within families
  #generateFprop2 will create an actual F matrix considering the rankings
  n<-sum(Finfo$change_F)
  F<-matrix(0,nrow=n-1,ncol=n-1)
  diag(F)<-seq(n,2)
  prob<-0
  F[2:(n-1),1]<-n-cumsum(Finfo$change_F)[-(n-1)]
  vintages<-1
  for (j in 3:(n-1)){
    if (Finfo$change_F[j-1]==2){F[j,2:(j-1)]<-F[j-1,2:(j-1)]-2; vintages<-c(vintages,max(vintages)+1)}
    if (Finfo$change_F[j-1]==1){
      who<-sample(seq(1,length(vintages)),1)
      whoval<-vintages[who]
      prob<-c(prob,1/length(vintages))
      vintages<-c(vintages,max(vintages)+1)
      vintages<-vintages[-who]
      if (whoval>1) {
        F[j,2:(whoval)]<-F[(j-1),2:(whoval)]-1
        F[j,(whoval+1):(j-1)]<-F[(j-1),(whoval+1):(j-1)]-2
        }else{
        F[j,(1+whoval):(j-1)]<-F[(j-1),(1+whoval):(j-1)]-2
      }
    }
    if (Finfo$change_F[j-1]==0){
      who<-sample(seq(1,length(vintages)),2)
      who1<-min(who)
      who2<-max(who)
      whoval1<-vintages[who1]
      whoval2<-vintages[who2]
      denom<-length(vintages)
      prob<-c(prob,2/(denom*(denom-1)))
      vintages<-c(vintages,max(vintages)+1)
      vintages<-vintages[-who]
      if (whoval1>1) {
        F[j,2:(whoval1)]<-F[(j-1),2:(whoval1)]
        F[j,(whoval1+1):(whoval2)]<-F[(j-1),(whoval1+1):(whoval2)]-1
        F[j,(whoval2+1):(j-1)]<-F[(j-1),(whoval2+1):(j-1)]-2
      }else{
        F[j,(whoval1+1):(whoval2)]<-F[(j-1),(whoval1+1):(whoval2)]-1
        F[j,(whoval2+1):(j-1)]<-F[(j-1),(whoval2+1):(j-1)]-2
      }
    }
  }
  F[F<0]<-0
  return(list(F=F,vintages=vintages[-1],prob=prob[-1]))
}


#' Simulate data according to infinite sites mutation model
#' 
#' @param mu overall mutation rate. The expected number of mutations
#'  will be the (mu)(L) where L is the tree length 
#' @param tree1 ape tree object
#'   
#' @export
simulate_data<-function(mu,tree1){
  n<-tree1$Nnode+1
  u<-coalescent.intervals(tree1)$interval.length
  tol<-.0000001
  u[u==0L]<-tol
  L<-sum(u*seq(n,2)) #total tree length
  # mu<-20 ##2*N*l, l is the number of loci
  mut<-rpois(1,L*mu)
  where<-runif(mut,0,L)
  wheresort<-sort(where)
  stretchedtree<-cumsum(rep(u,seq(n,2)))
  indicator<-rep(seq(n,2),seq(n,2))
  ind2<-matrix(0,nrow=n,ncol=n-1)
  ind3<-matrix(0,nrow=n,ncol=n-1)
  for (i in 1:(n-1)){
    node<-2*n-i
    #(parent,offspring)
    children<-tree1$edge[tree1$edge[,1]==node,2]
    for (k in children){
      if (k<=n){
        ind2[k,i]<-1
      }else{
        ind2[,i]<-ind2[,i]+ind2[,2*n-k]
      }
    }
    if (i==1){ind3[,i]<-seq(1,nrow(ind2))-cumsum(ind2[,i])+(cumsum(ind2[,i])>0)}
    else{
      val<-min(ind3[ind2[,i]==1,i-1])
      ind3[,i]<-ind3[,i-1]
      ind3[ind2[,i]==1,i]<-val
      if (sum(diff(ind3[,i])>1)>0){
        val<-diff(ind3[,i])[diff(ind3[,i])>1]
        where2<-seq(1,nrow(ind3)-1)[diff(ind3[,i])>1]+1
        ind3[where2:nrow(ind3),i]<-ind3[where2:nrow(ind3),i]-val+1
      }
    }
  }
  
  nbranches<-length(indicator)
  data<-matrix(0,nrow=n,ncol=mut)
  for (j in 1:mut){
    if (wheresort[j]<min(stretchedtree)) {who<-1}else{
      who<-max(seq(1,nbranches)[stretchedtree<wheresort[j]])+1}
    if (who<=n){data[who,j]<-1}
    if (who>n){
      epoch<-n-indicator[who]
      branch_order<-sum(indicator[1:who]==indicator[who])
      data[seq(1,n)[ind3[,epoch]==branch_order],j]<-1
    }
    if (j==1){data1<-as.character(paste(data[,j],collapse = '', sep = ''))}
    else{data1<-rbind(data1,as.character(paste(data[,j],collapse = '', sep = '')))}
  }
  
  return(data1)
  
}


getDistributionF0<-function(sufficient,Nsim){
  #update May 2
  #Function to validate marginal distribution of F
  sample1<-F_iter3(sufficient)$F
  Flist<-list(sample1)
  freqF<-1
  for (j in 1:(Nsim-1)){
    sample1<-F_iter3(sufficient)$F
    exists<-0
    for (i in 1:length(Flist)){
      if (sum(abs(Flist[[i]]-sample1))!=0) {
        exists<-exists+1
      }else{
        freqF[i]<-freqF[i]+1; break;
      }
      if (exists==length(Flist)){
        Flist[[length(Flist)+1]]<-sample1
        freqF<-c(freqF,1)
      }
    }
  }
  return(list(Flist=Flist,freqF=freqF))
}

coal_F<-function(F){
  #returns coalescent prior
  cherries<-sum(diff(F[,1])==-2)
  return(2^(F[1,1]-cherries-1)/factorial(F[1,1]-1))
}

getDistributionFA<-function(suf_ext,Nsim,Finfo=NULL){
  #update December 2017
  #Function to validate marginal distribution of F and allocation
  if (is.null(Finfo)){
    samp<-F_sample(suf_ext)
  }else{
    samp<-Finfo
  }
  samp
  sample1<-samp$F
  sample2<-samp$total_vintage
  Flist<-list(sample1)
  Flist2<-list(sample2)
  Flist3<-prod(samp$prob)
  freqF<-1
  for (j in 1:(Nsim-1)){
    samp<-F_sample(suf_ext)
    sample1<-samp$F
    sample2<-samp$total_vintage
    exists<-0
    for (i in 1:length(Flist)){
      if ((sum(abs(Flist[[i]]-sample1)>0)!=0) || (sum(abs(Flist2[[i]]-sample2))!=0)) { #The F matrixx is the same
        #if (sum(abs(Flist[[i]]-sample1)>0)!=0 || sum(abs(Flist2[[i]]-sample2))!=0) {
        exists<-exists+1
      }else{
        freqF[i]<-freqF[i]+1;
      }
      if (exists==length(Flist)){
        Flist[[length(Flist)+1]]<-sample1
        Flist2[[length(Flist2)+1]]<-sample2
        Flist3<-c(Flist3,prod(samp$prob))
        freqF<-c(freqF,1)
      }
    }
  }
  return(list(Flist=Flist,freqF=freqF,Flist2=Flist2,Flist3=Flist3))
  #return(list(Flist=Flist,Flist2=Flist2,Flist3=Flist3,freqF=freqF))
}

getDistributionPA<-function(suf_ext,Nsim,Finfo=NULL){
  #update December 2017
  #Function to validate marginal distribution of F and allocation
  if (is.null(Finfo)){
    samp<-python.call("F_sample", suf_ext)
  }else{
    samp<-Finfo
  }
  sample1<-matrix(unlist(samp$F_mat),nrow=length(samp$F_mat[[1]]),byrow=TRUE)
  sample2<-samp$total_vintage
  Flist<-list(sample1)
  Flist2<-list(sample2)
  Flist3<-prod(samp$probs)
  freqF<-1
  for (j in 1:(Nsim-1)){
    samp<-python.call("F_sample", suf_ext)
    sample1<-matrix(unlist(samp$F_mat),nrow=length(samp$F_mat[[1]]),byrow=TRUE)
    sample2<-samp$total_vintage
    exists<-0
    for (i in 1:length(Flist)){
      if ((sum(abs(Flist[[i]]-sample1)>0)!=0) || (sum(abs(Flist2[[i]]-sample2))!=0)) { #The F matrixx is the same
        #if (sum(abs(Flist[[i]]-sample1)>0)!=0 || sum(abs(Flist2[[i]]-sample2))!=0) {
        exists<-exists+1
      }else{
        freqF[i]<-freqF[i]+1;
      }
      if (exists==length(Flist)){
        Flist[[length(Flist)+1]]<-sample1
        Flist2[[length(Flist2)+1]]<-sample2
        Flist3<-c(Flist3,prod(samp$probs))
        freqF<-c(freqF,1)
      }
    }
  }
  return(list(Flist=Flist,freqF=freqF,Flist2=Flist2,Flist3=Flist3))
  #return(list(Flist=Flist,Flist2=Flist2,Flist3=Flist3,freqF=freqF))
}
getDistributionFAR<-function(sufficient,Nsim,Finfo=NULL){
  #update May 2
  #Function to validate marginal distribution of F and allocation
  if (is.null(Finfo)){
    samp<-F_iter5(sufficient)
  }else{
    samp<-Finfo
  }
  sample1<-samp$F
  sample2<-samp$total_vintage
  Flist<-list(sample1)
  Flist2<-list(sample2)
  Flist3<-prod(samp$prob)
  freqF<-1
  for (j in 1:(Nsim-1)){
    samp<-F_iter5(sufficient)
    sample1<-samp$F
    sample2<-samp$total_vintage
    exists<-0
    for (i in 1:length(Flist)){
    #  if (sum(abs(Flist[[i]]-sample1)>0)!=0) {
        if (sum(abs(Flist[[i]]-sample1)>0)!=0 || sum(abs(Flist2[[i]]-sample2)>0)!=0) {
        exists<-exists+1
      }else{
        freqF[i]<-freqF[i]+1; break;
      }
      if (exists==length(Flist)){
        Flist[[length(Flist)+1]]<-sample1
             Flist2[[length(Flist2)+1]]<-sample2
            Flist3<-c(Flist3,prod(samp$prob))
        freqF<-c(freqF,1)
      }
    }
  }
 # return(list(Flist=Flist,freqF=freqF))
  return(list(Flist=Flist,Flist2=Flist2,Flist3=Flist3,freqF=freqF))
}

getDistributionFA2<-function(sufficient,Nsim,Finfo){
  #update May 2 -- Not efficient
  #I want to the the distribution of allocations for a fixed F
  #Function to validate marginal distribution of F and allocation
    #samp<-F_iter3(sufficient)
  #  F<-Finfo$F
  F<-Finfo$F
  sample2<-Finfo$total_vintage
  #Flist<-list(sample1)
  Flist2<-list(sample2)
  #Flist3<-list(prod(samp$prob))
  freqF<-1
  for (j in 1:(Nsim-1)){
    samp<-F_iter3(sufficient)
    if (sum(abs(samp$F-F))==0){
      sample2<-samp$total_vintage
      exists<-0   
      for (i in 1:length(Flist2)){
        if (sum(abs(Flist2[[i]]-sample2))!=0) {
          exists<-exists+1
        }else{
          freqF[i]<-freqF[i]+1; break;
        }
        if (exists==length(Flist2)){
          Flist2[[length(Flist2)+1]]<-sample2
          freqF<-c(freqF,1)
        }
      }
    }
  }
  return(list(Flist2=Flist2,freqF=freqF))
}

getDistributionFA4<-function(sufficient,Nsim,Finfo){
  #samples with same F but different allocations - it is not efficient 3% acceptance with same Fs
  #I want to the the distribution of allocations for a fixed F
  #Function to validate marginal distribution of F and allocation
  #samp<-F_iter3(sufficient)
  #  F<-Finfo$F
  F<-Finfo$F
  sample2<-Finfo$total_vintage
  #Flist<-list(sample1)
  Flist2<-list(sample2)
  #Flist3<-list(prod(samp$prob))
  freqF<-1
  for (j in 1:(Nsim-1)){
    samp<-F_iter4(sufficient,Finfo$change_F)
    if (sum(abs(samp$F-F)>0)==0){
      sample2<-samp$total_vintage
      exists<-0   
      for (i in 1:length(Flist2)){
        if (sum(abs(Flist2[[i]]-sample2)>0)!=0) {
          exists<-exists+1
        }else{
          freqF[i]<-freqF[i]+1; break;
        }
        if (exists==length(Flist2)){
          Flist2[[length(Flist2)+1]]<-sample2
          freqF<-c(freqF,1)
        }
      }
    }
  }
  return(list(Flist2=Flist2,freqF=freqF))
}


getDistributionFA5<-function(sufficient,Nsim,Finfo){
  #samples with same F but different allocations - it is not efficient 3% acceptance with same Fs
  #I want to the the distribution of allocations for a fixed F
  #Function to validate marginal distribution of F and allocation
  #samp<-F_iter3(sufficient)
  #  F<-Finfo$F
  F<-Finfo$F
  sample2<-Finfo$total_vintage
  #Flist<-list(sample1)
  Flist2<-list(sample2)
  #Flist3<-list(prod(samp$prob))
  freqF<-1
  for (j in 1:(Nsim-1)){
    samp<-F_iter6(sufficient,Finfo$change_F)
    samp
    if (length(samp)>1){
    if (sum(abs(samp$F-F)>0)==0){
      sample2<-samp$total_vintage
      exists<-0   
      for (i in 1:length(Flist2)){
        if (sum(abs(Flist2[[i]]-sample2)>0)!=0) {
          exists<-exists+1
        }else{
          freqF[i]<-freqF[i]+1; break;
        }
        if (exists==length(Flist2)){
          Flist2[[length(Flist2)+1]]<-sample2
          freqF<-c(freqF,1)
        }
      }
    }
    }
  }
  return(list(Flist2=Flist2,freqF=freqF))
}


getDistributionF1<-function(sufficient,Nsim){
  #update April 29
  #Function to validate distribution of F first row given the data
  sample1<-F_iter2(sufficient)$change_F
  Flist<-list(sample1)
  freqF<-1
  for (j in 1:(Nsim-1)){
    sample1<-F_iter2(sufficient)$change_F
    exists<-0
    for (i in 1:length(Flist)){
      if (sum(abs(Flist[[i]]-sample1))!=0) {
        exists<-exists+1
      }else{
        freqF[i]<-freqF[i]+1; break;
      }
      if (exists==length(Flist)){
        Flist[[length(Flist)+1]]<-sample1
        freqF<-c(freqF,1)
      }
    }
  }
  return(list(Flist=Flist,freqF=freqF))
}
  

getDistributionF2<-function(sufficient,Nsim){
  #update April 29
  #Function to validate distribution of F first row given the data
  sample1<-F_iter2(sufficient)$change_F
  Flist<-list(sample1)
  freqF<-1
  for (j in 1:(Nsim-1)){
    sample1<-F_iter2(sufficient)$change_F
    exists<-0
    for (i in 1:length(Flist)){
      if (sum(abs(Flist[[i]]-sample1))!=0) {
        exists<-exists+1
      }else{
        freqF[i]<-freqF[i]+1; break;
      }
      if (exists==length(Flist)){
        Flist[[length(Flist)+1]]<-sample1
        freqF<-c(freqF,1)
      }
    }
  }
  return(list(Flist=Flist,freqF=freqF))
}

getDistribution<-function(Finfo,Nsim){
  #Update: April 29. Samples all Fs given the last row of F
  #Validation function
  F<-generateFprop(Finfo)
  Flist<-list(F$F)
  freqF<-1
  for (j in 1:(Nsim-1)){
    F<-generateFprop(Finfo)
    exists<-0
    for (i in 1:length(Flist)){
      if (sum(abs(Flist[[i]]-F$F))!=0) {
        exists<-exists+1
      }else{
        freqF[i]<-freqF[i]+1; break;
      }
      if (exists==length(Flist)){
        Flist[[length(Flist)+1]]<-F$F
        freqF<-c(freqF,1)
      }
    }
  }
  return(list(Flist=Flist,freqF=freqF))
  }
  
newickFmatrix<-function(F,br){
  #generates a newick format of an F matrix
  #x are intercoalescent times
  n <- F[1,1]
  D<-matrix(0,n-1,n-1)
  D[1:(n-2),1:(n-1)]<-F[1:(n-2),]-F[2:(n-1),]
  D[D<0]<-0
  D[n-1,]<-F[n-1,]
  nbr <- 2 * n - 2
  x <- if (is.character(br)) 
    (2 * rexp(n - 1)/(as.double(n:2) * as.double((n - 1):1)))
  else if (is.numeric(br)) 
    rep(br, length.out = n - 1)
  x<-c(0,x)
  vintage<-paste("(1:",x[2],",2:",x[2],")",sep="")
  labels<-seq(3,n,1)
  for (j in 3:(n-1)){
    if (F[j,1]-F[(j-1),1]==-2){
      vintage<-c(vintage,paste("(",labels[1],":",cumsum(x)[j],",",labels[2],":",cumsum(x)[j],")",sep=""))
      labels<-labels[-c(1,2)]
    }
    if (F[j,1]-F[(j-1),1]==-1){
      where=min(seq(1,n-1)[D[j-1,]==2])
      vintage<-c(vintage,paste("(",vintage[where-1],":",cumsum(x)[j]-cumsum(x)[where],",",labels[1],":",cumsum(x)[j],")",sep=""))
      labels<-labels[-1]
      #vintage<-vintage[-(where-1)]
    }
    if (F[j,1]-F[(j-1),1]==0){
      where2=min(seq(1,n-1)[D[j-1,]==2])
      where1=min(seq(1,n-1)[D[j-1,]==1])
      vintage<-c(vintage,paste("(",vintage[where1-1],":",cumsum(x)[j]-cumsum(x)[where1],",",vintage[where2-1],":",cumsum(x)[j]-cumsum(x)[where2],")",sep=""))
      
    }
  }
  j<-j+1
  where2=min(seq(1,n-1)[D[j-1,]==2])
  where1=min(seq(1,n-1)[D[j-1,]==1])
  if (where1==1){
    
    vintage<-c(vintage,paste("(",vintage[where2-1],":",cumsum(x)[j]-cumsum(x)[where2],",",labels[1],":",cumsum(x)[j],");",sep=""))
  # vintage<-c(vintage,paste("(",vintage[where2-1],":",cumsum(x)[j],",",labels[1],":",cumsum(x)[j],");",sep=""))
  }else{
    vintage<-c(vintage,paste("(",vintage[where1-1],":",cumsum(x)[j]-cumsum(x)[where1],",",vintage[where2-1],":",cumsum(x)[j]-cumsum(x)[where2],");",sep=""))
  }
  return(vintage[n-1])
}


#parent<-parent_size
#offspring<-offspring_sizes
#tajima<-bb
#
##I am here.
evidence<-function(parent,offspring,many,tajima,u,ancestors,ancestries,offspring_where,che=0){
  ##Bring the candidates with the same parent size
  candidates<-seq(1,length(tajima$familysize))[tajima$familysize==parent]
  valid_candidate<-rep(1,length(candidates))
  #Are those candidate nodes in the Tajima tree parents of the nodes that correspond to offspring from the data?
  
  
  ##Check constrains on each candidate, each candidate have to have offsprings in ancestors
  if (parent>2){
    for (j in 1:length(candidates)){
      descendents<-seq(1,nrow(ancestries))[ancestries[,candidates[j]]==1]
      condition<-0
      for (i in 1:length(offspring)){
        condition<-condition+ifelse(sum(ancestors[offspring_where[i],descendents])>=1, 1, 0)
      }
      if (condition<length(offspring)){valid_candidate[j]<-0}
    }
    candidates<-candidates[valid_candidate==1]
  }
  
  many_p<-rep(0,length(candidates))
  #many<-rep(1,length(offspring))
  allocations<-as.matrix(rep(0,sum(many)))
  for (i in 1:length(candidates)){
    potential<-seq(1,nrow(ancestries))[ancestries[,candidates[i]]==1]
    total_candidates<-vector("list",length(offspring))
    total_combinations<-vector("list",length(offspring))
    
    for (j in 1:length(offspring)){
      if (offspring[j]==1){
        posval1<-rep(0,candidates[i])
        posval2<-rep(0,candidates[i])
        potential<-c(potential,candidates[i])
        posval1[potential]<-1
        sing<-tajima$singletons[tajima$singletons<=candidates[i]]
        
        posval2[sing]<-1
        total_candidates[[j]]<-seq(1,candidates[i])[posval1==1 & posval2==1]
        #if (length(total_candidates[[j]])<many[j]){
          #find the cherries that need to be used as singletons only if parent is greater than 2
          if (parent>2){total_candidates[[j]]<-sort(c(total_candidates[[j]],total_candidates[[j]][tajima$familysize[total_candidates[[j]]]==2]))}
          
        #}
        if (length(total_candidates[[j]])==1){total_combinations[[j]]<-as.matrix(total_candidates[[j]])}else{
          total_combinations[[j]]<-combn(total_candidates[[j]],many[j])}
          ##correction for the cases the two branches of the cherries are not selected 
      #    cherr<-unique(total_candidates[[j]])[table(total_candidates[[j]])==2]
      #    if (many[j]>1){
      #      mylittle<-function(x,cc){return(sum(x==cc))}
      #    for (cc in cherr){
       #     total_combinations[[j]]<-total_combinations[[j]][,-seq(1,ncol(total_combinations[[j]]))[apply(total_combinations[[j]],2,mylittle,cc=cc)==1]]
      #    }}
      }else{
        sizes<-tajima$familysize[potential]
        cand2<-potential[sizes==offspring[j]]
        cand2<-cand2[ancestors[offspring_where[j],][cand2]==1]
        total_candidates[[j]]<-cand2
        if (length(cand2)==1){total_combinations[[j]]<-as.matrix(cand2)}else{
          total_combinations[[j]]<-combn(cand2,many[j])}
      }
    }
    ##Here is where we find all possible allocations
    if (length(offspring)>1){
      #if (ncol(total_combinations[[1]])<ncol(total_combinations[[2]])){
      #  allocations_1<-expand_grid_tajima(total_combinations[[2]],total_combinations[[1]],ancestries)
       # }else{
       allocations_1<-expand_grid_tajima(total_combinations[[1]],total_combinations[[2]],ancestries)
      #}
     # allocations_1<-expand_grid_tajima(total_combinations[[1]],total_combinations[[2]],ancestries)
      if (length(offspring)>2){
        for (j in 3:length(offspring)){
          #if (ncol(as.matrix(allocations_1))<ncol(total_combinations[[j]])){
           # allocations_1<-expand_grid_tajima(total_combinations[[j]],as.matrix(allocations_1),ancestries)}
            #else{
              allocations_1<-expand_grid_tajima(as.matrix(allocations_1),total_combinations[[j]],ancestries)
              }
      }
      }else{
    #  allocations_1<-cbind(rep(0,sum(many)),as.matrix(total_combinations[[1]]))
      allocations_1<-as.matrix(total_combinations[[1]])
    }
   #allocations<-allocations_1 
    allocations<-cbind(allocations,as.matrix(allocations_1))
    many_p[i]<-ncol(as.matrix(allocations_1))  
    
  }
  allocations<-as.matrix(allocations[,-1])
  if (sum(many)==1){allocations<-t(as.matrix(allocations))}
  allocations<-matrix(allocations,nrow=sum(many),ncol=sum(many_p))
  times<-matrix(0,nrow=sum(many),ncol=sum(many_p))
  flag<-rep(0,sum(many))
  if (sum(offspring==1)==1){
    cc<-sum(many[offspring>1])+1
    flag[cc:(cc+many[offspring==1]-1)]<-1
  }
  times[flag==0,]<-tajima$d[-1][allocations[flag==0,]]
  times[flag==1,]<-cumsum(u)[allocations[flag==1,]]
  
  #the list has the following hierarchy: candidates, offspring (descending order by size), many times
  return(list(candidates=candidates,many=many_p,times=times,allocations=allocations,flag=flag))
}


#parent<-parent_size
#offspring<-offspring_sizes
#tajima<-bb
#
##I am here.
evidence3<-function(parent,offspring,many,tajima,u,ancestors,ancestries,offspring_where,che=0){
  #special for when extra cherries
  
  ##Bring the candidates with the same parent size
  candidates<-seq(1,length(tajima$familysize))[tajima$familysize==parent]
  valid_candidate<-rep(1,length(candidates))
  #Are those candidate nodes in the Tajima tree parents of the nodes that correspond to offspring from the data?
  
  
  ##Check constrains on each candidate, each candidate have to have offsprings in ancestors
  if (parent>2){
    for (j in 1:length(candidates)){
      descendents<-seq(1,nrow(ancestries))[ancestries[,candidates[j]]==1]
      condition<-0
      for (i in 1:length(offspring)){
        condition<-condition+ifelse(sum(ancestors[offspring_where[i],descendents])>=1, 1, 0)
      }
      if (condition<length(offspring)){valid_candidate[j]<-0}
    }
    candidates<-candidates[valid_candidate==1]
  }
  ##Cherries needed in candidates?
  gost_cherries<-sum(tajima$familysize==2)
  available_sing<-floor(many[offspring==1]/2)
  if (sum(many[offspring==2])>0){num=many[offspring==2]}else{num=0}
  original2<-num
  if (gost_cherries>(num+available_sing) & num>0)
  { #then I can create up to available_sing

    many[offspring==2]<-original2+available_sing
    many[offspring==1]<-many[offspring==1]-available_sing*2
  }
 # permut<-factorial(many[offspring==2]) #the number of column
  flag<-rep(0,sum(many))
  offspring_where2<-rep(0,sum(many))
  if (sum(offspring==1)==1){
    cc<-sum(many[offspring>1])+1
    flag[cc:(cc+many[offspring==1]-1)]<-1
  }
  cc<-sum(many[offspring>2])+original2+1
  flag[cc:(cc+available_sing)]<-1
  newsizes<-rep(offspring,many)
 
  offspring_where2[newsizes>2]<-offspring_where[offspring>2]
  offspring_where2[newsizes==2 & flag==0]<-offspring_where[offspring==2]
  offspring_where2[newsizes==2 & flag==1]<-offspring_where[offspring==1]
  offspring_where2[newsizes==1 & flag==1]<-offspring_where[offspring==1]
  many_p<-rep(0,length(candidates))
  #many<-rep(1,length(offspring))
  allocations<-as.matrix(rep(0,sum(many)))
  for (i in 1:length(candidates)){
    potential<-seq(1,nrow(ancestries))[ancestries[,candidates[i]]==1]
    total_candidates<-vector("list",length(offspring))
    total_combinations<-vector("list",length(offspring))
    
    for (j in 1:length(offspring)){
      if (offspring[j]==1){
        #total<-many[j]
        #extra<-choose(many[j],2) this is not being used
        posval1<-rep(0,candidates[i])
        posval2<-rep(0,candidates[i])
        potential<-c(potential,candidates[i])
        posval1[potential]<-1
        sing<-tajima$singletons[tajima$singletons<=candidates[i]]
        
        posval2[sing]<-1
        total_candidates[[j]]<-seq(1,candidates[i])[posval1==1 & posval2==1]
       # if (length(total_candidates[[j]])<many[j]){
          #find the cherries that need to be used as singletons
        if (parent>2){ total_candidates[[j]]<-sort(c(total_candidates[[j]],total_candidates[[j]][tajima$familysize[total_candidates[[j]]]==2]))}
          
        #}
        if (length(total_candidates[[j]])==1){total_combinations[[j]]<-as.matrix(total_candidates[[j]])}else{
          total_combinations[[j]]<-combn(total_candidates[[j]],many[j])}
          ##correction for the cases the two branches of the cherries are not selected 
          # if (many[j]>1 & ncol(total_combinations[[j]])>1){
          #   cherr<-unique(total_candidates[[j]])[table(total_candidates[[j]])==2]
          # mylittle<-function(x,cc){return(sum(x==cc))}
          # for (cc in cherr){
          #   total_combinations[[j]]<-total_combinations[[j]][,-seq(1,ncol(total_combinations[[j]]))[apply(total_combinations[[j]],2,mylittle,cc=cc)==1]]
          # }
          # }
      }else{
        sizes<-tajima$familysize[potential]
        cand2<-potential[sizes==offspring[j]]
        cand2<-cand2[ancestors[offspring_where[j],][cand2]==1]
        total_candidates[[j]]<-cand2
        if (length(cand2)==1){total_combinations[[j]]<-as.matrix(cand2)}else{
          total_combinations[[j]]<-combn(cand2,many[j])}
        if (offspring[j]==2){##Permute the rows
          permut<-factorial(many[offspring==2])
          permuts<-allPerm(initMC(seq(1,many[j])))
          toreplace<-total_combinations[[j]][permuts[1,],]
          for (p in 2:permut){
            toreplace<-cbind(toreplace,total_combinations[[j]][permuts[p,],])
          }
          total_combinations[[j]]<-toreplace
      }
      }
    }
    ##Here is where we find all possible allocations
    if (length(offspring)>1){
      #if (ncol(total_combinations[[1]])<ncol(total_combinations[[2]])){
      #  allocations_1<-expand_grid_tajima(total_combinations[[2]],total_combinations[[1]],ancestries)
      # }else{
      allocations_1<-as.matrix(expand_grid_tajima(total_combinations[[1]],total_combinations[[2]],ancestries))
      #}
      # allocations_1<-expand_grid_tajima(total_combinations[[1]],total_combinations[[2]],ancestries)
      if (length(offspring)>2){
        for (j in 3:length(offspring)){
          #if (ncol(as.matrix(allocations_1))<ncol(total_combinations[[j]])){
          # allocations_1<-expand_grid_tajima(total_combinations[[j]],as.matrix(allocations_1),ancestries)}
          #else{
          allocations_1<-expand_grid_tajima(allocations_1,total_combinations[[j]],ancestries)
        }
      }
    }else{
      #  allocations_1<-cbind(rep(0,sum(many)),as.matrix(total_combinations[[1]]))
      allocations_1<-as.matrix(total_combinations[[1]])
    }
    #
    #original<-allocation_1
    #permuts<-allPerm(initMC(seq(1,many[offspring==2])))
    #newallocations<-rep(allocations_1,permuts)
    #for (j in 1:nrow(permuts)){
    #  allocations_1<-cbind(allocations_1,)
    #}
    #allocations_1
    #allocations<-allocations_1 
    allocations<-cbind(allocations,as.matrix(allocations_1))
   
    #many_p[i]<-ncol(as.matrix(allocations_1)*permut)  
    many_p[i]<-ncol(as.matrix(allocations_1))  
    }
    
  allocations<-as.matrix(allocations[,-1])
  if (sum(many)==1){allocations<-t(as.matrix(allocations))}
  allocations<-matrix(allocations,nrow=sum(many),ncol=sum(many_p))
  times<-matrix(0,nrow=sum(many),ncol=sum(many_p))
  
  times[flag==0,]<-tajima$d[-1][allocations[flag==0,]]
  times[flag==1,]<-cumsum(u)[allocations[flag==1,]]
  
  #the list has the following hierarchy: candidates, offspring (descending order by size), many times
  return(list(candidates=candidates,many=many_p,times=times,allocations=allocations,flag=flag,manyo=many,newsizes=newsizes,offspring_where2=offspring_where2))
}
#parent<-parent_size
#offspring<-offspring_sizes
#tajima<-bb_use
evidence2<-function(parent,offspring,many,tajima,u,ancestors,ancestries,offspring_where){
  ##Bring the candidates with the same parent size
  candidates<-seq(1,length(tajima$familysize))[tajima$familysize==parent]
  valid_candidate<-rep(1,length(candidates))
  #Are those candidate nodes in the Tajima tree parents of the nodes that correspond to offspring from the data?
  
  
  ##Check constrains on each candidate, each candidate have to have offsprings in ancestors
  if (parent>2){
    for (j in 1:length(candidates)){
      descendents<-seq(1,nrow(ancestries))[ancestries[,candidates[j]]==1]
      condition<-0
      for (i in 1:length(offspring)){
        condition<-condition+ifelse(sum(ancestors[offspring_where[i],descendents])>=1, 1, 0)
      }
      if (condition<length(offspring)){valid_candidate[j]<-0}
    }
    candidates<-candidates[valid_candidate==1]
  }
  
  many_p<-rep(0,length(candidates))
  #many<-rep(1,length(offspring))
  allocations<-as.matrix(rep(0,sum(many)))
  for (i in 1:length(candidates)){
    potential<-seq(1,nrow(ancestries))[ancestries[,candidates[i]]==1]
    total_candidates<-vector("list",length(offspring))
    total_combinations<-vector("list",length(offspring))
    
    for (j in 1:length(offspring)){
      if (offspring[j]==1){
        posval1<-rep(0,candidates[i])
        posval2<-rep(0,candidates[i])
        potential<-c(potential,candidates[i])
        posval1[potential]<-1
        sing<-tajima$singletons[tajima$singletons<=candidates[i]]
        posval2[sing]<-1
        total_candidates[[j]]<-seq(1,candidates[i])[posval1==1 & posval2==1]
          #find the cherries that need to be used as singletons
        if (parent>2){total_candidates[[j]]<-sort(c(total_candidates[[j]],total_candidates[[j]][tajima$familysize[total_candidates[[j]]]==2]))}
        
        if (length(total_candidates[[j]])==1){total_combinations[[j]]<-as.matrix(total_candidates[[j]])}else{
          total_combinations[[j]]<-combn(total_candidates[[j]],many[j])}
          ##correction for the cases the two branches of the cherries are not selected 
          # if (many[j]>1 & ncol(total_combinations[[j]])>1){
          #   cherr<-unique(total_candidates[[j]])[table(total_candidates[[j]])==2]
          # mylittle<-function(x,cc){return(sum(x==cc))}
          # for (cc in cherr){
          #   total_combinations[[j]]<-total_combinations[[j]][,-seq(1,ncol(total_combinations[[j]]))[apply(total_combinations[[j]],2,mylittle,cc=cc)==1]]
          # }}
      }else{
        sizes<-tajima$familysize[potential]
        cand2<-potential[sizes==offspring[j]]
        cand2<-cand2[ancestors[offspring_where[j],][cand2]==1]
        print(cand2)
        total_candidates[[j]]<-cand2
        if (length(cand2)==1){total_combinations[[j]]<-as.matrix(cand2)}else{
          total_combinations[[j]]<-combn(cand2,many[j])}
      }
    }
    ##Here is where we find all possible allocations
    if (length(offspring)>1){
      allocations_1<-expand_grid_tajima(total_combinations[[1]],total_combinations[[2]],ancestries)
      if (length(offspring)>2){
        for (j in 3:length(offspring)){
          allocations_1<-expand_grid_tajima(as.matrix(allocations_1),total_combinations[[j]],ancestries)
        }
      }
    }else{
      #  allocations_1<-cbind(rep(0,sum(many)),as.matrix(total_combinations[[1]]))
      allocations_1<-as.matrix(total_combinations[[1]])
    }
    #allocations<-allocations_1 
    allocations<-cbind(allocations,allocations_1)
    many_p[i]<-ncol(as.matrix(allocations_1))  
    
  }
  allocations<-allocations[,-1]
  if (sum(many)==1){allocations<-t(as.matrix(allocations))}
  allocations<-matrix(allocations,nrow=sum(many),ncol=sum(many_p))
  times<-matrix(0,nrow=sum(many),ncol=sum(many_p))
  flag<-rep(0,sum(many))
  if (sum(offspring==1)==1){
    cc<-sum(many[offspring>1])+1
    flag[cc:(cc+many[offspring==1]-1)]<-1
  }
  times[flag==0,]<-0
  times[flag==1,]<-cumsum(u)[allocations[flag==1,]]
  
  #the list has the following hierarchy: candidates, offspring (descending order by size), many times
  return(list(candidates=candidates,many=many_p,times=times,allocations=allocations,flag=flag,offspring_where2=offspring_where))
}

evidence4<-function(parent,offspring,many,tajima,u,ancestors,ancestries,offspring_where){
  ##Bring the candidates with the same parent size
  candidates<-seq(1,length(tajima$familysize))[tajima$familysize==parent]
  valid_candidate<-rep(1,length(candidates))
  #Are those candidate nodes in the Tajima tree parents of the nodes that correspond to offspring from the data?
  
  
  ##Check constrains on each candidate, each candidate have to have offsprings in ancestors
  if (parent>2){
    for (j in 1:length(candidates)){
      descendents<-seq(1,nrow(ancestries))[ancestries[,candidates[j]]==1]
      condition<-0
      for (i in 1:length(offspring)){
        condition<-condition+ifelse(sum(ancestors[offspring_where[i],descendents])>=1, 1, 0)
      }
      if (condition<length(offspring)){valid_candidate[j]<-0}
    }
    candidates<-candidates[valid_candidate==1]
  }
  ##Cherries needed in candidates?
  gost_cherries<-sum(tajima$familysize==2)
  available_sing<-floor(many[offspring==1]/2)
  if (sum(many[offspring==2])>0){num=many[offspring==2]}else{num=0}
  if (gost_cherries>(num+available_sing) & num>0)
  { #then I can create up to available_sing
    original2<-many[offspring==2]
    many[offspring==2]<-original2+available_sing
    many[offspring==1]<-many[offspring==1]-available_sing*2
  }
  flag<-rep(0,sum(many))
  offspring_where2<-rep(0,sum(many))
  if (sum(offspring==1)==1){
    cc<-sum(many[offspring>1])+1
    flag[cc:(cc+many[offspring==1]-1)]<-1
  }
  cc<-sum(many[offspring>2])+original2+1
  flag[cc:(cc+available_sing)]<-1
  newsizes<-rep(offspring,many)
  
  offspring_where2[newsizes>2]<-offspring_where[offspring>2]
  offspring_where2[newsizes==2 & flag==0]<-offspring_where[offspring==2]
  offspring_where2[newsizes==2 & flag==1]<-offspring_where[offspring==1]
  offspring_where2[newsizes==1 & flag==1]<-offspring_where[offspring==1]
  many_p<-rep(0,length(candidates))
  #many<-rep(1,length(offspring))
  allocations<-as.matrix(rep(0,sum(many)))
  for (i in 1:length(candidates)){
    potential<-seq(1,nrow(ancestries))[ancestries[,candidates[i]]==1]
    total_candidates<-vector("list",length(offspring))
    total_combinations<-vector("list",length(offspring))
    
    for (j in 1:length(offspring)){
      if (offspring[j]==1){
        posval1<-rep(0,candidates[i])
        posval2<-rep(0,candidates[i])
        potential<-c(potential,candidates[i])
        posval1[potential]<-1
        sing<-tajima$singletons[tajima$singletons<=candidates[i]]
        posval2[sing]<-1
        total_candidates[[j]]<-seq(1,candidates[i])[posval1==1 & posval2==1]
      #  if (length(total_candidates[[j]])<many[j]){
          #find the cherries that need to be used as singletons
        if (parent>2){ total_candidates[[j]]<-sort(c(total_candidates[[j]],total_candidates[[j]][tajima$familysize[total_candidates[[j]]]==2]))}
          
     #   }
        if (length(total_candidates[[j]])==1){total_combinations[[j]]<-as.matrix(total_candidates[[j]])}else{
          total_combinations[[j]]<-combn(total_candidates[[j]],many[j])}
          ##correction for the cases the two branches of the cherries are not selected 
          # if (many[j]>1 & ncol(total_combinations[[j]])>1){
          #   cherr<-unique(total_candidates[[j]])[table(total_candidates[[j]])==2]
          # mylittle<-function(x,cc){return(sum(x==cc))}
          # for (cc in cherr){
          #   total_combinations[[j]]<-total_combinations[[j]][,-seq(1,ncol(total_combinations[[j]]))[apply(total_combinations[[j]],2,mylittle,cc=cc)==1]]
          # }}
      }else{
        sizes<-tajima$familysize[potential]
        cand2<-potential[sizes==offspring[j]]
        cand2<-cand2[ancestors[offspring_where[j],][cand2]==1]
        total_candidates[[j]]<-cand2
        if (length(cand2)==1){total_combinations[[j]]<-as.matrix(cand2)}else{
          total_combinations[[j]]<-combn(cand2,many[j])}
       
         if (offspring[j]==2){##Permute the rows
          permut<-factorial(many[offspring==2])
          permuts<-allPerm(initMC(seq(1,many[j])))
          toreplace<-total_combinations[[j]][permuts[1,],]
          for (p in 2:permut){
            toreplace<-cbind(toreplace,total_combinations[[j]][permuts[p,],])
          }
          total_combinations[[j]]<-toreplace
        }
      }
    }
    ##Here is where we find all possible allocations
    if (length(offspring)>1){
      #if (ncol(total_combinations[[1]])<ncol(total_combinations[[2]])){
      #  allocations_1<-expand_grid_tajima(total_combinations[[2]],total_combinations[[1]],ancestries)
      # }else{
      allocations_1<-as.matrix(expand_grid_tajima(total_combinations[[1]],total_combinations[[2]],ancestries))
      if (length(offspring)>2){
        for (j in 3:length(offspring)){
          allocations_1<-expand_grid_tajima(allocations_1,total_combinations[[j]],ancestries)
        }
      }
    }else{
      allocations_1<-as.matrix(total_combinations[[1]])
    }
    #allocations<-allocations_1 
    allocations<-cbind(allocations,allocations_1)
    many_p[i]<-ncol(as.matrix(allocations_1))  
    
  }
  allocations<-allocations[,-1]
  if (sum(many)==1){allocations<-t(as.matrix(allocations))}
  allocations<-matrix(allocations,nrow=sum(many),ncol=sum(many_p))
  times<-matrix(0,nrow=sum(many),ncol=sum(many_p))
  
  times[flag==0,]<-0
  times[flag==1,]<-cumsum(u)[allocations[flag==1,]]
  
  
  #the list has the following hierarchy: candidates, offspring (descending order by size), many times
  return(list(candidates=candidates,many=many_p,times=times,allocations=allocations,flag=flag,manyo=many,newsizes=newsizes,offspring_where2=offspring_where2))
}

combine2<-function(new_alloc,old_alloc,ancestries){
  #they are either independent or connected
  found<-0
  for (j in 1:ncol(old_alloc)){
    if (found==0){
      combinehere<-old_alloc[,j]
      ant<-seq(1,ncol(ancestries))[ancestries[max(old_alloc[,j]),]==1]
      topaste<-rep(0,ncol(new_alloc))
      for (i in 1:ncol(new_alloc)){
        if ( (sum(ant==new_alloc[1,i])>0) || new_alloc[1,i]==old_alloc[nrow(old_alloc),j]){ #an ancestor or the last
          topaste[i]<-1
        }
      }
      if (sum(topaste)>0){
        combinehere<-matrix(rep(old_alloc[,j],sum(topaste)),byrow=FALSE,ncol=sum(topaste),nrow=nrow(old_alloc))
        combinehere<-rbind(combinehere,as.matrix(new_alloc[,topaste==1]))
        found<-1}else{
        found<-0
        }
    }else{
     # where<-ncol(combine)
      ant<-seq(1,ncol(ancestries))[ancestries[max(old_alloc[,j]),]==1]
      topaste<-rep(0,ncol(new_alloc))
      for (i in 1:ncol(new_alloc)){
        if ((sum(ant==new_alloc[1,i])>0) || new_alloc[1,i]==old_alloc[nrow(old_alloc),j]){
          topaste[i]<-1
        }
      }
      if (sum(topaste)>0){
      combine2<-matrix(rep(old_alloc[,j],sum(topaste)),byrow=FALSE,ncol=sum(topaste),nrow=nrow(old_alloc))
      combine2<-rbind(combine2,as.matrix(new_alloc[,topaste==1]))
      combinehere<-cbind(combinehere,combine2)
      }
    }
    
  }
  return(combinehere)
}
combine<-function(new_alloc,old_alloc){
  #rear<-sort(old_alloc[,1],index.return=T)$ix
  #old_alloc<-old_alloc[rear,]
  #rear2<-sort(new_alloc[,1],index.return=T)$ix
  #new_alloc<-new_alloc[rear2,]
  nr<-nrow(old_alloc)
  nc2<-ncol(new_alloc)
  nc<-ncol(old_alloc)
  match<-0
  for (co in 1:nc2){
    match<-match+sum(new_alloc[1,co]==old_alloc[nr,])
  }
  if (match==nc2){ #I keep the columns in the new table
    res<-matrix(0,nrow=nr+nrow(new_alloc)-1,ncol=nc2)
    for (co in 1:nc2){
      where<-seq(1,nc)[new_alloc[1,co]==old_alloc[nr,]]
      res[1:nr,co]<-old_alloc[,where]
      res[(nr+1):(nr+nrow(new_alloc)-1),co]<-new_alloc[-1,co]
    }
  }
  if (match>nc2){
    res<-matrix(0,nrow=nr+nrow(new_alloc)-1,ncol=nc)
    for (co in 1:nc2){
      where<-seq(1,nc)[new_alloc[1,co]==old_alloc[nr,]]
      res[1:nr,where]<-old_alloc[,where]
      res[(nr+1):(nr+nrow(new_alloc)-1),where]<-new_alloc[-1,co]
    }
  }
  if (match==0 & nc2==1){
    res<-matrix(0,nrow=nr+nrow(new_alloc),ncol=nc)
    res[1:nr,]<-old_alloc
    res[(nr+1):(nr+nrow(new_alloc)),seq(1,nc)]<-new_alloc
  }
  if (match==0 & nc==1){
    res<-matrix(0,nrow=nr+nrow(new_alloc),ncol=nc2)
    res[1:nr,]<-old_alloc
    res[(nr+1):(nr+nrow(new_alloc)),seq(1,nc2)]<-new_alloc
  }
  return(res)
}

qF<-function(suff,F,bb_use){
  #suff is the old suff
  #Calculation of qF for all possible allocations 
  #It still needs to be validated
  u<-rep(1,nrow(F))
  suff_use<-suff
  #bb_use<-tajima_info2(F)
  ancestries<-f_ancestries(bb_use)
  nnodes<-nrow(suff$nodes) 
  parent_node<-max(suff_use$nodes[,2])
  parent_list<-parent_node
  node_list<-NULL
  che<-0
  alloc<-0
  cherries<-sum(bb_use$familysize==2)
  needed_cherries<-sum(suff_use$nodes[,3]==2)
  if (cherries>needed_cherries){che<-1}
  
  flag_list<-1
  parent_where<-seq(1,nnodes)[suff_use$nodes[,1]==parent_node]
  m<-matrix(1,nrow=nnodes,ncol=length(u)) #the ncols is just to be consistent but I can't place mutations on n-1 node
  ancestors<-matrix(1,nrow=nnodes,ncol=length(u)) #It records the Tajima information that corresponds to the data
  while (parent_node>0){
    parent_size<-suff_use$nodes[parent_where,3] #relative to original list
    offspring_nodes<-suff_use$nodes[suff_use$nodes[,2]==parent_node,1]
    offspring_where<-rep(0,length(offspring_nodes))
    offspring_sizes<-offspring_where
    many<-offspring_where
    for (j in 1:length(offspring_nodes)){
      offspring_where[j]<-seq(1,nnodes)[suff_use$nodes[,1]==offspring_nodes[j]]
      many[j]<-suff_use$mylist[[offspring_where[j]]]$x
    }
    offspring_sizes<-suff_use$nodes[suff_use$nodes[,2]==parent_node,3]
    order_off<-sort(offspring_sizes,decreasing=TRUE,index.return=T)$ix
    offspring_where<-offspring_where[order_off]
    many<-many[order_off]
    offspring_sizes<-offspring_sizes[order_off]
    
    ##here the correction
    
    if (che==1 & many[offspring_sizes==1]>2){
      E<-evidence4(parent_size,offspring_sizes,many,bb_use,u,ancestors,ancestries,offspring_where)
      # many<-E$mayo
    }else{
      E<-evidence2(parent_size,offspring_sizes,many,bb_use,u,ancestors,ancestries,offspring_where)
    }
    
    if (parent_node==max(suff_use$nodes[,2])){alloc<-E$allocations}else{
      #  E<-evidence2(parent_size,offspring_sizes,many,bb_use,u,ancestors,ancestries,offspring_where)
      
      if (nrow(E$allocations)==1){alloc<-E$allocations
      }else{
        alloc<-combine2(E$allocations,alloc,ancestries)}
    }
    U<-update(E,parent_where,ancestors)
    ancestors<-U$ancestors
    #It's ok to repeat, because they are either nested or independent
    #If I remove the independent, I loose info, I need to remove repeats
    #in cal_q
    parent_list<-c(parent_list,rep(parent_node,nrow(E$allocations)))
    if (is.null(node_list)){node_list=suff_use$nodes[E$offspring_where2,1]}else{
      node_list<-c(node_list,suff_use$nodes[E$offspring_where2,1])}
    flag_list<-c(flag_list,E$flag)
    parent_node<-max(suff_use$nodes[suff_use$nodes[,2]<parent_node,2])
    parent_where<-seq(1,nnodes)[suff_use$nodes[,1]==parent_node]
    
  }
  
  parent_size<-length(u)+1 #relative to original list
  offspring_nodes<-suff_use$nodes[suff_use$nodes[,2]==parent_node,1]
  offspring_where<-rep(0,length(offspring_nodes))
  offspring_sizes<-offspring_where
  many<-offspring_where
  for (j in 1:length(offspring_nodes)){
    offspring_where[j]<-seq(1,nnodes)[suff_use$nodes[,1]==offspring_nodes[j]]
    many[j]<-suff_use$mylist[[offspring_where[j]]]$x
  }
  offspring_sizes<-suff_use$nodes[suff_use$nodes[,2]==parent_node,3]
  order_off<-sort(offspring_sizes,decreasing=TRUE,index.return=T)$ix
  offspring_where<-offspring_where[order_off]
  many<-many[order_off]
  offspring_sizes<-offspring_sizes[order_off]
  
  E<-evidence2(parent_size,offspring_sizes,many,bb_use,u,ancestors,ancestries,offspring_where)
  parent_list<-c(parent_list,rep(parent_node,nrow(E$allocations)))
  flag_list<-c(flag_list,E$flag)
  if (is.null(node_list)){node_list=suff_use$nodes[E$offspring_where2,1]}else{
    node_list<-c(node_list,suff_use$nodes[E$offspring_where2,1])}
  
  if (nrow(E$allocations)==1 || alloc==0){alloc<-E$allocations
  }else{
    alloc<-combine2(E$allocations,alloc,ancestries)}
  ##Note: the flag_list here indicates if the mutation happened at 
  ##the external branch or internal and not necessarily indicates true leaf
  #I need sufficient to use the leaf indicator, this will be done in cal_q
  # tofix<-unique(parent_list[flag_list==0])
  # for (s in tofix){
  #   if (sum(suff_use$nodes[,2]==s & suff_use$nodes[,3]>=2)==1){
  #  #   parent_list[parent_list==s & flag_list==0]<-suff_use$nodes[suff_use$nodes[,2]==s & suff_use$nodes[,3]==2,1]
  #     parent_list[parent_list==s & flag_list==0]<-suff_use$nodes[suff_use$nodes[,2]==s & suff_use$nodes[,3]>1,1]
  #   }
  # }
  parent_list[parent_list==0]<-max(suff_use$nodes[suff_use$nodes[,2]==0,1])
  
  return(list(alloc=alloc,parent_list=parent_list[-1],flag_list=flag_list[-1],ancestries=ancestries,node_list=node_list))  
} 


#path<-checkpath(path,E$candidates,ancestries)

checkpath<-function(path,candidates,ancestries){
  indicator<-rep(0,length(path))
  for (j in 1:length(path)){
    find<-0
    if (sum(candidates==path[j])==0){
      #then the path is not exactly the same node, but it can be descendent
      can2<-seq(1,nrow(ancestries))[ancestries[,path[j]]==1]
      while(find==0){
        if (sum(can2[length(can2)]==candidates)>0){
          find<-1
          path[j]<-max(can2)
        }
        else{
          if (length(can2)>2){
            can2<-can2[-length(can2)]}else{indicator[j]<-1; find<-1}
        }
      }
    }
      
    }
    return(list(indicator=indicator,path=path))
  }
  
  

TajimaLik<-function(u,suff,F){
  #Likelihood from a single F matrix
 
  suff_use<-suff
  bb<-tajima_info(u,F) #this is the only place I use u directly
  che<-0
  path=NULL
  cherries<-sum(bb$familysize==2) ##available
  needed_cherries<-sum(suff_use$nodes[,3]==2)
  if (cherries>needed_cherries){che<-1}
  ancestries<-f_ancestries(bb)
  suff_use<-suff
  nnodes<-nrow(suff$nodes) 
  parent_node<-max(suff_use$nodes[,2])
  parent_where<-seq(1,nnodes)[suff_use$nodes[,1]==parent_node]
  m<-matrix(1,nrow=nnodes,ncol=length(u)) #the ncols is just to be consistent but I can't place mutations on n-1 node
  ancestors<-matrix(1,nrow=nnodes,ncol=length(u)) #It records the Tajima information that corresponds to the data
  while (parent_node>0){
    parent_size<-suff_use$nodes[parent_where,3] #relative to original list
    offspring_nodes<-suff_use$nodes[suff_use$nodes[,2]==parent_node,1]
    offspring_where<-rep(0,length(offspring_nodes))
    offspring_sizes<-offspring_where
    many<-offspring_where
    for (j in 1:length(offspring_nodes)){
      offspring_where[j]<-seq(1,nnodes)[suff_use$nodes[,1]==offspring_nodes[j]]
      many[j]<-suff_use$mylist[[offspring_where[j]]]$x
    }
    offspring_sizes<-suff_use$nodes[suff_use$nodes[,2]==parent_node,3]
    order_off<-sort(offspring_sizes,decreasing=TRUE,index.return=TRUE)$ix
    offspring_where<-offspring_where[order_off]
    many<-many[order_off]
    offspring_sizes<-offspring_sizes[order_off]
    if (che==1 & many[offspring_sizes==1]>2){
      E<-evidence3(parent_size,offspring_sizes,many,bb,u,ancestors,ancestries,offspring_where)
      if (is.null(path)){
        path=E$candidates
        gd=list(matrix(0,nrow=nnodes,ncol=length(u)))
        if (length(path)>1){
          for (j in 2:length(path)){
          gd[[j]]<-matrix(0,nrow=nnodes,ncol=length(u))
        }}}else{
          tocheck<-checkpath(path,E$candidates,ancestries)
          toeliminate<-tocheck$indicator
          toeliminate<-rep(1,length(path))
          path<-tocheck$path[indicator==0]
          # for (j in 1:length(E$candidates)){
          #   if (sum(path==E$allocations[1,j])>0){
          #   where<-seq(1,length(path))[path==E$allocations[1,j]]
          #   path[where]<-E$candidates[j]
          #   toeliminate[where]<-0
          # }}
          if (sum(toeliminate)>0){
            whoto<-sort(seq(1,length(path))[toeliminate==1],decreasing=TRUE)   
            for (e in whoto){ 
              gd[[e]]<-NULL
            }
           # path<-path[-seq(1,length(path))[toeliminate==1]]
          }
        }
      # many<-E$mayo
      m_use<-CalculateConditional2(bb,parent_where,offspring_where,E,suff_use,many,m,gd,path,0)
      }else{
      E<-evidence(parent_size,offspring_sizes,many,bb,u,ancestors,ancestries,offspring_where)
      if (is.null(path)){
        path=E$candidates
        gd=list(matrix(0,nrow=nnodes,ncol=length(u)))
        if (length(path)>1){
          for (j in 2:length(path)){
          gd[[j]]<-matrix(0,nrow=nnodes,ncol=length(u))
        }}}else{
          toeliminate<-rep(1,length(path))
          for (j in 1:length(E$candidates)){
            where<-seq(1,length(path))[path==E$allocations[1,j]]
            path[where]<-E$candidates[j]
            toeliminate[where]<-0
          }
          if (sum(toeliminate)>0){
           whoto<-sort(seq(1,length(path))[toeliminate==1],decreasing=TRUE)   
           for (e in whoto){ 
           gd[[e]]<-NULL
           }
            path<-path[-seq(1,length(path))[toeliminate==1]]
          }
        }
      
        m_use<-CalculateConditional(bb,parent_where,offspring_where,E,suff_use,many,m,gd,path,0)
      }
    m<-m_use$m
    gd<-m_use$gd
    U<-update(E,parent_where,ancestors)
    ancestors<-U$ancestors
    parent_node<-max(suff_use$nodes[suff_use$nodes[,2]<parent_node,2])
    #parent_node3<-suff_use$nodes[suff_use$nodes[,2]==parent_node2 & suff_use$nodes[,3]==2,1]
    #if (parent_node2<parent_node3 & parent_node3<parent_node)
    #I should treat the cherries a bit different
    #if (suff_use$nodes[suff_use$nodes,2]==parent_node )
    parent_where<-seq(1,nnodes)[suff_use$nodes[,1]==parent_node]
  }
  
  parent_size<-length(u)+1 #relative to original list
  offspring_nodes<-suff_use$nodes[suff_use$nodes[,2]==parent_node,1]
  offspring_where<-rep(0,length(offspring_nodes))
  offspring_sizes<-offspring_where
  many<-offspring_where
  for (j in 1:length(offspring_nodes)){
    offspring_where[j]<-seq(1,nnodes)[suff_use$nodes[,1]==offspring_nodes[j]]
    many[j]<-suff_use$mylist[[offspring_where[j]]]$x
  }
  offspring_sizes<-suff_use$nodes[suff_use$nodes[,2]==parent_node,3]
  order_off<-sort(offspring_sizes,decreasing=TRUE,index.return=T)$ix
  offspring_where<-offspring_where[order_off]
  many<-many[order_off]
  offspring_sizes<-offspring_sizes[order_off]
  
  E<-evidence(parent_size,offspring_sizes,many,bb,u,ancestors,ancestries,offspring_where)
  # if (is.null(path)){
  #   path=E$candidates
  #   gd=list(matrix(0,nrow=nnodes,ncol=length(u)))
  #   if (length(path)>1){
  #     for (j in 2:length(path)){
  #       gd[[j]]<-matrix(0,nrow=nnodes,ncol=length(u))
  #     }}}
  if (is.null(path)){
    path=E$candidates
    gd=list(matrix(0,nrow=nnodes,ncol=length(u)))
    if (length(path)>1){
      for (j in 2:length(path)){
        gd[[j]]<-matrix(0,nrow=nnodes,ncol=length(u))
      }}}else{
        toeliminate<-rep(1,length(path))
        for (j in 1:length(E$candidates)){
          where<-seq(1,length(path))[path==E$allocations[1,j]]
          path[where]<-E$candidates[j]
          toeliminate[where]<-0
        }
        if (sum(toeliminate)>0){
          whoto<-sort(seq(1,length(path))[toeliminate==1],decreasing=TRUE)   
          for (e in whoto){ 
            gd[[e]]<-NULL
          }
          path<-path[-seq(1,length(path))[toeliminate==1]]
        }
      }
  optional<-CalculateConditional(bb,1,offspring_where,E,suff_use,many,m,gd,path,1)
  L<-sum(u*seq(length(u)+1,2))
  gd<-optional$gd
  res_gradient<-rep(0,length(u))
  likelihood<-log(sum(optional$optional))-L
  weight<-optional$optional/sum(optional$optional)
  if (length(gd)>1){
    res_gradient<-matrix(0,nrow=length(gd),ncol=length(u))
    for (j in 1:length(gd)){
      res_gradient<-res_gradient+weight[j]*apply(gd[[j]],2,sum)
    }
  }else{
    res_gradient<-apply(gd[[1]],2,sum)-seq(length(u)+1,2)
  }
  return(list(likelihood=likelihood,gradient=res_gradient))  
}  
  
listcond<-function(entry,value,list,wei){
  list2<-list()
  wei2<-c(0)
  nextv<-c(0)
  for (j in 1:length(list)){
    if (sum(abs(list[[j]][entry,(1:(entry-1))]-value))==0){
      list2[[length(list2)+1]]<-list[[j]]
      wei2<-c(wei2,wei[j])
      nextv<-c(nextv,list[[j]][(entry+1),1])
    }
  }
  return(list(wei2=wei2[-1],list2=list2,nextv=nextv[-1]))
}

listcond2<-function(entry,value,list,wei){
  list2<-list()
  wei2<-c(0)
  nextv<-rep(0,entry-2)
  for (j in 1:length(list)){
    if (list[[j]][entry,1]==value){
      list2[[length(list2)+1]]<-list[[j]]
      wei2<-c(wei2,wei[j])
      nextv<-rbind(nextv,list[[j]][entry,(2:(entry-1))])
    }
  }
  return(list(wei2=wei2[-1],list2=list2,nextv=nextv[-1,]))
}
listcond3<-function(entry,value,list,wei){
  list2<-list()
  wei2<-c(0)
  for (j in 1:length(list)){
    if (sum(abs(list[[j]][entry,2:(entry-1)]-value))==0){
      list2[[length(list2)+1]]<-list[[j]]
      wei2<-c(wei2,wei[j])
    }
  }
  return(list(wei2=wei2[-1],list2=list2))
}

FmarkovianMean<-function(n,test1){
  meanF<-matrix(0,nrow=n-1,ncol=n-1)
  diag(meanF)<-seq(n,2)
  meanF[2,1]<-n-2
  
  list1<-listcond(2,n-2,test1$Flist,exp(test1$Flist2)*test1$Flist4/sum(exp(test1$Flist2)*test1$Flist4))
  value<-round(sum(list1$nextv*list1$wei2)/sum(list1$wei2))
  meanF[3,1]<-value
  meanF[3,2]<-meanF[2,2]-2
  value<-meanF[3,1:2]
 
  for (i in 3:(n-2)){
    
    list1<-listcond(i,value,list1$list2,list1$wei2)
  #  list1<-listcond(i,value,test1$Flist,exp(test1$Flist2)*test1$Flist4/sum(exp(test1$Flist2)*test1$Flist4))
    
    table(list1$nextv)
    value<-round(sum(list1$nextv*list1$wei2)/sum(list1$wei2))
    meanF[(i+1),1]<-value
    if ((meanF[(i),1]-value[1])==2){
      meanF[(i+1),(2:i)]<-meanF[i,(2:i)]-2
      value<-meanF[(i+1),(1:i)]
    }
    if ((meanF[(i),1]-value[1])<2){
      list2<-listcond2((i+1),value,list1$list2,list1$wei2)
      if (length(list2$list2)>1){
        values<-apply(list2$nextv,1,paste,collapse="-")
        uniqval<-unique(values)
        freq<-rep(0,length(uniqval))
        for ( j in 1:length(uniqval)){
          freq[j]<-sum(list2$wei2[values==uniqval[j]])
        }
        who<-round(sum(seq(1,length(uniqval))*freq)/sum(freq))
        meanF[(i+1),(2:i)]<-as.numeric(strsplit(uniqval[who],"-")[[1]])
        list1<-listcond3(i+1,meanF[(i+1),(2:i)],list2$list2,list2$wei2)
      }else{
        meanF[(i+1),(2:i)]<-list2$nextv
      }
      value<-meanF[(i+1),(1:i)]
    }
    
  }
  return(meanF)
}

FposteriorMean<-function(n,test1){
  meanF<-matrix(0,nrow=n-1,ncol=n-1)
  diag(meanF)<-seq(n,2)
  meanF[2,1]<-n-2
  
  list1<-listcond(2,n-2,test1$Flist,test1$freqF/sum(test1$freqF))
  table(list1$nextv)
  value<-round(sum(list1$nextv*list1$wei2)/sum(list1$wei2))
  meanF[3,1]<-value
  meanF[3,2]<-meanF[2,2]-2
  value<-meanF[3,1:2]
  
  for (i in 3:(n-2)){
    
    list1<-listcond(i,value,list1$list2,list1$wei2)
    #  list1<-listcond(i,value,test1$Flist,exp(test1$Flist2)*test1$Flist4/sum(exp(test1$Flist2)*test1$Flist4))
    
    table(list1$nextv)
    value<-round(sum(list1$nextv*list1$wei2)/sum(list1$wei2))
    meanF[(i+1),1]<-value
    if ((meanF[(i),1]-value[1])==2){
      meanF[(i+1),(2:i)]<-meanF[i,(2:i)]-2
      value<-meanF[(i+1),(1:i)]
    }
    if ((meanF[(i),1]-value[1])<2){
      list2<-listcond2((i+1),value,list1$list2,list1$wei2)
      if (length(list2$list2)>1){
        values<-apply(list2$nextv,1,paste,collapse="-")
        uniqval<-unique(values)
        freq<-rep(0,length(uniqval))
        for ( j in 1:length(uniqval)){
          freq[j]<-sum(list2$wei2[values==uniqval[j]])
        }
        who<-round(sum(seq(1,length(uniqval))*freq)/sum(freq))
        meanF[(i+1),(2:i)]<-as.numeric(strsplit(uniqval[who],"-")[[1]])
        list1<-listcond3(i+1,meanF[(i+1),(2:i)],list2$list2,list2$wei2)
      }else{
        meanF[(i+1),(2:i)]<-list2$nextv
      }
      value<-meanF[(i+1),(1:i)]
    }
    
  }
  return(meanF)
}

FposteriorMode<-function(n,test1){
  meanF<-matrix(0,nrow=n-1,ncol=n-1)
  diag(meanF)<-seq(n,2)
  meanF[2,1]<-n-2
  
  list1<-listcond(2,n-2,test1$Flist,test1$freqF/sum(test1$freqF))
  mg <- aggregate(list1$wei2/sum(list1$wei2), by=list(list1$nextv), FUN=sum)
  value<-mg[mg[,2]==max(mg[,2]),1]
  meanF[3,1]<-value
  meanF[3,2]<-meanF[2,2]-2
  value<-meanF[3,1:2]
  
  for (i in 3:(n-2)){
    
    list1<-listcond(i,value,list1$list2,list1$wei2)
    #  list1<-listcond(i,value,test1$Flist,exp(test1$Flist2)*test1$Flist4/sum(exp(test1$Flist2)*test1$Flist4))
    mg <- aggregate(list1$wei2/sum(list1$wei2), by=list(list1$nextv), FUN=sum)
    #value<-round(sum(list1$nextv*list1$wei2)/sum(list1$wei2))
    value<-mg[mg[,2]==max(mg[,2]),1]
    
    #table(list1$nextv)
    #value<-round(sum(list1$nextv*list1$wei2)/sum(list1$wei2))
    meanF[(i+1),1]<-value
    
    if ((meanF[(i),1]-value[1])==2){
      meanF[(i+1),(2:i)]<-meanF[i,(2:i)]-2
      value<-meanF[(i+1),(1:i)]
    }
    if ((meanF[(i),1]-value[1])<2){
      list2<-listcond2((i+1),value,list1$list2,list1$wei2)
      if (length(list2$list2)>1){
        values<-apply(list2$nextv,1,paste,collapse="-")
        uniqval<-unique(values)
        freq<-rep(0,length(uniqval))
        for ( j in 1:length(uniqval)){
          freq[j]<-sum(list2$wei2[values==uniqval[j]])
        }
        # mg <- aggregate(freq/sum(freq), by=list(uniqval), FUN=sum)
        who<-seq(1,length(uniqval))[freq==max(freq)]
        meanF[(i+1),(2:i)]<-as.numeric(strsplit(uniqval[who],"-")[[1]])
        list1<-listcond3(i+1,meanF[(i+1),(2:i)],list2$list2,list2$wei2)
      }else{
        meanF[(i+1),(2:i)]<-list2$nextv
      }
      value=meanF[(i+1),(1:i)]
    }
  }
  return(meanF)
}



FmarkovianMode<-function(n,test1){
  meanF<-matrix(0,nrow=n-1,ncol=n-1)
  diag(meanF)<-seq(n,2)
  meanF[2,1]<-n-2
  
  list1<-listcond(2,n-2,test1$Flist,exp(test1$Flist2)*test1$Flist4)
  mg <- aggregate(list1$wei2/sum(list1$wei2), by=list(list1$nextv), FUN=sum)
  value<-mg[mg[,2]==max(mg[,2]),1]
  meanF[3,1]<-value
  meanF[3,2]<-meanF[2,2]-2
  
  
  for (i in 3:(n-2)){
    
    list1<-listcond(i,value,list1$list2,list1$wei2)
    table(list1$nextv)
    mg <- aggregate(list1$wei2/sum(list1$wei2), by=list(list1$nextv), FUN=sum)
    #value<-round(sum(list1$nextv*list1$wei2)/sum(list1$wei2))
    value<-mg[mg[,2]==max(mg[,2]),1]
    meanF[(i+1),1]<-value
    if ((meanF[(i),1]-value)==2){
      meanF[(i+1),(2:i)]<-meanF[i,(2:i)]-2
    }
    if ((meanF[(i),1]-value)<2){
      list2<-listcond2((i+1),value,list1$list2,list1$wei2)
      if (length(list2$list2)>1){
        values<-apply(list2$nextv,1,paste,collapse="-")
        uniqval<-unique(values)
        freq<-rep(0,length(uniqval))
        for ( j in 1:length(uniqval)){
          freq[j]<-sum(list2$wei2[values==uniqval[j]])
        }
       # mg <- aggregate(freq/sum(freq), by=list(uniqval), FUN=sum)
        who<-seq(1,length(uniqval))[freq==max(freq)]
        meanF[(i+1),(2:i)]<-as.numeric(strsplit(uniqval[who],"-")[[1]])
        list1<-listcond3(i+1,meanF[(i+1),(2:i)],list2$list2,list2$wei2)
      }else{
        meanF[(i+1),(2:i)]<-list2$nextv
      }
      
    }
    
  }
  return(meanF)
}

mat_shift_diff <- function(F1, F2) {
  dimf<-nrow(F1)
  diff1<- F1[2:(dimf),2:(dimf)] - F1[2:(dimf),1:(dimf-1)] 
  diff1[upper.tri(diff1)] <- 0
  diff2<-F2[2:(dimf),2:(dimf)]-F2[2:dimf,1:(dimf-1)]
  diff2[upper.tri(diff2)] <- 0
  mat_diff(diff1, diff2)
}

mat_diff <- function(F1, F2) {
  sum(abs(F1 - F2))
}



make_F <- function(tree) {
  
  # the number of leaves is number of internal nodes + 1
  n <- tree$Nnode + 1
  
  # initialize F matrix with correct diagonal entries
  F <- diag(seq(1,n))
  F[1,1] <- 0
  
  # initialize label for new nodes we'll create
  ind <- n + 1
  
  # fill in the entries one at a time; i is the row
  for (i in n:2) {
    
    # keep record of tips visited so far in this row (which have coalesced)
    alltips <- character()
    
    # j is the column
    for (j in (i-1):2) {
      
      # find all the tips that have been visited (have coalesced)
      alltips <- unique(c(alltips, tips(tree, i + j)))
      
      # the entry in F is the number that haven't coalesced
      F[i,j] = i - length(alltips)
    }
    
    # drop the two most recent tips and rename the resulting branch for next row
    if (i > 2) {
      tree <- drop.tip(tree, tips(tree, 2*i-1), trim.internal = FALSE)
      tree$tip.label[length(tree$tip.label)] <- paste("t",ind)
      ind <- ind + 1
    }
  }
  
  # for some reason this entry seems to get zeroed out
  F[2,2] <- 2
  
  # return F
  F
}


generateFprior<-function(n){
  #Only for n>=3 ##TODO: send a message of error when n<=1 and fix it to work for n=2 (just one entry=2) for n=3 (a 2x2 matrix)
  ##According to Genetics paper, Oct 2015
  ##Fixed the function, Feb 2019
  F<-matrix(0,nrow=n-1,ncol=n-1) #this one doesn't have the column of 0s
  diag(F)<-seq(n,2)
  F[2,1]<-n-2
  current.block.size<-c(2)
  vintage.size<-rep(1,n-1)
  vintage<-n
  c<-1
  for (j in 3:(n-1)){
    ran<-sum(rmultinom(1,1,p=c(choose(F[j-1,1],2),F[j-1,1]*(n+2-j-F[j-1,1]),choose(n+2-j-F[j-1,1],2)))*c(2,1,0))
    ran
    F[j,1]<-F[j-1,1]-ran
    if (ran==2){c<-c+1;vintage<-c(vintage,n+2-j);vintage.size[j-1]<-length(vintage);F[j,2:(j-1)]<-F[j-1,2:(j-1)]-2; current.block.size<-c(current.block.size,2)}
    if (ran==1){
      # F[j,1]<-F[j-1,1]-1
      who<-sample(1:vintage.size[j-2],1)
      vintage.who<-vintage[who]
      F[j,1:(n+1-vintage.who)]<-F[j-1,1:(n+1-vintage.who)]-1
      if ((n+2-vintage.who)<=(j-1)){
        F[j,(n+2-vintage.who):(j-1)]<-F[j-1,(n+2-vintage.who):(j-1)]-2}
      vintage<-c(vintage[-who],n+2-j)
      vintage.size[j-1]<-vintage.size[j-2]
    }
    if (ran==0){
      who<-sample(1:vintage.size[j-2],2)
      vintagewho<-vintage[who]
      whomin<-min(vintagewho)
      whomax<-max(vintagewho)
      F[j,1:(n+1-whomax)]<-F[j-1,1:(n+1-whomax)]
      if ((n+2-whomax)<=(n+1-whomin)){
        F[j,(n+2-whomax):(n+1-whomin)]<-F[j-1,(n+2-whomax):(n+1-whomin)]-1}
      if ((n+2-whomin)<=(j-1)){
        F[j,(n+2-whomin):(j-1)]<-F[j-1,(n+2-whomin):(j-1)]-2}
      vintage<-c(vintage[-who],n+2-j)
      vintage.size[j-1]<-vintage.size[j-2]-1
    }
  }
  return(F)
}



tree_from_F <- function(matF, coal_times){
  #generate an ape tree (Jaime's code)
  #F is the actual form used in the code that differs from the paper's notation
  
  n= dim(matF)[1]
  edge=matrix(rep(0,6*n-6),ncol=3)
  edge[,2]= (1:(2*n-2))
  vintages = c()
  times=c(rep(0,n),coal_times)
  for (j in n:2){
    new_node = 2*n-j+1
    next_leaf = intersect(which(edge[,1]==0),1:n)[1]
    F_difference = rev(matF[,j]-matF[,j-1])
    
    if (F_difference[1]==2){
      edge[next_leaf:(next_leaf+1),1]=new_node
      vintages=c(vintages, new_node)
    }
    else if (F_difference[1]==1){
      selected_vintage = which(F_difference == 2)[1]+n-1
      edge[selected_vintage,1]=new_node
      edge[next_leaf,1]=new_node
      vintages = c(vintages[vintages != selected_vintage],new_node)
    }
    else {
      selected_vintage1 =which(F_difference == 1)[1]+n-1
      selected_vintage2 =which(F_difference == 2)[1]+n-1
      edge[selected_vintage1,1]=new_node
      edge[selected_vintage2,1]=new_node
      vintages = vintages[vintages!=selected_vintage1]
      vintages = vintages[vintages!=selected_vintage2]
      vintages<-c(vintages,new_node)
    }
  }
  #edge[5:8,2]=c(6,7,8,5)
  edge[1:n,]=edge[order(edge[1:n,2]),]
  #edge=edge[order(edge[,1]),]
  
  for (j in 1:(2*n-2)) {
    #I don't understand this
    edge[j,3]=times[edge[j,1]]-times[edge[j,2]]
  }
  edge[,1]=3*n-edge[,1]
  edge[-(1:n),2]=3*n-edge[-(1:n),2]
  
  final_tree=rcoal(n,br=coal_times)
  final_tree$edge=edge[,-3]
  final_tree$edge.length=edge[,3]
  final_tree$Nnode=n-1
  class(final_tree) <- "phylo"
  final_tree <- reorder(final_tree,"postorder")
  final_tree$edge[final_tree$edge[, 2] <= n, 2] <- 1:n
  return(final_tree)
}


consolidateF<-function(Fju){
  ##generates an F matrix consistent with paper notation
  newF<-matrix(0,nrow=nrow(Fju),ncol=nrow(Fju))
  for (j in 1:nrow(newF)){
    newF[nrow(Fju)-j+1,]<-rev(Fju[,j])
  }
  newF2<-matrix(0,nrow=Fju[1,1],ncol=Fju[1,1])
  newF2[2:Fju[1,1],2:Fju[1,1]]<-newF
  return(newF2)
}

optimTimes<-function(times,nsites,iter=200,eps=5e-5,threshold=.1){
  score<-rep(0,iter)
  val1<-python.call("calcPF", result$change_F, result$F_nodes, result$family_size, oldsuff,times*nsites,"True")
  #times2<-times+eps*val1[[3]]*nsites
  #times2[times2<0]<-0.000001
  for (j in 1:iter){
    times2<-times+eps*val1[[3]]
    times2[times2<0]<-0.00000001
    val2<-python.call("calcPF", result$change_F, result$F_nodes, result$family_size, oldsuff,times2*nsites,"True")
    if ((val2[[2]]-val1[[2]])>threshold){
      times<-times2
      val1<-val2
      score[j]<-val1[[2]]
      print(score[j])
      if (score[j]>0){print(times)}
    }else{
      j<-iter
      print(paste("stop at",j,sep=""))
    }
  }
  return(list(times=times,score=score))
}
getDistributionF<-function(n,Nsim,Finfo=NULL){
  #update February 2018
  #Function to validate marginal distribution of F
  if (is.null(Finfo)){
    samp<-generateFprior(n)
  }else{
    samp<-Finfo
  }
  Flist<-list(samp)
  freqF<-1
  for (j in 1:(Nsim-1)){
    samp<-generateFprior(n)
    exists<-0
    for (i in 1:length(Flist)){
      if ((sum(abs(Flist[[i]]-samp)>0)!=0)) { #The F matrixx is the same
        exists<-exists+1
      }else{
        freqF[i]<-freqF[i]+1;
      }
      if (exists==length(Flist)){
        Flist[[length(Flist)+1]]<-samp
        freqF<-c(freqF,1)
      }
    }
  }
  return(list(Flist=Flist,freqF=freqF))
  #return(list(Flist=Flist,Flist2=Flist2,Flist3=Flist3,freqF=freqF))
}

plotall<-function(somelist,depthmax=NULL,main="DensiTree"){
  z<-length(somelist$Flist)
  treelist<-list()
  for (j in 1:z){
    treelist[[j]]<-ladderize(tree_from_F(consolidateF(somelist$Flist[[j]]),coal_times=cumsum(somelist$tList[[j]]/somelist$freqF[j])))
    #treelist[[j]]$tip.label<-sort(treelist[[j]]$tip.label)
    print(treelist[[j]]$tip.label)
  }
  
  class(treelist)<-"multiPhylo"
  x<-treelist 
  depth<-rep(0,length(x))
  for (j in 1:length(x)){
    depth[j]<-coalescent.intervals(x[[j]])$total.depth
  }
  if (is.null(depthmax)){depthmax<-max(depth)}
  print(depthmax)
  plot.new()
  plot.window(xlim = c(0, depthmax), ylim = c(0, n+1),main=main)
  
  for (j in 1:length(x)){
    width<-somelist$freqF[j]/max(somelist$freqF)
   # width<-.3
    treeindex<-j
    tmp <- reorder(x[[treeindex]])
    tmp <- ladderize(tmp)
    xy <- plotPhyloCoor(tmp)
    xx <- xy[, 1]+depthmax-max(xy[,1])
    #xx<-xx/max(xx)
    yy <- xy[, 2]
    #yy <- yy/max(yy)
    horizontal=TRUE
    e1 <- tmp$edge[, 1]
    Ntip <- min(e1) - 1L
    Nnode <- tmp$Nnode
    
    phylogram.plot(tmp$edge, Ntip, Nnode, xx, yy, horizontal, 
                   edge.color = adjustcolor("gray", alpha.f =width), 
                   edge.width = 1, edge.lty = 1)
    
  }
  #axisPhylo()
}

plotall_med<-function(somelist,depthmax=NULL){
  z<-length(somelist$Flist)
  treelist<-list()
  for (j in 1:z){
   if (is.vector(somelist$tlist2[[j]])){
     coal_time<-somelist$tlist2[[j]]
   }else{
    coal_time<-apply(somelist$tlist2[[j]],1,quantile,p=.5)
    }
    treelist[[j]]<-ladderize(tree_from_F(consolidateF(somelist$Flist[[j]]),coal_times=cumsum(coal_time)))
    #treelist[[j]]$tip.label<-sort(treelist[[j]]$tip.label)
    print(treelist[[j]]$tip.label)
  }
  
  class(treelist)<-"multiPhylo"
  x<-treelist 
  depth<-rep(0,length(x))
  for (j in 1:length(x)){
    depth[j]<-coalescent.intervals(x[[j]])$total.depth
  }
  if (is.null(depthmax)){depthmax<-max(depth)}
  print(depthmax)
  plot.new()
  plot.window(xlim = c(0, depthmax), ylim = c(0, n+1))
  
  for (j in 1:length(x)){
    width<-somelist$freqF[j]/max(somelist$freqF)
    # width<-.3
    treeindex<-j
    tmp <- reorder(x[[treeindex]])
    tmp <- ladderize(tmp)
    xy <- plotPhyloCoor(tmp)
    xx <- xy[, 1]+depthmax-max(xy[,1])
    #xx<-xx/max(xx)
    yy <- xy[, 2]
    #yy <- yy/max(yy)
    horizontal=TRUE
    e1 <- tmp$edge[, 1]
    Ntip <- min(e1) - 1L
    Nnode <- tmp$Nnode
    
    phylogram.plot(tmp$edge, Ntip, Nnode, xx, yy, horizontal, 
                   edge.color = adjustcolor("gray", alpha.f =width), 
                   edge.width = 1, edge.lty = 1)
    
  }
  #axisPhylo()
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



##This function takes a newick tree and generates the Tajima information needed for the python code  
createFpython<-function(tree){
  FMat<-create_F(tree)
  bb<-bring_branch_lengths(rep(1,nrow(FMat)),FMat)
  singletons<-bb$singletons[-length(bb$singletons)]
  result<-list(change_F=c(-1*diff(FMat[,1]),0),F_nodes=list(),family_size=c(bb$familysize[-1],FMat[1,1]))
    for (j in singletons){
      result$F_nodes[[j]]<-c(-1,j)
    }
    for ( j in unique(bb$nodes[,2])){
      #print(j)
      if ((sum(bb$nodes[,2]==j))==1){
        result$F_nodes[[j]]<-c(bb$nodes[bb$nodes[,2]==j,1],j)
      }else{
        result$F_nodes[[j]]<-list()
        result$F_nodes[[j]][[1]]<-c(min(bb$nodes[bb$nodes[,2]==j,1]),j)
        result$F_nodes[[j]][[2]]<-c(max(bb$nodes[bb$nodes[,2]==j,1]),j)
      }
    }
  changes<-seq(1,nrow(FMat))[FMat[nrow(FMat),]==1]
  if (min(changes)>1){
    result$F_nodes[[(FMat[1,1]-1)]]<-list()
    result$F_nodes[[(FMat[1,1]-1)]][[1]]<-c(min(changes)-1,FMat[1,1]-1)
    result$F_nodes[[(FMat[1,1]-1)]][[2]]<-c(max(changes),FMat[1,1]-1)
  }else{
    result$F_nodes[[(FMat[1,1]-1)]]<-c(max(changes),FMat[1,1]-1)
    result$change_F[(FMat[1,1]-1)]<-1
  }
    return(result)
  }
  