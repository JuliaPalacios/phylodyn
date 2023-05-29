# coding: utf-8
"""
    Function to be called from R. Adopted from notebook file
    
    Requirement:
        Install rPython
        
    Usage:
        library(rPython)
        ...create oldsuff...
        python.load("PATH/TO/Fsamp.py")
        result <- python.call("F_sample", oldsuff$nodes)    # result contains F_mat, total_vintage and change_F
"""
# # use this object to represent each node

# In[ ]:
import math
import copy
import collections
# from IPython.core.debugger import Tracer
class TreeNode:
    """
        self.lineage: list of int, keep track of constructed lineage under this node, 
        e.g. [0,0,1,1,2]
        0 means child node not constructed, 1 means available individual (for leaf node), 
        2 means child node constructed
        Initially it is [0]*(#descendents) for non-leaf nodes and [1]*(#descendents) for leaf node.
            - 12/17 change: compress to [#free, #constructed]   
            - 12/17 change2: to construct F mat, use dictionary {0:[list of free], 1:[list of constructed]}
                list of free is all -1, as initialization
                When combined, two -1s are removed from list_free. 
                One new element is added to constructed, with value as current index of operation (0,1,2...)
        
        self.parent: pointer to parent node
        
        self.childs: list of unconstructed child nodes. 
        When a child is constructed, it removes itself from parent.childs. 
        When childs is empty, add this parent to sample pool.
            - 12/17 change: record number of children only
            
        self.mutations: for computing probability
    """
    def __init__(self, name, isPF=False):
        self.name = int(name)
        self.isPF = isPF
        if isPF:
            self.lineage = {0:[],1:set()}
        else:
            self.lineage = {0:[],1:[]} 
        self.parent = None
        self.numChilds = 0  # use num of childeren left
        self.mutation = 0  # number of mutations
        self.mutFromSingleton = {}
        self.leaves = {-1:1}    # map of number of leaves (free individuals) from each lineage
        self.family_size = 0
        self.original_node = 1.0
    def __repr__(self):
        return "TreeNode " +' '+str(self.name)+': '+str(self.lineage)+' mutation: '+str(self.mutation)+' mut from sing: '+str(self.mutFromSingleton)+ 'original'+str(self.original_node)
        #return 'TreeNode '+str(self.name)+': '+str(self.lineage)
    def copy(self):
        """ return deep copy of current TreeNode"""
        newnode = TreeNode(self.name, self.isPF)
        newnode.lineage[0] = self.lineage[0][:]
        newnode.lineage[1] = copy.deepcopy(self.lineage[1])
        newnode.parent = self.parent
        newnode.numChilds = self.numChilds
        newnode.mutation = self.mutation
        newnode.mutFromSingleton = self.mutFromSingleton
        newnode.leaves = self.leaves # not used, no need to deep copy
        newnode.original_node=self.original_node#.......copy original node
        return newnode
    def restoreLineage(self, node):
        """ restore lineage of self from node"""
        self.lineage[0] = node.lineage[0][:]
        self.lineage[1] = copy.deepcopy(node.lineage[1])
        self.numChilds = node.numChilds


# # Initialization

# In[ ]:

import numpy as np
import random


# # initialize tree nodes

# In[ ]:

# load suf_ext matrix
# suf_ext = np.loadtxt(open("data.csv", "r"), delimiter=",",dtype=np.int)

def binomialCoeff(n, k):
    result = 1
    for i in range(1, k+1):
        result = result * (n-i+1) / i
    return result


# In[ ]:
def F_sample_het(oldsuff):
#def F_sample(oldsuff,seed):
    # current sample pool as list (set) of currently available nodes
        
    #np.random.seed(seed)
    oldsuff['nodes']=np.array(oldsuff['nodes'])
    poolN, nodesDict, Noldsuff,original = createTree(oldsuff)
    L=np.max(oldsuff['nodes'][:,0:2])
    pool, nodes = createTree_t0(oldsuff,L)
    #print(pool)
    #print(nodesDict)
    #print(Noldsuff)
    # In[ ]:
    # FIX: calculate actual individuals
    nno=np.array(oldsuff['nodes'])
    n = sum(nno[nno[:,1]==0,2]*nno[nno[:,1]==0,3])
    #print(n)
    #n=5
    # initialize var to store result
    total_vintage = []
    probs = []
    family_size = []
    F_nodes = []
    indicator = [0]*(n-1)
    #print(pool)
    sampleF_het(oldsuff,pool,nodes, total_vintage, probs, family_size, F_nodes, indicator,Noldsuff)

    if len(family_size)==n-1:
        family_group=family_groupGEN(Noldsuff,F_nodes)
        node_group=node_groupGEN(Noldsuff,family_group,total_vintage,indicator).tolist()
    else:
        node_group=[]
    #print(node_group[0])
    # return F_mat,total_vintage,change_F 
    # FIX:
    # return dict, return part of F_mat
    return {"total_vintage":total_vintage, "probs":probs, "family_size":family_size, "F_nodes":F_nodes, "indicator":indicator,"node_group":node_group}

# # Helper function

# In[ ]:
def createTree(data, isPF = False):

    # create tree bottom up
    data = copy.deepcopy(data)
    mutations = data['mylist']
    oldsuff = data['nodes']
    oldsuff = np.array(oldsuff,copy=True)
    suff_sorted = sorted([(olds,mut) for olds,mut in zip(oldsuff,mutations)], key=lambda row: row[0][1], reverse=True)
    #print("length mutations?"+str(len(mutations)))
    #original = np.zeros(oldsuff[:,0])+1.0
    original = np.zeros(len(mutations))+1.0
    #print("original creation"+str(original))
    # add another column at the end 
    #suff_sorted = np.concatenate(suff_sorted, suff_sorted[:,0],axis=1)
    # clean up the matrix (column 3 and column 4)
    newoldsuff = np.array([0]*oldsuff.shape[1])
    extraLines = []
    L = max(r[0][0] for r in suff_sorted)
    for row,mut in suff_sorted:
        #print("original row"+str(row[0]))
        if row[2]>1 and row[3]>1:
            #print("this is where we extend")
            for i in range(1,row[3]):
                newline = np.copy(row)
                newline[3] = 1
                newline[0] = L+1
                newoldsuff=np.vstack((newoldsuff,newline))
                newline = (newline,{'x':1,'y':mut['y'][i]})
                L += 1
                extraLines.append(newline)
                #print(len(original))
                if i==1: 
                    #print("value of x"+str(mut['x'])+"and factorial"+str(math.factorial(mut['x']))+"and permutations "+str(len(set(list(itertools.permutations(mut['y']))))))
                    #original=np.append(original,row[0])
                    original=np.append(original,math.factorial(mut['x'])/len(set(list(itertools.permutations(mut['y'])))))
                    #len(list(itertools.permutations(mut['y']))))
                else:   
                    original=np.append(original,1.0)
                #print("original extended"+str(original))
                # add 6th column
            row[3] = 1
            mut['x'] = 1
            mut['y'] = mut['y'][0]
    Noldsuff=np.vstack((newoldsuff,oldsuff))
    Noldsuff=np.delete(Noldsuff,(0),axis=0)
    suff_sorted += extraLines
    suff_sorted = sorted(suff_sorted, key=lambda row: row[0][1], reverse=True)
    pool = []
    nodes = {}
    i=0
    #print ("all good?")
    for row,mut in suff_sorted:
        #print("row is"+str(row[0]))
        #print("original after"+str(original))
        # assert mut['x'] == row[3], "wrong mutation['x']"
        i+=1
        #print("i is"+str(i))
        #newnode.original_node = original[i-1.0]
        if row[4] == 1 and row[2]==1: #leaf and no merged (?) L: and group of singletons I'd say
            # if is leaf, add parent node (but not strictly parent)
            # FIX: 
            # clamped individuals to parent when leaf and column 3 = 1
            if row[1] not in nodes:
                #print("this?"+str(row[0]))
                newnode = TreeNode(row[1], isPF)
                nodes[row[1]] = newnode 
                newnode.lineage[0] += [-row[6]]*row[3]
                newnode.family_size += row[3]
                #assert newnode.mutFromSingleton == [], 'mutFromSingleton initialized'
                newnode.mutFromSingleton[-row[6]]= [mut['y']] if isinstance(mut['y'],int) else mut['y']
                newnode.leaves[-row[6]]=1
                newnode.original_node=1.0
            else:
                #print("or this?")
                newnode = nodes[row[1]]
                newnode.lineage[0] += [-row[6]]*row[3]  
                newnode.family_size += row[3]
                #assert newnode.mutFromSingleton == [], 'mutFromSingleton initialized'
                newnode.mutFromSingleton[-row[6]] = [mut['y']] if isinstance(mut['y'],int) else mut['y']
                newnode.original_node=1.0
                newnode.leaves[-row[6]]=1
            if len(newnode.lineage[0])+len(newnode.lineage[1]) > 1  and newnode not in pool:
                # only add to pool if have more than one free
                pool.append(newnode)
        else:
            # leaf node with column 3 >1 or non-leaf
            if row[0] not in nodes:
                newnode = TreeNode(row[0], isPF)
                #print("this row is not in nodes"+str(row[0]))
                #newnode.original_node=original[i-1]
                #print("new node being created"+str(newnode.name)+"and original"+str(newnode.original_node))
            else:
                #print("check if this is called"+str(len(original))+" and row?"+str(row[0]-1))
                newnode = nodes[row[0]]
                #newnode.original_node=original[i-1]
            assert newnode.mutation == 0, 'mutation initialzed'
            newnode.leaves[-row[6]]=1 #TO CHECK!!!!!
            newnode.mutation = mut['y']
            newnode.original_node=original[i-1]
            if row[1] not in nodes:
                parent = TreeNode(row[1], isPF)
                parent.leaves[-row[6]]=1 #TO CHECK!!!!!
                nodes[row[1]] = parent
                
            else:
                parent = nodes[row[1]]
            parent.numChilds += 1
            newnode.parent = parent
            if row[4]==1:
                # leaf node with column 3 >1, has more than one free
                newnode.lineage[0] += [-row[6]]*row[2]
                newnode.family_size += row[2]
                pool.append(newnode)
            parent.family_size += newnode.family_size
            nodes[row[0]] = newnode
            #print(nodes[row[0]].original_node)
            #nodes[row[0]].original_node=original[row[0]]
        #print (nodes[i-1].name)
        #print(nodes[i-1].original_node)
    nodes[-1] = TreeNode(-1,isPF)
    nodes[0].parent = nodes[-1]
    nodes[-1].family_size = nodes[0].family_size
    
    #print("checking")
    return pool, nodes, Noldsuff,original




def createTree_t0(data,L, isPF = False):

    # create tree bottom up
    data = copy.deepcopy(data)
    mutations = data['mylist']
    oldsuff = data['nodes']
    oldsuff = np.array(oldsuff)
    idx = np.where(oldsuff[:,5]==0)[0] #1 is the first time event
    mutations=[mutations[i] for i in idx]
    oldsuff=oldsuff[idx]
    oldsuff = np.array(oldsuff,copy=True)
    suff_sorted = sorted([(olds,mut) for olds,mut in zip(oldsuff,mutations)], key=lambda row: row[0][1], reverse=True)
    #print("length mutations?"+str(len(mutations)))
    #original = np.zeros(oldsuff[:,0])+1.0
    original = np.zeros(len(mutations))+1.0
    #print("original creation"+str(original))
    # add another column at the end 
    #suff_sorted = np.concatenate(suff_sorted, suff_sorted[:,0],axis=1)
    # clean up the matrix (column 3 and column 4)
    
    #newoldsuff = np.array([0]*oldsuff.shape[1])
    extraLines = []
    #L = max(r[0][0] for r in suff_sorted)
    for row,mut in suff_sorted:
        #print("original row"+str(row[0]))
        if row[2]>1 and row[3]>1:
            #print("this is where we extend")
            for i in range(1,row[3]):
                newline = np.copy(row)
                newline[3] = 1
                newline[0] = L+1
                #newoldsuff=np.vstack((newoldsuff,newline))
                newline = (newline,{'x':1,'y':mut['y'][i]})
                L += 1
                extraLines.append(newline)
                #print(len(original))
                if i==1: 
                    #print("value of x"+str(mut['x'])+"and factorial"+str(math.factorial(mut['x']))+"and permutations "+str(len(set(list(itertools.permutations(mut['y']))))))
                    #original=np.append(original,row[0])
                    original=np.append(original,math.factorial(mut['x'])/len(set(list(itertools.permutations(mut['y'])))))
                    #len(list(itertools.permutations(mut['y']))))
                else:   
                    original=np.append(original,1.0)
                #print("original extended"+str(original))
                # add 6th column
            row[3] = 1
            mut['x'] = 1
            mut['y'] = mut['y'][0]
    #Noldsuff=np.vstack((newoldsuff,oldsuff))
    #Noldsuff=np.delete(Noldsuff,(0),axis=0)
    suff_sorted += extraLines
    suff_sorted = sorted(suff_sorted, key=lambda row: row[0][1], reverse=True)
    pool = []
    nodes = {}
    i=0
    #print ("all good?")
    
    for row,mut in suff_sorted:
        #print("row is"+str(row[0]))
        #print("original after"+str(original))
        # assert mut['x'] == row[3], "wrong mutation['x']"
        i+=1
        #print("i is"+str(i))
        #newnode.original_node = original[i-1.0]
        if row[4] == 1 and row[2]==1: #leaf and no merged (?) L: and group of singletons I'd say
            # if is leaf, add parent node (but not strictly parent)
            # FIX: 
            # clamped individuals to parent when leaf and column 3 = 1
            if row[1] not in nodes:
                newnode = TreeNode(row[1], isPF)
                nodes[row[1]] = newnode 
                newnode.lineage[0] += [-row[6]]*row[3]
                newnode.family_size += row[3]
                #assert newnode.mutFromSingleton == [], 'mutFromSingleton initialized'
                newnode.mutFromSingleton[-row[6]]= [mut['y']] if isinstance(mut['y'],int) else mut['y']
                newnode.leaves[-row[6]]=1
                newnode.original_node=1.0
            else:
                #print("or this?")
                newnode = nodes[row[1]]
                newnode.lineage[0] += [-row[6]]*row[3]  
                newnode.family_size += row[3]
                #assert newnode.mutFromSingleton == [], 'mutFromSingleton initialized'
                newnode.mutFromSingleton[-row[6]] = [mut['y']] if isinstance(mut['y'],int) else mut['y']
                newnode.original_node=1.0
                newnode.leaves[-row[6]]=1
            if len(newnode.lineage[0])+len(newnode.lineage[1]) > 1  and newnode not in pool:
                # only add to pool if have more than one free
                pool.append(newnode)
        else:
            # leaf node with column 3 >1 or non-leaf
            if row[0] not in nodes:
                newnode = TreeNode(row[0], isPF)
                #print("this row is not in nodes"+str(row[0]))
                #newnode.original_node=original[i-1]
                #print("new node being created"+str(newnode.name)+"and original"+str(newnode.original_node))
            else:
                #print("check if this is called"+str(len(original))+" and row?"+str(row[0]-1))
                newnode = nodes[row[0]]
                #newnode.original_node=original[i-1]
            assert newnode.mutation == 0, 'mutation initialzed'
            newnode.leaves[-row[6]]=1 #TO CHECK!!!!!
            newnode.mutation = mut['y']
            newnode.original_node=original[i-1]
            if row[1] not in nodes:
                parent = TreeNode(row[1], isPF)
                parent.leaves[-row[6]]=1 #TO CHECK!!!!!
                nodes[row[1]] = parent
                
            else:
                parent = nodes[row[1]]
            parent.numChilds += 1
            newnode.parent = parent
            if row[4]==1:
                # leaf node with column 3 >1, has more than one free
                newnode.lineage[0] += [-row[6]]*row[2]
                newnode.family_size += row[2]
                pool.append(newnode)
            parent.family_size += newnode.family_size
            nodes[row[0]] = newnode
            #print(nodes[row[0]].original_node)
            #nodes[row[0]].original_node=original[row[0]]
        #print (nodes[i-1].name)
        #print(nodes[i-1].original_node)
    nodes[-1] = TreeNode(-1,isPF)
    nodes[0].parent = nodes[-1]
    nodes[-1].family_size = nodes[0].family_size
    
    #print("checking")
    return pool, nodes


def updateTree(data,pool,nodes, cnt,L, isPF = False):

    # create tree bottom up
    data = copy.deepcopy(data)
    mutations = data['mylist']
    oldsuff = data['nodes']
    oldsuff=np.array(oldsuff,copy=True)
    idx = np.where(oldsuff[:,5]==(cnt+1))[0] #add the elements necessary at the next iteration.
    if len(idx)>0: #it does it only if there is sthg to update
        mutations=[mutations[i] for i in idx]
        oldsuff=oldsuff[idx]
        oldsuff = np.array(oldsuff,copy=True)
        suff_sorted = sorted([(olds,mut) for olds,mut in zip(oldsuff,mutations)], key=lambda row: row[0][1], reverse=True)
        #print("length mutations?"+str(len(mutations)))
        #original = np.zeros(oldsuff[:,0])+1.0
        original = np.zeros(len(mutations))+1.0
        #print("original creation"+str(original))
        # add another column at the end 
        #suff_sorted = np.concatenate(suff_sorted, suff_sorted[:,0],axis=1)
        # clean up the matrix (column 3 and column 4)
        extraLines = []
        L = max(r[0][0] for r in suff_sorted)
        for row,mut in suff_sorted:
            #print("original row"+str(row[0]))
            if row[2]>1 and row[3]>1:
                #print("this is where we extend")
                for i in range(1,row[3]):
                    newline = np.copy(row)
                    newline[3] = 1
                    newline[0] = L+1
                    newline = (newline,{'x':1,'y':mut['y'][i]})
                    L += 1
                    extraLines.append(newline)
                    #print(len(original))
                    if i==1: 
                        #print("value of x"+str(mut['x'])+"and factorial"+str(math.factorial(mut['x']))+"and permutations "+str(len(set(list(itertools.permutations(mut['y']))))))
                        #original=np.append(original,row[0])
                        original=np.append(original,math.factorial(mut['x'])/len(set(list(itertools.permutations(mut['y'])))))
                        #len(list(itertools.permutations(mut['y']))))
                    else:   
                        original=np.append(original,1.0)
                    #print("original extended"+str(original))
                    # add 6th column
                row[3] = 1
                mut['x'] = 1
                mut['y'] = mut['y'][0]
        suff_sorted += extraLines
        suff_sorted = sorted(suff_sorted, key=lambda row: row[0][1], reverse=True)
        i=0
        #print ("all good?")
        for row,mut in suff_sorted:
            #print("row is"+str(row[0]))
            #print("original after"+str(original))
            # assert mut['x'] == row[3], "wrong mutation['x']"
            i+=1
            #print("i is"+str(i))
            #newnode.original_node = original[i-1.0]
            if row[4] == 1 and row[2]==1: #leaf node and singletons
                # clamped individuals to parent when leaf and column 3 = 1
                if row[1] not in nodes:
                    #print("this?"+str(row[0]))
                    newnode = TreeNode(row[1], isPF)
                    nodes[row[1]] = newnode 
                    newnode.lineage[0] += [-row[6]]*row[3]
                    newnode.family_size += row[3]
                    #assert newnode.mutFromSingleton == [], 'mutFromSingleton initialized'
                    newnode.mutFromSingleton[-row[6]]= [mut['y']] if isinstance(mut['y'],int) else mut['y']
                    newnode.original_node=1.0
                    newnode.leaves[-row[6]]=1
                else:
                    #print("or this?")
                    newnode = nodes[row[1]]
                    newnode.lineage[0] += [-row[6]]*row[3]  
                    newnode.family_size += row[3]
                    #assert newnode.mutFromSingleton == [], 'mutFromSingleton initialized'
                    newnode.mutFromSingleton[-row[6]] = [mut['y']] if isinstance(mut['y'],int) else mut['y']
                    newnode.original_node=1.0
                    newnode.leaves[-row[6]]=1
                if len(newnode.lineage[0])+len(newnode.lineage[1]) > 1  and newnode not in pool:
                    # only add to pool if have more than one free
                    pool.append(newnode)
            else:
                # leaf node with column 3 >1 or non-leaf
                if row[0] not in nodes:
                    newnode = TreeNode(row[0], isPF)
                    #print("this row is not in nodes"+str(row[0]))
                    #newnode.original_node=original[i-1]
                    #print("new node being created"+str(newnode.name)+"and original"+str(newnode.original_node))
                else:
                    #print("check if this is called"+str(len(original))+" and row?"+str(row[0]-1))
                    newnode = nodes[row[0]]
                    #newnode.original_node=original[i-1]
                assert newnode.mutation == 0, 'mutation initialzed'
                newnode.mutation = mut['y']
                newnode.original_node=original[i-1]
                newnode.leaves[-row[6]]=1 #TO CHECK!!!!
                if row[1] not in nodes:
                    parent = TreeNode(row[1], isPF)
                    parent.leaves[-row[6]]=1 #TO CHECK!!!!!
                    nodes[row[1]] = parent
                else:
                    parent = nodes[row[1]]
                parent.numChilds += 1
                newnode.parent = parent
                if row[4]==1:
                    # leaf node with column 3 >1, has more than one free
                    newnode.lineage[0] += [-row[6]]*row[2]
                    newnode.family_size += row[2]
                    pool.append(newnode)
                parent.family_size += newnode.family_size
                nodes[row[0]] = newnode
                #print(nodes[row[0]].original_node)
                #nodes[row[0]].original_node=original[row[0]]
            #print (nodes[i-1].name)
            #print(nodes[i-1].original_node)
        nodes[-1] = TreeNode(-1,isPF)
        nodes[0].parent = nodes[-1]
        nodes[-1].family_size = nodes[0].family_size
    
    #print("checking")
    return pool, nodes
    
def calcFamilyNodes(nodes):
    """
        compute up-wards family_size list for each node
    """
    res = collections.defaultdict(set,set())
    def helper(node):
        if node is None:
            return set()
        if node.name not in res:
            res[node.name].add(node.family_size)
            res[node.name] = res[node.name].union(helper(node.parent))
        return res[node.name]
    for n in nodes:
        helper(nodes[n])
    return res
    
def calcFamilySample_het(F_nodes, family_size):
    """
        compute upward family_size list for each node in sample
    """
    res = collections.defaultdict(set,set())
    for ind,pair in reversed(list(enumerate(F_nodes))):
        fs = family_size[ind]
        res[ind].add(fs)
        res[pair[0]-1] = res[pair[0]-1].union(res[ind])
        res[pair[1]-1] = res[pair[1]-1].union(res[ind])
    #res.pop(-2) What's this for??
    return res

def comb2(total):
    # return number of choices of 2 from total
    return (total-1)*total/2.0

def sampleNodeFromList(nodeList):
    """
        Given a list of node, sample node according to number of non-zero choices in its lineage.
        
        Return sampled node.
        
        Time complexity: with a good data structure, O(1) to random pick an element. 
        Probably O(log(n)) to random pick with sepecified probability
        
        12/17: Use numpy.random.choice
        12/17: what to use for nodelist (aka pool)? - to calc prob, O(n) already
    """
    # p = [comb2(len(node.lineage[0])+len(node.lineage[1])) for node in nodeList]
    # if not any(p):
        # raise Exception("probability error!" +str(nodeList))
    # prob = [1.0*pp/sum(p) for pp in p]
    # choice = np.random.choice(nodeList, p=prob)
    # choice_prob = [pp for pp,node in (prob, nodeList) if node==choice]
    # if len(choice_prob)!=1:
        # raise Exception("choice prob error!")
    # choice = random.choice(nodeList) #FIX: uniform choice
    choice = nodeList[np.random.randint(len(nodeList))]
    return choice, 1.0/len(nodeList)


# In[ ]:

def sampleFromLineage_het(lineage,cnt):
    """
        Sample from lineage of one node, remove those two and combine
        
        Input: lineage, {0:[],1:[]}
        Return: types,ids
            types: two elements, each is among [0,1], denoting [free, constructed children or box].
            ids: index of two elements
    """
    nLineage0 = len(lineage[0])
    nLineage1 = len(lineage[1])
    chooseIds = random.sample(range(len(lineage[0])+len(lineage[1])),2)
    # chooseIds = np.random.choice(range(len(lineage[0])+len(lineage[1])),2,replace=False)
    chooseNums = [lineage[int(cid>=nLineage0)][cid-(cid>=nLineage0)*nLineage0] for cid in chooseIds] #L: some doubts about this first part 
    ids = []  # indexs of combined lineage
    types = []  # 0,1 indicating free or box
    #prob_choices = [comb2(len(lineage[0])), len(lineage[0])*len(lineage[1]), comb2(len(lineage[1]))]    # number of choices of [0,0],[0,1],[1,1]
    values=np.array(list(Counter(lineage[0]).keys()))
    counts_pre = np.array(list(Counter(lineage[0]).values()))
    for cid,cnum in zip(chooseIds,chooseNums):
        ids.append(cnum)
        types.append(int(cid>=nLineage0))
        lineage[int(cid>=nLineage0)].remove(cnum)
    lineage[1].append(cnt+1)  # add a new combined lineage
    counts_aft = [lineage[0].count(x) for x in set(values)]
    prob_num = np.prod([binomialCoeff(counts_pre[i],(counts_pre[i]-counts_aft[i])) for i in range(len(counts_pre))])
    prob_den = comb2(nLineage0+nLineage1)
    return types,ids, prob_num/prob_den, nLineage1


# In[ ]:

def combineAndUpdate_het(node, cnt, total_vintage, probs, family_size, F_nodes):
    """
        Given a node, sample two of lineages and combine. 
        Update total_vintage etc to reflect the combination.
        
        Return nothing. Update in-place.
        
        Time complexity: with a good data structure, O(1) remove, add time.
        To update F_mat and path, O(N)??
    """
    # sample two lineage ang update lineage
    types,ids,prob,origin_nLin1 = sampleFromLineage_het(node.lineage,cnt)
    probs.append(prob)
    ids.sort()
    
    # update F_mat and path
    total_vintage.append(node.name)
    new_family = node.leaves[ids[0]] + node.leaves[ids[1]]
    family_size.append(new_family)
    node.leaves[cnt+1] = new_family
    F_nodes.append((ids[0], ids[1],  cnt+1))

        
    


# # Sample function

# In[ ]:

def sampleF_het(oldsuff,pool,nodes, total_vintage, probs, family_size, F_nodes, indicator,Noldsuff):
    """
        sample available node from pool
        
        return F and path.
        
        In each iteration, one combination happens. Total nunmber of iterations = N.
    """

    cnt = 0
    while len(pool) > 0:
#         Tracer()()
        L=np.max(oldsuff['nodes'][:,0:2])
        node,node_prob = sampleNodeFromList(pool)
        probs.append(node_prob)
        combineAndUpdate_het(node, cnt, total_vintage, probs, family_size, F_nodes)
        updateTree(oldsuff,pool,nodes, cnt,L, isPF = False)
#         print (node.name, total_vintage, change_F)
        #print('ok')
        #print(cnt)
#         print (F_mat)
        if len(node.lineage[0])+len(node.lineage[1])==1:
            # this node has only one lineage left, merge to parent
            # remove node from pool, add parent to pool
            pool.remove(node)  # O(n) with list
            #if node.numChilds==0:
            if node.name != 0:
                if (Noldsuff[Noldsuff[:,0]==node.name,2]* Noldsuff[Noldsuff[:,0]==node.name,3]==family_size[-1]):    
                # if no child left 
                 #it could be possible that a node got free done, but still have children
                        # parent = node.parent
#                        parent.numChilds -= 1
#                        parent.lineage[1] += [cnt+1] # the children returns a box
#                        if len(node.lineage[1]) != 1:
#                            raise Exception('node not finished')
#                        parent.leaves[cnt+1] = node.leaves[node.lineage[1][0]]
#                        indicator[cnt] = 1
#                        # print ("node /{0}/ merge to parent /{1}/".format(node,parent))
#                        if len(parent.lineage[0])+len(parent.lineage[1])==2:
#                            pool.append(parent)
#                elif (oldsuff['nodes'][oldsuff['nodes'][:,0]==node.name,2]*oldsuff['nodes'][oldsuff['nodes'][:,0]==node.name,3])==family_size[-1]:
                    if node.name != -1:
                        parent = node.parent
                        parent.numChilds -= 1
                        parent.lineage[1] += [cnt+1] # the children returns a box
                        if len(node.lineage[1]) != 1:
                            raise Exception('node not finished')
                        parent.leaves[cnt+1] = node.leaves[node.lineage[1][0]]
                        indicator[cnt] = 1
                        # print ("node /{0}/ merge to parent /{1}/".format(node,parent))
                        if len(parent.lineage[0])+len(parent.lineage[1])==2:
                            pool.append(parent)
        cnt += 1

def testSpeed(mat):
    # test speed of R-Py interface
    mat = np.array(mat)
    return (mat+1).tolist()
# In[ ]:

# np.random.seed(2)
# F_mat = np.diag(list(range(n,0,-1)))
# total_vintage = []
# change_F = []
# sampleF(pool)


# In[ ]:

# print (total_vintage, change_F)

# print(F_mat)

def possibleUpdate_het(node, Fnode, ind):
    newprob = None
    oldnode = None
    compatible=True
    values=np.array(list(Counter(node.lineage[0]).keys()))
    counts_pre = np.array(list(Counter(node.lineage[0]).values()))
    nnode=len(node.lineage[0])+len(node.lineage[1])
    oldnode = node.copy()
    if ind+1==Fnode[2] and Fnode[0] in node.lineage[0]:
        node.lineage[0].remove(Fnode[0])
    elif ind+1==Fnode[2] and Fnode[0] in node.lineage[1]:
           # can produce the same Fnode (need to verify one Fnode - index of the box)
        node.lineage[1].remove(Fnode[0])
    else:
        compatible=False
    if ind+1==Fnode[2] and Fnode[1] in node.lineage[0]:
        node.lineage[0].remove(Fnode[1])
    elif ind+1==Fnode[2] and Fnode[1] in node.lineage[1]:
           # can produce the same Fnode (need to verify one Fnode - index of the box)
         node.lineage[1].remove(Fnode[1])
    else:
        compatible=False
    node.lineage[1].add(ind+1)
    counts_aft = [node.lineage[0].count(x) for x in set(values)]
    prob_num = np.prod([binomialCoeff(counts_pre[i],(counts_pre[i]-counts_aft[i])) for i in range(len(counts_pre))])
    prob_den = comb2(nnode)
    newprob=prob_num/prob_den
    if compatible==False:
        node.restoreLineage(oldnode)
    return compatible, oldnode, newprob, node
    
def update(pool, newnode, cnt):
    # update pool and parent for oldnode in pool
#     pool.remove(oldnode)
    isMarginalize = False
    node = newnode
    oldinfo = [None,None, None]
    if len(node.lineage[0])+len(node.lineage[1])==1:
        pool.remove(node)
        if node.numChilds==0 and node.name != 0:
            # node done, need to change its parent
            isMarginalize = True
            oldparent = node.parent.copy()
            oldinfo[1] = oldparent
            parent = node.parent
            parent.numChilds -= 1
            parent.lineage[1].add(cnt+1) # the children returns a box
            if len(node.lineage[1]) != 1:
                raise Exception('node not finished')
            # print ("node /{0}/ merge to parent /{1}/".format(node,parent))
            if len(parent.lineage[0])+len(parent.lineage[1])==2:
                pool.add(parent)
                oldinfo[2] = parent
    else:
        # add new node
        oldinfo[0] = True
    return oldinfo, isMarginalize or node.name==0

maps = {}
from collections import Counter
import itertools


def compare(sp,a):
    # compare if sp is distinct from a
    assert len(sp)==len(a)
    for i in range(len(sp)//2):
        if sorted(sp[i*2:i*2+2]) != sorted(a[i*2:i*2+2]):
            return True
    return False





def singFunc_het_block(muts, bls, inds, L,ids,idc):
    # modify grad in place
    mut=copy.copy(muts)
    perms,index=generatePerm(mut)
    sumpermLik=np.zeros(perms.shape,dtype=np.float64) 
    grad=np.zeros((perms.shape[0],L),dtype=np.float64)
    while len(bls)>0:
        idg=list(bls.keys())[0]
        mut=np.unique(muts[idg])
        for i in range(len(mut)):
            s=0
            for z in index[idg]:
                idm=np.where(perms[:,z]==mut[i])
                sumpermLik[idm,z]+=1.0*np.log(bls[idg][s])*mut[i]
                grad[idm,ids[-idg]:idc[inds[s]]+1]+=1.0*mut[i]/bls[idg][s]
                s+=1
        bls.pop(idg)        
    likBlock=np.exp(sumpermLik.sum(axis=1))
    gradBlock=np.multiply(np.reshape(np.repeat(likBlock,L),grad.shape),grad)
    lik=sum(likBlock)
    grad=gradBlock.sum(axis=0)
    
    return lik, grad
        
        
  
def generatePerm(muts):
    #generate permutations (for the heterochrnous less direct)
    index={}
    idg=list(muts.keys())[0]
    #perms = np.array(list((set(list(itertools.permutations(muts[idg]))))))
    perms=np.array(list(permutations(muts[idg]))) #new permutations
    index[idg]=list(range(0,len(muts[idg])))
    mId=max(index[idg])
    muts.pop(idg)
    while len(muts)>0:
        idg=list(muts.keys())[0]
        #perm=np.array(list((set(list(itertools.permutations(muts[idg]))))))
        perm=np.array(list(permutations(muts[idg]))) 
        perm1=np.tile(perms,(perm.shape[0],1))
        perm=np.repeat(perm,perms.shape[0],axis=0)
        perms=np.hstack((perm1,perm))
        index[idg]=list(range(mId+1,mId+1+len(muts[idg])))
        mId=max(index[idg])
        muts.pop(idg)
    
    return perms, index
        
        


def processMutation_het_singlenode(node,name, this_vintage, F_nodes, timeVec, tMatrix,ids,idc):
    # gradient wrt ts
    grad = np.zeros(tMatrix.shape[1])
    # sum log \sum_v lik_v
    loglik = 0


    # at what point the node got marginalized
    ind_marginalize = max((this_vintage==1)*range(len(this_vintage)))
#    if node.name==0: #It corrects for the special case of the zero node
#        ind_marginalize = np.where(this_vintage==1)[0][0]

    
    #print(node)
    if node.name != 0:
        #print("node name"+str(node.name)+" mut"+str(node.mutation))
        ind_nonZero = tMatrix[ind_marginalize]
        #Modification for the code check
        #loglik += node.mutation*np.log(timeVec[ind_marginalize][2])
        loglik += node.mutation*np.log(timeVec[ind_marginalize][2])
        
        #print("loglik other:"+str(loglik))
        #print ('extending mutation:{}, loglik:{}, len:{}'.format(node.mutation, node.mutation*np.log(timeVec[ind_marginalize][1]), timeVec[ind_marginalize][1]))
        grad[ind_nonZero==True] += 1.0*node.mutation/timeVec[ind_marginalize][2]
        #print("gradient other:"+str(grad))
    #print("mutations:"+str(node.mutation)+"time:"+str(timeVec[ind_marginalize][1])+str(ind_nonZero)+"gradient:"+str(grad))
    
#    if sum([len(node.mutFromSingleton[i]) for i in list(node.mutFromSingleton.keys())])==1: ##it is a dictionary now 
#        # branch from singleton: mutFromSingleton * log(branch length singleton)
#        loglik += (node.mutFromSingleton[list(node.mutFromSingleton.keys())[0]][0])*np.log(timeVec[ind_marginalize][0])-np.log(factorial(node.mutFromSingleton[list(node.mutFromSingleton.keys())[0]][0]))
#        grad[ids[-node.lineage[0][0]]:(ind_marginalize+1)] += 1.0*(node.mutFromSingleton[list(node.mutFromSingleton.keys())[0]][0])/timeVec[ind_marginalize][0]
#        #print("mutations sing:"+str(node.mutFromSingleton[0])+"time:"+str(timeVec[ind_marginalize][0])+"gradient"+str(grad))
#        #print("gradient from singletons:"+str(grad))
#        #print("loglik from singletons:"+str(loglik))

        
    if sum([len(node.mutFromSingleton[i]) for i in list(node.mutFromSingleton.keys())])>=1: ##CHECK IF THIS ONE WORKS
        # likelihood from singletons
        #branches = {1:[],2:[]}  #branches[cf] = [list of branches at change_F=cf]
        #cfs = []
        #inds = {1:[],2:[]}
        branch_vec = {}
        for i in list(node.mutFromSingleton.keys()):
            branch_vec[i]=[]
        inds_vec=[]
        #print("indnonzero"+str(tMatrix[ind_marginalize]))
        #print("ind_marginalize"+str(ind_marginalize))
        for i,tv in enumerate(this_vintage):
            #print("i index is"+str(i)+"tv is"+str(tv))
            if tv ==0 or (F_nodes[i][0]>0 and F_nodes[i][1]>0):
                continue
            if F_nodes[i][0]<=0:
               inds_vec+=[i] ###TO DO
               branch_vec[F_nodes[i][0]]+=[timeVec[i][0]]
            if F_nodes[i][1]<=0:#still in the first column of the vector
               inds_vec+=[i] ## TO DO
               branch_vec[F_nodes[i][1]]+=[timeVec[i][1]]
            #inds_vec +=[i]*change_F[i] ##????
            #print("inds_vec"+str(inds_vec))
            #print("branch_vec"+str(branch_vec))
            #cfs.append(change_F[i])
            #inds[change_F[i]].append(i)
            #branch_vec += [timeVec[i][0]]*change_F[i]
            #print("old inds"+str(inds))
        #print("inds_vec"+str(inds_vec))
        #print("branch_vec"+str(branch_vec))
        this_l,this_g = singFunc_het_block(node.mutFromSingleton,branch_vec,inds_vec,len(grad),ids,idc)
        loglik += np.log(this_l)
        #print("loglik is "+str(loglik))
        grad += 1.0*this_g/this_l
        #print("gradient is "+str(grad))
    #print(name)
    #print(loglik)
    return loglik,grad    # k start from 0?


def pF_het_paths(oldsuff,nodesDict, F_nodes, prob, res, paths,ids,idc, timeVec=[], tMatrix = [], no_output = False):
    """
        iterate through nodes to get probability of F_mat (defined by change_F+F_nodes)
    """
    res[0]+= prob
 
    #empty matrices to fill           
    likpaths=np.array([0]*paths.shape[0], dtype=np.float64)  
    gradpaths=np.tile(np.zeros(tMatrix.shape[1], dtype=np.float64),(paths.shape[0],1))
    for name,node in nodesDict.items():
        if name == -1:
            continue
        locpaths=copy.copy(paths)
        locpaths=(locpaths==name)*1
        while np.sum(locpaths)>0:
            ii=min(np.where(locpaths.sum(axis=1)>0)[0])
            row=locpaths[ii,:]
            diff=abs(locpaths[1:,:]-np.tile(row,(locpaths[1:,:].shape[0],1)))
            idx=np.where(diff.sum(axis=1)==0)[0]+1
            idx=np.append(idx,ii)
            loglik,grad = processMutation_het_singlenode(node,name, row, F_nodes, timeVec, tMatrix,ids,idc) #Julia: changed this part
            likpaths[idx]+=loglik
            gradpaths[idx,:]+=grad
            locpaths[idx,:]=0
            #print(node)
            #print(loglik)
    if not no_output:
        #print("likpaths")
        #print(likpaths)
        res[1]=np.log(sum(np.exp(likpaths)))
        #res[1]=sum(likpaths)
        #print(res[1])
        for i in range(len(likpaths)):
            gradpaths[i,:]=np.exp(likpaths[i])*gradpaths[i,:]
            #print(i)
            #print(gradpaths)
            #print(np.exp(likpaths[i]))
        res[2]=gradpaths.sum(axis=0)/sum(np.exp(likpaths))
       # print(res[2])
    
    return 




def calcTimeVec_het(F_nodes, times,indt):
    # calculate time vectors given F_nodes and change_F
    # vec[i][0] is branch length of singletons, vec[i][1] is extending length
    idc=[i for i, j in enumerate(indt) if j == 0]
    ids=[i for i, j in enumerate(indt) if j > 0]
    cumTimes=[sum(times[:i+1]) for i in range(len(times))]
    cumTimesc=[cumTimes[i] for i in idc]
    ctimes=[cumTimesc[i+1]-cumTimesc[i] for i in range(len(cumTimesc)-1)]
    ctimes.insert(0,sum(times[0:(idc[0]+1)]))
    stimes=[cumTimes[i] for i in ids]
    stimes.insert(0,0)
    vec = [[sum(ctimes[:ind+1]),sum(ctimes[:ind+1]),0] for ind in range(len(ctimes))]
    tMatrix = np.zeros((len(ctimes),len(times)),dtype=bool)
    #print(stimes)
    #print(ctimes)
    for fn in F_nodes:
        if fn[0]==fn[1]:
            vec[fn[2]-1][0]=vec[fn[2]-1][0]-stimes[-fn[0]]
            vec[fn[2]-1][1]=vec[fn[2]-1][1]-stimes[-fn[1]]
        else:
            if fn[0]<=0:
                vec[fn[2]-1][0]=vec[fn[2]-1][0]-stimes[-fn[0]]
            elif fn[0]>0:
                vec[fn[0]-1][2] += cumTimesc[fn[2]-1]-cumTimesc[fn[0]-1]
                tMatrix[fn[0]-1][(idc[fn[0]-1]+1):(idc[fn[2]-1]+1)] = True # these ts are involved for the extending branch
            if fn[1]<=0:
                vec[fn[2]-1][1]=vec[fn[2]-1][1]-stimes[-fn[1]]
            elif fn[1]>0:
                vec[fn[1]-1][2] += cumTimesc[fn[2]-1]-cumTimesc[fn[1]-1]
                tMatrix[fn[1]-1][(idc[fn[1]-1]+1):(idc[fn[2]-1]+1)] = True # these ts are involved for the extending branch
    return np.array(vec),tMatrix

def calcPF_het(F_nodes, family_size, oldsuff,node_group, times=None, indt=None, use_familySize = True, no_output = False):
    """ 
        Usage:
        pf <- python.call("calcPF", result$F_nodes, result$family_size, oldsuff) 
    """
    #print(oldsuff)
    #print(times)
#    print(F_nodes)
#    print(indt)
#    print(family_size)
    

    poolN, nodesDict, Noldsuff,original = createTree(oldsuff, isPF=True)
    L=np.max(np.array(oldsuff['nodes'])[:,0:2])
    pool, nodesDict0 = createTree_t0(oldsuff,L, isPF=True)
    
    n=sum(Noldsuff[Noldsuff[:,1]==0,2]*Noldsuff[Noldsuff[:,1]==0,3])
    leaf=Noldsuff[Noldsuff[:,4]==1,0]
    intern=Noldsuff[Noldsuff[:,4]==0,0]
    intern=np.append(intern,0)
    rel=vintParChil(F_nodes)
    child=childPerPhy(Noldsuff)
    nodeSize=NodeSize(Noldsuff)
    nodeInPath=NodePathFreq(Noldsuff,leaf,intern,nodeSize)
    family_group=family_groupGEN(Noldsuff,F_nodes)
    #node_group=node_groupGEN(Noldsuff,family_group,total_vintage,indicator)
    node_group=np.array(node_group)
    paths=createPath_het(n,child,rel,nodeInPath,leaf,nodeSize,family_size,Noldsuff,family_group,node_group)
    #print(paths)
    prob = 1
    correction = 0.0
    #print(paths)
    assert isinstance(times,list)
    #print("here")
    if times is not None:
        timeVec,tMatrix = calcTimeVec_het(F_nodes,times,indt)
    else:
        timeVec = []
        tMatrix = []
    res = [0,0,np.zeros(timeVec.shape[0])]
    ids=[i for i, j in enumerate(indt) if j > 0]
    ids.insert(0,0)
    idc=[i for i, j in enumerate(indt) if j == 0]
    treeLengthGrad(times,indt,oldsuff)

    if len(paths)>0:
        pF_het_paths(oldsuff,nodesDict, F_nodes, prob, res, paths, ids,idc, timeVec=timeVec, tMatrix = tMatrix, no_output = no_output)
        for ori in original: #CHECK!!! I DON'T KNOW WHAT THIS IS.
            #print("name"+str(name)+"and original"+str(node.original_node))
            correction+=np.log(ori)
        length,grad=treeLengthGrad(times,indt,oldsuff)
        #print(res[1])
        #print(length)
        res[1] -= length
        res[1] -= correction
        #res[2] -= grad
        #res[2] = list(res[2])
    else:
        res[1]= np.log(0)
        
    #print(res[1])
    #twin =twin_datasets_finder(oldsuff,nodesDict,paths[0],F_nodes)
    #print(twin)
    #corr2=2**twin
    res=res[0:2]
    #res[1]+=np.log(corr2)
    
    return res

def treeLengthGrad(times,indt,oldsuff):
    idc=[i for i, j in enumerate(indt) if j == 0]
    ids=[i for i, j in enumerate(indt) if j > 0]
    lbs=list(np.unique(indt))
    cumTimes=[sum(times[:i+1]) for i in range(len(times))]
    cumTimesc=[cumTimes[i] for i in idc]
    ctimes=[cumTimesc[i+1]-cumTimesc[i] for i in range(len(cumTimesc)-1)]
    ctimes.insert(0,sum(times[0:(idc[0]+1)]))
    stimes=[cumTimes[i] for i in ids]
    stimes.insert(0,0)
    isoLength=np.sum(np.arange(len(ctimes)+1, 1, -1) * np.array(ctimes))
    nodes=np.array(oldsuff['nodes'])
    corr=sum([nodes[i,2]*nodes[i,3]*stimes[nodes[i,6]] for i in range(len(nodes)) if nodes[i,4]==1])
    #corr=sum([nodes[i,2]*nodes[i,3]*stimes[nodes[i,6]] for i in range(len(nodes))])
    #print(isoLength)
    #print(corr)
    hetLength=isoLength-corr
    samp_sizes=[sum(nodes[(nodes[:,4]==1) & (nodes[:,6]==i),2]*nodes[(nodes[:,4]==1) & (nodes[:,6]==i),3]) for i in lbs]
    #ids.insert(0,0) #it is needed for the cycle below
    isoGrad=np.array([samp_sizes[0]]*len(times)) 
    for i in range(len(samp_sizes[1:])):
        isoGrad[(ids[i]+1):]=isoGrad[(ids[i]+1):]+samp_sizes[i+1]
    #isoGrad[0]=samp_sizes[0]
    lengthGrad=np.array([0]*len(times))
    lengthGrad[np.array(idc[:-1])+1]=-1*np.array(range(1,sum(samp_sizes)-1))
    #for i in ids[1:]:
    ids=np.array(ids)
    #print(ids) #GRADIENT HAS A BUG
    #for i in ids[ids>0]:
    #    idCS=max(np.where(np.array(idc)<i)[0])
    #    lengthGrad[i+1]=lengthGrad[idc[idCS]+1]
    corrGrad=isoGrad#+lengthGrad
    #print(corrGrad)
    return hetLength,corrGrad





def tuple_based(data):
    new_array = [tuple(row) for row in data]
    return np.unique(new_array)

def lexsort_based(data):                 
        sorted_data =  data[np.lexsort(data.T),:]
        row_mask = np.append([True],np.any(np.diff(sorted_data,axis=0),1))
        return sorted_data[row_mask]

#%% new functions to create the path. 

def vintParChil(F_nodes):
    
    #Define variables to fill
    rel={} #collect parent-child relationship in vintages

    for i in range(len(F_nodes)):
        if F_nodes[i][0]<=0:
            vint1=[]
        else:
            vint1=rel[F_nodes[i][0]]['vint']
            rel[F_nodes[i][0]]['par']=i+1
        if F_nodes[i][1]<=0:
            vint2=[]
        else:
            vint2=rel[F_nodes[i][1]]['vint']
            rel[F_nodes[i][1]]['par']=i+1
            
        vint=np.array(vint1+vint2)
        vint=list(vint[vint.argsort()])
        vint.append(i+1)
        rel[i+1]={'vint':vint,'par':[]}
        
        #Symmetry check
    return rel 

def childPerPhy(Noldsuff):     
    child={}
    parnode=np.unique(Noldsuff[:,1])
    for i in range(len(parnode)):
        node=parnode[i]
        ch=Noldsuff[Noldsuff[:,1]==node,0]
        if len(ch)>2 or sum(Noldsuff[Noldsuff[:,1]==node,3]>1)>=1:
            ch=np.insert(ch,0,node)
        child[node]=ch
    return child

def NodeSize(Noldsuff):

    nodeSize=[0]*(max(Noldsuff[:,0])+1)
    nodeSize[0]=sum(Noldsuff[Noldsuff[:,1]==0,2]*Noldsuff[Noldsuff[:,1]==0,3])
    for i in range(np.shape(Noldsuff)[0]):
        node=Noldsuff[i,0]
        nodeSize[node]=Noldsuff[i,2]
    nodeSize=[0 if i==1 else i for i in nodeSize]
    return np.array(nodeSize)




def NodePathFreq(Noldsuff,leaf,intern,nodeSize):
#Number of times a given node should appear in paths.
       nodeInPath=np.array([0]*(max(Noldsuff[:,0])+1))
       for i in intern:
          nodeInPath[i]=sum((Noldsuff[:,1]==i) & (Noldsuff[:,2]>1))-1
          if len(Noldsuff[(Noldsuff[:,1]==i) & (Noldsuff[:,2]==1),3])>0:
              nodeInPath[i]+= sum(Noldsuff[(Noldsuff[:,1]==i) & (Noldsuff[:,2]==1),3].astype(int))
       nodeInPath[leaf]=nodeSize[leaf]-1
       return nodeInPath

##### Start defininign the paths function

def createPath(n,child,rel,nodeInPath,leaf,nodeSize,family_size,Noldsuff):
    paths=np.array([-1]*(n-1))
    paths=np.reshape(paths,(1,n-1))
    paths[0,n-2]=0
    
    #for cycle here as well
    a=range(n-2)
    b=a[::-1]
    for i in b:
        paths,child=removeChild(Noldsuff,paths,child,nodeInPath,i+1)
        Ipath=paths[paths[:,i]==-1,:]
        if len(Ipath)==0:
            continue
        else:
            rIpath=paths[paths[:,i]!=-1,:]
            parI=rel[i+1]['par']-1 #inded, true parents +1
            parList=list(np.unique(Ipath[:,parI]))
            paths=np.array([-1]*(n-1))
            paths=np.reshape(paths,(1,n-1))
            for node in parList:
                SubIpath=Ipath[Ipath[:,parI]==node,:]
                parNode=child[node]
                PotNode=parNode[nodeSize[parNode]==family_size[i]]
                PotNode=np.append(PotNode,parNode[parNode==node])
                for j in range(np.shape(PotNode)[0]):
                        cIpath=copy.copy(SubIpath)
                        if PotNode[j] in leaf:
                            idSub=np.array(rel[i+1]['vint'])-1
                            for z in idSub:
                               cIpath[:,z]=np.array([PotNode[j]]*np.shape(SubIpath)[0])   
                        else:
                            cIpath[:,i]=np.array([PotNode[j]]*np.shape(SubIpath)[0])
                        paths=np.vstack((paths,cIpath))
                        #remove extra paths not true
                        l=np.equal(paths,PotNode[j]).sum(axis=1)
                        idc=np.where(l<=nodeInPath[PotNode[j]])[0]
                        paths=paths[idc,:]
                        #paths=tuple_based(paths)
                        paths=lexsort_based(paths)
                #in the for cycle above is missing the leaf condition at the moment.
            paths=np.vstack((paths,rIpath))
            paths=paths[1:,:]
    return paths

def createPath_het(n,child,rel,nodeInPath,leaf,nodeSize,family_size,Noldsuff,family_group,node_group):
    paths=np.array([-1]*(n-1))
    paths=np.reshape(paths,(1,n-1))
    paths[0,n-2]=0
    
    #for cycle here as well
    a=range(n-2)
    b=a[::-1]
    for i in b:
        paths,child=removeChild(Noldsuff,paths,child,nodeInPath,i+1)
        Ipath=paths[paths[:,i]==-1,:] # Not all rows are involved, only the one with -1 
        if len(Ipath)==0:
            continue
        else:
            rIpath=paths[paths[:,i]!=-1,:] #Paths that are not inolved
            parI=rel[i+1]['par']-1 #at position i, you are wondering which is the joining node
            parList=list(np.unique(Ipath[:,parI]))
            paths=np.array([-1]*(n-1))
            paths=np.reshape(paths,(1,n-1))
            for node in parList:
                SubIpath=Ipath[Ipath[:,parI]==node,:]
                parNode=child[node]
                #In the heterochrnous version is not enough the family size, I need to check groups compatibility. 
                ifam_gr=np.tile(family_group[:,i],(len(parNode),1)).transpose() #repeat the family group as many times as the par node length. The ,1) is used to create multple rows. 
                idd=np.sum(abs(node_group[:,parNode]-ifam_gr),axis=0) #substract from each family group the size to know if it is compatible
                PotNode=parNode[np.where(idd==0)[0]]
                PotNode=np.append(PotNode,parNode[parNode==node])
                for j in range(np.shape(PotNode)[0]):
                        cIpath=copy.copy(SubIpath)
                        if PotNode[j] in leaf:
                            idSub=np.array(rel[i+1]['vint'])-1
                            for z in idSub:
                               cIpath[:,z]=np.array([PotNode[j]]*np.shape(SubIpath)[0])   
                        else:
                            cIpath[:,i]=np.array([PotNode[j]]*np.shape(SubIpath)[0])
                        paths=np.vstack((paths,cIpath))
                        #remove extra paths not true
                        l=np.equal(paths,PotNode[j]).sum(axis=1)
                        idc=np.where(l<=nodeInPath[PotNode[j]])[0]
                        paths=paths[idc,:]
                        #paths=tuple_based(paths)
                        paths=lexsort_based(paths)
                #in the for cycle above is missing the leaf condition at the moment.
            paths=np.vstack((paths,rIpath))
            paths=paths[1:,:]
    return paths


def removeChild(Noldsuff,paths,child,nodeInPath,i):
    """
    Remove paths that are not completed for node that needs to be removed
    Remove node that needs to be marginalized from child so they are not considered anymore
    """ 
    ToDel=Noldsuff[Noldsuff[:,5]==i,:] # ???
    NodeDel=np.unique(ToDel[:,0])
    for j in range(len(NodeDel)):
        node=NodeDel[j]
        parnode=int(ToDel[ToDel[:,0]==node,1])
        child[parnode]=child[parnode][child[parnode]!=node] #remove from parent children list
        #keep only the path in which node has appeared the right amoount of time.
        if nodeInPath[node]>0:
            l=np.equal(paths,node).sum(axis=1)
            idc=np.where(l==nodeInPath[node])[0]
            paths=paths[idc,:]
    return paths, child


def family_groupGEN(Noldsuff,F_nodes):
    group=np.unique(Noldsuff[:,6])
    family_group=np.zeros([len(group),len(F_nodes)])
    for i in range(len(F_nodes)):
        if F_nodes[i][0]<=0:
            family_group[-F_nodes[i][0],i]+=1
        else:
            family_group[:,i]+=family_group[:,(F_nodes[i][0]-1)]
        if F_nodes[i][1]<=0:
            family_group[-F_nodes[i][1],i]+=1
        else:
            family_group[:,i]+=family_group[:,(F_nodes[i][1]-1)]
    return(family_group)
    
def node_groupGEN(Noldsuff,family_group,total_vintage,indicator):
    """
    
    """
    nodelist=np.unique(Noldsuff[:,0:2])
    node_group=np.zeros([family_group.shape[0],max(nodelist)+1])
    node_group[:,0]=family_group[:,-1]
    for i in np.unique(np.array(total_vintage)[np.array(total_vintage)>0]):
        idn=int(np.where((np.array(total_vintage)==i) & (np.array(indicator)==1))[0])
        node_group[:,i]=family_group[:,idn]
    return(node_group)

#%%  Permutation functions 
    

"""
This module encodes functions to generate the permutations of a multiset
following this algorithm:
Algorithm 1 
Visits the permutations of multiset E. The permutations are stored
in a singly-linked list pointed to by head pointer h. Each node in the linked
list has a value field v and a next field n. The init(E) call creates a
singly-linked list storing the elements of E in non-increasing order with h, i,
and j pointing to its first, second-last, and last nodes, respectively. The
null pointer is given by φ. Note: If E is empty, then init(E) should exit.
Also, if E contains only one element, then init(E) does not need to provide a
value for i.
[h, i, j] ← init(E) 
visit(h) 
while j.n ≠ φ orj.v <h.v do
    if j.n ≠    φ and i.v ≥ j.n.v then 
        s←j
    else
        s←i 
    end if
    t←s.n 
    s.n ← t.n 
    t.n ← h 
    if t.v < h.v then
        i←t 
    end if
    j←i.n 
    h←t 
    visit(h)
end while
... from "Loopless Generation of Multiset Permutations using a Constant Number
of Variables by Prefix Shifts."  Aaron Williams, 2009

"""

class ListElement:
    def __init__(self, value, next):
        self.value = value
        self.next = next
    def nth(self, n):
        o = self
        i = 0
        while i < n and o.next is not None:
            o = o.next
            i += 1
        return o

def init(multiset):
    multiset.sort() # ensures proper non-increasing order
    h = ListElement(multiset[0], None)
    for item in multiset[1:]:
        h = ListElement(item, h)
    return h, h.nth(len(multiset) - 2), h.nth(len(multiset) - 1)

def visit(h):
    """Converts our bespoke linked list to a python list."""
    o = h
    l = []
    while o is not None:
        l.append(o.value)
        o = o.next
    return l


def permutations(multiset):
    '''Generator providing all multiset permutations of a multiset.'''
    h, i, j = init(multiset)
    yield visit(h)
    while j.next is not None or j.value < h.value:
        if j.next is not None and i.value >= j.next.value:
            s = j
        else:
            s = i
        t = s.next
        s.next = t.next
        t.next = h
        if t.value < h.value:
            i = t
        j = i.next
        h = t
        yield visit(h)

#%% Counting functions
        
        '''
        
        These is a block a of function that are using in the sequential importance sampler
        
        '''
        

def g_prob_het(oldsuff,F_nodes,family_size,indicator,total_vintage,n):
#def F_sample(oldsuff,seed):
    # current sample pool as list (set) of currently available nodes
        
    #np.random.seed(seed)
    #print(oldsuff)
    oldsuff['nodes']=np.array(oldsuff['nodes'])
    poolN, nodesDict, Noldsuff, original = createTree(oldsuff)
    L=np.max(oldsuff['nodes'][:,0:2])
   
    #generate paths
    pool, nodesDict0 = createTree_t0(oldsuff,L, isPF=True)
    n=sum(Noldsuff[Noldsuff[:,1]==0,2]*Noldsuff[Noldsuff[:,1]==0,3])
    leaf=Noldsuff[Noldsuff[:,4]==1,0]
    intern=Noldsuff[Noldsuff[:,4]==0,0]
    intern=np.append(intern,0)
    rel=vintParChil(F_nodes)
    child=childPerPhy(Noldsuff)
    nodeSize=NodeSize(Noldsuff)
    nodeInPath=NodePathFreq(Noldsuff,leaf,intern,nodeSize)
    family_group=family_groupGEN(Noldsuff,F_nodes)
    node_group=node_groupGEN(Noldsuff,family_group,total_vintage,indicator)
    node_group=np.array(node_group)
    paths=createPath_het(n,child,rel,nodeInPath,leaf,nodeSize,family_size,Noldsuff,family_group,node_group)

    
    #print(paths)
    #print(n)
    #n=5
    # initialize var to store result
    probs_steps = []
    probs = []
    
    for path in paths:
        
    #print(pool)
    #print(nodes)
        pool, nodes = createTree_t0(oldsuff,L)
        prob_p=prob_path_het(oldsuff,pool,nodes,path, family_size, F_nodes,Noldsuff)
        probs.append(np.prod(prob_p))
        probs_steps.append(prob_p)

    #print(node_group[0])
    # return F_mat,total_vintage,change_F 
    # FIX:
    # return dict, return part of F_mat
    return probs,probs_steps
# # Helper function
    


def prob_path_het(oldsuff,pool,nodes,path, family_size, F_nodes,Noldsuff):
    """
        sample available node from pool
        
        return F and path.
        
        In each iteration, one combination happens. Total nunmber of iterations = N.
    """
    prob_p=[]
    cnt = 0
    while len(pool) > 0:
#         Tracer()()
        L=np.max(oldsuff['nodes'][:,0:2])
        id_node=path[cnt]
        id_node2=[nod for nod in range(len(pool)) if pool[nod].name==id_node][0] #the [0] is just to transform it into an integer
        node=pool[id_node2]
        #Option 1 

        #node_prob = 1.0/len(pool)
        #prob_p.append(node_prob)
        #Option 2
        nn=[len(nodew.lineage[0])+len(nodew.lineage[1]) for nodew in pool]
        #if not any(nn):
        #     raise Exception("probability error!" +str(pool))
        node_prob = [1.0*pp/sum(nn) for pp in nn]
        node_prob=node_prob[id_node2]
        prob_p.append(node_prob)
#
        lineageStepProb_het(node,cnt,F_nodes,prob_p)
        updateTree(oldsuff,pool,nodes, cnt,L, isPF = False)
##         print (node.name, total_vintage, change_F)
        #print(node.lineage)
#        print('ok')
        #print(cnt)
#         print (F_mat)
        if len(node.lineage[0])+len(node.lineage[1])==1:
            # this node has only one lineage left, merge to parent
            # remove node from pool, add parent to pool
          
            pool.remove(node)  # O(n) with list
            #if node.numChilds==0:
            if node.name != 0:
                if (Noldsuff[Noldsuff[:,0]==node.name,2]* Noldsuff[Noldsuff[:,0]==node.name,3]==family_size[cnt]):    
                # if no child left 
                 #it could be possible that a node got free done, but still have children
                        # parent = node.parent
#                        parent.numChilds -= 1
#                        parent.lineage[1] += [cnt+1] # the children returns a box
#                        if len(node.lineage[1]) != 1:
#                            raise Exception('node not finished')
#                        parent.leaves[cnt+1] = node.leaves[node.lineage[1][0]]
#                        indicator[cnt] = 1
#                        # print ("node /{0}/ merge to parent /{1}/".format(node,parent))
#                        if len(parent.lineage[0])+len(parent.lineage[1])==2:
#                            pool.append(parent)
#                elif (oldsuff['nodes'][oldsuff['nodes'][:,0]==node.name,2]*oldsuff['nodes'][oldsuff['nodes'][:,0]==node.name,3])==family_size[-1]:
                    if node.name != -1:
                        parent = node.parent
                        parent.numChilds -= 1
                        parent.lineage[1] += [cnt+1] # the children returns a box
                        if len(node.lineage[1]) != 1:
                            raise Exception('node not finished')
                        parent.leaves[cnt+1] = node.leaves[node.lineage[1][0]]
                        # print ("node /{0}/ merge to parent /{1}/".format(node,parent))
                        if len(parent.lineage[0])+len(parent.lineage[1])==2:
                            pool.append(parent)
        cnt += 1  

    return prob_p
        


def lineageStepProb_het(node,cnt,F_nodes,prob_p):
    """
        Sample from lineage of one node, remove those two and combine
        
        Input: lineage, {0:[],1:[]}
        Return: types,ids
            types: two elements, each is among [0,1], denoting [free, constructed children or box].
            ids: index of two elements
    """
    lineage=node.lineage
    nLineage0 = len(lineage[0])
    nLineage1 = len(lineage[1])
    # chooseIds = np.random.choice(range(len(lineage[0])+len(lineage[1])),2,replace=False)
    values=np.array(list(Counter(lineage[0]).keys()))
    counts_pre = np.array(list(Counter(lineage[0]).values()))

    ids = []  # indexs of combined lineage
    types = []  # 0,1 indicating free or box
    if F_nodes[cnt][0]<=0:
        types.append(0)
        ids.append(F_nodes[cnt][0])
        lineage[0].remove(F_nodes[cnt][0])
    else:
        types.append(1)
        ids.append(F_nodes[cnt][0])
        lineage[1].remove(F_nodes[cnt][0])
        
    if F_nodes[cnt][1]<=0:
        types.append(0)
        ids.append(F_nodes[cnt][1])
        lineage[0].remove(F_nodes[cnt][1])
    else:
        types.append(1)
        ids.append(F_nodes[cnt][1])
        lineage[1].remove(F_nodes[cnt][1])
    lineage[1].append(cnt+1)  # add a new combined lineage
    counts_aft = [lineage[0].count(x) for x in set(values)]
    prob_num = np.prod([binomialCoeff(counts_pre[i],(counts_pre[i]-counts_aft[i])) for i in range(len(counts_pre))])
    prob_den = comb2(nLineage0+nLineage1)
    
    prob_p.append(prob_num/prob_den)
    ids.sort()
    new_family = node.leaves[ids[0]] + node.leaves[ids[1]]

    node.leaves[cnt+1] = new_family
    

#%%Kingman counting
    
    ''' Function to count Kingman heterochronous trees '''
    
    
def king_prob_het(oldsuff,F_nodes,family_size,indicator,total_vintage,n):
#def F_sample(oldsuff,seed):
    # current sample pool as list (set) of currently available nodes
        
    #np.random.seed(seed)
    oldsuff['nodes']=np.array(oldsuff['nodes'])
    poolN, nodesDict, Noldsuff, original = createTree(oldsuff)
    L=np.max(oldsuff['nodes'][:,0:2])
    path=total_vintage
    pool, nodes = createTree_t0(oldsuff,L)
    probs=[]
    cnt = 0
    while len(pool) > 0:
#         Tracer()()
        L=np.max(oldsuff['nodes'][:,0:2])
        id_node=path[cnt]
        id_node2=[nod for nod in range(len(pool)) if pool[nod].name==id_node][0] #the [0] is just to transform it into an integer
        node=pool[id_node2]
        #Option 1 
        #node_prob = 1.0/len(pool)
        #prob_p.append(node_prob)
        #Option 2
        nn=[len(nodew.lineage[0])+len(nodew.lineage[1]) for nodew in pool]
        #if not any(nn):
        #     raise Exception("probability error!" +str(pool))
        node_prob = [1.0*pp/sum(nn) for pp in nn]
        node_prob=node_prob[id_node2]
        probs.append(node_prob)
        KinglineageStepProb_het(node,cnt,F_nodes,probs)
        updateTree(oldsuff,pool,nodes, cnt,L, isPF = False)
   
        if len(node.lineage[0])+len(node.lineage[1])==1:
            # this node has only one lineage left, merge to parent
            # remove node from pool, add parent to pool
          
            pool.remove(node)  # O(n) with list
            #if node.numChilds==0:
           
            if node.name != 0:
                if (Noldsuff[Noldsuff[:,0]==node.name,2]* Noldsuff[Noldsuff[:,0]==node.name,3]==family_size[cnt]):    
      
#                elif (oldsuff['nodes'][oldsuff['nodes'][:,0]==node.name,2]*oldsuff['nodes'][oldsuff['nodes'][:,0]==node.name,3])==family_size[-1]:
                    if node.name != -1:
                        parent = node.parent
                        parent.numChilds -= 1
                        parent.lineage[1] += [cnt+1] # the children returns a box
                        if len(node.lineage[1]) != 1:
                            raise Exception('node not finished')
                        parent.leaves[cnt+1] = node.leaves[node.lineage[1][0]]
                        # print ("node /{0}/ merge to parent /{1}/".format(node,parent))
                        if len(parent.lineage[0])+len(parent.lineage[1])==2:
                            pool.append(parent)
        cnt += 1  
        

    return np.prod(probs),probs

def KinglineageStepProb_het(node,cnt,F_nodes,prob_p):
    """
        Sample from lineage of one node, remove those two and combine
        
        Input: lineage, {0:[],1:[]}
        Return: types,ids
            types: two elements, each is among [0,1], denoting [free, constructed children or box].
            ids: index of two elements
    """
    lineage=node.lineage
    nLineage0 = len(lineage[0])
    nLineage1 = len(lineage[1])
    # chooseIds = np.random.choice(range(len(lineage[0])+len(lineage[1])),2,replace=False)

    ids = []  # indexs of combined lineage
    types = []  # 0,1 indicating free or box
    if F_nodes[cnt][0]<=0:
        types.append(0)
        ids.append(F_nodes[cnt][0])
        lineage[0].remove(F_nodes[cnt][0])
    else:
        types.append(1)
        ids.append(F_nodes[cnt][0])
        lineage[1].remove(F_nodes[cnt][0])
        
    if F_nodes[cnt][1]<=0:
        types.append(0)
        ids.append(F_nodes[cnt][1])
        lineage[0].remove(F_nodes[cnt][1])
    else:
        types.append(1)
        ids.append(F_nodes[cnt][1])
        lineage[1].remove(F_nodes[cnt][1])
    lineage[1].append(cnt+1)  # add a new combined lineage
    prob_num = 1
    prob_den = comb2(nLineage0+nLineage1)
    
    prob_p.append(prob_num/prob_den)
    ids.sort()
    new_family = node.leaves[ids[0]] + node.leaves[ids[1]]

    node.leaves[cnt+1] = new_family

#%% Code check function
    
    '''
    Here we put some utilities to be used to check the quality of the code 

    
    '''


    '''
    Given F_nodes, I would like to get the total vintage and indicator function
 
    '''

def auxiliary_info(oldsuff,F_nodes,family_size):
    
    '''
    Given F_nodes and and family_size and the heterochronous perfect phylogeny (check)
    fill in the missing information needed for the likelihood calculation
    
    '''
    
    print(oldsuff)
    print(F_nodes)
    print(family_size)

    oldsuff['nodes']=np.array(oldsuff['nodes'])
    poolN, nodesDict, Noldsuff,original = createTree(oldsuff)
    L=np.max(oldsuff['nodes'][:,0:2])
       
    #generate paths
    pool, nodes = createTree_t0(oldsuff,L)
    n=sum(Noldsuff[Noldsuff[:,1]==0,2]*Noldsuff[Noldsuff[:,1]==0,3])
    
    
    total_vintage = []
    indicator = [0]*(n-1) 
    probs=[]
    
    
    cnt = 0


    while len(pool) > 0:
    #         Tracer()()
        print(cnt)
        L=np.max(oldsuff['nodes'][:,0:2])
        #understand which node I have sampled from
        poolComp=generate_poolComp(F_nodes,cnt,pool)
        node,node_prob = sampleNodeFromList(poolComp)
        probs.append(node_prob)
        total_vintage.append(node.name)
        lineageStepProb_het(node,cnt,F_nodes,probs)
        updateTree(oldsuff,pool,nodes, cnt,L, isPF = False)
    #         print (node.name, total_vintage, change_F)
        print('ok')
        #print(cnt)
    #         print (F_mat)
        if len(node.lineage[0])+len(node.lineage[1])==1:
            # this node has only one lineage left, merge to parent
            # remove node from pool, add parent to pool
            pool.remove(node)  # O(n) with list
            #if node.numChilds==0:
            if node.name != 0:
                if (Noldsuff[Noldsuff[:,0]==node.name,2]* Noldsuff[Noldsuff[:,0]==node.name,3]==family_size[cnt]):    
                # if no child left 
                 #it could be possible that a node got free done, but still have children
                        # parent = node.parent
    #                        parent.numChilds -= 1
    #                        parent.lineage[1] += [cnt+1] # the children returns a box
    #                        if len(node.lineage[1]) != 1:
    #                            raise Exception('node not finished')
    #                        parent.leaves[cnt+1] = node.leaves[node.lineage[1][0]]
    #                        indicator[cnt] = 1
    #                        # print ("node /{0}/ merge to parent /{1}/".format(node,parent))
    #                        if len(parent.lineage[0])+len(parent.lineage[1])==2:
    #                            pool.append(parent)
    #                elif (oldsuff['nodes'][oldsuff['nodes'][:,0]==node.name,2]*oldsuff['nodes'][oldsuff['nodes'][:,0]==node.name,3])==family_size[-1]:
                    if node.name != -1:
                        parent = node.parent
                        parent.numChilds -= 1
                        parent.lineage[1] += [cnt+1] # the children returns a box
                        if len(node.lineage[1]) != 1:
                            raise Exception('node not finished')
                        parent.leaves[cnt+1] = node.leaves[node.lineage[1][0]]
                        indicator[cnt] = 1
                        # print ("node /{0}/ merge to parent /{1}/".format(node,parent))
                        if len(parent.lineage[0])+len(parent.lineage[1])==2:
                            pool.append(parent)
        cnt += 1
    
    return total_vintage,indicator    
    
    
def generate_poolComp(F_nodes,cnt,pool): 
    
    poolCopy = pool[:]
    poolComp=[]
    for f in poolCopy:
        tot_lin =f.lineage[0]+list(f.lineage[1])
        if F_nodes[cnt][0] in tot_lin:
            tot_lin.remove(F_nodes[cnt][0])
            if F_nodes[cnt][1] in tot_lin:
                poolComp.append(f)
    
    return poolComp          
    

def factorial(n):
    return float(np.math.factorial(n))



def twin_datasets_finder(oldsuff,nodesDict,path,F_nodes):

    twin = 0
    
    pooltwin=[]
    idx = [i for i in range(len(oldsuff['nodes'])) if oldsuff['nodes'][i][2]==1 and oldsuff['nodes'][i][3]==2]
    
    for j in idx:
        if len(idx)>0:
            if np.sum(oldsuff['mylist'][j]['y'])==2:
                par_node = oldsuff['nodes'][j][1]
                pooltwin.append(nodesDict[par_node])
                
    
    pooltwin=set(pooltwin)
    for nod in pooltwin:
        name=nod.name
        if name == -1:
            continue
        locpaths=copy.copy(path)
        locpaths=(locpaths==name)*1
        
    
        if sum([len(nod.mutFromSingleton[i]) for i in list(nod.mutFromSingleton.keys())])>=1:
            
            table = Counter(nod.lineage[0])
            potential_cherries=[t for t in table.keys() if table[t]==2 and np.sum(nod.mutFromSingleton[t])==2]
            #potential_cherries=[t for t in table.keys() if table[t]==2]
            
            if len(potential_cherries)>0:
                for pot in potential_cherries:
                    if sum(nod.mutFromSingleton[pot])==2:
                         idd =[ii for ii in range(len(path)) if path[ii]==name]
                         for iddd in idd:
                             Fnode=F_nodes[iddd]
                             if Fnode[0]==Fnode[1] and Fnode[0]==pot:
                                twin+=1
        
    return twin    


