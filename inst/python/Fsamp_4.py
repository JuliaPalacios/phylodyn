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
import timeit
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
        self.mutFromSingleton = []
        self.leaves = {-1:1}    # map of number of leaves (free individuals) from each lineage
        self.family_size = 0
        self.original_node = 1.0
    def __repr__(self):
        return "TreeNode "+str(self.__hash__())[-5:]+' '+str(self.name)+': '+str(self.lineage)+' mutation: '+str(self.mutation)+' mut from sing: '+str(self.mutFromSingleton)+ 'original'+str(self.original_node)
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


# In[ ]:
def F_sample(oldsuff):
#def F_sample(oldsuff,seed):
    # current sample pool as list (set) of currently available nodes
        
    #np.random.seed(seed)

    pool, nodes = createTree(oldsuff)

    # In[ ]:
    # FIX: calculate actual individuals
    n = sum([len(x.lineage[0]) for x in nodes.values()])
    
    # initialize var to store result
    F_mat = np.diag(list(range(n,0,-1)))
    total_vintage = []
    change_F = []
    probs = []
    family_size = []
    F_nodes = []
    indicator = [0]*(n-1)
    #print("here ok")
    sampleF(pool,total_vintage, change_F, F_mat, probs, family_size, F_nodes, indicator)

    # return F_mat,total_vintage,change_F 
    # FIX:
    # return dict, return part of F_mat
    return {"F_mat":F_mat[:-1,:-1].tolist(),"total_vintage":total_vintage,"change_F":change_F, "probs":probs, "family_size":family_size, "F_nodes":F_nodes, "indicator":indicator}

# # Helper function

# In[ ]:
def createTree(data, isPF = False):

    # create tree bottom up
    #isPF is to speed (changing list to set), when sampling it is set False
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
        if row[4] == 1 and row[2]==1: #leaf and no merged
            # if is leaf, add parent node (but not strictly parent)
            # FIX: 
            # clamped individuals to parent when leaf and column 3 = 1
            if row[1] not in nodes:
                #print("this?"+str(row[0]))
                newnode = TreeNode(row[1], isPF)
                nodes[row[1]] = newnode 
                newnode.lineage[0] += [-1]*row[3]
                newnode.family_size += row[3]
                assert newnode.mutFromSingleton == [], 'mutFromSingleton initialized'
                newnode.mutFromSingleton = [mut['y']] if isinstance(mut['y'],int) else mut['y']
                newnode.original_node=1.0
            else:
                #print("or this?")
                newnode = nodes[row[1]]
                newnode.lineage[0] += [-1]*row[3]  
                newnode.family_size += row[3]
                assert newnode.mutFromSingleton == [], 'mutFromSingleton initialized'
                newnode.mutFromSingleton = [mut['y']] if isinstance(mut['y'],int) else mut['y']
                newnode.original_node=1.0
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
            #print('DEBUG', i)
            newnode.original_node=original[i-1]
            if row[1] not in nodes:
                parent = TreeNode(row[1], isPF)
                nodes[row[1]] = parent
            else:
                parent = nodes[row[1]]
            parent.numChilds += 1
            newnode.parent = parent
            if row[4]==1:
                # leaf node with column 3 >1, has more than one free
                newnode.lineage[0] += [-1]*row[2]
                newnode.family_size += row[2]
                pool.append(newnode)
            parent.family_size += newnode.family_size
            nodes[row[0]] = newnode
            #print(nodes[row[0]].original_node)
            #nodes[row[0]].original_node=original[row[0]]
        #print (nodes[i-1].name)
        #print(nodes[i-1].original_node)
    nodes[-1] = TreeNode(-1,isPF)  ##should be back
    nodes[0].parent = nodes[-1] ##should be back
    nodes[-1].family_size = nodes[0].family_size ##should be back
    
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
    
def calcFamilySample(F_nodes,change_F, family_size):
    """
        compute upward family_size list for each node in sample
    """
    res = collections.defaultdict(set,set())
    for ind,pair in reversed(list(enumerate(F_nodes))):
        cf,fs = change_F[ind],family_size[ind]
        res[ind].add(fs)
        if cf != 0:
            res[pair[0]-1] = res[pair[0]-1].union(res[ind])
        else:
            res[pair[0][0]-1] = res[pair[0][0]-1].union(res[ind])
            res[pair[1][0]-1] = res[pair[1][0]-1].union(res[ind])
    res.pop(-2)
    return res

def comb2(total):
    # return number of choices of 2 from total
    return (total-1)*total/2.0



# In[ ]:



def sampleNodeFromList(nodeList, cnt, old_sample=None, start_changing=None):
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
    if old_sample is not None and cnt + 1 < start_changing:
        choice = None
        for i in range(len(nodeList)):
            if nodeList[i].name == old_sample['total_vintage'][cnt]:
                choice = nodeList[i]
                break
        assert choice is not None
    else:
         #Non uniform - proportional on the sample size
         p = [len(node.lineage[0])+len(node.lineage[1]) for node in nodeList]
         #print(p)
         if not any(p):
             raise Exception("probability error!" +str(nodeList))
         prob = [1.0*pp/sum(p) for pp in p]
         #choice = np.random.choice(nodeList, p=prob)
         idx=np.random.choice(len(nodeList), p=prob)
         choice=nodeList[idx]
         #choice_prob = [pp for pp,node in (prob, nodeList) if node==choice]
         choice_prob=prob[idx]
         #if len(int(choice_prob))!=1:
          #  raise Exception("choice prob error!")
       
    return choice, choice_prob
         #Uniform
#         choice = nodeList[np.random.randint(len(nodeList))]
#    return choice, 1.0/len(nodeList)


def sampleFromLineage(lineage,cnt):
    """
        Sample from lineage of one node, remove those two and combine
        
        Input: lineage, {0:[],1:[]}
        Return: types,ids
            types: two elements, each is among [0,1], denoting [free, constructed children or box].
            ids: index of two elements
    """
    choices = [0,1]
    nLineage0 = len(lineage[0])
    nLineage1 = len(lineage[1])
    chooseIds = random.sample(range(len(lineage[0])+len(lineage[1])),2)
    # chooseIds = np.random.choice(range(len(lineage[0])+len(lineage[1])),2,replace=False)
    chooseNums = [lineage[int(cid>=nLineage0)][cid-(cid>=nLineage0)*nLineage0] for cid in chooseIds]
    ids = []  # indexs of combined lineage
    types = []  # 0,1 indicating free or box
    prob_choices = [comb2(len(lineage[0])), len(lineage[0])*len(lineage[1]), comb2(len(lineage[1]))]    # number of choices of [0,0],[0,1],[1,1]
    for cid,cnum in zip(chooseIds,chooseNums):
        ids.append(cnum)
        types.append(int(cid>=nLineage0))
        lineage[int(cid>=nLineage0)].remove(cnum)
    lineage[1].append(cnt+1)  # add a new combined lineage
    return types,ids, 1.0*prob_choices[sum(types)]/sum(prob_choices), nLineage1


# In[ ]:

def combineAndUpdate(node, cnt, total_vintage, change_F, F_mat, probs, family_size, F_nodes):
    """
        Given a node, sample two of lineages and combine. 
        Update total_vintage etc to reflect the combination.
        
        Return nothing. Update in-place.
        
        Time complexity: with a good data structure, O(1) remove, add time.
        To update F_mat and path, O(N)??
    """
    # sample two lineage ang update lineage
    types,ids,prob,origin_nLin1 = sampleFromLineage(node.lineage,cnt)
    probs.append(prob)
    ids.sort()
    
    # update F_mat and path
    total_vintage.append(node.name)
    new_family = node.leaves[ids[0]] + node.leaves[ids[1]]
    family_size.append(new_family)
    node.leaves[cnt+1] = new_family
    if [0,0] == types:
        # all are free individual
        change_F.append(2)
        F_mat[cnt+1][:cnt+1] = F_mat[cnt][:cnt+1]-2
        F_nodes.append((-1, cnt+1))
    elif 0 in types:
        # one free and one box
        change_F.append(1)
        F_mat[cnt+1][:ids[1]] = F_mat[cnt][:ids[1]]-1
        F_mat[cnt+1][ids[1]:cnt+1] = F_mat[cnt][ids[1]:cnt+1]-2
        probs.append(1.0/origin_nLin1)
        F_nodes.append((ids[1], cnt+1))
    else:
        # two boxes
        change_F.append(0)
        F_mat[cnt+1][:cnt+1] = F_mat[cnt][:cnt+1]
        F_mat[cnt+1][ids[0]:ids[1]] -= 1
        F_mat[cnt+1][ids[1]:cnt+1] -= 2
        probs.append(1.0/comb2(origin_nLin1))
        F_nodes.append(((ids[0], cnt+1), (ids[1], cnt+1)))
        # F_nodes.append([ids[1], cnt+1])
    


# # Sample function

# In[ ]:

def sampleF(pool, total_vintage, change_F, F_mat, probs, family_size, F_nodes, indicator):
    """
        sample available node from pool
        
        return F and path.
        
        In each iteration, one combination happens. Total nunmber of iterations = N.
    """
    cnt = 0
    while len(pool) > 0:
#         Tracer()()
        node,node_prob = sampleNodeFromList(pool)
        probs.append(node_prob)
        combineAndUpdate(node, cnt, total_vintage, change_F, F_mat, probs, family_size, F_nodes)
        
#         print (node.name, total_vintage, change_F)
#         print (F_mat)
        if len(node.lineage[0])+len(node.lineage[1])==1:
            # this node has only one lineage left, merge to parent
            # remove node from pool, add parent to pool
            pool.remove(node)  # O(n) with list
            if node.numChilds==0:  
                # if no child left 
                # it could be possible that a node got free done, but still have children
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
#         print (pool)
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

def possibleUpdate(node, changeF, Fnode, ind):
    newprob = None
    oldnode = None
    prob_choices = [comb2(len(node.lineage[0])), \
                    len(node.lineage[0])*len(node.lineage[1]),\
                    comb2(len(node.lineage[1]))]
    prob_choices = [prob*1.0/sum(prob_choices) for prob in prob_choices]
    if changeF == 2:
        if len(node.lineage[0])>1:
            oldnode = node.copy()
            node.lineage[0].pop()
            node.lineage[0].pop()
            node.lineage[1].add(ind+1)
            newprob = prob_choices[0]
    elif changeF == 1:
        if len(node.lineage[0])>0 and len(node.lineage[1])>0 and \
            ind+1==Fnode[1] and Fnode[0] in node.lineage[1]:
            # can produce the same Fnode (need to verify one Fnode - index of the box)
            oldnode = node.copy()
            node.lineage[0].pop()
            node.lineage[1].remove(Fnode[0])
            node.lineage[1].add(ind+1)
            newprob = prob_choices[1]/len(oldnode.lineage[1])
    elif changeF == 0:
        if len(node.lineage[1])>1 and \
            Fnode[0][0] in node.lineage[1] and Fnode[1][0] in node.lineage[1] and \
            ind+1==Fnode[0][1]:
            # can produce the same Fnode (need to verify two Fnodes - index of two boxes)
            oldnode = node.copy()
            node.lineage[1].remove(Fnode[0][0])
            node.lineage[1].remove(Fnode[1][0])
            node.lineage[1].add(ind+1)
            newprob = prob_choices[2]/comb2(len(oldnode.lineage[1]))
    return newprob is not None, oldnode, newprob, node
    
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
def specialMap(muts,bls,cfs):
    key = tuple(sorted(cfs)+muts)
    #maps[key] will store the coefficient 
    if key not in maps:
        perms = set(list(itertools.permutations(muts)))
        perms = sorted(list(perms))
        n_cf2 = len(bls[2])*2
        n_cf1 = len(bls[1])
        coeff_1 = coeff_2 = 0
        B = []
        for perm in perms:
            if sorted(perm[n_cf2:]) not in B:
                B.append(sorted(perm[n_cf2:]))
                # permute perm[n_cf2:] for cf=1 to choose from
                permB = set(list(itertools.permutations(perm[n_cf2:])))
                n_repB = len(permB)
                # permute perm[:n_cf2] for cf=2 to choose from
                permA = set(list(itertools.permutations(perm[:n_cf2])))
                A = []
                if n_cf2>0:
                    for sp in permA:
                        assert len(sp)>=2, "permA length wrong"
                        if sorted(sp[:2]) not in A:
                            A.append(sorted(sp[:2]))
                n_repA = len(A)
                if n_cf1>0:
                    coeff_1 += n_repA*sum(pb[0] for pb in permB)
                if n_cf2>0:
                    coeff_2 += n_repB*sum(pa[0]+pa[1] for pa in A)
        maps[key] = [coeff_1, coeff_2]
    else:
        coeff_1,coeff_2 = maps[key]
    return coeff_1*sum(np.log(bls[1]))+coeff_2*sum(np.log(bls[2]))

def compare(sp,a):
    # compare if sp is distinct from a
    assert len(sp)==len(a)
    for i in range(len(sp)//2):
        if sorted(sp[i*2:i*2+2]) != sorted(a[i*2:i*2+2]):
            return True
    return False

def singFunc(muts, bls, inds, L):
    # modify grad in place
    
    sumpermLik = 0 # summation over likelihood of all permutations

    perms = set(list(itertools.permutations(muts)))
    perms = sorted(list(perms))
    #n_cf2 = len(bls[2])*2
    #print("perms"+str(len(perms)))
    #print("cherries"+str(n_cf2))
    #print("branches"+str(bls))
    #print("individuals"+str(inds))
    #n_cf1 = len(bls[1])
    #print("singletons"+str(n_cf1))
    #coeff_1 = coeff_2 = 0
    #B = []
    #A = []
    #lastCf1 = perms[0][:n_cf1]
    #print("lastCf1"+str(perms[0][:n_cf1]))
    grad = np.zeros(L)
    for perm in perms:
        #print("current perm:"+str(perm))
        s = 0
        grad1=np.zeros(len(grad))
        for i in range(len(bls)):
            grad1[:inds[i]+1]+=1.0*perm[i]/bls[i]
            #print("intermediate grad"+str(grad1)+"for i "+str(i)+"with mut"+str(perm[i])+"and individual"+str(inds[i])+"and branch"+str(bls[i]))
            s += 1.0*np.log(bls[i])*perm[i] 
        grad+=np.exp(s)*grad1
        sumpermLik += 1.0*np.exp(s)
        #print("supermLik is "+str(sumpermLik)) 
        #print("grad is"+str(grad))
        #for (bl, mut) in zip(bls,perm):
            #print("mut"+str(mut)+"bl"+str(bl))
            #s+=mut*np.log(bl)
            #grad1[inds[i]]=
         #grad+=s*grad1   


        # if muts for cf=1 get another permute, process the current batch of perms with same lastCf1
        # if lastCf1 != perm[:n_cf1]:
        #     print("this I need to understand, why the next assert?")
        #     # for all perms for cf=2
        #     assert len(lastCf1) == n_cf1
        #     assert len(inds[1]) == n_cf1
        #     s_cf1 = sum(xi*np.log(bi) for xi,bi in zip(lastCf1,bls[1]))
        #     print("xi"+str(lastCf1)+"branch:"+str(bls[1])+"scf1"+str(s_cf1)) ##singleton branch
        #     grad_cf1 = np.zeros(len(grad))
        #     for i in range(n_cf1):
        #         grad_cf1[:inds[1][i]+1] += 1.0*lastCf1[i]/bls[1][i]
        #         #print(str(i)+"grad_cf1:"+str(lastCf1[i])+"length:"+str(bls[1][i]))
        #         #print(str(i)+"grad_cf1:"+str(1.0*lastCf1[i]/bls[1][i]))
        #     for permA in A:
        #         print("permA first"+str(permA))
        #         print("s_cf1 first"+str(s_cf1))
        #         print("bls first"+str(bls[2])+"length:"+str(len(bls[2])))
        #         #s = np.exp(s_cf1+sum(np.log(bls[2][i])*(permA[2*i]+permA[2*i+1]) for i in range(len(bls[2]))))
        #         print("things are not being corrected here")
        #         s = np.exp(s_cf1+sum(np.log(bls[2][i])*(permA[2*i]+permA[2*i+1]) for i in range(len(bls[2])))) 
        #         ##print("s is:"+str(s))
        #         grad_this = np.zeros(len(grad))
        #         assert len(permA) == n_cf2
        #         for i in range(len(bls[2])):
        #             #grad_this[:inds[2][i]+1] += 1.0*(permA[2*i]+permA[2*i+1])/bls[2][i] 
        #             grad_this[:inds[2][i]+1] += 1.0*(permA[2*i]+permA[2*i+1])/bls[2][i] #this is my correction for cherries               
        #         grad += s*(grad_this+grad_cf1)
        #         sumpermLik += s
            
        #     A = []
        #     lastCf1 = perm[:n_cf1]
            
        
        # sp = perm[n_cf1:]
        # print("sp"+str(sp))
        # if len(A)==0 or all([compare(sp, A[i]) for i in range(len(A))]):
        #     print("this compare thing is used")
        #     A.append(sp)
        #     print("now A is"+str(A))
            
    # accumulate the last part of gradient
    #s_cf1 = sum(xi*np.log(bi) for xi,bi in zip(lastCf1,bls[1]))
    #s_cf1 = sum(xi for xi,bi in zip(lastCf1,bls[1]))
    #print("scf1 to see:"+str(s_cf1))
    # grad_cf1 = np.zeros(len(grad))
    # for i in range(n_cf1):
    #     grad_cf1[:inds[1][i]+1] += 1.0*lastCf1[i]/bls[1][i]
    # for permA in A:
    #     print("second perm"+str(permA)+"needs correction")
    #     s = np.exp(s_cf1+sum(np.log(bls[2][i])*(permA[2*i]+permA[2*i+1]) for i in range(len(bls[2]))))
    #     #s = s_cf1+sum(np.log(bls[2][i])*(permA[2*i]+permA[2*i+1]) for i in range(len(bls[2]))) #Julia: changed to log
    #     #print("s_cf1:"+str(s_cf1))
    #     #print("bls"+str(bls[2][i]))
    #     #print("s line 566: "+str(s))
    #     grad_this = np.zeros(len(grad))
    #     assert len(permA) == n_cf2
    #     for i in range(len(bls[2])):
    #         print("being corrected the gradient here")
    #         grad_this[:inds[2][i]+1] += 1.0*(permA[2*i]+permA[2*i+1])/bls[2][i]     
    #         #print("index"+str(i)+"grad is:"+str(grad_this)) 
    #         #print("grad_tis:"+str(grad_this))          
        
    #     # if supermLik==0: 
    #     #     supermLik=s
    #     #     grad=grad
    #     # else:
    #     #     grad=(1/(np.exp(supermLik)+np.exp(s)))*(np.exp(supermLik)*res[2]+np.exp(s)*grad)
    #     #     supermLik=np.log(np.exp(supermLik)+np.exp(s))

    #     grad += s*(grad_this+grad_cf1)
    #     #print("s is:"+str(s)+"gradient tot"+str(grad))
    #     sumpermLik += s
    #     #print("supermlink:"+str(sumpermLik))
    
    return sumpermLik,grad
        

    
def processMutation(nodesDict, this_vintage, change_F, timeVec, indicator, tMatrix):
    # gradient wrt ts
    grad = np.zeros(timeVec.shape[0])
    # sum log \sum_v lik_v
    loglik = 0
    #print("calling process Mut")



    for name,node in nodesDict.items():

        #print("the node right now is"+str(name))
        if name == -1:
            continue
            
        # at what point the node got marginalized
        #print('DEBUG', [tv==name and idc==1 for tv,idc in zip(this_vintage,indicator)])
        ind_marginalize = [tv==name and idc==1 for tv,idc in zip(this_vintage,indicator)].index(True)
        
        # extending branch: mutation * log(branch length extending)
        if node.name != 0:
            #print("node name"+str(node.name)+" mut"+str(node.mutation))
            ind_nonZero = tMatrix[ind_marginalize]
            loglik += node.mutation*np.log(timeVec[ind_marginalize][1])
            #print("loglik other:"+str(loglik))
            #print ('extending mutation:{}, loglik:{}, len:{}'.format(node.mutation, node.mutation*np.log(timeVec[ind_marginalize][1]), timeVec[ind_marginalize][1]))
            grad[ind_nonZero==True] += 1.0*node.mutation/timeVec[ind_marginalize][1]
            #print("gradient other:"+str(grad))
        #print("mutations:"+str(node.mutation)+"time:"+str(timeVec[ind_marginalize][1])+str(ind_nonZero)+"gradient:"+str(grad))
        
        if len(node.mutFromSingleton)==1:
            # branch from singleton: mutFromSingleton * log(branch length singleton)
            loglik += node.mutFromSingleton[0]*np.log(timeVec[ind_marginalize][0])
            #print("loglik from singletons:"+str(loglik))
            #print ('sing in parent mutation:{}, loglik:{}, len:{}'.format(node.mutFromSingleton[0], node.mutFromSingleton[0]*np.log(timeVec[ind_marginalize][0]), timeVec[ind_marginalize][0]))
            grad[:ind_marginalize+1] += 1.0*node.mutFromSingleton[0]/timeVec[ind_marginalize][0]
            #print("mutations sing:"+str(node.mutFromSingleton[0])+"time:"+str(timeVec[ind_marginalize][0])+"gradient"+str(grad))
            #print("gradient from singletons:"+str(grad))
            
        elif len(node.mutFromSingleton)>1:
            # likelihood from singletons
            branches = {1:[],2:[]}  #branches[cf] = [list of branches at change_F=cf]
            #cfs = []
            #inds = {1:[],2:[]}
            branch_vec = []
            inds_vec=[]
            #print("indnonzero"+str(tMatrix[ind_marginalize]))
            #print("ind_marginalize"+str(ind_marginalize))
            for i,tv in enumerate(this_vintage):
                #print("i index is"+str(i)+"tv is"+str(tv))
                if tv != name or change_F[i] == 0:
                    continue
                #branches[change_F[i]] += [timeVec[i][0]]
                inds_vec +=[i]*change_F[i]
                #print("inds_vec"+str(inds_vec))
                #cfs.append(change_F[i])
                #inds[change_F[i]].append(i)
                branch_vec += [timeVec[i][0]]*change_F[i]
                #print("old inds"+str(inds))
            #loglik += specialMap(node.mutFromSingleton,branches,cfs)
            #print("This is the part singFunc is called")
            this_l,this_g = singFunc(node.mutFromSingleton,branch_vec,inds_vec,len(grad))
            loglik += np.log(this_l)
            #print("loglik is "+str(loglik))
            grad += 1.0*this_g/this_l
            #print("gradient is "+str(grad))
    #return loglik,grad/loglik -np.arange(len(grad))  
    return loglik,grad    # k start from 0?
            
    
def pF(nodesDict, pool, change_F, F_nodes, ind, prob, res, this_vintage, indicator, timeVec = [], tMatrix = [], fs_sample=None, fs_data=None, no_output = False):
    """
        iterate through nodes to get probability of F_mat (defined by change_F+F_nodes)
    """
    #print("ind"+str(ind))
    # Tracer()()
    if ind >= len(change_F): #index of vintages, means end of a path
        # find one complete sample, accumulate probability, access indicator and timeVec
        res[0]+= prob
        #print("res 0 is"+str(prob))
        if not no_output:
            #print("computes here")
            #start = timeit.timeit()
            loglik,grad = processMutation(nodesDict, this_vintage, change_F, timeVec, indicator, tMatrix) #Julia: changed this part
            #end = timeit.timeit()
            #print end-start
            if res[1]==0: 
                res[1]=loglik
                res[2]=grad
            else:
                res[2]=(1/(np.exp(res[1])+np.exp(loglik)))*(np.exp(res[1])*res[2]+np.exp(loglik)*grad)
                res[1]=np.log(np.exp(res[1])+np.exp(loglik))

            #res[1]+=loglik
            #print("res is:"+str(res[1]))
            #print("res grad is:"+str(res[2]))
            #print ('found one!'+str(this_vintage)+" loglik from mutation: "+str(loglik)+"and the total"+str(res[1]))
            #res[2] += grad
        return
    if len(pool) == 0:
        # not found
        print ('not found! empty pool')
        return 
    L = len(pool)
    nodeProb = 1.0/L
    newpool = set([node for node in pool])  #something for speeding up
    n = len(change_F)
    for oldnode in pool:
        #print("oldnode name"+str(oldnode.name))
        # for all nodes in pool
        # check if fs_data[oldnode.name] is present in fs_sample[ind] fs=family size
        if fs_data is not None:
            compatible = True
            for fs in fs_data[oldnode.name]:
                if fs not in fs_sample[ind]:
                    compatible = False
                    break
            if not compatible:
                continue
        compatible, oldnode_copy, newprob, newnode = possibleUpdate(oldnode, change_F[ind], F_nodes[ind], ind)
        if compatible:
#             print (oldnode_copy, newnode, change_F[ind], F_nodes[ind])
#             print ('pool:'+ str(pool))
            # there is only one possible newnode&newprob
            oldinfo, isMarginalize = update(newpool, newnode, ind)
            #print("oldsinfo"+str(oldinfo))
            #print("the value of ind"+str(ind))
            pF(nodesDict, newpool, change_F, F_nodes, ind+1, prob*newprob*nodeProb, res, this_vintage+[oldnode.name], indicator+[isMarginalize],timeVec,tMatrix, fs_sample, fs_data, no_output)
            if not oldinfo[0]:
                # add old node (oldinfo[0])
                newpool.add(oldnode)
            # always restore old lineage
            oldnode.restoreLineage(oldnode_copy)
            if oldinfo[1]:
                # restore parent (oldinfo[1])
                oldnode.parent.restoreLineage(oldinfo[1])
            if oldinfo[2]:
               # remove newparent from pool
                newpool.remove(oldinfo[2])
    #if n > 15 and ind < 3 and not no_output:
    #    print ('search done @ dep='+str(ind)+' with '+str(this_vintage))
    return 

def calcTimeVec(F_nodes, change_F, times):
    # calculate time vectors given F_nodes and change_F
    # vec[i][0] is branch length of singletons, vec[i][1] is extending length
    vec = [[sum(times[:ind+1]),-sum(times[:ind+1])] for ind in range(len(times))]
    tMatrix = np.zeros((len(times),len(times)),dtype=bool)
    #print("defaul tMatrix"+str(tMatrix))

    for cf, fn in zip(change_F,F_nodes):
        if cf==2:
            continue
        elif cf==1:
            vec[fn[0]-1][1] += vec[fn[1]-1][0]
            tMatrix[fn[0]-1][fn[0]:fn[1]] = True # these ts are involved for the extending branch
        else:
            vec[fn[0][0]-1][1] += vec[fn[0][1]-1][0]
            tMatrix[fn[0][0]-1][fn[0][0]:fn[0][1]] = True
            vec[fn[1][0]-1][1] += vec[fn[0][1]-1][0]
            tMatrix[fn[1][0]-1][fn[1][0]:fn[0][1]] = True
        #print("tMatrix"+"fornode"+str(fn)+"and mat"+str(tMatrix))        
    return np.array(vec),tMatrix

def calcPF(change_F, F_nodes, family_size, oldsuff, times=None, use_familySize = True, no_output = False):
    """ 
        Usage:
        pf <- python.call("calcPF", result$change_F, result$F_nodes, result$family_size, oldsuff) 
    """

    pool, nodesDict = createTree(oldsuff, isPF=True)

    #print("ok?")
    #print("pool"+str(pool))
    ind = 0
    prob = 1
    correction = 0.0
    fs_data = calcFamilyNodes(nodesDict)
    fs_sample = calcFamilySample(F_nodes, change_F, family_size)
    
    assert isinstance(times,list)

    if times is not None:
        timeVec,tMatrix = calcTimeVec(F_nodes,change_F,times)
        #if not no_output:
        #print ("Time Vec: "+str(timeVec[0]))
        #print ("T Matrix: "+str(tMatrix[0]))
    
    else:
        timeVec = []
        tMatrix = []
    res = [0,0,np.zeros(timeVec.shape[0])]
    #print("what res"+str(res))
    if use_familySize:
        pF(nodesDict, pool, change_F, F_nodes, ind, prob, res, this_vintage=[], indicator=[], timeVec=timeVec, tMatrix = tMatrix, fs_sample=fs_sample, fs_data=fs_data, no_output = no_output)
    else:
        pF(pool, change_F, F_nodes, ind, prob, res, this_vintage=[], indicator=[], no_output=no_output)
    
    #print("number of nodes"+str(nodesDict[0].original_node)+"and length"+str(len(nodesDict)))
    for name,node in nodesDict.items():
        #print("name"+str(name)+"and original"+str(node.original_node))
        correction+=np.log(node.original_node)

    #print("number of nodes"+str(nodesDict[-1].original_node)+"and length"+str(nodesDict.items()))
    #print("number of nodes"+str(nodesDict[-2].original_node)+"and length"+str(len(nodesDict)))
    #print("number of nodes"+str(nodesDict[-1].value)+"and length"+str(len(nodesDict)))
    #correction = sum(np.log(nodesDict[x-1].original_node) for x in range(len(nodesDict)))
    #print("correction"+str(correction))
    
    #print("nodes"+str(nodes[2].original_node))
    
    # FIX: subtract terms of times
    #print ("what res 1 is here")
    res[1] -= np.sum(np.arange(len(res[2]) + 1, 1, -1) * np.array(times))
    res[1]-=correction
    res[2] -= np.arange(len(res[2]) + 1, 1, -1)

    res[2] = list(res[2])
    #if not no_output:
     #   print ('map of mutation: '+str(maps))
    #print("what res"+str(res))
    #res[1]=loglik
    return res
