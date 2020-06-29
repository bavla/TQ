# clustering temporal quantities
# by Vladimir Batagelj, June 2020

from numpy import random
import numpy as np
from copy import copy, deepcopy
import collections

def table(arr): return collections.Counter(arr)
def timer(): return datetime.datetime.now().ctime()
def unitTQ(nVar): return [ 'unit', [[[], 0]] * nVar ]

# unit = ['unit',[[[],0],[[],0]]]
infinity = float('inf')

# computes weighted squared Euclidean dissimilarity between TQs 
def distTQ(X,T,nVar,alpha):
  global infinity  
  d = 0
  for i in range(nVar):
    nX = X[i][1]; nT = T[i][1]
    if nT<=0: return infinity
    s = TQ.sum(TQ.prodConst(X[i][0],1/max(1,nX)),
               TQ.prodConst(T[i][0],-1/max(1,nT)))
    d = alpha[i]*nX*TQ.total(TQ.prod(s,s))
# if is.na(d): return(infinity)
  return d


# adapted leaders method
#---------------------------------------------
# VB, 16. julij 2010
#   dodaj omejitev na najmanjše število enot v skupini omejeni polmer
# VB, 1. February 2013: added: time, initial clustering
# Python version: VB, 26. June 2020
#  

def leaderTQ(TQs,maxL,nVar,alpha,clust=None,trace=2,tim=1,empty=0):
  global unit 
  numTQ = len(TQs)
  L = [ deepcopy(unit) for i in range(maxL)]
  Ro = np.zeros(maxL); K = np.zeros(maxL,dtype=int) 
# if not given, random partition into maxL clusters
  if clust is None:
    clust = np.random.choice(np.array(range(1,maxL-empty+1)),size=numTQ)  
  step = 0; print("LeaderTQ:",timer(),"\n\n")
  while tim > 0:
    step = step+1; K.fill(0)
  # new leaders - determine the leaders of clusters in current partition
    for j in range(maxL): L[j][0] = "L"+str(j+1)
    for i in range(numTQ):
      j = clust[i]-1
      for k in range(nVar):   
        L[j][1][k][0] = TQ.sum(L[j][1][k][0],TQs[i][1][k][0])
        L[j][1][k][1] = L[j][1][k][1]+TQs[i][1][k][1]
  # new partition - assign each unit to the nearest new leader
    clust.fill(0)
    R = np.zeros(maxL); p = np.zeros(maxL); d = np.zeros(maxL)
    for i in range(numTQ):     
      for k in range(maxL): d[k] = distTQ(TQs[i][1],L[k][1],nVar,alpha)
      j = np.argmin(d); r = d[j]
      if r == infinity:
        print("Infinite unit=",i,"\n"); print(TQs[i])
      clust[i] = j+1; p[j] = p[j] + r
      if R[j]<r: R[j] = r; K[j] = i
  # report
    if trace>0: print("Step",step,timer())
    delta = [ R[i]-Ro[i] for i in range(maxL) ]
    if trace>1:
      print(table(clust)); print(R); print(delta); print(p)
    if trace>0: print("P =",sum(p))
    if sum([abs(a) for a in delta])<0.0000001: break
    Ro = R; tim = tim-1
    if tim<1:
      print(table(clust)); print(R); print(delta); print(p); print("P =",sum(p)) 
      tim = int(input("Times repeat = :\n")); print(f'You entered {tim}')
  # TO DO: in the case of empty clusters use the most distant TQs as seeds
  return { 'proc':"leaderTQ", 'clust':clust.tolist(), 'leaders':L,
    'R':R.tolist(), 'p':p.tolist() }


# computes dissimilarity between clusters
def distCl(X,Y,nVar,alpha):
  d = 0
  for i in range(nVar):
    nX = X[i][1]; nY = Y[i][1]
    s = TQ.sum(TQ.prodConst(X[i][0],1/max(1,nX)),
               TQ.prodConst(Y[i][0],-1/max(1,nY)))
    d = d + alpha[i]*nX*nY/max(1,nX+nY)*TQ.total(TQ.prod(s,s))
  return d

# adapted Ward's hierarchical clustering
#---------------------------------------------
# VB, 16. julij 2010  Clamix in R
# Python version: VB, 23. junij 2020

def hclusTQ(L,nVar,alpha):
  global m
  def ordNode(S): return [ int(a) for a in S ]
  infinity = float('inf')  
  def orDendro(i):
    global m
    if i<0: return [-i]
    else: return orDendro(m[i-1,0])+orDendro(m[i-1,1])
  numL = len(L); numLm = numL-1
# each unit is a cluster; compute inter-cluster dissimilarity matrix
  D = np.zeros((numL,numL))
  for i in range(numL): D[i][i] = infinity
  for i in range(numL-1):
    for j in range(i+1,numL):      
      D[i][j] = distCl(L[i][1],L[j][1],nVar,alpha); D[j][i] = D[i][j]
  active = set(range(numL)); m = np.zeros((numLm,2),dtype=int) 
  node = np.zeros(numL,dtype=int); h = np.zeros(numLm)
  LL = []
  for k in range(numLm):
  # determine the closest pair of clusters (p,q)
    numA = len(active); dmin = infinity
    for a in active:
      for b in active:
        if D[a][b] < dmin:
          dmin = D[a][b]; p,q = a,b          
  # join the closest pair of clusters
    h[k] = dmin
    if node[p]==0: m[k][0] = -(p+1); Lp = L[p]
    else: m[k][0] = node[p]; Lp = LL[node[p]-1]
    if node[q]==0: m[k][1] = -(q+1); Lq = L[q]
    else: m[k][1] = node[q]; Lq = LL[node[q]-1]
    tq = []
    for t in range(nVar):
      tq.append([TQ.sum(Lp[1][t][0],Lq[1][t][0]), Lp[1][t][1]+Lq[1][t][1]])
    LL.append(['C'+str(k+1),tq[:]])  
    active.discard(p) 
  # determine dissimilarities to the new cluster
    for s in active:
      if s != q:
        if node[s]==0: Ls = L[s]
        else: Ls = LL[node[s]-1]
        D[q][s] = distCl(LL[k][1],Ls[1],nVar,alpha); D[s][q] = D[q][s]
    node[q] = k+1; 
#    orderNode = [ int(a) for a in order ]
  return {'proc':"hclusTQ", 'merge':m.tolist(), 'height':h.tolist(),
        'order': ordNode(orDendro(len(m))),
        'labels': [e[0] for e in L],
        'method':"hclusTQ", 'call':None, 'dist.method':"Ward", 'leaders':LL }
