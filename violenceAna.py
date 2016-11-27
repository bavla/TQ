import os, re, pickle, json
from TQ import *

os.chdir('C:/Users/batagelj/work/Python/temporal')

# with open("violence.p","rb") as p: V = pickle.load(p)
V = TQ.Ianus2Mat("./ten/violence.ten")
W = V['mat']
n = len(W); All = range(n)
names = V['nam']
times = V['tin']
# listMatrix(W,names,False)
for i,name in enumerate(names): print(i,name)

# All-All
ca = TQ.ActInt(W,All,All)
# police - 3
cb = TQ.sum(TQ.ActInt(W,[3],All),TQ.ActInt(W,All,[3]))
# fascists - 6
cc = TQ.sum(TQ.ActInt(W,[6],All),TQ.ActInt(W,All,[6]))

comp = {}
comp['All-All']    = ca
comp['Pol-All']    = cb
comp['Fas-All']    = cc
#with open('violComp.json',mode='w',encoding='utf-8') as f:
#   json.dump(comp,f,indent=2)

print("ca =",ca,"\nsummary =",TQ.VecSummary(ca))
print("cb =",cb,"\nsummary =",TQ.VecSummary(cb))
print("cc =",cc,"\nsummary =",TQ.VecSummary(cc))

print('-'*50)
cc2 = TQ.clusCoef(W,type=2)
print("Clustering coefficient 2"); TQ.VecList(cc2,names)

#with open('violCC2.json',mode='w',encoding='utf-8') as f:
#   json.dump(cc2,f,indent=2)

p = TQ.PartMaxVec(cc2,All)
