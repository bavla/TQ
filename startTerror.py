wdir = "C:/Users/batagelj/work/Python/graph/Nets/clusTQ"
gdir = 'c:/users/batagelj/work/python/graph/Nets'

import sys, os, re, datetime, json
sys.path = [gdir]+sys.path; os.chdir(wdir)
from TQ import *
from Nets import Network as N

import clusTQ as cl

G = N.loadNetsJSON("C:/Users/batagelj/work/Python/graph/JSON/terror/terror.json")
G.Info()
nVar = 1; alpha = [1]; unit = cl.unitTQ(nVar)
Ter = [[G._nodes[u][3]['lab'], [[G.TQnetSum(u), TQ.total(G.TQnetSum(u))]]]\
   for u in G.nodes() ]

# Rez = cl.leaderTQ(Ter,100,nVar,alpha,trace=1,tim=5)
# js = open("Terror100.json",'w'); json.dump(Rez, js, indent=1); js.close()

# Tot = [ (t[0],t[1][0][1]) for t in Ter ]
# js = open("Totals.json",'w'); json.dump(Tot, js, indent=1); js.close()

# HC = cl.hclusTQ(Rez['leaders'],nVar,alpha)
# js = open("TerrorHC.json",'w'); json.dump(HC, js, indent=1); js.close()

