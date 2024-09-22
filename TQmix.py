# support functions
# Sept 22, 2024

def width(tq):
   k = 0
   for i in tq:
      s,f,w = i; k = max(k,len(w[1]))
   return k

def span(tq):
   d = 0
   for i in tq: s,f,w = i; d = d + (f-s)
   return d

def third(t): return t[2]

def listNeighbors(R,u):
   v = R[u]["v"]; lab = R[u]["lab"]; tq = R[u]["tq"]
   print(u,". node:",lab)
   for i in tq:
      s,f,w = i; L = [int(x)-1 for x in w[1]]; L.sort()
      print(s,f,w[0],end="")
      for t in L: print(" ",R[t]["lab"],end="")
      print()

def stable(T,t=0):
   L = T._links
   S = [ (L[a][0],L[a][1],span(L[a][4]["tq"])) for a in L ]
   S.sort(reverse=True,key=third)
   for (u,v,d) in S:
      if d < t: break
      print(d,T._nodes[u][3]["lab"],">",T._nodes[v][3]["lab"]) 


