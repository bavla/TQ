# support functions

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
