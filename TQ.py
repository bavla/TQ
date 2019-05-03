'''
Created on October 29, 2015 / March 30, 2015 / August 15, 2014 / March 10, 2014
Updated on November 11, 2018
Apr 29 2019: added changeTime

@author: Vladimir Batagelj, Selena Praprotnik
'''
# operations with temporal quantities (tau = 0)

from copy import deepcopy
import pickle, random, operator, re, datetime, json
from math import sqrt

class TQ(object):

   class TQerror(RuntimeError): pass

   ''' constants '''
   inf = float("inf")
   sAdd = operator.add
   sMul = operator.mul
   sZero = 0; sOne = 1; sN = []; sE = [(1,inf,1)]
   rPF = 2

   ''' operations '''
   @staticmethod
   def pitagora(a,b): return(sqrt(a*a+b*b))

   @staticmethod
   def Minkowski():
      if TQ.rPF==TQ.inf: return(max)
      if TQ.rPF==1: return(operator.add)
      if TQ.rPF==2: return(TQ.pitagora)
      else: return(lambda a,b: (a**TQ.rPF+b**TQ.rPF)**(1/TQ.rPF))

   @staticmethod
   def geoAdd(a,b):
      (av,ac) = a; (bv,bc) = b
      return((min(av,bv), ac if av<bv else ac+bc if av==bv else bc))

   @staticmethod
   def geoMul(a,b):
      (av,ac) = a; (bv,bc) = b
      return((av+bv,ac*bc))

   ''' temporal semirings '''
   @staticmethod
   def combinatorial():
      TQ.sAdd  = operator.add;    TQ.sMul = operator.mul
      TQ.sZero = 0;               TQ.sOne = 1
      TQ.sN    = [];              TQ.sE   = [(1,TQ.inf,1)]
      TQ.semiring = TQ.combinatorial

   @staticmethod
   def path():
      TQ.sAdd  = min;             TQ.sMul = operator.add
      TQ.sZero = TQ.inf;          TQ.sOne = 0
      TQ.sN    = [];              TQ.sE   = [(1,TQ.inf,0)]
      TQ.semiring = TQ.path

   @staticmethod
   def maxmin():
      TQ.sAdd  = max;             TQ.sMul = min
      TQ.sZero = -TQ.inf;         TQ.sOne = TQ.inf
      TQ.sN    = [];              TQ.sE   = [(1,TQ.inf,TQ.inf)]
      TQ.semiring = TQ.maxmin

   @staticmethod
   def reach():
      TQ.sAdd  = operator.or_;    TQ.sMul = operator.and_
      TQ.sZero = 0;               TQ.sOne = 1
      TQ.sN    = [];              TQ.sE   = [(1,TQ.inf,1)]
      TQ.semiring = TQ.reach

   @staticmethod
   def PFsemi():
      TQ.sAdd  = min;             TQ.sMul = TQ.Minkowski()
      TQ.sZero = TQ.inf;          TQ.sOne = 0
      TQ.sN    = [];              TQ.sE   = [(1,TQ.inf,0)]
      TQ.semiring = TQ.PFsemi

   @staticmethod
   def geodetic():
      TQ.sAdd  = TQ.geoAdd;       TQ.sMul = TQ.geoMul
      TQ.sZero = (TQ.inf,0);      TQ.sOne = (0,1)
      TQ.sN    = [];              TQ.sE   = [(1,TQ.inf,TQ.sOne)]
      TQ.semiring = TQ.geodetic

   semiring = combinatorial

   @staticmethod
   def report():
      print("\nsemiring = ",TQ.semiring.__name__)
      print("add   = ",TQ.sAdd.__name__)
      print("mult  = ",TQ.sMul.__name__)
      print("sZero = ",TQ.sZero)
      print("sOne  = ",TQ.sOne)
      print("sN    = ",TQ.sN)
      print("sE    = ",TQ.sE)
      print("rPF   = ",TQ.rPF,"\n")

   ''' return the next triplet in a temporal quantity or a sentry value (inf,inf,0) '''
   @staticmethod
   def get(S):
      try: return(next(S))
      except StopIteration: return((TQ.inf,TQ.inf,TQ.sZero))

   ''' temporal quantities operations '''
   @staticmethod
   def standard(a):
      if len(a) == 0: return(a)
      fc = sc = TQ.inf; c = []
      vc = 0
      for (sa,fa,va) in a:
         if fc == sa:
            if vc == va: fc = fa
            else:
               if sc != fc: c.append((sc,fc,vc))
               vc = va; sc = sa; fc = fa
         else:
            if (fc != TQ.inf) and (sc != fc): c.append((sc,fc,vc))
            sc = sa; fc = fa; vc = va
      if sc != fc: c.append((sc,fc,vc))
      return(c)

   @staticmethod
   def isTq(a):
      if isinstance(a, list) and all(isinstance(e, tuple) for e in a): return True
      return False

   @staticmethod
   def TqSummary(a):
      vMin = tMin = TQ.inf; vMax = tMax = -TQ.inf
      for (sa,fa,va) in a:
         if sa < tMin: tMin = sa
         if fa > tMax: tMax = fa
         if va < vMin: vMin = va
         if va > vMax: vMax = va
      return((tMin,tMax,vMin,vMax))

   @staticmethod
   def TqMax(a):
      s = [ va for (sa,fa,va) in a ]
      return -TQ.inf if s == [] else max(s)

   @staticmethod
   def binary(a):
      return TQ.standard([ (sa,fa,1) for sa,fa,va in a if va > 0 ])

   @staticmethod
   def setConst(a,c):
      return TQ.standard([ (sa,fa,c) for sa,fa,va in a ])

   @staticmethod
   def fillGaps(a,s,f,const=inf):
      c = []; fo = s
      for (sa,fa,va) in a:
         if fo<sa: c.append((fo,sa,const))
         c.append((sa,fa,va)); fo = fa
      if fo<f: c.append((fo,f,const))
      return(c)

   @staticmethod
   def changeTime(a,p):
      i = 0; t = p[i]; v = 0; b = []; cut = False
      for sc,fc,vc in a:
         while t <= sc:
            if v > 0: 
               if not cut: b.append((i,i+1,v)); v = 0
            i = i+1; t = p[i]
         if t < fc:
            v = v + (t-sc)*vc; b.append((i,i+1,v))
            v = (fc-t)*vc; cut = True
         else: v = v + (fc-sc)*vc; cut = False
      if v > 0: b.append((i,i+1,v))
      return(b)

   @staticmethod
   def complement(a,s,f):
      c = []; fo = s
      for (sa,fa,va) in a:
         if fo<sa: c.append((fo,sa,1))
         fo = fa
      if fo<f: c.append((fo,f,1))
      return(c)

   @staticmethod
   def invert(a,vZero=0):
      return [ (sa,fa,1/va) if va!=0 else (sa,fa,vZero) for (sa,fa,va) in a ]

   @staticmethod
   def minus(a):
      return [(sa,fa,-va) for (sa,fa,va) in a]

   @staticmethod
   def prodConst(a,c):
      return [ (sa,fa,va*c) for sa,fa,va in a ]

   @staticmethod
   def cutGT(a,c):
      return [(sa,fa,va) for (sa,fa,va) in a if va > c]

   @staticmethod
   def cutGE(a,c):
      return [(sa,fa,va) for (sa,fa,va) in a if va >= c]

   @staticmethod
   def height(a):
      s = -TQ.inf
      for (sa,fa,va) in a: s = max(s,va)
      return(s)

   @staticmethod
   def total(a):
      return 0 if a == [] else sum([(fa-sa)*va for (sa,fa,va) in a])

   @staticmethod
   def sum(a,b):
      if len(a) == 0: return(b)
      if len(b) == 0: return(a)
      c = []; A = a.__iter__(); B = b.__iter__()
      (sa,fa,va) = TQ.get(A); (sb,fb,vb) = TQ.get(B)
      while (sa<TQ.inf) or (sb<TQ.inf):
         if sa < sb:
            sc = sa; vc = va
            if sb < fa: fc = sb; sa = sb
            else: fc = fa; (sa,fa,va) = TQ.get(A)
         elif sa == sb:
            sc = sa; fc = min(fa,fb); vc = TQ.sAdd(va,vb)
            sa = sb = fc; fA = fa
            if fA <= fb: (sa,fa,va) = TQ.get(A)
            if fb <= fA: (sb,fb,vb) = TQ.get(B)
         else:
            sc = sb; vc = vb
            if sa < fb: fc = sa; sb = sa
            else: fc = fb; (sb,fb,vb) = TQ.get(B)
         c.append((sc,fc,vc))
      return(TQ.standard(c))

   @staticmethod
   def prod(a,b):
      if len(a)*len(b) == 0: return([])
      c = []; A = a.__iter__(); B = b.__iter__()
      (sa,fa,va) = TQ.get(A); (sb,fb,vb) = TQ.get(B)
      while (sa<TQ.inf) or (sb<TQ.inf):
         if fa <= sb: (sa,fa,va) = TQ.get(A)
         elif fb <= sa: (sb,fb,vb) = TQ.get(B)
         else:
            sc = max(sa,sb); fc = min(fa,fb); vc = TQ.sMul(va,vb)
            c.append((sc,fc,vc))
            if fc == fa: (sa,fa,va) = TQ.get(A)
            if fc == fb: (sb,fb,vb) = TQ.get(B)
      return(TQ.standard(c))

   @staticmethod
   def proportion(a,b):
      if len(a) == 0: return([])
      c = []; A = a.__iter__(); B = b.__iter__()
      (sa,fa,va) = TQ.get(A); (sb,fb,vb) = TQ.get(B)
      while (sa<TQ.inf) or (sb<TQ.inf):
         if fa <= sb: (sa,fa,va) = TQ.get(A)
         elif fb <= sa: (sb,fb,vb) = TQ.get(B)
         else:
            sc = max(sa,sb); fc = min(fa,fb); vc = va/vb
            c.append((sc,fc,vc))
            if fc == fa: (sa,fa,va) = TQ.get(A)
            if fc == fb: (sb,fb,vb) = TQ.get(B)
      return(TQ.standard(c))

   @staticmethod
   def extract(a,b):
   # extract a from b
      if len(a)*len(b) == 0: return([])
      c = []; A = a.__iter__(); B = b.__iter__()
      (sa,fa,va) = TQ.get(A); (sb,fb,vb) = TQ.get(B)
      while (sa<TQ.inf) or (sb<TQ.inf):
         if fa <= sb: (sa,fa,va) = TQ.get(A)
         elif fb <= sa: (sb,fb,vb) = TQ.get(B)
         else:
            sc = max(sa,sb); fc = min(fa,fb); vc = vb
            c.append((sc,fc,vc))
            if fc == fa: (sa,fa,va) = TQ.get(A)
            if fc == fb: (sb,fb,vb) = TQ.get(B)
      return(TQ.standard(c))

   @staticmethod
   def PFcheck(a,b):
      if len(a) == 0: return(a)
      if len(b) == 0: return(a)
      c = []; A = a.__iter__(); B = b.__iter__()
      (sa,fa,va) = TQ.get(A); (sb,fb,vb) = TQ.get(B)
      while (sa<TQ.inf) or (sb<TQ.inf):
         if fa <= sb: (sa,fa,va) = TQ.get(A)
         elif fb <= sa: (sb,fb,vb) = TQ.get(B)
         else:
            sc = max(sa,sb); fc = min(fa,fb)
            if vb == va: c.append((sc,fc,va))
            if fc == fa: (sa,fa,va) = TQ.get(A)
            if fc == fb: (sb,fb,vb) = TQ.get(B)
      return(TQ.standard(c))

   @staticmethod
   def PartMaxAdd(a,g,b=None):
      if b==None: b = []
      if len(a)==0: return(b)
      if len(b)==0:
         p = []
         for (sa,fa,va) in a: p.append((sa,fa,(va,g)))
         return(p)
      p = []; A = a.__iter__(); B = b.__iter__()
      (sa,fa,va) = TQ.get(A); (sb,fb,eb) = TQ.get(B)
      if type(eb) is tuple: (vb,gb) = eb
      while (sa<TQ.inf) or (sb<TQ.inf):
         if fa <= sb:
            p.append((sa,fa,(va,g))); (sa,fa,va) = TQ.get(A)
         elif fb <= sa:
            p.append((sb,fb,(vb,gb))); (sb,fb,eb) = TQ.get(B)
            if type(eb) is tuple: (vb,gb) = eb
         else:
            sp = max(sa,sb); fp = min(fa,fb)
            if vb > va: vp = vb; gp = gb
            else: vp = va; gp = g
            p.append((sp,fp,(vp,gp)))
            if fp == fa: (sa,fa,va) = TQ.get(A)
            if fp == fb:
               (sb,fb,eb) = TQ.get(B);
               if type(eb) is tuple: (vb,gb) = eb
      return(TQ.standard(p))

   @staticmethod
   def PartMaxVec(v,clu):
      p = []
      for i in clu: p = TQ.PartMaxAdd(v[i],i+1,p)
      return(p)

   @staticmethod
   def VPart2Part(p):
      q = []
      for (sp,fp,ep) in p:
         (vp,gp) = ep
         q.append((sp,fp,gp))
      return(TQ.standard(q))

   @staticmethod
   def renumPart(p):
      C = {}; q = []
      for a in p:
         r = []
         for (sa,fa,ca) in a:
            if not(ca in C): C[ca] = 1+len(C)
            r.append((sa,fa,C[ca]))
         q.append(r)
      return(q)

   ''' temporal dictionaries operations '''
   @staticmethod
   def TQdictCut(D,level):
      CC = { u: TQ.cutGE(core,level) for u, core in D.items() }
      return { u: core for u, core in CC.items() if core!=[] }

   ''' temporal vectors operations '''
   @staticmethod
   def VecConst(n,const=sE):
      if not(TQ.isTq(const)): raise TQ.TQerror("VecConst: const must be a temporal quantity")
      return([ const for i in range(n) ])

   @staticmethod
   def VecList(a,names=None,all=True,skip=[]):
      for (i,e) in enumerate(a):
         if all or (e!=skip):
            if names==None: print (i+1, ":", e)
            else: print (i, names[i], "=" ,e)

   @staticmethod
   def VecSum(a,b):
      na = len(a); nb = len(b)
      if na!=nb: raise TQ.TQerror("VecSum: vectors of same size required")
      c = [ [] for i in range(na) ]
      for i in range(na): c[i] = TQ.sum(a[i],b[i])
      return(c)

   @staticmethod
   def VecProd(a,b):
      na = len(a); nb = len(b)
      if na!=nb: TQ.TQerror("VecProd: vectors of same size required")
      c = [ [] for i in range(na) ]
      for i in range(na): c[i] = TQ.prod(a[i],b[i])
      return(c)

   @staticmethod
   def VecInv(a):
      na = len(a)
      c = [ [] for i in range(na) ]
      for i in range(na): c[i] = TQ.invert(a[i])
      return(c)

   @staticmethod
   def VecTotal(V):
      n = len(V); T = [0]*n
      for u in range(n): T[u] = TQ.total(V[u])
      Q = sorted(enumerate(T),key=operator.itemgetter(1),reverse=True)
      return(Q)

   ''' temporal matrices operations '''
   @staticmethod
   def MatSummary(A):
      nr = len(A); nc = len(A[0])
      vMin = tMin = TQ.inf; vMax = tMax = -TQ.inf
      for u in range(nr):
         for v in range(nc):
            (ti,ta,vi,va) = TQ.TqSummary(A[u][v])
            tMin = min(tMin,ti); tMax = max(tMax,ta)
            vMin = min(vMin,vi); vMax = max(vMax,va)
      return((tMin,tMax,vMin,vMax))

   @staticmethod
   def MatList(M,names=None,all=True,skip=[]):
      nr = len(M); nc = len(M[0])
      for i in range(nr):
         for j in range(nc):
            if all or (M[i][j]!=skip):
               if names==None: print ('(',i+1,',',j+1,') =',M[i][j])
               else:
                  print ('(',names[i],',',names[j],') =\n  ',M[i][j])

   @staticmethod
   def MatClosure(R,strict=False):
      nr = len(R); nc = len(R[0])
      if nr!=nc: raise TQ.TQerror("MatClosure: square matrix required")
      C = deepcopy(R)
      for k in range(nr):
         for u in range(nr):
            for v in range(nr):
               C[u][v] = TQ.sum(C[u][v], TQ.prod(C[u][k],C[k][v]))
         if not strict: C[k][k] = TQ.sum(TQ.sE,C[k][k])
      return(C)

   @staticmethod
   def MatSetVal(A,val):
      nr = len(A); nc = len(A[0])
      B = [[ [] for i in range(nc)] for j in range(nr)]
      for u in range(nr):
         for v in range(nc):
            b = []
            for (sa,fa,va) in A[u][v]: b.append((sa,fa,val))
            B[u][v] = TQ.standard(b)
      return(B)

   @staticmethod
   def MatBin(A): return(TQ.MatSetVal(A,1))

   @staticmethod
   def MatConst(A,const):
      if not(TQ.isTq(const)): raise TQ.TQerror("MatConst: const must be a temporal quantity")
      nr = len(A); nc = len(A[0])
      B = [[ [] for i in range(nc)] for j in range(nr)]
      for u in range(nr):
         for v in range(nc):
            B[u][v] = TQ.prod(const,A[u][v])
      return(B)

   @staticmethod
   def MatSum(A,B):
      nra = len(A); nca = len(A[0]); nrb = len(B); ncb = len(B[0]);
      if (nra!=nrb) or (nca!=ncb):
         raise TQ.TQerror("MatSum: only matrices of the same size can be summed")
      C = [[ [] for i in range(nca)] for j in range(nra)]
      for u in range(nra):
         for v in range(nca):
            C[u][v] = TQ.sum(A[u][v],B[u][v])
      return(C)

   @staticmethod
   def MatInter(A,B):
      nra = len(A); nca = len(A[0]); nrb = len(B); ncb = len(B[0]);
      if (nra!=nrb) or (nca!=ncb):
         raise TQ.TQerror("MatInter: only matrices of the same size can be intersected")
      C = [[ [] for i in range(nca)] for j in range(nra)]
      old = TQ.semiring; TQ.reach();
      for u in range(nra):
         for v in range(nca):
            C[u][v] = TQ.prod(A[u][v],B[u][v])
      old()
      return(C)

   @staticmethod
   def MatSetDiag(A,const):
      nr = len(A); nc = len(A[0])
      if not(TQ.isTq(const)): raise TQ.TQerror("MatSetDiag: const must be a temporal quantity")
      if nr!=nc: raise TQ.TQerror("MatSetDiag: square matrix required")
      B = deepcopy(A)
      for u in range(nr): B[u][u] = const
      return(B)

   @staticmethod
   def MatPower(A,k):
      n = len(A); nc = len(A[0])
      if n!=nc: raise TQ.TQerror("MatPower: square matrix required")
      P = [[[] for u in range(n)] for v in range(n)]
      P = TQ.MatSetDiag(P,TQ.sE)
      if k > 0:
         i = k; Q = deepcopy(A)
         while True:
            if i%2 == 1: P = TQ.MatProd(P,Q)
            i = i // 2
            if i == 0: break
            Q = TQ.MatProd(Q,Q)
      return(P)

   @staticmethod
   def MatVecLeft(A,x):
      nr = len(A); nc = len(A[0]); nv = len(x)
      if nv!=nr:
         raise TQ.TQerror("MatVecLeft: vector length should be equal to number of rows")
      y = [ [] for i in range(nc)]
      for v in range(nc):
         s = []
         for u in range(nr):
            s = TQ.sum(s,TQ.prod(x[u],A[u][v]))
         y[v] = s
      return(y)

   @staticmethod
   def MatVecRight(A,x):
      nr = len(A); nc = len(A[0]); nv = len(x)
      if nv!=nc:
         raise TQ.TQerror("MatVecRight: vector length should be equal to number of columns")
      y = [ [] for i in range(nr)]
      for u in range(nr):
         s = []
         for v in range(nc):
            s = TQ.sum(s,TQ.prod(A[u][v],x[v]))
         y[u] = s
      return(y)

   @staticmethod
   def MatVecRow(A,x):
      nr = len(A); nc = len(A[0]); nv = len(x)
      if nv!=nr:
         raise TQ.TQerror("MatVecRow: vector length should be equal to number of rows")
      B = [[ [] for i in range(nc)] for j in range(nr)]
      for u in range(nr):
         for v in range(nc):
            B[u][v] = TQ.prod(x[u],A[u][v])
      return(B)

   @staticmethod
   def MatVecCol(A,x):
      nr = len(A); nc = len(A[0]); nv = len(x)
      if nv!=nr:
         raise TQ.TQerror("MatVecCol: vector length should be equal to number of columns")
      B = [[ [] for i in range(nc)] for j in range(nr)]
      for u in range(nr):
         for v in range(nc):
            B[u][v] = TQ.prod(A[u][v],x[v])
      return(B)

   @staticmethod
   def MatProd(A,B):
      nra = len(A); nca = len(A[0]); nrb = len(B); ncb = len(B[0]);
      if nca!=nrb:
         raise TQ.TQerror("MatProd: only compatible matrices can be multiplied")
      C = [[ [] for i in range(ncb)] for j in range(nra)]
      for u in range(nra):
         for v in range(ncb):
            s = []
            for t in range(nca):
               s = TQ.sum(s,TQ.prod(A[u][t],B[t][v]))
            C[u][v] = s
      return(C)

   @staticmethod
   def MatProdDiag(A,B):
      nra = len(A); nca = len(A[0]); nrb = len(B); ncb = len(B[0]);
      if (nca!=nrb) or (nra!=ncb):
         raise TQ.TQerror("MatProdDiag: only compatible matrices can be multiplied")
      x = [ [] for u in range(nra)]
      for u in range(nra):
         s = []
         for v in range(nca):
            s = TQ.sum(s,TQ.prod(A[u][v],B[v][u]))
         x[u] = s
      return(x)

   @staticmethod
   def MatTrans(A):
      nr = len(A); nc = len(A[0])
      B = [[[] for v in range(nr)] for u in range(nc)]
      for u in range(nr):
         for v in range(nc): B[v][u] = A[u][v]
      return(B)

   @staticmethod
   def MatSym(A):
      nr = len(A); nc = len(A[0])
      if nr!=nc: raise TQ.TQerror("MatSym: square matrix required")
      B = deepcopy(A)
      for u in range(nr-1):
         for v in range(u+1,nr):
            B[u][v] = TQ.sum(A[u][v],A[v][u])
            B[v][u] = B[u][v]
      return(B)

   @staticmethod
   def MatExtract(T,A):
      nt = len(T); nr = len(A); nc = len(A[0])
      if nr!=nt:
         raise TQ.TQerror("MatExtract: partition and matrix of same size required")
      B = [[[] for v in range(nc)] for u in range(nr)]
      for u in range(nr):
         for v in range(nc):
            B[u][v] = TQ.extract(T[u],A[u][v])
      return(B)

   ''' temporal networks methods '''
   @staticmethod
   def inDeg(A):
      return(TQ.MatVecLeft(TQ.MatBin(A),TQ.VecConst(len(A))))

   @staticmethod
   def outDeg(A):
      return(TQ.MatVecRight(TQ.MatBin(A),TQ.VecConst(len(A[0]))))

   @staticmethod
   def ActInt(C,Rows,Cols):
   # extend for Rows and Cols are temporal partitions
      s = []
      for u in Rows:
         for v in Cols: s = TQ.sum(s,C[v][u])
      return(s)

   @staticmethod
   def clusCoef(A,type=1):
   # type = 1 - standard clustering coefficient
   # type = 2 - corrected clustering coefficient / temporal degMax
   # type = 3 - corrected clustering coefficient / overall degMax
      nr = len(A); nc = len(A[0])
      if nr!=nc: raise TQ.TQerror("clusCoef: square matrix required")
      if (type<1)or(type>3): raise TQ.TQerror("clusCoef: unsupported type")
      Ab = TQ.MatSetDiag(TQ.MatBin(A),TQ.sN)
      S = TQ.MatBin(TQ.MatSym(Ab))
      deg = TQ.MatVecRight(S,TQ.VecConst(nr))
      if type == 1:
         fac = TQ.VecProd(deg,TQ.VecSum(deg,TQ.VecConst(nr,const=[(1,TQ.inf,-1)])))
      else:
         old = TQ.semiring(); TQ.maxmin(); delta = []
         for d in deg: delta = TQ.sum(delta,d)
         if type == 3:
            Delta = max([v for (s,f,v) in delta])
            delta = [(1,TQ.inf,Delta)]
         TQ.combinatorial(); fac = []
         degm = TQ.VecSum(deg,TQ.VecConst(nr,const=[(1,TQ.inf,-1)]))
         for d in degm: fac.append(TQ.prod(delta,d))
      tri = TQ.MatProdDiag(TQ.MatProd(S,Ab),S)
      cc = TQ.VecProd(TQ.VecInv(fac),tri)
      old()
      return(cc)

   @staticmethod
   def pathFinder(W,r=1,q=inf,closure=False):
   # W is a dissimilarity matrix !!!
      nr = len(W); nc = len(W[0])
      if nr!=nc: raise TQ.TQerror("pathFinder: square matrix required")
      PF = [[[] for v in range(nr)] for u in range(nr)]
      old = TQ.semiring; rold = TQ.rPF; TQ.rPF = r; TQ.PFsemi()
      Z = TQ.MatClosure(W) if q>nr else TQ.MatPower(TQ.MatSetDiag(W,TQ.sE),q)
      if closure:
         print("Matrix/Closure"); TQ.MatList(Z,all=False)
      TQ.rPF = rold; old()
      for u in range(nr):
         for v in range(nr):
   #        pf = []
   ##       for t in T:
   ##          if Z[u][v][t] == W[u][v][t]: pf.append(W[u][v][t])
           PF[u][v] = TQ.PFcheck(W[u][v],Z[u][v])
      return(PF)

   @staticmethod
   def eqMat2Part(E):
      old = TQ.semiring; TQ.combinatorial()
      v = [ [(1,TQ.inf,i+1)] for i in range(len(E))]
      random.shuffle(v)
      p = TQ.MatVecRight(E,v)
      old()
      return(TQ.renumPart(p))

   @staticmethod
   def weakConnMat(A):
      old = TQ.semiring; TQ.reach()
      W = TQ.MatClosure(TQ.MatSym(TQ.MatBin(A)),strict=True)
      old()
      return(W)

   @staticmethod
   def strongConnMat(A):
      old = TQ.semiring; TQ.reach()
      R = TQ.MatClosure(TQ.MatBin(A),strict=True)
      S = TQ.MatInter(R,TQ.MatTrans(R))
      old()
      return(S)

   @staticmethod
   def attraction(A,act=None):
      n = len(A); att = [[]]*n
      old = TQ.semiring; TQ.combinatorial()
      if act==None:
         act = [[]]*n;
         for u in range(n):
            a = []
            for v in range(n):
               a = TQ.sum(a,A[u][v])
            act[u] = a
      for u in range(n):
         s = []
         for v in range(n):
            if (A[u][v]!=[]) and (u!=v):
               s = TQ.sum(s,TQ.proportion(A[v][u],act[v]))
         att[u] = TQ.proportion(s,[(1,TQ.inf,n-1)])
# alternative - replace n-1 with indeg
      old()
      return(att)

   @staticmethod
   def closeness(A,type=2):
# type: 1 - output, 2 - all, 3 - input
      (s,f,x,y) = TQ.MatSummary(A)
      old = TQ.semiring; TQ.path(); n = len(A)
      C = TQ.MatClosure(A,strict=True)
#      print("Matrix Dist");  TQ.MatList(C,all=False)
      cl = [ [] for i in range(n) ]
      TQ.combinatorial()
      fac = [(1,TQ.inf,(n-1)*(2 if type==2 else 1))]
      for v in range(n):
         d = []
         for u in range(n):
            if u!=v:
               if type<3: d = TQ.sum(d,TQ.fillGaps(C[v][u],s,f))
               if type>1: d = TQ.sum(d,TQ.fillGaps(C[u][v],s,f))
         cl[v] = TQ.prod(fac,TQ.invert(d))
      old()
      return(cl)

   @staticmethod
   def between(uv,vw,uw):
      if len(uv)==0: return([])
      if len(vw)==0: return([])
      if len(uw)==0: return([])
      r = []; A = uv.__iter__(); B = vw.__iter__(); C = uw.__iter__()
      (sa,fa,va) = TQ.get(A); (sb,fb,vb) = TQ.get(B); (sc,fc,vc) = TQ.get(C)
      if type(va) is tuple: (da,ca) = va
      if type(vb) is tuple: (db,cb) = vb
      if type(vc) is tuple: (dc,cc) = vc
      while (sa<TQ.inf) or (sb<TQ.inf) or (sc<TQ.inf):
         sr = max(sa,sb,sc); fr = min(fa,fb,fc)
         if fa <= sr:
            (sa,fa,va) = TQ.get(A)
            if type(va) is tuple: (da,ca) = va
         elif fb <= sr:
            (sb,fb,vb) = TQ.get(B)
            if type(vb) is tuple: (db,cb) = vb
         elif fc <= sr:
            (sc,fc,vc) = TQ.get(C)
            if type(vc) is tuple: (dc,cc) = vc
         else:
            if da+db == dc: r.append((sr,fr,ca*cb/cc))
            if fr == fa:
               (sa,fa,va) = TQ.get(A)
               if type(va) is tuple: (da,ca) = va
            if fr == fb:
               (sb,fb,vb) = TQ.get(B)
               if type(vb) is tuple: (db,cb) = vb
            if fr == fc:
               (sc,fc,vc) = TQ.get(C)
               if type(vc) is tuple: (dc,cc) = vc
      return(TQ.standard(r))

   @staticmethod
   def betweenness(A):
      G = TQ.MatSetVal(A,(1,1))
      old = TQ.semiring; TQ.geodetic(); n = len(G)
      C = TQ.MatClosure(G,strict=True)
#      print("Matrix C1");  TQ.MatList(C,all=False)
      bw = [ [] for i in range(n) ]
      TQ.combinatorial()
      fac = [(1,TQ.inf,1/(n-1)/(n-2))]
      for v in range(n):
         b = []
         for u in range(n):
            for w in range(n):
               if (C[u][w]!=[]) and (u!=w) and (u!=v) and (v!=w):
                  b = TQ.sum(b,TQ.between(C[u][v],C[v][w],C[u][w]))
         bw[v] = TQ.prod(b,fac)
      old()
      return(bw)

   ''' saving and loading operations '''
   @staticmethod
   def VecSave(a,file="test.vet"):
      nr = len(a)
      with open(file, "w") as vet:
         vet.write("*nodes "+str(nr)+"\n")
         for i in range(nr):
            vet.write(str(a[i])+'\n')
         vet.close()

   @staticmethod
   def Ianus2Mat(tenFile):
      net = open(tenFile,'r',encoding='utf-8')
#      net = open(tenFile,encoding='utf-8-sig','r')
      line = net.readline()
      if not line: raise TQ.TQerror("readTen: empty file")
      if line[0:6]!="%Ianus": raise TQ.TQerror("readTen: not an Ianus file")
      control = ""; k = 1; typ = ['simple', 'directed']
      while True:
         line = net.readline()
         if not line: break
         k = k+1
         if line[0]=='*':
            control = line[1:4].upper()
            if control=="NET":
               L = re.split(':',line.strip(),1)
               tit = L[1]
               if tit[0]=='"': tit = tit[1:-1]
            elif control=="TYP":
               L = re.split(':',line.strip(),1)
               typ = eval(L[1])
            elif control=="NOD":
               L = re.split('\s+',line.strip())
               nc = eval(L[1])
               names = ["v"+str(i+1) for i in range(nc)]
               Tnode = [[(tMin,tMax,1)]]*nc
               twoMode = len(L)>2
               if twoMode:
                  nr = eval(L[2]); nc = nc-nr
               else:
                  nr = nc
            elif control=="TIM":
               L = re.split('\s+',line.strip())
               tMin = eval(L[1]); tMax = eval(L[2])
               Tlabs = []
            elif control=="MET":
               meta =""
            elif control=="ARC":
               W = [[ [] for i in range(nc)] for j in range(nr)]
            elif control=="EDG":
               W = [[ [] for i in range(nc)] for j in range(nr)]
            else: raise TQ.TQerror("readTen: unsupported keyword: "+control)
         else:
            if control=="NOD":
               L = re.split('\s+',line.strip(),1)
               i = eval(L[0]); bf = L[1]
               if bf[0]=='"':
                  j = bf.find('"',1); name = bf[1:j]
                  bf = bf[j+1:].lstrip()
               else:
                  L = re.split('\s+',bf,1)
                  name = L[0]; bf = L[1]
               names[i-1] = name
               if len(bf)>0:
                  Tnode[i-1] = eval(bf)
            if control=="ARC":
               L = re.split('\s+',line.strip(),2)
               u = eval(L[0])-1; v = eval(L[1])-1; w = eval(L[2])
               if twoMode:
                  W[u][v-nr] = TQ.sum(W[u][v-nr],w)
               else:
                  W[u][v] = TQ.sum(W[u][v],w)
            if control=="EDG":
               L = re.split('\s+',line.strip(),2)
               u = eval(L[0])-1; v = eval(L[1])-1; w = eval(L[2])
               if twoMode:
                  W[u][v-nr] = TQ.sum(W[u][v-nr],w)
               else:
                  W[u][v] = TQ.sum(W[u][v],w)
            elif control=="TIM":
               L = re.split('\s+',line.strip())
               i = eval(L[0]); lab = L[1]
               if lab[0]=='"': lab = lab[1:-1]
               Tlabs += [ (i, lab) ]
            elif control=="MET":
               meta += line
      net.close()
      if twoMode:
         if not('twomode' in typ): typ += ['twomode']
         if 'onemode' in typ: typ.remove('onemode')
      else:
         if 'twomode' in typ: typ.remove('twomode')
         if not('onemode' in typ): typ += ['onemode']
      dim = (nr,nc,tMin,tMax)
      rez = {}
      rez['dim'] = dim
      rez['met'] = meta
      rez['typ'] = typ
      rez['nam'] = names
      rez['mat'] = W
      rez['tit'] = tit
      rez['tin'] = Tnode
      rez['til'] = Tlabs
      return(rez)

   @staticmethod
   def MatSave(N,file="test.ten"):
      K = N.keys()
      if not('dim' in K): raise TQ.TQerror("MatSave: missing DIM")
      if not('mat' in K): raise TQ.TQerror("MatSave: missing MAT")
      (nr, nc, minT, maxT) = N['dim']; M = N['mat']
      with open(file, "w",encoding='utf-8') as net:
         net.write("%Ianus\n*metadata\n")
         if 'met' in K: net.write(N['met'])
         net.write("\nre\ndt "+datetime.datetime.now().ctime()+
            "\nti saved from Ianus\ner\n")
         tit = "no title"
         if 'tit' in K: tit = N['tit']
         net.write('*network:"'+tit+'"\n')
         typ = ["simple", "onemode", "directed"]
         if 'typ' in K: typ = N['typ']
         net.write('*type:'+str(typ)+'\n')
         twoMode = (nr != nc) or ('twomode' in typ)
         net.write('*timescale '+str(minT)+' '+str(maxT)+'\n')
         if 'til' in K:
            for (i,l) in N['til']:
               net.write(str(i)+' "'+l+'"\n')
         names = None; times = None
         if 'nam' in K: names = N['nam']
         if 'tin' in K: times = N['tin']
         if twoMode:
            net.write("*nodes "+str(nr+nc)+" "+str(nr)+"\n")
            if names != None:
               for i in range(nr+nc):
                  if times==None:
                     net.write(str(i+1)+' "'+names[i]+'"\n')
                  else:
                     net.write(str(i+1)+' "'+names[i]+'" '+str(times[i])+'\n')
            net.write("*arcs\n")
            for i in range(nr):
               for j in range(nc):
                  if len(M[i][j])>0:
                     net.write(str(i+1)+' '+str(nr+j+1)+' '+str(M[i][j])+'\n')
         else:
            net.write("*nodes "+str(nr)+"\n")
            if names != None:
               for i in range(nr):
                  if times==None:
                     net.write(str(i+1)+' "'+names[i]+'"\n')
                  else:
                     net.write(str(i+1)+' "'+names[i]+'" '+str(times[i])+'\n')
            net.write("*arcs\n")
            for i in range(nr):
               for j in range(nc):
                  if len(M[i][j])>0:
                     net.write(str(i+1)+' '+str(j+1)+' '+str(M[i][j])+'\n')
         # raise TQ.TQerror("MatSave: not implemented yet")
         net.close()

   @staticmethod
   def Ianus2netJSON(N,fileJSON="test.json",indent=None):
      K = N.keys()
      if not('dim' in K): raise TQ.TQerror("Ianus2netJSON: missing DIM")
      if not('mat' in K): raise TQ.TQerror("Ianus2netJSON: missing MAT")
      (nr, nc, minT, maxT) = N['dim']; M = N['mat']; org = 1
      info = {}; types = N.get('typ',[])
      info['network'] = N.get('tit',"test");
      info['org'] = org; info['temporal'] = True
      info['simple'] = 'simple' in types
      info['directed'] = 'directed' in types
      info['multirel'] = 'multirel' in types
      info['mode'] = 1 if nr==nc else 2
      info['meta'] = [ N.get('met',{}) ]
      info['meta'].append({"date": datetime.datetime.now().ctime(),\
         "title": "saved from Ianus to netJSON" })
      info['nNodes'] = nr; time = { "Tmin": minT, "Tmax": maxT }
#      if 'til' in K: Tlabs = { str(k): v for k,v in N['til'] }
#      else:
      Tlabs = { str(y):str(y) for y in range(minT,maxT+1)}
      time['Tlabs'] = Tlabs
      if nr==nc:
         nodes = []; names = N.get('nam', ['v'+str(v+org) for v in range(nr)])
         nodeAct = N.get('tin', [[(minT, maxT+1, 1)] for v in range(nr)])
         for v in range(nr):
            Node = { 'id': v+1, 'lab': names[v], 'tq': nodeAct[v] }
            nodes.append(Node)
         links = []; ltype = "arc"
         for u in range(nr):
            for v in range(nc):
               if M[u][v]!=[]:
                  Link = {"type": ltype, "n1": u+1, "n2": v+1, # "rel": link[3]
                     "tq": M[u][v] }
                  links.append(Link)
      else:
#         raise TQ.TQerror("Ianus2netJSON: two-mode not implemented yet")
         n = nr+nc; info['nNodes'] = n; info['dim'] = [nr, nc]
         nodes = []; names = N.get('nam', ['v'+str(v+org) for v in range(n)])
         nodeAct = N.get('tin', [[(minT, maxT+1, 1)] for v in range(n)])
         for v in range(nr):
            Node = { 'id': v+1, 'lab': names[v], 'mode':1, 'tq': nodeAct[v] }
            nodes.append(Node)
         for v in range(nr,n):
            Node = { 'id': v+1, 'lab': names[v], 'mode':2, 'tq': nodeAct[v] }
            nodes.append(Node)
         links = []; ltype = "arc"
         for u in range(nr):
            for v in range(nc):
               if M[u][v]!=[]:
                  Link = {"type": ltype, "n1": u+1, "n2": nr+v+1, # "rel": link[3]
                     "tq": M[u][v] }
                  links.append(Link)
      info['nArcs'] = len(links); info['nEdges'] = 0; info['time'] = time
      net = {"netJSON": "basic", "info": info, "nodes": nodes, "links": links}
#      js = open(info['network']+'.json','w')
      js = open(fileJSON,'w')
      json.dump(net, js, ensure_ascii=False, indent=indent)
      js.close()

   @staticmethod
   def Mat2Pajek(M,file="test.net",names=None,labels=False):
      nr = len(M); nc = len(M[0])
      with open(file, "w") as net:
         if nr == nc: # one-mode
            net.write("*vertices "+str(nr)+"\n")
            if names != None:
               for i in range(nr):
                  net.write(str(i+1)+' "'+names[i]+'"\n')
            net.write("*arcs\n")
            for i in range(nr):
               for j in range(nc):
                  if len(M[i][j])>0:
                     net.write(str(i+1)+' '+str(j+1))
                     if labels: net.write(' 1 l "'+str(M[i][j])+'"\n')
                     else: net.write('\n')
         else: raise TQ.TQerror("Mat2Pajek: not implemented yet")
         net.close()

   ''' additional unasigned operations '''
   @staticmethod
   def project(a,ind):
      c = []
      for (sa,fa,va) in a: c.append((sa,fa,va[ind]))
      return(TQ.standard(c))

   @staticmethod
   def filter(a,test):
      c = []
      for (sa,fa,va) in a:
         if test(va): c.append((sa,fa,va))
      return(c)

   @staticmethod
   def central(v,C):
      cc = c = []; n = len(C)
      for t in range(n):
         if t!=v:
            c = TQ.sum(c,TQ.invert(C[v][t]))
      for (s,f,v) in c: cc.append((s,f,v/(n-1)))
      return(cc)

   @staticmethod
   def minTime(A):
      old = TQ.semiring; TQ.reach(); n = len(A)
      T = [ [] for i in range(n) ]
      for u in range(n):
         d = []
         for v in range(n):
            if A[u][v]!=[]: d = TQ.sum(d,TQ.binary(A[u][v]))
         T[u] = d
      old()
      return(T)

   @staticmethod
   def MatBidir(A):
      nr = len(A); nc = len(A[0])
      if nr!=nc: raise TQ.TQerror("MatBidir: square matrix required")
      B = TQ.MatBin(A)
      C = [[[] for v in range(nr)] for u in range(nr)]
      old = TQ.semiring; TQ.reach()
      for u in range(nr-1):
         for v in range(u+1,nr):
            C[u][v] = TQ.prod(B[u][v],B[v][u])
            C[v][u] = C[u][v]
      old()
      return(C)

