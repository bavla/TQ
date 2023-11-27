# trajectoR
# Vladimir Batagelj, november 2023

traj2Pajek <- function(E,kMax,Net){
  I <- order(E$s,E$f)
  n <- length(I); k <- 0; r <- 0
  cn <- rep(0,kMax); cs <- rep("",kMax)
  N <- data.frame(u=cs,v=cs,s=cn,f=cn,rel=cs)
  cat("% traj2Pajek",date(),"\nevents",n,"\n")
  for(p in 1:(n-1)){i <- I[p]; tm <- E$f[i]
    cat("."); if(p%%50==0) {cat(p,k,date(),"\n")}; flush.console()
    for(q in (p+1):n){j <- I[q]; r <- r+1
      if(E$s[j]>tm) break
      if(E$S[i]==E$S[j]){
        fm <- min(E$f[i],E$f[j]); sM <- max(E$s[i],E$s[j]); T <- fm-sM
        if(T>0){k <- k+1; if(k>kMax) stop("kMax too small")
          N[k,] <- list(u=E$ID[i],v=E$ID[j],s=sM,f=fm,rel=E$S[i])}
      }
    }
  }
  cat("\n",date(),"\ndensity R =",2*r/n/(n-1),"  tests =",r,
    "\ndensity E =",2*k/n/(n-1),"  edges =",k,"\n"); flush.console()
  N <- N[1:k,]; sf <- as.matrix(N[,c("s","f")])
  uvrwt2net(N$u,N$v,w=N$f-N$s,r=N$rel,t=sf,directed=FALSE,Net=Net)
  cat("% finished",date(),"\n")
}
