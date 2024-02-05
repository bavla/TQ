# trajectoR / version 4
# Vladimir Batagelj, December 2023 / February 2024
# Network of time interval intersections / equal part_of_cv and institution_si
# different units; union of intervals

datum <- function(d,m,y) {
  dd <- ifelse(d!="",d,ifelse(is.na(m)||(m==""),paste("01/01/",y,sep=""),
    paste("01/",m,"/",y,sep="")))
  di <- as.integer(as.Date(dd,format="%d/%m/%Y",origin="1970-01-01"))
  if(is.na(di)) {OK <<- FALSE
    di <- as.integer(as.Date(dd,format="%m/%d/%Y",origin="1970-01-01"))}
  return(di)
}

raw2events <- function(D,Tday="",CSV="events.csv"){
  if(Tday=="") Tday <- format(Sys.time(),"%d/%m/%Y")
  csv <- file(CSV,"w",encoding="UTF-8")
  cat("R  ID  s  f  S  T\n",file=csv)
  n <- nrow(D); OK <- TRUE; k <- 0; td <- datum(Tday,0,0) 
  for(i in 1:n){add <- TRUE
    ID <- trimws(D$person_name[i]); rel <- trimws(tolower(D$part_of_cv[i]))
    test <- trimws(tolower(D$institution_si[i]))
    ds <- D$start_day[i]; ms <- D$start_month[i]; ys <- D$start_year[i]
    sd <- datum(ds,ms,ys)
    if(sd>=td) add <- FALSE
    if(!OK) {cat(i,":",ID,ds,ms,ys,rel,'*** wrong date\n')
      flush.console(); OK <- TRUE; add <- FALSE}
    de <- D$end_day[i]; me <- D$end_month[i]; ye <- D$end_year[i]
    ed <- if(ye==2100) datum(Tday,NA,0) else datum(de,me,ye)
    ed <- min(ed,td)
    if(!OK) {cat(i,":",ID,de,me,ye,rel,'*** wrong date\n')
      flush.console(); OK <- TRUE; add <- FALSE}
    if(add) cat(i,' "',ID,'" ',sd,' ',ed,' "',rel,'" "',test,'"\n',sep='',file=csv) else
      k <- k+1
  }
  close(csv)
  if(k>0) {cat(k," events skipped\n"); flush.console()}
}

intUnion <- function(N){
  k <- nrow(N); z <- rep(0,k); L <- data.frame(s=z,f=z)
  i <- 1; r <- 0; s <- N$s[i]; f <- N$f[i]
  while(i<k){i <- i+1
    if(N$s[i]<=(f+1)){ if(N$f[i]>f) f <- N$f[i]
    } else { r <- r+1; L[r,] <- list(s=s,f=f); s <- N$s[i]; f <- N$f[i] }
  }
  r <- r+1; L[r,] <- list(s=s,f=f)
  return(L[1:r,])
}

traj2Pajek <- function(E,kMax,Net){
  I <- order(E$s,E$f)
  n <- length(I); k <- 0; r <- 0; sent <- "?????"
  cn <- rep(0,kMax); cs <- rep("",kMax)
  N <- data.frame(u=cs,v=cs,s=cn,f=cn,rel=cs,test=cs)
  cat("% traj2Pajek",date(),"\nevents",n,"\n")
  for(p in 1:(n-1)){i <- I[p]; tm <- E$f[i]
    cat("."); if(p%%50==0) {cat(p,k,date(),"\n")}; flush.console()
    for(q in (p+1):n){j <- I[q]; r <- r+1
      if(E$s[j]>tm) break
      uID <- E$ID[i]; vID <- E$ID[j]
      if(uID!=vID) {
        if(uID > vID) {vID <- E$ID[i]; uID <- E$ID[j]}
        if(E$S[i]==E$S[j]) if(E$T[i]==E$T[j]) {
          fm <- min(E$f[i],E$f[j]); sM <- max(E$s[i],E$s[j]); T <- fm-sM+1
          if(T>0){k <- k+1; if(k>kMax) stop("kMax too small")
            N[k,] <- list(u=uID,v=vID,s=sM,f=fm,rel=E$S[i],test=E$T[i])}
        }
      }
    }
  }
  cat("\n",date(),"\ndensity R =",2*r/n/(n-1),"  tests =",r,
    "\ndensity E =",2*k/n/(n-1),"  edges =",k,"\n"); flush.console()
  N <- N[1:k,]; I <- order(N$u,N$v,N$rel,N$test,N$s,N$f); N <- N[I,]
  N[k+1,] <- list(u=sent,v=sent,s=0,f=0,rel="",test="")
  cn <- rep(0,k); cs <- rep("",k); r <- 0; q <- 0; k <- k+1
  M <- data.frame(u=cs,v=cs,s=cn,f=cn,rel=cs,test=cs)
  Q <- data.frame(s=cn,f=cn); i <- 1; H <- N[i,c("u","v","rel","test")]
  for(p in 1:k){
    G <- N[p,c("u","v","rel","test")]
    if(all(H==G)) {q <- q+1; Q[q,] <- list(s=N$s[p],f=N$f[p])} else {
      L <- intUnion(Q[1:q,])
      for(h in 1:nrow(L)){ r <- r+1
        M[r,] <- list(u=H$u,v=H$v,s=L$s[h],f=L$f[h],rel=H$rel,test=H$test)
      }
      if(p==k) break
      H <- G; i <- p; q <- 0
    }
  }
  M <- M[1:r,]; sf <- as.matrix(M[,c("s","f")])
  cat("density E' =",2*r/n/(n-1),"  edges' =",r,"\n"); flush.console()
  uvrwt2net(M$u,M$v,w=M$f-M$s+1,r=M$rel,t=sf,directed=FALSE,Net=Net)
  cat("% finished",date(),"\n")
}

explain <- function(E,p1,p2,kMax=500){
  i12 <- c(which(E$ID == p1),which(E$ID == p2))
  C <- E[i12,]; n <- nrow(C); I <- order(C$s,C$f); k <- 0; r <- 0
  cn <- rep(0,kMax); cs <- rep("",kMax); sent <- "?????"
  N <- data.frame(u=cs,v=cs,s=cn,f=cn,d=cn,rel=cs,test=cs,Ru=cs,Rv=cs)
  for(p in 1:(n-1)){i <- I[p]; tm <- C$f[i]
    for(q in (p+1):n){j <- I[q]; r <- r+1
      if(C$s[j]>tm) break
      uID <- C$ID[i]; vID <- C$ID[j]
      if(uID!=vID) 
      if(C$S[i]==C$S[j]) if(C$T[i]==C$T[j]) {
        fm <- min(C$f[i],C$f[j]); sM <- max(C$s[i],C$s[j]); T <- fm-sM+1
        if(T>0){k <- k+1; if(k>kMax) stop("kMax too small")
          if(uID > vID) {vID <- C$ID[i]; uID <- C$ID[j]}
          N[k,] <- list(u=uID,v=vID,s=sM,f=fm,d=T,
            rel=C$S[i],test=C$T[i],Ru=C$R[i],Rv=C$R[j])}
      }
    }
  }
  tA <- sum(N$d); tS <- sum(N[N$u!=N$v,"d"])
  cat(p1,":",p2,"\n",date(),"\n density R =",2*r/n/(n-1),
    "  tests =",r,"\n density E =",2*k/n/(n-1),"  edges =",k,
    "\n time all =",tA,"  time strict =",tS,"\n")
  N <- N[1:k,]; I <- order(N$u,N$v,N$rel,N$test,N$s,N$f); N <- N[I,]
  k <- k+1; N[k,] <- list(u=sent,v=sent,s=0,f=0,rel="",test="")
  cn <- rep(0,k); cs <- rep("",k); r <- 0; q <- 0
  M <- data.frame(u=cs,v=cs,s=cn,f=cn,d=cn,rel=cs,test=cs)
  Q <- data.frame(s=cn,f=cn); i <- 1; H <- N[i,c("u","v","rel","test")]
  for(p in 1:k){
    G <- N[p,c("u","v","rel","test")]
    if(all(H==G)) {q <- q+1; Q[q,] <- list(s=N$s[p],f=N$f[p])} else {
      L <- intUnion(Q[1:q,])
      for(h in 1:nrow(L)){ r <- r+1
        M[r,] <- list(u=H$u,v=H$v,s=L$s[h],f=L$f[h],d=0,rel=H$rel,test=H$test)
      }
      if(p==k) break
      H <- G; i <- p; q <- 0
    }
  }
  M <- M[1:r,]; M$d <- M$f-M$s+1; tR <- sum(M$d)
  cat(" time reduced =",tR,"\n"); flush.console()
  M$s <- as.Date(M$s,origin='1970-01-01')
  M$f <- as.Date(M$f,origin='1970-01-01'); return(M)
}
 

