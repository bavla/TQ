# trajectoR
# Vladimir Batagelj, December 2023 / January 2024
# Network of time interval intersections / equal part_of_cv and institution_si

datum <- function(d,m,y) {
  dd <- ifelse(d!="",d,ifelse(is.na(m),paste("01/01/",y,sep=""),
    paste("01/",m,"/",y,sep="")))
  di <- as.integer(as.Date(dd,format="%d/%m/%Y",origin="1970-01-01"))
  if(is.na(di)) {OK <<- FALSE
    di <- as.integer(as.Date(dd,format="%m/%d/%Y",origin="1970-01-01"))}
  return(di)
}

raw2events <- function(D,CSV="events.csv"){
  csv <- file(CSV,"w",encoding="UTF-8")
  cat("R  ID  s  f  S  T\n",file=csv)
  n <- nrow(D); OK <- TRUE 
  for(i in 1:n){
    ID <- D$person_name[i]; rel <- tolower(D$part_of_cv[i])
    test <- trimws(tolower(D$institution_si[i]))
    ds <- D$start_day[i]; ms <- D$start_month[i]; ys <- D$start_year[i]
    sd <- datum(ds,ms,ys)
    if(!OK) {cat(i,":",ID,ds,ms,ys,rel,'*** wrong date\n')
      flush.console(); OK <- TRUE}
    de <- D$end_day[i]; me <- D$end_month[i]; ye <- D$end_year[i]
    ed <- if(ye==2100) datum("01/01/2024",NA,2024) else datum(de,me,ye)
    if(!OK) {cat(i,":",ID,de,me,ye,rel,'*** wrong date\n')
      flush.console(); OK <- TRUE}
    cat(i,' "',ID,'" ',sd,' ',ed,' "',rel,'" "',test,'"\n',sep='',file=csv)
  }
  close(csv)
}

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
      if(E$S[i]==E$S[j]) if(E$T[i]==E$T[j]) {
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

explain <- function(E,p1,p2,kMax=500){
  i12 <- c(which(E$ID == p1),which(E$ID == p2))
  C <- E[i12,]; n <- nrow(C); I <- order(C$s,C$f); k <- 0; r <- 0
  cn <- rep(0,kMax); cs <- rep("",kMax)
  N <- data.frame(u=cs,v=cs,s=cn,f=cn,rel=cs,d=cn,Ru=cs,Rv=cs)
  for(p in 1:(n-1)){i <- I[p]; tm <- C$f[i]
    # cat(p,k,date(),"\n"); flush.console()
    for(q in (p+1):n){j <- I[q]; r <- r+1
      if(C$s[j]>tm) break
      # if(C$R[i]!=C$R[j]) 
      if(C$S[i]==C$S[j]) if(C$T[i]==C$T[j]) {
        fm <- min(C$f[i],C$f[j]); sM <- max(C$s[i],C$s[j]); T <- fm-sM
        if(T>0){k <- k+1; if(k>kMax) stop("kMax too small")
          N[k,] <- list(u=C$ID[i],v=C$ID[j],s=sM,f=fm,rel=C$S[i],d=T,
                        Ru=C$R[i],Rv=C$R[j])}
      }
    }
  }
  tA <- sum(N$d); tS <- sum(N[N$u!=N$v,"d"])
  cat(p1,":",p2,"\n",date(),"\n density R =",2*r/n/(n-1),
    "  tests =",r,"\n density E =",2*k/n/(n-1),"  edges =",k,
    "\n time all =",tA,"  time strict =",tS,"\n")
  return(N[1:k,])
}

