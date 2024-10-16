# TQvis - visualization of temporal quantities
# October 2024 by Vladimir Batagelj
# http://localhost:8800/doku.php?id=work:alg:icon
# https://github.com/bavla/TQ/tree/master/R
# http://localhost:8800/doku.php?id=work:alg:cutq#python_to_r
# http://localhost:8800/doku.php?id=work:alg:cutq#importing_in_r
# C:/Users/vlado/work/OpenAlex/vis

TQlistGet <- function(D,i){
  n <- length(D[[i]][[2]])
  if(n==0) return(NULL)
  if(n==1) tq <- as.data.frame(t(D[[i]][[2]][[1]])) else
    tq <- data.frame(Reduce(rbind,D[[i]][[2]]))
  colnames(tq) <- c("s","f","v"); row.names(tq) <- paste("e",1:nrow(tq),sep="")
  return(tq)
}

TQlistLab <- function(D,i) return(D[[i]][[1]])

TQmaxVal <- function(D){
  maxV <- 0
  for(i in 1:length(D)){
    n <- length(D[[i]][[2]])
    if(n==0) next
    if(n==1) tq <- as.data.frame(t(D[[i]][[2]][[1]])) else
      tq <- data.frame(Reduce(rbind,D[[i]][[2]]))
    maxV <- max(c(maxV,tq[,3]),na.rm=TRUE)
  }
  return(maxV)
}

TQspan <- function(tq){
  span <- 0
  for(i in 1:nrow(tq)) span <- span + (tq[i,2]-tq[i,1])
  return(span)
}

TQtotal <- function(tq){
  tot <- 0
  for(i in 1:nrow(tq)) tot <- tot + (tq[i,2]-tq[i,1])*tq[i,3]
  return(tot)
}

TQmax <- function(tq) return(ifelse(is.null(tq),0,max(tq[,3])))

tqComps <- function(tq,tMin,tMax,sent=NULL,trans){
  n <- nrow(tq)
  if(is.null(sent)) sent <- max(tq[,3],na.rm=TRUE)+1 
  tq <- rbind(tq,c(tq[n,2],tMax,sent)); n <- n+1
  x <- h <- w <- c(); last <- tMin
  for(i in 1:n){
    s <- tq[i,1]; f <- tq[i,2]; v <- ifelse(i<n,round(trans(tq[i,3])),sent)  
    if(last < s){x <- c(x,last-tMin); w <- c(w,s-last); h <- c(h,sent)}
    x <- c(x,s-tMin); w <- c(w,f-s); h <- c(h,v); last <- f   
  }
  return(list(x=x,w=w,h=h))
}

TQicons <- function(D,I,tMin,tMax,sent,col,type=1,g=0.15,pts=100,f=function(x) x){
  grid.newpage()
  dt <- tMax-tMin
  for(i in 1:length(I)){
    y0 <- 1.005-0.1*i; j <- I[i]; lab <- TQlistLab(D,j); tq <- TQlistGet(D,j)
    xwh <- tqComps(tq,tMin,tMax,sent=sent,trans=f)
    # cat(i,j,lab,"\n"); cat(xwh$x,"\n",xwh$w,"\n",xwh$h,"\n")
    tqLab <- textbox_grob(lab,x=0,y=y0,width=unit(pts,"pt"),
      hjust=0,vjust=0,halign=0,valign=0,gp=gpar(fontsize=10),
      box_gp = gpar(col="black",fill="lightcyan1", 
        just=c("left","bottom"),lty="blank"),
      padding=unit(c(2,2,2,5),"pt"),margin=unit(c(0,5,0,5),"pt"),
      r = unit(5,"pt"),name = paste0("T",j)
    )
    m <- length(xwh$x); xx <- g + xwh$x*(0.99-g)/dt; ww <- xwh$w*(0.99-g)/dt
    if(type==1){ C <- col[xwh$h]; h <- 0.09 } else {
      C <- rep(col[1],m); C[xwh$h==sent] <- col[2]; h <- xwh$h*0.09/sent }
    tqIcon <- rectGrob(x=xx,y=y0,width=ww,height=h,
      gp=gpar(lty="blank",fill=C),just=c("left","bottom"))    
    grid.draw(gList(tqLab,tqIcon))
  }
}


 
