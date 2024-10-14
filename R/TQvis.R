# TQvis - visualization of temporal quantities
# October 2024 by Vladimir Batagelj
# http://localhost:8800/doku.php?id=work:alg:icon
# https://github.com/bavla/TQ/tree/master/R

TQmaxVal <- function(D){
  maxV <- 0
  for(i in 1:length(D)){
    tq <- data.frame(Reduce(rbind,D[[i]][[2]]))
    maxV <- max(c(maxV,tq[,3]),na.rm=TRUE)
  }
  return(maxV)
}

tqComps <- function(tq,tMin,tMax,sent=NULL,trans){
  n <- nrow(tq)
  if(is.null(sent)) sent <- max(tq[,3],na.rm=TRUE)+1 
  tq <- rbind(tq,c(tq[n,2],tMax,sent)); n <- n+1
  x <- h <- w <- c(); last <- tMin
  for(i in 1:n){
    s <- tq[i,1]; f <- tq[i,2]; v <- round(trans(tq[i,3]))  
    if(last < s){x <- c(x,last-tMin); w <- c(w,s-last); h <- c(h,sent)}
    x <- c(x,s-tMin); w <- c(w,f-s); h <- c(h,v); last <- f   
  }
  return(list(x=x,w=w,h=h))
}

TQicons <- function(D,I,sent,col,type=1,g=0.15,pts=100,f=function(x) x){
  grid.newpage()
  for(i in 1:length(I)){
    y0 <- 1.005-0.1*i; j <- I[i]; lab <- D[[j]][[1]]
    tq <- data.frame(Reduce(rbind,D[[j]][[2]]))
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


 
