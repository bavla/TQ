# R TQ

```
> setwd(wdir <- "C:/work/OpenAlex/vis")
> library(grid); library(rjson); library(gridtext)
> source("C:/work/OpenAlex/vis/TQvis.R")
> ddir <- "C:/work/python/pqCores"
> Jfile <- paste0(ddir,"/T50wdegNodes.json")
> DN <- fromJSON(file=Jfile)
> 
> # fun <- function(x) 1+log(x)
> # fun <- function(x) x
> fun <- sqrt
> tMin <- 1; tMax <- 67; dt <- tMax-tMin
> vMax <- round(fun(TQmaxVal(DN))); sent <- vMax+1
> Fire <- c(heat.colors(vMax)[vMax:1],"#00FF00")
>
> TQicons(DN,1:10,tMin=tMin,tMax=tMax,sent=sent,col=Fire,type=1,f=fun)
```
![T50_1-10a](https://github.com/user-attachments/assets/3a1d6ca3-8bbf-4be9-bee1-461620b521b3)


```
> TQicons(DN,1:10,tMin=tMin,tMax=tMax,sent=sent,col=c("red","lightcyan1"),type=2,f=sqrt)
```


![T50_1-10b](https://github.com/user-attachments/assets/77dde1c1-0ee7-441f-a6d0-25526500c0c2)
