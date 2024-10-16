# R TQ

```
> setwd(wdir <- "C:/Users/vlado/work/OpenAlex/vis")
> library(grid); library(rjson); library(gridtext)
> source("C:/Users/vlado/work/OpenAlex/vis/TQvis.R")
> ddir <- "C:/Users/vlado/work2/python/pqCores"
> Jfile <- paste0(ddir,"/T50wdegNodes.json")
> DN <- fromJSON(file=Jfile)
> 
> # fun <- function(x) 1+log(x)
> # fun <- function(x) x
> fun <- sqrt
> tMin <- 1; tMax <- 67; dt <- tMax-tMin
> vmax <- TQmaxVal(DN)
> vMax <- round(fun(vmax)); sent <- vMax+1
> Fire <- c(heat.colors(vMax)[vMax:1],"#00FF00")
>
> TQicons(DN,1:10,tMin=tMin,tMax=tMax,sent=sent,col=Fire,type=1,f=fun)
```
