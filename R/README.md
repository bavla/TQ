# R TQ

```
> setwd(wdir <- "C:/work/OpenAlex/vis")
> library(grid); library(rjson); library(gridtext)
> source("C:/work/OpenAlex/vis/TQvis.R")
> ddir <- "C:/work/python/pqCores"
> Jfile <- paste0(ddir,"/T50wdegN.json")
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
![T50_1-10a](https://github.com/user-attachments/assets/99e43f15-2907-40fc-aa08-644c99ece6d5)
