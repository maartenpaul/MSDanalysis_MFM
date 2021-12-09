#required packages
library(reticulate)
use_condaenv("R",required = T)
library(plyr)
library(lattice)
library(stats)
library(MSDtracking)
library(ggplot2)#required packages
library(ggpol)
library(doParallel)


source('R/MSD.R')
source('R/MSD_fit.R')
source('R/analysis functions.R')
source('python/ML_py.R')

directory <- "D:/Stack/Genetics/data3/"

load(file.path(directory,"msd_fit_all.Rdata"))
load(file.path(directory,"segments_all_angles.Rdata"))

segments_all2 <- as_data_frame(ldply(segments_all))
segments_all <- segments_all2
rm(segments_all2)

segments_all <- ddply(segments_all,.variables = c("condition","cellID","track"),.parallel = F,function(x){
  segment <- 1
  tracklets <- vector(length=nrow(x))
  tracklets[1] <- segment
  
  for (i in 2:nrow(x)){
    if(x$state[i]==x$state[i-1]){
      tracklets[i] <-segment
    } else {
      segment <- segment+1
      tracklets[i] <-segment
    }
  }
  x$tracklet <- paste0(x$track,".",tracklets)
  
  return(x)
  
  })
    

numPmsd <- 4
numPmss <- 4
minLen <- 10

p <- seq(from=0.5,to=6,length.out=12)

segments_all$`WT MMC` <- ddply(segments_all$`WT MMC`,.variables = c("cellID","tracklet"),.parallel = F,function(x){
  if(nrow(x)>minLen){
  getMSDandMSS(x$X,x$Y,numPmsd,numPmss,p)
  }
  return(x)
  
})


stopCluster(cl)
