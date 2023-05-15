library(tidyverse)
plot(segments_all$`WT MMC`$radial_displacement1,segments_all$`WT MMC`$distanceToNearest3D)
ggplot(segments_all$`WT MMC`,aes(distanceToNearest3D,radial_displacement1))+geom_line(stat="summary_bin",binwidth=0.01)

binsize<-0.01
bins<- seq(0,0.5,binsize)
results <- rep(0,length(bins))
i <- 0
for(bin in bins){
  i<-i+1
  data <- dplyr::filter(segments_all$`WT MMC`,distanceToNearest3D>bin&distanceToNearest3D<(bin+binsize))
  results[i] <- mean(data$radial_displacement1,na.rm = T)
  
}
results <- results*(bins+0.5*binsize)
plot(x = (bins+0.5*binsize),y = results,type = "l")

bins<- seq(0,0.5,binsize)
results <- rep(0,length(bins))
i <- 0
for(bin in bins){
  i<-i+1
  data <- dplyr::filter(segments_all$`WT MMC`,distanceToNearest3D_gt1>bin&distanceToNearest3D_gt1<(bin+binsize))
  results[i] <- mean(data$radial_displacement1,na.rm = T)
  
}

plot(x = (bins+0.5*binsize),y = results,type = "l")
