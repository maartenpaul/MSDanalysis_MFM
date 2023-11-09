library(tidyverse)

binsize<-0.05
bins<- seq(0,3,binsize)
results <- rep(0,length(bins))
i <- 0
for(bin in bins){
  i<-i+1
  data <- dplyr::filter(segments_all$`WT MMC`,distanceToNearest2D>=bin&distanceToNearest2D<(bin+binsize)&inMask==F)
  results[i] <- mean(data$radial_displacement1,na.rm = T)
  
}
results <- results#*(bins+0.5*binsize)
plot(x = (bins+0.5*binsize),y = results,type = "l",xlab="Radial position (um)",ylab=("r*dr (um2)"),col="red",ylim=c(-0.20,0.20))
abline(h=0,lty=2)

#gt1
bins<- seq(0,3,binsize)
results <- rep(0,length(bins))
i <- 0
for(bin in bins){
  i<-i+1
  data <- dplyr::filter(segments_all$`WT MMC`,distanceToNearest2D_gt1>bin&distanceToNearest2D_gt1<(bin+binsize))
  results[i] <- mean(data$radial_displacement1_gt1,na.rm = T)
  
}
results <- results*(bins+0.5*binsize)
lines(x = (bins+0.5*binsize),y = results,,col="blue",add=T)

##only centrioid in foci
bins<- seq(0,2,binsize)
results <- rep(0,length(bins))
i <- 0
for(bin in bins){
  i<-i+1
  data <- dplyr::filter(segments_all$`WT MMC`,distanceCentroid2D>bin&distanceCentroid2D<(bin+binsize))
  results[i] <- mean(data$radial_displacement1,na.rm = T)*mean(data$distanceToNearest2D,na.rm = T)
  
}
results <- results*(bins+0.5*binsize)
plot(x = (bins+0.5*binsize),y = results,col="blue",typ="l",xlim=c(0,1),ylim=c(-0.01,0.01))



##only centrioid in foci
bins<- seq(0,2,binsize)
results <- rep(0,length(bins))
i <- 0
for(bin in bins){
  i<-i+1
  data <- dplyr::filter(segments_all$`WT MMC`,distanceCentroid2D_gt1>bin&distanceCentroid2D_gt1<(bin+binsize))
  results[i] <- mean(data$radial_displacement1_gt1,na.rm = T)*bins[i]
  
}
results <- results*(bins+0.5*binsize)
plot(x = (bins+0.5*binsize),y = results,col="blue",typ="l",xlim=c(0,1),ylim=c(-0.01,0.01))
