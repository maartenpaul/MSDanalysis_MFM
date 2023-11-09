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
<<<<<<< HEAD
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
=======
results <- results
plot(x = (bins+0.5*binsize),y = results,type = "l",xlab="Radial position (um)",ylab=("<dr> (um)"),col="red",ylim=c(-0.020,0.02))
abline(h=0,lty=2)

#gt1
bins<- seq(0,3,binsize)
>>>>>>> d0489fe33dd544f3c894d8d2a7a7b4466d576435
results <- rep(0,length(bins))
i <- 0
for(bin in bins){
  i<-i+1
<<<<<<< HEAD
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
=======
  data <- dplyr::filter(segments_all$`WT MMC`,distanceToNearest2D_gt1>bin&distanceToNearest2D_gt1<(bin+binsize))
  results[i] <- mean(data$radial_displacement1_gt1,na.rm = T)
  
}
results <- results
lines(x = (bins+0.5*binsize),y = results,,col="blue")

legend(2., 0.017, legend=c("Data", "Random"),
       col=c("red", "blue"), lty=1:2, cex=0.8)


#####
library(tidyverse)

binsize<-0.05
bins<- seq(0,3,binsize)
results <- array(dim = c(length(bins),2),data = 0)
i <- 0
for(bin in bins){
  i<-i+1
  data <- dplyr::filter(segments_all$`WT MMC`,distanceToNearest2D>=bin&distanceToNearest2D<(bin+binsize)&inMask==F)
  results[i,1] <- mean(data$radial_displacement1,na.rm = T)
  results[i,2] <- mean(data$distanceToNearest2D,na.rm = T)
  
  
}
results[,1] <- results[,1]*results[,2]
plot(x = (bins+0.5*binsize),y = results[,1],type = "l",xlab="Radial position (um)",ylab=("<dr>r (um2)"),col="red",ylim=c(-0.020,0.02))
abline(h=0,lty=2)

#gt1
bins<- seq(0,3,binsize)
results <- array(dim = c(length(bins),2),data = 0)
i <- 0
for(bin in bins){
  i<-i+1
  data <- dplyr::filter(segments_all$`WT MMC`,distanceToNearest2D_gt1>=bin&distanceToNearest2D_gt1<(bin+binsize)&inMask==F)
  results[i,1] <- mean(data$radial_displacement1_gt1,na.rm = T)
  results[i,2] <- mean(data$distanceToNearest2D_gt1,na.rm = T)
  
  
}
results[,1] <- results[,1]*results[,2]
lines(x = (bins+0.5*binsize),y = results[,1],,col="blue")

legend(2., 0.017, legend=c("Data", "Random"),
       col=c("red", "blue"), lty=1:2, cex=0.8)




>>>>>>> d0489fe33dd544f3c894d8d2a7a7b4466d576435
