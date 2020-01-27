#required packages
library(plyr)
library(lattice)
library(stats)
library(MSDtracking)
library(ggplot2)
library(ggpol)

#input variables
framerate <- 1/52 #1/200 #1/ms
n <- 4 #number of timepoints taken into account for MSD fit
fitzero <- TRUE #should fit go through origin (0,0)

min_length <- 6 #minimum length track
#pixelsize <- c(1000,1000,1000) #nm
pixelsize = 1000
fitMSD <- T
offset <- 4*(0.01)^2 #experimentally determined
max_tracks <- 500 #maximum number of tracks per frame else exclude tracks from dataset, avoids mislinking of tracks
dim <- 2 #number of dimensions of tracking

directory <- "D:/OneDrive/Data2/slow_track/"
condition_list <- list.dirs(directory,full.names = F,recursive = F)
#condition_list <- condition_list[c(1,2)]


#msd_analyze_data_mosaic_mask(directory,condition_list,framerate,n,fitzero,min_length,pixelsize,fitMSD,offset,max_tracks,dim=dim)

segments_all <- list()
msd_fit_all <- list()
track_stats_all <- list()
library(readr)

for (i in 1:length(condition_list)){
  segments <- list()
  dir <- file.path(directory,condition_list[i])
  filelist <- list.files(dir,full.names = T,recursive = F,pattern = "^Traj_.*.\\mask.csv$")
  total <- length(filelist)
  # create progress bar
  for(j in 1:total){
    Sys.sleep(0.1)
    tracks_simple <- read_csv(filelist[j])
    names(tracks_simple) <- c("track","X","Y","Z","time","frame","step_x","step_y","step_z","inMask")
    #trackID,pos_x,pos_y,pos_z,time,frame,step_x,step_y,step_z,inMask
    tracks_simple <- tracks_simple[,c(6,2,3,1,4,5,7,8,9,10)]
    #3,4,5,2,6,7,8,9,10,11,12)]
    tracks_simple$frame <- tracks_simple$frame
    segments[[j]] <- data.frame(SEGMENT_STAT(tracks_simple),"cellID"=basename(filelist[j]))
  }
  segments_all[[basename(dir)]] <- data.frame(ldply(segments),"condition"=basename(dir))


}

#save data to the folder
save(segments_all,file=file.path(directory,"segments_all.Rdata"))
load(file=file.path(directory,"segments_all.Rdata"))
head(segments_all$WT_200)

stats <- llply(segments_all,function(x) {
  ddply(x,.variables = c("track","cellID"), function(x) {
    data.frame("TrackN"=nrow(x),"inMask"=any(x$inMask==TRUE))
  })
})

stats$WT_50$length <- stats$WT_50$TrackN*0.05
stats$WT_50 <- stats$WT_50[stats$WT_50$length>=0.010,]

stats$WT_200$length <- stats$WT_200$TrackN*0.25
stats$WT_200 <- stats$WT_200[stats$WT_200$length>0.4,]
stats$WT_1000$length <- stats$WT_1000$TrackN*1.05
stats$WT_1000 <- stats$WT_1000[stats$WT_1000$length>2.0,]
stats$WT_3000$length <- stats$WT_3000$TrackN*3.05
stats$WT_3000 <- stats$WT_3000[stats$WT_3000$length>6.0,]

stats$dDBD_200$length <- stats$dDBD_200$TrackN*0.25
stats$dDBD_200 <- stats$dDBD_200[stats$dDBD_200$length>0.4,]
stats$dDBD_1000$length <- stats$dDBD_1000$TrackN*1.05
stats$dDBD_1000 <- stats$dDBD_1000[stats$dDBD_1000$length>2.0,]
stats$dDBD_3000$length <- stats$dDBD_3000$TrackN*3.05
stats$dDBD_3000 <- stats$dDBD_3000[stats$dDBD_3000$length>6.0,]

stats$dCTD_200$length <- stats$dCTD_200$TrackN*0.25
stats$dCTD_200 <- stats$dCTD_200[stats$dCTD_200$length>0.4,]
stats$dCTD_1000$length <- stats$dCTD_1000$TrackN*1.05
stats$dCTD_1000 <- stats$dCTD_1000[stats$dCTD_1000$length>2.0,]
stats$dCTD_3000$length <- stats$dCTD_3000$TrackN*3.05
stats$dCTD_3000 <- stats$dCTD_3000[stats$dCTD_3000$length>6.0,]



x50 <- hist(subset(stats$WT_50,inMask==TRUE)$length,breaks=seq(0,100,0.05),col = "black")
x200 <- hist(subset(stats$dDBD_200,inMask==TRUE)$length,breaks=seq(0,100,0.05),add=TRUE,col = "red")
x1000 <- hist(subset(stats$dDBD_1000,inMask==TRUE)$length,breaks=seq(0,100,0.05),add=TRUE,col="blue")
x3000 <- hist(subset(stats$dDBD_3000,inMask==TRUE)$length,breaks=seq(0,100,0.05),add=TRUE,col="green")

data <- rbind(data.frame("time"=x50$mids+0.05,"density"=x50$counts,tl=0.05),
              data.frame("time"=x200$mids+0.05,"density"=x200$counts,tl=0.25),
              data.frame("time"=x1000$mids+0.05,"density"=x1000$counts,tl=1.05),
              data.frame("time"=x3000$mids+0.05,"density"=x3000$counts,tl=3.05))

data$interval <- as.character(data$tl)
data <- data[data$density>10,]
data <- ddply(data,.variables="interval", function(x){
  x$weight <- 1/nrow(x)
  return(x)
})
tint = 0.05


library(reticulate)
#use_virtualenv("base")
scipy <- import("scipy")
lmfit <- import('lmfit')

py$time <- as.vector(data$time)
py$density <- as.vector(data$density)
py$tl <- as.vector(data$tl)
py$weight <- as.vector(data$weight)

library(minpack.lm)
result <- nlsLM(formula=density~A*exp(-((kb*(tint/tl)+koff1)*time)),
                start=list(koff1=5,kb=1,A=1),
                lower=c(0.00001,0,0),upper = c(100,100,10),data = data,trace = T,weights = data$weight)
result


#formula=density~A*(B*(kb*(tint/tl)+koff1)*exp(-((kb*(tint/tl)+koff1)*time))+
#                     C*(kb*(tint/tl)+koff2)*exp(-((kb*(tint/tl)+koff2)*time))+
#                     (1-B-C)*(kb*(tint/tl)+koff3)*exp(-(kb*((tint/tl)+koff3)*time)))


#formula=density~A*((B*(kb*(tint/tl)+koff1)*exp(-(kb*(tint/tl)+koff1)*time))+
#((1-B)*(kb*(tint/tl)+koff2)*exp(-(kb*(tint/tl)+koff2)*time))),

#result <- nlsLM(formula=density~A*(B*exp(-(kb*(tint/tl)+koff1)*time)+(1-B)*exp(-(kb*(tint/tl)+koff2)*time)),
#                start=list(koff1=0.0001,koff2=1,kb=1,A=2,B=0.5),
#                lower=c(0.00001,0.5,0.01,0,0),upper = c(30,30,10,5,1),data = data)

source_python("python/lsqfit.py")
plot(data$time,data$density,col="black",xlab="track length (seconds)",ylab="frequency",log="xy")
lines((data$time[data$tl==.05]),(data$density[data$tl==.05]), main = "data",col="red",ylab="log10 density")
lines((data$time[data$tl==.25]),(data$density[data$tl==.25]), main = "data",col="green",xlab="log10 track length (seconds)",ylab="log10 density")
lines((data$time[data$tl==1.05]),(data$density[data$tl==1.05]), main = "data",col="blue",xlab="log10 track length (seconds)",ylab="log10 density")
lines((data$time[data$tl==3.05]),(data$density[data$tl==3.05]), main = "data",col="orange",xlab="log10 track length (seconds)",ylab="log10 density")

plot(data$time,data$density,col="black",xlab="track length (seconds)",ylab="frequency",log="")
lines((data$time[data$tl==.05]),(data$density[data$tl==.05]), main = "data",col="red",ylab="log10 density")
lines((data$time[data$tl==.25]),(data$density[data$tl==.25]), main = "data",col="green",xlab="log10 track length (seconds)",ylab="log10 density")
lines((data$time[data$tl==1.05]),(data$density[data$tl==1.05]), main = "data",col="blue",xlab="log10 track length (seconds)",ylab="log10 density")
lines((data$time[data$tl==3.05]),(data$density[data$tl==3.05]), main = "data",col="orange",xlab="log10 track length (seconds)",ylab="log10 density")



lines(log10(data$time[data$tl==.05]), log10(py$model[data$tl==0.05]), col = 1, lwd = 2)
lines(log10(data$time[data$tl==0.25]), log10(py$model[data$tl==0.25]), col = 2, lwd = 2)
lines(log10(data$time[data$tl==1.05]), log10(py$model[data$tl==1.05]), col = 3, lwd = 2)
lines(log10(data$time[data$tl==3.05]), log10(py$model[data$tl==3.05]), col = 4, lwd = 2)





plot(data$time,data$density, main = "data",col="black")


lines(data$time[data$tl==.05], fitted(result)[data$tl==0.05], col = 1, lwd = 2)
lines(data$time[data$tl==0.25], fitted(result)[data$tl==0.25], col = 2, lwd = 2)
lines(data$time[data$tl==1.05], fitted(result)[data$tl==1.05], col = 3, lwd = 2)
lines(data$time[data$tl==3.05], fitted(result)[data$tl==3.05], col = 4, lwd = 2)

lines(log10(data$time[data$tl==.05]), log10(fitted(result)[data$tl==0.05]), col = 1, lwd = 2)
lines(log10(data$time[data$tl==0.25]), log10(fitted(result)[data$tl==0.25]), col = 2, lwd = 2)
lines(log10(data$time[data$tl==1.05]), log10(fitted(result)[data$tl==1.05]), col = 3, lwd = 2)
lines(log10(data$time[data$tl==3.05]), log10(fitted(result)[data$tl==3.05]), col = 4, lwd = 2)


lines(data$time[data$tl==.05], py$model[data$tl==0.05], col = 1, lwd = 2)
lines(data$time[data$tl==0.25], py$model[data$tl==0.25], col = 2, lwd = 2)
lines(data$time[data$tl==1.05], py$model[data$tl==1.05], col = 3, lwd = 2)
lines(data$time[data$tl==3.05], py$model[data$tl==3.05], col = 4, lwd = 2)




resFun <- f

nls.lm()
