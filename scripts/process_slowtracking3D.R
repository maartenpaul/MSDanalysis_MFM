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

directory <- "/media/DATA/Maarten/slow_track/"
#directory <- "/media/DATA/Maarten/slow_track2/"

condition_list <- list.dirs(directory,full.names = F,recursive = F)
#condition_list <- condition_list[c(1,2)]


#msd_analyze_data_mosaic_mask(directory,condition_list,framerate,n,fitzero,min_length,pixelsize,fitMSD,offset,max_tracks,dim=dim)

segments_all <- list()
msd_fit_all <- list()
track_stats_all <- list()
library(readr)

#load data
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
write_delim(ldply(segments_all),file = file.path(directory,"segments_all.txt"))
load(file=file.path(directory,"segments_all.Rdata"))
head(segments_all$WT_200)

#calculate track stats; mark in mask if any in segment inside
stats <- llply(segments_all,function(x) {
  ddply(x,.variables = c("track","cellID"), function(x) {
    data.frame("TrackN"=nrow(x),"inMask"=any(x$inMask==TRUE))
  })
})

#convert number of frames to time (s), filter out short tracks
stats$WT_50$length <- stats$WT_50$TrackN*0.05
stats$WT_50 <- stats$WT_50[stats$WT_50$length>=0.010,]
stats$WT_200$length <- stats$WT_200$TrackN*0.25
stats$WT_200 <- stats$WT_200[stats$WT_200$length>0.5,]
stats$WT_1000$length <- stats$WT_1000$TrackN*1.05
stats$WT_1000 <- stats$WT_1000[stats$WT_1000$length>2.0,]
stats$WT_3000$length <- stats$WT_3000$TrackN*3.05
stats$WT_3000 <- stats$WT_3000[stats$WT_3000$length>5.0,]

# stats$dDBD_200$length <- stats$dDBD_200$TrackN*0.25
# stats$dDBD_200 <- stats$dDBD_200[stats$dDBD_200$length>0.4,]
# stats$dDBD_1000$length <- stats$dDBD_1000$TrackN*1.05
# stats$dDBD_1000 <- stats$dDBD_1000[stats$dDBD_1000$length>2.0,]
# stats$dDBD_3000$length <- stats$dDBD_3000$TrackN*3.05
# stats$dDBD_3000 <- stats$dDBD_3000[stats$dDBD_3000$length>6.0,]
# 
# stats$dCTD_200$length <- stats$dCTD_200$TrackN*0.25
# stats$dCTD_200 <- stats$dCTD_200[stats$dCTD_200$length>0.4,]
# stats$dCTD_1000$length <- stats$dCTD_1000$TrackN*1.05
# stats$dCTD_1000 <- stats$dCTD_1000[stats$dCTD_1000$length>2.0,]
# stats$dCTD_3000$length <- stats$dCTD_3000$TrackN*3.05
# stats$dCTD_3000 <- stats$dCTD_3000[stats$dCTD_3000$length>6.0,]

save(stats,file=file.path(directory,"slow_track_stats.Rdata"))
write_delim(ldply(stats),file = file.path(directory,"slow_track_stats.txt"))

#x50 <- hist(subset(stats$WT_50,inMask==FALSE)$length,breaks=seq(0,100,0.05),col = "black")
# x200 <- hist(subset(stats$WT_200,inMask==FALSE)$length,breaks=seq(0,100,0.2),col = "red")
# x1000 <- hist(subset(stats$WT_1000,inMask==FALSE)$length,breaks=seq(0,100,1),add=TRUE,col="blue")
# x3000 <- hist(subset(stats$WT_3000,inMask==FALSE)$length,breaks=seq(0,100,3),add=TRUE,col="green")

####all
#create frequency histogram
x50 <- hist((stats$WT_50)$length,breaks=seq(0,100,0.05),col = "red")
x200 <- hist((stats$WT_200)$length,breaks=seq(0,100,0.2),col = "red")
x1000 <- hist((stats$WT_1000)$length,breaks=seq(0,200,1),add=TRUE,col="blue")
x3000 <- hist((stats$WT_3000)$length,breaks=seq(0,400,3),add=TRUE,col="green")


# data <- rbind(data.frame("time"=x200$mids+0.05,"density"=rev(cumsum(rev(x200$counts))),tl=0.25),
#               data.frame("time"=x1000$mids+0.05,"density"=rev(cumsum(rev(x1000$counts))),tl=1.05),
#               data.frame("time"=x3000$mids+0.05,"density"=rev(cumsum(rev(x3000$counts))),tl=3.05))

#make cumulative distriibution, filter out low density
s50  <- data.frame("time"=x50$mids+0.05,"density"=rev(cumsum(rev(x50$counts))))[-c(1,2,3,4,5),]
s50 <- s50[s50$density>1,]
s200  <- data.frame("time"=x200$mids+0.05,"density"=rev(cumsum(rev(x200$counts))))[-c(1,2,3),]
s200 <- s200[s200$density>1,]
s1000 <-  data.frame("time"=x1000$mids+0.05,"density"=rev(cumsum(rev(x1000$counts))))[-c(1,2),]
s1000 <- s1000[s1000$density>1,]
s3000 <-  data.frame("time"=x3000$mids+0.05,"density"=rev(cumsum(rev(x3000$counts))))[-c(1,2),]
s3000 <- s3000[s3000$density>1,]

#make matrix to fill data
survival_matrix <- matrix(data = 0,nrow=max(c(nrow(s50),nrow(s200),nrow(s1000),nrow(s3000))),ncol=8)
survival_matrix[1:nrow(s50),1:2] <- as.matrix(s50)
survival_matrix[1:nrow(s200),3:4] <- as.matrix(s200)
survival_matrix[1:nrow(s1000),5:6] <- as.matrix(s1000)

survival_matrix[1:nrow(s3000),7:8] <- as.matrix(s3000)

#write to file
write_delim(as.data.frame(survival_matrix),file = file.path(directory,"survival functions_all.txt"),col_names = FALSE)

#inside
x50 <- hist(subset(stats$WT_50,inMask==TRUE)$length,breaks=seq(0,100,0.05),col = "red")
x200 <- hist(subset(stats$WT_200,inMask==TRUE)$length,breaks=seq(0,100,0.2),col = "red")
x1000 <- hist(subset(stats$WT_1000,inMask==TRUE)$length,breaks=seq(0,300,1),add=TRUE,col="blue")
x3000 <- hist(subset(stats$WT_3000,inMask==TRUE)$length,breaks=seq(0,400,3),add=TRUE,col="green")

s50  <- data.frame("time"=x50$mids+0.05,"density"=rev(cumsum(rev(x50$counts))))[-c(1,2,3,4,5),]
s50 <- s50[s50$density>1,]
s200  <- data.frame("time"=x200$mids+0.05,"density"=rev(cumsum(rev(x200$counts))))[-c(1,2,3),]
s200 <- s200[s200$density>1,]
s1000 <-  data.frame("time"=x1000$mids+0.05,"density"=rev(cumsum(rev(x1000$counts))))[-c(1,2),]
s1000 <- s1000[s1000$density>1,]
s3000 <-  data.frame("time"=x3000$mids+0.05,"density"=rev(cumsum(rev(x3000$counts))))[-c(1,2),]
s3000 <- s3000[s3000$density>1,]

survival_matrix <- matrix(data = 0,nrow=max(c(nrow(s50),nrow(s200),nrow(s1000),nrow(s3000))),ncol=8)
survival_matrix[1:nrow(s50),1:2] <- as.matrix(s50)
survival_matrix[1:nrow(s200),3:4] <- as.matrix(s200)
survival_matrix[1:nrow(s1000),5:6] <- as.matrix(s1000)

survival_matrix[1:nrow(s3000),7:8] <- as.matrix(s3000)


write_delim(as.data.frame(survival_matrix),file = file.path(directory,"survival functions_inside.txt"),col_names = FALSE)

x50 <- hist(subset(stats$WT_50,inMask==FALSE)$length,breaks=seq(0,100,0.05),col = "red")
x200 <- hist(subset(stats$WT_200,inMask==FALSE)$length,breaks=seq(0,100,0.2),col = "red")
x1000 <- hist(subset(stats$WT_1000,inMask==FALSE)$length,breaks=seq(0,300,1),add=TRUE,col="blue")
x3000 <- hist(subset(stats$WT_3000,inMask==FALSE)$length,breaks=seq(0,400,3),add=TRUE,col="green")


data <- rbind(data.frame("time"=x200$mids+0.05,"density"=rev(cumsum(rev(x200$counts))),tl=0.25),
              data.frame("time"=x1000$mids+0.05,"density"=rev(cumsum(rev(x1000$counts))),tl=1.05),
              data.frame("time"=x3000$mids+0.05,"density"=rev(cumsum(rev(x3000$counts))),tl=3.05))

s50  <- data.frame("time"=x50$mids+0.05,"density"=rev(cumsum(rev(x50$counts))))[-c(1,2,3,4,5),]
s50 <- s50[s50$density>1,]
s200  <- data.frame("time"=x200$mids+0.05,"density"=rev(cumsum(rev(x200$counts))))[-c(1,2,3),]
s200 <- s200[s200$density>1,]
s1000 <-  data.frame("time"=x1000$mids+0.05,"density"=rev(cumsum(rev(x1000$counts))))[-c(1,2),]
s1000 <- s1000[s1000$density>1,]
s3000 <-  data.frame("time"=x3000$mids+0.05,"density"=rev(cumsum(rev(x3000$counts))))[-c(1,2),]
s3000 <- s3000[s3000$density>1,]

survival_matrix <- matrix(data = 0,nrow=max(c(nrow(s50),nrow(s200),nrow(s1000),nrow(s3000))),ncol=8)
survival_matrix[1:nrow(s50),1:2] <- as.matrix(s50)
survival_matrix[1:nrow(s200),3:4] <- as.matrix(s200)
survival_matrix[1:nrow(s1000),5:6] <- as.matrix(s1000)

survival_matrix[1:nrow(s3000),7:8] <- as.matrix(s3000)


write_delim(as.data.frame(survival_matrix),file = file.path(directory,"survival functions_outside.txt"),col_names = FALSE)
