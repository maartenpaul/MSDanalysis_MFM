#required packages
library(reticulate)
#use_condaenv("R",required = T)
library(tidyverse)
library(plyr)
library(lattice)
library(stats)
library(MSDtracking) #https://github.com/maartenpaul/MSDtracking
library(ggplot2)
library(ggpol)
library(doParallel)

source('R/MSD.R')
source('R/MSD_fit.R')
source('R/analysis functions.R')
source('python/ML_py.R')

#input variables
framerate <- 1/52 #1/200 #1/ms
n <- 4 #number of timepoints taken into account for MSD fit
fitzero <- TRUE #should fit go through origin (0,0)
pixelsize = 1000
fitMSD <- T
offset <- 4*(0.01)^2 #experimentally determined
max_tracks <- 500 #maximum number of tracks per frame else exclude tracks from dataset, avoids mislinking of tracks
dim <- 2 #number of dimensions of tracking

directory <- "/media/DATA/Maarten/MFM/data_gtv2/"
condition_list <- list.dirs(directory,full.names = F,recursive = F)

#import data organize in lists of the different data sets
msd_analyze_data_mosaic_mask_parallel_intensity(directory,condition_list,framerate,n,fitzero,min_length,pixelsize,fitMSD,offset,max_tracks,dim=dim,groundtruth=TRUE)

load(file.path(directory,"msd_fit_all.Rdata"))
load(file.path(directory,"segments_all.Rdata"))

#save segments to text files
for (i in 1:length(msd_fit_all)){
  inner_join(segments_all[[i]],select(msd_fit_all[[i]],c(.id,track,D)),by=c(".id","track"),keep=FALSE) %>%
    write_tsv(file.path(directory,paste0(names(msd_fit_all)[i],"_segmentsD.txt")))
}

#estimate states using MLMSS
segments_all <- llply(segments_all,function(x){
  ddply(x, .variables = "cellID", function(x){
    ML_segment_tracks(x)
  })
})

#calculate angle displacements
ptm <- proc.time()
#initialize cluster
nodes <- detectCores()
cl <- makeCluster(nodes-10)
registerDoParallel(cl)
for (j in 1:length(segments_all)){
  segments_all[[j]] <- ddply(segments_all[[j]],.variables = c("cellID"), .parallel = T, function(x){
    get_angle_3D <- function(A,B,C){
      seg_angle <- vector()
      AB <- B[1:2]-A[1:2]
      CB <- C[1:2]-B[1:2]
      
      #dAB <- sqrt((B[1]-A[1])^2+(B[2]-A[2])^2)
      #dBC <- sqrt((C[1]-B[1])^2+(C[2]-B[2])^2)
      #Formula obtained from https://gitlab.com/anders.sejr.hansen/anisotropy
      angle <- abs(atan2(det(cbind(AB,CB)),AB%*%CB))
      angle <- angle/pi*180
      return(angle)
      
      
    } 
    get_angles <- function(x,n){
      x$frame  <- x$frame-x$frame[1]+1
      angles <- rep(x=-1,nrow(x))
      if(nrow(x)>=(3*n+n-1)){
        
        #loop over all steps of tracks
        for (k in 1:(nrow(x)-2*n)){
          
          if (is.element(x[k,2]+n,x$frame)&&is.element(x[k,2]+2*n,x$frame)){ #check if point is present in track, this deals with gaps
            x1 <- as.numeric(x[k,c(3,4,6)])
            which_point1 <- which(x[,2]==x[k,2]+n)
            x2 <- as.numeric(x[which_point1,c(3,4,6)])
            which_point2 <- which(x[,2]==x[k,2]+2*n)
            x3 <- as.numeric(x[which_point2,c(3,4,6)])
            angles[which_point1] <- get_angle_3D(x1,x2,x3)
            
          }}
      }
      return(angles)
    }
    get_displacements <- function(x,n){
      x$frame  <- x$frame-x$frame[1]+1
      displ <- cbind(rep(x=-1,nrow(x)),rep(x=-1,nrow(x)))
      if(nrow(x)>=(3*n+n-1)){
        #loop over all steps of tracks
        for (k in 1:(nrow(x)-2*n)){
          
          if (is.element(x[k,2]+n,x$frame)&&is.element(x[k,2]+2*n,x$frame)){ #check if point is present in track, this deals with gaps
            x1 <- as.numeric(x[k,c(3,4,6)])
            which_point1 <- which(x[,2]==x[k,2]+n)
            x2 <- as.numeric(x[which_point1,c(3,4,6)])
            which_point2 <- which(x[,2]==x[k,2]+2*n)
            x3 <- as.numeric(x[which_point2,c(3,4,6)])
            displ[which_point1,] <- c( sqrt((x1[1]-x2[1])^2+(x1[2]-x2[2])^2) ,sqrt((x2[1]-x3[1])^2+(x2[2]-x3[2])^2) )
            
          }}
      }
      return(displ)
    }
    
  
    ddply(x,.variables = c("track"), .parallel = F, function(x){
      
      
      
      x$angle1 <- get_angles(x,1)
      x[c("displacement1","displacement2")] <- get_displacements(x,1)
      x$angle2 <- get_angles(x,2)
      x$angle3 <- get_angles(x,3)
      x$angle4 <- get_angles(x,4)
      x$angle5 <- get_angles(x,5)
      x$angle6 <- get_angles(x,6)
      x$angle7 <- get_angles(x,7)
      x$angle8 <- get_angles(x,8)
      x$angle9 <- get_angles(x,9)
      x$angle10 <- get_angles(x,10)
      x$angle11 <- get_angles(x,11)
      x$angle12 <- get_angles(x,12)
      x$angle13 <- get_angles(x,13)
      x$angle14 <- get_angles(x,14)
      x$angle15 <- get_angles(x,15)
      x$angle16 <- get_angles(x,16)
      x$angle17 <- get_angles(x,17)
      x$angle18 <- get_angles(x,18)
      x$angle19 <- get_angles(x,19)
      x$angle20 <- get_angles(x,20)
      
      
      return(x)
      
    })
    
  })
}
stopCluster(cl)
proc.time() - ptm

save(segments_all,file=file.path(directory,"segments_all_angles.Rdata"))
load(file=file.path(directory,"segments_all_angles.Rdata"))


#######calculate MSD and MSS
####add tracklet column
segments_all <- as_data_frame(ldply(segments_all))

for (j in 1:length(segments_all)){
  
  segments_all[[j]] <- ddply(segments_all[[j]],.variables = c("cellID"),.parallel = T,function(x){
    ddply(x,.variables = c("track"), .parallel = F, function(x){
      
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
      
    })})
}

stopCluster(cl)
proc.time() - ptm

save(segments_all,file=file.path(directory,"segments_all_angles.Rdata"))

####split tracks in inside outside column
for (j in 1:length(segments_all)){
  
  segments_all[[j]] <- ddply(segments_all[[j]],.variables = c("cellID"),.parallel = T,function(x){
    ddply(x,.variables = c("track"), .parallel = F, function(x){
      
      segment <- 1
      tracklets <- vector(length=nrow(x))
      tracklets[1] <- segment
      
      for (i in 2:nrow(x)){
        if(x$inMask[i]==x$inMask[i-1]){
          tracklets[i] <-segment
        } else {
          segment <- segment+1
          tracklets[i] <-segment
        }
      }
      x$focus_tracklet <- paste0(x$track,".",tracklets)
      
      return(x)
      
    })})
}

####split tracks in inside outside gt1 column
for (j in 1:length(segments_all)){
  
  segments_all[[j]] <- ddply(segments_all[[j]],.variables = c("cellID"),.parallel = T,function(x){
    ddply(x,.variables = c("track"), .parallel = F, function(x){
      
      segment <- 1
      tracklets <- vector(length=nrow(x))
      tracklets[1] <- segment
      
      for (i in 2:nrow(x)){
        if(x$inMask_gt1[i]==x$inMask_gt1[i-1]){
          tracklets[i] <-segment
        } else {
          segment <- segment+1
          tracklets[i] <-segment
        }
      }
      x$focus_tracklet_gt1 <- paste0(x$track,".",tracklets)
      
      return(x)
      
    })})
}

save(segments_all,file=file.path(directory,"segments_all_angles.Rdata"))
load(file=file.path(directory,"segments_all_angles.Rdata"))

###get msd and mss from tracklets
numPmsd <- 4
numPmss <- 4
minLen <- 5
p <- seq(from=0.5,to=6,length.out=12)
py$pixSize <- 0.100
py$t <- 0.052
source_python('python/getMSDandMSS_R.py')

MSD_MSS <- function(x){
  if(nrow(x)>minLen){
    out <- getMSDandMSS_R(x$X,x$Y)
    return(tibble("D_ML"=out[[1]],"D_Smmss"=out[[2]]))
  } else {
    return(tibble("D_ML"=-1.0,"D_Smmss"=-1.0))
    return(NA)
  }
}

MSD_MSS_focus <- function(x){
  if(nrow(x)>minLen){
    out <- getMSDandMSS_R(x$X,x$Y)
    return(tibble("D_ML_focus"=out[[1]],"D_Smmss_focus"=out[[2]]))
  } else {
    return(tibble("D_ML_focus"=-1.0,"D_Smmss_focus"=-1.0))
    return(NA)
  }
}

MSD_MSS_focus_gt1 <- function(x){
  if(nrow(x)>minLen){
    out <- getMSDandMSS_R(x$X,x$Y)
    return(tibble("D_ML_focus_gt1"=out[[1]],"D_Smmss_focus_gt1"=out[[2]]))
  } else {
    return(tibble("D_ML_focus_gt1"=-1.0,"D_Smmss_focus_gt1"=-1.0))
    return(NA)
  }
}

MSD_only <- function(x){
  if(nrow(x)>minLen){
    out <- getMSDandMSS_R(x$X,x$Y)
    return(out[[1]])
  } else {
    return(NA)
  }
}

segs_nest <- ldply(segments_all)
segs_nest <-segs_nest%>%
  filter(condition=="WT MMC")%>%
  select(condition,cellID,focus_tracklet,X,Y) %>%
  group_by(condition,cellID,focus_tracklet) %>%
  group_modify(~MSD_MSS_focus(.x)) %>%
  inner_join(y=segs_nest,by=c("condition","cellID","focus_tracklet")) %>%
  ungroup()

#gt1
segs_nest <-segs_nest%>%
  select(condition,cellID,focus_tracklet_gt1,X,Y) %>%
  group_by(condition,cellID,focus_tracklet_gt1) %>%
  group_modify(~MSD_MSS_focus_gt1(.x)) %>%
  inner_join(y=segs_nest,by=c("condition","cellID","focus_tracklet_gt1")) %>%
  ungroup()

segs_nest <- segs_nest %>%
  select(condition,cellID,tracklet,X,Y) %>%
  group_by(condition,cellID,tracklet) %>%
  group_modify(~MSD_MSS(.x),.keep=T) %>%
  inner_join(y=segs_nest,by=c("condition","cellID","tracklet")) 

segs_nest$condition <- droplevels(segs_nest$condition)

save(segs_nest,file=file.path(directory,"segs_nest.Rdata"))
segs_nest %>%
  nest(-.id) %>%
  pwalk(~write_delim(x = .y, file = file.path(directory,paste0(.x, ".txt") )) )
  
write_delim(segs_nest,file = file.path(directory,"segs_nest.txt"))
load(file=file.path(directory,"segs_nest.Rdata"))


