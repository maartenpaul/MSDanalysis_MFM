#required packages
library(reticulate)
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
use_condaenv("r-reticulate")
source('python/ML_py.R')

#input variables
framerate <- 1/52 #1/200 #1/ms
n <- 4 #number of timepoints taken into account for MSD fit
fitzero <- TRUE #should fit go through origin (0,0)
pixelsize = 1000
fitMSD <- F
max_tracks <- 500 #maximum number of tracks per frame else exclude tracks from dataset, avoids mislinking of tracks
dim <- 2 #number of dimensions of tracking

directory <- "/media/DATA/Maarten/MFM/data_2023/"
condition_list <- list.dirs(directory,full.names = F,recursive = F)

#import data organize in list of the different data sets
msd_analyze_data_mosaic_mask_parallel_intensity(directory,condition_list[5],framerate,n,fitzero,min_length,pixelsize,fitMSD,offset,max_tracks,dim=dim,groundtruth=TRUE)

load(file.path(directory,"msd_fit_all.Rdata"))
load(file.path(directory,"segments_all.Rdata"))

#estimate states using MLMSS
segments_all <- llply(segments_all,function(x){
  ddply(x, .variables = "cellID", function(x){
    ML_segment_tracks(x)
  })
})

save(segments_all,file=file.path(directory,"segments_all_MLSS.Rdata"))
load(file=file.path(directory,"segments_all_MLSS.Rdata"))

#calculate angle displacements
ptm <- proc.time()
#initialize cluster
nodes <- detectCores()
cl <- makeCluster(nodes-15)
registerDoParallel(cl)
for (j in 1:length(segments_all)){
  segments_all[[j]] <- ddply(segments_all[[j]],.variables = c("cellID"), .parallel = T, function(x){
    get_angle <- function(A,B,C,dim=2){
      if (dim==2){
        AB <- B[1:2]-A[1:2]
        CB <- C[1:2]-B[1:2]      
      }else if (dim==3){
        AB <- B[1:3]-A[1:3]
        CB <- C[1:3]-B[1:3]
      }
      #dAB <- sqrt((B[1]-A[1])^2+(B[2]-A[2])^2)
      #dBC <- sqrt((C[1]-B[1])^2+(C[2]-B[2])^2)
      #Formula obtained from https://gitlab.com/anders.sejr.hansen/anisotropy
     # angle <- abs(atan2(det(cbind(AB,CB)),AB%*%CB))
      library(pracma)
      angle<- abs(acos((AB%*%CB)/(Norm(AB)* Norm(CB) )))
      angle <- angle/pi*180
      angle
      return(angle)
    } 
    
    get_angles <- function(x,n,dim=2){
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
            angles[which_point1] <- get_angle(x1,x2,x3,dim=dim)
            
          }}
      }
      return(angles)
    }
    get_focus_angles <- function(x,n,gt="",dim=2){
      x$frame  <- x$frame-x$frame[1]+1
      angles <- rep(x=-1,nrow(x))
      if(nrow(x)>=(3*n+n-1)){
        
        #loop over all steps of tracks
        for (k in 1:(nrow(x)-2*n)){
        
          if (is.element(x[k,"frame"]+n,x$frame)&&x[1,paste0("center_x",gt)]!=-1){ #check if point is present in track
            x1 <- as.numeric(c(x[1,paste0("center_y",gt)],x[1,paste0("center_x",gt)],x[1,paste0("center_z",gt)]))
            x2 <- as.numeric(x[k,c(3,4,6)])
            x2[1:2] <- x2[1:2]+0.96
            which_point1 <- which(x[,"frame"]==x[k,"frame"]+n)
            x3 <- as.numeric(x[which_point1,c(3,4,6)])
            #which_point2 <- which(x[,2]==x[k,2]+2*n)
            #x3 <- as.numeric(x[which_point2,c(3,4,6)])
            x3[1:2] <- x3[1:2]+0.96
            
            angles[k] <- get_angle(x1,x2,x3,dim)
            
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
    get_radial_displacements <- function(x,n,gt=""){
      x$frame  <- x$frame-x$frame[1]+1
      displ <- cbind(rep(x=NA,nrow(x)),rep(x=NA,nrow(x)))
      if(nrow(x)>=(3*n+n-1)){
        #loop over all steps of tracks
        for (k in 1:(nrow(x)-2*n)){
          
          if (is.element(x[k,"frame"]+n,x$frame)){ #check if point is present in track, this deals with gaps
            x1 <- as.numeric(x[,paste0("distanceToNearest3D",gt)][k])
            which_point1 <- which(x[,2]==x[k,2]+n)
            x2 <- as.numeric(x[,paste0("distanceToNearest3D",gt)][which_point1])
            which_point2 <- which(x[,2]==x[k,2]+2*n)
            x3 <- as.numeric(x[,paste0("distanceToNearest3D",gt)][which_point2])
            displ[k,] <- c((x2[1]-x1[1]))
            
          }}
      }
      return(displ)
    }
    
  
    ddply(x,.variables = c("track"), .parallel = F, function(x){
      
      
      x$angle1 <- get_angles(x,1)
      x$angle1_3D <- get_angles(x,1,dim=3)
      x[c("displacement1","displacement2")] <- get_displacements(x,1)
      x$radial_displacement1 <- get_radial_displacements(x,1,"")
      x$radial_displacement1_gt1 <- get_radial_displacements(x,1,"_gt1")
      x$focus_angle <- get_focus_angles(x,1,gt="",dim=2)
      x$focus_angle_3D <- get_focus_angles(x,1,gt="",dim=3)
      x$focus_angle_gt1 <- get_focus_angles(x,1,"_gt1",dim=2)
      x$focus_angle_3D_gt1 <- get_focus_angles(x,1,"_gt1",dim=3)
      
      return(x)
      
    })
    
  })
}
stopCluster(cl)
proc.time() - ptm

save(segments_all,file=file.path(directory,"segments_all_angles.Rdata"))
load(file=file.path(directory,"segments_all_angles.Rdata"))

####add tracklet column
ptm <- proc.time()
#initialize cluster
nodes <- detectCores()
cl <- makeCluster(nodes-10)
registerDoParallel(cl)
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


####split tracks in inside outside column
ptm <- proc.time()
#initialize cluster
nodes <- detectCores()
cl <- makeCluster(nodes-10)
registerDoParallel(cl)

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
stopCluster(cl)
proc.time() - ptm
####split tracks in inside outside gt1 column
ptm <- proc.time()
#initialize cluster
nodes <- detectCores()
cl <- makeCluster(nodes-10)
registerDoParallel(cl)
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
stopCluster(cl)
proc.time() - ptm
save(segments_all,file=file.path(directory,"segments_all_angles_msd.Rdata"))
load(file=file.path(directory,"segments_all_angles_msd.Rdata"))

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
#  filter(condition=="WT MMC")%>%
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

#add tracklet in Mask column
segs_nest <- segs_nest %>%
  group_by(.id,tracklet,condition)%>%
  dplyr::mutate(trackletInMask=any(inMask==T))

#save csv files
segs_nest$condition <- droplevels(segs_nest$condition)

save(segs_nest,file=file.path(directory,"segs_nest.Rdata"))
segs_nest %>%
  nest(-.id) %>%
  pwalk(~write_delim(x = .y, file = file.path(directory,paste0(.x, ".txt") )) )

write_delim(segs_nest,file = file.path(directory,"segs_nest.txt"))
load(file=file.path(directory,"segs_nest.Rdata"))


