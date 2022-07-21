#this script skips the import of data and can be used for subsequent steps of analysis

#required packages
library(plyr)
library(lattice)
library(stats)
library(MSDtracking)
library(ggplot2)#required packages
library(ggpol)
library(doParallel)
library(reticulate)
source('R/MSD.R')
source('R/MSD_fit.R')
source('R/analysis functions.R')
source('python/ML_py.R')

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
directory <- "/media/DATA/Maarten/data3/"
directory <- "/media/DATA/Maarten/data_gtv2/"

condition_list <- list.dirs(directory,full.names = F,recursive = F)
#condition_list <- condition_list[c(2,3,5)]

load(file.path(directory,"msd_fit_all.Rdata"))
#load(file.path(directory,"segments_all.Rdata"))
load(file=file.path(directory,"segments_all_angles.Rdata"))

ptm <- proc.time()
#initialize cluster
nodes <- detectCores()
cl <- makeCluster(nodes-12)
registerDoParallel(cl)

####make tracklet column
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

numPmsd <- 4
numPmss <- 4
minLen <- 10

p <- seq(from=0.5,to=6,length.out=12)

for (j in 1:length(segments_all)){
  
  segments_all[[j]] <- ddply(segments_all[[j]],.variables = c("cellID"),.parallel = T,function(x){
    ddply(x,.variables = c("track"), .parallel = F, function(x){
      if(nrow(x)>minLen){
        getMSDandMSS(x$X,x$Y,numPmsd,numPmss,p)
        
      }
      
      
      return(x)
      
    })})
}

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

save(segments_all,file=file.path(directory,"segments_all_angles.Rdata"))
load(file=file.path(directory,"segments_all_angles.Rdata"))

###get msd and mss from tracklets
numPmsd <- 4
numPmss <- 4
minLen <- 10
p <- seq(from=0.5,to=6,length.out=12)
py$pixSize <- 0.100
py$t <- 0.032
source_python('python/getMSDandMSS_R.py')

MSD_MSS <- function(x){
  if(nrow(x)>10){
    out <- getMSDandMSS_R(x$X,x$Y)
   return(tibble("D_ML"=out[[1]],"D_Smmss"=out[[2]]))
  } else {
  return(tibble("D_ML"=-1.0,"D_Smmss"=-1.0))
    return(NA)
  }
}

MSD_MSS_focus <- function(x){
  if(nrow(x)>10){
    out <- getMSDandMSS_R(x$X,x$Y)
    return(tibble("D_ML_focus"=out[[1]],"D_Smmss_focus"=out[[2]]))
  } else {
    return(tibble("D_ML_focus"=-1.0,"D_Smmss_focus"=-1.0))
    return(NA)
  }
}

MSD_only <- function(x){
  if(nrow(x)>10){
    out <- getMSDandMSS_R(x$X,x$Y)
    return(out[[1]])
  } else {
    return(NA)
  }
}

library(tidyverse)
segments <- ldply(segments_all)
segs_nest <-segments%>%
  select(condition,cellID,focus_tracklet,X,Y) %>%
  group_by(condition,cellID,focus_tracklet) %>%
  group_modify(~MSD_MSS_focus(.x)) %>%
  inner_join(y=segments,by=c("condition","cellID","focus_tracklet")) %>%
  ungroup()


segs_nest <- segs_nest %>%
  select(condition,cellID,tracklet,X,Y) %>%
  group_by(condition,cellID,tracklet) %>%
  group_modify(~MSD_MSS(.x),keep=T) %>%
  inner_join(y=segs_nest,by=c("condition","cellID","tracklet")) 
  
save(segs_nest,file=file.path(directory,"segs_nest.Rdata"))
load(file=file.path(directory,"segs_nest.Rdata"))

###make D histograms
k <- segs_nest %>%
  filter(D_ML_focus>0)%>%
  dplyr::distinct(condition,cellID,tracklet,.keep_all=T)%>%
  mutate(state_str=as.character(inMask))

x <- group_by(k,condition) %>%
  group_by(condition,inMask)%>%
  dplyr::summarise(number=n())%>%
  group_by(condition) %>%
  dplyr::mutate(fraction=number/sum(number))

y <- group_by(k,condition) %>%
  group_by(condition,cellID,inMask)%>%
  dplyr::summarise(number=n())%>%
  group_by(condition,cellID) %>%
  dplyr::mutate(fraction=number/sum(number)) %>%
  group_by(condition,inMask) %>%
  dplyr::summarise(mean_fraction=round(mean(fraction),digits = 2),sd_fraction=round(sd(fraction),digits = 2))


p <- segs_nest %>%
  filter(D_ML_focus>0)%>%
  dplyr::distinct(condition,cellID,tracklet,.keep_all=T)%>%
  mutate(state_str=as.character(inMask))%>%
  ggplot(aes(x=D_ML_focus*100,fill=state_str,y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]))+geom_histogram(position="identity",alpha=0.5)+scale_x_log10(limits=c(0.0001,10))+facet_wrap(.~condition,ncol=2)+
  scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=12)+ theme(legend.position = "none")+ylab("relative frequency")+ xlab(expression(D[app]~mu~m^{2}/s)) +
  geom_text(data=subset(y,inMask==0 ), aes(x=1., y=0.12, label=mean_fraction), colour="#c00000", inherit.aes=FALSE, parse=FALSE)+
  geom_text(data=subset(y,inMask==0 ), aes(x=1., y=0.10, label=paste0("+/-",sd_fraction)), colour="#c00000", inherit.aes=FALSE, parse=FALSE)+
  
  geom_text(data=subset(y,inMask==1 ), aes(x=0.08, y=0.06, label=mean_fraction), colour="#fdae61", inherit.aes=FALSE, parse=FALSE)+
  geom_text(data=subset(y,inMask==1), aes(x=.08, y=0.04, label=paste0("+/-",sd_fraction)), colour="#fdae61", inherit.aes=FALSE, parse=FALSE)
  #geom_text(data=subset(y,state==2 ), aes(x=0.003, y=0.09, label=mean_fraction), colour="#1f497d", inherit.aes=FALSE, parse=FALSE)+
  #geom_text(data=subset(y,state==2 ), aes(x=.003, y=0.07, label=paste0("+/-",sd_fraction)), colour="#1f497d", inherit.aes=FALSE, parse=FALSE)

p


# head(segs_nest)
# flatten(segs_nest[2:4])
# 
# segs_nest2 <- segs_nest %>%
#   map_dfr(MSD_MSS(data))

##use group_by
segs_nest <- segments_all[[1]][1:1000,] %>%
  select(cellID,tracklet,X,Y) %>%
  group_by(cellID,tracklet) %>%
  dplyr::summarise(MSD_MSS(.))


###write to feather format for fast import in Python
library(feather)

for (i in condition_list){
  segments%>%
    filter(condition==i)%>%
    write_feather(file.path(directory,paste0(i,"_segmentsD_ML.feather")))
}
