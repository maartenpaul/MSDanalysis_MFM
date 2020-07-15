#required packages
library(plyr)
library(lattice)
library(stats)
library(MSDtracking)
library(ggplot2)
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
pixelsize <- c(1000,1000,1000) #nm
#pixelsize = 1000
fitMSD <- T
offset <- 4*(0.01)^2 #experimentally determined
max_tracks <- 500 #maximum number of tracks per frame else exclude tracks from dataset, avoids mislinking of tracks
dim <- 3 #number of dimensions of tracking

#directory <- "/media/DATA/Maarten/data3/"
directory <- "/media/DATA/Maarten/data_gtv2/"

condition_list <- list.dirs(directory,full.names = F,recursive = F)
#condition_list <- condition_list[c(2,3,5)] #use if only certain conditions should be analyzed, otherwise will take all conditions from folder

#import all data and estimate MSD for all tracks
msd_analyze_data_mosaic_mask_parallel_intensity(directory,condition_list,framerate,n,fitzero,min_length,pixelsize,fitMSD,offset,max_tracks,dim=dim,groundtruth=TRUE)

load(file.path(directory,"msd_fit_all.Rdata"))
load(file.path(directory,"segments_all.Rdata"))

#save segments to text files
for (i in 1:length(msd_fit_all)){
   inner_join(segments_all[[i]],select(msd_fit_all[[i]],c(.id,track,D)),by=c(".id","track"),keep=FALSE) %>%
  write_tsv(file.path(directory,paste0(names(msd_fit_all)[i],"_segmentsD.txt")))
}

#apply machine learning segmentation
segments_all <- llply(segments_all,function(x){
  ddply(x, .variables = "cellID", function(x){
  ML_segment_tracks(x)
  })
})

save(segments_all,file=file.path(directory,"segments_all_ML.Rdata"))
load(file=file.path(directory,"segments_all_ML.Rdata"))

#initialize cluster
nodes <- detectCores()
cl <- makeCluster(nodes-16)
registerDoParallel(cl)

ptm <- proc.time()
#make tracklet column, based on ML results
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

save(segments_all,file=file.path(directory,"segments_all_ML.Rdata"))
load(file=file.path(directory,"segments_all_ML.Rdata"))


###calculate angles
#initialize cluster
nodes <- detectCores()
cl <- makeCluster(nodes-20)
registerDoParallel(cl)

ptm <- proc.time()
#calculate angles and displacements
for (j in 1:length(segments_all)){
  segments_all[[j]] <- ddply(segments_all[[j]],.variables = c("cellID"), .parallel = T, function(x){
    get_angles <- function(x,n){
      x$frame  <- x$frame-x$frame[1]+1
      angles <- rep(x=-1,nrow(x))
      if(nrow(x)>=(3*n+n-1)){
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
      # x$angle6 <- get_angles(x,6)
      # x$angle7 <- get_angles(x,7)
      # x$angle8 <- get_angles(x,8)
      # x$angle9 <- get_angles(x,9)
      # x$angle10 <- get_angles(x,10)
      # x$angle11 <- get_angles(x,11)
      # x$angle12 <- get_angles(x,12)
      # x$angle13 <- get_angles(x,13)
      # x$angle14 <- get_angles(x,14)
      # x$angle15 <- get_angles(x,15)
      # x$angle16 <- get_angles(x,16)
      # x$angle17 <- get_angles(x,17)
      # x$angle18 <- get_angles(x,18)
      # x$angle19 <- get_angles(x,19)
      # x$angle20 <- get_angles(x,20)
      # 
      
      return(x)
      
    })
    
  })
}
stopCluster(cl)
proc.time() - ptm

save(segments_all,file=file.path(directory,"segments_all_angles.Rdata"))
load(file=file.path(directory,"segments_all_angles.Rdata"))


#calculate MSD and MSS from tracklets
source_python('python/getMSDandMSS_R.py')
py$pixSize <- 0.12
py$t <- 0.03

# MSD_only2 <- function(x,y){
#   if(length(x)>10){
#     out <- getMSDandMSS_R(x*10,y*10)
#     return(out[[1]])
#   } else {
#     return(NA)
#   }
# }

MSD_MSS2 <- function(x,y){
  if(length(x)>10){
    out <- getMSDandMSS_R(x*10,y*10)
    return(tibble("MSD_ML"=out[[1]],"MSS_ML"=out[[2]]))
  } else {
    return(tibble("MSD_ML"=NA,"MSS_ML"=NA))
  }
}

#get msd and mss 
ptm <- proc.time()
all <- as_tibble(bind_rows(segments_all,.id = "condition")) %>%
  group_by(condition,.id,tracklet)%>%
  mutate("data"=list(MSD_MSS2(X,Y)))
save(all,file=file.path(directory,"all.rdata"))


all <- all %>%
   unnest_legacy(cols=data) %>%
  ungroup()

save(all,file=file.path(directory,"all.rdata"))
load(file=file.path(directory,"all.rdata"))

segments_all <- all %>%
  group_by(condition) %>%
  group_split() %>%
  setNames(unique(all$condition))

save(segments_all,file=file.path(directory,"segments_nest.Rdata"))
load(file=file.path(directory,"segments_nest.Rdata"))

for (i in 1:length(segments_all)){
  segments_all[[i]]%>%
    select(-data)%>%
    write_tsv(file.path(directory,paste0(names(segments_all)[i],"_segmentsD_ML.txt")))
}

###write to feather format for fast import in Python
library(feather)

for (i in 1:length(segments_all)){
  segments_all[[i]]%>%
    select(-data)%>%
    write_feather(file.path(directory,paste0(names(segments_all)[i],"_segmentsD_ML.feather")))
}

### Plots -------------------------------------------------------------------
theme_Publication <- function(base_size=14, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line.x = element_line(colour="black"),
            axis.line.y = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#ffffff"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.4, "cm"),
            legend.margin = unit(0.2, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#ffffff",fill="#ffffff"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#c00000","#6599d9","#1f497d","#542788","#de77ae","#217d68","#6dc5aa")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#c00000","#6599d9","#1f497d","#542788","#de77ae","#217d68","#6dc5aa")), ...)
  
}

p1 <- ggplot(data = ldply(msd_fit_all),aes(x=D,y=..density..,color=inMask,fill=inMask)) +
  geom_density(position="identity",alpha=0.5,adjust = 1) +  geom_histogram(position="identity",alpha=0.5)+
  scale_y_continuous(limits = c(0,1),expand = c(0, 0))+  scale_x_log10(limits=c(0.00005,2))+ scale_colour_Publication()+scale_fill_Publication()+facet_wrap(.~condition)
p1
+theme_Publication(base_size=25)

save(p1,file = file.path(directory,"all_histograms.pdf"))

ggplot(data = msd_fit_all$`WT untr`,aes(x=D,y=..density..,color=inMask,fill=inMask)) +
  geom_density(position="identity",alpha=0.5,adjust = 1) +  geom_histogram(position="identity",alpha=0.5)+
  scale_y_continuous(limits = c(0,0.8),expand = c(0, 0))+  scale_x_log10(limits=c(0.00005,2))+ scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=25)

ggplot(data = msd_fit_all$var,aes(x=D,y=..density..)) +
  geom_density(position="identity",alpha=0.5,adjust = 1) +  geom_histogram(position="identity",alpha=0.5)+
  scale_y_continuous(limits = c(0,0.8),expand = c(0, 0))+  scale_x_log10(limits=c(0.00005,2))+ scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=25)

ggplot(angles,aes(x = mids,y=density,fill=density))+geom_bar(stat='identity',width=16) +ylim(c(0,.01))+scale_x_continuous(breaks=seq(0, 350, 45))+coord_polar(start = 0.5*pi,direction = 1)+
  theme(legend.position = "none",text = element_text(size=15))

msd_fit_all <- llply(msd_fit_all,function(x){
  x$cellID <- x$.id
  return(x)
})
msd_histogram(msd_fit_all[c(8,5,7,4,2,6,3,1)],directory = directory,merge_var = "cellID",threshold = 0.05)
