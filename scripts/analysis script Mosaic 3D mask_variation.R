#required packages
library(plyr)
library(lattice)
library(stats)
library(MSDtracking)
library(ggplot2)#required packages
library(plyr)
library(lattice)
library(stats)
library(MSDtracking)
library(ggplot2)
library(ggpol)
library(doParallel)

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

directory <- "D:/Maarten/variation check/"
directory <- "/media/DATA/data/"

condition_list <- list.dirs(directory,full.names = F,recursive = F)
condition_list <- condition_list[c(1,2)]


msd_analyze_data_mosaic_mask_parallel(directory,condition_list,framerate,n,fitzero,min_length,pixelsize,fitMSD,offset,max_tracks,dim=dim,groundtruth=TRUE)




load(file.path(directory,"msd_fit_all.Rdata"))
#load(file.path(directory,"track_stats_all.Rdata"))
load(file.path(directory,"segments_all.Rdata"))

#msd_histogram(msd_fit_all,directory,threshold = 0.05)

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

+theme_Publication(base_size=25)

save(p1,file = file.path(directory,"all_histograms.pdf"))

ggplot(data = msd_fit_all$var,aes(x=D,y=..density..,color=inMask_gt1,fill=inMask_gt1)) +
  geom_density(position="identity",alpha=0.5,adjust = 1) +  geom_histogram(position="identity",alpha=0.5)+
  scale_y_continuous(limits = c(0,0.8),expand = c(0, 0))+  scale_x_log10(limits=c(0.00005,2))+ scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=25)

ggplot(data = msd_fit_all$var,aes(x=D,y=..density..)) +
  geom_density(position="identity",alpha=0.5,adjust = 1) +  geom_histogram(position="identity",alpha=0.5)+
  scale_y_continuous(limits = c(0,0.8),expand = c(0, 0))+  scale_x_log10(limits=c(0.00005,2))+ scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=25)

p1 <- 


msd_fit_all <- llply(msd_fit_all,function(x){
  x$cellID <- x$.id
  return(x)
})
msd_histogram(msd_fit_all[c(8,5,7,4,2,6,3,1)],directory = directory,merge_var = "cellID",threshold = 0.05)

filter <- llply(msd_fit_all[c(6,4,2,5,3,1)],function(x) return(x[x$inMask==T]))

msd_histogram(llply(msd_fit_all[c(6,4,2,5,3,1)],function(x) return(x[x$inMask==F,])),directory = directory,merge_var = "cellID",threshold = 0.05)


#angle
library(reticulate)

ML_load()



segments_all <- llply(segments_all,function(x) {ddply(x,.variables = c("cellID","track"),function(x){
  if(nrow(x)>5){
    if(any(x$inMask==TRUE)){
      x$trackInMask <-TRUE
    } else {
      x$trackInMask <-FALSE

    }
    return(x)
  }
})})



#ML_segment_tracks(subset(segments_all$`dCTD HU`,cellID=="Traj_190316exp3_53bp1_GFP_B2dCTDA2_HU_50ms_100_f488int_0001__Ch1_preprocessed_tracks_mask.csv"))

segments_all <- llply(segments_all,function(x){
  ML_segment_tracks(x,directory=directory)
}
)

segs <- subset(segments_all$`dCTD MMC`,segments_all$`dCTD MMC`$cellID==segments_all$`dCTD MMC`$cellID[1])

save(segments_all,file=file.path(directory,"segments_all_ML.Rdata"))

name <- 'WT MMC'

segments_all[[name]]$displacement <- sqrt(segments_all[[name]]$step_x^2+segments_all[[name]]$step_y^2)
inputdata <- subset(segments_all[[name]],displacement>0.2&angle>0&inMask==F&state<2)


angles <- hist(c(inputdata$angle,abs(360-inputdata$angle)),breaks = seq(0,360,15),plot=F)

angles <- data.frame('mids'=angles$mids,'density'=angles$density)

ggplot(angles,aes(x = mids,y=density,fill=density))+geom_bar(stat='identity',width=16) +ylim(c(0,.01))+scale_x_continuous(breaks=seq(0, 350, 45))+coord_polar(start = 0.5*pi,direction = 1)+
  theme(legend.position = "none",text = element_text(size=15))
angles
