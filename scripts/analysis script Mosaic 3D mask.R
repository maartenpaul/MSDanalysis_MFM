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

directory <- "E:/track_mask/"
directory <- "D:/OneDrive/Data2/MFM/track_mask"
condition_list <- list.dirs(directory,full.names = F,recursive = F)
#condition_list <- condition_list[c(1,2)]


# MSD analysis ------------------------------------------------------------


msd_analyze_data_mosaic_mask(directory,condition_list,framerate,n,fitzero,min_length,pixelsize,fitMSD,offset,max_tracks,dim=dim)


load(file.path(directory,"msd_fit_all_  .Rdata"))
#load(file.path(directory,"track_stats_all.Rdata"))
load(file.path(directory,"segments_all_  .Rdata"))

#msd_histogram(msd_fit_all,directory,threshold = 0.05)

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#c00000","#6599d9","#1f497d","#542788","#de77ae","#217d68","#6dc5aa")), ...)

}


scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#c00000","#6599d9","#1f497d","#542788","#de77ae","#217d68","#6dc5aa")), ...)

}

ggplot(data = msd_fit_all$`WT MMC`,aes(x=D,y=..density..,color=inMask,fill=inMask)) +
  geom_density(position="identity",alpha=0.5,adjust = 1) +  geom_histogram(position="identity",alpha=0.5)+
  scale_y_continuous(limits = c(0,0.8),expand = c(0, 0))+  scale_x_log10(limits=c(0.00005,2))+ scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=25)

#for all data sets
for (i in 1:length(condition_list)){
  msdplot <- ggplot(data = msd_fit_all[[i]],aes(x=D,y=..density..,color=inMask,fill=inMask)) +
    geom_density(position="identity",alpha=0.5,adjust = 1) +  geom_histogram(position="identity",alpha=0.5)+
    scale_y_continuous(limits = c(0,0.9),expand = c(0, 0))+  scale_x_log10(limits=c(0.00005,2))+ scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=25)
  ggsave(msdplot,filename = file.path(directory,condition_list[i],"MSD_plot.png"))
}



segments_all[1] <- llply(segments_all[1],function(x) {ddply(x,.variables = c("cellID","track"),function(x){
  if(nrow(x)>5){
    if(any(x$inMask==TRUE)){
      x$trackInMask <-TRUE
    } else {
      x$trackInMask <-FALSE

    }

    if(x$inMask[1]==0){
      currState <- "0"
    } else {
      currState <- "1"
    }
    stateSeq <- ""
    lengthState <- ""
    lengthCurState <- 1
    for(i in 1:(nrow(x)-1)){
      if(x$inMask[i]==x$inMask[i+1]){
        lengthCurState <- lengthCurState+1
      } else {
        lengthState <- paste0(lengthState, ".",lengthCurState)
        lengthCurState <- 1
        stateSeq <- paste0(stateSeq,currState)
        if(x$inMask[i+1]==0){
          currState <- "0"
        } else {
          currState <- "1"
        }
      }
    }
    if(x$inMask[i]==x$inMask[i+1]){
      lengthState <- paste0(lengthState, ".",lengthCurState)
      stateSeq <- paste0(stateSeq,currState)
    } else {
      lengthState <- paste0(lengthState, ".",1)
      stateSeq <- paste0(stateSeq,x$inMask[i+1])
    }
    x$stateLength <- substring(lengthState, 2)
    x$stateSeq <- stateSeq
    return(x)
  }
})})

uniqie


# Further analysis --------------------------------------------------------
#load machine learning scripts,remove py variable to avoid errors

rm(py)
ML_load()
py$pixSize <- 0.12
py$t <- 0.03
source_python("python/SMMsplot.py")


`# Process ML --------------------------------------------------------------


#determine in and outside tracks
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

#calculate MLMSS data from the tracks
segments_all[c(6)] <- llply(segments_all[c(6)],function(x){
  ML_segment_tracks(x,directory=directory)
  }
  )
ML_segment_tracks(segments_all[c(6)],directory=directory)


x <- py_load_object(file.path(directory,"WT MMC","x.pydata"))
y <- py_load_object(file.path(directory,"WT MMC","y.pydata"))
allStates <- py_load_object(file.path(directory,"WT MMC","allStates.pydata"))


out <- SMMsplot(x,y,allStates)
repl_python()

segs <- subset(segments_all$`dCTD MMC`,segments_all$`dCTD MMC`$cellID==segments_all$`dCTD MMC`$cellID[1])

save(segments_all,file=file.path(directory,"segments_all_ML.Rdata"))


#Load segmentsdata with ML data
load(file=file.path(directory,"segments_all_ML.Rdata"))
name <- 'WT MMC'

segments_all[[name]]$displacement <- sqrt(segments_all[[name]]$step_x^2+segments_all[[name]]$step_y^2)
inputdata <- subset(segments_all[[name]],displacement>0.2&angle>0&inMask==T&state<2)


angles <- hist(c(inputdata$angle,abs(360-inputdata$angle)),breaks = seq(0,360,15),plot=F)

angles <- data.frame('mids'=angles$mids,'density'=angles$density)

ggplot(angles,aes(x = mids,y=density,fill=density))+geom_bar(stat='identity',width=16) +ylim(c(0,.01))+scale_x_continuous(breaks=seq(0, 350, 45))+coord_polar(start = 0.5*pi,direction = 1)+
  theme(legend.position = "none",text = element_text(size=15))
angles


# loop to make all plots --------------------------------------------------


for (i in 1:length(condition_list)){
  #SMMsplot all
  #segments_all[[i]] <- ML_segment_tracks(segments_all[[i]],directory=directory)


  x <- dlply(segments_all[[i]],.variables = c("cellID","track"), function(x){

    return(x$X*10)
  })
  names(x) <- NULL



  y <- dlply(segments_all[[i]],.variables = c("cellID","track"), function(x){
    return(x$Y*10)
  })
  names(y) <- NULL

  allStates <- dlply(segments_all[[i]],.variables = c("cellID","track"), function(x){
    return(x$state)
  })
  names(allStates) <- NULL


 # x <- py_load_object(file.path(directory,condition_list[i],"x.pydata"))
#  y <- py_load_object(file.path(directory,condition_list[i],"y.pydata"))
#  allStates <- py_load_object(file.path(directory,condition_list[i],"allStates.pydata"))
  SMMsplot(x,y,allStates)
  py_run_string(paste0("plt.savefig('",file.path(directory,condition_list[i],"SMMsplot_all.png"),"')"))

  list_tracks_inside <- unique(data.frame("cellID"=segments_all[[i]]$cellID,"track" = segments_all[[i]]$track,"trackInMask"=segments_all[[i]]$trackInMask))
  SMMsplot(x[list_tracks_inside$trackInMask],y[list_tracks_inside$trackInMask],allStates[list_tracks_inside$trackInMask])
  py_run_string(paste0("plt.savefig('",file.path(directory,condition_list[i],"SMMsplot_inside.png"),"')"))
  SMMsplot(x[!list_tracks_inside$trackInMask],y[!list_tracks_inside$trackInMask],allStates[!list_tracks_inside$trackInMask])
  py_run_string(paste0("plt.savefig('",file.path(directory,condition_list[i],"SMMsplot_outside.png"),"')"))
  py$plt$close('all')

  #angle plot ALL
  segments_all[[i]]$displacement <- sqrt(segments_all[[i]]$step_x^2+segments_all[[i]]$step_y^2)
  inputdata <- subset(segments_all[[i]],displacement>0.2&angle>0&state<2)
  angles <- hist(c(inputdata$angle,abs(360-inputdata$angle)),breaks = seq(0,360,15),plot=F)
  angles <- data.frame('mids'=angles$mids,'density'=angles$density)

  plot1 <- ggplot(angles,aes(x = mids,y=density,fill=density))+geom_bar(stat='identity',width=16) +ylim(c(0,.01))+scale_x_continuous(breaks=seq(0, 350, 45))+coord_polar(start = 0.5*pi,direction = 1)+
    theme(legend.position = "none",text = element_text(size=15))
  ggsave(plot1,filename = file.path(directory,condition_list[i],"angleplot_all.png"))

  #angle plot inside
  segments_all[[i]]$displacement <- sqrt(segments_all[[i]]$step_x^2+segments_all[[i]]$step_y^2)
  inputdata <- subset(segments_all[[i]],displacement>0.2&angle>0&inMask==T&state<2)
  angles <- hist(c(inputdata$angle,abs(360-inputdata$angle)),breaks = seq(0,360,15),plot=F)
  angles <- data.frame('mids'=angles$mids,'density'=angles$density)

  plot1 <- ggplot(angles,aes(x = mids,y=density,fill=density))+geom_bar(stat='identity',width=16) +ylim(c(0,.01))+scale_x_continuous(breaks=seq(0, 350, 45))+coord_polar(start = 0.5*pi,direction = 1)+
    theme(legend.position = "none",text = element_text(size=15))
  ggsave(plot1,filename = file.path(directory,condition_list[i],"angleplot_inside.png"))

  #angle plot outside
  inputdata <- subset(segments_all[[i]],displacement>0.2&angle>0&inMask==F&state<2)
  angles <- hist(c(inputdata$angle,abs(360-inputdata$angle)),breaks = seq(0,360,15),plot=F)
  angles <- data.frame('mids'=angles$mids,'density'=angles$density)

  plot1 <- ggplot(angles,aes(x = mids,y=density,fill=density))+geom_bar(stat='identity',width=16) +ylim(c(0,.01))+scale_x_continuous(breaks=seq(0, 350, 45))+coord_polar(start = 0.5*pi,direction = 1)+
    theme(legend.position = "none",text = element_text(size=15))
  ggsave(plot1,filename = file.path(directory,condition_list[i],"angleplot_outside.png"))

}

#transition matrices
data <- subset(segments_all$`WT MMC` ,trackInMask==TRUE)
data$trackID <- paste0(data$cellID,"_",data$track)
data_states <- dlply(.data = data,.variables = "trackID",function(x){
  return(x$state)
})
names(data_states) <- NULL
py$allStates <- data_states
