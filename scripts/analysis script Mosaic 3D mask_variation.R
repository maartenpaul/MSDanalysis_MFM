<<<<<<< HEAD
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

condition_list <- list.dirs(directory,full.names = F,recursive = F)
#condition_list <- condition_list[c(2,3,5)]

msd_analyze_data_mosaic_mask_parallel_intensity(directory,condition_list,framerate,n,fitzero,min_length,pixelsize,fitMSD,offset,max_tracks,dim=dim,groundtruth=TRUE)

load(file.path(directory,"msd_fit_all.Rdata"))
load(file.path(directory,"segments_all.Rdata"))

#save segments to text files
for (i in 1:length(msd_fit_all)){
   inner_join(segments_all[[i]],select(msd_fit_all[[i]],c(.id,track,D)),by=c(".id","track"),keep=FALSE) %>%
  write_tsv(file.path(directory,paste0(names(msd_fit_all)[i],"_segmentsD.txt")))
}

#msd_histogram(msd_fit_all,directory,threshold = 0.05)



#filter <- llply(msd_fit_all[c(6,4,2,5,3,1)],function(x) return(x[x$inMask==T]))

#msd_histogram(llply(msd_fit_all[c(6,4,2,5,3,1)],function(x) return(x[x$inMask==F,])),directory = directory,merge_var = "cellID",threshold = 0.05)


# nodes <- detectCores()
# cl <- makeCluster(nodes-6)
# registerDoParallel(cl)
# segments_all <- llply(segments_all,.parallel = TRUE,function(x) {ddply(x,.variables = c("cellID","track"),function(x){
#   if(nrow(x)>5){
#     if(any(x$inMask==TRUE)){
#       x$trackInMask <-TRUE
#     } else {
#       x$trackInMask <-FALSE
# 
#     }
#     return(x)
#   }
# })})
# stopCluster(cl)



segments_all <- llply(segments_all,function(x){
  ddply(x, .variables = "cellID", function(x){
  ML_segment_tracks(x)
  })
})

save(segments_all,file=file.path(directory,"segments_all_ML.Rdata"))
load(file=file.path(directory,"segments_all_ML.Rdata"))

i <- 1
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


indices <- daply(segments_all[[i]],.variables = c("cellID"), function(x){
  return(length(unique(x$track)))
})
names(indices) <- NULL

indices <- cumsum(indices)


py$x <- x
py$y <- y 
py$allStates <- allStates
py$indices <- indices
py$t <- 0.03
numPmsd<- 4
numPmss <- 4
p <- as.array(seq(from=0.5,to = 6,length.out = 12))


source_python("python/make_table.py")

segments_all[[1]] <- ddply(segments_all[[1]],.variables = c("cellID","track"), function(x) {
  if (nrow(x)>4) {
    out <- getMSDandMSS(x$X*10,x$Y*10,numPmsd,numPmss,p )
    
  }
})

source_python("python/SMMsplot_as_script.py")
py_run_string(paste0("plt.savefig('",file.path(directory,condition_list[i],"SMMsplot_all.png"),"')"))






for (i in 1:length(segments_all)){
  write_tsv(segments_all[[i]],file.path(directory,paste0(names(segments_all)[i],"_segments.txt")))
}



segments_all2 <- segments_all

name <- 'WT MMC'

segments_all[[name]]$displacement <- sqrt(segments_all[[name]]$step_x^2+segments_all[[name]]$step_y^2)
inputdata <- subset(segments_all[[name]],displacement>0.2&angle>0&inMask==F&state<2)


angles <- hist(c(inputdata$angle,abs(360-inputdata$angle)),breaks = seq(0,360,15),plot=F)

angles <- data.frame('mids'=angles$mids,'density'=angles$density)


# Plots -------------------------------------------------------------------
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
=======
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

condition_list <- list.dirs(directory,full.names = F,recursive = F)
#condition_list <- condition_list[c(2,3,5)]

msd_analyze_data_mosaic_mask_parallel_intensity(directory,condition_list,framerate,n,fitzero,min_length,pixelsize,fitMSD,offset,max_tracks,dim=dim,groundtruth=TRUE)

load(file.path(directory,"msd_fit_all.Rdata"))
load(file.path(directory,"segments_all.Rdata"))

#save segments to text files
for (i in 1:length(msd_fit_all)){
   inner_join(segments_all[[i]],select(msd_fit_all[[i]],c(.id,track,D)),by=c(".id","track"),keep=FALSE) %>%
  write_tsv(file.path(directory,paste0(names(msd_fit_all)[i],"_segmentsD.txt")))
}

#msd_histogram(msd_fit_all,directory,threshold = 0.05)



#filter <- llply(msd_fit_all[c(6,4,2,5,3,1)],function(x) return(x[x$inMask==T]))

#msd_histogram(llply(msd_fit_all[c(6,4,2,5,3,1)],function(x) return(x[x$inMask==F,])),directory = directory,merge_var = "cellID",threshold = 0.05)


# nodes <- detectCores()
# cl <- makeCluster(nodes-6)
# registerDoParallel(cl)
# segments_all <- llply(segments_all,.parallel = TRUE,function(x) {ddply(x,.variables = c("cellID","track"),function(x){
#   if(nrow(x)>5){
#     if(any(x$inMask==TRUE)){
#       x$trackInMask <-TRUE
#     } else {
#       x$trackInMask <-FALSE
# 
#     }
#     return(x)
#   }
# })})
# stopCluster(cl)



segments_all <- llply(segments_all,function(x){
  ddply(x, .variables = "cellID", function(x){
  ML_segment_tracks(x)
  })
})

save(segments_all,file=file.path(directory,"segments_all_ML.Rdata"))
load(file=file.path(directory,"segments_all_ML.Rdata"))

i <- 1
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


indices <- daply(segments_all[[i]],.variables = c("cellID"), function(x){
  return(length(unique(x$track)))
})
names(indices) <- NULL

indices <- cumsum(indices)


py$x <- x
py$y <- y 
py$allStates <- allStates
py$indices <- indices
py$t <- 0.03
numPmsd<- 4
numPmss <- 4
p <- as.array(seq(from=0.5,to = 6,length.out = 12))


source_python("python/make_table.py")

segments_all[[1]] <- ddply(segments_all[[1]],.variables = c("cellID","track"), function(x) {
  if (nrow(x)>4) {
    out <- getMSDandMSS(x$X*10,x$Y*10,numPmsd,numPmss,p )
    
  }
})

source_python("python/SMMsplot_as_script.py")
py_run_string(paste0("plt.savefig('",file.path(directory,condition_list[i],"SMMsplot_all.png"),"')"))






for (i in 1:length(segments_all)){
  write_tsv(segments_all[[i]],file.path(directory,paste0(names(segments_all)[i],"_segments.txt")))
}



segments_all2 <- segments_all

name <- 'WT MMC'

segments_all[[name]]$displacement <- sqrt(segments_all[[name]]$step_x^2+segments_all[[name]]$step_y^2)
inputdata <- subset(segments_all[[name]],displacement>0.2&angle>0&inMask==F&state<2)


angles <- hist(c(inputdata$angle,abs(360-inputdata$angle)),breaks = seq(0,360,15),plot=F)

angles <- data.frame('mids'=angles$mids,'density'=angles$density)


# Plots -------------------------------------------------------------------
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
>>>>>>> c2032d864c6deb6bcf4ec51b540031195cc977a3
