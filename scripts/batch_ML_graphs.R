<<<<<<< HEAD
#required packages
library(plyr)
library(lattice)
library(stats)
library(MSDtracking)
library(ggplot2)
library(ggpol)
library(reticulate)
library(tidyverse)
library(reshape2)
#rm(py)
source("python/ML_py.R")
source("scripts/importdata.R")
py$pixSize <- 0.12
py$t <- 0.03
py$directory <- directory
k <- 7

for (k in 1:length(condition_list)){
  #SMMsplot all
  #segments_all[[k]] <- ML_segment_tracks(segments_all[[k]],directory=directory)
  py$x <- NULL
  py$y <- NULL
  py$allStates <- NULL
  py$condition <- condition_list[k]

  x <- dlply(segments_all[[k]],.variables = c("cellID","track"), function(x){

    return(x$X*10)
  })
  names(x) <- NULL


  y <- dlply(segments_all[[k]],.variables = c("cellID","track"), function(x){
    return(x$Y*10)
  })
  names(y) <- NULL
  
  allStates <- dlply(segments_all[[k]],.variables = c("cellID","track"), function(x){
    return(x$state)
  })
  names(allStates) <- NULL
  
  list_tracks_inside <- unique(data.frame("cellID"=segments_all[[k]]$cellID,"track" = segments_all[[k]]$track,"trackInMask"=segments_all[[k]]$trackInMask))
  list_tracks_inside_GT <- unique(data.frame("cellID"=segments_all[[k]]$cellID,"track" = segments_all[[k]]$track,"trackInMask"=segments_all[[k]]$trackInMask_gt1))
  
  inx <- x[list_tracks_inside$trackInMask]
  iny <- y[list_tracks_inside$trackInMask]
  inStates <- allStates[list_tracks_inside$trackInMask]
  
  inGTx <- x[list_tracks_inside_GT$trackInMask]
  inGTy <- y[list_tracks_inside_GT$trackInMask]
  inGTStates <- allStates[list_tracks_inside_GT$trackInMask] 

  outx <- x[!list_tracks_inside$trackInMask]
  outy <- y[!list_tracks_inside$trackInMask]
  outStates <- allStates[!list_tracks_inside$trackInMask]
  py$plt <- NULL
  
  py$x <- x
  py$y <- y
  py$allStates <- allStates
  
  py$state <- 'all'
  py_run_file("python/SMMsplot_as_script_batch.py")
  #py_run_string(paste0("plt.savefig('",file.path(directory,condition_list[k],"SMMsplot_all.png"),"')"))
  #py_run_string(paste0("plt.savefig('",file.path(directory,condition_list[k],"SMMsplot_all.pdf"),"')"))
  py$plt <- NULL
  py$x <- NULL
  py$y <- NULL
  py$allStates <- NULL
 
   #inside
  py$x <- inx
  py$y <- iny
  py$allStates <- inStates
  py$state <- 'inside'
  
  py_run_file("python/SMMsplot_as_script_batch.py")
  #py_run_string(paste0("plt.savefig('",file.path(directory,condition_list[k],"SMMsplot_inside.png"),"')"))
  #py_run_string(paste0("plt.savefig('",file.path(directory,condition_list[k],"SMMsplot_inside.pdf"),"')"))
    
  py$x <- NULL
  py$y <- NULL
  py$allStates <- NULL
  #outside
  py$x <- outx
  py$y <- outy
  py$allStates <- outStates
  py$state <- 'outside'
  py$plt <- NULL
  
  py_run_file("python/SMMsplot_as_script_batch.py")
  #py_run_string(paste0("plt.savefig('",file.path(directory,condition_list[k],"SMMsplot_outside.png"),"')"))
  #py_run_string(paste0("plt.savefig('",file.path(directory,condition_list[k],"SMMsplot_outside.pdf"),"')"))

  py$x <- NULL
  py$y <- NULL
  py$allStates <- NULL
  #inside GT1
  py$x <- inGTx
  py$y <- inGTy
  py$allStates <- inGTStates
  py$state <- 'insideGT1'
  
  py_run_file("python/SMMsplot_as_script_batch.py")
  #py_run_string(paste0("plt.savefig('",file.path(directory,condition_list[k],"SMMsplot_insideGT1.png"),"')"))
  #py_run_string(paste0("plt.savefig('",file.path(directory,condition_list[k],"SMMsplot_insideGT1.pdf"),"')"))
  
  
  #angle plot ALL
  segments_all[[k]]$displacement <- sqrt(segments_all[[k]]$step_x^2+segments_all[[k]]$step_y^2)
  inputdata <- subset(segments_all[[k]],displacement>0.2&angle>0&state<2)
  angles <- hist(c(inputdata$angle,abs(360-inputdata$angle)),breaks = seq(0,360,15),plot=F)
  angles <- data.frame('mids'=angles$mids,'density'=angles$density)

  plot1 <- ggplot(angles,aes(x = mids,y=density,fill=density))+geom_bar(stat='identity',width=16) +ylim(c(0,.01))+scale_x_continuous(breaks=seq(0, 350, 45))+coord_polar(start = 0.5*pi,direction = 1)+
    theme(legend.position = "none",text = element_text(size=15))
  ggsave(plot1,filename = file.path(directory,condition_list[k],"angleplot_all.png"))
  ggsave(plot1,filename = file.path(directory,condition_list[k],"angleplot_all.pdf"))
  

  #angle plot outside
  segments_all[[k]]$displacement <- sqrt(segments_all[[k]]$step_x^2+segments_all[[k]]$step_y^2)
  inputdata <- subset(segments_all[[k]],displacement>0.2&angle>0&inMask==T&state<2)
  angles <- hist(c(inputdata$angle,abs(360-inputdata$angle)),breaks = seq(0,360,15),plot=F)
  angles <- data.frame('mids'=angles$mids,'density'=angles$density)

  plot1 <- ggplot(angles,aes(x = mids,y=density,fill=density))+geom_bar(stat='identity',width=16) +ylim(c(0,.01))+scale_x_continuous(breaks=seq(0, 350, 45))+coord_polar(start = 0.5*pi,direction = 1)+
    theme(legend.position = "none",text = element_text(size=15))
  ggsave(plot1,filename = file.path(directory,condition_list[k],"angleplot_inside.png"))
  ggsave(plot1,filename = file.path(directory,condition_list[k],"angleplot_inside.png"))
  
  #angle plot inside
  segments_all[[k]]$displacement <- sqrt(segments_all[[k]]$step_x^2+segments_all[[k]]$step_y^2)
  inputdata <- subset(segments_all[[k]],displacement>0.2&angle>0&inMask==T&state<2)
  angles <- hist(c(inputdata$angle,abs(360-inputdata$angle)),breaks = seq(0,360,15),plot=F)
  angles <- data.frame('mids'=angles$mids,'density'=angles$density)
  
  plot1 <- ggplot(angles,aes(x = mids,y=density,fill=density))+geom_bar(stat='identity',width=16) +ylim(c(0,.01))+scale_x_continuous(breaks=seq(0, 350, 45))+coord_polar(start = 0.5*pi,direction = 1)+
    theme(legend.position = "none",text = element_text(size=15))
  ggsave(plot1,filename = file.path(directory,condition_list[k],"angleplot_inside.png"))
  ggsave(plot1,filename = file.path(directory,condition_list[k],"angleplot_inside.png"))
  

  #angle plot inside GT1
  inputdata <- subset(segments_all[[k]],displacement>0.2&angle>0&inMask_gt1==T&state<2)
  angles <- hist(c(inputdata$angle,abs(360-inputdata$angle)),breaks = seq(0,360,15),plot=F)
  angles <- data.frame('mids'=angles$mids,'density'=angles$density)

  plot1 <- ggplot(angles,aes(x = mids,y=density,fill=density))+geom_bar(stat='identity',width=16) +ylim(c(0,.01))+scale_x_continuous(breaks=seq(0, 350, 45))+coord_polar(start = 0.5*pi,direction = 1)+
    theme(legend.position = "none",text = element_text(size=15))
  ggsave(plot1,filename = file.path(directory,condition_list[k],"angleplot_inside_GT1.png"))
  ggsave(plot1,filename = file.path(directory,condition_list[k],"angleplot_inside_GT1.pdf"))
  

}

cellstats <- ldply(msd_fit_all,.id="condition", function(x) ddply(x,.variables = ".id", function(x){
    data.frame("fraction"=length(x$inMask[x$inMask])/length(x$inMask),"Ninside"=length(x$inMask[x$inMask]),"fraction_gt1"=length(x$inMask_gt1[x$inMask_gt1])/length(x$inMask_gt1),"Ninside_gt1"=length(x$inMask_gt1[x$inMask_gt1]),"Ntotal"=length(x$inMask_gt1))
}))


ggplot(cellstats,aes(x=condition,y=fraction))+geom_bar()

barplot(cellstats$Ninside)



=======
#required packages
library(plyr)
library(lattice)
library(stats)
library(MSDtracking)
library(ggplot2)
library(ggpol)
library(reticulate)
library(tidyverse)
library(reshape2)
#rm(py)
source("python/ML_py.R")
source("scripts/importdata.R")
py$pixSize <- 0.12
py$t <- 0.03
py$directory <- directory
k <- 7

for (k in 1:length(condition_list)){
  #SMMsplot all
  #segments_all[[k]] <- ML_segment_tracks(segments_all[[k]],directory=directory)
  py$x <- NULL
  py$y <- NULL
  py$allStates <- NULL
  py$condition <- condition_list[k]

  x <- dlply(segments_all[[k]],.variables = c("cellID","track"), function(x){

    return(x$X*10)
  })
  names(x) <- NULL


  y <- dlply(segments_all[[k]],.variables = c("cellID","track"), function(x){
    return(x$Y*10)
  })
  names(y) <- NULL
  
  allStates <- dlply(segments_all[[k]],.variables = c("cellID","track"), function(x){
    return(x$state)
  })
  names(allStates) <- NULL
  
  list_tracks_inside <- unique(data.frame("cellID"=segments_all[[k]]$cellID,"track" = segments_all[[k]]$track,"trackInMask"=segments_all[[k]]$trackInMask))
  list_tracks_inside_GT <- unique(data.frame("cellID"=segments_all[[k]]$cellID,"track" = segments_all[[k]]$track,"trackInMask"=segments_all[[k]]$trackInMask_gt1))
  
  inx <- x[list_tracks_inside$trackInMask]
  iny <- y[list_tracks_inside$trackInMask]
  inStates <- allStates[list_tracks_inside$trackInMask]
  
  inGTx <- x[list_tracks_inside_GT$trackInMask]
  inGTy <- y[list_tracks_inside_GT$trackInMask]
  inGTStates <- allStates[list_tracks_inside_GT$trackInMask] 

  outx <- x[!list_tracks_inside$trackInMask]
  outy <- y[!list_tracks_inside$trackInMask]
  outStates <- allStates[!list_tracks_inside$trackInMask]
  py$plt <- NULL
  
  py$x <- x
  py$y <- y
  py$allStates <- allStates
  
  py$state <- 'all'
  py_run_file("python/SMMsplot_as_script_batch.py")
  #py_run_string(paste0("plt.savefig('",file.path(directory,condition_list[k],"SMMsplot_all.png"),"')"))
  #py_run_string(paste0("plt.savefig('",file.path(directory,condition_list[k],"SMMsplot_all.pdf"),"')"))
  py$plt <- NULL
  py$x <- NULL
  py$y <- NULL
  py$allStates <- NULL
 
   #inside
  py$x <- inx
  py$y <- iny
  py$allStates <- inStates
  py$state <- 'inside'
  
  py_run_file("python/SMMsplot_as_script_batch.py")
  #py_run_string(paste0("plt.savefig('",file.path(directory,condition_list[k],"SMMsplot_inside.png"),"')"))
  #py_run_string(paste0("plt.savefig('",file.path(directory,condition_list[k],"SMMsplot_inside.pdf"),"')"))
    
  py$x <- NULL
  py$y <- NULL
  py$allStates <- NULL
  #outside
  py$x <- outx
  py$y <- outy
  py$allStates <- outStates
  py$state <- 'outside'
  py$plt <- NULL
  
  py_run_file("python/SMMsplot_as_script_batch.py")
  #py_run_string(paste0("plt.savefig('",file.path(directory,condition_list[k],"SMMsplot_outside.png"),"')"))
  #py_run_string(paste0("plt.savefig('",file.path(directory,condition_list[k],"SMMsplot_outside.pdf"),"')"))

  py$x <- NULL
  py$y <- NULL
  py$allStates <- NULL
  #inside GT1
  py$x <- inGTx
  py$y <- inGTy
  py$allStates <- inGTStates
  py$state <- 'insideGT1'
  
  py_run_file("python/SMMsplot_as_script_batch.py")
  #py_run_string(paste0("plt.savefig('",file.path(directory,condition_list[k],"SMMsplot_insideGT1.png"),"')"))
  #py_run_string(paste0("plt.savefig('",file.path(directory,condition_list[k],"SMMsplot_insideGT1.pdf"),"')"))
  
  
  #angle plot ALL
  segments_all[[k]]$displacement <- sqrt(segments_all[[k]]$step_x^2+segments_all[[k]]$step_y^2)
  inputdata <- subset(segments_all[[k]],displacement>0.2&angle>0&state>0)
  angles <- hist(c(inputdata$angle,abs(360-inputdata$angle)),breaks = seq(0,360,15),plot=F)
  angles <- data.frame('mids'=angles$mids,'density'=angles$density)

  plot1 <- ggplot(angles,aes(x = mids,y=density,fill=density))+geom_bar(stat='identity',width=16) +ylim(c(0,.01))+scale_x_continuous(breaks=seq(0, 350, 45))+coord_polar(start = 0.5*pi,direction = 1)+
    theme(legend.position = "none",text = element_text(size=15))
  ggsave(plot1,filename = file.path(directory,condition_list[k],"angleplot_all.png"))
  ggsave(plot1,filename = file.path(directory,condition_list[k],"angleplot_all.pdf"))
  

  #angle plot outside
  segments_all[[k]]$displacement <- sqrt(segments_all[[k]]$step_x^2+segments_all[[k]]$step_y^2)
  inputdata <- subset(segments_all[[k]],displacement>0.2&angle>0&inMask==F&state>0)
  angles <- hist(c(inputdata$angle,abs(360-inputdata$angle)),breaks = seq(0,360,15),plot=F)
  angles <- data.frame('mids'=angles$mids,'density'=angles$density)

  plot1 <- ggplot(angles,aes(x = mids,y=density,fill=density))+geom_bar(stat='identity',width=16) +ylim(c(0,.01))+scale_x_continuous(breaks=seq(0, 350, 45))+coord_polar(start = 0.5*pi,direction = 1)+
    theme(legend.position = "none",text = element_text(size=15))
  ggsave(plot1,filename = file.path(directory,condition_list[k],"angleplot_inside.png"))
  ggsave(plot1,filename = file.path(directory,condition_list[k],"angleplot_inside.png"))
  
  #angle plot inside
  segments_all[[k]]$displacement <- sqrt(segments_all[[k]]$step_x^2+segments_all[[k]]$step_y^2)
  inputdata <- subset(segments_all[[k]],displacement>0.2&angle>0&inMask==T&state>0)
  angles <- hist(c(inputdata$angle,abs(360-inputdata$angle)),breaks = seq(0,360,15),plot=F)
  angles <- data.frame('mids'=angles$mids,'density'=angles$density)
  
  plot1 <- ggplot(angles,aes(x = mids,y=density,fill=density))+geom_bar(stat='identity',width=16) +ylim(c(0,.01))+scale_x_continuous(breaks=seq(0, 350, 45))+coord_polar(start = 0.5*pi,direction = 1)+
    theme(legend.position = "none",text = element_text(size=15))
  plot1
  ggsave(plot1,filename = file.path(directory,condition_list[k],"angleplot_inside.png"))
  ggsave(plot1,filename = file.path(directory,condition_list[k],"angleplot_inside.png"))
  

  #angle plot inside GT1
  inputdata <- subset(segments_all[[k]],displacement>0.2&angle>0&inMask_gt1==T&state<2)
  angles <- hist(c(inputdata$angle,abs(360-inputdata$angle)),breaks = seq(0,360,15),plot=F)
  angles <- data.frame('mids'=angles$mids,'density'=angles$density)

  plot1 <- ggplot(angles,aes(x = mids,y=density,fill=density))+geom_bar(stat='identity',width=16) +ylim(c(0,.01))+scale_x_continuous(breaks=seq(0, 350, 45))+coord_polar(start = 0.5*pi,direction = 1)+
    theme(legend.position = "none",text = element_text(size=15))
  ggsave(plot1,filename = file.path(directory,condition_list[k],"angleplot_inside_GT1.png"))
  ggsave(plot1,filename = file.path(directory,condition_list[k],"angleplot_inside_GT1.pdf"))
  

}

cellstats <- ldply(msd_fit_all,.id="condition", function(x) ddply(x,.variables = ".id", function(x){
    data.frame("fraction"=length(x$inMask[x$inMask])/length(x$inMask),"Ninside"=length(x$inMask[x$inMask]),"fraction_gt1"=length(x$inMask_gt1[x$inMask_gt1])/length(x$inMask_gt1),"Ninside_gt1"=length(x$inMask_gt1[x$inMask_gt1]),"Ntotal"=length(x$inMask_gt1))
}))




ggplot(cellstats,aes(x=condition,y=fraction))+geom_bar()

barplot(cellstats$Ninside)



>>>>>>> c2032d864c6deb6bcf4ec51b540031195cc977a3
