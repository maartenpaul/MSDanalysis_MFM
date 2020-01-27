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
rm(py)
source("scripts/ML_py.R")
ML_load()
py$pixSize <- 0.12
py$t <- 0.03
source_python("python/SMMsplot.py")

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
  py_run_string(paste0("plt.savefig('",file.path(directory,condition_list[i],"SMMsplot_inside.png"),"')"))


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
