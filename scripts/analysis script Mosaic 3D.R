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
fitzero <- FALSE #should fit go through origin (0,0)

min_length <- 6 #minimum length track
pixelsize <- c(120,120,412) #nm
fitMSD <- T
offset <- 0 #experimentally determined
max_tracks <- 500 #maximum number of tracks per frame else exclude tracks from dataset, avoids mislinking of tracks
dim <- 3 #number of dimensions of tracking

directory <- "D:/Maarten/MFMdata/Tracking/"
condition_list <- list.dirs(directory,full.names = F,recursive = F)
#condition_list <- condition_list[c(1,2)]


msd_analyze_data_mosaic(directory,condition_list[c(2,4,6)],framerate,n,fitzero,min_length,pixelsize,fitMSD,offset,max_tracks,dim=dim)

load(file.path(directory,"msd_fit_all_  .Rdata"))
#load(file.path(directory,"track_stats_all.Rdata"))
load(file.path(directory,"segments_all.Rdata"))

msd_histogram(msd_fit_all,directory,threshold = 0.05)

