library(plotly)
library(tidyverse)
library(plyr)

file <- "E:/track_mask/WT MMC/Traj_190312_53bp1_GFP_B2WTG10_MMC_50ms_100_f488int_0001__Ch1_preprocessed_tracks_mask.csv"
load("E:/track_mask/WT MMC/msd_fit.Rdata")
#load("E:/track_mask/WT MMC/track_stats.Rdata")
tracks_cell$ <- read.csv(file,sep = ",",header = T)
dats <- subset(ldply(tracks),cellID==basename(file))

names(tracks_cell) <- c("track", "pos_x", "pos_y","pos_z", "time", "frame", "step_x", "step_y", "step_z", "inMask")

#tracks <- tracks[,c(3,4,5,2,6,7,8,9,10,11,12)]
head(tracks_cell)

tracks <- subset(segments_all$`WT MMC`,cellID==segments_all$`WT MMC`$cellID[1])


pixelsize <- c(120,120,410) #nm

tracks$X <- tracks$X*pixelsize[1]/1000
tracks$Y <- tracks$Y*pixelsize[2]/1000
tracks$Z <- tracks$Z*pixelsize[3]/1000

tracks$logD <- log10(tracks$D)


minD <- min(tracks$logD)
maxD <- max(tracks$logD)

tracks2 <- subset(tracks,inMask.x==TRUE)
trackIDs <- unique(tracks2$track)


p=plot_ly(data=subset(tracks2,track==trackIDs[1]), x = ~X, y = ~Y, z = ~Z, type = 'scatter3d', mode = 'lines',
          line = list(color = ~logD, width = 1,colorscale="Rainbow",cmin=minD,cmax=maxD),showlegend=FALSE)
for(i in 2:length(trackIDs)){
  p <- add_trace(p,data=subset(tracks2,track==trackIDs[i]),x = ~X, y = ~Y, z = ~Z,
                 line = list(color = ~logD, width = 1,colorscale="Rainbow",cmin=minD,cmax=maxD))
}

p %>%
  layout(scene = list(aspectmode = "manual", aspectratio = list(x=1, y=1, z=0.2)))


%>%
  layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
