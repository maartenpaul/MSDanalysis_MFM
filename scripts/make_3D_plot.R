<<<<<<< HEAD
library(plotly)
library(tidyverse)
library(dplyr)
library(plyr)


file <- "/media/DATA/Maarten/data/WT MMC/Traj_190312_53bp1_GFP_B2WTG10_MMC_50ms_100_f488int_0001__Ch1_preprocessed_tracks_mask.csv"
file <- "/media/DATA/Maarten/data/WT MMC/Traj_190313exp1_53bp1_GFP_B2WTG10_MMC_50ms_100_f488int_0003__Ch1_preprocessed_tracks_mask.csv"
load("/media/DATA/Maarten/data/msd_fit_all.Rdata")

#load("E:/track_mask/WT MMC/track_stats.Rdata")
tracks_cell <- read.csv(file,sep = ",",header = T)

#names(tracks_cell) <- c("track", "pos_x", "pos_y","pos_z", "time", "frame", "step_x", "step_y", "step_z", "inMask")

#tracks <- tracks[,c(3,4,5,2,6,7,8,9,10,11,12)]
head(tracks_cell)

trackstats <- subset(msd_fit_all$`WT MMC`,.id==basename(file))
trackstats$trackID <- trackstats$track

tracks <- inner_join(trackstats,tracks_cell,by="trackID")

pixelsize <- c(120,120,410) #nm

tracks$X <- tracks$pos_x*pixelsize[1]/1000
tracks$Y <- tracks$pos_y*pixelsize[2]/1000
tracks$Z <- tracks$pos_z*pixelsize[3]/1000

tracks$logD <- log10(tracks$D)


minD <- min(tracks$logD,na.rm = T)
maxD <- max(tracks$logD,na.rm = T)

tracks2 <- subset(tracks,inMask.x==TRUE)
tracks2 <- tracks

trackIDs <- unique(tracks2$track)

tracks2$Z <- tracks2$Z*10

p=plot_ly(data=subset(tracks2,track==trackIDs[1]), x = ~X, y = ~Y, z = ~Z, type = 'scatter3d', mode = 'lines',
          line = list(color = ~logD, width = 1,colorscale="Rainbow",cmin=minD,cmax=maxD),showlegend=FALSE)
for(i in 2:length(trackIDs)){
  p <- add_trace(p,data=subset(tracks2,track==trackIDs[i]),x = ~X, y = ~Y, z = ~Z,
                 line = list(color = ~logD, width = 1,colorscale="Rainbow",cmin=minD,cmax=maxD))
}

p %>%
  layout(scene = list(aspectmode = "manual", aspectratio = list(x=1, y=1, z=0.2)))


tracks2$col <- 0
tracks2$col[tracks2$inMask.x] <- 1
p=plot_ly(data=subset(tracks2,track==trackIDs[1]), x = ~X, y = ~Y, z = ~Z, type = 'scatter3d', mode = 'lines',
          line = list(color = ~col, width = 1,colorscale="Rainbow",cmin=0,cmax=4),showlegend=FALSE)
for(i in 2:length(trackIDs)){
  p <- add_trace(p,data=subset(tracks2,track==trackIDs[i]),x = ~X, y = ~Y, z = ~Z,
                 line = list(color = ~col, width = 1,colorscale="Rainbow",cmin=0,cmax=4))
}

p %>%
  layout(scene = list(aspectmode = "manual", aspectratio = list(x=1, y=1, z=0.2)))
minZ <- min(tracks2$Z)
maxZ <- max(tracks2$Z)

p=plot_ly(data=subset(tracks2,track==trackIDs[1]), x = ~X, y = ~Y, z = ~Z, type = 'scatter3d', mode = 'lines',
          line = list(color = ~Z, width = 1,colorscale="Rainbow",cmin=minZ,cmax=maxZ),showlegend=FALSE)
for(i in 2:length(trackIDs)){
  p <- add_trace(p,data=subset(tracks2,track==trackIDs[i]),x = ~X, y = ~Y, z = ~Z,
                 line = list(color = ~Z, width = 1,colorscale="Rainbow",cmin=minZ,cmax=maxZ))
}

p %>%
  layout(scene = list(aspectmode = "manual", aspectratio = list(x=1, y=1, z=0.2)))
=======
library(plotly)
library(tidyverse)
library(dplyr)
library(plyr)


file <- "/media/DATA/Maarten/data/WT MMC/Traj_190312_53bp1_GFP_B2WTG10_MMC_50ms_100_f488int_0001__Ch1_preprocessed_tracks_mask2.csv"
file <- "/media/DATA/Maarten/data/WT MMC/Traj_190313exp1_53bp1_GFP_B2WTG10_MMC_50ms_100_f488int_0003__Ch1_preprocessed_tracks_mask.csv"
load("/media/DATA/Maarten/data/msd_fit_all.Rdata")

#load("E:/track_mask/WT MMC/track_stats.Rdata")
tracks_cell <- read.csv(file,sep = ",",header = T)
tracks_cell <- subset(segments_all$`WT MMC`,.id==basename(file))


#names(tracks_cell) <- c("track", "pos_x", "pos_y","pos_z", "time", "frame", "step_x", "step_y", "step_z", "inMask")

#tracks <- tracks[,c(3,4,5,2,6,7,8,9,10,11,12)]
head(tracks_cell)

trackstats <- subset(msd_fit_all$`WT MMC`,.id==basename(file))

trackstats$trackID <- trackstats$track

tracks <- inner_join(trackstats,tracks_cell,by="track")

out <- tracks %>%
  mutate(X_pix=round(X/20*50),Y_pix=round(Y/20*50)) %>%
  group_by(X_pix,Y_pix) %>%
  summarise(mean_logD=mean(logD),n=n(),fraction=sum(state==2)/n())
library(viridis)
out %>%
ggplot(aes(X_pix, Y_pix, fill= mean_logD)) + 
  geom_tile()+scale_fill_viridis(option = "B")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                         panel.background = element_blank(), axis.line = element_line(colour = "black"))

out %>%
  ggplot(aes(X_pix, Y_pix, fill= fraction)) + 
  geom_tile()+scale_fill_viridis(option = "B")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                     panel.background = element_blank(), axis.line = element_line(colour = "black"))

out %>%
  ggplot(aes(X_pix, Y_pix, fill= n)) + 
  geom_tile()+scale_fill_viridis(option = "B")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                     panel.background = element_blank(), axis.line = element_line(colour = "black"))


pixelsize <- c(120,120,410) #nm

tracks$X <- tracks$pos_x*pixelsize[1]/1000
tracks$Y <- tracks$pos_y*pixelsize[2]/1000
tracks$Z <- tracks$pos_z*pixelsize[3]/1000

tracks$logD <- log10(tracks$D)

tracks




minD <- min(tracks$logD,na.rm = T)
maxD <- max(tracks$logD,na.rm = T)

tracks2 <- subset(tracks,inMask.x==TRUE)
tracks2 <- tracks

trackIDs <- unique(tracks2$track)

tracks2$Z <- tracks2$Z*10

p=plot_ly(data=subset(tracks2,track==trackIDs[1]), x = ~X, y = ~Y, z = ~Z, type = 'scatter3d', mode = 'lines',
          line = list(color = ~logD, width = 1,colorscale="Rainbow",cmin=minD,cmax=maxD),showlegend=FALSE)
for(i in 2:length(trackIDs)){
  p <- add_trace(p,data=subset(tracks2,track==trackIDs[i]),x = ~X, y = ~Y, z = ~Z,
                 line = list(color = ~logD, width = 1,colorscale="Rainbow",cmin=minD,cmax=maxD))
}

p %>%
  layout(scene = list(aspectmode = "manual", aspectratio = list(x=1, y=1, z=0.2)))


tracks2$col <- 0
tracks2$col[tracks2$inMask.x] <- 1
p=plot_ly(data=subset(tracks2,track==trackIDs[1]), x = ~X, y = ~Y, z = ~Z, type = 'scatter3d', mode = 'lines',
          line = list(color = ~col, width = 1,colorscale="Rainbow",cmin=0,cmax=4),showlegend=FALSE)
for(i in 2:length(trackIDs)){
  p <- add_trace(p,data=subset(tracks2,track==trackIDs[i]),x = ~X, y = ~Y, z = ~Z,
                 line = list(color = ~col, width = 1,colorscale="Rainbow",cmin=0,cmax=4))
}

p %>%
  layout(scene = list(aspectmode = "manual", aspectratio = list(x=1, y=1, z=0.2)))
minZ <- min(tracks2$Z)
maxZ <- max(tracks2$Z)

p=plot_ly(data=subset(tracks2,track==trackIDs[1]), x = ~X, y = ~Y, z = ~Z, type = 'scatter3d', mode = 'lines',
          line = list(color = ~Z, width = 1,colorscale="Rainbow",cmin=minZ,cmax=maxZ),showlegend=FALSE)
for(i in 2:length(trackIDs)){
  p <- add_trace(p,data=subset(tracks2,track==trackIDs[i]),x = ~X, y = ~Y, z = ~Z,
                 line = list(color = ~Z, width = 1,colorscale="Rainbow",cmin=minZ,cmax=maxZ))
}

p %>%
  layout(scene = list(aspectmode = "manual", aspectratio = list(x=1, y=1, z=0.2)))
>>>>>>> c2032d864c6deb6bcf4ec51b540031195cc977a3
