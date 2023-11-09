#this script skips the import of data and can be used for subsequent steps of analysis

#required packages
library(plyr)
library(tidyverse)
library(lattice)
library(stats)
library(MSDtracking)
library(ggplot2)#required packages
library(ggpol)
library(doParallel)
library(reticulate)
library(ggbeeswarm)
library(ggsignif)
source('R/MSD.R')
source('R/MSD_fit.R')
source('R/analysis functions.R')

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
  discrete_scale("fill","Publication",manual_pal(values = c("#1F497D","#c00000","#542788","#217D68","#386cb0","#fdb462","#386cb0","#fdb462","#386cb0","#fdb462","#386cb0","#fdb462")), ...)
  
}
scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#1F497D","#c00000","#542788","#217D68","#386cb0","#fdb462","#386cb0","#fdb462","#386cb0","#fdb462","#386cb0","#fdb462")), ...)
  
}

#input variables
directory <- "/media/OIC-station2/MFM/TM/"

#condition_list <- list.dirs(directory,full.names = F,recursive = F)
#condition_list <- condition_list[c(2,3,5)]

load(file=file.path(directory,"segs_nest.Rdata"))

condition_id <- "WT MMC"


# plots -------------------------------------------------------------------


#histogram in-out foci
p <- segs_nest %>%
  filter(D_ML_focus>0,condition==condition_id)%>%
  dplyr::distinct(condition,cellID,focus_tracklet,.keep_all=T)%>%
  mutate(state_str=as.character(inMask))%>%
  ggplot(aes(x=D_ML_focus*100,fill=state_str,y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]))+geom_histogram(position="identity",alpha=0.5)+scale_x_log10(limits=c(0.0001,10))+facet_wrap(.~state_str,ncol=2)+
  scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=8)+ theme(legend.position = "none")+ylab("relative frequency")+ xlab(expression(D[app]~mu~m^{2}/s)) + 
  scale_y_continuous(expand=c(0,0))
  
p
ggsave(p,filename = paste0("/media/DATA/Maarten/OneDrive/Data2/MFM/R plots/",condition_id,"_msd_histogram_inout_mask.pdf"),width = 10,height = 5,units = "cm")
ggsave(p,filename = paste0("/media/DATA/Maarten/OneDrive/Data2/MFM/R plots/",condition_id,"_msd_histogram_inout_mask.png"),width = 10,height = 5,units = "cm")

#histogram in-out foci gt1
p <- segs_nest %>%
  filter(D_ML_focus_gt1>0&condition==condition_id)%>%
  dplyr::distinct(condition,cellID,focus_tracklet_gt1,.keep_all=T)%>%
  mutate(state_str=as.factor(inMask_gt1))%>%
  ggplot(aes(x=D_ML_focus_gt1*100,fill=state_str,y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]))+
  geom_histogram(position="identity",alpha=0.5)+scale_x_log10(limits=c(0.0001,10))+facet_wrap(.~state_str,ncol=2)+
  scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=8)+ theme(legend.position = "none")+ylab("relative frequency")+ xlab(expression(D[app]~mu~m^{2}/s)) 
p
ggsave(p,filename = paste0("/media/DATA/Maarten/OneDrive/Data2/MFM/R plots/",condition_id,"_msd_histogram_inout_mask_gt1.pdf"),width = 10,height = 5,units = "cm",bg = "transparent")
ggsave(p,filename = paste0("/media/DATA/Maarten/OneDrive/Data2/MFM/R plots/",condition_id,"_msd_histogram_inout_mask_gt1.png"),width = 10,height = 5,units = "cm",bg = "transparent")

#plotting of displacements in out
p <- segs_nest %>%
  filter(displacement1>0&condition==condition_id)%>%
  ggplot(aes(x=displacement1*1000,fill=as.character(inMask)))+geom_histogram(bins = 60,aes(y=..density..),alpha=0.5,position='identity')+xlim(0,750)+
  geom_density(alpha=0.5,aes(color=as.character(inMask),fill=NULL))+scale_colour_Publication()+
  scale_fill_Publication()+theme_Publication(base_size=12)+ xlab("2D displacement (nm)")+ theme(legend.position = "none")+scale_y_continuous(expand=c(0,0))+scale_x_continuous(expand=c(0,0),limits=c(0,500))
p
ggsave(p,filename = paste0("/media/DATA/Maarten/OneDrive/Data2/MFM/R plots/",condition_id,"_displacements_inout_mask.pdf"),width = 10,height = 5,units = "cm",bg = "transparent")
ggsave(p,filename = paste0("/media/DATA/Maarten/OneDrive/Data2/MFM/R plots/",condition_id,"_displacements_inout_mask.png"),width = 10,height = 5,units = "cm",bg = "transparent")


p <- segs_nest %>%
  filter(disp2d>0&condition==condition_id)%>%
  ggplot(aes(x=disp2d*1000,fill=as.character(inMask)))+geom_histogram(bins = 60,aes(y=..density..),alpha=0.5,position='identity')+xlim(0,750)+
  geom_density(alpha=0.5,aes(color=as.character(inMask),fill=NULL))+scale_colour_Publication()+
  scale_fill_Publication()+theme_Publication(base_size=12)+ xlab("2D displacement (nm)")+ theme(legend.position = "none")+scale_y_continuous(expand=c(0,0))+scale_x_continuous(expand=c(0,0),limits=c(0,500))
p
ggsave(p,filename = paste0("/media/DATA/Maarten/OneDrive/Data2/MFM/R plots/",condition_id,"_displacements_inout_2d_mask.pdf"),width = 10,height = 5,units = "cm",bg = "transparent")
ggsave(p,filename = paste0("/media/DATA/Maarten/OneDrive/Data2/MFM/R plots/",condition_id,"_displacements_inout_2d_mask.png"),width = 10,height = 5,units = "cm",bg = "transparent")

p <- segs_nest %>%
  filter(disp2d>0&condition==condition_id)%>%
  ggplot(aes(x=disp2d*1000,fill=as.character(inMask_gt1)))+geom_histogram(bins = 60,aes(y=..density..),alpha=0.5,position='identity')+xlim(0,750)+
  geom_density(alpha=0.5,aes(color=as.character(inMask_gt1),fill=NULL))+scale_colour_Publication()+
  scale_fill_Publication()+theme_Publication(base_size=12)+ xlab("2D displacement (nm)")+ theme(legend.position = "none")+scale_y_continuous(expand=c(0,0))+scale_x_continuous(expand=c(0,0),limits=c(0,500))
p
ggsave(p,filename = paste0("/media/DATA/Maarten/OneDrive/Data2/MFM/R plots/",condition_id,"_displacements_inout_2d_mask_gt1.pdf"),width = 10,height = 5,units = "cm",bg = "transparent")
ggsave(p,filename = paste0("/media/DATA/Maarten/OneDrive/Data2/MFM/R plots/",condition_id,"_displacements_inout_2d_mask_gt1.png"),width = 10,height = 5,units = "cm",bg = "transparent")

p <- segs_nest %>%
  filter(disp3d>0&condition==condition_id)%>%
  ggplot(aes(x=disp3d*1000,fill=as.character(inMask)))+geom_histogram(bins = 60,aes(y=..density..),alpha=0.5,position='identity')+xlim(0,750)+
  geom_density(alpha=0.5,aes(color=as.character(inMask),fill=NULL))+scale_colour_Publication()+
  scale_fill_Publication()+theme_Publication(base_size=12)+ xlab("3D displacement (nm)")+ theme(legend.position = "none")+scale_y_continuous(expand=c(0,0))+scale_x_continuous(expand=c(0,0),limits=c(0,500))
p
ggsave(p,filename = paste0("/media/DATA/Maarten/OneDrive/Data2/MFM/R plots/",condition_id,"_displacements_inout_3d_mask.pdf"),width = 10,height = 5,units = "cm",bg = "transparent")
ggsave(p,filename = paste0("/media/DATA/Maarten/OneDrive/Data2/MFM/R plots/",condition_id,"_displacements_inout_3d_mask.png"),width = 10,height = 5,units = "cm",bg = "transparent")

p <- segs_nest %>%
  filter(disp3d>0&condition==condition_id)%>%
  ggplot(aes(x=disp3d*1000,fill=as.character(inMask_gt1)))+geom_histogram(bins = 60,aes(y=..density..),alpha=0.5,position='identity')+xlim(0,750)+
  geom_density(alpha=0.5,aes(color=as.character(inMask_gt1),fill=NULL))+scale_colour_Publication()+
  scale_fill_Publication()+theme_Publication(base_size=12)+ xlab("3D displacement (nm)")+ theme(legend.position = "none")+scale_y_continuous(expand=c(0,0))+scale_x_continuous(expand=c(0,0),limits=c(0,500))
p
ggsave(p,filename = paste0("/media/DATA/Maarten/OneDrive/Data2/MFM/R plots/",condition_id,"_displacements_inout_3d_mask_gt1.pdf"),width = 10,height = 5,units = "cm",bg = "transparent")
ggsave(p,filename = paste0("/media/DATA/Maarten/OneDrive/Data2/MFM/R plots/",condition_id,"_displacements_inout_3d_mask_gt1.png"),width = 10,height = 5,units = "cm",bg = "transparent")

#plotting of displacements in out gt1
p <- segs_nest %>%
  filter(displacement1>0&condition==condition_id)%>%
  ggplot(aes(x=displacement1*1000,fill=as.character(inMask_gt1)))+geom_histogram(bins = 60,aes(y=..density..),alpha=0.5,position='identity')+
  geom_density(aes(alpha=0.5,color=as.character(inMask_gt1),fill=NULL))+scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=12)+
  xlab("2D displacement (nm)")+ theme(legend.position = "none")+scale_y_continuous(expand=c(0,0))+scale_x_continuous(expand=c(0,0),limits=c(0,500))
p
ggsave(p,filename = paste0( "/media/DATA/Maarten/OneDrive/Data2/MFM/R plots/",condition_id,"_displacements_inout_mask_gt1.pdf"),width = 10,height = 5,units = "cm",bg = "transparent")
ggsave(p,filename =  paste0("/media/DATA/Maarten/OneDrive/Data2/MFM/R plots/",condition_id,"_displacements_inout_mask_gt1.png"),width = 10,height = 5,units = "cm",bg = "transparent")


plot_data <-   dplyr::filter(segs_nest,condition==condition_id)
  
p <-  rbind(plot_data %>%
             group_by(.id,inMask,state)%>%
             dplyr::summarise(n=n()) %>%
             group_by(.id,inMask) %>%
             dplyr:: mutate(fraction=n/sum(n))%>%
             dplyr::mutate(labels=paste0(c("outside","inside")[inMask+1],"_",c("fast","slow",'immobile')[state+1]))%>%
             dplyr::mutate(labels=factor(labels,levels=c("outside_immobile","outside_slow","outside_fast","inside_immobile","inside_slow","inside_fast")))%>%
             group_by(labels),
           plot_data%>%
             group_by(.id,inMask_gt1,state)%>%
             dplyr::summarise(n=n()) %>%
             group_by(.id,inMask_gt1) %>%
             dplyr::mutate(fraction=n/sum(n))%>%
             dplyr::mutate(labels=paste0(c("outside_GT","inside_GT")[inMask_gt1+1],"_",c("fast","slow",'immobile')[state+1]))%>%
             dplyr::mutate(labels=factor(labels,levels=c("outside_GT_immobile","outside_GT_slow","outside_GT_fast","inside_GT_immobile","inside_GT_slow","inside_GT_fast")))%>%
             group_by(labels)) %>%
  dplyr::filter(labels!="outside_GT_immobile"&labels!="outside_GT_slow"&labels!="outside_GT_fast")%>%
  ggplot(aes(x=labels,y=fraction, fill=labels))+geom_boxplot()+geom_quasirandom(size=0.1)+xlab("")+
  theme_Publication(base_size=30)+scale_colour_Publication()+scale_fill_Publication()+
  theme(legend.position = "none",text = element_text(size=20),axis.text.x = element_text(angle = 90))+
  scale_y_continuous(expand=c(0,0), limits = c(0,1.3))+
    geom_signif(comparisons = list(c("outside_immobile", "inside_immobile"),c("inside_immobile", "inside_GT_immobile")), 
              map_signif_level=F,     y_position = c(1,1.1))
p

ggsave(p,filename = paste0("/media/DATA/Maarten/OneDrive/Data2/MFM/R plots/",condition_id,"_WT_MMC_mean_fraction_localizations_beeswarm.pdf"),width = 20,height = 20,units = "cm")
ggsave(p,filename = paste0("/media/DATA/Maarten/OneDrive/Data2/MFM/R plots/",condition_id,"_WT_MMC_mean_fraction_localizations_beeswarm.png"),width = 20,height = 20,units = "cm")
plot_data <-   dplyr::filter(segs_nest,condition==condition_id)

#plot fraction of tracklets
plot_data <-   dplyr::filter(segs_nest,condition==condition_id,D_ML>0)

p <- rbind(plot_data %>%
        group_by(.id,trackletInMask,state)%>%
        dplyr::summarise(n=n_distinct(tracklet)) %>%
        group_by(.id,trackletInMask) %>%
        dplyr:: mutate(fraction=n/sum(n))%>%
        dplyr::mutate(labels=paste0(c("outside","inside")[trackletInMask+1],"_",c("fast","slow",'immobile')[state+1]))%>%
        dplyr::mutate(labels=factor(labels,levels=c("outside_immobile","outside_slow","outside_fast","inside_immobile","inside_slow","inside_fast")))%>%
        group_by(labels),
      plot_data%>%
         group_by(.id,trackletInMask_gt1,state)%>%
         dplyr::summarise(n=n_distinct(tracklet)) %>%
         group_by(.id,trackletInMask_gt1) %>%
         dplyr::mutate(fraction=n/sum(n))%>%
         dplyr::mutate(labels=paste0(c("outside_GT","inside_GT")[trackletInMask_gt1+1],"_",c("fast","slow",'immobile')[state+1]))%>%
         dplyr::mutate(labels=factor(labels,levels=c("outside_GT_immobile","outside_GT_slow","outside_GT_fast","inside_GT_immobile","inside_GT_slow","inside_GT_fast")))%>%
         group_by(labels)) %>%
        dplyr::filter(labels!="outside_GT_immobile"&labels!="outside_GT_slow"&labels!="outside_GT_fast")%>%
  ggplot(aes(x=labels,y=fraction, fill=labels))+geom_boxplot()+geom_quasirandom(size=0.1)+xlab("")+
  theme_Publication(base_size=30)+scale_colour_Publication()+scale_fill_Publication()+
  theme(legend.position = "none",text = element_text(size=20),axis.text.x = element_text(angle = 90))+
  scale_y_continuous(expand=c(0,0), limits = c(0,1.3))+
    geom_signif(comparisons = list(c("outside_immobile", "inside_immobile"),c("inside_immobile", "inside_GT_immobile")), 
              map_signif_level=F,     y_position = c(1,1.1))
p

ggsave(p,filename = paste0("/media/DATA/Maarten/OneDrive/Data2/MFM/R plots/",condition_id,"_WT_MMC_mean_fraction_tracklets_beeswarm.pdf"),width = 20,height = 20,units = "cm")
ggsave(p,filename = paste0("/media/DATA/Maarten/OneDrive/Data2/MFM/R plots/",condition_id,"_WT_MMC_mean_fraction_tracklets_beeswarm.png"),width = 20,height = 20,units = "cm")

rbind(plot_data %>%
        group_by(.id,trackletInMask,state)%>%
        dplyr::summarise(n=n_distinct(tracklet)) %>%
        group_by(.id,trackletInMask) %>%
        dplyr:: mutate(fraction=n/sum(n))%>%
        dplyr::mutate(labels=paste0(c("outside","inside")[trackletInMask+1],"_",c("fast","slow",'immobile')[state+1]))%>%
        dplyr::mutate(labels=factor(labels,levels=c("outside_immobile","outside_slow","outside_fast","inside_immobile","inside_slow","inside_fast")))%>%
        group_by(labels),
      plot_data%>%
        group_by(.id,trackletInMask_gt1,state)%>%
        dplyr::summarise(n=n_distinct(tracklet)) %>%
        group_by(.id,trackletInMask_gt1) %>%
        dplyr::mutate(fraction=n/sum(n))%>%
        dplyr::mutate(labels=paste0(c("outside_GT","inside_GT")[trackletInMask_gt1+1],"_",c("fast","slow",'immobile')[state+1]))%>%
        dplyr::mutate(labels=factor(labels,levels=c("outside_GT_immobile","outside_GT_slow","outside_GT_fast","inside_GT_immobile","inside_GT_slow","inside_GT_fast")))%>%
        group_by(labels)) %>%
  dplyr::filter(labels!="outside_GT_immobile"&labels!="outside_GT_slow"&labels!="outside_GT_fast") %>%
  group_by(labels)%>%
  dplyr::summarize(mean=mean(fraction),median=median(fraction))


# plot fraction of localizations in the mask ------------------------------


  
p <- segs_nest %>%
  filter(condition==condition_id&D_ML>0)%>%
  group_by(cellID)%>%
  dplyr::summarise(n=n(),data=sum(inMask)/n(),random=sum(inMask_gt1)/n())%>%
  melt(id.vars="cellID",measure.vars=c("data","random"))%>%
  ggplot(aes(y=value,x=variable,fill=variable))+ geom_boxplot(outlier.shape = NA,notch = T)+geom_quasirandom()+
  scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=8)+ theme(legend.position = "none")+
  xlab("")+ylab("fraction localizations inside mask")+  
  scale_y_continuous(expand=c(0,0),limits = c(0,0.3)) +
geom_signif(comparisons = list(c("data", "random")), 
            map_signif_level=F,     y_position = c(0.25))
p

segs_nest %>%
filter(condition==condition_id&D_ML>0)%>%
  group_by(cellID)%>%
  dplyr::summarise(n=n(),data=sum(inMask)/n(),random=sum(inMask_gt1)/n())%>%
  melt(id.vars="cellID",measure.vars=c("data","random")) %>% group_by(variable) %>% dplyr::summarize(mean(value),median(value))

ggsave(p,filename = paste0("/media/DATA/Maarten/OneDrive/Data2/MFM/R plots/",condition_id,"_fraction_localizations_inside.pdf"),width = 10,height = 10,units = "cm")
ggsave(p,filename = paste0("/media/DATA/Maarten/OneDrive/Data2/MFM/R plots/",condition_id,"_fraction_localizations_inside.png"),width = 10,height = 10,units = "cm")


stats <- segs_nest %>%
  filter(condition==condition_id)%>%
  group_by(cellID)%>%
  dplyr::summarise(n=n(),data=sum(inMask)/n(),random=sum(inMask_gt1)/n())%>%
  melt(id.vars="cellID",measure.vars=c("data","random"))%>%
  group_by(variable,cellID) %>%
  dplyr::summarize(x=mean(value)) %>%
  group_by(variable)%>%
  dplyr::summarize(out=mean(x))

stats$out[1]/stats$out[2]




# angle plots -------------------------------------------------------------
#angle histogram in out

inputdata <- segs_nest %>%
  dplyr::filter(displacement1>0.1,angle1>0,state!=2,condition==condition_id)
p <-inputdata %>%
  dplyr::mutate(angle1=360-angle1)%>%
  rbind(inputdata)%>%
 ggplot(aes(x=angle1,fill=as.character(inMask)))+geom_histogram(breaks=seq(0, 360, 10),aes(y=..density..),alpha=0.4,position='identity')+
  scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=16)+theme(axis.line = element_line(size = 0.5),axis.ticks = element_line(size = 0.5)) +
  scale_y_continuous(expand=c(0,0),limits = c(0,0.010)) + scale_x_continuous(expand=c(0,0),limits = c(0,360),breaks=c(0,90,180,270,360))+xlab("")+ theme(legend.position = "none")

p 
ggsave(p,filename = paste0("/media/DATA/Maarten/OneDrive/Data2/MFM/R plots/",condition_id,"_angle_histogram_in_out.pdf"),width = 10,height = 10,units = "cm")
ggsave(p,filename = paste0("/media/DATA/Maarten/OneDrive/Data2/MFM/R plots/",condition_id,"_angle_histogram_in_out.png"),width = 10,height = 10,units = "cm")

#angle histogram gt_1
inputdata <- segs_nest %>%
  dplyr::filter(displacement1>0.1,angle1>0,state!=2,condition==condition_id)
p <-inputdata %>%
  dplyr::mutate(angle1=360-angle1)%>%
  rbind(inputdata)%>%
  ggplot(aes(x=angle1,fill=as.character(inMask_gt1)))+geom_histogram(breaks=seq(0, 360, 10),aes(y=..density..),alpha=0.4,position='identity')+
  scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=16)+ 
  scale_y_continuous(expand=c(0,0),limits = c(0,0.010)) + scale_x_continuous(expand=c(0,0),limits = c(0,360),breaks=c(0,90,180,270,360))+xlab("")+
  theme(legend.position = "none")
p
ggsave(p,filename = paste0("/media/DATA/Maarten/OneDrive/Data2/MFM/R plots/",condition_id,"_angle_histogram_in_out_gt1.pdf"),width = 10,height = 10,units = "cm")
ggsave(p,filename = paste0("/media/DATA/Maarten/OneDrive/Data2/MFM/R plots/",condition_id,"_angle_histogram_in_out_gt1.png"),width = 10,height = 10,units = "cm")


#angle histogram stat1 slow
inputdata <- segs_nest %>%
  filter(displacement1>0.1,angle1>0,state==1,condition==condition_id)
p <- inputdata %>%
  dplyr::mutate(angle1=360-angle1)%>%
  rbind(inputdata)%>%
  ggplot(aes(x=angle1,fill=as.character(inMask)))+geom_histogram(breaks=seq(0, 360, 10),aes(y=..density..),
                                                                alpha=0.4,position='identity')+scale_colour_Publication()+
  scale_fill_Publication()+theme_Publication(base_size=16)+ 
  scale_y_continuous(expand=c(0,0),limits = c(0,0.010)) + scale_x_continuous(expand=c(0,0),limits = c(0,360),breaks=c(0,90,180,270,360))+xlab("")+
  theme(legend.position = "none")

p
ggsave(p,filename = paste0("/media/DATA/Maarten/OneDrive/Data2/MFM/R plots/",condition_id,"_angle_histogram_in_out_state1.pdf"),width = 10,height = 10,units = "cm")
ggsave(p,filename = paste0("/media/DATA/Maarten/OneDrive/Data2/MFM/R plots/",condition_id,"_angle_histogram_in_out_state1.png"),width = 10,height = 10,units = "cm")

#angle histogram stat0 fast
inputdata <- segs_nest %>%
  filter(displacement1>0.1,angle1>0,state!=2,condition==condition_id,inMask==T)
p <- inputdata %>%
  dplyr::mutate(angle1=360-angle1)%>%
  rbind(inputdata)%>%
  ggplot(aes(x=angle1,fill=as.character(state)))+geom_histogram(breaks=seq(0, 360, 10),aes(y=..density..),
                                                                alpha=0.4,position='identity')+
  scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=16)+ 
  scale_y_continuous(expand=c(0,0),limits = c(0,0.010)) + scale_x_continuous(expand=c(0,0),limits = c(0,360),breaks=c(0,90,180,270,360))+xlab("")+ theme(legend.position = "none")
p
ggsave(p,filename = paste0("/media/DATA/Maarten/OneDrive/Data2/MFM/R plots/",condition_id,"_angle_histogram_in_fast_slow.pdf"),width = 10,height = 10,units = "cm")
ggsave(p,filename = paste0("/media/DATA/Maarten/OneDrive/Data2/MFM/R plots/",condition_id,"_angle_histogram_in_fast_slow.png"),width = 10,height = 10,units = "cm")

#angle histogram stat0 fast
inputdata <- segs_nest %>%
  filter(displacement1>0.1,angle1>0,state==0,condition==condition_id)
p <- inputdata %>%
  dplyr::mutate(angle1=360-angle1)%>%
  rbind(inputdata)%>%
  ggplot(aes(x=angle1,fill=as.character(inMask)))+geom_histogram(breaks=seq(0, 360, 10),aes(y=..density..),
                                                                 alpha=0.4,position='identity')+
  scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=16)+ 
  scale_y_continuous(expand=c(0,0),limits = c(0,0.010)) + scale_x_continuous(expand=c(0,0),limits = c(0,360),breaks=c(0,90,180,270,360))+xlab("")+ theme(legend.position = "none")
p
ggsave(p,filename = paste0("/media/DATA/Maarten/OneDrive/Data2/MFM/R plots/",condition_id,"_angle_histogram_in_out_state0_fast.pdf"),width = 10,height = 10,units = "cm")
ggsave(p,filename = paste0("/media/DATA/Maarten/OneDrive/Data2/MFM/R plots/",condition_id,"_angle_histogram_in_out_state0_fast.png"),width = 10,height = 10,units = "cm")



# inputdata <- segs_nest %>%
#   filter(displacement1>0.1,angle1>0,inMask==T,state<2,condition==condition)
# p <- inputdata %>%
#   dplyr::mutate(angle1=360-angle1)%>%
#   rbind(inputdata)%>%
#   ggplot(aes(x=angle1,fill=as.character(state)))+geom_histogram(breaks=seq(0, 350, 10),aes(y=..density..),
#                                                                 alpha=0.5,position='identity')+ylim(0,0.020)+
#   scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=8)
# 
# p

#making linear line plot of fold anisotropy vs mean displacement

break_size <-0.025
data <- segs_nest %>%
  dplyr::filter(angle1>0,state<2,condition==condition_id,inMask==T,displacement1>0,displacement2>0)
data$mean_disp <- (data$displacement1+data$displacement2)/2
bins <- seq(0,0.5,break_size)

results_inside <- data.frame("mean_disp"=bins,"fold"=-1)
results_inside <- ddply(results_inside,.variables = c("mean_disp"),function(x){
  
  datasub <- subset(data,mean_disp>=x$mean_disp&mean_disp<(x$mean_disp+0.01))
  x$fold <- length(datasub$angle1[datasub$angle1>165])/length(datasub$angle1[datasub$angle1<15])
  return(x)
})
results_inside$location <- "inside"

#outside
data <- segs_nest %>%
  filter(angle1>0,state<2,condition==condition_id,inMask==F,displacement1>0,displacement2>0)
data$mean_disp <- (data$displacement1+data$displacement2)/2
bins <- seq(0,0.5,break_size)

results_outside <- data.frame("mean_disp"=bins,"fold"=-1)
results_outside <- ddply(results_outside,.variables = c("mean_disp"),function(x){
  
  datasub <- subset(data,mean_disp>=x$mean_disp&mean_disp<(x$mean_disp+break_size))
  x$fold <- length(datasub$angle1[datasub$angle1>165])/length(datasub$angle1[datasub$angle1<15])
  return(x)
})

results_outside$location <- "outside"
#gt1
data <- segs_nest %>%
  filter(angle>0,state<2,condition==condition_id,inMask_gt1==T,displacement1>0,displacement2>0)
data$mean_disp <- (data$displacement1+data$displacement2)/2
bins <- seq(0,0.5,break_size)

results_gt1 <- data.frame("mean_disp"=bins,"fold"=-1)
results_gt1 <- ddply(results_outside,.variables = c("mean_disp"),function(x){
  
  datasub <- subset(data,mean_disp>x$mean_disp&mean_disp<(x$mean_disp+break_size))
  x$fold <- length(datasub$angle1[datasub$angle1>165])/length(datasub$angle1[datasub$angle1<15])
  return(x)
})

#gt2
data <- segs_nest %>%
  filter(angle>0,state<2,condition==condition_id,inMask_gt2==T,displacement1>0,displacement2>0)
data$mean_disp <- (data$displacement1+data$displacement2)/2
bins <- seq(0,0.5,break_size)

results_gt2 <- data.frame("mean_disp"=bins,"fold"=-1)
results_gt2 <- ddply(results_outside,.variables = c("mean_disp"),function(x){
  
  datasub <- subset(data,mean_disp>x$mean_disp&mean_disp<(x$mean_disp+break_size))
  x$fold <- length(datasub$angle1[datasub$angle1>165])/length(datasub$angle1[datasub$angle1<15])
  return(x)
})

#gt3
data <- segs_nest %>%
  filter(angle1>0,state<2,condition==condition_id,inMask_gt3==T,displacement1>0,displacement2>0)
data$mean_disp <- (data$displacement1+data$displacement2)/2
bins <- seq(0,0.5,break_size)

results_gt3 <- data.frame("mean_disp"=bins,"fold"=-1)
results_gt3 <- ddply(results_outside,.variables = c("mean_disp"),function(x){
  
  datasub <- subset(data,mean_disp>x$mean_disp&mean_disp<(x$mean_disp+break_size))
  x$fold <- length(datasub$angle1[datasub$angle1>165])/length(datasub$angle1[datasub$angle1<15])
  return(x)
})


results_gt1$fold <- (results_gt1$fold+results_gt2$fold+results_gt3$fold)/3
results_gt1$location <- "mean_gt"
results <- rbind(results_inside,results_outside,results_gt1)
p <- ggplot(results,aes(x=mean_disp,y=fold,color=location))+geom_line()+geom_point()+xlab("Mean displacement (um)")+ylab("Fold anisotropy")+
  scale_colour_Publication()+theme_Publication(base_size=16)+
  scale_y_continuous(expand=c(0,0),limits = c(0,8)) + scale_x_continuous(expand=c(0,0),limits = c(0,0.5))
p


ggsave(p,filename = paste0("/media/DATA/Maarten/OneDrive/Data2/MFM/R plots/",condition_id,"_angle_lineplot.pdf"),width = 10,height = 7,units = "cm")
ggsave(p,filename = paste0("/media/DATA/Maarten/OneDrive/Data2/MFM/R plots/",condition_id,"_angle_lineplot.png"),width = 10,height = 7,units = "cm")

