#this script skips the import of data and can be used for subsequent steps of analysis

#required packages
library(tidyverse)
library(plyr)
library(lattice)
library(stats)
library(MSDtracking)
library(ggplot2)#required packages
library(ggpol)
library(doParallel)
library(reticulate)
library(ggbeeswarm)
source('R/MSD.R')
source('R/MSD_fit.R')
source('R/analysis functions.R')
source('python/ML_py.R')

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
  discrete_scale("fill","Publication",manual_pal(values = c("#1f497d","#6599d9","#c00000","#1f497d","#6599d9","#c00000","#1f497d","#6599d9","#c00000","#1f497d","#6599d9","#c00000")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#1f497d","#6599d9","#c00000","#1f497d","#6599d9","#c00000","#1f497d","#6599d9","#c00000","#1f497d","#6599d9","#c00000")), ...)
  
}

#input variables
directory <- "/media/DATA/Maarten/MFM/data_gtv2/"

condition_list <- list.dirs(directory,full.names = F,recursive = F)
#condition_list <- condition_list[c(2,3,5)]

load(file=file.path(directory,"segs_nest.Rdata"))

segs_nest <- dplyr::filter(segs_nest,condition=="WT MMC")
segs_nest <- segs_nest %>%
  group_by(.id,tracklet,condition)%>%
  dplyr::mutate(trackletInMask=any(inMask==T))

segs_nest <- segs_nest %>%
  group_by(.id,tracklet,condition)%>%
  dplyr::mutate(trackletInMask_gt1=any(inMask_gt1==T))


# Diffusion histograms in-out mask ----------------------------------------
# k <- segs_nest %>%
#   filter(D_ML_focus>0)%>%
#   dplyr::distinct(condition,cellID,tracklet,.keep_all=T)%>%
#   mutate(state_str=as.character(inMask))
# 
# x <- group_by(k,condition) %>%
#   group_by(condition,inMask)%>%
#   dplyr::summarise(number=n())%>%
#   group_by(condition) %>%
#   dplyr::mutate(fraction=number/sum(number))
# 
# y <- group_by(k,condition) %>%
#   group_by(condition,cellID,inMask)%>%
#   dplyr::summarise(number=n())%>%
#   group_by(condition,cellID) %>%
#   dplyr::mutate(fraction=number/sum(number)) %>%
#   group_by(condition,inMask) %>%
#   dplyr::summarise(mean_fraction=round(mean(fraction),digits = 2),sd_fraction=round(sd(fraction),digits = 2))

#histogram in-out foci
p <- segs_nest %>%
  filter(D_ML_focus>0,condition=="WT MMC")%>%
  dplyr::distinct(condition,cellID,focus_tracklet,.keep_all=T)%>%
  mutate(state_str=as.character(inMask))%>%
  ggplot(aes(x=D_ML_focus*100,fill=state_str,y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]))+geom_histogram(position="identity",alpha=0.5)+scale_x_log10(limits=c(0.0001,10))+facet_wrap(.~state_str,ncol=2)+
  scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=8)+ theme(legend.position = "none")+ylab("relative frequency")+ xlab(expression(D[app]~mu~m^{2}/s)) 
p
ggsave(p,filename = "/media/DATA/Maarten/OneDrive/Documents/Manuscripts/in preparation - MFM BRCA2 tracking/Figure 2/msd_histogram_inout_mask.pdf",width = 10,height = 5,units = "cm")
#histogram in-out foci gt1
p <- segs_nest %>%
  filter(D_ML_focus_gt1>0,condition=="WT MMC")%>%
  dplyr::distinct(condition,cellID,focus_tracklet_gt1,.keep_all=T)%>%
  mutate(state_str=as.character(inMask_gt1))%>%
  ggplot(aes(x=D_ML_focus_gt1*100,fill=state_str,y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]))+geom_histogram(position="identity",alpha=0.5)+scale_x_log10(limits=c(0.0001,10))+facet_wrap(.~state_str,ncol=2)+
  scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=8)+ theme(legend.position = "none")+ylab("relative frequency")+ xlab(expression(D[app]~mu~m^{2}/s)) 
p
ggsave(p,filename = "/media/DATA/Maarten/OneDrive/Documents/Manuscripts/in preparation - MFM BRCA2 tracking/Figure 2/msd_histogram_inout_mask_gt1.pdf",width = 10,height = 5,units = "cm")

#plotting of displacements in out
p <- segs_nest %>%
  filter(displacement1>0)%>%
  ggplot(aes(x=displacement1*1000,fill=as.character(inMask)))+geom_histogram(bins = 60,aes(y=..density..),alpha=0.5,position='identity')+
  geom_density(alpha=0.5,aes(color=as.character(inMask),fill=NULL))+scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=12)+ xlab("x-y displacement (nm)")+ theme(legend.position = "none")
p
ggsave(p,filename = "/media/DATA/Maarten/OneDrive/Documents/Manuscripts/in preparation - MFM BRCA2 tracking/Figure 2/displacements_inout_mask.pdf",width = 10,height = 5,units = "cm")

#plotting of displacements in out gt1

p <- segs_nest %>%
  filter(displacement1>0)%>%
  ggplot(aes(x=displacement1*1000,fill=as.character(inMask_gt1)))+geom_histogram(bins = 60,aes(y=..density..),alpha=0.5,position='identity')+
  geom_density(aes(alpha=0.5,color=as.character(inMask_gt1),fill=NULL))+scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=12)+ xlab("x-y displacement (nm)")+ theme(legend.position = "none")
p
ggsave(p,filename = "/media/DATA/Maarten/OneDrive/Documents/Manuscripts/in preparation - MFM BRCA2 tracking/Figure 2/displacements_inout_mask_gt1.pdf",width = 10,height = 5,units = "cm")


###make D histograms fractions
#MSD
msd <- segs_nest %>%
  dplyr::filter(D_ML>0,condition=="WT MMC")%>%
  dplyr::distinct(condition,cellID,tracklet,.keep_all=T)%>%
  dplyr::mutate(state_str=as.character(state))%>%
  ggplot(aes(x=D_ML*100,y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..],fill=state_str))+geom_histogram(position="identity",alpha=0.5)+scale_x_log10(limits=c(0.0001,10))+
  scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=12)+ theme(legend.position = "none")+ylab("relative frequency")+ xlab(expression(D[app]~mu~m^{2}/s))
msd
ggsave(msd,filename = "/media/DATA/Maarten/OneDrive/Documents/Manuscripts/in preparation - MFM BRCA2 tracking/Figure 2/msd_plot.pdf",width = 10,height = 5,units = "cm")

#MSS
mss <- segs_nest %>%
  dplyr::filter(D_ML>0,condition=="WT MMC")%>%
  dplyr::distinct(condition,cellID,tracklet,.keep_all=T)%>%
  dplyr::mutate(state_str=as.character(state))%>%
  ggplot(aes(x=D_Smmss,y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..],fill=state_str))+geom_histogram(position="identity",alpha=0.5)+xlim(-0.5,1)+
  scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=12)+ theme(legend.position = "none")+ylab("relative frequency")+ xlab(expression(S[MSS]))
mss
ggsave(mss,filename = "/media/DATA/Maarten/OneDrive/Documents/Manuscripts/in preparation - MFM BRCA2 tracking/Figure 2/mss_plot.pdf",width = 10,height = 5,units = "cm")


#scatter plot
scatter <- segs_nest %>%
  dplyr::filter(D_ML>0,condition=="WT MMC")%>%
  dplyr::distinct(condition,cellID,tracklet,.keep_all=T)%>%
  dplyr::mutate(state_str=as.character(state))%>%
  ggplot(aes(x=D_ML*100,y=D_Smmss,color=state_str))+geom_point(alpha=0.05,size=0.5)+scale_x_log10(limits=c(0.0001,10))+ylim(-0.5,1)+
  scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=12)+ theme(legend.position = "none")+ylab(expression(S[MSS]))+ xlab(expression(D[app]~mu~m^{2}/s))
scatter
ggsave(scatter,filename = "/media/DATA/Maarten/OneDrive/Documents/Manuscripts/in preparation - MFM BRCA2 tracking/Figure 2/scatter_plot.pdf",width = 10,height = 10,units = "cm")
ggsave(scatter,filename = "/media/DATA/Maarten/OneDrive/Documents/Manuscripts/in preparation - MFM BRCA2 tracking/Figure 2/scatter_plot.png",width = 10,height = 10,units = "cm")


#inside MSD
msd <- segs_nest %>%
  dplyr::filter(D_ML>0,condition=="WT MMC",trackletInMask==T)%>%
  dplyr::distinct(condition,cellID,tracklet,.keep_all=T)%>%
  dplyr::mutate(state_str=as.character(state))%>%
  ggplot(aes(x=D_ML*100,y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..],fill=state_str))+geom_histogram(position="identity",alpha=0.5)+scale_x_log10(limits=c(0.0001,10))+
  scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=12)+ theme(legend.position = "none")+ylab("relative frequency")+ xlab(expression(D[app]~mu~m^{2}/s))+
  scale_y_continuous(expand=c(0,0))

msd
ggsave(msd,filename = "/media/DATA/Maarten/OneDrive/Documents/Manuscripts/in preparation - MFM BRCA2 tracking/Figure 2/msd_plot_inside.pdf",width = 10,height = 5,units = "cm")

#outside MSD
msd <- segs_nest %>%
  dplyr::filter(D_ML>0,condition=="WT MMC",trackletInMask==F)%>%
  dplyr::distinct(condition,cellID,tracklet,.keep_all=T)%>%
  dplyr::mutate(state_str=as.character(state))%>%
  ggplot(aes(x=D_ML*100,y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..],fill=state_str))+geom_histogram(position="identity",alpha=0.5)+scale_x_log10(limits=c(0.0001,10))+
  scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=12)+ theme(legend.position = "none")+ylab("relative frequency")+ xlab(expression(D[app]~mu~m^{2}/s))+  scale_y_continuous(expand=c(0,0))

msd
ggsave(msd,filename = "/media/DATA/Maarten/OneDrive/Documents/Manuscripts/in preparation - MFM BRCA2 tracking/Figure 2/msd_plot_outside.pdf",width = 10,height = 5,units = "cm")

#random MSD
msd <- segs_nest %>%
  dplyr::filter(D_ML>0,condition=="WT MMC",trackletInMask_gt1==T)%>%
  dplyr::distinct(condition,cellID,tracklet,.keep_all=T)%>%
  dplyr::mutate(state_str=as.character(state))%>%
  ggplot(aes(x=D_ML*100,y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..],fill=state_str))+geom_histogram(position="identity",alpha=0.5)+scale_x_log10(limits=c(0.0001,10))+
  scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=12)+ theme(legend.position = "none")+ylab("relative frequency")+ xlab(expression(D[app]~mu~m^{2}/s))+  scale_y_continuous(expand=c(0,0))

msd
ggsave(msd,filename = "/media/DATA/Maarten/OneDrive/Documents/Manuscripts/in preparation - MFM BRCA2 tracking/Figure 2/msd_plot_random.pdf",width = 10,height = 5,units = "cm")


# plot fraction of segments per cell --------------------------------------
segments <- segs_nest %>%
  filter(condition=="WT MMC"&D_ML>0)
p <- rbind(segments %>%
             group_by(.id,inMask,state)%>%
             dplyr::summarise(n=n()) %>%
             group_by(.id,inMask) %>%
             dplyr::mutate(fraction=n/sum(n))%>%
             dplyr::mutate(labels=paste0(c("outside","inside")[inMask+1],"_",c("fast","slow",'immobile')[state+1]))%>%
             mutate(labels=factor(labels,levels=c("outside_immobile","outside_slow","outside_fast","inside_immobile","inside_slow","inside_fast")))%>%
             group_by(labels,inMask)%>%
             dplyr::summarise(var=sd(fraction,na.rm = T)/sqrt(n()),fraction=mean(fraction)),
           segments %>%
             group_by(.id,inMask_gt1,state)%>%
             dplyr::summarise(n=n()) %>%
             group_by(.id,inMask_gt1) %>%
             dplyr:: mutate(fraction=n/sum(n))%>%
             dplyr::mutate(labels=paste0(c("outside_gt1","inside_gt1")[inMask_gt1+1],"_",c("fast","slow",'immobile')[state+1]))%>%
             dplyr::mutate(labels=factor(labels,levels=c("outside_gt1_immobile","outside_gt1_slow","outside_gt1_fast","inside_gt1_immobile","inside_gt1_slow","inside_gt1_fast")))%>%
             group_by(labels,inMask_gt1)%>%
             dplyr::summarise(var=sd(fraction,na.rm = T)/sqrt(n()),fraction=mean(fraction))%>%
             filter(inMask_gt1==1))%>%
  ggplot(aes(y=fraction, x=labels, fill=labels))+geom_bar(stat="identity")+ylim(0,1)+geom_errorbar(stat="identity",
                                                                                                   aes(ymax=fraction+var,ymin=fraction-var))+theme_Publication(base_size=25)+xlab("")+
  scale_colour_Publication()+scale_fill_Publication()+
  theme(legend.position = "none",text = element_text(size=15),axis.text.x = element_text(angle = 90,size=6,hjust=1,vjust=0.5))

p


ggsave(p,filename = file.path("/media/DATA/Maarten/OneDrive/Documents/Manuscripts/in preparation - MFM BRCA2 tracking/Figure 2/","WT_MMC_mean_fraction_segments.pdf"))
ggsave(p,filename = file.path("/media/DATA/Maarten/OneDrive/Documents/Manuscripts/in preparation - MFM BRCA2 tracking/Figure 2/","WT_MM_mean_fraction_segments.png"))



data <- segs_nest %>%
  filter(condition=="WT MMC"&D_ML>0) 

p <- rbind(data %>%
        group_by(.id,trackletInMask,state)%>%
        dplyr::summarise(n=n_distinct(tracklet)) %>%
        group_by(.id,trackletInMask) %>%
        dplyr:: mutate(fraction=n/sum(n))%>%
        dplyr::mutate(labels=paste0(c("outside","inside")[trackletInMask+1],"_",c("fast","slow",'immobile')[state+1]))%>%
        dplyr::mutate(labels=factor(labels,levels=c("outside_immobile","outside_slow","outside_fast","inside_immobile","inside_slow","inside_fast")))%>%
        group_by(labels)%>%
        dplyr::summarise(var=sd(fraction,na.rm = T)/sqrt(n()),fraction=mean(fraction),state=state[1],trackInMask=trackletInMask[1]),
      (data%>%
         group_by(.id,trackletInMask_gt1,state)%>%
         dplyr::summarise(n=n_distinct(tracklet)) %>%
         group_by(.id,trackletInMask_gt1) %>%
         dplyr::mutate(fraction=n/sum(n))%>%
         dplyr::mutate(labels=paste0(c("outside","inside_GT")[trackletInMask_gt1+1],"_",c("fast","slow",'immobile')[state+1]))%>%
         dplyr::mutate(labels=factor(labels,levels=c("outside_immobile","outside_slow","outside_fast","inside_GT_immobile","inside_GT_slow","inside_GT_fast")))%>%
         group_by(labels)%>%
         dplyr::summarise(var=sd(fraction,na.rm = T)/sqrt(n()),fraction=mean(fraction),state=state[1],trackInMask=trackletInMask_gt1[1]))[4:6,] ) %>%
  ggplot(aes(x=labels, fill=labels))+geom_bar(stat="identity",aes(y=fraction))+xlab("")+
  ylim(0,1)+geom_errorbar(stat="identity",aes(ymax=fraction+var,ymin=fraction-var)) +theme_Publication(base_size=25)+scale_colour_Publication()+scale_fill_Publication()+theme(legend.position = "none",text = element_text(size=15),axis.text.x = element_text(angle = 90))+
  scale_y_continuous(expand=c(0,0))


p

trackletInMask

ggsave(p,filename = file.path("/media/DATA/Maarten/OneDrive/Documents/Manuscripts/in preparation - MFM BRCA2 tracking/Figure 2/","WT_MMC_mean_fraction_tracklets.pdf"))
ggsave(p,filename = file.path("/media/DATA/Maarten/OneDrive/Documents/Manuscripts/in preparation - MFM BRCA2 tracking/Figure 2/","WT_MM_mean_fraction_tracklets.png"))

p <- rbind(data %>%
        group_by(.id,trackInMask,state)%>%
        dplyr::summarise(n=n_distinct(tracklet)) %>%
        group_by(.id,trackInMask) %>%
        dplyr:: mutate(fraction=n/sum(n))%>%
        dplyr::mutate(labels=paste0(c("outside","inside")[trackInMask+1],"_",c("fast","slow",'immobile')[state+1]))%>%
        dplyr::mutate(labels=factor(labels,levels=c("outside_immobile","outside_slow","outside_fast","inside_immobile","inside_slow","inside_fast")))%>%
        group_by(labels),
      data%>%
         group_by(.id,trackInMask_gt1,state)%>%
         dplyr::summarise(n=n_distinct(tracklet)) %>%
         group_by(.id,trackInMask_gt1) %>%
         dplyr::mutate(fraction=n/sum(n))%>%
         dplyr::mutate(labels=paste0(c("outside_GT","inside_GT")[trackInMask_gt1+1],"_",c("fast","slow",'immobile')[state+1]))%>%
         dplyr::mutate(labels=factor(labels,levels=c("outside_GT_immobile","outside_GT_slow","outside_GT_fast","inside_GT_immobile","inside_GT_slow","inside_GT_fast")))%>%
         group_by(labels)) %>%
        dplyr::filter(labels!="outside_GT_immobile"&labels!="outside_GT_slow"&labels!="outside_GT_fast")%>%
  ggplot(aes(x=labels,y=fraction, fill=labels))+geom_boxplot()+geom_quasirandom(size=0.1)+xlab("")+
  ylim(0,1) +theme_Publication(base_size=25)+scale_colour_Publication()+scale_fill_Publication()+
  theme(legend.position = "none",text = element_text(size=15),axis.text.x = element_text(angle = 90))+
  scale_y_continuous(expand=c(0,0))

  p

  ggsave(p,filename = file.path("/media/DATA/Maarten/OneDrive/Documents/Manuscripts/in preparation - MFM BRCA2 tracking/Figure 2/","WT_MMC_mean_fraction_tracklets_beeswarm.pdf"))
  ggsave(p,filename = file.path("/media/DATA/Maarten/OneDrive/Documents/Manuscripts/in preparation - MFM BRCA2 tracking/Figure 2/","WT_MM_mean_fraction_tracklets_beeswarm.png"))
# plot fraction of localizations in the mask ------------------------------

library(ggbeeswarm)

p <- segs_nest %>%
  filter(condition=="WT MMC")%>%
  group_by(cellID)%>%
  dplyr::summarise(n=n(),data=sum(inMask)/n(),random=sum(inMask_gt1)/n())%>%
  melt(id.vars="cellID",measure.vars=c("data","random"))%>%
  ggplot(aes(y=value,x=variable,fill=variable))+ geom_boxplot(outlier.shape = NA,notch = T)+geom_quasirandom()+ylim(0,0.5)+
  scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=8)+ theme(legend.position = "none")+
  xlab("")+ylab("fraction localizations inside mask")
p
ggsave(p,filename = "/media/DATA/Maarten/OneDrive/Documents/Manuscripts/in preparation - MFM BRCA2 tracking/Figure 2/fraction_localizations_inside.pdf",width = 10,height = 10,units = "cm")

stats <- segs_nest %>%
  filter(condition=="WT MMC")%>%
  group_by(cellID)%>%
  dplyr::summarise(n=n(),data=sum(inMask)/n(),random=sum(inMask_gt1)/n())%>%
  melt(id.vars="cellID",measure.vars=c("data","random"))%>%
  group_by(variable,cellID) %>%
  dplyr::summarize(x=mean(value)) %>%
  group_by(variable)%>%
  dplyr::summarize(out=mean(x))


# angle plots -------------------------------------------------------------
#angle histogram in out

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#1F497D","#c00000","#542788","#217D68","#386cb0","#fdb462","#386cb0","#fdb462","#386cb0","#fdb462","#386cb0","#fdb462")), ...)
  
}
scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#c00000","#6599D9","#542788","#217D68","#386cb0","#fdb462","#386cb0","#fdb462","#386cb0","#fdb462","#386cb0","#fdb462")), ...)
  
}
inputdata <- segs_nest %>%
  dplyr::filter(displacement1>0.1,angle>0,state<2,condition=="WT MMC")
p <-inputdata %>%
  dplyr::mutate(angle=360-angle)%>%
  rbind(inputdata)%>%
 ggplot(aes(x=angle,fill=as.character(inMask)))+geom_histogram(breaks=seq(0, 350, 10),aes(y=..density..),alpha=0.5,position='identity')+
  ylim(0,0.007)+scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=8)
p 
ggsave(p,filename = "/media/DATA/Maarten/OneDrive/Documents/Manuscripts/in preparation - MFM BRCA2 tracking/Figure 2/angle_histogram_in_out.pdf",width = 10,height = 10,units = "cm")

#angle histogram gt_1
inputdata <- segs_nest %>%
  filter(displacement1>0.1,angle>0,state<2,condition=="WT MMC")
p <- inputdata %>%
  dplyr::mutate(angle=360-angle)%>%
  rbind(inputdata)%>%
  ggplot(aes(x=angle,fill=as.character(inMask_gt1)))+geom_histogram(breaks=seq(0, 350, 10),aes(y=..density..),
                                                                alpha=0.5,position='identity')+ylim(0,0.007)+scale_colour_Publication()+scale_fill_Publication()+
  theme_Publication(base_size=8)
p
ggsave(p,filename = "/media/DATA/Maarten/OneDrive/Documents/Manuscripts/in preparation - MFM BRCA2 tracking/Figure 2/angle_histogram_in_out_gt1.pdf",width = 10,height = 10,units = "cm")

#angle histogram stat1 slow
inputdata <- segs_nest %>%
  filter(displacement1>0.1,angle>0,state==1,condition=="WT MMC")
p <- inputdata %>%
  dplyr::mutate(angle=360-angle)%>%
  rbind(inputdata)%>%
  ggplot(aes(x=angle,fill=as.character(inMask)))+geom_histogram(breaks=seq(0, 350, 10),aes(y=..density..),
                                                                alpha=0.5,position='identity')+ylim(0,0.008)+scale_colour_Publication()+
  scale_fill_Publication()+theme_Publication(base_size=8)
p
ggsave(p,filename = "/media/DATA/Maarten/OneDrive/Documents/Manuscripts/in preparation - MFM BRCA2 tracking/Figure 2/angle_histogram_in_out_state1.pdf",width = 10,height = 10,units = "cm")

#angle histogram stat0 fast

inputdata <- segs_nest %>%
  filter(displacement1>0.1,angle>0,state==0,condition=="WT MMC")
p <- inputdata %>%
  dplyr::mutate(angle=360-angle)%>%
  rbind(inputdata)%>%
  ggplot(aes(x=angle,fill=as.character(inMask)))+geom_histogram(breaks=seq(0, 350, 10),aes(y=..density..),
                                                                alpha=0.5,position='identity')+ylim(0,0.008)+
  scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=8)
p
ggsave(p,filename = "/media/DATA/Maarten/OneDrive/Documents/Manuscripts/in preparation - MFM BRCA2 tracking/Figure 2/angle_histogram_in_out_state0_fast.pdf",width = 10,height = 10,units = "cm")

inputdata <- segs_nest %>%
  filter(displacement1>0.1,angle>0,inMask==T,state<2,condition=="WT MMC")
p <- inputdata %>%
  dplyr::mutate(angle=360-angle)%>%
  rbind(inputdata)%>%
  ggplot(aes(x=angle,fill=as.character(state)))+geom_histogram(breaks=seq(0, 350, 10),aes(y=..density..),
                                                                alpha=0.5,position='identity')+ylim(0,0.008)+
  scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=8)

p

#making linear line plot of fold anisotropy vs mean displacement

break_size <-0.015
data <- segs_nest %>%
  filter(angle>0,state<2,condition=="WT MMC",inMask==T)
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
  filter(angle>0,state<2,condition=="WT MMC",inMask==F)
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
  filter(angle>0,state<2,condition=="WT MMC",inMask_gt1==T)
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
  filter(angle>0,state<2,condition=="WT MMC",inMask_gt2==T)
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
  filter(angle>0,state<2,condition=="WT MMC",inMask_gt3==T)
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
p <- ggplot(results,aes(x=mean_disp,y=fold,color=location))+geom_line()+geom_point()+xlab("Mean displacement (um)")+ylab("Fold anisotropy")+ylim(0,7)+xlim(0,0.5)+scale_colour_Publication()+theme_Publication(base_size=12)+
  scale_y_continuous(expand=c(0,0))
p

ggsave(p,filename = "/media/DATA/Maarten/OneDrive/Documents/Manuscripts/in preparation - MFM BRCA2 tracking/Figure 2/angle_lineplot.pdf",width = 10,height = 7,units = "cm")
