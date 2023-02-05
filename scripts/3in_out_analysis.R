#in out analysis
library(tidyverse)

#state0=fast;state1=slow;state2=immobile
#input variables
directory <- "/media/DATA/Maarten/data_gtv2/"

condition_list <- list.dirs(directory,full.names = F,recursive = F)
#condition_list <- condition_list[c(2,3,5)]

load(file=file.path(directory,"segs_nest.Rdata"))

segments <- dplyr::filter(segs_nest,condition=="WT MMC")

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
  discrete_scale("fill","Publication",manual_pal(values = c("#c00000","#1F497D","#542788","#217D68","#386cb0","#fdb462","#386cb0","#fdb462","#386cb0","#fdb462","#386cb0","#fdb462")), ...)
  
}
scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#c00000","#6599D9","#542788","#217D68","#386cb0","#fdb462","#386cb0","#fdb462","#386cb0","#fdb462","#386cb0","#fdb462")), ...)
  
}
data <- segments %>%
  filter(condition=="WT MMC")%>%
  group_by(cellID,track)%>%
  dplyr::summarise(n=n())%>%
  filter(n>5) %>%
  left_join(.,segments,by = c("cellID","track"))%>%
  group_by(cellID,track)%>%
  dplyr::mutate(inMask_shift=paste0(c(inMask[2:length(inMask)],-1),c(inMask[3:length(inMask)],-1,-1)),
         inMask_gt1_shift=paste0(c(inMask_gt1[2:length(inMask_gt1)],-1),c(inMask_gt1[3:length(inMask_gt1)],-1,-1)),
         inMask_gt2_shift=paste0(c(inMask_gt2[2:length(inMask_gt2)],-1),c(inMask_gt2[3:length(inMask_gt2)],-1,-1)),
         inMask_gt3_shift=paste0(c(inMask_gt3[2:length(inMask_gt3)],-1),c(inMask_gt3[3:length(inMask_gt3)],-1,-1)))

hist(data$distMask)

#out -> in
real <- data%>% 
  ungroup() %>%
  filter(distMask!=0,distMask<0.5,state!=2,inMask_shift!="-1-1"&inMask_shift!="0-1"&inMask_shift!="1-1")%>%
  group_by(cellID)%>%
  dplyr::summarise(prob=length(inMask_shift[inMask_shift==11])/length(inMask_shift))

real$data <- "data"

gt1 <- data%>% 
  ungroup() %>%
  filter(distMask_gt1!=0,distMask_gt1<0.5,state!=2,inMask_gt1_shift!="-1-1"&inMask_gt1_shift!="0-1"&inMask_gt1_shift!="1-1")%>%
  group_by(cellID)%>%
  #at least two frames inside
  dplyr::summarise(prob=length(inMask_gt1_shift[inMask_gt1_shift==11])/length(inMask_gt1_shift))

gt1$data <- "gt1"

gt2 <- data%>% 
  ungroup() %>%
  filter(distMask_gt2!=0,distMask_gt2<0.5,state!=2,inMask_gt2_shift!="-1-1"&inMask_gt2_shift!="0-1"&inMask_gt2_shift!="1-1")%>%
  group_by(cellID)%>%
  dplyr::summarise(prob=length(inMask_gt2_shift[inMask_gt2_shift==11])/length(inMask_gt2_shift))
gt2$data <- "gt2"


gt3 <- data%>% 
  ungroup() %>%
  filter(distMask_gt3!=0,distMask_gt3<0.5,state!=2,inMask_gt3_shift!="-1-1"&inMask_gt3_shift!="0-1"&inMask_gt3_shift!="1-1")%>%
  group_by(cellID)%>%
  dplyr::summarise(prob=length(inMask_gt3_shift[inMask_gt3_shift==11])/length(inMask_gt3_shift))
gt3$data <- "gt3"

plotdata <- rbind(real,gt1,gt2,gt3)
plotdata$data <- factor(x = plotdata$data,levels=c("data","gt1","gt2","gt3","random"))
plotdata$data[plotdata$data!="data"] <- "random"

p <-plotdata%>%
 ggplot(aes(y=prob,x=data,fill=data))+geom_boxplot() + 
  scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=20)+ ylim(0,0.1)+
  theme(legend.position = "none")+xlab("")+ylab("probability")
p
ggsave(p,filename = "/media/DATA/Maarten/OneDrive/Documents/Manuscripts/in preparation - MFM BRCA2 tracking/Figure 2/out_in_probability.pdf",width = 10,height = 10,units = "cm")
(real$prob-mean(c(gt1$prob,gt2$prob,gt3$prob)))/sd(c(gt1$prob,gt2$prob,gt3$prob))

#in -> out all

real <- data%>% 
  ungroup() %>%
  filter(distMask==0,state!=2,inMask_shift!="-1-1"&inMask_shift!="0-1"&inMask_shift!="1-1")%>%
  group_by(cellID)%>%
  dplyr::summarise(prob=length(inMask_shift[inMask_shift=="00"])/length(inMask_shift))
real$data <- "data"

gt1 <- data%>% 
  ungroup() %>%
  filter(distMask_gt1==0,state!=2,inMask_gt1_shift!="-1-1"&inMask_gt1_shift!="0-1"&inMask_gt1_shift!="1-1")%>%
  group_by(cellID)%>%
  dplyr::summarise(prob=length(inMask_gt1_shift[inMask_gt1_shift=="00"])/length(inMask_gt1_shift))
gt1$data <- "gt1"

gt2 <- data%>% 
  ungroup() %>%
  filter(distMask_gt2==0,state!=2,inMask_gt2_shift!="-1-1"&inMask_gt2_shift!="0-1"&inMask_gt2_shift!="1-1")%>%
  group_by(cellID)%>%
  dplyr::summarise(prob=length(inMask_gt2_shift[inMask_gt2_shift=="00"])/length(inMask_gt2_shift))
gt2$data <- "gt2"

gt3 <- data%>% 
  ungroup() %>%
  filter(distMask_gt3==0,state!=2,inMask_gt3_shift!="-1-1"&inMask_gt3_shift!="0-1"&inMask_gt3_shift!="1-1")%>%
  group_by(cellID)%>%
  dplyr::summarise(prob=length(inMask_gt3_shift[inMask_gt3_shift=="00"])/length(inMask_gt3_shift))
gt3$data <- "gt3"


plotdata <- rbind(real,gt1,gt2,gt3)
plotdata$data <- factor(x = plotdata$data,levels=c("data","gt1","gt2","gt3","random"))
plotdata$data[plotdata$data!="data"] <- "random"

p <- ggplot(plotdata,aes(y=prob,x=data,fill=data))+geom_boxplot()+ 
  scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=20)+ylim(0,1)+
  theme(legend.position = "none")+xlab("")+ylab("probability")
p
ggsave(p,filename = "/media/DATA/Maarten/OneDrive/Documents/Manuscripts/in preparation - MFM BRCA2 tracking/Figure 2/in_out_probability.pdf",width = 10,height = 10,units = "cm")


#in -> out slow

real <- data%>% 
  ungroup() %>%
  filter(distMask==0,state==1,inMask_shift!="-1-1"&inMask_shift!="0-1"&inMask_shift!="1-1")%>%
  group_by(cellID)%>%
  dplyr::summarise(prob=length(inMask_shift[inMask_shift=="00"])/length(inMask_shift))
real$data <- "data"

gt1 <- data%>% 
  ungroup() %>%
  filter(distMask_gt1==0,state==1,inMask_gt1_shift!="-1-1"&inMask_gt1_shift!="0-1"&inMask_gt1_shift!="1-1")%>%
  group_by(cellID)%>%
  dplyr::summarise(prob=length(inMask_gt1_shift[inMask_gt1_shift=="00"])/length(inMask_gt1_shift))
gt1$data <- "gt1"

gt2 <- data%>% 
  ungroup() %>%
  filter(distMask_gt2==0,state==1,inMask_gt2_shift!="-1-1"&inMask_gt2_shift!="0-1"&inMask_gt2_shift!="1-1")%>%
  group_by(cellID)%>%
  dplyr::summarise(prob=length(inMask_gt2_shift[inMask_gt2_shift=="00"])/length(inMask_gt2_shift))
gt2$data <- "gt2"

gt3 <- data%>% 
  ungroup() %>%
  filter(distMask_gt3==0,state==1,inMask_gt3_shift!="-1-1"&inMask_gt3_shift!="0-1"&inMask_gt3_shift!="1-1")%>%
  group_by(cellID)%>%
  dplyr::summarise(prob=length(inMask_gt3_shift[inMask_gt3_shift=="00"])/length(inMask_gt3_shift))
gt3$data <- "gt3"

plotdata <- rbind(real,gt1,gt2,gt3)
plotdata$data <- factor(x = plotdata$data,levels=c("data","gt1","gt2","gt3","random"))
plotdata$data[plotdata$data!="data"] <- "random"


p <- ggplot(plotdata,aes(y=prob,x=data,fill=data))+geom_boxplot()+ 
  scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=20)+ theme(legend.position = "none")+ylim(0,1)+xlab("")+ylab("probability")
p
ggsave(p,filename = "/media/DATA/Maarten/OneDrive/Documents/Manuscripts/in preparation - MFM BRCA2 tracking/Figure 2/in_out_probability_slow.pdf",width = 10,height = 10,units = "cm")

#in -> out fast

real <- data%>% 
  ungroup() %>%
  filter(distMask==0,state==0,inMask_shift!="-1-1"&inMask_shift!="0-1"&inMask_shift!="1-1")%>%
  group_by(cellID)%>%
  dplyr::summarise(prob=length(inMask_shift[inMask_shift=="00"])/length(inMask_shift))
real$data <- "data"

gt1 <- data%>% 
  ungroup() %>%
  filter(distMask_gt1==0,state==0,inMask_gt1_shift!="-1-1"&inMask_gt1_shift!="0-1"&inMask_gt1_shift!="1-1")%>%
  group_by(cellID)%>%
  dplyr::summarise(prob=length(inMask_gt1_shift[inMask_gt1_shift=="00"])/length(inMask_gt1_shift))
gt1$data <- "gt1"


gt2 <- data%>% 
  ungroup() %>%
  filter(distMask_gt2==0,state==0,inMask_gt2_shift!="-1-1"&inMask_gt2_shift!="0-1"&inMask_gt2_shift!="1-1")%>%
  group_by(cellID)%>%
  dplyr::summarise(prob=length(inMask_gt2_shift[inMask_gt2_shift=="00"])/length(inMask_gt2_shift))
gt2$data <- "gt2"

gt3 <- data%>% 
  ungroup() %>%
  filter(distMask_gt3==0,state==0,inMask_gt3_shift!="-1-1"&inMask_gt3_shift!="0-1"&inMask_gt3_shift!="1-1")%>%
  group_by(cellID)%>%
  dplyr::summarise(prob=length(inMask_gt3_shift[inMask_gt3_shift=="00"])/length(inMask_gt3_shift))
gt3$data <- "gt3"

plotdata <- rbind(real,gt1,gt2,gt3)
plotdata$data <- factor(x = plotdata$data,levels=c("data","gt1","gt2","gt3","random"))
plotdata$data[plotdata$data!="data"] <- "random"

p <- ggplot(plotdata,aes(y=prob,x=data,fill=data))+geom_boxplot()+ 
  scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=20)+ theme(legend.position = "none")+xlab("")+ylab("probability")+ylim(0,1)
p
ggsave(p,filename = "/media/DATA/Maarten/OneDrive/Documents/Manuscripts/in preparation - MFM BRCA2 tracking/Figure 2/in_out_probability_fast.pdf",width = 10,height = 10,units = "cm")





