library(dplyr)
library(gridExtra)
segments <- as_tibble(segments)

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


# plot fractions of all segments --------------------------------------------------------------------
for (cond in condition_list) {
  p <- rbind(segments %>%
        filter(condition==cond) %>%
        group_by(cellID)%>%
        summarise(n=length(unique(track))) %>%
        filter(n<10000) %>%
        left_join(.,segments,by = "cellID") %>%
        group_by(cellID,track)%>%
        summarise(n=n())%>%
        filter(n>4) %>%
        left_join(.,segments,by = c("cellID","track")) %>%
        group_by(inMask,state)%>%
        dplyr::summarise(n=n()) %>%
        group_by(inMask) %>%
        mutate(fraction=n/sum(n))%>%
        mutate(labels=paste0(c("outside","inside")[inMask+1],"_",c("fast","slow",'immobile')[state+1]))%>%
        mutate(labels=factor(labels,levels=c("outside_immobile","outside_slow","outside_fast","inside_immobile","inside_slow","inside_fast"))),
      segments %>%
        filter(condition==cond) %>%
        group_by(cellID)%>%
        summarise(n=length(unique(track))) %>%
        filter(n<10000) %>%
        left_join(.,segments,by = "cellID") %>%
        group_by(cellID,track)%>%
        summarise(n=n())%>%
        filter(n>4) %>%
        left_join(.,segments,by = c("cellID","track")) %>%
        group_by(inMask_gt1,state)%>%
        dplyr::summarise(n=n()) %>%
        group_by(inMask_gt1) %>%
        mutate(fraction=n/sum(n))%>%
        mutate(labels=paste0(c("outside_gt1","inside_gt1")[inMask_gt1+1],"_",c("fast","slow",'immobile')[state+1]))%>%
        mutate(labels=factor(labels,levels=c("outside_gt1_immobile","outside_gt1_slow","outside_gt1_fast","inside_gt1_immobile","inside_gt1_slow","inside_gt1_fast"))) %>%
        filter(inMask_gt1==1))%>%
   ggplot(aes(y=fraction, x=labels, fill=labels))+geom_bar(stat="identity")+ylim(0,1)+xlab("")+ggtitle(cond)+
   theme_Publication(base_size=25)+scale_colour_Publication()+scale_fill_Publication()+theme(legend.position = "none",text = element_text(size=15),axis.text.x = element_text(angle = 90))
  ggsave(p,filename = file.path(directory,paste0(cond,"_fraction_segments.pdf")))
  ggsave(p,filename = file.path(directory,paste0(cond,"_fraction_segments.png")))
  
} 
# plot fractions of all segments per cell --------------------------------------------------------------------
allplots <- list()
for (cond in condition_list) {
p <- rbind(segments %>%
  filter(condition==cond) %>%
  #filter out cells with more than 20k tracks
    group_by(cellID)%>%
  summarise(n=length(unique(track))) %>%
  filter(n<20000) %>%
  left_join(.,segments,by = "cellID") %>%
  #filter out short tracks
  # group_by(cellID,track)%>%
  # summarise(n=n())%>%
  # filter(n>4) %>%
  # left_join(.,segments,by = c("cellID","track")) %>%
  #calculate statistics for plot
    group_by(cellID,track)%>%
    summarise(n=n())%>%
    filter(n>4) %>%
    left_join(.,segments,by = c("cellID","track")) %>%
  group_by(.id,inMask,state)%>%
  dplyr::summarise(n=n()) %>%
  group_by(.id,inMask) %>%
  mutate(fraction=n/sum(n))%>%
  mutate(labels=paste0(c("outside","inside")[inMask+1],"_",c("fast","slow",'immobile')[state+1]))%>%
  mutate(labels=factor(labels,levels=c("outside_immobile","outside_slow","outside_fast","inside_immobile","inside_slow","inside_fast")))%>%
  group_by(labels,inMask)%>%
  summarise(var=sd(fraction,na.rm = T)/sqrt(n()),fraction=mean(fraction)),
  segments %>%
    filter(condition==cond) %>%
    #filter out cells with more than 20k tracks
    group_by(cellID)%>%
    summarise(n=length(unique(track))) %>%
    filter(n<20000) %>%
    left_join(.,segments,by = "cellID") %>%
    #filter out short tracks
    # group_by(cellID,track)%>%
    # summarise(n=n())%>%
    # filter(n>4) %>%
    # left_join(.,segments,by = c("cellID","track")) %>%
    #calculate statistics for plot
    group_by(cellID,track)%>%
    summarise(n=n())%>%
    filter(n>4) %>%
    left_join(.,segments,by = c("cellID","track")) %>%
    group_by(.id,inMask_gt1,state)%>%
    dplyr::summarise(n=n()) %>%
    group_by(.id,inMask_gt1) %>%
    mutate(fraction=n/sum(n))%>%
    mutate(labels=paste0(c("outside_gt1","inside_gt1")[inMask_gt1+1],"_",c("fast","slow",'immobile')[state+1]))%>%
    mutate(labels=factor(labels,levels=c("outside_gt1_immobile","outside_gt1_slow","outside_gt1_fast","inside_gt1_immobile","inside_gt1_slow","inside_gt1_fast")))%>%
    group_by(labels,inMask_gt1)%>%
    summarise(var=sd(fraction,na.rm = T)/sqrt(n()),fraction=mean(fraction))%>%
    filter(inMask_gt1==1))%>%
  ggplot(aes(y=fraction, x=labels, fill=labels))+geom_bar(stat="identity")+ylim(0,1)+geom_errorbar(stat="identity",
          aes(ymax=fraction+var,ymin=fraction-var))+theme_Publication(base_size=25)+xlab("")+
      scale_colour_Publication()+scale_fill_Publication()+
  theme(legend.position = "none",text = element_text(size=15),axis.text.x = element_text(angle = 90,size=6,hjust=1,vjust=0.5))
  ggsave(p,filename = file.path(directory,paste0(cond,"_mean_fraction_segments.pdf")))
  ggsave(p,filename = file.path(directory,paste0(cond,"_mean_fraction_segments.png")))
  allplots[[cond]] <- p
  }

grid.arrange(allplots[[7]]+ggtitle(names(allplots)[7]),allplots[[6]],allplots[[4]],allplots[[3]],allplots[[2]],allplots[[1]],allplots[[8]],nrow=2)
ggsave(allplots,filename = file.path(directory,"200825_fractions all.pdf"))

# plot fractions of all tracklets --------------------------------------------------------------------
segments %>%
  filter(condition=="WT MMC") %>%
  #filter out cells with more than 20k tracks
  group_by(cellID)%>%
  summarise(n=length(unique(track))) %>%
  filter(n<20000) %>%
  left_join(.,segments,by = "cellID") %>%
  #filter out short tracks
  group_by(cellID,track)%>%
  summarise(n=n())%>%
  filter(n>4) %>%
  left_join(.,segments,by = c("cellID","track")) %>%
  #calculate statistics for plot
  group_by(trackInMask,state)%>%
  dplyr::summarise(n=n_distinct(tracklet)) %>%
  group_by(trackInMask) %>%
  mutate(fraction=n/sum(n))%>%
  mutate(labels=paste0(c("outside","inside")[trackInMask+1],"_",c("fast","slow",'immobile')[state+1]))%>%
  mutate(labels=factor(labels,levels=c("outside_immobile","outside_slow","outside_fast","inside_immobile","inside_slow","inside_fast")))%>%
  ggplot(aes(y=fraction, x=labels, fill=labels))+geom_bar(stat="identity")+ylim(0,1)

install.packages("Rmisc")
library(Rmisc)

# plot fractions of all tracklet per cell per tracklet --------------------------------------------------------------------
data <- segments %>%
  filter(condition=="WT MMC") %>%
  #filter out cells with more than 20k tracks
  group_by(cellID)%>%
  summarise(n=length(unique(track))) %>%
  filter(n<20000) %>%
  left_join(.,segments,by = "cellID") %>%
  #filter out short tracks
  group_by(cellID,track)%>%
  summarise(n=n())%>%
  filter(n>5) %>%
  left_join(.,segments,by = c("cellID","track"))
 

rbind(data %>%
  group_by(.id,trackInMask,state)%>%
  dplyr::summarise(n=n_distinct(tracklet)) %>%
  group_by(.id,trackInMask) %>%
  mutate(fraction=n/sum(n))%>%
  mutate(labels=paste0(c("outside","inside")[trackInMask+1],"_",c("fast","slow",'immobile')[state+1]))%>%
  mutate(labels=factor(labels,levels=c("outside_immobile","outside_slow","outside_fast","inside_immobile","inside_slow","inside_fast")))%>%
  group_by(labels)%>%
  summarise(var=sd(fraction,na.rm = T)/sqrt(n()),fraction=mean(fraction),state=state[1],trackInMask=trackInMask[1]),
  (data%>%
  group_by(.id,trackInMask_gt1,state)%>%
  dplyr::summarise(n=n_distinct(tracklet)) %>%
  group_by(.id,trackInMask_gt1) %>%
  mutate(fraction=n/sum(n))%>%
  mutate(labels=paste0(c("outside","inside_GT")[trackInMask_gt1+1],"_",c("fast","slow",'immobile')[state+1]))%>%
  mutate(labels=factor(labels,levels=c("outside_immobile","outside_slow","outside_fast","inside_GT_immobile","inside_GT_slow","inside_GT_fast")))%>%
  group_by(labels)%>%
  summarise(var=sd(fraction,na.rm = T)/sqrt(n()),fraction=mean(fraction),state=state[1],trackInMask=trackInMask_gt1[1]))[4:6,] ) %>%
  ggplot(aes(x=labels, fill=labels))+geom_bar(stat="identity",aes(y=fraction))+xlab("")+
  ylim(0,1)+geom_errorbar(stat="identity",aes(ymax=fraction+var,ymin=fraction-var)) +theme_Publication(base_size=25)+scale_colour_Publication()+scale_fill_Publication()+theme(legend.position = "none",text = element_text(size=15),axis.text.x = element_text(angle = 90))

#plot average lengths of the different states?

#plot fraction of tracklets inside mask vs outside
library(reshape2)
library(ggbeeswarm)


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
x<- data %>%
  select(cellID,track,trackInMask)%>%
  distinct()%>%
  group_by(cellID)%>%
  summarise(n=n(),fraction=sum(trackInMask)/n())
mean(x$fraction)
x %>%
  ggplot(aes(y=fraction))+ geom_boxplot()+ylim(0,1)



data %>%
  group_by(cellID)%>%
  summarise(n=n(),data=sum(inMask)/n(),random=sum(inMask_gt1)/n())%>%
  melt(id.vars="cellID",measure.vars=c("data","random"))%>%
  ggplot(aes(y=value,x=variable,fill=variable))+ geom_boxplot(outlier.shape = NA,notch = T)+geom_quasirandom()+ylim(0,0.5)+
  scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=16)+ theme(legend.position = "none")+
  xlab("")+ylab("fraction localizations inside mask")

data %>%
  group_by(cellID)%>%
  summarise(n=n(),fraction=sum(inMask)/n(),fraction_gt1=sum(inMask_gt1)/n())%>%
  ggplot(aes(y=fraction_gt1))+ geom_boxplot()+ylim(0,0.5)


data %>%
  summarise()
