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
  dplyr::summarise(n=n(),fraction=sum(trackInMask)/n())
mean(x$fraction)
x %>%
  ggplot(aes(y=fraction))+ geom_boxplot()+ylim(0,1)



data %>%
  group_by(cellID)%>%
  dplyr::summarise(n=n(),data=sum(inMask)/n(),random=sum(inMask_gt1)/n())%>%
  melt(id.vars="cellID",measure.vars=c("data","random"))%>%
  ggplot(aes(y=value,x=variable,fill=variable))+ geom_boxplot(outlier.shape = NA,notch = T)+geom_quasirandom()+ylim(0,0.5)+
  scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=16)+ theme(legend.position = "none")+
  xlab("")+ylab("fraction localizations inside mask")

data %>%
  group_by(cellID)%>%
  dplyr::summarise(n=n(),fraction=sum(inMask)/n(),fraction_gt1=sum(inMask_gt1)/n())%>%
  ggplot(aes(y=fraction_gt1))+ geom_boxplot()+ylim(0,0.5)

