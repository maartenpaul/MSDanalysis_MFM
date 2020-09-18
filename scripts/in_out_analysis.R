#in out analysis

data <- segments %>%
  filter(condition=="WT MMC")%>%
  group_by(cellID,track)%>%
  summarise(n=n())%>%
  filter(n>5) %>%
  left_join(.,segments,by = c("cellID","track"))%>%
  group_by(cellID,track)%>%
  mutate(inMask_shift=paste0(c(inMask[2:length(inMask)],-1),c(inMask[3:length(inMask)],-1,-1)),
         inMask_gt1_shift=paste0(c(inMask_gt1[2:length(inMask_gt1)],-1),c(inMask_gt1[3:length(inMask_gt1)],-1,-1)),
         inMask_gt2_shift=paste0(c(inMask_gt2[2:length(inMask_gt2)],-1),c(inMask_gt2[3:length(inMask_gt2)],-1,-1)),
         inMask_gt3_shift=paste0(c(inMask_gt3[2:length(inMask_gt3)],-1),c(inMask_gt3[3:length(inMask_gt3)],-1,-1)))

hist(data$distMask)

#out -> in
real <- data%>% 
  ungroup() %>%
  filter(distMask!=0,distMask<0.5,state!=2,inMask_shift!="-1-1"&inMask_shift!="0-1"&inMask_shift!="1-1")%>%
  group_by(cellID)%>%
  summarise(prob=length(inMask_shift[inMask_shift==11])/length(inMask_shift))

real$data <- "real"

gt1 <- data%>% 
  ungroup() %>%
  filter(distMask_gt1!=0,distMask_gt1<0.5,state!=2,inMask_gt1_shift!="-1-1"&inMask_gt1_shift!="0-1"&inMask_gt1_shift!="1-1")%>%
  group_by(cellID)%>%
  #at least two frames inside
  summarise(prob=length(inMask_gt1_shift[inMask_gt1_shift==11])/length(inMask_gt1_shift))

gt1$data <- "gt1"

gt2 <- data%>% 
  ungroup() %>%
  filter(distMask_gt2!=0,distMask_gt2<0.5,state!=2,inMask_gt2_shift!="-1-1"&inMask_gt2_shift!="0-1"&inMask_gt2_shift!="1-1")%>%
  group_by(cellID)%>%
  summarise(prob=length(inMask_gt2_shift[inMask_gt2_shift==11])/length(inMask_gt2_shift))
gt2$data <- "gt2"


gt3 <- data%>% 
  ungroup() %>%
  filter(distMask_gt3!=0,distMask_gt3<0.5,state!=2,inMask_gt3_shift!="-1-1"&inMask_gt3_shift!="0-1"&inMask_gt3_shift!="1-1")%>%
  group_by(cellID)%>%
  summarise(prob=length(inMask_gt3_shift[inMask_gt3_shift==11])/length(inMask_gt3_shift))
gt3$data <- "gt3"

plotdata <- rbind(real,gt1,gt2,gt3)
plotdata$data <- factor(x = plotdata$data,levels=c("real","gt1","gt2","gt3"))

ggplot(plotdata,aes(y=prob,fill=data))+geom_boxplot()

(real$prob-mean(c(gt1$prob,gt2$prob,gt3$prob)))/sd(c(gt1$prob,gt2$prob,gt3$prob))

#in -> out

real <- data%>% 
  ungroup() %>%
  filter(distMask==0,state!=2,inMask_shift!="-1-1"&inMask_shift!="0-1"&inMask_shift!="1-1")%>%
  group_by(cellID)%>%
  summarise(prob=length(inMask_shift[inMask_shift=="00"])/length(inMask_shift))
real$data <- "real"

gt1 <- data%>% 
  ungroup() %>%
  filter(distMask_gt1==0,state!=2,inMask_gt1_shift!="-1-1"&inMask_gt1_shift!="0-1"&inMask_gt1_shift!="1-1")%>%
  group_by(cellID)%>%
  summarise(prob=length(inMask_gt1_shift[inMask_gt1_shift=="00"])/length(inMask_gt1_shift))
gt1$data <- "gt1"


plotdata <- rbind(real,gt1)
plotdata$data <- factor(x = plotdata$data,levels=c("real","gt1","gt2","gt3"))

ggplot(plotdata,aes(y=prob,fill=data))+geom_boxplot()

