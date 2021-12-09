data <- segments_all$`WT MMC`
library(doParallel)
library(reticulate)
library(plyr)
ptm <- proc.time()
#initialize cluster
nodes <- detectCores()
cl <- makeCluster(nodes-4)
registerDoParallel(cl)

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

for (j in 1:length(segments_all)){
  segments_all[[j]] <- ddply(segments_all[[j]],.variables = c("cellID"), .parallel = T, function(x){
  get_angles <- function(x,n){
    x$frame  <- x$frame-x$frame[1]+1
    angles <- rep(x=-1,nrow(x))
    if(nrow(x)>=(3*n+n-1)){
      get_angle_3D <- function(A,B,C){
        seg_angle <- vector()
        AB <- B[1:2]-A[1:2]
        CB <- C[1:2]-B[1:2]
        
        #dAB <- sqrt((B[1]-A[1])^2+(B[2]-A[2])^2)
        #dBC <- sqrt((C[1]-B[1])^2+(C[2]-B[2])^2)
        #Formula obtained from https://gitlab.com/anders.sejr.hansen/anisotropy
        angle <- abs(atan2(det(cbind(AB,CB)),AB%*%CB))
        angle <- angle/pi*180
        return(angle)
        
        
      } 
      #loop over all steps of tracks
      for (k in 1:(nrow(x)-2*n)){
        
        if (is.element(x[k,2]+n,x$frame)&&is.element(x[k,2]+2*n,x$frame)){ #check if point is present in track, this deals with gaps
          x1 <- as.numeric(x[k,c(3,4,6)])
          which_point1 <- which(x[,2]==x[k,2]+n)
          x2 <- as.numeric(x[which_point1,c(3,4,6)])
          which_point2 <- which(x[,2]==x[k,2]+2*n)
          x3 <- as.numeric(x[which_point2,c(3,4,6)])
          angles[which_point1] <- get_angle_3D(x1,x2,x3)
          
        }}
    }
    return(angles)
  }
  
  ddply(x,.variables = c("track"), .parallel = F, function(x){
  
 
  
  x$angle1 <- get_angles(x,1)
  x$angle2 <- get_angles(x,2)
  x$angle3 <- get_angles(x,3)
  x$angle4 <- get_angles(x,4)
  x$angle5 <- get_angles(x,5)
  x$angle6 <- get_angles(x,6)
  x$angle7 <- get_angles(x,7)
  x$angle8 <- get_angles(x,8)
  x$angle9 <- get_angles(x,9)
  x$angle10 <- get_angles(x,10)
  x$angle11 <- get_angles(x,11)
  x$angle12 <- get_angles(x,12)
  x$angle13 <- get_angles(x,13)
  x$angle14 <- get_angles(x,14)
  x$angle15 <- get_angles(x,15)
  x$angle16 <- get_angles(x,16)
  x$angle17 <- get_angles(x,17)
  x$angle18 <- get_angles(x,18)
  x$angle19 <- get_angles(x,19)
  x$angle20 <- get_angles(x,20)
  
    
    return(x)
    
})
  
})
}
stopCluster(cl)
proc.time() - ptm

save(segments_all,file=file.path(directory,"segments_all_angles.Rdata"))

#make a plot
library(tidyverse)
 
for (k in 1:length(segments_all)){
  segments_all[[k]]$displacement <- sqrt(segments_all[[k]]$step_x^2+segments_all[[k]]$step_y^2)
  inputdata <- subset(segments_all[[k]],displacement>0.2&angle>0&state<2)
  
  fold <- tibble("time"=numeric(),"location"=character(),"fold"=numeric())
  
                
  for (i in 1:10){
    data3 <- inputdata[inputdata$inMask==1,]
    datasub <- data3[paste0("angle",i)][data3[paste0("angle",i)]>0]
     fold <- add_row(fold,"time"=0.05*i,location='inside',"fold"= length(datasub[datasub>165])/length(datasub[datasub<15]))
    
     
     data3 <- inputdata[inputdata$inMask==0,]
    datasub <- data3[paste0("angle",i)][data3[paste0("angle",i)]>0]
    fold <- add_row(fold,"time"=0.05*i,location='outside',"fold"= length(datasub[datasub>165])/length(datasub[datasub<15]))  
    
    if(i==1){
      angles <- hist(c(datasub,abs(360-datasub)),breaks = seq(0,360,15),plot=F)
      
      angles <- data.frame('mids'=angles$mids,'density'=angles$density)
      
      ggplot(angles,aes(x = mids,y=density,fill=density))+
        geom_hline(yintercept = seq(0, 0.01, by = 0.002), colour = "black", size = 0.2) +
        geom_vline(xintercept = seq(0, 360-1, by = 45), colour = "black", size = 0.2) +
        scale_x_continuous(breaks=seq(0, 350, 45))+coord_polar(start = 0.5*pi,direction = 1)+
        geom_bar(stat='identity',width=16) +ylim(c(0,.01))+
        theme_Publication(base_size=18)+ 
        theme(legend.position = "none",text = element_text(size=15),axis.text.x = element_text(hjust = -2),axis.text.y=element_blank(),axis.ticks=element_blank())+scale_colour_Publication()
      #scale_fill_Publication()
    }
    
  }


  p <- ggplot(fold,aes(x=time,y=fold,col=location))+geom_line()+xlab("Interval (s)")+ylab("fold difference 180/0 degrees")+ 
    ylim(0,6)+scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=18)
 print(p)
   ggsave(p,filename = file.path(directory,names(segments_all)[k],"anglefoldplot.png"))
   Sys.sleep(10)
  ggsave(p,filename = file.path(directory,names(segments_all)[k],"anglefoldplot.pdf"))
  
  
  }
      