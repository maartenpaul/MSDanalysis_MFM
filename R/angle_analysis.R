library(doParallel)
library(reticulate)
library(plyr)


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

ptm <- proc.time()
#initialize cluster
nodes <- detectCores()
cl <- makeCluster(nodes-4)
registerDoParallel(cl)
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
    get_displacements <- function(x,n){
      x$frame  <- x$frame-x$frame[1]+1
      displ <- cbind(rep(x=-1,nrow(x)),rep(x=-1,nrow(x)))
      if(nrow(x)>=(3*n+n-1)){
        #loop over all steps of tracks
        for (k in 1:(nrow(x)-2*n)){
          
          if (is.element(x[k,2]+n,x$frame)&&is.element(x[k,2]+2*n,x$frame)){ #check if point is present in track, this deals with gaps
            x1 <- as.numeric(x[k,c(3,4,6)])
            which_point1 <- which(x[,2]==x[k,2]+n)
            x2 <- as.numeric(x[which_point1,c(3,4,6)])
            which_point2 <- which(x[,2]==x[k,2]+2*n)
            x3 <- as.numeric(x[which_point2,c(3,4,6)])
            displ[which_point1,] <- c( sqrt((x1[1]-x2[1])^2+(x1[2]-x2[2])^2) ,sqrt((x2[1]-x3[1])^2+(x2[2]-x3[2])^2) )
            
          }}
      }
      return(displ)
    }
    
    
    
    ddply(x,.variables = c("track"), .parallel = F, function(x){
      
      
      
      x$angle1 <- get_angles(x,1)
      x[c("displacement1","displacement2")] <- get_displacements(x,1)
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
  inputdata <- subset(segments_all[[k]],displacement>0.1&angle>0&state<2)
  
  fold <- tibble("time"=numeric(),"location"=character(),"fold"=numeric())
  
                
  for (i in 1:10){
    data3 <- inputdata[inputdata$inMask==1,]
    datasub <- data3[paste0("angle",i)][data3[paste0("angle",i)]>0]
     fold <- add_row(fold,"time"=0.05*i,location='inside',"fold"= length(datasub[datasub>165])/length(datasub[datasub<15]))
    
     data3 <- inputdata[inputdata$inMask_gt1==1,]
     datasub <- data3[paste0("angle",i)][data3[paste0("angle",i)]>0]
     fold <- add_row(fold,"time"=0.05*i,location='inside_gt1',"fold"= length(datasub[datasub>165])/length(datasub[datasub<15]))
     
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
        theme_Publication(base_size=18)+ + labs(x=NULL, y=NULL)+
        theme(legend.position = "none",text = element_text(size=15),axis.text.x = element_text(hjust = -2),axis.text.y=element_blank(),axis.ticks=element_blank(),line = element_blank(),panel.border=element_blank())+
        scale_colour_Publication()
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

#spatial angle analysis


segments_all$`WT MMC`$displacement <- sqrt(segments_all$`WT MMC`$step_x^2+segments_all$`WT MMC`$step_y^2)

data <- subset(segments_all$`WT MMC`,angle>0&state<2)
bins <- seq(0,0.5,0.02)

#results <- matrix(nrow = length(bins),ncol = length(bins))
results <- data.frame("disp1"=rep(bins,each=length(bins)),"disp2"=rep(bins,times=length(bins)),"fold"=0)
results <- ddply(results,.variables = c("disp1","disp2"),function(x){
    datasub <- subset(data,displacement1>=x$disp1&displacement1<(x$disp1+0.01)&displacement2>=x$disp2&displacement2<(x$disp2+0.01))
    x$fold <- length(datasub$angle1[datasub$angle1>165])/length(datasub$angle1[datasub$angle1<15])
    return(x)
})
results$fold[is.na(results$fold)] <- 0
results$fold[is.infinite(results$fold)] <- 0


ggplot(data = results, aes(x=disp1, y=disp2, fill=fold)) + 
  geom_tile()+
  scale_fill_distiller(palette="YlOrBr",direction = 1) + theme_bw()

#repeat for different ML states
#repeat for filtered in and out

#making linear line plot of mean displacement
data <- subset(segments_all$`WT MMC`,angle>0&state<2)
data$mean_disp <- (data$displacement1+data$displacement2)/2

bins <- seq(0,0.5,0.02)

results <- data.frame("mean_disp"=bins,"fold"=2)
results <- ddply(results,.variables = c("mean_disp"),function(x){
  
  datasub <- subset(data,mean_disp>=x$mean_disp&mean_disp<(x$mean_disp+0.01))
  x$fold <- length(datasub$angle1[datasub$angle1>165])/length(datasub$angle1[datasub$angle1<15])
  return(x)
})

ggplot(results,aes(x=mean_disp,y=fold))+geom_line()+geom_point()+xlab("Mean displacement (um)")+ylab("Fold anisotropy")+ylim(0,5)

segments_all$`WT MMC`$logD_ML <- log10(segments_all$`WT MMC`$D_ML*100)
data <- subset(segments_all$`WT MMC`,angle>0&state!=2&inMask==F)

bins_logD <- seq(-4,1,length.out = 25)
binsize_logD <- bins_logD[2]-bins_logD[1]
bins_Smss <- seq(0,1,length.out = 25)
binsize_Smss <- bins_Smss[2]-bins_Smss[1]
results <- data.frame("D_ML"=rep(bins_logD,each=length(bins)),"Smss_ML"=rep(bins_Smss,times=length(bins)),"fold"=0)
results <- ddply(results,.variables = c("D_ML","Smss_ML"),function(x){
  datasub <- subset(data,logD_ML>=x$D_ML&logD_ML<(x$D_ML+binsize_logD)&Smss_ML>=x$Smss_ML&Smss_ML<(x$Smss_ML+binsize_Smss))
  x$fold <- length(datasub$angle1[datasub$angle1>165])/length(datasub$angle1[datasub$angle1<15])
  return(x)
})
results$fold[is.na(results$fold)] <- 0
results$fold[is.infinite(results$fold)] <- 0


ggplot(data = results, aes(x=D_ML, y=Smss_ML, fill=fold)) + 
  geom_tile()+
  scale_fill_distiller(palette="YlOrBr",direction = 1) + theme_bw()+xlab("log10 Dapp")+ylab("Smss")
