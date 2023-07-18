msd_analyze_data_mosaic_mask_parallel_intensity <- function(directory,condition_list,framerate,n,fitzero,pixelsize,fitMSD,offset,max_tracks,track_file_name="tracks.simple.filtered.txt",extension="",dim,groundtruth=TRUE){
  segments_all <- list()
  msd_fit_all <- list()
  track_stats_all <- list()
  library(readr)
  for (i in 1:length(condition_list)){
    dirs <- file.path(directory,condition_list[i])
    filelist <- list.files(dirs,full.names = T,recursive = F,pattern = "^Traj_.*.\\mask2.csv$")
    filelist <- filelist[-grep(filelist,pattern = "^.*transformed.*$")]
    total <- seq(1:length(filelist))
    
    nodes <- detectCores()
    cl <- makeCluster(nodes-12)
    registerDoParallel(cl)
    dimensions <- dim
    offst <- offset
    pixelsize <- pixelsize
    n <- n
    fitMSD<- fitMSD
    fitzero <-fitzero
    framerate<-framerate
    groundtruth <- groundtruth
    cat(paste0("Analyzing ",condition_list[i],"\n"))
    output <-alply(filelist,.margins = 1,.parallel = TRUE,.paropts = list(.export=c("filelist","dirs","groundtruth","n","framerate","offst","pixelsize","dimensions","fitMSD","fitzero"),.packages=c("MSDtracking","readr","stats")),function(j){
      Sys.sleep(0.1)
      tracks_simple <- read_csv(j)
      names(tracks_simple) <- c("track","X","Y","Z","time","frame","step_x","step_y","step_z","inMask","rawInt","meanInt","normInt", "sdInt","maxInt","minInt","distMask","invDistMask", "label", "center_x", "center_y", "center_z", "distanceCentroid3D", "distanceCentroid2D", "distanceToNearest3D", "distanceToNearest2D")
      
      tracks_simple <- tracks_simple[,c(6,2,3,1,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26)]
      tracks_simple$frame <- tracks_simple$frame
      segments <- data.frame(SEGMENT_STAT(tracks_simple),"cellID"=basename(j))
      segments <- ddply(segments,.variables = "track",function(x){
          if(any(x$inMask==TRUE)){
            x$trackInMask <-TRUE
          } else {
            x$trackInMask <-FALSE
          }
          return(x)
        })
      
      
      if(groundtruth){
        filelist_gt <- list.files(dirs,full.names = T,recursive = F,pattern = paste0("^Traj_transformed.*.",basename(j)))
        #cat(paste0("number groundtruths found:", k))
        for(k in 1:length(filelist_gt)){
          tracks_simple_gt <- read_csv(filelist_gt[k])
          names(tracks_simple) <- c("track","X","Y","Z","time","frame","step_x","step_y","step_z","inMask","rawInt","meanInt","normInt", "sdInt","maxInt","minInt","distMask","invDistMask", "label", "center_x", "center_y", "center_z", "distanceCentroid3D", "distanceCentroid2D", "distanceToNearest3D", "distanceToNearest2D")
          segments[[paste0("inMask_gt",k)]] <- tracks_simple_gt$inMask
          segments[[paste0("rawInt_gt",k)]] <- tracks_simple_gt$rawInt
          segments[[paste0("meanInt_gt",k)]] <- tracks_simple_gt$meanInt
          segments[[paste0("sdInt_gt",k)]] <- tracks_simple_gt$sdInt
          segments[[paste0("maxInt_gt",k)]] <- tracks_simple_gt$maxInt
          segments[[paste0("minInt_gt",k)]] <- tracks_simple_gt$minInt
          segments[[paste0("distMask_gt",k)]] <- tracks_simple_gt$distMask
          segments[[paste0("invDistMask_gt",k)]] <- tracks_simple_gt$invDistMask
          segments[[paste0("center_x_gt",k)]] <- tracks_simple_gt$center_x
          segments[[paste0("center_y_gt",k)]] <- tracks_simple_gt$center_y
          segments[[paste0("center_z_gt",k)]] <- tracks_simple_gt$center_z
          segments[[paste0("distanceCentroid3D_gt",k)]] <- tracks_simple_gt$distanceCentroid3D
          segments[[paste0("distanceCentroid2D_gt",k)]] <- tracks_simple_gt$distanceCentroid2D
          segments[[paste0("distanceToNearest3D_gt",k)]] <- tracks_simple_gt$distanceToNearest3D
          segments[[paste0("distanceToNearest2D_gt",k)]] <- tracks_simple_gt$distanceToNearest2D
          segments <- ddply(segments,.variables = "track",function(x){
            if(any(x[[paste0("inMask_gt",k)]]==TRUE)){
              x[[paste0("trackInMask_gt",k)]] <-TRUE
            } else {
              x[[paste0("trackInMask_gt",k)]] <-FALSE
            }
            return(x)
          })
          
        }
        
      }
      track_msd <- TRACK_MSD(segments,n = n,framerate=framerate,pxsize = pixelsize,dim=dimensions)
      if(fitMSD==TRUE){
        #tracks <-  TRACK_MSD_fit(track_msd,n = n,fitzero = fitzero,framerate=framerate,pxsize = pixelsize,offset=offset,dim=dimensions)
        tracks <-  TRACK_MSD_fit(track_msd,n = n,fitzero = fitzero,framerate=framerate,pxsize = pixelsize,dim=dimensions)
        
        for (k in 1:length(tracks$track)){
          tracks$inMask = FALSE
          if(groundtruth){
            for(l in 1:length(filelist_gt)){
              tracks[[paste0("inMask_gt",l)]] = FALSE
            }
          }
        }
        
        for(l in 1:length(tracks$track)){
          tracks$inMask[l] <- any(segments$inMask[segments$track==tracks$track[l]]==1)
        }
        if(groundtruth){
          for(l in 1:length(filelist_gt)){
            
            for(m in 1:length(tracks$track)){
              tracks[[paste0("inMask_gt",l)]][m] <- any(segments[[paste0("inMask_gt",l)]][segments$track==tracks$track[m]]==1)
            }
          }
        }
      }
      if(fitMSD==TRUE){
        return(list("segments"=segments,"tracks"=tracks))
      }
      else {
        return(segments)
      }
    })
    stopCluster(cl)
    segments <- list()
   # tracks <- list()
    for(m in 1:length(output)){
      segments[[m]] <- output[[m]]
    #  tracks[[m]] <- output[[m]]$tracks
    }
    names(segments) <- basename(filelist)
    #names(tracks) <- basename(filelist)
    output <- NULL
    # stats <- TRACK_STAT(x=segments)
    #save(stats,file=file.path(dir,"track_stats.Rdata"))
    segments <- data.frame(ldply(segments),"condition"=condition_list[i])
   # tracks <- data.frame(ldply(tracks),"condition"=condition_list[i])
    save(segments,file = file.path(directory,condition_list[i],"segments.Rdata"))
   # save(tracks,file = file.path(directory,condition_list[i],"tracks.Rdata"))
    rm(segments)
   # rm(tracks)
    #track_stats_all[[basename(dir)]] <- data.frame(ldply(stats),"condition"=basename(dir))
    
    
  }
  #save data to the folder
  cat("Finishing off, putting all data together")
  for (i in condition_list){
    load(file=file.path(directory,i,"segments.Rdata"))
    #load(file = file.path(directory,i,"tracks.Rdata"))
    segments_all[[i]] <- segments
  #  msd_fit_all[[i]] <- tracks
  }
  
  
  save(segments_all,file=file.path(directory,paste("segments_all.Rdata")))
 # save(msd_fit_all,file=file.path(directory,paste("msd_fit_all.Rdata")))
#  save(track_stats_all,file=file.path(directory,paste("track_stats_all.Rdata")))
}

msd_analyze_data_mosaic_mask_parallel_intensityTM <- function(directory,condition_list,framerate,n,fitzero,pixelsize,fitMSD,offset,groundtruth=TRUE){
  segments_all <- list()
  msd_fit_all <- list()
  track_stats_all <- list()
  library(readr)
  for (i in 1:length(condition_list)){
    dirs <- file.path(directory,condition_list[i])
    filelist <- list.files(dirs,full.names = T,recursive = F,pattern = "^Traj_.*.\\mask2.csv$")
    filelist <- filelist[-grep(filelist,pattern = "^.*transformed.*$")]
    total <- seq(1:length(filelist))
    
    nodes <- detectCores()
    cl <- makeCluster(nodes-12)
    registerDoParallel(cl)
    dimensions <- dim
    offst <- offset
    pixelsize <- pixelsize
    n <- n
    fitMSD<- fitMSD
    fitzero <-fitzero
    framerate<-framerate
    groundtruth <- groundtruth
    cat(paste0("Analyzing ",condition_list[i],"\n"))
    output <-alply(filelist,.margins = 1,.parallel = T,.paropts = list(.export=c("filelist","dirs","groundtruth","n","framerate","offst","pixelsize","dimensions","fitMSD","fitzero"),.packages=c("MSDtracking","readr","stats","stringr")),function(j){
      Sys.sleep(0.1)
      tracks_simple <- read_csv(j)
      names(tracks_simple) <- c("track","X","Y","Z","time","frame","step_x","step_y","step_z","inMask","rawInt","meanInt","normInt", "sdInt","maxInt","minInt","distMask","invDistMask", "label", "center_x", "center_y", "center_z", "distanceCentroid3D", "distanceCentroid2D", "distanceToNearest3D", "distanceToNearest2D")
      
      tracks_simple <- tracks_simple[,c(6,2,3,1,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26)]
      tracks_simple$frame <- tracks_simple$frame
      segments <- data.frame(SEGMENT_STAT(tracks_simple),"cellID"=basename(j))
      segments <- ddply(segments,.variables = "track",function(x){
        if(any(x$inMask==TRUE)){
          x$trackInMask <-TRUE
        } else {
          x$trackInMask <-FALSE
        }
        return(x)
      })
      
      
      if(groundtruth){
        filelist_gt <- list.files(dirs,full.names = T,recursive = F,pattern = paste0("^Traj_transformed.*.",str_sub(basename(j),6)))
        #cat(paste0("number groundtruths found:", k))
        for(k in 1:length(filelist_gt)){
          tracks_simple_gt <- read_csv(filelist_gt[k])
          names(tracks_simple) <- c("track","X","Y","Z","time","frame","step_x","step_y","step_z","inMask","rawInt","meanInt","normInt", "sdInt","maxInt","minInt","distMask","invDistMask", "label", "center_x", "center_y", "center_z", "distanceCentroid3D", "distanceCentroid2D", "distanceToNearest3D", "distanceToNearest2D")
          segments[[paste0("inMask_gt",k)]] <- tracks_simple_gt$inMask
          segments[[paste0("rawInt_gt",k)]] <- tracks_simple_gt$rawInt
          segments[[paste0("meanInt_gt",k)]] <- tracks_simple_gt$meanInt
          segments[[paste0("sdInt_gt",k)]] <- tracks_simple_gt$sdInt
          segments[[paste0("maxInt_gt",k)]] <- tracks_simple_gt$maxInt
          segments[[paste0("minInt_gt",k)]] <- tracks_simple_gt$minInt
          segments[[paste0("distMask_gt",k)]] <- tracks_simple_gt$distMask
          segments[[paste0("invDistMask_gt",k)]] <- tracks_simple_gt$invDistMask
          segments[[paste0("center_x_gt",k)]] <- tracks_simple_gt$center_x
          segments[[paste0("center_y_gt",k)]] <- tracks_simple_gt$center_y
          segments[[paste0("center_z_gt",k)]] <- tracks_simple_gt$center_z
          segments[[paste0("distanceCentroid3D_gt",k)]] <- tracks_simple_gt$distanceCentroid3D
          segments[[paste0("distanceCentroid2D_gt",k)]] <- tracks_simple_gt$distanceCentroid2D
          segments[[paste0("distanceToNearest3D_gt",k)]] <- tracks_simple_gt$distanceToNearest3D
          segments[[paste0("distanceToNearest2D_gt",k)]] <- tracks_simple_gt$distanceToNearest2D
          segments <- ddply(segments,.variables = "track",function(x){
            if(any(x[[paste0("inMask_gt",k)]]==TRUE)){
              x[[paste0("trackInMask_gt",k)]] <-TRUE
            } else {
              x[[paste0("trackInMask_gt",k)]] <-FALSE
            }
            return(x)
          })
          
        }
        
      }
      track_msd <- TRACK_MSD(segments,n = n,framerate=framerate,pxsize = pixelsize,dim=dimensions)
      if(fitMSD==TRUE){
        #tracks <-  TRACK_MSD_fit(track_msd,n = n,fitzero = fitzero,framerate=framerate,pxsize = pixelsize,offset=offset,dim=dimensions)
        tracks <-  TRACK_MSD_fit(track_msd,n = n,fitzero = fitzero,framerate=framerate,pxsize = pixelsize,dim=dimensions)
        
        for (k in 1:length(tracks$track)){
          tracks$inMask = FALSE
          if(groundtruth){
            for(l in 1:length(filelist_gt)){
              tracks[[paste0("inMask_gt",l)]] = FALSE
            }
          }
        }
        
        for(l in 1:length(tracks$track)){
          tracks$inMask[l] <- any(segments$inMask[segments$track==tracks$track[l]]==1)
        }
        if(groundtruth){
          for(l in 1:length(filelist_gt)){
            
            for(m in 1:length(tracks$track)){
              tracks[[paste0("inMask_gt",l)]][m] <- any(segments[[paste0("inMask_gt",l)]][segments$track==tracks$track[m]]==1)
            }
          }
        }
      }
      if(fitMSD==TRUE){
        return(list("segments"=segments,"tracks"=tracks))
      }
      else {
        return(segments)
      }
    })
    stopCluster(cl)
    segments <- list()
    # tracks <- list()
    for(m in 1:length(output)){
      segments[[m]] <- output[[m]]
      #  tracks[[m]] <- output[[m]]$tracks
    }
    names(segments) <- basename(filelist)
    #names(tracks) <- basename(filelist)
    output <- NULL
    # stats <- TRACK_STAT(x=segments)
    #save(stats,file=file.path(dir,"track_stats.Rdata"))
    segments <- data.frame(ldply(segments),"condition"=condition_list[i])
    # tracks <- data.frame(ldply(tracks),"condition"=condition_list[i])
    save(segments,file = file.path(directory,condition_list[i],"segments.Rdata"))
    # save(tracks,file = file.path(directory,condition_list[i],"tracks.Rdata"))
    rm(segments)
    # rm(tracks)
    #track_stats_all[[basename(dir)]] <- data.frame(ldply(stats),"condition"=basename(dir))
    
    
  }
  #save data to the folder
  cat("Finishing off, putting all data together")
  for (i in condition_list){
    load(file=file.path(directory,i,"segments.Rdata"))
    #load(file = file.path(directory,i,"tracks.Rdata"))
    segments_all[[i]] <- segments
    #  msd_fit_all[[i]] <- tracks
  }
  
  
  save(segments_all,file=file.path(directory,paste("segments_all.Rdata")))
  # save(msd_fit_all,file=file.path(directory,paste("msd_fit_all.Rdata")))
  #  save(track_stats_all,file=file.path(directory,paste("track_stats_all.Rdata")))
}


sos_to_spoton_csv <- function(condition_list,directory,framerate,pixelsize){
  condition_list <- list.dirs(directory,full.names = F,recursive = F)
  #condition_list <- condition_list[c(4)]
  #condition_list <- "SNAP-SiR"
  segments_all <- list()
  msd_fit_all <- list()
  track_stats_all <- list()
  #load data
  for (i in 1:length(condition_list)){
    segments <- list()
    dir <- file.path(directory,condition_list[i])
    filelist <- list.dirs(dir,full.names = T,recursive = F)
    #filelist <- filelist[-grep("skip",x = filelist)]
    total <- length(filelist)
    # create progress bar
    for(j in 1:total){
      #  Sys.sleep(0.1)
      tracks_simple <- read.csv(file.path(filelist[j],"tracks.simple.filtered.txt"),sep = "\t",header = F)
      spoton <- tracks_simple[c(1,1,4,2,3)]
      spoton[,2] <- spoton[,2]/framerate/1000 #seconds
      spoton[,4:5] <- spoton[,4:5]*pixelsize/1000 #um
      spoton <- cbind(seq(0,nrow(spoton)-1),spoton)
      names(spoton) <- c("","frame","t","trajectory","x","y")
      write.csv(spoton,file = file.path(dirname(filelist[j]),paste(basename(dirname(filelist[j])),j,"tracks.spoton.txt",sep = "_")),row.names=FALSE,quote = FALSE)
    }
  }
}

