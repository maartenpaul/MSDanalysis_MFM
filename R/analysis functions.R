msd_analyze_data <- function(directory,condition_list,framerate,n,fitzero,min_length,pixelsize,fitMSD,offset,max_tracks,track_file_name="tracks.simple.filtered.txt",extension=""){
  segments_all <- list()
  msd_fit_all <- list()
  track_stats_all <- list()

  for (i in 1:length(condition_list)){
    segments <- list()
    dir <- file.path(directory,condition_list[i])
    filelist <- list.dirs(dir,full.names = T,recursive = F)
    total <- length(filelist)
    # create progress bar
    for(j in 1:total){
      Sys.sleep(0.1)
      tracks_simple <- read.csv(file.path(filelist[j],track_file_name),sep = "\t",header = F)
      names(tracks_simple) <- c("frame","X","Y","track","displacement","intensity","sigma","fit_error")
      segments[[j]] <- data.frame(SEGMENT_STAT(tracks_simple),"cellID"=basename(dirname(file.path(filelist[j],track_file_name))))
    }
    track_msd <- TRACK_MSD(segments,n = n,framerate=framerate,truncate = FALSE,pxsize = pixelsize)
    save(track_msd,file=file.path(dir,"track_msd.Rdata"))

    if(fitMSD){
      tracks <-  TRACK_MSD_fit(track_msd,n = n,fitzero = fitzero,framerate=framerate,offset=offset,pxsize = pixelsize)
      save(tracks,file=file.path(dir,"msd_fit.Rdata"))

      stats <- TRACK_STAT(x=segments)
      save(stats,file=file.path(dir,"track_stats.Rdata"))

    }
    segments_all[[basename(dir)]] <- data.frame(ldply(segments),"condition"=basename(dir))
    msd_fit_all[[basename(dir)]] <- data.frame(ldply(tracks),"condition"=basename(dir))
    track_stats_all[[basename(dir)]] <- data.frame(ldply(stats),"condition"=basename(dir))


  }
  #save data to the folder
  save(segments_all,file=file.path(directory,paste("segments_all_",extension,".Rdata")))
  save(msd_fit_all,file=file.path(directory,paste("msd_fit_all_",extension,".Rdata")))
  save(track_stats_all,file=file.path(directory,paste("track_stats_all_",extension,".Rdata")))
}

msd_analyze_data_mosaic <- function(directory,condition_list,framerate,n,fitzero,min_length,pixelsize,fitMSD,offset,max_tracks,track_file_name="tracks.simple.filtered.txt",extension="",dim){
  segments_all <- list()
  msd_fit_all <- list()
  track_stats_all <- list()

  for (i in 1:length(condition_list)){
    segments <- list()
    dir <- file.path(directory,condition_list[i])
    filelist <- list.files(dir,full.names = T,recursive = F,pattern = "^Traj_.*.\\.csv$")
    total <- length(filelist)
    # create progress bar
    for(j in 1:total){
      Sys.sleep(0.1)
      tracks_simple <- read.csv(filelist[i],sep = ",",header = T)
      names(tracks_simple) <- c("id","track","frame","X","Y","Z","m0","m1","m2","m3","m4","NPscore")
#trackID,pos_x,pos_y,pos_z,time,frame,step_x,step_y,step_z,inMask
      tracks_simple <- tracks_simple[,c(3,4,5,2,6,7,8,9,10,11,12)]
      tracks_simple$frame <- tracks_simple$frame-1
      segments[[j]] <- data.frame(SEGMENT_STAT(tracks_simple),"cellID"=basename(filelist[j]))
    }
    track_msd <- TRACK_MSD(segments,n = n,framerate=framerate,truncate = FALSE,pxsize = pixelsize,dim=dim)
    save(track_msd,file=file.path(dir,"track_msd.Rdata"))

    if(fitMSD){
      tracks <-  TRACK_MSD_fit(track_msd,n = n,fitzero = fitzero,framerate=framerate,pxsize = pixelsize,offset=offset,dim=dim)
      save(tracks,file=file.path(dir,"msd_fit.Rdata"))

      stats <- TRACK_STAT(x=segments)
      save(stats,file=file.path(dir,"track_stats.Rdata"))

    }
    segments_all[[basename(dir)]] <- data.frame(ldply(segments),"condition"=basename(dir))
    msd_fit_all[[basename(dir)]] <- data.frame(ldply(tracks),"condition"=basename(dir))
    track_stats_all[[basename(dir)]] <- data.frame(ldply(stats),"condition"=basename(dir))


  }
  #save data to the folder
  save(segments_all,file=file.path(directory,paste("segments_all_",extension,".Rdata")))
  save(msd_fit_all,file=file.path(directory,paste("msd_fit_all_",extension,".Rdata")))
  save(track_stats_all,file=file.path(directory,paste("track_stats_all_",extension,".Rdata")))
}

msd_analyze_data_mosaic_mask <- function(directory,condition_list,framerate,n,fitzero,min_length,pixelsize,fitMSD,offset,max_tracks,track_file_name="tracks.simple.filtered.txt",extension="",dim,groundtruth=FALSE){
  segments_all <- list()
  msd_fit_all <- list()
  track_stats_all <- list()

  for (i in 1:length(condition_list)){
    segments <- list()
    dir <- file.path(directory,condition_list[i])
    filelist <- list.files(dir,full.names = T,recursive = F,pattern = "^Traj_.*.\\mask.csv$")
    filelist <- filelist[-grep(filelist,pattern = "^.*transformed.*$")]
    total <- length(filelist)
    # create progress bar
    for(j in 1:total){
      Sys.sleep(0.1)
      library(readr)
      tracks_simple <- read_csv(filelist[j])
      names(tracks_simple) <- c("track","X","Y","Z","time","frame","step_x","step_y","step_z","inMask")
      #trackID,pos_x,pos_y,pos_z,time,frame,step_x,step_y,step_z,inMask
      tracks_simple <- tracks_simple[,c(6,2,3,1,4,5,7,8,9,10)]
                                        #3,4,5,2,6,7,8,9,10,11,12)]
      tracks_simple$frame <- tracks_simple$frame
      segments[[j]] <- data.frame(SEGMENT_STAT(tracks_simple),"cellID"=basename(filelist[j]))
      if(groundtruth){
        filelist_gt <- list.files(dir,full.names = T,recursive = F,pattern = paste0("^Traj_transformed.*.",basename(filelist[j])))
        cat(paste0("number groundtruths found:", k))
        for(k in 1:length(filelist_gt)){
          tracks_simple_gt <- read_csv(filelist_gt[k])
          names(tracks_simple_gt) <- c("track","X","Y","Z","time","frame","step_x","step_y","step_z","inMask")
          segments[[j]][[paste0("inMask_gt",k)]] <- tracks_simple_gt$inMask
        }

      }
    }
    track_msd <- TRACK_MSD(segments,n = n,framerate=framerate,truncate = FALSE,pxsize = pixelsize,dim=dim)
    save(track_msd,file=file.path(dir,"track_msd.Rdata"))

    if(fitMSD){
      tracks <-  TRACK_MSD_fit(track_msd,n = n,fitzero = fitzero,framerate=framerate,pxsize = pixelsize,offset=offset,dim=dim)
      for (k in 1:length(tracks)){
        tracks[[k]]$inMask = FALSE
        for(l in 1:length(filelist_gt)){
          tracks[[k]][[paste0("inMask_gt",l)]] = FALSE
        }


        for(l in 1:length(tracks[[k]]$track)){
          tracks[[k]]$inMask[l] <- any(segments[[k]]$inMask[segments[[k]]$track==tracks[[k]]$track[l]]==1)
        }
        if(groundtruth){
          for(l in 1:length(filelist_gt)){

            for(m in 1:length(tracks[[k]]$track)){
              tracks[[k]][[paste0("inMask_gt",l)]][m] <- any(segments[[k]][[paste0("inMask_gt",l)]][segments[[k]]$track==tracks[[k]]$track[m]]==1)
            }
        }
        }
        }
      }

      save(tracks,file=file.path(dir,"msd_fit.Rdata"))

     # stats <- TRACK_STAT(x=segments)
      #save(stats,file=file.path(dir,"track_stats.Rdata"))


    segments_all[[basename(dir)]] <- data.frame(ldply(segments),"condition"=basename(dir))
    msd_fit_all[[basename(dir)]] <- data.frame(ldply(tracks),"condition"=basename(dir))
    #track_stats_all[[basename(dir)]] <- data.frame(ldply(stats),"condition"=basename(dir))


  }
  #save data to the folder
  save(segments_all,file=file.path(directory,paste("segments_all_",extension,".Rdata")))
  save(msd_fit_all,file=file.path(directory,paste("msd_fit_all_",extension,".Rdata")))
  save(track_stats_all,file=file.path(directory,paste("track_stats_all_",extension,".Rdata")))
}

msd_analyze_data_mosaic_mask_parallel <- function(directory,condition_list,framerate,n,fitzero,min_length,pixelsize,fitMSD,offset,max_tracks,track_file_name="tracks.simple.filtered.txt",extension="",dim,groundtruth=TRUE){
  segments_all <- list()
  msd_fit_all <- list()
  track_stats_all <- list()
  library(readr)
  for (i in 1:length(condition_list)){
    dirs <- file.path(directory,condition_list[i])
    filelist <- list.files(dirs,full.names = T,recursive = F,pattern = "^Traj_.*.\\mask.csv$")
    filelist <- filelist[-grep(filelist,pattern = "^.*transformed.*$")]
    total <- seq(1:length(filelist))
    
    nodes <- detectCores()
    cl <- makeCluster(nodes-6)
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
      names(tracks_simple) <- c("track","X","Y","Z","time","frame","step_x","step_y","step_z","inMask")
      #trackID,pos_x,pos_y,pos_z,time,frame,step_x,step_y,step_z,inMask
      tracks_simple <- tracks_simple[,c(6,2,3,1,4,5,7,8,9,10)]
      #3,4,5,2,6,7,8,9,10,11,12)]
      tracks_simple$frame <- tracks_simple$frame
      segments <- data.frame(SEGMENT_STAT(tracks_simple),"cellID"=basename(j))
      if(groundtruth){
        filelist_gt <- list.files(dirs,full.names = T,recursive = F,pattern = paste0("^Traj_transformed.*.",basename(j)))
        #cat(paste0("number groundtruths found:", k))
        for(k in 1:length(filelist_gt)){
          tracks_simple_gt <- read_csv(filelist_gt[k])
          names(tracks_simple_gt) <- c("track","X","Y","Z","time","frame","step_x","step_y","step_z","inMask")
          segments[[paste0("inMask_gt",k)]] <- tracks_simple_gt$inMask
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
    tracks <- list()
    for(m in 1:length(output)){
      segments[[m]] <- output[[m]]$segments
      tracks[[m]] <- output[[m]]$tracks
    }
    names(segments) <- basename(filelist)
    names(tracks) <- basename(filelist)
    output <- NULL
    # stats <- TRACK_STAT(x=segments)
    #save(stats,file=file.path(dir,"track_stats.Rdata"))
    
    
    segments_all[[basename(dirs)]] <- data.frame(ldply(segments),"condition"=basename(dirs))
    msd_fit_all[[basename(dirs)]] <- data.frame(ldply(tracks),"condition"=basename(dirs))
    #track_stats_all[[basename(dir)]] <- data.frame(ldply(stats),"condition"=basename(dir))
    
    
  }
  #save data to the folder
  save(segments_all,file=file.path(directory,paste("segments_all.Rdata")))
  save(msd_fit_all,file=file.path(directory,paste("msd_fit_all.Rdata")))
  save(track_stats_all,file=file.path(directory,paste("track_stats_all.Rdata")))
}

msd_analyze_data_mosaic_mask_parallel_intensity <- function(directory,condition_list,framerate,n,fitzero,min_length,pixelsize,fitMSD,offset,max_tracks,track_file_name="tracks.simple.filtered.txt",extension="",dim,groundtruth=TRUE){
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
      if(length(names(tracks_simple))==16){
        names(tracks_simple) <- c("track","X","Y","Z","time","frame","step_x","step_y","step_z","inMask","rawInt","meanInt","sdInt","maxInt","minInt","distMask")
        tracks_simple <- tracks_simple[,c(6,2,3,1,4,5,7,8,9,10,11,12,13,14,15,16)]
        
      } else if (length(names(tracks_simple))==17) {
        names(tracks_simple) <- c("track","X","Y","Z","time","frame","step_x","step_y","step_z","inMask","rawInt","normInt","meanInt","sdInt","maxInt","minInt","distMask")
        tracks_simple <- tracks_simple[,c(6,2,3,1,4,5,7,8,9,10,11,12,13,14,15,16,17)]
        
        
      }
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
          if(length(names(tracks_simple))==16){
            names(tracks_simple) <- c("track","X","Y","Z","time","frame","step_x","step_y","step_z","inMask","rawInt","meanInt","sdInt","maxInt","minInt","distMask")
          } else if (length(names(tracks_simple))==17) {
            names(tracks_simple) <- c("track","X","Y","Z","time","frame","step_x","step_y","step_z","inMask","rawInt","normInt","meanInt","sdInt","maxInt","minInt","distMask")

          }          
          segments[[paste0("inMask_gt",k)]] <- tracks_simple_gt$inMask
          segments[[paste0("rawInt_gt",k)]] <- tracks_simple_gt$rawInt
          segments[[paste0("meanInt_gt",k)]] <- tracks_simple_gt$meanInt
          segments[[paste0("sdInt_gt",k)]] <- tracks_simple_gt$sdInt
          segments[[paste0("maxInt_gt",k)]] <- tracks_simple_gt$maxInt
          segments[[paste0("minInt_gt",k)]] <- tracks_simple_gt$minInt
          segments[[paste0("distMask_gt",k)]] <- tracks_simple_gt$distMask
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
    tracks <- list()
    for(m in 1:length(output)){
      segments[[m]] <- output[[m]]$segments
      tracks[[m]] <- output[[m]]$tracks
    }
    names(segments) <- basename(filelist)
    names(tracks) <- basename(filelist)
    output <- NULL
    # stats <- TRACK_STAT(x=segments)
    #save(stats,file=file.path(dir,"track_stats.Rdata"))
    segments <- data.frame(ldply(segments),"condition"=condition_list[i])
    tracks <- data.frame(ldply(tracks),"condition"=condition_list[i])
    save(segments,file = file.path(directory,condition_list[i],"segments.Rdata"))
    save(tracks,file = file.path(directory,condition_list[i],"tracks.Rdata"))
    rm(segments)
    rm(tracks)
    #track_stats_all[[basename(dir)]] <- data.frame(ldply(stats),"condition"=basename(dir))
    
    
  }
  #save data to the folder
  cat("Finishing off, putting all data together")
  for (i in condition_list){
    load(file=file.path(directory,i,"segments.Rdata"))
    load(file = file.path(directory,i,"tracks.Rdata"))
    segments_all[[i]] <- segments
    msd_fit_all[[i]] <- tracks
  }
  
  
  save(segments_all,file=file.path(directory,paste("segments_all.Rdata")))
  save(msd_fit_all,file=file.path(directory,paste("msd_fit_all.Rdata")))
  save(track_stats_all,file=file.path(directory,paste("track_stats_all.Rdata")))
}


msd_analyze_data_submask <- function(directory,condition_list,framerate,n,fitzero,min_length,pixelsize,fitMSD,offset,max_tracks,track_file_name="tracks.simple.filtered.txt",extension=""){
  library(RImageJROI)
  library(spatstat.utils)
  library(spatstat)
  segments_inside_all <- list()
  segments_outside_all <- list()

  msd_fit_inside_all <- list()
  msd_fit_outside_all <- list()
  track_stats_inside_all <- list()
  track_stats_outside_all <- list()


  for (i in 1:length(condition_list)){
    segments_inside <- list()
    segments_outside <- list()
    dir <- file.path(directory,condition_list[i])
    filelist <- list.dirs(dir,full.names = T,recursive = F)
    total <- length(filelist)
    # create progress bar
    for(j in 1:total){
      Sys.sleep(0.1)
      tracks_simple <- read.csv(file.path(filelist[j],track_file_name),sep = "\t",header = F)
      tracks_simple$inside <- FALSE
      roi <- read.ijzip(file.path(filelist[j],"mask_a.zip"))
      roi <- llply(roi, function(x){
        area <- Area.xypolygon( list(x = x$coords[, 1], y = x$coords[, 2]))
        if (area<0){
          x$coords <- x$coords[ nrow(x$coords):1, ]
        }
        return(x)
      })
      roi <- llply(roi, function(x){
        ij2spatstat(x)
      })




      for(i in 1:length(roi)){
        tracks_simple$inside<-tracks_simple$inside|inside.owin(x = tracks_simple$V2,y = tracks_simple$V3,w=roi[[i]])
      }
      tracks_simple <- ddply(tracks_simple,"V4",function(x){
        if(any(x$inside)){
          x$inside <- TRUE
        }
        return(x)
      })
      inside <-tracks_simple[tracks_simple$inside==TRUE,1:8]
      track_id <- unique(inside$V4)
      for(k in 1:length(track_id)){
        inside$V4[inside$V4==track_id[k]] <- (k-1)
      }
        outside <- tracks_simple[tracks_simple$inside==FALSE,1:8]
        track_id <- unique(outside$V4)

        for(k in 1:length(track_id)){
          outside$V4[outside$V4==track_id[k]] <- (k-1)
        }
      write.table(x = cbind(sprintf("%.2f",inside$V1),sprintf("%.2f",inside$V2),sprintf("%.2f",inside$V3),inside$V4,
                            sprintf("%.2f",inside$V5),sprintf("%.2f",inside$V6),sprintf("%.2f",inside$V7),sprintf("%.2f",inside$V8)),
                            file = file.path(filelist[j],"tracks.simple.filtered_a.txt"),sep = "\t",col.names = F, row.names = F,quote = F)
      write.table(x = cbind(sprintf("%.2f",outside$V1),sprintf("%.2f",outside$V2),sprintf("%.2f",outside$V3),outside$V4,
                            sprintf("%.2f",outside$V5),sprintf("%.2f",outside$V6),sprintf("%.2f",outside$V7),sprintf("%.2f",outside$V8)),
                  file = file.path(filelist[j],"tracks.simple.filtered_b.txt"),sep = "\t",col.names = F, row.names = F,quote = F)
    names(tracks_simple) <- c("frame","X","Y","track","displacement","intensity","sigma","fit_error","inside")
      segments_inside[[j]] <- data.frame(SEGMENT_STAT(tracks_simple[tracks_simple$inside==TRUE,]),"cellID"=basename(dirname(file.path(filelist[j],track_file_name))))
      segments_outside[[j]] <- data.frame(SEGMENT_STAT(tracks_simple[tracks_simple$inside==FALSE,]),"cellID"=basename(dirname(file.path(filelist[j],track_file_name))))

    }
    track_msd_inside <- TRACK_MSD(llply(segments_inside,function(x) {
      x[x$inside,]
    }),n = n,framerate=framerate,truncate = FALSE)
    track_msd_outside <- TRACK_MSD(llply(segments_outside,function(x) {
      x[!x$inside,]
    }),n = n,framerate=framerate,truncate = FALSE)


      tracks_inside <-  TRACK_MSD_fit(track_msd_inside,n = n,fitzero = fitzero,framerate=framerate,offset)
      tracks_outside <-  TRACK_MSD_fit(track_msd_outside,n = n,fitzero = fitzero,framerate=framerate,offset)
      save(tracks_inside,file=file.path(dir,"msd_inside_fit.Rdata"))
      save(tracks_outside,file=file.path(dir,"msd_outside_fit.Rdata"))

      stats_inside <- TRACK_STAT(x=segments_inside)
      stats_outside <- TRACK_STAT(x=segments_outside)



    segments_inside_all[[basename(dir)]] <- data.frame(ldply(segments_inside),"condition"=basename(dir))
    segments_outside_all[[basename(dir)]] <- data.frame(ldply(segments_outside),"condition"=basename(dir))

    msd_fit_inside_all[[basename(dir)]] <- data.frame(ldply(tracks_inside),"condition"=basename(dir))
    msd_fit_outside_all[[basename(dir)]] <- data.frame(ldply(tracks_outside),"condition"=basename(dir))
    track_stats_inside_all[[basename(dir)]] <- data.frame(ldply(stats_inside),"condition"=basename(dir))
    track_stats_outside_all[[basename(dir)]] <- data.frame(ldply(stats_outside),"condition"=basename(dir))



  }
  #save data to the folder
  save(segments_inside_all,file=file.path(directory,paste0("segments_inside_all_",extension,".Rdata")))
  save(segments_outside_all,file=file.path(directory,paste0("segments_outside_all_",extension,".Rdata")))

  save(msd_fit_inside_all,file=file.path(directory,paste0("msd_fit_inside_all_",extension,".Rdata")))
  save(msd_fit_outside_all,file=file.path(directory,paste0("msd_fit_outside_all_",extension,".Rdata")))

  save(track_stats_all,file=file.path(directory,paste0("track_stats_all_",extension,".Rdata")))
  save(track_stats_inside_all,file=file.path(directory,paste0("track_stats_inside_all_",extension,".Rdata")))

}

msd_analyze_data_submask <- function(directory,condition_list,framerate,n,fitzero,min_length,pixelsize,fitMSD,offset,max_tracks,track_file_name="tracks.simple.filtered.txt",extension=""){
  library(RImageJROI)
  library(spatstat.utils)
  library(spatstat)
  segments_inside_all <- list()
  segments_outside_all <- list()
  
  msd_fit_inside_all <- list()
  msd_fit_outside_all <- list()
  track_stats_inside_all <- list()
  track_stats_outside_all <- list()
  
  
  for (i in 1:length(condition_list)){
    segments_inside <- list()
    segments_outside <- list()
    dir <- file.path(directory,condition_list[i])
    filelist <- list.dirs(dir,full.names = T,recursive = F)
    # create progress bar
    
    nodes <- detectCores()
    cl <- makeCluster(nodes)
    registerDoParallel(cl)
    for(j in 1:total){
      Sys.sleep(0.1)
      tracks_simple <- read.csv(file.path(filelist[j],track_file_name),sep = "\t",header = F)
      tracks_simple$inside <- FALSE
      roi <- read.ijzip(file.path(filelist[j],"mask_a.zip"))
      roi <- llply(roi, function(x){
        area <- Area.xypolygon( list(x = x$coords[, 1], y = x$coords[, 2]))
        if (area<0){
          x$coords <- x$coords[ nrow(x$coords):1, ]
        }
        return(x)
      })
      roi <- llply(roi, function(x){
        ij2spatstat(x)
      })
      
      
      for(i in 1:length(roi)){
        tracks_simple$inside<-tracks_simple$inside|inside.owin(x = tracks_simple$V2,y = tracks_simple$V3,w=roi[[i]])
      }
      tracks_simple <- ddply(tracks_simple,"V4",function(x){
        if(any(x$inside)){
          x$inside <- TRUE
        }
        return(x)
      })
      inside <-tracks_simple[tracks_simple$inside==TRUE,1:8]
      track_id <- unique(inside$V4)
      for(k in 1:length(track_id)){
        inside$V4[inside$V4==track_id[k]] <- (k-1)
      }
      outside <- tracks_simple[tracks_simple$inside==FALSE,1:8]
      track_id <- unique(outside$V4)
      
      for(k in 1:length(track_id)){
        outside$V4[outside$V4==track_id[k]] <- (k-1)
      }
      write.table(x = cbind(sprintf("%.2f",inside$V1),sprintf("%.2f",inside$V2),sprintf("%.2f",inside$V3),inside$V4,
                            sprintf("%.2f",inside$V5),sprintf("%.2f",inside$V6),sprintf("%.2f",inside$V7),sprintf("%.2f",inside$V8)),
                  file = file.path(filelist[j],"tracks.simple.filtered_a.txt"),sep = "\t",col.names = F, row.names = F,quote = F)
      write.table(x = cbind(sprintf("%.2f",outside$V1),sprintf("%.2f",outside$V2),sprintf("%.2f",outside$V3),outside$V4,
                            sprintf("%.2f",outside$V5),sprintf("%.2f",outside$V6),sprintf("%.2f",outside$V7),sprintf("%.2f",outside$V8)),
                  file = file.path(filelist[j],"tracks.simple.filtered_b.txt"),sep = "\t",col.names = F, row.names = F,quote = F)
      names(tracks_simple) <- c("frame","X","Y","track","displacement","intensity","sigma","fit_error","inside")
      segments_inside[[j]] <- data.frame(SEGMENT_STAT(tracks_simple[tracks_simple$inside==TRUE,]),"cellID"=basename(dirname(file.path(filelist[j],track_file_name))))
      segments_outside[[j]] <- data.frame(SEGMENT_STAT(tracks_simple[tracks_simple$inside==FALSE,]),"cellID"=basename(dirname(file.path(filelist[j],track_file_name))))
      
    }
    track_msd_inside <- TRACK_MSD(llply(segments_inside,function(x) {
      x[x$inside,]
    }),n = n,framerate=framerate,truncate = FALSE)
    track_msd_outside <- TRACK_MSD(llply(segments_outside,function(x) {
      x[!x$inside,]
    }),n = n,framerate=framerate,truncate = FALSE)
    
    
    tracks_inside <-  TRACK_MSD_fit(track_msd_inside,n = n,fitzero = fitzero,framerate=framerate,offset)
    tracks_outside <-  TRACK_MSD_fit(track_msd_outside,n = n,fitzero = fitzero,framerate=framerate,offset)
    save(tracks_inside,file=file.path(dir,"msd_inside_fit.Rdata"))
    save(tracks_outside,file=file.path(dir,"msd_outside_fit.Rdata"))
    
    stats_inside <- TRACK_STAT(x=segments_inside)
    stats_outside <- TRACK_STAT(x=segments_outside)
    
    
    
    segments_inside_all[[basename(dir)]] <- data.frame(ldply(segments_inside),"condition"=basename(dir))
    segments_outside_all[[basename(dir)]] <- data.frame(ldply(segments_outside),"condition"=basename(dir))
    
    msd_fit_inside_all[[basename(dir)]] <- data.frame(ldply(tracks_inside),"condition"=basename(dir))
    msd_fit_outside_all[[basename(dir)]] <- data.frame(ldply(tracks_outside),"condition"=basename(dir))
    track_stats_inside_all[[basename(dir)]] <- data.frame(ldply(stats_inside),"condition"=basename(dir))
    track_stats_outside_all[[basename(dir)]] <- data.frame(ldply(stats_outside),"condition"=basename(dir))
    
    
    
  }
  #save data to the folder
  save(segments_inside_all,file=file.path(directory,paste0("segments_inside_all_",extension,".Rdata")))
  save(segments_outside_all,file=file.path(directory,paste0("segments_outside_all_",extension,".Rdata")))
  
  save(msd_fit_inside_all,file=file.path(directory,paste0("msd_fit_inside_all_",extension,".Rdata")))
  save(msd_fit_outside_all,file=file.path(directory,paste0("msd_fit_outside_all_",extension,".Rdata")))
  
  save(track_stats_all,file=file.path(directory,paste0("track_stats_all_",extension,".Rdata")))
  save(track_stats_inside_all,file=file.path(directory,paste0("track_stats_inside_all_",extension,".Rdata")))
  
}


msd_histogram <- function(msd_fit_all,directory,name="",threshold=0.05,order=NULL,ymax=0.12,merge_var="cellID"){
 library(plyr)
  library(ggplot2)
  library(reshape2)
  library(ggpol)
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
    discrete_scale("fill","Publication",manual_pal(values = c("#c00000","#fdae61","#1f497d","#6599d9","#542788","#de77ae","#271d68","#6dc5aa")), ...)

  }

  scale_colour_Publication <- function(...){
    library(scales)
    discrete_scale("colour","Publication",manual_pal(values = c("#c00000","#fdae61","#1f497d","#6599d9","#542788","#de77ae","#217d68","#6dc5aa")), ...)

  }

  tracklength <- ldply(msd_fit_all,function(x){
    nrow(x)
  })
  tracks <- ldply(msd_fit_all)
  tracks$condition <-factor(tracks$condition, levels=names(msd_fit_all))
  out <- ddply(tracks,.variables = c("condition",merge_var)  ,function(x){
    out <- hist(log10(x$D),breaks = seq(-7,2,0.1),plot = F)
    return(out$counts/sum(out$counts))
  })

  n_tracks <- ddply(tracks,.variables = "condition",function(x){
    nrow(x)
  })

  n_exp <- ddply(tracks,.variables = "condition",function(x){
    length(unique(x$experiment))
  })

  n_cells <- ddply(tracks,.variables = "condition",function(x){
    length(unique(x$cellID))
  })

  n_statistics <- data.frame(n_tracks$condition,"n_tracks" = n_tracks$V1,"n_cells" = n_cells$V1,"n_exp"=n_exp$V1)
  write.table(n_statistics,file=file.path(directory,paste0(name,"_N_statistics.txt")),row.names = F,col.names = T)

  means <- ddply(out,.variables = "condition",function(x){
    colMeans(x[,-c(1,2)])
  })
  sds <- ddply(out,.variables = "condition",function(x){
    apply(x[,-c(1,2)], 2, sd)#/sqrt(length(table(x$cellID))) #if se instead of stdev
  })

  means <- melt(means)
  sds <- melt(sds)
  mids <- seq(-6.9,2,0.1)
  mids <- 10^mids

  histdata <- data.frame(".id"=means$condition,"x"=rep(mids,each=length(msd_fit_all)),'mean'=means$value,'se'=sds$value)
  if(!is.null(order)){
    histdata$.id <- factor(histdata$.id,levels(histdata$.id)[order])
  }

  #histdata$se <- 0
  histdata <- na.omit(histdata)

  data <- data.frame(histdata[,1:2])
  names(data) <- c("x","y")

  #levels(histdata$.id) <- c(levels(histdata$.id)[2],levels(histdata$.id)[1])

  if(length(msd_fit_all)>2){
  q1 <- ggplot(histdata, aes(x=x, y=mean,colour=.id,fill=.id)) +
    geom_line(stat="identity")+
    geom_ribbon(aes(ymin=mean-abs(se), ymax=mean+abs(se),colour=NULL), alpha=0.4)+
    ylim(-0.005,ymax)+
    # geom_smooth(stat="smooth",span = 0.2)+
    theme(text = element_text(size=12),legend.position = "none")+
    facet_wrap(~.id,nrow=2,dir="v") +scale_x_continuous(trans="log10")+
    scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=12)+
    theme(legend.position = "none")+ geom_vline(aes(xintercept=threshold),color = "red",linetype="dashed")+
    xlab(expression(D[app]~mu~m^{2}/s))+ylab("")+theme(axis.line.x = element_line(color="black"),
                                                       axis.line.y = element_line(color="black"))
  } else {
    q1 <- ggplot(histdata, aes(x=x, y=mean,colour=.id,fill=.id)) +
      geom_line(stat="identity")+
      geom_ribbon(aes(ymin=mean-se, ymax=mean+se,colour=NULL), alpha=0.4)+

      # geom_smooth(stat="smooth",span = 0.2)+
      ylim(-0.005,ymax)+
      theme(text = element_text(size=12),legend.position = "none")+
      facet_wrap(~.id,nrow=1,dir="h") +scale_x_continuous(trans="log10")+
      scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=12)+
      theme(legend.position = "none")+ geom_vline(aes(xintercept=threshold),color = "red",linetype="dashed")+
      xlab(expression(D[app]~mu~m^{2}/s))+ylab("")+theme(axis.line.x = element_line(color="black"),
                                                               axis.line.y = element_line(color="black"))
  }

  print(q1)

  ggsave(filename = file.path(directory,paste0(name,"_D_histogram.pdf")),plot = q1,dpi=300,units="mm",height=100, width=150)
  ggsave(filename = file.path(directory,paste0(name,"_D_histogram.png")),plot = q1+
  theme(
    plot.background = element_rect(fill = "transparent"),axis.line.x = element_line(colour="grey"),
    axis.line.y = element_line(colour="grey"),
    axis.ticks = element_line(colour="grey"),
    axis.title.y = element_text(angle=90,vjust =2,colour="grey"),
    axis.title.x = element_text(vjust = -0.2,colour="grey"),
    axis.text = element_text(colour="grey")),bg = "transparent",dpi=400,limitsize = F,units="mm",height=100, width=150)

  results <- llply(msd_fit_all,function(x) {
    ddply(x,.variables = merge_var,function(x){
      out <- table(x$D>threshold)/length(x$D)
      if(length(out)==1){
        if(names(out)=="FALSE"){
          out <- c(1,0)
        } else {
          out <- c(0,1)
        }

      }
      # out <- c(out,x$cellID[1])
      return(out)
    })
  }
  )

  results2 <- ldply(results)
  results2$.id <- factor(results2$.id,names(msd_fit_all))
  names(results2) <- c(".id","condition","immobile","mobile")
  write.table(results2,file=file.path(directory,paste0(name,"_imm_fractions.txt")),row.names = F,col.names = T)
  write.table(tracklength,file=file.path(directory,paste0(name,"_N_tracks.txt")),row.names = F,col.names = T)
  if(!is.null(order)){
    results2$.id <- factor(results2$.id,levels(histdata$.id)[order])
  }
  results2 <- na.omit(results2)
  out <- ldply(results, function(x){
    means <- mean(x[,2])
    sems <- sd(x[,2])/sqrt(length(x[,2]))
    return(data.frame("mean"=means,"se"=sems))
  }
  )
  out$.id <- factor(out$.id,names(msd_fit_all))

  #out$.id <- factor(out$.id,out$.id[c(1,2)])


  #out$.id <- factor(out$.id,levels(histdata$.id)[c(1,2,3,4)])
  out <- na.omit(out)
  q2 <- ggplot(out, aes(x=.id, y=mean,fill=.id)) +
    geom_bar(position=position_dodge(width = 1.1), stat="identity")+  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.2)+xlab("")+ylab("immobile fraction")+
    scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=16)+theme(legend.position = "none")+
    ylim(0,1.1)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+ylab("immobile fraction")
  print(q2)

  ggsave(filename = file.path(directory,paste0(name,"_immobile_bar_graph.pdf")),plot = q2,units="mm",height=100, width=150)

  ggsave(filename = file.path(directory,paste0(name,"_immobile_bar_graph.png")),plot = q2+
           theme(
             plot.background = element_rect(fill = "transparent"),axis.line.x = element_line(colour="grey"),
             axis.line.y = element_line(colour="grey"),
             axis.ticks = element_line(colour="grey"),
             axis.title.y = element_text(angle=90,vjust =2,colour="grey"),
             axis.title.x = element_text(vjust = -0.2,colour="grey"),
             axis.text = element_text(colour="grey")),bg = "transparent",dpi=400,limitsize = F,units="mm",height=100, width=150)


  if(!is.null(order)&&length(msd_fit_all)%%2==0){
    datapoints <- list()
    for (i in 1:length(msd_fit_all)){
      if (i%%2==1){
        one <- data.frame(results[[i]],'date'=ldply(strsplit(results[[i]]$cellID,split = "_"))[,1])
        one_stat <- aggregate(one$'FALSE.'~one$date, data=one, FUN=function(x) c(mean=mean(x)))
        two <- data.frame(results[[i+1]],'date'=ldply(strsplit(results[[i+1]]$cellID,split = "_"))[,1])
        two_stat <- aggregate(two$'FALSE.'~two$date, data=two, FUN=function(x) c(mean=mean(x)))
        merged <- (two_stat[,2]-one_stat[,2])/one_stat[,2]*100
        datapoints[[ldply(strsplit(names(results),split = " "))[,1][i]]] <- merged
      }
    }
    datapoints <- melt(datapoints)
    datapoints <- aggregate(value~L1,data=datapoints, FUN=function(x) c(mean=mean(x)),simplify=T)
    qx <- ggplot(datapoints,aes(x=L1,y=value,fill=L1))+geom_bar(stat = "identity")+
      scale_colour_Publication()+scale_fill_Publication()+theme(legend.position = "none")+
      xlab("")+ylab("% increase in immobilization")+ylim(c(0,100))
    print(qx)
    ggsave(filename = file.path(directory,paste0(name,"_immobile_increase.pdf")),plot = qx)


  }

#dotplot
  q3 <- ggplot(results2,aes(x=factor(.id),y=immobile,fill=.id))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    geom_dotplot(binaxis = "y", stackdir = "center",dotsize = 0.7)+ xlab("")+
    scale_colour_Publication()+scale_fill_Publication()+theme(legend.position = "none")+theme_Publication(base_size=12)+
    ylim(0,1.1)+stat_summary(fun.data=mean_cl_normal,
                             geom="errorbar", color="black", width=0.3,size=1) +
    stat_summary(fun.y=mean, geom="point", color="black")
  print(q3)
  ggsave(filename = file.path(directory,paste0(name,"_immobile_dot_plot.pdf")),plot = q3,units="mm",height=100, width=150)
  ggsave(filename = file.path(directory,paste0(name,"_immobile_dot_plot.png")),plot = q3+
           theme(
             plot.background = element_rect(fill = "transparent"),axis.line.x = element_line(colour="grey"),
             axis.line.y = element_line(colour="grey"),
             axis.ticks = element_line(colour="grey"),
             axis.title.y = element_text(angle=90,vjust =2,colour="grey"),
             axis.title.x = element_text(vjust = -0.2,colour="grey"),
             axis.text = element_text(colour="grey")),bg = "transparent",dpi=400,limitsize = F,units="mm",height=100, width=150)


#boxplot
  q4 <- ggplot(results2,aes(x=factor(.id),y=immobile,fill=.id))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    geom_boxjitter(errorbar.draw = TRUE,jitter.height = 0, jitter.width = 0.04)+ xlab("")+
    scale_colour_Publication()+scale_fill_Publication()+theme(legend.position = "none")+theme_Publication(base_size=12)+
    ylim(0,1.1)
  print(q4)
  ggsave(filename = file.path(directory,paste0(name,"_immobile_boxplot.pdf")),plot = q4,units="mm",height=100, width=150)
  ggsave(filename = file.path(directory,paste0(name,"_immobile_boxplot.png")),plot = q4,bg = "transparent",dpi=400,limitsize = F,units="mm",height=100, width=75)

  #qa <- grid.arrange(q1,q2,ncol=2)
}

merge_dataset <- function(dirs,conds,con_name){
  if(length(dirs)!=length(conds)){
    stop("dirs not equal to conds")
  }
  for(i in 1:length(dirs)){
    load(file.path(dirs[i],"msd_fit_all.Rdata"))
    date <- basename(dirs[i])
    date <- strsplit(date," ")
    date <- date[[1]][1]
    msd_fit_all[[conds[i]]]$cellID <- paste0(date,"_",msd_fit_all[[conds[i]]]$cellID)
    msd_fit_all[[conds[i]]]$experiment <- date
    if (i==1) {
      all <- msd_fit_all[[conds[i]]]
    } else {
      all <- rbind(all,msd_fit_all[[conds[i]]])
    }
  }
  all$condition <- con_name
  return(all)
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
