track_msd <- function(x,n=5,framerate=0.1,pxsize=100,truncate=F,dim=2){
  if (dim==2){
    x$X <- (x$X*pxsize)/1000 #convert pixelsize to um
    x$Y <- (x$Y*pxsize)/1000
  } else if (dim==3){
    x$X <- (x$X*pxsize[1])/1000 #convert pixelsize to um
    x$Y <- (x$Y*pxsize[2])/1000
    x$Z <- (x$Z*pxsize[3])/1000
  }
  #calculate MSD for n time intervals
  coef <- ddply(x,.variables = "track",.fun= function(x) {
    result <- vector(length = n)
    #put first frame of track at zero
    x$frame  <- x$frame-x$frame[1]+1
    if(nrow(x)>n+1){
      if(truncate){ #test effect of truncating tracks to first number of steps, skipped by default
        x <- x[1:6,]
      }
      #loop over different time intervals
      sapply(1:n,function(j){
        #loop over all steps of tracks
        sum <- unlist(sapply(1:nrow(x),simplify = F,function(k){
          if (is.element(x[k,1]+j,x[,1])){ #check if point is present in track, this deals with gaps
            which_point <- which(x[,1]==x[k,1]+j)
            if (dim==2){
              sum <- ((x[k,2]-x[which_point,2])^2+(x[k,3]-x[which_point,3])^2) #determine squared displacement between points
            } else if (dim==3){
              sum <- ((x[k,2]-x[which_point,2])^2+(x[k,3]-x[which_point,3])^2+(x[k,5]-x[which_point,5])^2)
            }
            return(sum)

          }}))

        result[j] <<- mean(sum) #calculate MSD for j time interval

      })
      return(result) #return MSD for different time intervals as vector

    }
  })
  return(coef)
}


TRACK_MSD <- function(x,n=5,framerate=0.1,pxsize=100,truncate=F,dim=2){
  UseMethod("TRACK_MSD")
}

TRACK_MSD.default <- function(x,n=5,framerate=0.1,pxsize=100,truncate=F,dim=2){
  stop("MSD requires data frame")
}

TRACK_MSD.data.frame <-  function(x,n=5,framerate=0.1,pxsize=100,truncate=F,dim=2){
  track_msd(x,n,framerate,pxsize,truncate,dim)

}

TRACK_MSD.list <-  function(x,n=5,framerate=0.1,pxsize=100,truncate=F,dim=2){

  llply(x,function(x){
    out <- TRACK_MSD(x,n,framerate,pxsize,truncate,dim)
    out$cellID <- x$cellID[1]
    return(out)
  })

}

