

track_stat <- function(x,framerate=30,pxsize=100){

  x$X <- (x$X*pxsize)/1000
  x$Y <- (x$Y*pxsize)/1000
  out <- ddply(x,.variables = "track",.fun= function(x) {
    speed <- 0
    for (i in 2:nrow(x)){
      speed <- speed + (x$X[i]-x$X[i-1])^2+(x$Y[i]-x$Y[i-1])^2/((x$frame[i]-x$frame[i-1])*framerate)
    }
    speed <- speed/nrow(x)

    coord <- cbind(x$X,x$Y)
    # use pricipal component analysis on X and Y coordinates to get eigenvectors: major and minor axis
    D <- princomp(coord)
    angle <- atan2(D$loadings[2,1],D$loadings[1,1])
    #calculate convex hull and futher statistics
    y <- chull(coord)
    area <- pracma::polyarea(coord[rev(y),1], coord[rev(y),2])
    perimeter <- pracma::poly_length(coord[rev(y),1], coord[rev(y),2])
    D_chull <- princomp(coord[y,])


    #return(data.frame("sd"=((sd(x$X)+sd(x$Y))/2)*2.35,"N"=nrow(x),"channel"=1))
    # return(data.frame(,"N"=nrow(x),"channel"=1))
    return(data.frame("N"=nrow(x),"meanX"=mean(x$X),"meanY"=mean(x$Y),"meanspeed"=speed ,
                      "sd"=((sd(x$X)+sd(x$Y))/2),"sdpri"=((D$sdev[1]+D$sdev[2])/2),"major"=D$sdev[1],"minor"=D$sdev[2],
                      "width"=(max(D$scores[,1])-min(D$scores[,1])),"ratio"=(D$sdev[1]/D$sdev[2]),"angle"=angle,
                      "chull_area"=area,"chull_perimeter"=perimeter,"chull_major"=D_chull$sdev[1],"chull_minor"=D_chull$sdev[2]))
  })

  return(out)
}

TRACK_STAT <- function(x,framerate=30,pxsize=100){
  UseMethod("TRACK_STAT")
}

TRACK_STAT.default <- function(x,framerate=30,pxsize=100){
  stop("TRACK_STAT requires data frame")
}

TRACK_STAT.data.frame <-  function(x,framerate=30,pxsize=100){
  track_stat(x,framerate,pxsize)

}

TRACK_STAT.list <-  function(x,framerate=30,pxsize=100){
  llply(x,function(x){
    TRACK_STAT(x,framerate,pxsize)
  })
}


segment_stat <- function(x){

  get_angle <- function(x){
    seg_angle <- vector()
    for(i in 1:(nrow(x)-2)){
      A <- as.numeric(x[i,2:3])
      B <- as.numeric(x[i+1,2:3])
      C <- as.numeric(x[i+2,2:3])
      AB <- B-A
      CB <- C-B

      #dAB <- sqrt((B[1]-A[1])^2+(B[2]-A[2])^2)
      #dBC <- sqrt((C[1]-B[1])^2+(C[2]-B[2])^2)
      #Formula obtained from https://gitlab.com/anders.sejr.hansen/anisotropy
      angle <- abs(atan2(det(cbind(AB,CB)),AB%*%CB))
      angle <- angle/pi*180
      seg_angle <- c(seg_angle,angle)


    }
    seg_angle <- c(-1,seg_angle,-1)
    return(seg_angle)
  }
  get_angle_3D <- function(x){
    seg_angle <- vector()
    for(i in 1:(nrow(x)-2)){
      A <- as.numeric(x[i,c(2,3,5)])
      B <- as.numeric(x[i+1,c(2,3,5)])
      C <- as.numeric(x[i+2,c(2,3,5)])
      AB <- B-A
      CB <- C-B
      
      #dAB <- sqrt((B[1]-A[1])^2+(B[2]-A[2])^2)
      #dBC <- sqrt((C[1]-B[1])^2+(C[2]-B[2])^2)
      #Formula obtained from https://gitlab.com/anders.sejr.hansen/anisotropy
      angle <- abs(atan2(det(cbind(AB,CB)),AB%*%CB))
      angle <- angle/pi*180
      seg_angle <- c(seg_angle,angle)
      
      
    }
    seg_angle <- c(-1,seg_angle,-1)
    return(seg_angle)
  }
  

  result <- ddply(x,.variables = "track",function(x){
  if(nrow(x)>4){x$angle <- get_angle(x)}else{
    x$angle<- 0
  }
    if(nrow(x)>4){x$angle <- get_angle_3D(x)}else{
      x$angle3D<- 0
    } 
  if(!is.null(x$displacement)){
    x$disp_squared <- x$displacement^2
    x$disp_squared[x$displacement==-1] <- -1
  }
  return(x)
  }
  )
  return(result)
}

SEGMENT_STAT <- function(x){
  UseMethod("SEGMENT_STAT")
}

SEGMENT_STAT.default <- function(x){
  stop("SEGMENT_STAT requires data frame")
}

SEGMENT_STAT.data.frame <-  function(x){
  segment_stat(x)

}

SEGMENT_STAT.list <-  function(x){
  result <- llply(x,function(x){
    SEGMENT_STAT(x)
  })
  return(result)
}
