

track_msd_fit <- function(x,n=5,fitzero=T,framerate=1/33,pxsize=100,offset=4*0.01^2,dim=2){

  #apply function for every track (column V4)
  coef <- ddply(x[,-6],.variables = "track", function(x) {
    mst <- t(x[,-1])
    time <- seq(1,length(mst),1)
    time <- time/framerate/1000 #convert frame to second
    #plot(mst,ylim=c(0,0.5),xlim=c(0,10))

    if(fitzero){

      mst <- mst-offset

      fit <- tryCatch(lm(formula = mst ~ time-1), error=function(e) NULL)

    } else if (fitzero==FALSE) {

      fit <- tryCatch(nls(formula = mst[,1] ~ time*x1+x2,lower = list(x1=0,x2=0),start=list(x1=0,x2=0),algorithm = "port"), error=function(e) NULL)
    }
    if (!is.null(fit)){
      RSS.p <- sum(residuals(fit)^2)
      TSS <- sum((mst - mean(mst,na.rm = T))^2,na.rm = T)
      rsq <- 1-(RSS.p/TSS)

      coef <- c(coef(fit)[1],offset,rsq)
      coef[1] <- coef[1]/(2*dim)
      return(coef)

    }

    # MSD=4*D*t


  })
  names(coef) <-  c("track","D","intercept","Rsquared")#,paste("dx_",1:n,sep=""),paste("N_",1:n,sep=""))

  return (coef)
}

TRACK_MSD_fit <- function(x,n,fitzero,framerate,pxsize,offset,dim){
  UseMethod("TRACK_MSD_fit")
}

TRACK_MSD_fit.default <- function(x,n=5,fitzero=T,framerate=1/33,pxsize=100,offset=4*0.01^2,dim=2){
  stop("MSD_fit requires data frame")
}

TRACK_MSD_fit.data.frame <-  function(x,n=5,fitzero=T,framerate=1/33,pxsize=100,offset=4*0.01^2,dim=2){
  track_msd_fit(x,n,fitzero,framerate,pxsize,offset,dim)

}

TRACK_MSD_fit.list <-  function(x,n=5,fitzero=T,framerate=1/33,pxsize=100,offset=4*0.01^2,dim=2){
  llply(x,.progress = "text",function(x){
    out <- TRACK_MSD_fit(x,n,fitzero,framerate,pxsize,offset,dim)
    out$cellID <- x$cellID[1]
    return(out)
  })
}

