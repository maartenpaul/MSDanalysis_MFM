library(reticulate)
library(tidyverse)
library(reshape2)

#use with tensorflow==1.14.0 and keras==2.2.0

#py_install(c("statsmodels", "seaborn", "pandas", "matplotlib"))
#py_install("scikit-learn")

path_to_model<-"python/Model_Bidirectional_NoShape_3state_Tr10000"
  rm(py)
  #py$MLfolder <- 'D:/Imaging_data/ML_test/'
  cat("Loading Python functions and ML model")
  source_python(file.path(getwd(),"python/MLE_functions.py"))
  py$path_to_model <- path_to_model
  source_python(file.path(getwd(),"python/Load_model.py"))


ML_segment_tracks <- function(tracks,directory=NULL){

  #source_python("python/Load_data.py")

  x <- dlply(tracks,.variables = c("cellID","track"), function(x){

    return(x$X*10)
  })
  names(x) <- NULL

  py$x <- x


  y <- dlply(tracks,.variables = c("cellID","track"), function(x){
    return(x$Y*10)
  })
  names(y) <- NULL
  py$y <- y

  indices <- cumsum(daply(tracks,.variables = "cellID", function(x){
    return(length(unique(x$track)))
  }))

  indices <- c(0,indices-1)
  cat("import R data")
  source_python("python/Load_Rdata.py")
  cat("Get ML predictions")

  source_python("python/Get_predictions.py")
  cat("Return data")
  #source_python('python/make_table.py')

  if(!is.null(directory)){
    py_save_object(object = x,filename = file.path(directory,tracks$condition[1],"x.pydata"))
    py_save_object(object = y,filename = file.path(directory,tracks$condition[1],"y.pydata"))
    py_save_object(object = allStates,filename = file.path(directory,tracks$condition[1],"allStates.pydata"))



  }

out <-  py$getTrackPiecesForInfo(x, y, allStates)
x0 <- out[[1]]
y0 <- out[[2]]
x1 <- out[[3]]
y1 <- out[[4]]
x2 <- out[[5]]
y2 <- out[[6]]
trnums0 <- out[[7]]
trnums1 <- out[[8]]
trnums2  <- out[[9]]

  names(x0) <- trnums0
  names(x1) <- trnums1
  names(x2) <- trnums2
  names(y0) <- trnums0
  names(y1) <- trnums1
  names(y2) <- trnums2


  x0 <- melt(x0)
  y0 <- melt(y0)
  x0$y <- y0$value
  x0$state <- 0

  x1 <- melt(x1)
  y1 <- melt(y1)
  x1$y <- y1$value
  x1$state <- 1

  x2 <- melt(x2)
  y2 <- melt(y2)
  x2$y <- y2$value
  x2$state <- 2

  all <- rbind(x0,x1,x2)
  all$L1 <- as.numeric(all$L1)
  all <- plyr::arrange(df = all, L1)

  tracks$state <- all$state

  return(tracks)
}
#plotting
#source_python('python/mss.py')
source_python('python/SMMsplot.py')



