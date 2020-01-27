directory <- "D:/Stack/Genetics/171024 PALB2-Halo tracking/"

condition_list <- list.dirs(directory,full.names = F,recursive = F)
#condition_list <- condition_list[c(2,3)]
#condition_list <- "SNAP-SiR"
area_all <- list()
#load data
for (i in 1:length(condition_list)){
  segments <- list()
  dir <- file.path(directory,condition_list[i])
  filelist <- list.dirs(dir,full.names = T,recursive = F)
  #filelist <- filelist[-grep("skip",x = filelist)]
  total <- length(filelist)
  # create progress bar
  area <- vector(length=length(total))
  for(j in 1:total){
    #  Sys.sleep(0.1)
  rois <- read.csv(file.path(filelist[j],"Mask_results.txt"),sep = "\t",header = T)
  area[j] <- rois$Area[1]
  }
  area_all[[condition_list[i]]] <- area

}

