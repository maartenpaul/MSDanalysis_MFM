directory <- "D:/Stack/Genetics/171024 PALB2-Halo tracking/"
directory <- "D:/Stack/Genetics/170809 tracking MMC WT_dDBD_+CTD/"

framerate <- 1/32 #1/200 #ms
pixelsize <- 100 #nm

condition_list <- list.dirs(directory,full.names = F,recursive = F)

sos_to_spoton_csv(condition_list,directory,framerate,pixelsize)

