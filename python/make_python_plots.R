library(reticulate)
source("scripts/ML_py.R")
source_python('python/mss.py')
source_python('python/SMMsplot.py')
ML_load()



condition = "dCTD MMC"
x <- py_load_object(file.path("D:/Imaging_data/MFM_test_data/x.pydata"))
y <- py_load_object(file.path("D:/Imaging_data/MFM_test_data/y.pydata"))
allStates <- py_load_object(file.path("D:/Imaging_data/MFM_test_data/allStates.pydata"))

out <- dlply(segments_all[[condition]],.variables = c("cellID","track"))


inmask <- laply(out,function(x){
  any(x$inMask)
})


x2 <- x[inmask]
y2 <- y[inmask]
allStates2 <- allStates[inmask]
x3 <- x[!inmask]
y3 <- y[!inmask]
allStates3 <- allStates[!inmask]

py_run_string("pixSize= 0.1")
py_run_string("t = 0.05")



#mss(x,y,allStates)
plt <- SMMsplot(x2,y2,allStates2)

py_run_string("plt.show()")
matplt <- import("matplotlib")
test <- plot
test$show

##
inside <- getTrackPieces(x2,y2,allStates2)
lengths <- laply(inside[[5]],length)/10
lenghts_df <- data.frame("lengths"=lengths)
ggplot(lenghts_df, aes(x = lengths))+ geom_histogram()  +xlim(c(-1,20))+xlab("Track Length (s)")
outside <- getTrackPieces(x3,y3,allStates3)
lengths_out <- laply(outside[[5]],length)
lenghts_out_df <- data.frame("lengths"=lengths_out/10,"inMask"=FALSE)
ggplot(rbind(lenghts_df,lenghts_out_df), aes(x = lengths,y=..density..,fill=inMask)) +geom_histogram(position="identity",alpha=0.35)  +xlab("Track Length (s)")+xlim(c(0,10))



ggplot(data = msd_fit_all$`dDBD HU`,aes(x=D,y=..density..,color=inMask)) +
  geom_histogram(position="identity",fill="white",alpha=0.1) + scale_x_log10(limits=c(0.00005,2))

