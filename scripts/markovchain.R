library(markovchain)
library(plyr)
markovchain::markovchainFit()

statesequence

statesequence <- llply(segments_all,function(x) {
  dlply(subset(x,trackInMask==TRUE),.variables = c("cellID","track"),function(x){

    return(x$inMask)

  })
})


statesequence <- llply(segments_all[c(4,6)],function(x) {
  dlply(x,.variables = c("cellID","track"),function(x){

    return(x$inMask)

  })
})


markovmatrix <- llply(statesequence, function(x){
  markovchainFit(x,confint=FALSE)
})

save(markovmatrix,file="markovchain.Rdata")
load(file="markovchain.Rdata")



save(segments_all,file="ML_segments_all.Rdata")


load(file="ML_segments_all.Rdata")
