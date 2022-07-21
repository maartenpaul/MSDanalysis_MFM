#200713 make summary table from all data

#from: https://stackoverflow.com/questions/24704344/copy-an-r-data-frame-to-an-excel-spreadsheet
library(clipr)
write.excel <- function(x,row.names=FALSE,col.names=TRUE,...) {
  write_clip(x)
}



segments <- ldply(segments_all)

segments %>%
  mutate(trackID=paste(cellID, track))%>%
  group_by(condition)%>%
  dplyr::summarise(number_of_cells=length(unique(cellID)),number_of_tracks=length(unique(trackID)),number_of_segments=n()) %>%
  write.excel()


cell_statistics <- segments %>%
  group_by(condition,cellID)%>%
  dplyr::summarise(number_of_tracks=length(unique(track)),number_of_segments=n(),trackInside=length(unique(track[trackInMask])),fractionTrackInside=length(unique(track[trackInMask]))/length(unique(track)),
                   segmentsInside=length(inMask[inMask]), fractionSegmentsInside=length(inMask[inMask])/n(), trackInside_gt1=length(unique(track[trackInMask_gt1])),fractionTrackInside_gt1=length(unique(track[trackInMask_gt1]))/length(unique(track)),
                   segmentsInside_gt1=length(inMask_gt1[inMask_gt1]), fractionSegmentsInside_gt1=length(inMask[inMask_gt1])/n(),number_segments_immobile=length(state[state==2]),number_segments_slow=length(state[state==1]),
                   number_segments_fast=length(state[state==0]),
                   medianTracklength=median(table(track)),
                   meanTracklength=mean(table(track)))

merge <- segments %>%
  group_by(condition,cellID,track)%>%
  dplyr::summarise(n_tracks=n()) %>%
 group_by(condition,cellID) %>%
  dplyr::summarise(n1=n(),n5=length(n_tracks[n_tracks>4]),n10=length(n_tracks[n_tracks>9]),n50=length(n_tracks[n_tracks>49]))

cell_statistics <- left_join(cell_statistics,merge,by=c("condition","cellID"))

cell_statistics %>% write.excel()

  cell_statistics%>%  group_by(condition) %>%
  dplyr::summarise(number_of_cells=length(unique(cellID)),total_number_of_tracks=sum(number_of_tracks),total_number_of_segments=sum(number_of_segments),meanNtracks = mean(number_of_tracks),maxNtracks = max(number_of_tracks),minNtracks = min(number_of_tracks),
                   fractionIn=mean(trackInside/number_of_tracks),fractionIn_gt1=mean(trackInside_gt1/number_of_tracks),fractionSegmentsInside=mean(segmentsInside/number_of_segments)) %>% write.excel()
  
  
  
  dplyr::summarise(number_of_cells=length(unique(cellID)))
  
  
  library(plotly)
  x <- ggplot(cell_statistics,aes(x=fractionSegmentsInside,y=number_of_tracks,color=condition))+geom_point()+theme(legend.position = "none")
  ggplotly(x)
  
  x <- ggplot(cell_statistics,aes(x=meanTracklength,y=number_of_tracks,color=condition))+geom_point()+theme(legend.position = "none")
  ggplotly(x)
  