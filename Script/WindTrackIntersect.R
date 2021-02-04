library(sf)
library(raster)
library(ncdf4)
library(parallel)


### Tracks
tracks <- read.csv("~/Google Drive/Science/ProjectsData/MovSim/Tracks/trackMovebank.csv")
tracks$timestamp <- as.POSIXct(tracks$timestamp, tz = "GMT")


### ERA5 dataset
fls <- list.files("~/Google Drive/Science/ProjectsData/MovSim/ERA5/", "adaptor", full.names = T)

fileList <- do.call("rbind", lapply(fls, function(x) {
  id = sapply(strsplit(x, "//"), function(y) y[[2]])
  nf <- nc_open(x)
  tms    <- as.POSIXct(nf$var[[1]]$dim[[4]]$vals*60*60, "1900-01-01", tz = "GMT")
  nc_close(nf)
  
  data.frame(path = x, id = sapply(strsplit(x, "//"), function(y) y[[2]]), tms = tms)
}))


### Track interpolation & movement vs. residency
interpTracks <- do.call("rbind", lapply(split(tracks, tracks$tag.local.identifier), function(x) {
  
  tms = seq(min(fileList$tms), max(fileList$tms), by = 1.5*60*60)
  
    tt      <- merge(data.frame(timestamp = seq(min(x$timestamp), max(x$timestamp), by = 1)), x[,c("timestamp", "location.long", "location.lat", "argos.altitude")], all.x = T)
    tt$location.long <- ifelse(tt$location.long<0, 180 + (tt$location.long + 180), tt$location.long)
    ttAll   <- cbind(tt$timestamp, apply(tt[,-1], 2, zoo::na.approx, rule = 2)) 
  
  hourTrack                    <- data.frame(tms = tms, ttAll[unlist(mclapply(as.numeric(tms), function(y) which.min(abs(y-ttAll[,1])), mc.cores = 3)),-1])
  hourTrack$location.long.wrap <- ifelse(hourTrack$location.long>180, -180 + (hourTrack$location.long - 180), hourTrack$location.long)
  hourTrack$dist               <- c(sapply(1:(nrow(hourTrack)-1), function(z) geosphere::distVincentySphere(hourTrack[z, c(5,3)], hourTrack[z+1, c(5,3)])/1000), 0)
  hourTrack$bearing            <- c(sapply(1:(nrow(hourTrack)-1), function(z) geosphere::bearing(hourTrack[z, c(5,3)], hourTrack[z+1, c(5,3)])/1000), 0)
  hourTrack$speed              <- (hourTrack$dist*1000)/c(as.numeric(diff(hourTrack$tms)*60*60), NA)
  hourTrack$mov <- hourTrack$dist>25
  
  hourTrack <- hourTrack[,c(1,2,5,3,4,6,7,8,9)]
  
  # plot(hourTrack$location.long, hourTrack$location.lat, pch = 16, col = ifelse(hourTrack$dist>25, "firebrick", "grey90"))
  # plot(wrld_simpl, add = T)
  # plot(maptools::elide(wrld_simpl, shift = c(360, 0)), add = T)
  
  hourTrack$ERAind <- mapply(function(i) {
    ind <- which.min(abs(i - fileList$tms))
    if(as.numeric(abs(difftime(i, fileList$tms[ind], units = "hours")))>5) {
      NA
    } else ind
  }, i = hourTrack$tms)
  
  hourTrack
  
}))

save(interpTracks, file = "~/Google Drive/Science/ProjectsData/MovSim/Tracks/interpTracks.rda")










