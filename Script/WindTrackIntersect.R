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


### Track era5 index
tracks$era5ind <- unlist(mclapply(tracks$timestamp, function(x) {
  ind <- which.min(abs(x-fileList$tms))
  if(abs(difftime(x, fileList$tms[ind], units = "hours"))<3) ind else NA
}, mc.cores = 5))

### Track movement vs. residency
lapply(split(tracks, tracks$tag.local.identifier) function(x) {
  
  tm   <- diff(x$timestamp)/60/60
  dist <- geosphere::distm(x[,c("location.long", "location.lat")])[cbind(2:nrow(x), 1:(nrow(x)-1))]*0.001
  
  
  
  hist(tracks$argos.altitude[tracks$argos.altitude>10])
  
  
  plot(x$location.long, x$location.lat)
  plot(x$timestamp, x$argos.altitude)
  
}










