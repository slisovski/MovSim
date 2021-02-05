library(sf)
library(raster)
library(ncdf4)
library(parallel)
library(rgdal)
library(rgeos)

### Tracks
tracks <- read.csv("~/Google Drive/Science/ProjectsData/MovSim/Tracks/trackMovebank.csv")
tracks$timestamp <- as.POSIXct(tracks$timestamp, tz = "GMT")


### Map
land   <- read_sf("~/Google Drive/GeoDat/NaturalEarth/50m_physical/ne_50m_land/ne_50m_land.shp") %>% st_geometry()
ocean  <- st_sym_difference(st_as_sfc(st_bbox(c(xmin = -180, xmax = 180, ymax = -90, ymin = 90), crs = st_crs(4326))), st_union(land))


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
  
  hourTrack                    <- data.frame(id = x$tag.local.identifier[1], tms = tms, ttAll[unlist(mclapply(as.numeric(tms), function(y) which.min(abs(y-ttAll[,1])), mc.cores = 3)),-1])
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


movFreq <- aggregate(interpTracks$mov, by = list(interpTracks$ERAind), FUN = function(x) sum(x)/length(x))





## snow Cover (24km)
fls.gz <- list.files("/Volumes/slisovski/RemoteSensedData/IMS_DailyNHSnowIceAnalysis/24km/2020/", pattern = ".asc.gz", recursive = T,  full.names = T)
dates  <- as.Date(as.POSIXct(unlist(lapply(strsplit(fls.gz, "ims"), function(x) strsplit(x[[2]], "_24km")))[c(TRUE, FALSE)], format = "%Y%j"))
prj <- "+proj=stere +lat_0=90 +lat_ts=60 +lon_0=-80 +k=1 +x_0=0 +y_0=0 +a=6378137 +b=6356257 +units=m +no_defs"


## mkTiles and wind properties per track

# m10 dates
m10path  <- "~/Google Drive/Science/ProjectsData/MovSim/ERA5/era5_wind_10m.nc"
nf       <- nc_open(m10path)
tms10    <- as.POSIXct(nf$var[[1]]$dim[[3]]$vals*60*60, "1900-01-01", tz = "GMT")
nc_close(nf)
tms10ind <- mapply(function(x) which.min(abs(x-tms10)), x = interpTracks$tms[!duplicated(interpTracks$ERAind)])

interpTracksWind <- cbind(interpTracks, matrix(ncol = 8, nrow = nrow(interpTracks)))
names(interpTracksWind) <- c(names(interpTracks), c("m10u", "m10v", "p700u", "p700v", "p850u", "p850v", "p925u", "p925v"))

Date     <- NA

for(i in 1:max(interpTracks$ERAind)) {
  
  cat(sprintf('\rDate %d of %d',
              i, max(interpTracks$ERAind)))
  
  path <- as.character(fileList$path[i])
  
  windR  <- rotate(stack(brick(m10path, varname = "u10", level = 1)[[i]],
                         brick(m10path, varname = "v10", level = 1)[[i]],
                         brick(path, varname =  "u", level = 1)[[i]],
                         brick(path, varname =  "v", level = 1)[[i]],
                         brick(path, varname =  "u", level = 2)[[i]],
                         brick(path, varname =  "v", level = 2)[[i]],
                         brick(path, varname =  "u", level = 3)[[i]],
                         brick(path, varname =  "v", level = 3)[[i]]))
  names(windR) <- c("m10u", "m10v", "p700u", "p700v", "p850u", "p850v", "p925u", "p925v")
  
  trackInd <- which(interpTracks$ERAind==i)
  if(length(trackInd)>0) {
  interpTracksWind[trackInd, which(names(interpTracksWind)%in%names(windR))] <- raster::extract(windR, interpTracks[trackInd, c("location.long", "location.lat")])
  }
  
  if(is.na(Date) | Date !=as.Date(median(interpTracks$tms[trackInd]))) {
    Date <- as.Date(median(interpTracks$tms[trackInd]))
    indSnow <- which(dates==Date)
    
    tab0 <- readLines(fls.gz[indSnow])
    ind <- unlist(suppressWarnings(parallel::mclapply(tab0, function(x) is.na(as.numeric(gsub(" ", "", x))), mc.cores = 5)))
    tab <- tab0[-which(ind)]
    
    z = do.call("rbind", parallel::mclapply(tab, function(.line) as.numeric(strsplit(.line, '')[[1]]), mc.cores = 5))
    r0 <- raster(z[nrow(z):1,], crs = CRS(prj))
    extent(r0) <- c(-12126597.0, -12126597.0 + 1024*23684.997, -12126840.0, -12126597.0 + 1024*23684.997)
  }
  
  uOut <- windR[[1]]
  uOut[] <- apply(windR[[c(1,3,5)]][], 1, function(x) median(x, na.rm = T))
  uOut <- mask(uOut, as(ocean, "Spatial"))
  
  vOut <- windR[[1]]
  vOut[] <- apply(windR[[c(2,4,6)]][], 1, function(x) median(x, na.rm = T))
  vOut <- mask(vOut, as(ocean, "Spatial"))
  
  snowe <- raster::extract(r0, project(coordinates(snow)[is.na(vOut[]),], prj))
  snow  <- vOut; snow[] <- NA; snow[is.na(vOut[])] <- ifelse(snowe==4, 1, NA)

  dateBrick <- brick(uOut, vOut, snow)
  
  writeRaster(dateBrick, paste0("~/Google Drive/Science/ProjectsData/MovSim/geoTiffs/WindSnow_", i, ".tif"), options = c('TFW=YES'), overwrite=TRUE)
  
}

save(interpTracksWind, file = "~/Google Drive/Science/ProjectsData/MovSim/Tracks/interpTracksWind.rda")

