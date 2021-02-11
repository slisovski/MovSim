library(sf)
library(raster)
library(ncdf4)
library(parallel)
library(rgdal)
library(rgeos)
library(future.apply)
plan(multisession)

### Tracks
tracks <- read.csv("~/Google Drive/Science/ProjectsData/MovSim/Tracks/trackMovebank.csv")
tracks$timestamp <- as.POSIXct(tracks$timestamp, tz = "GMT")
tracks <- subset(tracks, !tag.local.identifier%in%c(195831, 195835, 195840, 195842))

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
  
  data.frame(path = x, id = sapply(strsplit(x, "//"), function(y) y[[2]]), tms = tms, file = which(fls==x), fileID = 1:length(tms))
}))


### Track interpolation & movement vs. residency
interpTracks <- do.call("rbind", lapply(split(tracks, tracks$tag.local.identifier), function(x) {
  
  tms = seq(min(fileList$tms), max(fileList$tms), by = 1.5*60*60)
  
  indTms <- do.call("rbind", lapply(tms, function(y) {
    dfx <- abs(difftime(y, x$timestamp, units = "hours"))
    data.frame(tms = x$timestamp[which.min(dfx)], diff = as.numeric(dfx[which.min(dfx)]), ind = which.min(dfx))
  }))
  
  indTms$min <- unlist(sapply(unique(indTms$ind), function(y) {
    tmp <- indTms[indTms$ind==y,]
    out <- rep(NA, nrow(tmp))
    out[which.min(tmp$diff)] <- y
    out
  }))
  
  hourTrack <- data.frame(timestamp = tms, indTrack = indTms$min, location.long = NA, location.lat = NA, argos.altitude = NA) 
  hourTrack[which(!is.na(hourTrack$indTrack)),c("location.long", "location.lat", "argos.altitude")] <- x[hourTrack$indTrack[!is.na(hourTrack$indTrack)], c("location.long", "location.lat", "argos.altitude")]
  hourTrack$location.long      <- ifelse(hourTrack$location.long<0, 180 + (hourTrack$location.long + 180), hourTrack$location.long)
  
  ### interpolation
  hourTrackInterp <- do.call("rbind", lapply(unique(hourTrack$indTrack[!is.na(hourTrack$indTrack)]), function(z) {
    
    if(sum(hourTrack$indTrack>z, na.rm = T)){
      end   <- min(which(hourTrack$indTrack>z))
      tmp   <- hourTrack[which(!is.na(hourTrack$indTrack) & hourTrack$indTrack==z):end,]
      diff  <- as.numeric(difftime(tmp$timestamp[nrow(tmp)], tmp$timestamp[1], units = "hours"))
      tmp[,c("location.long", "location.lat")] <- apply(tmp[,c("location.long", "location.lat")], 2, function(a) zoo::na.approx(a, method = ifelse(diff>24*3, "constant", "linear")))
      cbind(tmp[-nrow(tmp),], type = ifelse(diff>24*3, 2, 1))
    } else {
      cbind(hourTrack[which(!is.na(hourTrack$indTrack) & hourTrack$indTrack==z):nrow(hourTrack),], type = 3)
    }
    
  }))
  
  
  # plot(land, xlim = c(160, 200))
  # lines(hourTrackInterp$location.long, hourTrackInterp$location.lat, col = "darkgreen")
  # points(hourTrackInterp$location.long, hourTrackInterp$location.lat, pch = 16, cex = ifelse(hourTrackInterp$type==1, 0.4, 1), 
  #        col = ifelse(hourTrackInterp$type==1, "cornflowerblue", "purple"))  
  
  hourTrackInterp$location.long.wrap <- ifelse(hourTrackInterp$location.long>180, -180 + (hourTrackInterp$location.long - 180), hourTrackInterp$location.long)
  hourTrackInterp$dist               <- c(sapply(1:(nrow(hourTrackInterp)-1), function(z) geosphere::distVincentySphere(hourTrackInterp[z, c(7,4)], hourTrackInterp[z+1, c(7,4)])/1000), 0)
  hourTrackInterp$bearing            <- c(sapply(1:(nrow(hourTrackInterp)-1), function(z) geosphere::bearing(hourTrackInterp[z, c(7,4)], hourTrackInterp[z+1, c(7,4)])), 0)
  hourTrackInterp$speed              <- (hourTrackInterp$dist*1000)/c(as.numeric(diff(hourTrackInterp$timestamp)*60*60), NA)
  
  hourTrackInterp <- hourTrackInterp[,c(1,2,6,3,7,4,5,8,9,10)]

  hourTrackInterp
  
}))

interpTracks <- cbind(interpTracks, t(mapply(function(t) {
  ind <- which.min(abs(t - fileList$tms))
  td  <- abs(difftime(t, fileList$tms, units = "hours"))[ind]
  if(td>5) {
    c(NA, NA)
  } else unlist(fileList[ind, c("file", "fileID")])
}, t = interpTracks$timestamp)))

# save(interpTracks, file = "~/Google Drive/Science/ProjectsData/MovSim/Tracks/interpTracks.rda")
load("~/Google Drive/Science/ProjectsData/MovSim/Tracks/interpTracks.rda")

interpTracksWind <- cbind(interpTracks, matrix(ncol = 8, nrow = nrow(interpTracks)))
names(interpTracksWind) <- c(names(interpTracks), c("m10u", "m10v", "p700u", "p700v", "p850u", "p850v", "p925u", "p925v"))


## snow Cover (24km)
fls.gz <- list.files("/Volumes/slisovski/RemoteSensedData/IMS_DailyNHSnowIceAnalysis/24km/2020/", pattern = ".asc.gz", recursive = T,  full.names = T)
datesSnow  <- as.POSIXct(sapply(strsplit(fls.gz, "ims"), function(x) sapply(strsplit(x[[2]], "_24km"), function(z) z[[1]])), format = "%Y%j", tz = "GMT")
prj <- "+proj=stere +lat_0=90 +lat_ts=60 +lon_0=-80 +k=1 +x_0=0 +y_0=0 +a=6378137 +b=6356257 +units=m +no_defs"


############################################
## mkTiles and wind properties per track ###
############################################

# m10 dates
m10path  <- "~/Google Drive/Science/ProjectsData/MovSim/ERA5/era5_wind_10m.nc"
nf       <- nc_open(m10path)
tms10    <- as.POSIXct(nf$var[[1]]$dim[[3]]$vals*60*60, "1900-01-01", tz = "GMT")
nc_close(nf)
tms10ind <- mapply(function(x) which.min(abs(x-tms10)), x = interpTracks$tms[!duplicated(interpTracks$ERAind)])

ERAindex      <- data.frame(file = interpTracksWind$file, ind = interpTracksWind$fileID)[!duplicated(paste(interpTracksWind$file, interpTracksWind$fileID)),]
ERAindex$tms  <- as.POSIXct(apply(ERAindex, 1, function(x) fileList$tms[fileList$file==x[1] & fileList$fileID==x[2]]), origin = "1970-01-01", tz = "GMT")
ERAindex$snow <- sapply(ERAindex$tms, function(x) which.min(abs(x-datesSnow)))

for(i in i:nrow(ERAindex)) {
  
  cat(sprintf('\rDate %d of %d',
              i, nrow(ERAindex)))
  
  path <- fls[ERAindex[i,1]]
  
  trackInd <- which(interpTracksWind$file==ERAindex[i,1] & interpTracksWind$fileID==ERAindex[i,2])
  
  windR  <- rotate(stack(brick(m10path, varname = "u10", level = 1)[[i]],
                         brick(m10path, varname = "v10", level = 1)[[i]],
                         brick(path, varname =  "u", level = 1)[[ERAindex[i,2]]],
                         brick(path, varname =  "v", level = 1)[[ERAindex[i,2]]],
                         brick(path, varname =  "u", level = 2)[[ERAindex[i,2]]],
                         brick(path, varname =  "v", level = 2)[[ERAindex[i,2]]],
                         brick(path, varname =  "u", level = 3)[[ERAindex[i,2]]],
                         brick(path, varname =  "v", level = 3)[[ERAindex[i,2]]]))
  names(windR) <- c("m10u", "m10v", "p700u", "p700v", "p850u", "p850v", "p925u", "p925v")
  
  
  if(length(trackInd)>0) {
    interpTracksWind[trackInd, which(names(interpTracksWind)%in%names(windR))] <- raster::extract(windR, interpTracks[trackInd, c("location.long", "location.lat")])
  }

    tab0    <- readLines(fls.gz[ERAindex$snow[i]])
    ind     <- unlist(suppressWarnings(parallel::mclapply(tab0, function(x) is.na(as.numeric(gsub(" ", "", x))), mc.cores = 5)))
    tab     <- tab0[-which(ind)]

    z   <-  do.call("rbind", parallel::mclapply(tab, function(.line) as.numeric(strsplit(.line, '')[[1]]), mc.cores = 5))
    r0  <- raster(z[nrow(z):1,], crs = CRS(prj))
    extent(r0) <- c(-12126597.0, -12126597.0 + 1024*23684.997, -12126840.0, -12126597.0 + 1024*23684.997)

  msk  <- mask(windR, as(ocean, "Spatial"))
  uOut <- msk[[1]]; uOut[] <- future_apply(msk[[c(1,3,5)]][], MARGIN = 1L, FUN = median)
  vOut <- msk[[1]]; vOut[] <- future_apply(msk[[c(2,4,6)]][], MARGIN = 1L, FUN = median)

  snowe <- raster::extract(r0, project(coordinates(snow)[is.na(vOut[]),], prj))
  snow  <- vOut; snow[] <- NA; snow[is.na(vOut[])] <- ifelse(snowe==4, 1, NA)

  dateBrick <- brick(uOut, vOut, snow)

  writeRaster(dateBrick, paste0("~/Google Drive/Science/ProjectsData/MovSim/geoTiffs/WindSnow_", i, ".tif"), options = c('TFW=YES'), overwrite=TRUE)

  if(i/200 == floor(i/200) | i==nrow(ERAindex)) save(interpTracksWind, file = "~/Google Drive/Science/ProjectsData/MovSim/Tracks/interpTracksWind.rda")

   
}


