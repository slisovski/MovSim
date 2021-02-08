library(raster)
library(unikn)
library(abind)
library(dichromat)
library(rayrender)
library(sphereplot)
library(sf)
library(rayimage)

### global elevation model
topo <- raster("~/Google Drive/GeoDat/ETOPO1_Ice_g_geotiff.tif")

### Map
land_map   <- read_sf("~/Google Drive/GeoDat/NaturalEarth/50m_physical/ne_50m_land/ne_50m_land.shp") %>% st_geometry()
ocean_map  <- st_sym_difference(st_as_sfc(st_bbox(c(xmin = -180, xmax = 180, ymax = -90, ymin = 90), crs = st_crs(4326))), st_union(land_map))


### tracks
load("~/Google Drive/Science/ProjectsData/MovSim/Tracks/interpTracksWind.rda")
head(interpTracksWind)



### time sequence
start <- as.POSIXct("2020-03-01", tz = "GMT")
end   <- as.POSIXct("2020-03-15", tz = "GMT")

tmTab <- subset(interpTracks, tms>=start & tms<=end & !duplicated(tms), select = c("tms", "file", "fileID"))

ERAindex <- cbind(interpTracksWind$file, interpTracksWind$fileID)[!duplicated(paste(interpTracksWind$file, interpTracksWind$fileID)),]
ERAindex <- ERAindex[order(ERAindex[,1]),]

TIFFindex <- data.frame(tms = tmTab$tms, ind = apply(tmTab[,c("file", "fileID")], 1, function(x) which(x[1]==ERAindex[,1] & x[2]==ERAindex[,2]))) 

### sun zenith
rSun    <- raster(xmn = -180, xmx = 180, ymn = -90, ymx = 90, res = 1)
TIFFindex <- cbind(TIFFindex, t(apply(TIFFindex, 1, function(x) {
  solar   <- SGAT::solar(as.POSIXct(paste0(format(as.POSIXct(x[1], origin = "1970-01-01", tz = "GMT"), "%Y-%m-%d"), " 12:00:00"), tz = "GMT"))
  crds    <- coordinates(rSun)
  angle   <- unlist(parallel::mclapply(1:nrow(crds), function(x) SGAT::zenith(solar, crds[x,1], crds[x,2]), mc.cores = 4))
  sun     <- crds[which.min(abs(angle-0)),]
  sinCrds <- (sph2car(cbind(-sun[1], sun[2]),1)[,c(1,3,2)]*9)
})))



############
### loop ###
############

## globals
blues  <- colorRampPalette(pal_seeblau)
greens <- colorRampPalette(pal_seegruen)

## anlge init
rot <- cbind(seq(0, -40, length = nrow(TIFFindex)), 
             seq(0, 125, length = nrow(TIFFindex)))

## wind track init
nrTracks <-  1000

for(i in 1:nrow(TIFFindex)) {
  
  if(i==1 || TIFFindex$ind[i-1]!=TIFFindex$ind[i]) {
    
    tmpBrick <- brick(paste0("~/Google Drive/Science/ProjectsData/MovSim/geoTiffs/WindSnow_", i, ".tif"))
    spdR <- sqrt(tmpBrick[[1]]^2 + tmpBrick[[2]]^2)
    
    top  <- spdR; top[] <- NA
    crds <- coordinates(top)[is.na(tmpBrick[[1]][]) & is.na(tmpBrick[[3]][]),]
    top[is.na(tmpBrick[[1]][]) & is.na(tmpBrick[[3]][])] <- raster::extract(topo, crds)
    top[coordinates(top)[,2]< -57] <- NA
    
    ocean  <- as.array(RGB(spdR, col = blues(100)))
    land   <- as.array(RGB(top,  col = greens(100)))
    
    ocean[land!=255]  <- land[land!=255]
    img <- abind(lapply(1:3, function(x) t(ocean[,dim(ocean)[2]:1,x]/255)), along = 3)
    
  }
  
  ## tracks
  new <- array(dim = c(8, 2, nrTracks))
  probWind <- spdR[!is.na(spdR[])]
  new[1,,] <- t(coordinates(spdR)[!is.na(spdR[]),][sample(1:sum(!is.na(spdR[])), nrTracks),]) #, 
                                                          # prob = ((probWind-min(probWind)) / (max(probWind) - min(probWind))) * (1-0.75)*0.75),])
  
  if(i==1) tracks <- new else tracks <- abind(tracks, new, along = 3)
  
  ind <- do.call("rbind", parallel::mclapply(1:dim(tracks)[3], function(x) {
    ind   <- max(which(!is.na(tracks[,1,x])))
    matrix(c(ind, tracks[ind,,x]), ncol = 3)
  }, mc.cores = parallel::detectCores()))
  
  del <- which(ind[,1]>=(dim(tracks)[1]-1))
  if(length(del)>0) {
    tracks <- tracks[,,-del]
    ind <- ind[-del,]
  }
  
  wnd  <- extract(tmpBrick[[1:2]], ind[,2:3])
  dir  <- atan2(wnd[,1], wnd[,2]) * (180/pi)
  spd  <- sqrt(wnd[,1]^2 + wnd[,2]^2)
  
  dest <- geosphere::destPoint(ind[,2:3], dir, spd*60*60*5)
  
  invisible(mapply(function(x,y) {tracks[x,,y] <<- dest[y,]}, ind[,1]+1, 1:dim(tracks)[3]))
  
  ls <- lapply(1:dim(tracks)[3], function(x) {
    tmp <- tracks[!is.na(tracks[,1,x]),,x]
    if(is.matrix(tmp))   st_linestring(tmp)
  })
  
  t <- suppressMessages(st_wrap_dateline(st_sf(st_sfc(st_multilinestring(ls[!sapply(ls, is.null)]), crs = 4326))) %>% st_intersection(ocean_map))
  
  grTracks <- dplyr::bind_rows(parallel::mclapply(split(as.data.frame(st_coordinates(t)[,1:2]), st_coordinates(t)[,3]), function(x) {
    x[,2] <- x[,2]*-1
    sph2car(cbind(x,1))[,c(1,3,2)] %>%
      path(x = 0, y = 7, material = diffuse(color="grey80"), width = 0.002, angle=c(180+rot[i,1],rot[i,2],0))
  }, mc.cores = 4))
  
  
  generate_studio(material=diffuse(color = "grey10")) %>%
    add_object(sphere(x = 0, y = 7, radius=0.9999, angle=c(rot[i,1],rot[i,2],0),
                      material = glossy(gloss=0.3, image_texture =  img))) %>%
    add_object(group_objects(grTracks)) %>%
    add_object(sphere(y = 9, z = 10, x = 20, radius = 6, material=light(intensity=10))) %>%
    render_scene(width=450, height=450, aperture=0, fov=14, sample_method = "random", parallel = TRUE,
                 samples= 200, clamp_value=10, lookfrom=c(1,7,10), lookat=c(0,7,0), camera_up = c(0,1,0), filename=glue::glue("~/Desktop/tmp/globe_{i}.png"))
  
  
  title_mat = matrix(0,60,450) %>%
    add_title(title_text = glue::glue("Bar-tailed godwit migration          {format(TIFFindex[i,1], '%Y-%m-%d')}"), 
              title_bar_alpha = 1, title_bar_color = "grey30", title_size = 20,
              title_color = "white", filename=glue::glue("~/Desktop/tmp/title_{i}.png"))
  
  
  world_image <- magick::image_read(glue::glue("~/Desktop/tmp/globe_{i}.png"))
  titel_image <- magick::image_read(glue::glue("~/Desktop/tmp/title_{i}.png"))
  
  magick::image_append(c(titel_image, world_image), stack = TRUE) %>% 
    magick::image_write(glue::glue("~/Desktop/tmp/full_image_{i}.png"))
  
  file.remove(c(glue::glue("~/Desktop/tmp/globe_{i}.png"), glue::glue("~/Desktop/tmp/title_{i}.png")))
}

av::av_encode_video(glue::glue("~/Desktop/tmp/full_image_{1:225}.png"), 
                    output = "~/Desktop/globe_viz.mp4", framerate = 15)


# generate_studio(material=diffuse(color = "grey10")) %>%
#  add_object(sphere(x = 0, y = 1, radius=0.9999, angle=c(0,0,0),
#                    material = glossy(gloss=0.3, image_texture =  img))) %>%
#   add_object(sphere(y = 9, z = 10, x = 20, radius = 6, material=light(intensity=10))) %>%
#   render_scene(width=1200, height=1200, aperture = 0, fov = 15, sample_method = "stratified", parallel = TRUE,
#               samples= 4000, clamp_value=10, lookat = c(0,1,0), lookfrom = c(0,1,10), camera_up = c(0,1,0))



