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
land_map     <- read_sf("~/Google Drive/GeoDat/NaturalEarth/50m_physical/ne_50m_land/ne_50m_land.shp") %>% st_geometry()
ocean_map    <- st_sym_difference(st_as_sfc(st_bbox(c(xmin = -180, xmax = 180, ymax = -90, ymin = 90), crs = st_crs(4326))), st_union(land_map))
ocean_sample <- st_sfc(st_polygon(list(matrix(c(-80, -80, 89, 89, -80,-90, 90, 90, -90,-90), ncol = 2))), crs = 4326)

### tracks
load("~/Google Drive/Science/ProjectsData/MovSim/Tracks/interpTracksWind.rda")
head(interpTracksWind)

  ### wind support
  tracksMap <- cbind(interpTracksWind, t(apply(interpTracksWind[,c(8,9,10,13:20)], 1, function(x) {
    # x <- interpTracksWind[ 262,c(8,9,10,13:20)]
    if(!is.na(x[1]) & x[1]>20) {
      out <- c(RNCEP::NCEP.Airspeed(as.numeric(x[4]), as.numeric(x[5]), as.numeric(x[2]), as.numeric(x[3]))$fa,
               RNCEP::NCEP.Airspeed(as.numeric(x[6]), as.numeric(x[7]), as.numeric(x[2]), as.numeric(x[3]))$fa,
               RNCEP::NCEP.Airspeed(as.numeric(x[8]), as.numeric(x[9]), as.numeric(x[2]), as.numeric(x[3]))$fa,
               RNCEP::NCEP.Airspeed(as.numeric(x[10]), as.numeric(x[11]), as.numeric(x[2]), as.numeric(x[3]))$fa)
      t(data.frame(c(out, max(out, na.rm = T))))
    } else rep(NA, 5)
  })))
  tracksMap$id <- sapply(strsplit(row.names(tracksMap), "[.]"), function(x) x[[1]])

  
# plot(land_map, xlim = c(160, 200))  
# points(tracksMap$location.long, tracksMap$location.lat, pch = 16, col = ifelse(is.infinite(abs(tracksMap$'5')), "grey90", "orange"))

############
### loop ###
############
  
### time sequence
start <- as.POSIXct("2020-03-01", tz = "GMT")
end   <- as.POSIXct("2020-04-14", tz = "GMT")
  
tmTab <- subset(tracksMap, timestamp>=start & timestamp<=end & !duplicated(timestamp), select = c("timestamp", "file", "fileID"))
  
ERAindex <- cbind(tracksMap$file, tracksMap$fileID)[!duplicated(paste(tracksMap$file, tracksMap$fileID)),]
TIFFindex <- data.frame(tms = tmTab$timestamp, ind = apply(tmTab[,c("file", "fileID")], 1, function(x) which(x[1]==ERAindex[,1] & x[2]==ERAindex[,2]))) 
  
  
## col inds
cls <- brewer.pal(11, "RdYlGn")
tracksMap$colID <- cut(tracksMap$'5', seq(-26, 26, length = length(cls)), labels = FALSE)


## global colors
blues  <- colorRampPalette(pal_seeblau)
greens <- colorRampPalette(pal_seegruen)

## anlge init
rot <- data.frame(long  = seq(130, 145, length = nrow(TIFFindex)),
                  fromx = seq(3,     0, length = nrow(TIFFindex)),
                  fromy = seq(7,     9, length = nrow(TIFFindex)),
                  fromz = seq(10,  8.5, length = nrow(TIFFindex)))


## wind track init
nrTracks <-  50

for(i in 1:nrow(TIFFindex)) {
  
  if(i==1 || TIFFindex$ind[i-1]!=TIFFindex$ind[i]) {
    
    tmpBrick <- brick(paste0("~/Google Drive/Science/ProjectsData/MovSim/geoTiffs/WindSnow_", TIFFindex[i,2], ".tif"))
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
  new[1,,] <- t(st_coordinates(suppressMessages(st_sample(ocean_map %>% st_difference(ocean_sample), size = dim(new)[3]))))
   
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
      path(x = 0, y = 7, material = diffuse(color="grey80"), width = 0.002, angle=c(180,rot[i,1],0))
  }, mc.cores = 4))

  
  ### Bar-tailed godwit tracks
  realTracks <- subset(tracksMap, timestamp<=TIFFindex[i,1])
  last       <- lapply(unique(realTracks$id), function(x) realTracks[realTracks$id==x,][sum(realTracks$id==x),])

  t_real <- st_wrap_dateline(st_sf(st_sfc(st_multipoint(as.matrix(do.call("rbind", last)[,c("location.long", "location.lat")]))), crs = 4326))
   
  grPoints <- dplyr::bind_rows(parallel::mclapply(split(as.data.frame(cbind(st_coordinates(t_real)[,1:2], id = 1:nrow(st_coordinates(t_real)), colID = do.call("rbind", last)$colID)),
                                                        1:nrow(st_coordinates(t_real))),
          function(x) {
              sp <- sph2car(matrix(c((x$X+rot[i,1])*-1, x$Y, 1), ncol = 3))[,c(1,3,2)]
              sphere(x = sp[1], y = 7+sp[2], z = sp[3], material = diffuse(color=ifelse(is.na(x$colID), "grey30", cls[x$colID])), radius = 0.02)
            }, mc.cores = 4))
  
  
  generate_studio(material=diffuse(color = "grey10")) %>%
    add_object(sphere(x = 0, y = 7, radius=0.9999, angle = c(0, rot[i,1], 0),
                      material = glossy(gloss=0.3, image_texture =  img))) %>%
    add_object(group_objects(grTracks)) %>%
    add_object(group_objects(grPoints)) %>%
    add_object(sphere(y = 9, z = 10, x = 20, radius = 6, material=light(intensity=10))) %>%
    render_scene(width=350, height=350, aperture=0, fov=14, sample_method = "random", parallel = TRUE,
                 samples= 200, clamp_value=10, lookfrom=c(rot[i,2],rot[i,3],rot[i,4]), lookat=c(0,7,0), camera_up = c(0,1,0), filename=glue::glue("~/Desktop/tmp/globe_{i}.png"))
  
  
  title_mat = matrix(0,60,350) %>%
    add_title(title_text = glue::glue("Title:  {format(TIFFindex[i,1], '%Y-%m-%d')}"), 
              title_bar_alpha = 1, title_bar_color = "grey30", title_size = 20,
              title_color = "white", filename=glue::glue("~/Desktop/tmp/title_{i}.png"))
  
  
  world_image <- magick::image_read(glue::glue("~/Desktop/tmp/globe_{i}.png"))
  titel_image <- magick::image_read(glue::glue("~/Desktop/tmp/title_{i}.png"))
  
  magick::image_append(c(titel_image, world_image), stack = TRUE) %>% 
    magick::image_write(glue::glue("~/Google Drive/tmp/full_image_{i}.png"))
  
  file.remove(c(glue::glue("~/Desktop/tmp/globe_{i}.png"), glue::glue("~/Desktop/tmp/title_{i}.png")))
}

# av::av_encode_video(glue::glue("~/Google Drive//tmp/full_image_{1:121}.png"), 
                    # output = "~/Google Drive//globe_viz.mp4", framerate = 15)




