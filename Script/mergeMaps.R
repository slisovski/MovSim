library(raster)
library(unikn)
library(abind)
library(dichromat)
library(rayrender)
library(sphereplot)
library(sf)
library(rayimage)
library(RColorBrewer)

### global elevation model
topo <- raster("GeoDat/ETOPO1_Ice_g_geotiff.tif")

### Map
land_map     <- read_sf("GeoDat/NaturalEarth/50m_physical/ne_50m_land/ne_50m_land.shp") %>% st_geometry()
ocean_map    <- st_sym_difference(st_as_sfc(st_bbox(c(xmin = -180, xmax = 180, ymax = -90, ymin = 90), crs = st_crs(4326))), st_union(land_map))
ocean_sample <- st_sfc(st_polygon(list(matrix(c(-80, -80, 89, 89, -80,-90, 90, 90, -90,-90), ncol = 2))), crs = 4326)

### tracks
load("Tracks/interpTracksWind.rda")

  ### wind support
  tracksMap <- cbind(interpTracksWind, id = sapply(strsplit(row.names(interpTracksWind), "[.]"), function(x) x[[1]]),
  t(apply(interpTracksWind[,c(8,9,10,13:20)], 1, function(x) {
    # x <- interpTracksWind[ 262,c(8,9,10,13:20)]
    if(!is.na(x[1]) & x[1]>0) {
      out <- c(RNCEP::NCEP.Tailwind(as.numeric(x[4]),  as.numeric(x[5]),  as.numeric(x[2]), as.numeric(x[3]))$fa,
               RNCEP::NCEP.Tailwind(as.numeric(x[6]),  as.numeric(x[7]),  as.numeric(x[2]), as.numeric(x[3]))$fa,
               RNCEP::NCEP.Tailwind(as.numeric(x[8]),  as.numeric(x[9]),  as.numeric(x[2]), as.numeric(x[3]))$fa,
               RNCEP::NCEP.Tailwind(as.numeric(x[10]), as.numeric(x[11]), as.numeric(x[2]), as.numeric(x[3]))$fa)
      t(data.frame(c(out, ifelse(any(!is.na(out)), max(out, na.rm = T), NA))))
    } else rep(NA, 5)
  })))
  # save(tracksMap, file = "/Tracks/tracksMap.rda")

# plot(land_map, xlim = c(160, 200))
# points(tracksMap$location.long, tracksMap$location.lat, pch = 16, 
#        col = ifelse(is.na(tracksMap$`5`), "grey90", "transparent"))

############
### loop ###
############
  
loopTracks <- cbind(tracksMap, leg = unlist(parallel::mclapply(as.numeric(tracksMap$timestamp), function(x) {
  x <- as.POSIXct(x, origin = "1970-01-01", tz = "GMT")
  leg <- c(1,2,3)[c(x>=as.POSIXct("2020-03-01", tz = "GMT") & x<=as.POSIXct("2020-04-14", tz = "GMT"),
                    x>=as.POSIXct("2020-04-25", tz = "GMT") & x<=as.POSIXct("2020-06-03", tz = "GMT"),
                    x>=as.POSIXct("2020-09-08", tz = "GMT") & x<=as.POSIXct("2020-10-04", tz = "GMT"))]
  if(length(leg)>0) leg else NA
}, mc.cores = 5)))
  
  
tmTab <- subset(loopTracks, !is.na(leg) & !duplicated(timestamp), select = c("timestamp", "leg", "file", "fileID"))
  
ERAindex <- cbind(tmTab$file, tmTab$fileID)[!duplicated(paste(tmTab$file, tmTab$fileID)),]
tmTab$TIFFindex <- apply(tmTab[,c("file", "fileID")], 1, function(x) which(x[1]==ERAindex[,1] & x[2]==ERAindex[,2]))


## col inds
clsF <- scales::colour_ramp(brewer.pal(11, "RdYlGn"))
cls  <- clsF(seq(0, 1, length = 100))
loopTracks$colID <- cut(loopTracks$'5', seq(-42, 42, length = length(cls)), labels = FALSE)


## global colors
blues  <- colorRampPalette(pal_seeblau)
greens <- colorRampPalette(pal_seegruen)

## anlge init
rot <- data.frame(long  = seq(130, 145, length = nrow(tmTab)),
                  fromx = seq(3,     0, length = nrow(tmTab)),
                  fromy = seq(7,     9, length = nrow(tmTab)),
                  fromz = seq(10,  8.5, length = nrow(tmTab)))


rot <- data.frame(long = rep(NA, nrow(tmTab)), fromx = NA, fromy = NA, fromz = NA)
rot[min(which(tmTab$leg==1)),] <- c(130, 3,  7, 10)
rot[max(which(tmTab$leg==1)),] <- c(145, 0,  9, 8.5)
rot[max(which(tmTab$leg==2)),] <- c(105, 0, 12, 6)
rot[max(which(tmTab$leg==3)),] <- c(115, 0,  8, 10)

rot <- apply(rot, 2, function(x) zoo::na.approx(x, rule = 3))

## wind track init
nrTracks <-  150

for(i in 1:nrow(tmTab)) {
  
  if(i==1 || tmTab$TIFFindex[i-1]!=tmTab$TIFFindex[i]) {

    tmpBrick <- brick(paste0("/bioing/user/slisovsk/MovSimData/geoTiffs/WindSnow_", tmTab$TIFFindex[i], ".tif"))
    spdR     <- sqrt(tmpBrick[[1]]^2 + tmpBrick[[2]]^2)

    top  <- spdR; top[] <- NA
    crds <- coordinates(top)[is.na(tmpBrick[[1]][]) & is.na(tmpBrick[[3]][]),]
    top[is.na(tmpBrick[[1]][]) & is.na(tmpBrick[[3]][])] <- raster::extract(topo, crds)
    top[coordinates(top)[,2]< -57] <- NA

    ocean  <- as.array(RGB(spdR, col = blues(100)))
    land   <- as.array(RGB(top,  col = greens(100)))

    ocean[land!=255]  <- land[land!=255]
    img <- abind(lapply(1:3, function(x) t(ocean[,dim(ocean)[2]:1,x]/255)), along = 3)

  }
  
  ## tracks ----
  new <- array(dim = c(10, 2, nrTracks))
  new[1,,] <- t(st_coordinates(suppressMessages(st_sample(ocean_map %>% st_difference(ocean_sample), size = dim(new)[3]))))

  if(i %in% which(!duplicated(tmTab$leg))) tracks <- new else tracks <- abind(tracks, new, along = 3)

  ind <- do.call("rbind", parallel::mclapply(1:dim(tracks)[3], function(x) {
    ind   <- max(which(!is.na(tracks[,1,x])))
    matrix(c(ind, tracks[ind,,x]), ncol = 3)
  }, mc.cores = parallel::detectCores()))

  wnd  <- extract(tmpBrick[[1:2]], ind[,2:3])
  dir  <- atan2(wnd[,1], wnd[,2]) * (180/pi)
  spd  <- sqrt(wnd[,1]^2 + wnd[,2]^2)

  dest <- geosphere::destPoint(ind[,2:3], dir, spd*60*60*5)

  del <- which(ind[,1]>=(dim(tracks)[1]-1) | is.na(dest[,1]))
  if(length(del)>0) {
    tracks <- tracks[,,-del]
    ind    <- ind[-del,]
    dest   <- dest[-del,]
  }

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

  
  
  ### Bar-tailed godwit tracks ----
  realTracks <- subset(loopTracks, timestamp<=tmTab[i,1] & leg==tmTab[i,2])
  last       <- lapply(unique(realTracks$id), function(x) realTracks[realTracks$id==x,][sum(realTracks$id==x),])
    lastTab  <- do.call("rbind", last)[!is.na(do.call("rbind", last)[,"location.long"]),]

  t_real <- st_wrap_dateline(st_sf(st_sfc(st_multipoint(as.matrix(lastTab[,c("location.long", "location.lat")]))), crs = 4326))

  grPoints <- dplyr::bind_rows(parallel::mclapply(split(as.data.frame(cbind(st_coordinates(t_real)[,1:2],
                                                                            colID = lastTab$colID, dist = lastTab$dist,
                                                                            support = lastTab$'5')),
                                                        1:nrow(st_coordinates(t_real))),
          function(x) {
              sp <- sph2car(matrix(c((x$X+rot[i,1])*-1, x$Y, 1), ncol = 3))[,c(1,3,2)]
              sphere(x = sp[1], y = 7+sp[2], z = sp[3], material = diffuse(color=ifelse(is.na(x$colID) | x$dist<5, "grey30", cls[x$colID])),
                     radius = ifelse(is.na(x$colID) | x$dist<5, 0.015, approx(c(0, 40), c(0.015, 0.045), abs(x$support))$y))
            }, mc.cores = 4))


  ls_lines <- lapply(unique(realTracks$id), function(x) {
    tmp <- subset(realTracks, id==x & !is.na(location.long), select = c("location.long", "location.lat"))
    if(nrow(tmp)>2)   st_linestring(as.matrix(tmp))
  })

  t_lines <- st_wrap_dateline(st_sf(st_sfc(st_multilinestring(ls_lines[!sapply(ls_lines, is.null)]), crs = 4326)))

  grLines <- dplyr::bind_rows(parallel::mclapply(split(as.data.frame(st_coordinates(t_lines)[,1:2]), st_coordinates(t_lines)[,3]), function(x) {
    x[,2] <- x[,2]*-1
    sph2car(cbind(x,1))[,c(1,3,2)] %>%
      path(x = 0, y = 7, material = diffuse(color="grey50"), width = 0.005, angle=c(180,rot[i,1],0))
  }, mc.cores = 4))

  
  
  # ### Plot ----
  if(i==1) {
      png(glue::glue("{tmp_path}/plot_color.png"), width = 300, height = 100, res = 100)
      opar1 <- par(mar = c(3,2,1,1), bty = "n", cex.axis = 0.5)

      plot(NA, ylim = c(0,1), xlim = c(-40, 40), yaxt = "n", xlab = "",
           ylab = "", mgp = c(0.1,0.5,0.1))

      plotrix::color.legend(-40, 0.1, 40, 0.9, legend = " ",
                            rect.col = clsF(seq(0,1,length = 40)), gradient="x")

      par(opar1)
      dev.off()
  }


  png(glue::glue("{tmp_path}/plot_winds_{i}.png"), width = 300, height = 500, res = 100)
  opar2 <- par(mfrow = c(3,1), mar = c(1,3,0.5,1), oma = c(3,0,0,0), bty = "n", cex.axis = 0.7)

  ## leg 1
  leg1 <- subset(loopTracks, timestamp<=tmTab[i,1] & leg==1 & !is.na(`5`) & dist>15)
  if(nrow(leg1)>10) {

    d <- density(leg1$`5`)

    plot(NA, ylim = c(0, max(d$y)), xlim = c(-40, 40), yaxt = "n", xlab = "", ylab = "", mgp = c(0.1,0.5,0.1))
    polygon(c(d$x, rev(d$x)), c(d$y, rep(0, length(d$y))), col = cls[cut(median(leg1$`5`, na.rm = T), seq(-42, 42, length = length(cls)), labels = FALSE)],
            border = NA)
    lines(d$x, d$y, col = "grey40")

    mtext("Leg 1", 3, line = -2, at = -35, cex = 0.7)

    if(tmTab[i,2]==1) {
      oparXPD <- par(xpd = NA)
      rect(-49, -0, -46, max(d$y), col = "grey80", border = NA)
      par(oparXPD)
    }

  } else {
    plot(NA, xlim = c(-40, 40), ylim = c(0,1), yaxt = "n", mgp = c(0.1,0.5,0.1), xlab = "", ylab = "")
    mtext("Leg 1", 3, line = -2, at = -35, cex = 0.7)
    if(tmTab[i,2]==1) {
      oparXPD <- par(xpd = NA)
      rect(-49, -0, -46, 1, col = "grey80", border = NA)
      par(oparXPD)
    }
  }

  leg2 <- subset(loopTracks, timestamp<=tmTab[i,1] & leg==2 & !is.na(`5`) & dist>15)
  if(nrow(leg2)>10) {

    d <- density(leg2$`5`)

    plot(NA, ylim = c(0, max(d$y)), xlim = c(-40, 40), yaxt = "n", xlab = "", ylab = "", mgp = c(0.1,0.5,0.1))
    polygon(c(d$x, rev(d$x)), c(d$y, rep(0, length(d$y))), col = cls[cut(median(leg2$`5`, na.rm = T), seq(-42, 42, length = length(cls)), labels = FALSE)],
            border = NA)
    lines(d$x, d$y, col = "grey40")

    mtext("Leg 2", 3, line = -2, at = -35, cex = 0.7)

    if(tmTab[i,2]==2) {
      oparXPD <- par(xpd = NA)
      rect(-49, -0, -46, max(d$y), col = "grey80", border = NA)
      par(oparXPD)
    }

  } else {
    plot(NA, xlim = c(-40, 40), ylim = c(0,1), yaxt = "n", mgp = c(0.1,0.5,0.1), xlab = "", ylab = "")
    mtext("Leg 2", 3, line = -2, at = -35, cex = 0.7)

    if(tmTab[i,2]==2) {
      oparXPD <- par(xpd = NA)
      rect(-49, -0, -46, 1, col = "grey80", border = NA)
      par(oparXPD)
    }
  }

  leg3 <- subset(loopTracks, timestamp<=tmTab[i,1] & leg==3 & !is.na(`5`) & dist>15)
  if(nrow(leg3)>10) {

    d <- density(leg3$`5`)

    plot(NA, ylim = c(0, max(d$y)), xlim = c(-40, 40), yaxt = "n", xlab = "", ylab = "", mgp = c(0.1,0.5,0.1))
    polygon(c(d$x, rev(d$x)), c(d$y, rep(0, length(d$y))), col = cls[cut(median(leg3$`5`, na.rm = T), seq(-42, 42, length = length(cls)), labels = FALSE)],
            border = NA)
    lines(d$x, d$y, col = "grey40")

    mtext("Leg 3", 3, line = -2, at = -35, cex = 0.7)
    mtext("Wind support (m/s)", 1, line = 2.5, cex = 1)

    if(tmTab[i,2]==3) {
      oparXPD <- par(xpd = NA)
      rect(-49, -0, -46, max(d$y), col = "grey80", border = NA)
      par(oparXPD)
    }

  } else {
    plot(NA, xlim = c(-40, 40), ylim = c(0,1), yaxt = "n", mgp = c(0.1,0.5,0.1), xlab = "", ylab = "")
    mtext("Leg 3", 3, line = -2, at = -35, cex = 0.7)
    mtext("Wind support (m/s)", 1, line = 2.5, cex = 1)

    if(tmTab[i,2]==3) {
      oparXPD <- par(xpd = NA)
      rect(-49, -0, -46, 1, col = "grey80", border = NA)
      par(oparXPD)
    }
  }

  par(opar2)
  dev.off()

  colScale <- magick::image_read(glue::glue("~/Desktop/tmp/plot_color.png"))
  plot1    <- magick::image_read(glue::glue("~/Desktop/tmp/plot_winds_{i}.png"))

  magick::image_append(c(magick::image_read(glue::glue("~/Desktop/tmp/plot_color.png")), 
                         magick::image_read(glue::glue("~/Desktop/tmp/plot_winds_{i}.png"))), stack = TRUE) %>% 
    magick::image_write(glue::glue("{tmp_path}/sidePlot_{i}.png"))
   
  #### render globe ----
  
  generate_studio(material=diffuse(color = "grey10")) %>%
    add_object(sphere(x = 0, y = 7, radius=0.9999, angle = c(0, rot[i,1], 0),
                      material = glossy(gloss=0.3, image_texture =  img))) %>%
    add_object(group_objects(grTracks)) %>%
    add_object(group_objects(grLines)) %>%
    add_object(group_objects(grPoints)) %>%
    add_object(sphere(y = 9, z = 10, x = 20, radius = 6, material=light(intensity=10))) %>%
    render_scene(width=550, height=600, aperture=0, fov=14, sample_method = "random", parallel = TRUE,
                 samples= 250, clamp_value=10, lookfrom=c(rot[i,2],rot[i,3],rot[i,4]), lookat=c(0,7,0), camera_up = c(0,1,0),
                 filename=glue::glue("{tmp_path}/globe_{i}.png"))
  
  title_mat = matrix(0,60,850) %>%
    add_title(title_text = glue::glue("{format(tmTab[i,1], '%Y-%m-%d')}        Bar-tailed godwit migration along the East Asian-Australasian Flyway"),
              title_bar_alpha = 1, title_bar_color = "grey30", title_size = 21,
              title_color = "white", filename=glue::glue("~/Desktop/tmp/title_{i}.png"))


  titel_image <- magick::image_read(glue::glue("{tmp_path}/title_{i}.png"))
  world_image <- magick::image_read(glue::glue("{tmp_path}/globe_{i}.png"))
  side_imge   <- magick::image_read(glue::glue("{tmp_path}/sidePlot_{i}.png"))
  foot_image  <- magick::image_read(glue::glue("MovSim/img/footer.png"))

    magick::image_append(c(
      titel_image,
        magick::image_append(
          c(world_image,
            side_imge), stack = FALSE),
        foot_image), stack = TRUE) %>% # magick::image_resize(magick::geometry_size_pixels(500, 500, preserve_aspect = F)) %>%
    magick::image_write(glue::glue("{tmp_path}/sim_image_{i}.png"))

}

av::av_encode_video(glue::glue("{tmp_path}/sim_image_{1:1690}.png"),
                    output = "Video/BTGviz.mov", framerate = 20)
 

