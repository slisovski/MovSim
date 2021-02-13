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
topo <- raster("~/Google Drive/GeoDat/ETOPO1_Ice_g_geotiff.tif")

### Map
land_map     <- read_sf("~/Google Drive/GeoDat/NaturalEarth/50m_physical/ne_50m_land/ne_50m_land.shp") %>% st_geometry()
ocean_map    <- st_sym_difference(st_as_sfc(st_bbox(c(xmin = -180, xmax = 180, ymax = -90, ymin = 90), crs = st_crs(4326))), st_union(land_map))
ocean_sample <- st_sfc(st_polygon(list(matrix(c(-80, -80, 89, 89, -80,-90, 90, 90, -90,-90), ncol = 2))), crs = 4326)

### tracks
load("~/Google Drive/Science/ProjectsData/MovSim/Tracks/interpTracksWind.rda")
# head(interpTracksWind)

  ### wind support
  tracksMap <- cbind(interpTracksWind, t(apply(interpTracksWind[,c(8,9,10,13:20)], 1, function(x) {
    # x <- interpTracksWind[ 262,c(8,9,10,13:20)]
    if(!is.na(x[1]) & x[1]>0) {
      out <- c(RNCEP::NCEP.Tailwind(as.numeric(x[4]),  as.numeric(x[5]),  as.numeric(x[2]), as.numeric(x[3]))$fa,
               RNCEP::NCEP.Tailwind(as.numeric(x[6]),  as.numeric(x[7]),  as.numeric(x[2]), as.numeric(x[3]))$fa,
               RNCEP::NCEP.Tailwind(as.numeric(x[8]),  as.numeric(x[9]),  as.numeric(x[2]), as.numeric(x[3]))$fa,
               RNCEP::NCEP.Tailwind(as.numeric(x[10]), as.numeric(x[11]), as.numeric(x[2]), as.numeric(x[3]))$fa)
      t(data.frame(c(out, max(out, na.rm = T))))
    } else rep(NA, 5)
  })))
  tracksMap$id <- sapply(strsplit(row.names(tracksMap), "[.]"), function(x) x[[1]])


# plot(land_map, xlim = c(160, 200))
# points(tracksMap$location.long, tracksMap$location.lat, pch = 16, col = ifelse(tracksMap$dist>0, "grey90", "orange"))

############
### loop ###
############
  
loopTracks <- cbind(tracksMap, leg = unlist(parallel::mclapply(as.numeric(tracksMap$timestamp), function(x) {
  x <- as.POSIXct(x, origin = "1970-01-01", tz = "GMT")
  leg <- c(1,2,3)[c(x>=as.POSIXct("2020-03-01", tz = "GMT") & x<=as.POSIXct("2020-04-14", tz = "GMT"),
                    x>=as.POSIXct("2020-04-25", tz = "GMT") & x<=as.POSIXct("2020-06-05", tz = "GMT"),
                    x>=as.POSIXct("2020-09-08", tz = "GMT") & x<=as.POSIXct("2020-10-04", tz = "GMT"))]
  if(length(leg)>0) leg else NA
}, mc.cores = 5)))
  
  
tmTab <- subset(loopTracks, !is.na(leg) & !duplicated(timestamp), select = c("timestamp", "leg", "file", "fileID"))
  
ERAindex <- cbind(loopTab$file, loopTab$fileID)[!duplicated(paste(loopTab$file, loopTab$fileID)),]
tmTab$TIFFindex <- apply(tmTab[,c("file", "fileID")], 1, function(x) which(x[1]==ERAindex[,1] & x[2]==ERAindex[,2]))


## col inds
cls <- brewer.pal(11, "RdYlGn")
loopTracks$colID <- cut(loopTab$'5', seq(-40, 40, length = length(cls)), labels = FALSE)


## global colors
blues  <- colorRampPalette(pal_seeblau)
greens <- colorRampPalette(pal_seegruen)

## anlge init
rot <- data.frame(long  = seq(130, 145, length = nrow(TIFFindex)),
                  fromx = seq(3,     0, length = nrow(TIFFindex)),
                  fromy = seq(7,     9, length = nrow(TIFFindex)),
                  fromz = seq(10,  8.5, length = nrow(TIFFindex)))


rot <- data.frame(long = rep(NA, nrow(tmTab)), fromx = NA, fromy = NA, fromz = NA)
rot[min(which(tmTab$leg==1)),] <- c(130, 3,  7, 10)
rot[max(which(tmTab$leg==1)),] <- c(145, 0,  9, 8.5)
rot[max(which(tmTab$leg==2)),] <- c(105, 0, 12, 6)
rot[max(which(tmTab$leg==3)),] <- c(115, 0,  8, 10)

rot <- apply(rot, 2, function(x) zoo::na.approx(x, rule = 3))

## wind track init
nrTracks <-  75

for(i in 1:nrow(tmTab)) {
  
  if(i==1 || tmTab$TIFFindex[i-1]!=tmTab$TIFFindex[i]) {
    
    tmpBrick <- brick(paste0("~/Google Drive/Science/ProjectsData/MovSim/geoTiffs/WindSnow_", tmTab$TIFFindex[i], ".tif"))
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
  
  if(i==1) tracks <- new else tracks <- abind(tracks, new, along = 3)

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
  realTracks <- subset(loopTracks, timestamp<=tmTab[i,1])
  last       <- lapply(unique(realTracks$id), function(x) realTracks[realTracks$id==x,][sum(realTracks$id==x),])

  t_real <- st_wrap_dateline(st_sf(st_sfc(st_multipoint(as.matrix(do.call("rbind", last)[!is.na(do.call("rbind", last)[,"location.long"]),c("location.long", "location.lat")]))), crs = 4326))
   
  grPoints <- dplyr::bind_rows(parallel::mclapply(split(as.data.frame(cbind(st_coordinates(t_real)[,1:2], id = 1:nrow(st_coordinates(t_real)), colID = do.call("rbind", last)$colID, dist = do.call("rbind", last)$dist)),
                                                        1:nrow(st_coordinates(t_real))),
          function(x) {
              sp <- sph2car(matrix(c((x$X+rot[i,1])*-1, x$Y, 1), ncol = 3))[,c(1,3,2)]
              sphere(x = sp[1], y = 7+sp[2], z = sp[3], material = diffuse(color=ifelse(is.na(x$colID) | x$dist<10, "grey30", cls[x$colID])), radius = 0.02)
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
  
  
  
  # ### Plot1 ----
  # png(glue::glue("~/Desktop/tmp/plot1_{i}.png"), width = 700, height = 450, res = 100)
  # opar1 <- par(mfrow = c(1,3), mar = c(3,1,1,1), oma = c(1,8,1,0), bty = "n")
  # 
  # plot(NA, ylim = c(-8,length(unique(tracksMap$id))), xlim = c(-40, 40), yaxt = "n", xlab = "",
  #      ylab = "")
  # mtext("Leg 1", 3, cex = 0.8)
  # 
  # out1 <- do.call("rbind", lapply(unique(tracksMap$id), function(x) {
  #   tmp <- subset(realTracks, id == x)
  #   if(nrow(tmp)>0) {
  #    points(tmp$`5`, jitter(rep(which(tmp$id[1]==unique(tracksMap$id)), nrow(tmp)), 0.2), pch = 16,
  #           col = adjustcolor(cls[tmp$colID], alpha.f = 0.6))
  #    data.frame(s = tmp$`5`, id = which(tmp$id[1]==unique(tracksMap$id)), col = cls[tmp$colID], type = c(rep(1, nrow(tmp)-1),2))[!is.na(tmp$`5`) & !is.infinite(abs(tmp$`5`)),]
  #   }
  # }))
  # 
  # with(subset(out1, type == 2), points(s, id, pch = 21, bg = col,  col = "grey30", lwd = 1.6))
  # 
  # oparXPD <- par(xpd = NA, mar = par()$mar + c(2.5, 0, 1, 0))
  # rect(-44, -8.5, 44, 16, border = "grey70")
  # par(oparXPD)
  # 
  # if(nrow(out1)>10) {
  # d <- density(out1$s)
  # 
  # par(new = TRUE)
  # plot(NA, ylim = c(0,0.35), xlim = c(-40, 40), yaxt = "n", xlab = "", ylab = "")
  # polygon(c(d$x, rev(d$x)), c(d$y, rep(0, length(d$y))), col = "burlywood", border = NA)
  # lines(d$x, d$y, col = "grey40")
  # }
  # 
  # plot(NA, ylim = c(1,length(unique(tracksMap$id))), xlim = c(-40, 40), yaxt = "n", 
  #      xlab = "",  ylab = "")
  # mtext("Wind support", 1, cex = 1.2, line = 3)
  # mtext("Leg 2", 3, cex = 0.8)
  # 
  # plot(NA, ylim = c(1,length(unique(tracksMap$id))), xlim = c(-40, 40), yaxt = "n", 
  #      xlab = "",  ylab = "")
  # mtext("Leg 3", 3, cex = 0.8)
  # 
  # 
  # par(opar1)
  # dev.off()
  # 
  # ### Plot2 ----
  # png(glue::glue("~/Desktop/tmp/plot2_{i}.png"), width = 700, height = 450, res = 100)
  # opar2 <- par(mar = c(3,6,1,1))
  # 
  # plot(NA, xlim = range(datesTab$date), ylim = c(0, 3200), las = 1, ylab = "", xlab = "", xaxt = "n", 
  #      yaxt = "n")
  # axis(2, at = c(10, 779, 1502, 3130), las = 1)
  # mtext("Altitude with max. support (a.s.l.)", 2, line = 4, cex = 1.2)
  # 
  # brks <- as.numeric(datesTab$date)[suppressWarnings(which(datesTab$leg!=datesTab$leg[-1]))][1:2]
  # abline(v = brks, lty = 3, col = "grey80")
  # 
  # axis(1, at = as.numeric(datesTab$date[c(TRUE, rep(FALSE, 10))]), labels = as.Date(datesTab$indDate)[c(TRUE, rep(FALSE, 10))])
  # plotrix::axis.break(1, brks[1], style="slash")
  # plotrix::axis.break(1, brks[2], style="slash")
  # plotrix::axis.break(3, brks[1], style="slash")
  # plotrix::axis.break(3, brks[2], style="slash")
  # 
  # 
  # datAlt <- data.frame(t = as.numeric(realTracks[!is.na(realTracks$`5`) & realTracks$dist > 25, "timestamp"]),
  #                      a = unlist(apply(realTracks[!is.na(realTracks$`5`) & realTracks$dist > 25, c('1', '4', '3', '2')], 1, function(x) {
  #                           if(any(!is.na(x) | !is.nan(x))) c(10, 779, 1502, 3130)[which.max(x)] else NA })))
  # 
  # points(datAlt$t, jitter(datAlt$a, 0.2), pch = 16, col = adjustcolor("grey40", alpha.f = 0.4))
  # 
  # if(nrow(datAlt)>10) {
  # mod <- mgcv::gam(a~s(t), data = datAlt)
  # fitt <- data.frame(t = seq(min(datAlt$t), max(datAlt$t), by = 12*60*60))
  # fit  <- predict(mod, newdata = fitt, type="response", se=T)$fit
  # se   <- predict(mod, newdata=fitt, type="response", se=T)$se.fit
  # 
  # polygon(c(fitt$t, rev(fitt$t)), c(fit-1.96*se, rev(fit+1.96*se)), border = NA, col = adjustcolor("grey50", alpha.f = 0.5))
  # lines(fitt$t, fit, lwd = 2, col = "orange")
  # }
  # par(opar2)
  # dev.off()
  # 
  
  
  #### render ----
  
  generate_studio(material=diffuse(color = "grey10")) %>%
    add_object(sphere(x = 0, y = 7, radius=0.9999, angle = c(0, rot[i,1], 0),
                      material = glossy(gloss=0.3, image_texture =  img))) %>%
    add_object(group_objects(grTracks)) %>%
    add_object(group_objects(grLines)) %>%
    add_object(group_objects(grPoints)) %>%
    add_object(sphere(y = 9, z = 10, x = 20, radius = 6, material=light(intensity=10))) %>%
    render_scene(width=400, height=400, aperture=0, fov=14, sample_method = "random", parallel = TRUE,
                 samples= 200, clamp_value=10, lookfrom=c(rot[i,2],rot[i,3],rot[i,4]), lookat=c(0,7,0), camera_up = c(0,1,0), filename=glue::glue("~/Desktop/tmp/globe_{i}.png"))
  
  
  # title_mat = matrix(0,55,1400) %>%
  #   add_title(title_text = glue::glue("Bar-tailed godwit migration         {format(TIFFindex[i,1], '%Y-%m-%d')}"), 
  #             title_bar_alpha = 1, title_bar_color = "grey30", title_size = 25,
  #             title_color = "white", filename=glue::glue("~/Desktop/tmp/title_{i}.png"))
  # 
  # # title_mat = matrix(1, 45, 1400) %>%
  # #   add_title(title_text = glue::glue("by: Simeon Lisovski; traking data from: Jesse Conklin, Phil Battley; wind data: ECWMF ERA5; Land topography: ETOPO2; Snow Cover: IMS.                            Rcode = {'https://github.com/slisovski/MovSim'}"), 
  # #             title_bar_alpha = 1, title_bar_color = "grey30", title_size = 10,
  # #             title_color = "white", filename=glue::glue("~/Desktop/tmp/footer_{i}.png"))
  # 
  # titel_image <- magick::image_read(glue::glue("~/Desktop/tmp/title_{i}.png"))
  # world_image <- magick::image_read(glue::glue("~/Desktop/tmp/globe_{i}.png"))
  # plt1_imge   <- magick::image_read(glue::glue("~/Desktop/tmp/plot1_{i}.png"))
  # plt2_imge   <- magick::image_read(glue::glue("~/Desktop/tmp/plot2_{i}.png"))
  # foot_image  <- magick::image_read(glue::glue("~/Desktop/tmp/footer.png"))
  # 
  # magick::image_append(c(titel_image,
  #     magick::image_append(
  #       c(world_image,
  #        magick::image_append(c(plt1_imge, plt2_imge), stack = TRUE)),
  #     ),
  #     foot_image), stack = TRUE) %>% 
  #   magick::image_write(glue::glue("~/Desktop/tmp/full_image_{i}.png"))
  # 
  # suppressWarnings({
  # file.remove(c(glue::glue("~/Desktop/tmp/title_{i}.png"),
  #               glue::glue("~/Desktop/tmp/globe_{i}.png"), 
  #               glue::glue("~/Desktop/tmp/plot1_{i}.png"),
  #               glue::glue("~/Desktop/tmp/plot2_{i}.png"),
  #               glue::glue("~/Desktop/tmp/footer_{i}.png")
  #               ))})
}

# av::av_encode_video(glue::glue("~/Google Drive//tmp/full_image_{1:121}.png"), 
                    # output = "~/Google Drive//globe_viz.mp4", framerate = 15)




