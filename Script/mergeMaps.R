library(raster)
library(unikn)
library(abind)
library(dichromat)
library(rayrender)

### global elevation model
topo <- raster("/Users/slisovsk/Google Drive/GeoDat/ETOPO1_Ice_g_geotiff.tif")



### example raster
tmpBrick <- brick("~/Google Drive/Science/ProjectsData/MovSim/geoTiffs/WindSnow_1.tif")



## spd
spdR <- sqrt(tmpBrick[[1]]^2 + tmpBrick[[2]]^2)

top  <- spdR; top[] <- NA
crds <- coordinates(top)[is.na(tmpBrick[[1]][]) & is.na(tmpBrick[[3]][]),]
top[is.na(tmpBrick[[1]][]) & is.na(tmpBrick[[3]][])] <- raster::extract(topo, crds)
top[coordinates(top)[,2]< -57] <- NA

blues  <- colorRampPalette(pal_seeblau)
greens <- colorRampPalette(pal_seegruen)

ocean  <- as.array(RGB(spdR, col = blues(100)))
land   <- as.array(RGB(top,  col = greens(100)))

ocean[land!=255]  <- land[land!=255]

img <- abind(lapply(1:3, function(x) t(ocean[,dim(ocean)[2]:1,x]/255)), along = 3)


generate_studio(material=dielectric(color = "white")) %>%
  add_object(sphere(x = 0, y = 1, radius=0.9999, angle=c(0,0,0),
                    material = glossy(gloss=0.3, image_texture =  img))) %>%
  add_object(sphere(y = 9, z = 10, x = 20, radius = 6, material=light(intensity=10))) %>%
  render_scene(width=200, height=200, aperture = 0, fov = 15, sample_method = "stratified", parallel = TRUE,
               samples= 200, clamp_value=10, lookfrom=c(0,1,10), camera_up = c(0,1,0))
