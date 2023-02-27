

##### GEOLOGICAL AND ENVIRONMENTAL FACTORS SHAPE REEF #####
##### HABITAT STRUCTURE IN THE MAIN HAWAIIAN ISLANDS #####
################### ASBURY ET AL. 2023 ###################

library(rgeos)
library(sp)
source("R/functions.R")

# DEMs & DEM extent 
files <- list.files("data")

# allocate space
complexity <- data.frame()

for (f in seq(1, length(files), 2)){
  
  borderfilename = files[f]
  infilename= files[f+1]
  
  # Load rasters (DEM & bounding extent)
  temp.ras <- raster(paste0("data/", infilename))
  data.b <- raster(paste0("data/", borderfilename))
  plot(data)
  
  # Bounding box
  a <- aggregate(data.b, fact = 100)
  r <- a > -Inf
  pp <- rasterToPolygons(r, dissolve = TRUE)
  buff <- gBuffer(pp, width = -0.5)
  ex <- extent(buff)
  plot(buff, add = TRUE)
  
  # Scope (extent), scales of variation, and resolution (grain)
  L <- 2 #2x2m box
  scl <- L / c(1, 2, 3.125, 6.25, 12.5, 25, 50, 100, 200) #2m-1cm
  L0 <- min(scl) 
  resolution <- 0.01 
  rep <- 1
  i <- 0
  j <- 0
  site = strsplit(names(data), split = "_") 

  # Delineate 5 boxes within a single plot
  while (i < 5) {
    
    j <- j + 1
    x0 <- runif(1, ex[1], ex[2])
    y0 <- runif(1, ex[3], ex[4])
    rect(x0, y0, x0 + L, y0 + L, border="red")
    
    data <- crop(temp.ras, extent(x0, x0 + L, y0, y0 + L))
    
    # Ensure box is fully inside plot
    x.crd <- c(x0,x0,x0,x0+L/2,x0+L,x0+L,x0+L,x0+L/2)
    y.crd <- c(y0,y0+L/2,y0+L,y0+L,y0+L,y0+L/2,y0,y0)
    
    overlap <- point.in.polygon(x.crd, y.crd, buff@polygons[[1]]@Polygons[[1]]@coords[,1], buff@polygons[[1]]@Polygons[[1]]@coords[,2])
    
    if(sum(overlap) == 8){
      rect(x0, y0, x0 + L, y0 + L, border="blue")
      
      # Correct slopes
      x <- seq(data@extent[1]+res(data)[1]/2, data@extent[2]-res(data)[1]/2, res(data)[1])
      y <- seq(data@extent[3]+res(data)[2]/2, data@extent[4]-res(data)[2]/2, res(data)[2])
      z <- data@data@values
      xy <- expand.grid(y, x)
      
      mod <- lm(z ~ xy$Var1 + xy$Var2)
      
      data@data@values <- mod$residuals
      
      #plot(data)
      #writeRaster(data, (paste0("output/rasters/", site, "_", rep, "_fixed.tif"))) 
      
      temp <- crop(data, extent(x0, x0 + L, y0, y0 + L))
      i <- i + 1
       
      # Ensure no NAs in plot
      na.fix <- table(is.na(temp[])) }}
      if(is.na(na.fix[2])){
      
        # Extract habitat metrics
        start.time <- Sys.time()
        example <- height_variation(write=TRUE, return=TRUE)
        out <- rdh(example)
  
        rep <- rep + 1
        i <- i + 1
        
        complexity <- rbind(complexity, data.frame(site = site[[1]][1],
                           box = i,
                           x = x0,
                           y = y0,
                           rugosity = out$R,
                           rugosity_theory = out$R_theory,
                           fractal = out$D,
                           fractal_theory = out$D_theory,
                           height = out$H))
        
        end.time <- Sys.time()
        print(end.time - start.time)

      }
  
  # Start on next plot if 5 boxes cannot fit 
  if (j > 200000){
    i = 5
    print(paste0(site[[1]][1], "  too thin"))
  }
}
  
write.csv(complexity, "output/MHI2019_Complexity.csv")
 



