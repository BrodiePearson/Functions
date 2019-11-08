# ### Dependencies
# 
# # converting to spatial data
# library(spatstat)
# library(maptools)
# library(sf)
# library(rgdal)
# 
# 
# ### Data
# path <- '~/Google Drive (jenna_pearson@brown.edu)/Baylor & Jenna/Projects/Eulerian-Lagrangian SFs (Obs)/LASER/Data/Structure Functions/'
# filename <-'LASER_All_SPOT_Drifters_latlon.txt'
# 
# 
# df <- read.delim(paste(path,filename, sep=""), header=FALSE)
# colnames(df)<- c("time", "lon","lat")

beg <-100
fin <-9643
# loop through timesteps
for(tt in beg:fin) {
  print(tt)
  #subset data for timestep
  dftt <- df[df$time == tt, ]

  #remove nans
  dftt <- na.omit(dftt)
  
  #remove time id column
  dftt$time <- NULL

    
  #convert lat lon to local coordinates
  coordinates(dftt) <- c("lon", "lat")
  proj4string(dftt) <- CRS("+init=epsg:4326") # WGS 84
  CRS.new <- CRS("+init=epsg:4489")
  dftt <- spTransform(dftt, CRS.new)
  
  # Convert the dataframe to a spatial object. Note that the
  # crs= 4326 parameter assigns a WGS84 coordinate system to the
  # spatial object
  p.sf <- st_as_sf(dftt, coords = c("lon", "lat"), crs = 4326)

  # Convert to ppp object
  p.sp  <- as(p.sf, "Spatial")  # Create Spatial* object
  p.ppp <- as(p.sp, "ppp")      # Create ppp object
  
  # par(mfrow=c(1,2))
  # plot(dftt, axes=TRUE, main="Original lat-lon", cex.axis=.95)
  # plot(dftt.ch1903, axes=TRUE, main="Projected", cex.axis=.95)
  # unclass(dftt.ch1903)

  
  #K <- Kest(p.ppp, correction="Ripley", rmax = '10')

  L <- Lest(p.ppp, correction="Ripley")
  #Correct to set to zero
  L$iso  <- L$iso  - L$r
  L$theo <- L$theo - L$r
  p <- 0.05
  n <- 200
  # E <- envelope(p.ppp, Lest, nsim = n,rank=(p * (n + 1)))
  # E$hi  <- E$hi  - E$r
  # E$lo <- E$lo - E$r

  #--------------------------------------------------
  # Plot Results
  #--------------------------------------------------
  
  plot(L$r/1000, L$iso,type = 'l', col = 'cornflowerblue',xlab="r (km)", 
       ylab="L(r)", main = 'All LASER Drifters at All times', cex = 2, ylim = c(-2000,150000),
       xlim = c(0,500), cex.lab=1.25, cex.axis=1.25, cex.main=1.25, cex.sub=1.5, lwd = 1)
  # polygon(c(E$r/1000,rev(E$r/1000)),c(E$hi ,rev(E$lo)),col="thistle",border=NA)
  # lines(L$r/1000, L$theo,type = 'l', col = 'black', lwd = 4)
  # 
  # legend(50,100000,legend=c("LASER Drifters", "Simulation Envelope"), col=c("cornflowerblue", "thistle"), 
  #        lty=c(1,1), lwd = 4, cex=1, bty = "n")
  par(new = 'T')
}