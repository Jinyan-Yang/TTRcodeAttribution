sm.2000.2010 <- raster::brick('//fs1-cbr.nexus.csiro.au/{mmrg}/work/users/yan190/repo/ERA5_monthly/sm/sm_2000_2010.nc')
sm.2010.2020 <- raster::brick('//fs1-cbr.nexus.csiro.au/{mmrg}/work/users/yan190/repo/ERA5_monthly/sm/sm_2010_2020.nc')
tair.2000.2010 <- raster::brick('//fs1-cbr.nexus.csiro.au/{mmrg}/work/users/yan190/repo/ERA5_monthly/tair/tair_2000_2010.nc')
tair.2010.2020 <- raster::brick('//fs1-cbr.nexus.csiro.au/{mmrg}/work/users/yan190/repo/ERA5_monthly/tair/tair_2010_2020.nc')

# library(rasterVis)
# levelplot(sm.2000.2010)
# levelplot(sm.2000.2010[[1]])
# names(sm.2000.2010)
# names(sm.2010.2020)

# plot(sm.2000.2010[[1]])
# points(x = ls.point.df$lon,ls.point.df$lat)
# crs(sm.2000.2010)
