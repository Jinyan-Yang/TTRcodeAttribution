library(raster)
library(sf)
# df.biome <- readRDS('//fs1-cbr.nexus.csiro.au/{mmrg}/work/users/yan190/repo/delta_n_15/cache/ls.d15n.slope.global.rds')
# df.biome.au <- df.biome[df.biome$lon > 110 & 
#                           df.biome$lon <155&
#                           df.biome$lat > -45 & 
#                           df.biome$lat < -10,]
# sample havplot
hav.df <- read.csv('//fs1-cbr.nexus.csiro.au/{lw-bio-digi}/work/nativeness/plot_data/HAVPlot_nativeness_2021-02-05.csv')
hav.df.good <- hav.df[hav.df$proportionNativeCover > 0.8,]
hav.df.good$lon <- hav.df.good$decimalLongitude
hav.df.good$lat <- hav.df.good$decimalLatitude

kpn.ra <- raster('//fs1-cbr.nexus.csiro.au/{mmrg}/work/users/yan190/repo/dynamics_HC/data/kpn/kpnall.txt')

hav.df.good$kpn <- extract(kpn.ra,cbind(hav.df.good$lon,
                                        hav.df.good$lat))

get.sample.func <- function(dat.in,sample.size = 800){
  
  dat.in$x <- dat.in$lon
  dat.in$y <- dat.in$lat
  
  if(sample.size >= nrow(dat.in)){
    coord.bad <- dat.in
    coord.bad$x <- coord.bad$lon
    coord.bad$y <- coord.bad$lat
    coord.bad <- coord.bad[,c('x','y')]
  }else{
    st.obj.bad <- st_multipoint(as.matrix(dat.in[,c('x','y')]))#cbind(bad.df[,c('x')],bad.df[,c('y')])
    
    smp.bad <- st_sample(st.obj.bad,size = sample.size,type='regular')
    
    coord.bad <- data.frame(x = smp.bad[[1]][,1],y = smp.bad[[1]][,2])
  }
  return(coord.bad)
}
# #####
sample.func <- function(hav.df.good){
  # hav.df.good$lon <- hav.df.good$decimalLongitude
  # hav.df.good$lat <- hav.df.good$decimalLatitude
  # hav.df.good <- get.kpn.func(hav.df.good)
  # 
  hav.df.good.sub <- hav.df.good[!is.na(hav.df.good$kpn),]
  
  
  kpm.vec <- unique(hav.df.good.sub$kpn)
  
  hav.df.good.sub$x <- hav.df.good.sub$lon
  hav.df.good.sub$y <- hav.df.good.sub$lat

  hav.good.ls <- lapply(kpm.vec, function(k.nm){
    
    tmp.df <- get.sample.func(hav.df.good.sub[hav.df.good.sub$kpn == k.nm,],
                              sample.size = 100)

    tmp.df$kpn <- k.nm
    return(tmp.df)
  })
  names(hav.good.ls) <- kpm.vec
  coord.hav.good <-  do.call(rbind, Map(cbind, 
                                        Name = names(hav.good.ls), 
                                        hav.good.ls))
  return(coord.hav.good)
}
#########
sample.df <- sample.func(hav.df.good)
sample.df.sub <- sample.df[,c("x",'y')]
names(sample.df.sub) <- c('lon','lat')
xxxxx <- hav.df.good[!is.na(hav.df.good$decimalLongitude),]
xxxxx[duplicated(xxxxx),]
xxxxx[xxxxx$decimalLongitude == 149.3349,]
df.biome.sub <- hav.df.good[!is.na(hav.df.good$kpn),]

# both.df <- rbind(sample.df.sub,df.biome.sub)
both.df <- dplyr::left_join(sample.df.sub,df.biome.sub)
both.df$yr <- lubridate::year(as.Date(both.df$obsStartDate))
both.df <- both.df[both.df$yr > 1983,]
both.df <- both.df[order(both.df$yr,decreasing = T),]
both.df <- both.df[!duplicated(both.df[,c("lon",'lat')]),]

# both.df <- split(both.df,both.df$plotID)
# sample.index <- which(duplicated(both.df))
# sample.index <- sample.index - nrow(sample.df.sub)
# df.biome.sub[sample.index[1],]
# saveRDS(sample.index,'cache/auSites.rds')
write.csv(both.df,'cache/auSitesGPS_havPlot.csv',row.names=F)


xxx.df <- read.csv('cache/auSitesGPS_havPlot.csv')
xxx.df$lon
