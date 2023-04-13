library(raster)
df.biome <- readRDS('//fs1-cbr.nexus.csiro.au/{mmrg}/work/users/yan190/repo/delta_n_15/cache/ls.d15n.slope.global.rds')
df.biome.au <- df.biome[df.biome$lon > 110 & 
                          df.biome$lon <155&
                          df.biome$lat > -45 & 
                          df.biome$lat < -10,]
kpn.ra <- raster('//fs1-cbr.nexus.csiro.au/{mmrg}/work/users/yan190/repo/dynamics_HC/data/kpn/kpnall.txt')

df.biome.au$kpn <- extract(kpn.ra,cbind(df.biome.au$lon,df.biome.au$lat))

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
sample.df <- sample.func(df.biome)
sample.df.sub <- sample.df[,c("x",'y')]
names(sample.df.sub) <- c('lon','lat')
df.biome.sub <- df.biome[,c('lon','lat')]

both.df <- rbind(sample.df.sub,df.biome.sub)
sample.index <- which(duplicated(both.df))
sample.index <- sample.index - nrow(sample.df.sub)
# df.biome.sub[sample.index[1],]
saveRDS(sample.index,'cache/auSites.rds')
write.csv(df.biome.sub[sample.index,],'cache/auSitesGPS.csv',row.names=F)
