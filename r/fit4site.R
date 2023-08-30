library(DEoptim)
library(foreach)
library(doParallel)
library(plante)
library(lubridate)
library(doBy)
library(randomForest)
library(raster)
library(parallel)
library(coda)
# 
source('r/readInputs.R')
source('r/function_abc.R')
source('r/fitDeoptim.R')
source('r/readERA5.R')
source('r/function_fitting.R')
# 

# read data
landsat.ts.ls.sample <- readRDS('cache/sampleAU_havPlot.ls.rds')
landsat.ts.ls.sample[[1]]


#fitting even sample 0.1 degree####
s.time <- Sys.time()
out.ls <- lapply(landsat.ts.ls.sample, fit.func)
e.time <- Sys.time() - s.time

saveRDS(out.ls,'cache/fittedParAuSites.rds')

#fitting havplot####

s.time <- Sys.time()
out.ls.hav <- lapply(landsat.n15.ls, fit.func)
e.time <- Sys.time() - s.time

# ls.point.df <- landsat.n15.ls[[1]]
saveRDS(out.ls.hav,'cache/fittedParHavPlot.rds')

out.hav.median.ls <- lapply(out.ls.hav, function(df){
  # print(head(df))
  if(!is.na(df$lon[1])){
    require(doBy)
    # print(df)
    df.sum <- summaryBy(.~lon+lat,
                        data = df,
                        FUN=c(median),
                        na.rm=T,keep.names = T)
    
    return(df.sum)
  }else{
    return(data.frame(lon = NA,lat = NA))
  }
  
})

out.hav.median.df <- dplyr::bind_rows(out.hav.median.ls)

# analysis#####
out.ls <- readRDS('cache/fittedParHavPlot.rds')
# out.ls <- lapply(out.ls, function(x) if(is.null(x)) data.frame(lon = NA) else x)
# 

out.median.ls <- lapply(out.ls, function(df){
  # print(head(df))
  df <- df[complete.cases(df),]
  if(!is.null(nrow(df))){

      # require(doBy)
      # df.sum <- summaryBy(.~lon+lat,
      #                     data = df,
      #                     FUN=c(median),
      #                     na.rm=T,keep.names = T)
      df.sum <- df[1,]
    
    return(df.sum)
  }else{
    return(data.frame(lon = NA,lat=NA))
  }
 
})

out.median.df <- dplyr::bind_rows(out.median.ls)
out.median.df$site.no <- row.names(out.median.df)

out.median.df <- out.median.df[order(out.median.df$nrmse.ndvi),]

good.df <- out.median.df[1:10,]
good.index <- as.numeric(good.df$site.no)

good.input.ls <- landsat.ts.ls.sample[good.index]
# good.input.ls[[6]]

for (i in seq_along(good.input.ls)){
  print(paste0(i,': ',nrow(good.input.ls[[i]])))
}

good.input.ls <- good.input.ls[c(5,7,8,9,10)]
# get met data
input.val.ls <- lapply(good.input.ls, function(dat){
  
  dat$lon <- as.numeric(dat$lon)
  dat$lat <- as.numeric(dat$lat)
  met.df <- data.frame(sm = extract(sm.ra,
                                     cbind(dat$lon[1],
                                           dat$lat[1]))[1,],
                       tair = extract(tair.ra,
                                       cbind(dat$lon[1],
                                             dat$lat[1]))[1,],
                       date = names(sm.ra))
  
  met.df$date <- as.Date(met.df$date,'X%Y.%m.%d')
  met.df$yr <- year(met.df$date)
  met.df$mon <- month(met.df$date)
  dat$yr <- year(dat$date)
  dat$mon <- month(dat$date)
  
  out.df <- merge(dat,
                  met.df[,which(colnames(met.df) != 'date')],
                  by=c('yr','mon'),all=T)
  
  out.df <- out.df[order(out.df$date),]
  out.df <- out.df[!is.na(out.df$sm),]
  out.df$tair <- out.df$tair - 272.15
  out.df$tsoil <- out.df$tair
  
  out.df <- out.df[order(out.df$yr,out.df$mon),]
  return(out.df)
})
# xxx <- input.val.ls[[1]]
# get par value
good.df <- good.df[c(5,7,8,9,10),]
# put together
save.ls <- list(fittedPar = good.df,
                inputList = input.val.ls)
# xxx <-save.ls$inputList[[4]]
saveRDS(save.ls,'cache/data4testing.rds')
save.ls <- readRDS('cache/data4testing.rds')

# plot##########
library(geodata)
clim.ra <- worldclim_country("Australia", var="bio", res = 10,path='wc10/')

out.median.df$map <- extract(clim.ra[[12]],cbind(out.median.df$lon,
                                                 out.median.df$lat))[,1]
out.median.df$mat <- extract(clim.ra[[1]],cbind(out.median.df$lon,
                                                out.median.df$lat))[,1]

out.median.df.narm <- out.median.df[complete.cases(out.median.df),]
out.median.df.narm$nrmse.d15n[out.median.df.narm$nrmse.d15n>1] <- 1

out.median.df$col.d15n <- cut(out.median.df$nrmse.d15n,
                              breaks = c(seq(0,1,by=0.2),2.5))
# range(out.median.df$nrmse.d15n)
palette(c(hcl.colors(5)[1:4],'coral','red'))
plot(mat~log(map,base = 10),
     data = out.median.df[order(out.median.df$nrmse.d15n),],
     pch=16,
     col = col.d15n,cex = 1)
legend('bottomleft',legend = levels(out.median.df$col.d15n),
       col = palette(),pch=16)


out.median.df$col.ndvi <- cut(out.median.df$nrmse.ndvi,
                              breaks = c(seq(0,1,by=0.2),2.5))

out.median.df$col.dn15 <- cut(out.median.df$nrmse.d15n,
                              breaks = c(seq(0,1,by=0.2),2.5))


out.median.df$col.k <- cut(out.median.df$k.slope,
                           breaks = c(seq(0,100,by=20),500))


hist(out.median.df$k.slope)
# range(out.median.df$k.slope)
palette(c(hcl.colors(5)[1:4],'coral','red'))
plot(mat~log(map,base = 10),
     data = out.median.df[order(out.median.df$nrmse.ndvi),],
     pch=16,
     col = col.d15n,cex = 1)
legend('bottomleft',legend = levels(out.median.df$col.ndvi),
       col = palette(),pch=16)
# 
plot(mat~map,
     data = out.median.df[order(out.median.df$nrmse.ndvi),],
     pch=16,
     col = col.dn15,cex = 1)
legend('bottomleft',legend = levels(out.median.df$col.dn15),
       col = palette(),pch=16)
# 
plot(mat~map,
     data = out.median.df[order(out.median.df$nrmse.ndvi),],
     pch=16,
     col = col.dn15,cex = 1)
legend('bottomleft',legend = levels(out.median.df$col.dn15),
       col = palette(),pch=16)

# 
plot(mat~map,
     data = out.median.df[order(out.median.df$nrmse.ndvi),],
     pch=16,
     col = col.k,cex = 1)
legend('bottomleft',legend = levels(out.median.df$col.k),
       col = palette(),pch=16)

plot(log(k.slope)~mat,
     data = out.median.df,
     pch=16,cex = 1)

# 
plot(clim.ra[[1]],legend=F,xlim=c(110,155),ylim=c(-45,-10),col='grey80')
points(x = out.median.df$lon, 
       y = out.median.df$lat,
       col = out.median.df$col.k,pch=15,cex=0.5)


plot(clim.ra[[1]],legend=F,
     xlim=c(110,155),
     ylim=c(-45,-10),col='grey80')
points(x = out.median.df$lon, 
       y = out.median.df$lat,
       col = out.median.df$col.dn15,pch=15,cex=0.5)
# ma1.ra <- rasterFromXYZ(out.median.df[,c('lon','lat','ma1')])
# plot(ma1.ra)
# k.slope.ra <- rasterFromXYZ(out.median.df[,c('lon','lat','k.slope')])
# plot(k.slope.ra)
# 
# c2c.ra <- rasterFromXYZ(out.median.df[,c('lon','lat','c2c')])
# plot(c2c.ra)
# 
# N0.ra <- rasterFromXYZ(out.median.df[,c('lon','lat','N0')])
# plot(N0.ra)
# 
# tr2.ra <- rasterFromXYZ(out.median.df[,c('lon','lat','tr2')])
# plot(tr2.ra)

col.vec <- names(out.median.df.narm)[-c(1,2)]

col.vec.t <- col.vec[grep('t',col.vec)]
col.vec.tn <- col.vec[grep('tn',col.vec)]
col.vec.tr <- col.vec[grep('tr',col.vec)]
col.vec.tg <- col.vec[grep('tg',col.vec)]
col.vec.mn <- col.vec[c(grep('ma',col.vec))]
col.vec.ma <- col.vec[c(grep('mn',col.vec))]
col.vec.mg <- col.vec[c(grep('mg',col.vec))]

col.vec.cn <- c('k.slope','c2c')

col.vec.error <- c('nrmse.d15n','nrmse.ndvi')


# 
r.t = stack(lapply(col.vec.t, 
                 function(j) rasterFromXYZ(out.median.df.narm[,c('lon','lat',j)])))
# 
r.tn.all = stack(lapply(col.vec.tn, 
                   function(j) rasterFromXYZ(out.median.df[,c('lon','lat',j)])))
r.tn <- calc(r.tn.all,sum)
# plot(r.tn)
# 
r.tr.all = stack(lapply(col.vec.tr, 
                        function(j) rasterFromXYZ(out.median.df[,c('lon','lat',j)])))
r.tr<- calc(r.tr.all,sum)
# plot(r.tr)
# 
r.nrmse = stack(lapply(col.vec.error, 
                   function(j) rasterFromXYZ(out.median.df[,c('lon','lat',j)])))
# 
r.tg.all = stack(lapply(col.vec.tg, 
                        function(j) rasterFromXYZ(out.median.df[,c('lon','lat',j)])))
r.tg.peak <- calc(r.tg.all[[1:2]],sum)
r.tg.stop <- calc(r.tg.all,sum)
# plot(r.tg.stop)
# 
r.mn.all = stack(lapply(col.vec.mn, 
                   function(j) rasterFromXYZ(out.median.df[,c('lon','lat',j)])))
r.mn.peak <- calc(r.mn.all[[1:2]],sum)
r.mn.stop <- calc(r.mn.all,sum)
# 
r.ma.all = stack(lapply(col.vec.ma, 
                        function(j) rasterFromXYZ(out.median.df[,c('lon','lat',j)])))
r.ma <- calc(r.ma.all,sum)
# 
r.mg.all = stack(lapply(col.vec.mg, 
                        function(j) rasterFromXYZ(out.median.df[,c('lon','lat',j)])))
r.mg <- calc(r.mg.all,sum)

# 
r.cn = stack(lapply(col.vec.cn, 
                   function(j) rasterFromXYZ(out.median.df[,c('lon','lat',j)])))
r.cn[[1]] <- r.cn[[1]]/10
# 
library(rasterVis)
levelplot(stack(list(tn = r.tn,
                     tr = r.tr,
                     tg.peak = r.tg.peak,
                     tg.stop = r.tg.stop)))
levelplot(stack(list(mg = r.mg,
                     ma = r.ma,
                     mn.peak = r.mn.peak,
                     mn.stop = r.mn.stop)))
levelplot(r.cn,pretty =T,margin=F)

levelplot(r.nrmse,at = c(seq(0,1,by=0.25),2), contour=TRUE)

bad.df <- out.median.df[out.median.df$nrmse.d15n>0.8 |
                          out.median.df$nrmse.ndvi>0.8,]

plot(r.tn)
points(bad.df$lon,bad.df$lat)






# plot.col.func <- function(){
#   tr2.ra <- rasterFromXYZ(out.median.df[,c('lon','lat','tr2')])
#   plot(tr2.ra)
#   legend('top')
# }
# # run mcmc abc  ########
# s.time <- Sys.time()
# fit.value <- get.ini.func(par.df = par.df,
#                           iter.max = 500,
#                           n.core = 25)
# par.fit.best <- unname(fit.value$optim$bestmem)
# mcmc.out <- run_MCMC_ABC(startvalue = par.fit.best,
#                          inputDF = inputs.ls,
#                          iterations = 100,
#                          thinning_interval = 1,
#                          thresh.in = 0.4,
#                          size_of_move = 1e-2)
# e.time <- Sys.time() - s.time
# # chain.thin <- mcmc(mcmc.out, thin=100)
# mcmc.out <- as.data.frame(mcmc.out)
# 
# names(mcmc.out) <- c(names(par.df),'nrmse.d15n','nrmse.ndvi')


# plot#####
mcmc.out$nrmse.sum <- mcmc.out$nrmse.d15n + mcmc.out$nrmse.ndvi
mcmc.out <- mcmc.out[complete.cases(mcmc.out),]
mcmc.out.ordered <- mcmc.out[order(mcmc.out$nrmse.sum,decreasing = F),]
# nrmse.in <- get.nmse.func(obs = in.out.df$dn15.pred,prd = in.out.df$d15n)
# nrmse.in <- get.nmse.func(obs = in.out.df$ndvi,prd = in.out.df$ndvi.pred)
# out.df <- sample_by_rejection(dat = inputs.ls,
#                               n_iterations = 100,
#                               accept_or_reject_function = 
#                                 accept_or_reject_with_squared_distance,
#                               acceptance_threshold = nrmse.in)

chain.thin <- mcmc(mcmc.out[,1:25], thin=1)
effectiveSize(chain.thin)
# autocorr.plot(chain.thin)
# plot(chain.thin)
par.mcmc <- par.fit.best#unname(apply(as.matrix(mcmc.out), 2, median))
# par.mcmc <- unname(as.matrix(mcmc.out.ordered[1,]))

out.mcmc.df <- model.de.func(parset = par.mcmc,
                             dat = inputs.ls,
                             is.evalue = F)
in.out.mcmc.df <- cbind(inputs.ls$obs,out.mcmc.df)

in.out.mcmc.df$days <- 1:nrow(in.out.mcmc.df)
in.out.mcmc.df$ndvi.pred <- 1-exp(-in.out.mcmc.df$plant.s.biomass * par.mcmc[22])

# 
get.nmse.func(obs = in.out.mcmc.df$dn15.pred,
              prd = in.out.mcmc.df$d15n)
get.nmse.func(obs = in.out.mcmc.df$ndvi,
              prd = in.out.mcmc.df$ndvi.pred)
# ndvi
plot <- ggplot(in.out.mcmc.df, aes(days,ndvi)) + 
  geom_point()+ ylim(0, 1)+ geom_smooth(method = 'gam',formula = y ~ s(x, k = 5))


plot  + 
  geom_line(aes(in.out.mcmc.df$days,
                in.out.mcmc.df$ndvi.pred),
            colour = "red")

# dn15
plot <- ggplot(in.out.mcmc.df, aes(days,dn15.pred)) + 
  geom_point()+ geom_smooth(method = 'gam',formula = y ~ s(x, k = 5))


plot  + 
  geom_line(aes(in.out.mcmc.df$days,
                in.out.mcmc.df$d15n),
            colour = "red")

# 
plot <- ggplot(in.out.mcmc.df[in.out.mcmc.df$n.up >0,], aes(days,f.value)) + 
  geom_point()+ geom_smooth(method = 'gam',formula = y ~ s(x, k = 5))
# 
plot <- ggplot(in.out.mcmc.df, aes(days,n.mass.inhib)) + 
  geom_point()+ geom_smooth()
