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
# 

landsat.ts.ls.sample <- readRDS('cache/sampleAU.ls.rds')

length(landsat.ts.ls.sample)
# func#####
get.input.func <- function(ls.point.df){#,
                           # rf.fit = readRDS('//fs1-cbr.nexus.csiro.au/{mmrg}/work/users/yan190/repo/delta_n_15/cache/rf.kFold.n15.rds')){
  # ls.point.df <- landsat.ts.ls[[42627]]#landsat.ts.ls[[4113]]
  # 
  ls.point.df$mon <- month(ls.point.df$date)
  ls.point.df$yr <- year(ls.point.df$date)
  
  ls.point.df <- as.data.frame(apply(ls.point.df,2,as.numeric))
  # ls.point.df$dn15.pred <- predict(rf.fit,ls.point.df)
  
  ls.point.sum.df <- summaryBy( dn15.pred + ndvi +
                                  blue+ green +red+
                                  nir+swir1+swir2 ~  lon+lat+mon+yr,
                                data = ls.point.df,
                                FUN = median,na.rm=T,keep.names = T)
  
  # extract met data
  met.df <- data.frame(sm = (extract(sm.ra,
                                    cbind(ls.point.sum.df$lon[1],
                                          ls.point.sum.df$lat[1])))[1,],
                       tair = (extract(tair.ra,
                                      cbind(ls.point.sum.df$lon[1],
                                            ls.point.sum.df$lat[1])))[1,],
                       date = names(sm.ra))
  
  met.df$date <- as.Date(met.df$date,'X%Y.%m.%d')
  met.df$yr <- year(met.df$date)
  met.df$mon <- month(met.df$date)
  # 
  ls.met.df <- merge(ls.point.sum.df,met.df,by = c('mon','yr'),all=T)
  ls.met.df <- ls.met.df[order(ls.met.df$date),]
  ls.met.df <- ls.met.df[!is.na(ls.met.df$sm),]
  ls.met.df$tair.c <- ls.met.df$tair - 272.15
  # # month()
  inputs.ls <- readRDS(file="MyData_SA-BFA-BNP.rds")
  inputs.ls$X <- data.frame(PHOTO = seq(350,400,length.out = nrow(ls.met.df)),
                            TAIR = ls.met.df$tair.c,
                            TSOIL = ls.met.df$tair.c,
                            M = ls.met.df$sm)

  inputs.ls$obs <- data.frame(ls.met.df[,c("date","dn15.pred","ndvi")])
  inputs.ls$steps <- nrow(ls.met.df)
  
  inputs.ls <- inputs.ls 
  
  ndvi.mean <- mean(inputs.ls$obs$ndvi,na.rm=T)
  ndvi.sd <- sd(inputs.ls$obs$ndvi,na.rm=T)
  
  d15n.mean <- mean(inputs.ls$obs$dn15.pred,na.rm=T)
  d15n.sd <- sd(inputs.ls$obs$dn15.pred,na.rm=T)
  # make input par list####
  par.df <- data.frame(ma1 = c(0.01,0.15),
                       ma2 = c(0.05,0.25),
                       tn1 = c(5,15),
                       tn2 = c(5,20),
                       mn1 = c(0.01,0.1),
                       mn2 = c(0.01,0.1),
                       mn3 = c(0.01,0.1),
                       mn4 = c(0.01,0.1),
                       tg1 = c(5,10),
                       tg2 = c(1,10),
                       tg3 = c(1,10),
                       tg4 = c(1,10),
                       mg1 = c(0.01,0.1),
                       mg2 = c(0.05,0.2),
                       tr1 = c(1,10),     
                       tr2 = c(5,20),
                       f1 = c(0.01,0.1),     
                       f2 = c(0.05,0.2),
                       A0 = c(0.05,0.4),
                       N0 = c(0.005,0.04),
                       # uMs = c(1,10),#sigma[1]/1e3,
                       # uMr = c(1,10),#sigma[2]/1e3,
                       # uCs = c(1,10),#sigma[3]/1e3,
                       # uCr = c(1,10),#sigma[4]/1e3,
                       # uNs = c(1,10),#sigma[5]/1e3,
                       # uNr = c(1,10),#sigma[6]/1e3,
                       k.slope = c(10,500),
                       # ud15N = c(1,10),
                       c2c = c(1,30),
                       ms = c(0.1,5),
                       mr = c(0.1,15),
                       d15n = c(-10,10),#c(d15n.mean - d15n.sd,d15n.mean + d15n.sd),
                       row.names = c('min','max'))
  
  return(list(inputs.ls = inputs.ls,par.df = par.df))
}

fit.func <- function(ls.point.df){
  input.val.ls <- get.input.func(ls.point.df)
  
  fit.value <- try(get.ini.func(par.df = input.val.ls$par.df,
                                inputs.ls = input.val.ls$inputs.ls,
                                iter.max = 500,
                                n.core = 30))
  if(class(fit.value) != 'try-error'){
    par.fit.best <- unname(fit.value$optim$bestmem)
    mcmc.out <- run_MCMC_ABC(startvalue = par.fit.best,
                             inputDF = input.val.ls$inputs.ls,
                             iterations = 100,
                             thinning_interval = 1,
                             thresh.in = 0.4,
                             size_of_move = 1e-2)
    
    mcmc.out <- as.data.frame(mcmc.out)
    
    names(mcmc.out) <- c(names(input.val.ls$par.df),
                         'nrmse.d15n',
                         'nrmse.ndvi')
    mcmc.out$lon <- as.numeric(ls.point.df$lon[1])
    mcmc.out$lat <- as.numeric(ls.point.df$lat[1])
    return(mcmc.out)
  }else{
    return(data.frame(lon=NA))
  }
 
}

#fitting even sample 0.1 degree####
s.time <- Sys.time()
out.ls <- lapply(landsat.ts.ls.sample, fit.func)
e.time <- Sys.time() - s.time

saveRDS(out.ls,'cache/fittedParAuSites.rds')

#fitting havplot####
s.time <- Sys.time()
out.ls <- lapply(landsat.n15.ls[2], fit.func)
e.time <- Sys.time() - s.time
# ls.point.df <- landsat.n15.ls[[1]]
saveRDS(out.ls,'cache/fittedParHavPlot.rds')

# length(landsat.n15.ls)
# length(out.ls)
# nrow(out.median.df)
# get the poorly fitted sites
vec.1 <- which(out.median.df$nrmse.ndvi>0.6)
vec.2 <- which(out.median.df$nrmse.d15n>0.6)
test <- intersect(vec.1,vec.2)

# redo fitting with more data
s.time <- Sys.time()
out.ls.bad <- lapply(landsat.ts.ls.sample[test], fit.func)
e.time <- Sys.time() - s.time
# analysis#####
out.ls <- readRDS('cache/fittedParAuSites.rds')
# out.ls <- lapply(out.ls, function(x) if(is.null(x)) data.frame(lon = NA) else x)
# 
out.median.ls <- lapply(out.ls, function(df){
  # print(head(df))
  if(!is.null(nrow(df))){
    require(doBy)
    df.sum <- summaryBy(.~lon+lat,
                        data = df,
                        FUN=c(median),
                        na.rm=T,keep.names = T)
    
    return(df.sum)
  }else{
    return(data.frame(lon = NA))
  }
 
})
out.median.df <- dplyr::bind_rows(out.median.ls)
# 
out.bad.ls <- lapply(out.ls.bad, function(df){
  # print(head(df))
  if(!is.null(nrow(df))){
    require(doBy)
    df.sum <- summaryBy(.~lon+lat,
                        data = df,
                        FUN=c(median),
                        na.rm=T,keep.names = T)
    
    return(df.sum)
  }else{
    return(data.frame(lon = NA))
  }
  
})
out.bad.df <- dplyr::bind_rows(out.bad.ls)
out.bad.df.narm <- out.bad.df[complete.cases(out.bad.df),]

# 
library(geodata)
clim.ra <- worldclim_country("Australia", var="bio", res = 10,path='wc10/')

out.median.df$map <- extract(clim.ra[[12]],cbind(out.median.df$lon,
                                                 out.median.df$lat))[,1]
out.median.df$mat <- extract(clim.ra[[1]],cbind(out.median.df$lon,
                                                out.median.df$lat))[,1]
out.median.df$nrmse.d15n
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

col.vec <- names(out.median.df)[-c(1,2)]

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
                 function(j) rasterFromXYZ(out.median.df[,c('lon','lat',j)])))
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
