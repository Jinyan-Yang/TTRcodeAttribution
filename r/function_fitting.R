get.input.func <- function(ls.point.df){
  # ########
  # This function prepares the inputs to the format and time step of the model
  # the input is ls.point.df
  # it need to have the following columns:
  # Date, lon, lat
  # dn15.pred, ndvi,blue, green ,red,nir,swir1,swir2
  # tair.c, sm
  
  #########
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
  met.df <- data.frame(sm = extract(sm.ra,
                                     cbind(ls.point.sum.df$lon[1],
                                           ls.point.sum.df$lat[1]))[1,],
                       tair = extract(tair.ra,
                                       cbind(ls.point.sum.df$lon[1],
                                             ls.point.sum.df$lat[1]))[1,],
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
  
  inputs.ls$
    inputs.ls$X <- data.frame(PHOTO = seq(350,400,length.out = nrow(ls.met.df)),
                              TAIR = ls.met.df$tair.c,
                              TSOIL = ls.met.df$tair.c,
                              M = ls.met.df$sm)
  
  inputs.ls$obs <- data.frame(ls.met.df[,c("date",'lon','lat',"dn15.pred","ndvi")])
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
  if(sum(met.df$sm,na.rm=T)!=0){
    return(list(inputs.ls = inputs.ls,par.df = par.df))
  }
  
}

# func#####
fit.func <- function(ls.point.df){
  # This function performa a two-step fit:
  # 1. it fits the TTR model to data using deoptim
  # 2. it used the ABC method to get a CI for the fitting
  # 
  # read climate
  input.val.ls <- try(get.input.func(ls.point.df))
  # use deopt to fit
  fit.value <- try(get.ini.func(par.df = input.val.ls$par.df,
                                inputs.ls = input.val.ls$inputs.ls,
                                iter.max = 500,
                                n.core = 29))
  # get the fit results
  par.fit.best <- try(unname(fit.value$optim$bestmem))
  # use abc to get uncertainty
  mcmc.out <- try(run_MCMC_ABC(startvalue = par.fit.best,
                               inputDF = input.val.ls$inputs.ls,
                               iterations = 100,
                               thinning_interval = 1,
                               thresh.in = 0.5,
                               size_of_move = 1e-2))
  
  if(class(fit.value) == 'try-error' |
     class(input.val.ls) == 'try-error'){
    
    class(mcmc.out) = 'try-error'
  }
  # give results
  if(class(mcmc.out)[1] != 'try-error'){
    mcmc.out <- as.data.frame(mcmc.out)
    
    names(mcmc.out) <- c(names(input.val.ls$par.df),
                         'nrmse.d15n',
                         'nrmse.ndvi')
    mcmc.out$lon <- as.numeric(ls.point.df$lon[1])
    mcmc.out$lat <- as.numeric(ls.point.df$lat[1])
    mcmc.out <- mcmc.out[!duplicated(mcmc.out),]
    return(mcmc.out)
  }else{
    return(data.frame(lon=NA))
  }
}