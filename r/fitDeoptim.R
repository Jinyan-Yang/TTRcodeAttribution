library(DEoptim)
library(foreach)
library(doParallel)
library(plante)
library(lubridate)
library(doBy)
library(randomForest)
library(raster)
library(parallel)

source('r/readERA5.r')
# read rf model fit fr dn15
rf.fit <- readRDS('\\\\fs1-cbr.nexus.csiro.au/{mmrg}/work/users/yan190/repo/delta_n_15/cache/rf.kFold.n15.rds')
# read landsat
inputs.ls <- readRDS(file="MyData_SA-BFA-BNP.rds")

landsat.ts.ls <- readRDS('\\\\fs1-cbr.nexus.csiro.au/{mmrg}/work/users/yan190/repo/delta_n_15/cache/landsat.ts.rds')

ls.point.df <- landsat.ts.ls[[42627]]#landsat.ts.ls[[4113]]
# 
ls.point.df$mon <- month(ls.point.df$date)
ls.point.df$yr <- year(ls.point.df$date)

ls.point.df <- as.data.frame(apply(ls.point.df,2,as.numeric))
ls.point.df$dn15.pred <- predict(rf.fit,ls.point.df)

ls.point.sum.df <- summaryBy( dn15.pred + ndvi +
                                blue+ green +red+
                                nir+swir1+swir2 ~  lon+lat+mon+yr,
                              data = ls.point.df,
                              FUN = median,na.rm=T,keep.names = T)

# extract met data
met.df <- data.frame(sm = c(extract(sm.2000.2010,cbind(ls.point.sum.df$lon[1],
                                                     ls.point.sum.df$lat[1])),
                            extract(sm.2010.2020,cbind(ls.point.sum.df$lon[1],
                                                       ls.point.sum.df$lat[1]))),
                     tair = c(extract(tair.2000.2010,cbind(ls.point.sum.df$lon[1],
                                                         ls.point.sum.df$lat[1])),
                              extract(tair.2010.2020,cbind(ls.point.sum.df$lon[1],
                                                         ls.point.sum.df$lat[1]))),
                     date = c(names(sm.2000.2010),names(sm.2010.2020)))

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
# inputs.ls$initial <- data.frame(ms = 1,
#                                 mr = 1,
#                                 d15n = 0)
# # dat = inputs.ls
parm.df <- data.frame(nm = inputs.ls$parm.names,
                      value = inputs.ls$pars[c(inputs.ls$pos.parset,
                                               inputs.ls$pos.sigma,
                                               inputs.ls$pos.ou,
                                               inputs.ls$pos.x0)])
model.de.func <- function(parset,
                          dat,is.evalue=TRUE){
  
  # parm = dat$pars
  # parset <- pars#parm[dat$pos.parset]
  # sigma <- parm[dat$pos.sigma]
  # ou <- parm[dat$pos.ou] # not used in this function
  # x0 <- parm[dat$pos.x0]
  
  # initial <- data.frame(ms = parset[23],
  #                       mr = parset[24],
  #                       d15n = parset[25])
  
  initials <- data.frame(ms = parset[23],
                        mr = parset[24],
                        d15n = parset[25])
  #old model####
  # d15n.pred.df <- ttr_d15n(steps = dat$steps,
  #                          initials = initial,
  #                          TAIR = dat$X[,"TAIR"],
  #                          TSOIL = dat$X[,"TSOIL"],
  #                          M = dat$X[,"M"],
  #                          # N = dat$X[,"N"],
  #                          # FIRE = dat$X[,"FIRE"],
  #                          A = dat$X[,"PHOTO"],
  #                          Kl = dat$Kl,
  #                          gs = dat$gs,
  #                          gr= dat$gr, 
  #                          KM= dat$KM,
  #                          KA= dat$KA,
  #                          Jc= dat$Jc,
  #                          Jn = dat$Jn,
  #                          q = dat$q,
  #                          RHOc= dat$RHOc,
  #                          RHOn=dat$RHOn,
  #                          Fc= dat$Fc,
  #                          Fn= dat$Fn,
  #                          ma1= parset[1],    
  #                          ma2= parset[1]+parset[2],  
  #                          
  #                          tn1= parset[3],     
  #                          tn2= parset[3]+parset[4],    
  #                          
  #                          mn1= parset[5],   
  #                          mn2= parset[5]+parset[6],     
  #                          mn3= parset[5]+parset[6]+parset[7],     
  #                          mn4= parset[5]+parset[6]+parset[7]+parset[8], 
  #                          
  #                          tg1= parset[9],     
  #                          tg2= parset[9]+parset[10],    
  #                          tg3= parset[9]+parset[10]+parset[11],    
  #                          tg4= parset[9]+parset[10]+parset[11]+parset[12],    
  #                          
  #                          mg1= parset[1],    
  #                          mg2= parset[1]+parset[2],     
  #                          
  #                          tr1= parset[3],     
  #                          tr2= parset[3] + parset[4],    
  #                          
  #                          f1= parset[1],     
  #                          f2= parset[1] + parset[2], 
  #                          
  #                          A0= parset[13],
  #                          N0= parset[14],
  #                          uMs= parset[15],#sigma[1]/1e3,
  #                          uMr= parset[16],#sigma[2]/1e3,
  #                          uCs= parset[17],#sigma[3]/1e3,
  #                          uCr= parset[18],#sigma[4]/1e3,
  #                          uNs= parset[19],#sigma[5]/1e3,
  #                          uNr= parset[20],#sigma[6]/1e3,
  #                          # new param for d15N
  #                          n15.rich = 6.93,
  #                          ndepleted = -5.96,
  #                          k.slope = parset[21],#5e4,
  #                          ud15N = parset[22],#sigma[6]/1e3,
  #                          fc.700 = 1.4
  # )
  #model######
  d15n.pred.df <- ttr_d15n(steps = dat$steps,
                           initials = initials,
                           TAIR = dat$X[,"TAIR"],
                           TSOIL = dat$X[,"TSOIL"],
                           M = dat$X[,"M"],
                           # N = dat$X[,"N"],
                           # FIRE = dat$X[,"FIRE"],
                           A = dat$X[,"PHOTO"],
                           Kl = dat$Kl,
                           gs = dat$gs/7,
                           gr= dat$gr/7,
                           KM= dat$KM,
                           KA= dat$KA,
                           Jc= dat$Jc,
                           Jn = dat$Jn,
                           q = dat$q,
                           RHOc= dat$RHOc,
                           RHOn=dat$RHOn,
                           Fc= dat$Fc,
                           Fn= dat$Fn,
                           ma1= parset[1],
                           ma2= parset[1]+parset[2],
                           
                           tn1= parset[3],
                           tn2= parset[3]+parset[4],
                           
                           mn1= parset[5],
                           mn2= parset[5]+parset[6],
                           
                           mn3= parset[5]+parset[6]+parset[7],
                           mn4= parset[5]+parset[6]+parset[7]+parset[8],
                           
                           tg1= parset[9],
                           tg2= parset[9]+parset[10],
                           tg3= parset[9]+parset[10]+parset[11],
                           tg4= parset[9]+parset[10]+parset[11]+parset[12],
                           
                           mg1= parset[13],
                           mg2= parset[13]+parset[14],
                           
                           tr1= parset[15],
                           tr2= parset[15] + parset[16],
                           
                           f1= parset[17],
                           f2= parset[17] + parset[18],
                           
                           A0= parset[19],
                           N0= parset[20],
                           uMs= 0,#parset[21],#sigma[1]/1e3,
                           uMr= 0,#parset[22],#sigma[2]/1e3,
                           uCs= 0,#parset[23],#sigma[3]/1e3,
                           uCr= 0,#parset[24],#sigma[4]/1e3,
                           uNs= 0,#parset[25],#sigma[5]/1e3,
                           uNr= 0,#parset[26],#sigma[6]/1e3,
                           # new param for d15N
                           n15.rich = 6.93,
                           ndepleted = -5.96,
                           k.slope = parset[21],#5e4,
                           ud15N = 0,#parset[28],#sigma[6]/1e3,
                           fc.700 = 1.4
  )
  
  #output #####
  #standardise with sd
  sd.cover <- sd(dat$obs$ndvi,na.rm = T)
  
  dat$obs$cover.pred <- 1 - exp(-parset[22] * d15n.pred.df$c.s)
  dat$obs$d15n <- d15n.pred.df$d15n 
  
  # d15n.pred.df$pred.cover <- d15n.pred.df$c.s * parset[29]  
  
  df.tmp <- dat$obs[!is.na(dat$obs$ndvi),]
  
  resid.cover <- ((df.tmp$cover.pred - df.tmp$ndvi)/sd.cover)^2
  # resid.cover[is.na(resid.cover)] <- 1e3
  # 
  sd.d15n <- sd(df.tmp$dn15.pred,na.rm = T)
  
  resid.swc <- ((df.tmp$dn15.pred - df.tmp$d15n)/sd.d15n)^2
  # resid.swc[is.na(resid.swc)] <- 1e3
  # 
  if(is.evalue){
    return(sum(resid.cover) + sum(resid.swc))
  }else{
    return(d15n.pred.df)
  }

}
# # inputs.ls$parm.names
# inputs.ls$pars[5] <- 17
# inputs.ls$pars[6:8] <- 5
# inputs.ls$pars[9] <- 17
# inputs.ls$pars[10:12] <- 5
# inputs.ls$pars[3] <- 17
# inputs.ls$pars[4] <- 15
# inputs.ls$pars[15] <- 17
# inputs.ls$pars[16] <- 15
# model.de.func(parset = c(inputs.ls$pars[1:20],rep(1e-5,6),100,1e-5),
#               dat = inputs.ls)


# func to use DEoptim to calaulate initial values
get.ini.func <- function(par.df,...){
  # on.exit(stopCluster(cl = NULL))
  # setting control parameters and limits to values
  lower <- unlist(par.df['min',])
  upper <- unlist(par.df['max',]) 
  NPmax <- 500
  maxiter <- 1000
  # 
  set.seed(1935)
  OptBB.de.fit <- DEoptim(fn= model.de.func,
                          lower = lower,
                          upper = upper,
                          dat = inputs.ls,
                          DEoptim.control(VTR = 1,
                                          NP = NPmax,
                                          itermax = maxiter,
                                          trace = 10,
                                          parallelType = 'foreach',
                                          parVar = list('ttr_d15n'),
                                          cluster = makeCluster(20)))
  
  # OptBB.de.fit <- DEoptim(fn= model.de.func,
  #                         lower = lower,
  #                         upper = upper,
  #                         dat = inputs.ls,
  #                         DEoptim.control(VTR = 1,
  #                                         NP = NPmax,
  #                                         itermax = maxiter,
  #                                         parallelType ="none"
  #                                         ))
  # Sys.sleep(10)
  initial.vec <- OptBB.de.fit#unname(OptBB.de.fit$optim$bestmem)
  return(initial.vec)
}

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
                     c2c = c(1,20),
                     ms = c(0.1,5),
                     mr = c(0.1,15),
                     d15n = c(d15n.mean - d15n.sd,d15n.mean + d15n.sd)
                     )


# par.df <- data.frame(ma1 = c(0.05,0.5),
#                      ma2 = c(0.05,0.5),
#                      
#                      tn1 = c(0,20),
#                      tn2 = c(1,20),
#                      
#                      mn1 = c(0.05,0.25),
#                      mn2 = c(0.1,0.25),
#                      mn3 = c(0.1,0.25),
#                      mn4 = c(0.1,0.25),
#                      
#                      tg1 = c(0,20),
#                      tg2 = c(1,10),
#                      tg3 = c(1,10),
#                      tg4 = c(1,10),
#                      
#                      A0 = c(0.01,0.4),
#                      N0 = c(0.005,0.04),
#                      
#                      uMs = c(1e-6,1e-5),#sigma[1]/1e3,
#                      uMr = c(1e-6,1e-5),#sigma[2]/1e3,
#                      uCs = c(1e-6,1e-5),#sigma[3]/1e3,
#                      uCr = c(1e-6,1e-5),#sigma[4]/1e3,
#                      uNs = c(1e-6,1e-5),#sigma[5]/1e3,
#                      uNr = c(1e-6,1e-5),#sigma[6]/1e3,
#                      
#                      k.slope = c(10,500),
#                      ud15N = c(1e-6,1e-5),
#                      c2c = c(0.1,2),
#                      ms = c(1,500),
#                      mr = c(1,50),
#                      d15n = c(-10,10)
# )

row.names(par.df) <- c('min','max')
# ######
fit.value <- get.ini.func(par.df = par.df)
par.fit.best <- unname(fit.value$optim$bestmem)

out.df <- model.de.func(parset = par.fit.best,
                        dat = inputs.ls,is.evalue = F)
in.out.df <- cbind(ls.met.df,out.df)

in.out.df$days <- 1:nrow(in.out.df)

# plot(c(c.s* par.fit.best[29])~ndvi,data = in.out.df)
# plot(in.out.df$c.s * par.fit.best[29],pch=16,col='red')
# points(in.out.df$ndvi,pch=16,col='grey')
# 
# plot(in.out.df$dn15.pred,pch=16,col='red')
# points(in.out.df$d15n,pch=16,col='grey')
# 
# plot(ls.met.df$sm)

library(ggplot2)
plot <- ggplot(in.out.df, aes(days,ndvi)) + 
  geom_point()+ ylim(0, 1)


plot + geom_smooth(method = 'gam',formula = y ~ s(x, k = 80)) + 
  geom_point(aes(in.out.df$days,
                 1-exp(-in.out.df$c.s * par.fit.best[22])),
             colour = "red")

# x <- in.out.df[!is.na(in.out.df$ndvi),]
in.out.df.sub <- in.out.df#[in.out.df$yr >2010,]
plot.d15n <- ggplot(in.out.df.sub, aes(days,dn15.pred)) + geom_point()


plot.d15n + geom_smooth(method = 'gam',formula = y ~ s(x, k = 80)) + 
  geom_point(aes(in.out.df.sub$days,
                 in.out.df.sub$d15n),
             colour = "red")
