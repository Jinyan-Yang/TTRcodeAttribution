# functions to fit deoptim #####
#opt model######
model.de.func <- function(parset,
                          dat,is.evalue=TRUE){
  # anything in the parset is to be fitted
  # the order in parset is obviously critical
  
  # initial values
  initials <- data.frame(ms = parset[23],
                         mr = parset[24],
                         d15n = parset[25])
  
  #model######
  d15n.pred.df <- ttr_d15n(steps = dat$steps,
                           initials = initials,
                           TAIR = unname(unlist(dat$X[,"TAIR"])),
                           TSOIL = unname(unlist(dat$X[,"TSOIL"])),
                           M = unname(unlist(dat$X[,"M"])),
                           # N = dat$X[,"N"],
                           # FIRE = dat$X[,"FIRE"],
                           A = unname(unlist(dat$X[,"PHOTO"])),
                           Kl = 0.5/30,
                           gs = 1/0.025 / 60,# dat$gs/7,
                           gr = 1/0.025 / 60,#dat$gr/7,#
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
                           
                           mn3= 1,#parset[5]+parset[6]+parset[7],
                           mn4= 1, #parset[5]+parset[6]+parset[7]+parset[8],
                           
                           n15.rich = parset[7],#6.93,#
                           ndepleted = parset[8],#-5.96,#
                           
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
                           # n15.rich = parset[26],#6.93,#
                           # ndepleted = parset[27],#-5.96,#
                           
                           # n15.rich = 6.93,#
                           # ndepleted = -5.96,#
                           
                           k.slope = parset[21],#5e4,
                           ud15N = 0,#parset[28],#sigma[6]/1e3,
                           fc.700 = 1.4
  )
  
  #output #####
  # #standardise with max

  # 
  dat$obs$cover.pred <- 1 - exp(-parset[22] * d15n.pred.df$plant.s.biomass)
  dat$obs$d15n <- d15n.pred.df$d15n
  df.tmp <- dat$obs[!is.na(dat$obs$ndvi),]
  
  dat$obs$n.up <- d15n.pred.df$n.up
  dat$obs$n.mass.inhib <- d15n.pred.df$n.mass.inhib
  
  dat$obs$c.s <- d15n.pred.df$c.s
  dat$obs$n.s <- d15n.pred.df$n.s
  dat$obs$n.r <- d15n.pred.df$n.r
  dat$obs$plant.biomass <- d15n.pred.df$plant.biomass
  dat$obs$plant.s.biomass <- d15n.pred.df$plant.s.biomass
  
  
  # #sd based metric
  # d15n.pred.df$pred.cover <- d15n.pred.df$c.s * parset[29]
  # 
  
  # 
  sd.cover <- 2*IQR(dat$obs$ndvi,na.rm = T)
  resid.cover <- ((df.tmp$cover.pred - df.tmp$ndvi)/sd.cover)^2

  resid.cover[is.na(resid.cover)] <- 1e3

  #
  sd.d15n <- 2*IQR(df.tmp$dn15.pred,na.rm = T)

  resid.swc <- ((df.tmp$dn15.pred - df.tmp$d15n)/sd.d15n)^2
  resid.swc[is.na(resid.swc)] <- 1e3
  # 
  get.nmse.func <- function(obs,prd){
    
    val.m <- mean((obs - prd)^2,na.rm=T)
    val.quantile <- #unname(quantile(obs,probs = c(0.05,0.95),na.rm=T))
    2*IQR(obs,na.rm = T)
    # sqrt( / 
    #   (quantile(obs,probs = 0.95,na.rm=T) - 
    #      quantile(obs,probs = 0.05,na.rm=T))
    
    # print(val.m)
    # print(val.quantile)
    
    return(sqrt(val.m) / val.quantile)#(val.quantile[2] - val.quantile[1]))
  }
  # # # nrmse based sum error
  resid.swc <- get.nmse.func(obs = df.tmp$dn15.pred,
                             prd = df.tmp$d15n)
  resid.cover <- get.nmse.func(obs = df.tmp$ndvi,
                               prd = df.tmp$cover.pred)
  # 
  dat$obs$nrmse.ndvi <- resid.cover
  dat$obs$nrmse.d15n <- resid.swc
  # 
  
  out.info <- sum(resid.swc) + sum(resid.cover)
  if(is.evalue){
    return(out.info)
  }else{
    return(dat$obs)
  }
  
}

# func to use DEoptim to calculate initial values
get.ini.func <- function(par.df,inputs.ls,iter.max=50,n.core = 20){

  # setting control parameters and limits to values
  lower <- unlist(par.df['min',])
  upper <- unlist(par.df['max',]) 
  NPmax <- 300#$population in each iteration
  # 
  set.seed(1935)
  OptBB.de.fit <- DEoptim(fn= model.de.func,
                          lower = lower,
                          upper = upper,
                          dat = inputs.ls,
                          DEoptim.control(VTR = 0.05,
                                          reltol = 0.1,
                                          NP = NPmax,
                                          itermax = iter.max,
                                          trace = 50,
                                          parallelType = 'foreach',
                                          parVar = list('ttr_d15n'),
                                          cluster = makeCluster(n.core)))
  
  
  initial.vec <- OptBB.de.fit
  return(initial.vec)
}



