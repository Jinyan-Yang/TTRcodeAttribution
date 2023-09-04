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
                           TAIR = dat$X[,"TAIR"],
                           TSOIL = dat$X[,"TSOIL"],
                           M = dat$X[,"M"],
                           # N = dat$X[,"N"],
                           # FIRE = dat$X[,"FIRE"],
                           A = dat$X[,"PHOTO"],
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
  # #standardise with sd
  # sd.cover <- sd(dat$obs$ndvi,na.rm = T)
  # 
  dat$obs$cover.pred <- 1 - exp(-parset[22] * d15n.pred.df$plant.s.biomass)
  dat$obs$d15n <- d15n.pred.df$d15n
  df.tmp <- dat$obs[!is.na(dat$obs$ndvi),]
  
  # #sd based metric
  # d15n.pred.df$pred.cover <- d15n.pred.df$c.s * parset[29]
  # 
  
  # 
  # resid.cover <- ((df.tmp$cover.pred - df.tmp$ndvi)/sd.cover)^2
  # 
  # resid.cover[is.na(resid.cover)] <- 1e3
  # 
  # # 
  # sd.d15n <- sd(df.tmp$dn15.pred,na.rm = T)
  # 
  # resid.swc <- ((df.tmp$dn15.pred - df.tmp$d15n)/sd.d15n)^2
  # resid.swc[is.na(resid.swc)] <- 1e3
  
  # nrmse based sum error
  diff.n15 <- get.nmse.func(obs = df.tmp$dn15.pred,
                            prd = df.tmp$d15n)
  diff.ndvi <- get.nmse.func(obs = df.tmp$ndvi,
                             prd = df.tmp$cover.pred)
  
  # 
  if(is.evalue){
    return(diff.n15 + diff.ndvi)
  }else{
    
    return(d15n.pred.df)
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
                          DEoptim.control(VTR = 1,
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



