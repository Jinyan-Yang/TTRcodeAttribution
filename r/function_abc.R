
# source('r/fitDeoptim.R')
# modified from https://rpubs.com/boussau/BasicABC
get.nmse.func <- function(obs,prd){
  
  val.m <- mean((obs - prd)^2,na.rm=T)
  val.quantile <- unname(quantile(obs,probs = c(0.05,0.95),na.rm=T))
  # sqrt( / 
  #   (quantile(obs,probs = 0.95,na.rm=T) - 
  #      quantile(obs,probs = 0.05,na.rm=T))
  
  # print(val.m)
  # print(val.quantile)
  
  return(sqrt(val.m) / (val.quantile[2] - val.quantile[1]))
}

compute_quantiles <- function(data) {
  return (quantile(data, probs=c(0.1, 0.5, 0.9),na.rm=T))
}
# First method to compare a simulated sample to the observed data
compare_quantiles_with_squared_distance <- function (true, simulated) {
  distance = sqrt(sum(mapply(function(x,y) (x-y)^2, true, simulated),na.rm=T))
  return(distance)
}

# Accept or reject based on the first method to compare a simulated sample to the observed data
accept_or_reject_with_squared_distance <- function (true, simulated, acceptance_threshold) {
  distance = #compare_quantiles_with_squared_distance(compute_quantiles(true), compute_quantiles(simulated))
    get.nmse.func(true,simulated)
  # print(distance)
  # print(distance)
  # print(acceptance_threshold)
  if((distance < acceptance_threshold) ) return(T) else return(F)
}

draw_no <- function(min.in,max.in){
  runif(1,min = min.in,max = max.in)
}

#sample func####
sample_by_rejection <- function (dat, n_iterations, 
                                 acceptance_threshold = 0.2, 
                                 accept_or_reject_function) {
  true_data <- dat$obs
  number_of_data_points = length(true_data)
  accepted_or_rejected <- vector(length = n_iterations)
  para.full <- list()

  # get par list######
  d15n.mean <- mean(true_data$dn15.pred,na.rm=T)
  d15n.sd <- sd(true_data$dn15.pred,na.rm=T)
  ndvi.sd <- sd(true_data$ndvi,na.rm=T)

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
                       d15n = c(d15n.mean - d15n.sd,d15n.mean + d15n.sd),
                       row.names = c('min','max'))
  

  #iteration######
  for (i in 1:n_iterations){
    
    # param.val <- apply(par.df, 2,
    #                    function(vec.in)draw_no(min.in = vec.in[1],
    #                                            max.in = vec.in[2]))
    param.val <- sapply(par.fit.best, function(vals){
      rnorm(1,mean = vals,sd = vals*1e-2)
    })
    # 
    para.full[[i]] <- as.data.frame(t(param.val))
    parset <- unname(param.val)
    # get inital state var values
    initials <- data.frame(ms = parset[23],
                           mr = parset[24],
                           d15n = parset[25])
    # start model####
    simulated_data <- ttr_d15n(steps = dat$steps,
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
    
    
    simulated_data$ndvi.pred <- 1-exp(-simulated_data$plant.s.biomass * 
                                        parset[22])
    # get error####
    ndvi.accept = accept_or_reject_function(true_data$ndvi, 
                                            simulated_data$ndvi.pred, 
                                            acceptance_threshold)
    n.accept = accept_or_reject_function(true_data$dn15.pred, 
                                         simulated_data$d15n, 
                                         acceptance_threshold)
    
    accepted_or_rejected[i] <- all(c(ndvi.accept,n.accept))
  }
  
  out.df <- do.call(rbind,para.full)
  out.df$accept <- accepted_or_rejected
  
  return(out.df)
}
# #####
ABC_acceptance <- function(parset,dat,threshold){
  out.ls <- list()
  # prior to avoid negative standard deviation
  if (any(parset[-length(parset)] <= 0)){
    out.ls$diff.n15 <- NA
    out.ls$diff.ndvi <- NA
    out.ls$trueFalse <- F
    return(out.ls) 
  } 
  
  # start model####
  simulated_data <- ttr_d15n(steps = dat$steps,
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
  
  ######
  simulated_data$ndvi.pred <- 1-exp(-simulated_data$plant.s.biomass * 
                                      parset[22])
  # print(simulated_data)
  simulated_data <- cbind(simulated_data,dat$obs)
  # comparison with the observed summary statistics
  diff.n15 <- get.nmse.func(obs = simulated_data$dn15.pred,
                            prd = simulated_data$d15n)
  diff.ndvi <- get.nmse.func(obs = simulated_data$ndvi,
                             prd = simulated_data$ndvi.pred)
  if((diff.n15 < threshold) & (diff.ndvi < threshold)){
    out.ls$diff.n15 <- diff.n15
    out.ls$diff.ndvi <- diff.ndvi
    out.ls$trueFalse <- T
    return(out.ls)
  }  else{
    out.ls$diff.n15 <- diff.n15
    out.ls$diff.ndvi <- diff.ndvi
    out.ls$trueFalse <- F
    return(out.ls)
  } 
}
# simulated_data = in.out.df

# we plug this in a standard MCMC, 
# with the metropolis acceptance exchanged for the ABC acceptance
#MCMC####
run_MCMC_ABC <- function(startvalue, 
                         inputDF,
                         iterations, 
                         thinning_interval=1, 
                         size_of_move=0.7,
                         thresh.in = 0.5){
  n.par <- length(startvalue)
  chain = array(dim = c(iterations+1,n.par))
  chain[1,] = startvalue
  out.accpt <- ABC_acceptance(startvalue,
                              dat = inputDF,
                              threshold = thresh.in)
  d15n.nrmse <- c(out.accpt$diff.n15)
  ndvi.nrmse <- c(out.accpt$diff.ndvi)
  for (i in 1:iterations){
    # proposalfunction
    
    proposal = rnorm(n.par,
                     mean = chain[i,], 
                     sd= abs(size_of_move * chain[i,]))
    
    out.accpt.it <- ABC_acceptance(proposal,
                                   dat = inputDF,
                                   threshold = thresh.in)
    
    if(out.accpt.it$trueFalse){
      chain[i+1,] = proposal
      d15n.nrmse[i+1] <- out.accpt.it$diff.n15
      ndvi.nrmse[i+1] <- out.accpt.it$diff.ndvi
    }else{
      chain[i+1,] = chain[i,]
      d15n.nrmse[i+1] <- d15n.nrmse[i]
      ndvi.nrmse[i+1] <- ndvi.nrmse[i]
    }
    

  }
  return(do.call(cbind,list(chain,
                            as.matrix(d15n.nrmse,ncol=1),
                            as.matrix(ndvi.nrmse,ncol=1))))
}
# mcmc(chain, thin=thinning_interval)
