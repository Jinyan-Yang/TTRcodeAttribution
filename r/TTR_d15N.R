# paper
# https://www.nature.com/articles/s41561-022-01114-x

# 
# make sure the parameter values are in per day so it can be converted properly
# 
# acutal function#####
ttr_d15n <- function( steps, #number of total steps
                      nDay=30,#number of days in each steps
                      initials, #intial values of the pools
                      TAIR, # air temperature
                      TSOIL, #soil temperature
                      M, #soil moisture
                      # N, 
                      # FIRE, 
                      A,  #co2 concentration
                      Kl, #carbon senescence rate
                      gs, # shoot growth rate
                      gr, #root growth rate
                      KM, #mass based growth limitation factor
                      A0, #carbon assimilation rate
                      N0,# nitrogen uptake rate
                      KA, #mass based yotake limitation factor
                      Jc, #inhibition of source activity when substrate concentrations are high
                      Jn, #inhibition of source activity when substrate concentrations are high
                      q, #resistent power not used: q=1
                      RHOc, # specifc transport resistances for carbon
                      RHOn,#specifc transport resistances for nitrogen
                      Fc, #fraction of structural C in structural dry matter
                      Fn, #fraction of structural N in structural dry matter
                      # //  to be fitted environmental dependencies
                      # soil moisture control on carbon assimilation
                      ma1,
                      ma2,
                      # temperature control on nitrogen uptake
                      tn1,
                      tn2, 
                      # soil moisture control on nitrogen uptake
                      mn1,
                      mn2, 
                      mn3, 
                      mn4,
                      # temperature control on growth
                      tg1,
                      tg2, 
                      tg3,
                      tg4,
                      # soil moisture control on
                      mg1,
                      mg2,
                      # temperature control on root growth
                      tr1,
                      tr2,
                      # fire related parameters ; not used
                      f1, 
                      f2,
                      # Jinyan adding d15N params
                      n15.rich = 6.93, #d15N of 15N rich source
                      ndepleted = -5.96,#d15N of 15N deplete source
                      k.slope, #mass based fraction of n from rich source
                      ud15N, #d15n observation error
                      
                      fc.700,# co2 response rate
                      # 
                      # // uncertainty parameters
                      uMs, uMr, uCs, uCr,
                      uNs, uNr)
{
  # utilities#####
  # water and t effect; eqn13
  trap1 <- function( x,  a,  b) {
    ret = max(min((x-a)/(b-a),1),0)
    return(ret) 
  }
  # water and t effect;eqn 14
  trap2<- function(x, a, b,c,d){
    
    ret = max( min( c((x-a)/(b-a), 1, (d-x)/(d-c))) , 0 )
    return (ret)
  }
  # control the c to biomass ratio which is used to calculate mass dependent transport, tau, eqn 7
  # not really used as parameter are assumed to be 1
  F_RsC <- function(RHOc, Ms, q){
    ret = RHOc / (Ms^q)
    return(ret)
  }
  # uptake of c and N
  F_P <- function(A0, Ms, KA, Cs, Jc){
    ret = ( A0*Ms ) / ( (1.0 + Ms / KA) * (1.0 + Cs / (Jc*Ms))  )
    return(ret)
  }
  # growth
  F_Gs <- function(gs, Ms, Cs, Ns){
    ret = gs * Ms * Cs/ Ms * Ns/ Ms 
    return (ret)
  }
  # mass affected resistance between root and shoot for n and C
  F_TAUc <- function(Cs, Cr, Ms, Mr, RsC, RrC){
    ret = ( Cs/Ms - Cr/Mr ) / ( RsC + RrC)
    return (ret)
  }
  # biomass dynamics
  F_dMs_dt<- function(Gs, Kl, KM, Ms){
    ret = Gs - ( Kl*Ms ) / ( 1.0 + KM / Ms )
    return (ret)
  }
  # shoot c pool dynamics
  F_dCs_dt<- function(P, Fc, Gs, TAUc){
    
    ret  = P - Fc * Gs - TAUc
    return (ret)
  }
  # root c pool dynamics
  F_dCr_dt<- function(Fc, Gr, TAUc)
  {
    ret  = TAUc - Fc * Gr
    return (ret)
  }
  # plant n uptake changes n sources
  # higher f mean plant take more from 15n rich source
  f.func <- function(x,k=5) 1-x^k #1 / (1 + exp(k*(x - 0.5)))#exp(-k*x) # 1 / (1 + k*x ) # 
  # plant response to co2s
  f_co2 <- function(x,fc.700=1.4){
    #the 1.4 is taken from 
    # 10.1111/j.1365-3040.2007.01641.x
    # fig 3
    # also the same in 3pg
    # https://www.sciencedirect.com/science/article/pii/S0378112797000261
    fC = fc.700 / (2 - fc.700)
    falpha = fC * x / (350* (fC - 1) + x)
    return(falpha)
  }
  # control the c to biomass ratio which is used to calculate mass dependent transport, tau, eqn 7
  # not really used as parameter are assumed to be 1
  F_RsC <- function(RHOc, Ms, q){
    ret = RHOc / (Ms^q)
    return(ret)
  }
  # uptake of c and N
  F_P <- function(A0, Ms, KA, Cs, Jc){
    ret = ( A0*Ms ) / ( (1.0 + Ms / KA) * (1.0 + Cs / (Jc*Ms))  )
    return(ret)
  }
  # growth
  F_Gs <- function(gs, Ms, Cs, Ns){
    ret = gs * (Cs * Ns) / Ms
    return (ret)
  }
  # mass affected resistance between root and shoot for n and C
  F_TAUc <- function(Cs, Cr, Ms, Mr, RsC, RrC){
    ret = ( Cs/Ms - Cr/Mr ) / ( RsC + RrC)
    return (ret)
  }
  # biomass dynamics
  F_dMs_dt<- function(Gs, Kl, KM, Ms){
    ret = Gs - ( Kl*Ms ) / ( 1.0 + KM / Ms )
    return (ret)
  }
  # shoot c pool dynamics
  F_dCs_dt<- function(P, Fc, Gs, TAUc){
    
    ret  = P - Fc * Gs - TAUc
    return (ret)
  }
  # root c pool dynamics
  F_dCr_dt<- function(Fc, Gr, TAUc)
  {
    ret  = TAUc - Fc * Gr
    return (ret)
  }
  # plant n uptake changes n sources
  # higher f mean plant take more from 15n rich source
  f.func <- function(x,k=1) exp(-k*x)#1 / (1 + k*x )
  # plant response to co2s
  f_co2 <- function(x,fc.700=1.4){
    #the 1.4 is taken from 
    # 10.1111/j.1365-3040.2007.01641.x
    # fig 3
    # also the same in 3pg
    # https://www.sciencedirect.com/science/article/pii/S0378112797000261
    fC = fc.700 / (2 - fc.700)
    falpha = fC * x / (350* (fC - 1) + x)
    return(falpha)
  }
  # controls the number of decimal place of state variables
  persision.val <- 1e-6
  #start of model####
  sumb = 0
  # collect output
  cs.vec <- c()
  ns.vec <- c()
  nr.vec <- c()
  ABUTIME <- c()
  n.up.vec <- c()
  plant.d15n.vec <- c()
  f.value <- c()
  n.mass.inhib <- c()
  shoot.biomass <- c()
  photo <- c()
  # setup initial conditions
  Ms_ = initials$ms 
  Mr_ = initials$mr
  d15n <- initials$d15n
  
  Cs_ = Ms_ * Fc
  Cr_ = Mr_ * Fc
  Ns_ = Ms_ * Fn
  Nr_ = Ms_ * Fn
  
  # loop over the time steps
  for(index in c(1:steps)) { 
    # c uptake photosynthesis
    # A0_E = A0* A[index]* trap1(M[index],ma1,ma2) 
    A0_E = A0 * nDay * 
      f_co2(A[index],fc.700 = fc.700) * 
      trap1(M[index],ma1,ma2) 
    # N uptake
    N0_E = N0 * nDay * trap1(TSOIL[index],tn1,tn2) *
      trap2(M[index],mn1,mn2,mn3,mn4)  
    # N0_E_15N <- N0_E*
    # growth factor detemined by t and mass
    GF = trap2(TSOIL[index],tg1,tg2,tg3,tg4)*trap1(M[index],mg1,mg2)             
    # growth rate modiffed
    gs_E = gs * nDay * GF
    gr_E = gr * nDay * GF 
    # size dependent
    RsC = F_RsC(RHOc,Ms_,q)
    RrC = F_RsC(RHOc,Mr_,q)
    RsN = F_RsC(RHOn,Ms_,q)
    RrN = F_RsC(RHOn,Mr_,q)
    # actual uptake rates
    P  = F_P(A0_E,Ms_,KA,Cs_,Jc) * (Ms_ - persision.val) / Ms_
    Un = F_P(N0_E,Mr_,KA,Nr_,Jn) * (Mr_ - persision.val) / Mr_
    Gs = F_Gs(gs = gs_E,Ms = Ms_,Cs = Cs_,Ns = Ns_) * (Ns_ - persision.val) / Ns_
    Gr = F_Gs(gr_E,Mr_,Cr_,Nr_) * (Nr_ - persision.val) / Nr_
    # size dependent resistance
    TAUc = max(F_TAUc(Cs_,Cr_,Ms_,Mr_,RsC,RrC),0) * (Ms_ - persision.val) / Ms_
    TAUn = max(F_TAUc(Nr_,Ns_,Mr_,Ms_,RsN,RrN),0) * (Nr_ - persision.val) / Nr_
    # senecience
    # Kl = A0 *0.1

    LOSS = (Kl*0.5 + Kl*0.5*trap1(TAIR[index],tr1,tr2))*nDay

    # change in pools
    Ms_ = max(persision.val, 
              Ms_ + F_dMs_dt(Gs = Gs,Kl = LOSS,KM = KM,Ms = Ms_) + rnorm(1,0.0, (1e-10)*uMs) )
    Mr_ = max(persision.val, 
              Mr_ + F_dMs_dt(Gr,LOSS,KM,Mr_) + rnorm(1,0.0, (1e-10)*uMr) )
    Cs_ = max(persision.val, 
              Cs_ + F_dCs_dt(P = P,Fc = Fc,Gs = Gs,TAUc = TAUc)   + rnorm(1,0.0, (1e-10)*uCs) )
    Cr_ = max(persision.val, 
              Cr_ + F_dCr_dt(Fc,Gr,TAUc)     + rnorm(1,0.0, (1e-10)*uCr) )
    Ns_ = max(persision.val, 
              Ns_ + F_dCr_dt(Fn,Gs,TAUn)     + rnorm(1,0.0, (1e-10)*uNs) )
    Nr_ = max(persision.val, 
              Nr_ + F_dCs_dt(Un,Fn,Gr,TAUn)  + rnorm(1,0.0, (1e-10)*uNr) )
    
    # 
    # Ms_ = ( 1.0 - trap1(FIRE[index],f1,f2) ) * Ms_ 
    
    cs.vec[index] <- Cs_
    ns.vec[index] <- Ns_
    nr.vec[index] <- Nr_
    n.up.vec[index] <- max(0,Un)
    photo[index] <- P
    # sumb = max(0.0, Ms_ + Mr_)
    ABUTIME[index]=Ms_ + Mr_
    shoot.biomass[index]=Ms_
    #use the mass reduction on uptake here eqn 6 
    if(N0_E>0){
      n.mass.inhib[index] <- Un / N0_E
    }else{
      n.mass.inhib[index] <- 0
    }
    
    # if(N0_E>0){
    f <- max(f.func(n.mass.inhib[index],k=k.slope),0.0001)
      
    # }else{
    #   f <- NA
    # }
    
    #caclualte d15n based on masss wieghting
    f.value[index] <- f
    # n15.rich <- 15
    # ndepleted <- -10
    plant.d15n.new.f <- n15.rich * f + ndepleted *(1-f)
    
    d15n <- (d15n * (Nr_ - n.up.vec[index]) + 
               plant.d15n.new.f * n.up.vec[index]) / Nr_
    
    plant.d15n.vec[index]  <- d15n  + rnorm(1,0.0, (1e-10)*ud15N) 
  }
  return(data.frame(c.s = cs.vec,
                    n.s = ns.vec,
                    n.r = nr.vec,
                    n.up = n.up.vec,
                    plant.biomass = ABUTIME,
                    plant.s.biomass = shoot.biomass,
                    photo = photo,
                    f.value = f.value,
                    n.mass.inhib = n.mass.inhib,
                    d15n = plant.d15n.vec))
}

# plot(sim~c(ABUTIME*parset[21]*1000))
# abline(a=0,b=1)


