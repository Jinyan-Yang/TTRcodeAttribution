# paper
# https://www.nature.com/articles/s41561-022-01114-x
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
# acutal function#####
ThornTimeARU <- function( steps, 
                          initials, 
                          TAIR, 
                          TSOIL, 
                          M, 
                          N, 
                          FIRE, 
                          A,
                          Kl, 
                          gs, 
                          gr, 
                          KM, 
                          A0, 
                          N0,
                          KA, 
                          Jc, 
                          Jn, 
                          q, 
                          RHOc, 
                          RHOn,
                          Fc, 
                          Fn, 
                          # //  to be fitted environmental dependencies
                          ma1,
                          ma2,
                          tn1,
                          tn2, 
                          mn1,
                          mn2, 
                          mn3, 
                          mn4,
                          tg1,
                          tg2, 
                          tg3,
                          tg4,
                          mg1,
                          mg2,
                          tr1,
                          tr2,
                          f1, 
                          f2,
                          # // uncertainty parameters
                          uMs, uMr, uCs, uCr,
                          uNs, uNr)
{

  sumb=0
  # collect output
  cs.vec <- c()
  ns.vec <- c()
  ABUTIME <- c()
  
  # setup initial conditions
  Ms_ = initials 
  Mr_ = initials 
  Cs_ = Ms_ * 0.05
  Cr_ = Mr_ * 0.05
  Ns_ = Ms_ * 0.01
  Nr_ = Ms_ * 0.01
  
  # loop over the time steps
  for(index in 1:steps) { 
    # c uptake photosynthesis
    A0_E = A0* A[index]*  trap1(M[index],ma1,ma2) 
    # N uptake
    N0_E = N0* trap1(TSOIL[index],tn1,tn2) *
      trap2(M[index],mn1,mn2,mn3,mn4)  
    # N0_E_15N <- N0_E*
    # growth factor detemined by t and mass
    GF = trap2(TSOIL[index],tg1,tg2,tg3,tg4)*trap1(M[index],mg1,mg2)             
    # growth rate modiffed
    gs_E = gs* GF
    gr_E = gr* GF 
    # size dependent
    RsC = F_RsC(RHOc,Ms_,q)
    RrC = F_RsC(RHOc,Mr_,q)
    RsN = F_RsC(RHOn,Ms_,q)
    RrN = F_RsC(RHOn,Mr_,q)
    # actual uptake rates
    P  = F_P(A0_E,Ms_,KA,Cs_,Jc)
    Un = F_P(N0_E,Mr_,KA,Nr_,Jn)
    Gs = F_Gs(gs_E,Ms_,Cs_,Ns_)
    Gr = F_Gs(gr_E,Mr_,Cr_,Nr_)
    # size dependent resistence
    TAUc = F_TAUc(Cs_,Cr_,Ms_,Mr_,RsC,RrC)
    TAUn = F_TAUc(Nr_,Ns_,Mr_,Ms_,RsN,RrN)
    # senecience
    LOSS = Kl*0.5 + Kl*0.5*trap1(TAIR[index],tr1,tr2)
    # change in pools
    Ms_ = max(0.0, Ms_ + F_dMs_dt(Gs,LOSS,KM,Ms_) + rnorm(1,0.0, uMs) )
    Mr_ = max(0.0, Mr_ + F_dMs_dt(Gr,LOSS,KM,Mr_) + rnorm(1,0.0, uMr) )
    Cs_ = max(0.0, Cs_ + F_dCs_dt(P,Fc,Gs,TAUc)   + rnorm(1,0.0, uCs) )
    Cr_ = max(0.0, Cr_ + F_dCr_dt(Fc,Gr,TAUc)     + rnorm(1,0.0, uCr) )
    Ns_ = max(0.0, Ns_ + F_dCr_dt(Fn,Gs,TAUn)     + rnorm(1,0.0, uNs) )
    Nr_ = max(0.0, Nr_ + F_dCs_dt(Un,Fn,Gr,TAUn)  + rnorm(1,0.0, uNr) )
    
    Ms_ = ( 1.0 - trap1(FIRE[index],f1,f2) ) * Ms_ 
    
    cs.vec[index] <- Cs_
    ns.vec[index] <- Ns_
    
    sumb = max(0.0, Ms_ + Mr_)
    ABUTIME[index]=sumb
   
    # d15nplant = 
  }
  return(data.frame(cs = cs.vec,
                    ns = ns.vec,
                    ABUTIME = ABUTIME))
}

# plot(sim~c(ABUTIME*parset[21]*1000))
# abline(a=0,b=1)
