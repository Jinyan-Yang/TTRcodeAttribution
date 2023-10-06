source('r/readERA5.r')
source('r/TTR_d15N.R')
# inputs#####
# read rf model fit fr dn15
# rf.fit <- readRDS('\\\\fs1-cbr.nexus.csiro.au/{mmrg}/work/users/yan190/repo/delta_n_15/cache/rf.kFold.n15.rds')
# read landsat
inputs.ls <- readRDS(file="MyData_SA-BFA-BNP.rds")
# inputs.ls$X
# landsat.ts.ls <- readRDS('//fs1-cbr.nexus.csiro.au/{mmrg}/work/users/yan190/repo/delta_n_15/cache/ls.ts.0.1.part2.rds')


# landsat.ls <- readRDS('')

if(file.exists('cache/sampleAU_havPlot.ls.rds')){
  landsat.n15.ls <- readRDS('cache/sampleAU_havPlot.ls.rds')
}else{
  library(randomForest)
  library(caret)
  source('//fs1-cbr.nexus.csiro.au/{mmrg}/work/users/yan190/repo/delta_n_15/r/functions_json.R')
  # 
  # fit.rf.n15 <- readRDS('//fs1-cbr.nexus.csiro.au/{mmrg}/work/users/yan190/repo/delta_n_15/cache/rf.kFold.n15.rds')
  # tmp.df <- read.csv('siteSampleAU.csv')
  # tmp.ls <- apply(tmp.df,1, get.TS.func,
  #                 lat.col = 3,lon.col=2,json.col=4,
  #                 date.col=100,n15.col=100)
  # landsat.n15.ls <- lapply(tmp.ls, get.dn154ts.new.func)
  # # landsat.n15.ls[[1]]
  # saveRDS(landsat.n15.ls,'cache/sampleAU.ls.rds')
  fit.rf.n15 <- readRDS('//fs1-cbr.nexus.csiro.au/{mmrg}/work/users/yan190/repo/delta_n_15/cache/rf.kFold.n15.rds')
  tmp.df <- read.csv('cache/siteSampleAU_havPlot.csv')
  tmp.ls <- apply(tmp.df,1, get.TS.func,
                  lat.col = 3,lon.col=2,json.col=18,
                  date.col=15,n15.col=4)
  landsat.n15.ls <- lapply(tmp.ls, get.dn154ts.new.func)
  # landsat.n15.ls[[1]]
  saveRDS(landsat.n15.ls,'cache/sampleAU_havPlot.ls.rds')
  
}


# landsat.ls <- readRDS('')
#########
if(file.exists('cache/sampleAU_60k.rds')){
  landsat.n15.ls <- readRDS('cache/sampleAU_60k.rds')
}else{
  library(randomForest)
  library(caret)
  source('//fs1-cbr.nexus.csiro.au/{mmrg}/work/users/yan190/repo/delta_n_15/r/functions_json.R')
  # 
  # fit.rf.n15 <- readRDS('//fs1-cbr.nexus.csiro.au/{mmrg}/work/users/yan190/repo/delta_n_15/cache/rf.kFold.n15.rds')
  # tmp.df <- read.csv('siteSampleAU.csv')
  # tmp.ls <- apply(tmp.df,1, get.TS.func,
  #                 lat.col = 3,lon.col=2,json.col=4,
  #                 date.col=100,n15.col=100)
  # landsat.n15.ls <- lapply(tmp.ls, get.dn154ts.new.func)
  # # landsat.n15.ls[[1]]
  # saveRDS(landsat.n15.ls,'cache/sampleAU.ls.rds')
  fit.rf.n15 <- readRDS('//fs1-cbr.nexus.csiro.au/{mmrg}/work/users/yan190/repo/delta_n_15/cache/rf.kFold.n15.rds')
  tmp.df <- read.csv('cache/auByClass.csv')
  tmp.ls <- apply(tmp.df,1, 
                  get.TS.func,
                  lat.col = 3,
                  lon.col = 2,
                  json.col=5,
                  date.col=4,
                  n15.col=10)
  landsat.n15.ls <- lapply(tmp.ls, 
                           get.dn154ts.new.func)
  # landsat.n15.ls[[1]]
  saveRDS(landsat.n15.ls,'cache/sampleAU_60k.rds')
  
}
