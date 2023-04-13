source('r/readERA5.r')
source('r/TTR_d15N.R')
# inputs#####
# read rf model fit fr dn15
# rf.fit <- readRDS('\\\\fs1-cbr.nexus.csiro.au/{mmrg}/work/users/yan190/repo/delta_n_15/cache/rf.kFold.n15.rds')
# read landsat
inputs.ls <- readRDS(file="MyData_SA-BFA-BNP.rds")

# landsat.ts.ls <- readRDS('//fs1-cbr.nexus.csiro.au/{mmrg}/work/users/yan190/repo/delta_n_15/cache/ls.ts.0.1.part2.rds')

library(randomForest)
library(caret)
source('//fs1-cbr.nexus.csiro.au/{mmrg}/work/users/yan190/repo/delta_n_15/r/functions_json.R')

# landsat.ls <- readRDS('')
fit.rf.n15 <- readRDS('//fs1-cbr.nexus.csiro.au/{mmrg}/work/users/yan190/repo/delta_n_15/cache/rf.kFold.n15.rds')
tmp.df <- read.csv('siteSampleAU.csv')
tmp.ls <- apply(tmp.df,1, get.TS.func,
                lat.col = 3,lon.col=2,json.col=4,
                date.col=100,n15.col=100)

landsat.n15.ls <- lapply(tmp.ls, get.dn154ts.new.func)
# landsat.n15.ls[[1]]
saveRDS(landsat.n15.ls,'cache/sampleAU.ls.rds')
