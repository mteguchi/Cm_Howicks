#Closed model 1

# using Mt

rm(list=ls())

library(rjags)
SAVE <- TRUE
runDate <- Sys.Date()
tBegin <- Sys.time()
sysInfo <- Sys.info()

# Change the directory structure according to the computer
ifelse(sysInfo[1] == 'Linux',
       source('~/Documents/R/TomosFunctions.R'),
       source('~/R/TomosFunctions.R'))
D00 <- dirSelector()

## set MCMC parameters
n.adapt <- 50000
n.update <- 10000
n.iter <- 10000
n.chains <- 5

modelName <- paste(D00$Rdir, "Howicks2014/Model_Mt.txt", sep = "")
Ddata <- paste(D00$Dtomo, "turtles/Australia/Howicks2014/", sep = "")
# first get data:
file01 <- paste(Ddata, 'HW14KAPT_01SEP2014A_01_AGECL_J.csv', sep = "")

data01 <- read.table(file = file01, header = FALSE, sep = ",")
#dataID <- read.table(file = fileID, header = FALSE, sep = ",")
if (SAVE == TRUE) saveFname <- paste(D00$Rdir, 
                                     "Howicks2014/RData/Mt_", 
                                     runDate, ".RData", 
                                     sep = "")

# Augment data by X
nz <- 600

yobs <- as.matrix(data01)
yaug <- rbind(yobs, array(0, dim = c(nz, dim(yobs)[2])))

bugs.data <- list(yaug = yaug, M = nrow(yaug), T = ncol(yaug))

inits <- function() list(z = rep(1, nrow(yaug)), 
                         p = runif(ncol(yaug), 0, 1))

params <- c("N", "p", "Omega", "deviance")
params2 <- c("N", "p", "Omega", "deviance", "pD")

jm <- jags.model(modelName,
                 data = bugs.data,
                 inits,
                 n.chains = n.chains,
                 n.adapt = n.adapt)

if (SAVE == TRUE) save(list = ls(all = TRUE), file = saveFname)

update(jm, n.iter = n.update)
if (SAVE == TRUE) save(list = ls(all = TRUE), file = saveFname)

load.module("dic")
zm <- coda.samples(jm, 
                   variable.names = params,
                   n.iter = n.iter)
if (SAVE == TRUE) save(list = ls(all = TRUE), file = saveFname)

dicOut1 <- dic.samples(jm, 
                      n.iter=n.iter,
                      type = "pD")
                      
dicOut2 <- dic.samples(jm, 
                      n.iter=n.iter,
                      type = "popt")

g.diag <- gelman.diag(zm)
h.diag <- heidel.diag(zm)
r.diag <- raftery.diag(zm)

sum_zm <- summary(zm)

meanDev <- sum_zm$statistics[rownames(sum_zm$statistics) == "deviance", 
                             colnames(sum_zm$statistics) == "Mean"]
sdDev <- sum_zm$statistics[rownames(sum_zm$statistics) == "deviance", 
                           colnames(sum_zm$statistics) == "SD"]
DIC <- 0.5*(sdDev^2) + meanDev

tEnd <- Sys.time()

if (SAVE == TRUE) save(list = ls(all = TRUE), file = saveFname)


