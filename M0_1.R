#Closed model 1

# using M0

rm(list=ls())

library(rjags)
SAVE <- F
runDate <- Sys.Date()
tBegin <- Sys.time()
sysInfo <- Sys.info()

# Change the directory structure according to the computer
ifelse(sysInfo[1] == 'Linux',
       source('~/Documents/R/tools/TomosFunctions.R'),
       source('~/R/tools/TomosFunctions.R'))
D00 <- dirSelector()

## set MCMC parameters
n.adapt <- 5000
n.update <- 1000
n.iter <- 1000
n.chains <- 5

modelName <- "models/Model_M0.txt"

# first get data:
file01 <- 'data/HW14KAPT_01SEP2014A_01_AGECL_J.csv'
#fileID <- paste(Ddata, 'Resights01ID_04102014.csv', sep = "")
data01 <- read.table(file = file01, header = FALSE, sep = ",")
#dataID <- read.table(file = fileID, header = FALSE, sep = ",")

if (SAVE == TRUE) saveFname <- paste("RData/M0_", 
                                     runDate, ".RData", 
                                     sep = "")


# Augment data by X
nz <- 900

yobs <- as.matrix(data01)
yaug <- rbind(yobs, array(0, dim = c(nz, dim(yobs)[2])))

bugs.data <- list(yaug = yaug, M = nrow(yaug), T = ncol(yaug))

inits <- function() list(z = rep(1, nrow(yaug)), p = runif(1, 0, 1))

params <- c("N", "p", "Omega", "deviance")

jm <- jags.model(modelName,
                 data = bugs.data,
                 inits,
                 n.chains = n.chains,
                 n.adapt = n.adapt)

if (SAVE == TRUE) save(list = ls(all = TRUE), file = saveFname)

update(jm, n.iter = n.update)
if (SAVE == TRUE) save(list = ls(all = TRUE), file = saveFname)

load.module("dic")
zm <- coda.samples(jm, variable.names = params,
                  n.iter = n.iter)
if (SAVE == TRUE) save(list = ls(all = TRUE), file = saveFname)

dicOut <- dic.samples(jm, 
                      n.iter = n.iter,
                      type = "pD")
dicOut2 <- dic.samples(jm, 
                      n.iter = n.iter,
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

