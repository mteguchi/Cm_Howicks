#Closed model 1

# using M0

rm(list=ls())

library(jagsUI)
library(bayesplot)

SAVE <- TRUE

model.file <- "models/Model_M0.txt"

MCMC.n.chains <- 5
MCMC.n.samples <- 5000
MCMC.n.burnin <- 2000
MCMC.n.thin <- 5

# first get data:
file01 <- 'data/HW14KAPT_01SEP2014A_01_AGECL_J.csv'
#fileID <- paste(Ddata, 'Resights01ID_04102014.csv', sep = "")
data01 <- read.table(file = file01, header = FALSE, sep = ",")
#dataID <- read.table(file = fileID, header = FALSE, sep = ",")

if (SAVE == TRUE) saveFname <- paste("RData/M0_", 
                                     Sys.Date(), ".rds", 
                                     sep = "")

inits <- function() list(z = rep(1, nrow(yaug)), p = runif(1, 0, 1))

jags.parameters <- c("N", "p", "Omega", "deviance")#, "log.likelihood")

# Augment data by X
nz <- 100

yobs <- as.matrix(data01)
yaug <- rbind(yobs, 
              array(0, dim = c(nz, dim(yobs)[2])))
jags.data <- list(yaug = yaug, 
                  M = nrow(yaug), 
                  T = ncol(yaug))

t.Begin <- Sys.time()

# initialize the model with all other stuff
jm <- jags(data = jags.data,
           inits = inits,
           parameters.to.save= jags.parameters,
           model.file = model.file,
           n.chains = MCMC.n.chains,
           n.burnin = MCMC.n.burnin,
           n.thin = MCMC.n.thin,
           n.iter = MCMC.n.samples,
           DIC = T, 
           parallel=T)

t.End <- Sys.time()

# look at trace and posteriors
mcmc_trace(jm$samples, c("N", "p"))

mcmc_dens(jm$samples, c("N", "p"))

out.list <- list(jm = jm,
                 data = jags.data,
                 MCMC.params = list(MCMC.n.chains = MCMC.n.chains,
                                    MCMC.n.samples = MCMC.n.samples,
                                    MCMC.n.burnin = MCMC.n.burnin,
                                    MCMC.n.thin = MCMC.n.thin),
                 model.file = model.file,
                 elapsed.time = t.End - t.Begin)

if (SAVE == TRUE) saveRDS(out.list, file = saveFname)



