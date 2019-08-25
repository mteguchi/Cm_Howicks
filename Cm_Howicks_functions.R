
get.data.clack <- function(filename){
  col.def <- cols(PTAG = col_character(),
                  TAG = col_integer(),
                  TAGST = col_character(),
                  SP = col_character(),
                  SEX = col_character(),
                  AGECL = col_character(),
                  DATE = col_integer(),
                  MONTH = col_integer(),
                  YEAR = col_integer(),
                  LOCAL = col_character(),
                  SECT = col_character(),
                  HABIT = col_character(),
                  PRIAC = col_character(),
                  SECAC = col_character(),
                  TERAC = col_character(),
                  EXPER = col_character(),
                  CCL = col_double())
  
  dat.1 <- read_csv(file = filename, col_types = col.def)
  
  dat.1 %>% mutate(PTAG = as.factor(PTAG),
                   TAGST = as.factor(TAGST),
                   SP = as.factor(SP),
                   SEX = as.factor(SEX),
                   AGECL = as.factor(AGECL),
                   LOCAL = as.factor(LOCAL),
                   HABIT = as.factor(HABIT),
                   YEAR = YEAR + 1900,
                   CDATE = as.Date(paste0(YEAR, "-", MONTH, "-", DATE))) %>%
    filter(AGECL != "A" & AGECL != "SA") %>% 
    transmute(ID = paste(PTAG, TAG, sep = "-"),
              detect = 1,
              TAGST = TAGST,
              SEX = SEX,
              AGECL = AGECL,
              DATE = CDATE,
              LOCAL = LOCAL,
              HABIT = HABIT,
              YEAR = YEAR,
              CCL = CCL)-> dat.1
  
  return(dat.1)
}
get.data <- function(filename){
  col.def <- cols(PTAG = col_character(),
                  TAG = col_integer(),
                  TAGST = col_character(),
                  SP = col_character(),
                  SEX = col_character(),
                  AGECL = col_character(),
                  CDATE = col_date(format = "%m/%d/%Y"),
                  LOCAL = col_character(),
                  YEAR = col_integer(),
                  HABIT = col_character(),
                  PRIAC = col_character(),
                  SECAC = col_character(),
                  CCL = col_double())
  
  dat.1 <- read_csv(file = filename, col_types = col.def)
  
  dat.1 %>% mutate(PTAG = as.factor(PTAG),
                   TAGST = as.factor(TAGST),
                   SP = as.factor(SP),
                   SEX = as.factor(SEX),
                   AGECL = as.factor(AGECL),
                   LOCAL = as.factor(LOCAL),
                   HABIT = as.factor(HABIT)) %>%
    filter(AGECL != "A" & AGECL != "SA") %>% 
    transmute(ID = paste(PTAG, TAG, sep = "-"),
              detect = 1,
              TAGST = TAGST,
              SEX = SEX,
              AGECL = AGECL,
              DATE = CDATE,
              LOCAL = LOCAL,
              HABIT = HABIT,
              YEAR = YEAR,
              CCL = CCL)-> dat.1
  
  return(dat.1)
}


extract.posterior <- function(varname, samples){
  
  all.samples <- do.call(rbind, samples)
  var.names <- dimnames(all.samples)[[2]]
  var.idx <- grep(varname, var.names)
  var.samples <- all.samples[, var.idx]

  out.df <- data.frame(var.samples)
  colnames(out.df) <- dimnames(var.samples)[[2]]
  
  return(out.df)
  
}

dat2CJS <- function(dat.1, save.file = FALSE){
  
  # Create ID by Date and assign 1s
  tmp <-melt(dat.1, 
             id.var = c("ID", "YEAR"), 
             measure.var = "detect")
  
  # make a table with ID by year
  dat.01 <- cast(tmp, ID ~ YEAR)
  
  # replace > 1 with ones
  dat.01 <- as.data.frame(dat.01) %>%
    remove_rownames() %>%
    column_to_rownames(var = "ID")
  
  dat.01[(dat.01 > 1)] <- 1
  
  # save file for later
  if (save.file){
    out.name <- "data/juve_Cm_01_CJS_pt.csv"
    write.csv(dat.01, 
              file = out.name, 
              row.names = T,
              quote = F)
    
  } else {
    out.name <- "not.saved"
  }
  
  out <- list(filename = out.name,
              data = dat.01)
  
  return(out)
  
}


known.state.cjs <- function(ch){
  state <- ch
  for (i in 1:dim(ch)[1]){
    n1 <- min(which(ch[i,]==1))
    n2 <- max(which(ch[i,]==1))
    state[i,n1:n2] <- 1
    state[i,n1] <- NA
  }
  state[state==0] <- NA
  return(state)
}

dat2dat01 <- function(dat.1, year, save.file = FALSE){
  dat.year <- filter(dat.1, YEAR == year)
  
  # Create ID by Date and assign 1s
  tmp <-melt(dat.year, 
             id.var = c("ID", "DATE"), 
             measure.var = "detect")
  
  # make a table with ID by Date
  dat.01.year <- cast(tmp, ID ~ DATE)
  
  # replace NAs with zeros
  dat.01.year[is.na(dat.01.year)] <- 0
  dat.01.year <- as.data.frame(dat.01.year)
  
  # save file for later
  if (save.file){
    out.name = paste0("data/juve_Cm_01_", year, ".csv")
    write.csv(dat.01.year, 
              file = out.name, 
              row.names = F,
              quote = F)
    
  } else {
    out.name <- "not.saved"
  }
  
  out <- list(filename = out.name,
              data = dat.01.year)
  
  return(out)
  
}


dat2data01 <- function(dat.1, save.file = FALSE){
  
  # Create ID by Date and assign 1s
  tmp <-melt(dat.1, 
             id.var = c("ID", "DATE"), 
             measure.var = "detect")
  
  # make a table with ID by year
  dat.01 <- cast(tmp, ID ~ DATE)
  
  # replace > 1 with ones
  dat.01 <- as.data.frame(dat.01) %>%
    remove_rownames() %>%
    column_to_rownames(var = "ID")
  
  dat.01[(dat.01 > 1)] <- 1
  
  # save file for later
  if (save.file){
    out.name <- "data/juve_Cm_01_all.csv"
    write.csv(dat.01, 
              file = out.name, 
              row.names = T,
              quote = F)
    
  } else {
    out.name <- "not.saved"
  }
  
  out <- list(filename = out.name,
              data = dat.01)
  
  return(out)
  
}


jags.M0 <- function(data01, nz = 200, MCMC.params){
  inits <- function() list(z = rep(1, nrow(yaug)), 
                           p = runif(1, 0, 1))
  
  jags.parameters <- c("N", "p", "Omega", "deviance")#, "log.likelihood")
  
  # Augment data by X
  nz <- nz
  
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
             model.file = MCMC.params$model.file,
             n.chains = MCMC.params$n.chains,
             n.burnin = MCMC.params$n.burnin,
             n.thin = MCMC.params$n.thin,
             n.iter = MCMC.params$n.samples,
             DIC = T, 
             parallel=T)
  
  t.End <- Sys.time()
  
  out.list <- list(jags.out = jm,
                   jags.data = jags.data,
                   elapsed.time = t.Begin - t.End)
  
  return(out.list)
}

jags.Mt <- function(data01, nz = 200, MCMC.params){
  
  inits <- function() list(z = rep(1, nrow(yaug)), 
                           p = runif(ncol(yaug), 0, 1))
  
  jags.parameters <- c("N", "p", "Omega", "deviance")#, "log.likelihood")
  
  # Augment data by X
  nz <- nz
  
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
             model.file = MCMC.params$model.file,
             n.chains = MCMC.params$n.chains,
             n.burnin = MCMC.params$n.burnin,
             n.thin = MCMC.params$n.thin,
             n.iter = MCMC.params$n.samples,
             DIC = T, 
             parallel=T)
  
  t.End <- Sys.time()
  
  out.list <- list(jags.out = jm,
                   jags.data = jags.data,
                   elapsed.time = t.Begin - t.End)
  
  return(out.list)
}

# extract first CCL measurement for each individual
dat2CCL <- function(dat.1, save.file = FALSE){
  dat.1 %>% select(ID, CCL, DATE) %>%
    group_by(ID) %>%
    summarise(CCL_min = min(CCL, na.rm = T),
              CCL_mean = mean(CCL, na.rm = T),
              CCL_max = max(CCL, na.rm = T),
              Date1 = min(DATE)) -> dat.CCL1
  
  # save file for later
  if (save.file){
    out.name <- "data/juve_Cm_CCL.csv"
    write.csv(dat.CCL1, 
              file = out.name, 
              row.names = F,
              quote = F)
    
  } else {
    out.name <- "not.saved"
  }
  
  out <- list(filename = out.name,
              data = dat.CCL1)
  
  return(out)
}

get.first <- function(x) min(which(x != 0))
