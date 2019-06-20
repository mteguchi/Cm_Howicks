# Summarize Howicks data

# Tomo EGuchi
#22 August 2014

rm(list = ls())
runDate <- Sys.Date()

source("~/R/TomosFunctions.R")
D00 <- dirSelector()

fig <- 1
#runDate <- Sys.Date()

fname.base <- 'HW14KAPT_01SEP2014C'
datafile <- paste(D00$Dtomo, 
                  'turtles/Australia/Howicks2014/', 
                  fname.base, '.csv',
                  sep = "")

dat0 <- read.table(file = datafile, sep=",", header = TRUE,
                  na.strings = "")

# extract certain columns 
cols <- c("TAG", "CCL", "WEIGHT", "SEX", 
          "AGECL", "HABIT", "EXPER1", "CDATE",
          "ISRP", "HORMONES","GENETICS", "SP")

dat0_all <- dat0[!is.na(dat0$TAG) & 
                   dat0$SP == "G", cols]

dat0_juve <- dat0[!is.na(dat0$TAG) & 
                    dat0$SP == "G" & 
                    dat0$AGECL == "J", cols]

dat0_subadult <- dat0[!is.na(dat0$TAG) & 
                        dat0$SP == "G" & 
                        (dat0$AGECL == "SA" | 
                           dat0$AGECL == "SP"), cols]

dat0_adult <- dat0[!is.na(dat0$TAG) & 
                    dat0$SP == "G" & 
                    dat0$AGECL == "A", cols]

write.table(dat0_juve, 
            file = paste(D00$Dtomo,
                         'turtles/Australia/Howicks2014/',
                         fname.base, '_juve_', runDate, '.csv',
                         sep = ""),
            sep = ",",
            quote = FALSE,
            na = "NA",
            row.names = FALSE)

dat_gen <- dat0[!is.na(dat0$GENETICS) & 
                  dat0$SP == "G", cols]

dat_hormone_H <- dat0[!is.na(dat0$HORMONES) & 
                      dat0$SP == "H", cols]

dat_hormone <- dat0[!is.na(dat0$HORMONES) & 
                      dat0$SP == "G", cols]

dat_hormone_all <- dat0[!is.na(dat0$HORMONES), cols]

dat_lap <- dat0[!is.na(dat0$EXPER1) & 
                      dat0$SP == "G", cols]

dat_WSR <- dat0[dat0$ISRP == "WSR", cols]

dat_ISR <- dat0[dat0$ISRP == "ISR", cols]

dat_Ei <- dat0[dat0$SP == "H", cols]

# only greens in dat0_ data frames 
uniq_all <- rbind(dat0_juve[dat0_juve$ISRP != 'WSR', ],
                  dat0_subadult[dat0_subadult$ISRP != 'WSR',],
                  dat0_adult[dat0_adult$ISRP != 'WSR', ])

uniq_all$CCL2 <- uniq_all$CCL^2
fit1_CCL_MASS <- lm(log(WEIGHT) ~ CCL, data = uniq_all)

# poly doesn't give the same coefficients... 
#fit2_CCL_MASS <- lm(log(WEIGHT) ~ poly(CCL, degree = 2), data = uniq_all)

fit2_CCL_MASS <- lm(log(WEIGHT) ~ CCL + CCL2, data = uniq_all)
AICs <- AIC(fit1_CCL_MASS, fit2_CCL_MASS)

fit_used <- fit2_CCL_MASS

# AIC indicates the 2nd degree polynomial seems to fit a lot better
CCLseq <- seq(min(uniq_all$CCL, na.rm = TRUE), 
              max(uniq_all$CCL, na.rm = TRUE), 0.5)

# sex ratios:
#cols2 <- c("TAG", "CCL", "WEIGHT", "SEX", "AGECL", "HABIT")
dat_juve <- uniq_all[uniq_all$AGECL == "J", cols]
summary(dat_juve)

# subadult - SA and Subadult pubescent
dat_subadult <- uniq_all[uniq_all$AGECL == "SA" | 
                           uniq_all$AGECL == "SP", 
                         cols]
summary(dat_subadult)

dat_adult <- uniq_all[uniq_all$AGECL == "A" | 
                        uniq_all$AGECL == "AT", 
                      cols]
summary(dat_adult)

if (fig == 1){
  par(mfrow = c(2,1))
  
  # legend box upper left corner
  leg_x <- 90
  leg_y <- 40
  
  # histogram x axis min and max, adjust with observed min and max
  # %% is modulus in R
  byVal <- 5
  minVal <- min(dat0_juve$CCL, na.rm = TRUE) - 
    (min(dat0_juve$CCL, na.rm = TRUE) %% byVal)
  maxVal <- (byVal - max(dat0_adult$CCL, na.rm = TRUE) %% byVal) + 
    max(dat0_adult$CCL,na.rm = TRUE)
  
  breaks <- seq(from = minVal, to = maxVal, by = byVal)
  hist(uniq_all[uniq_all$HABIT == "RF",]$CCL, breaks = breaks, 
       main = paste("CCL-Reef (n=", 
                    dim(uniq_all[uniq_all$HABIT == "RF",])[1], ")", 
                    sep = ""),
       xlab = "",
       ylab = "Frequency",
       bty = "l")
  hist(uniq_all[uniq_all$HABIT == "M", ]$CCL, breaks = breaks, 
       main = paste("CCL-Mangroves (n=", 
                    dim(uniq_all[uniq_all$HABIT == "M", ])[1], ")", 
                    sep = ""), 
       xlab = "CCL (cm)",
       ylab = "Frequency",
       bty = "l")
  
  par(mfrow = c(1,1))
  minVal <- min(uniq_all$CCL, na.rm = TRUE) - 
    (min(uniq_all$CCL, na.rm = TRUE) %% byVal)
  maxVal <- (byVal - max(uniq_all$CCL, na.rm = TRUE) %% byVal) + 
    max(uniq_all$CCL,na.rm = TRUE)
  
  breaks <- seq(from = minVal, to = maxVal, by = byVal)
  hist(uniq_all$CCL, breaks = breaks, 
       main = paste("CCL (n=", 
                    length(na.omit(uniq_all$CCL)), ")", 
                    sep = ""),
       xlab = "CCL (cm)",
       ylab = "Frequency",
       bty = "l")
    
  byVal <- 10
  minVal <- min(uniq_all$WEIGHT, na.rm = TRUE) - 
    (min(uniq_all$WEIGHT, na.rm = TRUE) %% byVal)
  maxVal <- (byVal - max(uniq_all$WEIGHT, na.rm = TRUE) %% byVal) + 
    max(uniq_all$WEIGHT,na.rm = TRUE)
  
  breaks <- seq(from = minVal, to = maxVal, by = byVal)
  hist(uniq_all$WEIGHT, breaks = breaks, 
       main = paste("Mass (n=", 
                    length(na.omit(uniq_all$WEIGHT)), ")", 
                    sep = ""),
       xlab = "Mass (kg)",
       ylab = "Frequency",
       bty = "l")
  
  par(mfrow = c(1,1))
  plot(uniq_all$CCL, 
       uniq_all$WEIGHT, 
       type = "n",
       xlab = "CCL (cm)", ylab = "Mass (kg)",
       bty = "l")
  points(uniq_all[uniq_all$HABIT == "RF",]$CCL, 
         uniq_all[uniq_all$HABIT == "RF",]$WEIGHT, pch = 16,
         col = "red")
  points(uniq_all[uniq_all$HABIT == "M",]$CCL, 
         uniq_all[uniq_all$HABIT == "M",]$WEIGHT, pch = 1,
         col = "blue")
  legend(x=leg_x,y=leg_y,
         legend=c("Reef", "Mangroves"), 
         col = c("red", "blue"),
         pch = c(16, 1))
  
  lines(CCLseq, exp(predict.lm(fit_used, 
                               data.frame(CCL = CCLseq, 
                                          CCL2 = CCLseq^2))),
        col = 'black')
  
  par(mfrow = c(1,1))
  plot(uniq_all$CCL, log(uniq_all$WEIGHT), type = "n",
       xlab = "CCL (cm)", ylab = "ln(Mass)",
       bty = "l")
  points(uniq_all[uniq_all$HABIT == "RF", ]$CCL, 
         log(uniq_all[uniq_all$HABIT == "RF", ]$WEIGHT), 
         pch = 16,
         col = "red")
  points(uniq_all[uniq_all$HABIT == "M", ]$CCL, 
         log(uniq_all[uniq_all$HABIT == "M", ]$WEIGHT), 
         pch = 1,
         col = "blue")
  legend(x=leg_x,y=log(leg_y),
         legend=c("Reef", "Mangroves"), 
         col = c("red", "blue"),
         pch = c(16, 1))
  lines(CCLseq, predict.lm(fit_used, 
                           data.frame(CCL = CCLseq, 
                                      CCL2 = CCLseq^2)),
        col = 'black')
}

