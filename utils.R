##################################################
# File containing functions to use in analyses
##################################################

##################################################
# Function to load data, merge, and return a dataset for analyses
##################################################

load_data <- function(path) {
  
  # install.packages("data.table")
  library(data.table)

  # get data
  total1 = fread(file = paste0(path, 'hesin_registry_assess_cent_all_v2.csv'), sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  
  # polygenic risk score
  prs = fread(file = paste0(path, 'prs_cad_ukbb.tsv'), sep = '\t', header = TRUE, stringsAsFactors = FALSE, na.strings = "")
  
  # SNP matrix of genome-wide significant SNPs for MI
  snp_mat <- fread(file = paste0(path, 'variant_gw_cad_allT.tsv'), sep = '\t', header = TRUE, stringsAsFactors = FALSE)

  # merge phenotype data, polygenic risk score, and SNP matrix
  totalo = merge(total1, prs, by.x = "eid", by.y = "ID", all.x = TRUE)
  total = merge(totalo, snp_mat, by = "eid", all.x = TRUE)
  
  return(data.frame(total))

}

##################################################
# Function to plot average polygenic risk score for MI, no MI, ICD and no ICD
##################################################
        
plot_icd_mi_prs <- function(icdcode, alldata, allmi,grouptype="MI") {
  
  alldataU <- alldata[!duplicated(alldata$eid), c("eid", "sex", "age", "death", "death_date", "endfollowup", "startfollowup","PRS_0.5")]
  
  alldataT <- alldata[alldata$icd10 == icdcode & !is.na(alldata$icd10), ]
  alldata.sorted <- alldataT[order(alldataT$eid, as.numeric(alldataT$icd10_date)), ]
  alldata.sortedU <- alldata.sorted[!duplicated(alldata.sorted$eid), c("eid", "icd10_date")]
  
  temp1 <- merge(alldataU, alldata.sortedU, all.x = T, by = "eid")
  
  total.comorb <- merge(temp1, allmi, all.x = T, by = "eid")
  
  total.comorb$pred <- ifelse(is.na(total.comorb$icd10_date), 0, 1)
  total.comorb$mi <- ifelse(is.na(total.comorb$mi_date), 0, 1)
  
  icd.mi <- total.comorb[which(total.comorb$pred == 1 & total.comorb$mi == 1),]
  icd.no.mi <- total.comorb[which(total.comorb$pred == 1 & total.comorb$mi == 0),]
  no.icd.mi <- total.comorb[which(total.comorb$pred == 0 & total.comorb$mi == 1),]
  no.icd.no.mi <- total.comorb[which(total.comorb$pred == 0 & total.comorb$mi == 0),]
  
  df <- data.frame(prs_mean = c(mean(icd.mi$PRS_0.5, na.rm = T), mean(icd.no.mi$PRS_0.5, na.rm = T), mean(no.icd.mi$PRS_0.5, na.rm = T), mean(no.icd.no.mi$PRS_0.5, na.rm = T)), prs_se = c(sd(icd.mi$PRS_0.5, na.rm = T) / sqrt(nrow(icd.mi)), sd(icd.no.mi$PRS_0.5, na.rm = T) / sqrt(nrow(icd.no.mi)), sd(no.icd.mi$PRS_0.5, na.rm = T) / sqrt(nrow(no.icd.mi)), sd(no.icd.no.mi$PRS_0.5, na.rm = T) / sqrt(nrow(no.icd.no.mi))), mi = c("MI", "NO MI", "MI", "NO MI"), ICD = c(icdcode, icdcode, paste(".NO", icdcode), paste(".NO", icdcode)))
  df$CIU <- df$prs_mean + 1.96*df$prs_se
  df$CID <- df$prs_mean - 1.96*df$prs_se
  
  if (grouptype == "ICD") {
    ggplot(aes(y = prs_mean, x = mi), data = df) + geom_bar(aes(fill = ICD), stat = "identity", position = "dodge") + geom_errorbar(width = 0.2, aes(ymin = CID, ymax = CIU, fill = ICD), position = position_dodge(.9)) + theme_bw() + ylab("Mean Polygenic Score for MI")
  }
  else if (grouptype == "MI") {
    ggplot(aes(y = prs_mean, x = ICD), data = df) + geom_bar(aes(fill = mi), stat = "identity", position = "dodge") + geom_errorbar(width = 0.2, aes(ymin = CID, ymax = CIU, fill = mi), position = position_dodge(.9)) + theme_bw() + ylab("Mean Polygenic Score for MI")
  }

}

##################################################
# Function to plot association between ICD code and MI as a function of polygenic risk score
##################################################

plot_icd_mi_prs_int <- function(icdcode, alldata, allmi) {
  
  alldataU <- alldata[!duplicated(alldata$eid), c("eid", "sex", "age", "death", "death_date", "endfollowup", "startfollowup","PRS_0.5")]
  
  alldataT <- alldata[alldata$icd10 == icdcode & !is.na(alldata$icd10), ]
  alldata.sorted <- alldataT[order(alldataT$eid, as.numeric(alldataT$icd10_date)), ]
  alldata.sortedU <- alldata.sorted[!duplicated(alldata.sorted$eid), c("eid", "icd10_date")]
  temp1 <- merge(alldataU, alldata.sortedU, all.x = T, by = "eid")
  total.comorb <- merge(temp1, allmi, all.x = T, by = "eid")
  total.comorb$pred <- ifelse(is.na(total.comorb$icd10_date), 0, 1)
  total.comorb$mi <- ifelse(is.na(total.comorb$mi_date), 0, 1)
  
  data1 <- subset(total.comorb, pred == 1)
  data1$tstart <- 0
  data1$tstop <- as.numeric(data1$icd10_date - data1$startfollowup)
  data1$outP <- 0
  data1$trt <- 0
  
  data2 <- subset(total.comorb, pred == 1)
  data2$tstart <- as.numeric(data2$icd10_date - data2$startfollowup)
  data2$tstop <- as.numeric(data2$endfollowup - data2$startfollowup)
  data2$outP <- data2$mi
  data2$trt <- 1
  
  data3 <- subset(total.comorb, pred == 0 & mi == 1)
  data3$tstart <- 0
  data3$tstop <- as.numeric(data3$endfollowup - data3$startfollowup)
  data3$outP <- 1
  data3$trt <- 0
  
  dataF <- rbind(data1, data2, data3)
  
  otherdata <- subset(total.comorb, pred != 1 & mi != 1)
  
  otherdata$tstart <- 0
  otherdata$tstop <- as.numeric(otherdata$endfollowup - otherdata$startfollowup)
  otherdata$outP <- 0
  otherdata$trt <- 0
  
  findata <- rbind(dataF, otherdata)
  
  findata <- findata[findata$tstart < findata$tstop, ]
  findata$PRS_quartiles <- as.numeric(cut(findata$PRS_0.5, quantile(findata$PRS_0.5, seq(0, 1, 0.25), na.rm = T)))
  
  mod1 <- coxph(Surv(tstart, tstop, outP) ~ age + trt + sex, data = findata[findata$PRS_quartiles == 1, ])
  mod2 <- coxph(Surv(tstart, tstop, outP) ~ age + trt + sex, data = findata[findata$PRS_quartiles == 2, ])
  mod3 <- coxph(Surv(tstart, tstop, outP) ~ age + trt + sex, data = findata[findata$PRS_quartiles == 3, ])
  mod4 <- coxph(Surv(tstart, tstop, outP) ~ age + trt + sex, data = findata[findata$PRS_quartiles == 4, ])
  
  df <- data.frame(BETA = c(coef(mod1)[2], coef(mod2)[2], coef(mod3)[2], coef(mod4)[2]),
                   SE = c(coef(summary(mod1))[2, 3], coef(summary(mod2))[2, 3], coef(summary(mod3))[2, 3], coef(summary(mod4))[2, 3]),
                   P = c(coef(summary(mod1))[2, 5], coef(summary(mod2))[2, 5], coef(summary(mod3))[2, 5], coef(summary(mod4))[2, 5]),
                   X = c("1st", "2nd", "3rd", "4th"))
  
  df$CIU <- exp(df$BETA+1.96*df$SE)
  df$CID <- exp(df$BETA-1.96*df$SE)
  
  ggplot(aes(y = exp(BETA), x = X), data = df) + geom_errorbar(width = .1, aes(ymin = CID, ymax = CIU), colour="red") + geom_point(shape = 21, size = 3, fill = "white") + theme_bw() + ylab(paste("HR for Association Between", icdcode, "and MI")) + xlab("Quartiles of PRS for MI")
  
}

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}