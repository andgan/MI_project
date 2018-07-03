##################################################
# Goal: Find association between MI and ICD and check how association relates to PRS
# Sensitivity analysis: Keep only ICD10 codes after baseline
##################################################

# load libraries 
# install.packages("survival")
# install.packages("ggplot2")
# install.packages("ggrepel")
# install.packages("reshape2")
# install.packages("caret")
# install.packages("pROC")
# install.packages("plotROC")
library(survival)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(caret)
library(pROC)
library(plotROC)

# identify data path
path <- "C:/Jiwoo Lee/Myocardial Infarction Research Project 2017/"
path <- "/Users/andreaganna/Documents/Work/Post_doc/jiwoo/"

source(paste0(path,"MI_project/utils.R"))

# load data with function from utils.R
total <- load_data(path)

##################################################
# Goal: Clean data
##################################################

# change mi_date from string to date
total$mi_date <- as.Date(total$mi_date, "%Y-%m-%d")

# change icd_date from string to date
total$icd10_date <- as.Date(total$icd10_date, "%Y-%m-%d")

# create end-of-follow-up analysis, which is first occurence of 2015-03-31 or MI or death
total$endfollowup <- ifelse(!is.na(total$death_date) & is.na(total$mi_date), as.Date(total$death_date, "%Y-%m-%d"), ifelse(!is.na(total$mi_date), as.Date(total$mi_date, "%Y-%m-%d"), as.Date("2015-03-01", "%Y-%m-%d")))
class(total$endfollowup) <- "Date"

# create start-of-follow-up, which is entrance in the study (for sensitivity analysis)
total$startfollowup <- as.Date(total$ac_date, "%Y-%m-%d")
nrow(total)
#2864652 rows
length(unique(total$eid))
#498828 individuals
length(unique(total$icd10))
#1472 ICD codes
length(unique(total[which(!is.na(total$mi_date)),]$eid))
#17959 individuals with MI

# set ICD10 codes and dates to NA if they are not related to disease (e.g., car crash) 
total.new = total
total.new$icd10[grepl("^V|^Z", total.new$icd10)] <- NA
total.new$icd10_date[grepl("^V|^Z", total.new$icd10)] <- NA
nrow(total.new)
#2864652 rows
length(unique(total.new$eid))
#498828 individuals
length(unique(total.new$icd10))
#1393 ICD codes
length(unique(total.new[which(!is.na(total.new$mi_date)),]$eid))
#17959 individuals with MI

# remove individuals with MI outside of timeframe
to_remove <- total.new$eid[total.new$mi_date > as.Date("2015-03-01", "%Y-%m-%d") | total.new$mi_date < total.new$startfollowup]
to_remove <- unique(to_remove[!is.na(to_remove)])
total.new = total.new[!total.new$eid %in% to_remove,]
nrow(total.new)
#2681645 rows
length(unique(total.new$eid))
#486626 individuals
length(unique(total.new$icd10))
#1389 ICD codes
length(unique(total.new[which(!is.na(total.new$mi_date)),]$eid))
#5757 individuals with MI

# set ICD10 codes and dates to NA if they are outside of timeframe or less than 7 days from MI
total.new$icd10[total.new$icd10_date > as.Date("2015-03-01", "%Y-%m-%d") | total.new$icd10_date < total.new$startfollowup] <- NA
total.new$icd10_date[total.new$icd10_date > as.Date("2015-03-01", "%Y-%m-%d") | total.new$icd10_date < total.new$startfollowup] <- NA
# this is different between pre- and post-analyses
total.new$icd10[(total.new$mi_date - total.new$icd10_date) < 8] <- NA
total.new$icd10_date[(total.new$mi_date - total.new$icd10_date) < 8] <- NA
nrow(total.new)
#2681645 rows
length(unique(total.new$eid))
#486626 individuals
length(unique(total.new$icd10))
#1274 ICD codes
length(unique(total.new[which(!is.na(total.new$mi_date)),]$eid))
#5757 individuals with MI

# remove individuals with I25 before MI
# LET'S CHECK WITH ERIK IF THIS EXCLUSION IS NEED IT
to_remove <- total.new$eid[total.new$icd10 == "I25"]
to_remove <- unique(to_remove[!is.na(to_remove)])
total.new = total.new[!total.new$eid %in% to_remove,]
nrow(total.new)
#2590659 rows
length(unique(total.new$eid))
#479990 individuals
length(unique(total.new$icd10))
#1268 ICD codes
length(unique(total.new[which(!is.na(total.new$mi_date)),]$eid))
#5245 individuals with MI

# keep only unique individuals with MI
total.newT <- total.new[!is.na(total.new$mi_date),]
total.newMI <- total.newT[!duplicated(total.newT$eid), c("eid","mi_date")]

# keep only unique individuals and a subset of variables from total.new
total.newU <- total.new[!duplicated(total.new$eid), c("eid", "sex", "age", "death", "death_date", "startfollowup", "endfollowup", "PRS_0.5", "sbp", "smoke", "bmi", "diabetes", "lipid_lowering")]

##################################################
# Goal: Examine PRS
##################################################

# keep only unique individuals and a subset of variables from total.new
total.new.prs <- total.new[!duplicated(total.new$eid), c("eid", "sex", "age", "PRS_0.5", "mi_date")]
total.new.prs$mi <- ifelse(is.na(total.new.prs$mi_date), 0, 1)
total.new.prs$mi <- factor(total.new.prs$mi)
# plot frequency of PRS stratified for MI
ggplot(total.new.prs, aes(x = PRS_0.5, fill = mi)) + 
  geom_histogram(binwidth = 1) + 
  labs(x = "Polygenic Risk Score", y = "Frequency") + 
  scale_fill_discrete(name = "", labels = c("No MI", "MI"))
# plot average PRS for MI and no MI
prs.sum <- summarySE(total.new.prs[!is.na(total.new.prs$PRS_0.5),], measurevar = "PRS_0.5", groupvars = "mi")
prs.sum2 <- prs.sum
prs.sum2$mi <- factor(prs.sum2$mi)
ggplot(prs.sum2, aes(x = mi, y = PRS_0.5)) + 
  geom_bar(position = position_dodge(), stat = "identity") +
  geom_errorbar(aes(ymin = PRS_0.5-ci, ymax = PRS_0.5+ci), width = 0.2, position = position_dodge(1)) +
  xlab("MI Status") +
  ylab("Average Polygenic Risk Score")
# plot prevalence of MI for different bins of PRS
total.new.prs$bin <- round(total.new.prs$PRS_0.5, digits = 0)
mi.prev <- NULL
mi.prev <- data.frame(mi.prev)
for (i in seq(-38, 38, 3)) {
  mi.prev.temp <- length(which(total.new.prs$bin == i & total.new.prs$mi == 1)) / length(which(total.new.prs$bin == i))
  mi.prev <- rbind(mi.prev, mi.prev.temp)
  print(i)
}
mi.prev$bin <- seq(-38, 38, 3)
colnames(mi.prev) <- c("prev", "bin")
ggplot(mi.prev, aes(x = bin, y = prev)) +
  geom_point() +
  xlim(-38, 38) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))  

##################################################
# Goal: Calculate statistics for all ICD10 codes in a forloop
##################################################

# loop across ICD10 codes associated with MI and n > 10
icd.mi <- total.new[!is.na(total.new$mi_date), c("eid", "icd10")]
icd.mi <- unique(icd.mi[, c('eid','icd10')])
icd.mi <- icd.mi[!is.na(icd.mi$icd10), "icd10"]
icd.all <- names(table(icd.mi))[table(icd.mi) > 10]
res <- NULL
icd.final <- NULL

for (i in 1:length(icd.all)) {
  
  # keep unique record for comorbity examined
  total.newT <- total.new[total.new$icd10 == icd.all[i] & !is.na(total.new$icd10),]
  
  # sort and keep first event for each individual
  total.new.sorted <- total.newT[order(total.newT$eid, as.numeric(total.newT$icd10_date)),]
  total.new.sortedU <- total.new.sorted[!duplicated(total.new.sorted$eid), c("eid", "icd10_date")]
  
  # merge ICD10-only dataset with  main dataset
  temp1 <- merge(total.newU, total.new.sortedU, all.x = T, by = "eid")
  
  # add MI-only dataset
  total.comorb <- merge(temp1, total.newMI, all.x = T, by = "eid")
  
  # assign pred and MI
  total.comorb$pred <- ifelse(is.na(total.comorb$icd10_date), 0, 1)
  total.comorb$mi <- ifelse(is.na(total.comorb$mi_date), 0, 1)
  
  # calculate mean time distance between ICD10 code and MI
  total.comorbT <- total.comorb[total.comorb$pred == 1 & total.comorb$mi == 1,]
  mean_day_dist <- mean(total.comorbT$mi_date - total.comorbT$icd10_date)
  sd_day_dist <- sd(total.comorbT$mi_date - total.comorbT$icd10_date)
  
  # if you have the ICD10 code, you contribute as unexposed (trt=0) from when you enter to when you get the ICD10 code
  data1 <- subset(total.comorb, pred == 1)
  data1$tstart <- 0
  data1$tstop <- as.numeric(data1$icd10_date - data1$startfollowup)
  data1$outP <- 0
  data1$trt <- 0
  
  # If you have the ICD10 code, you contribute as exposed (trt=1) from when you get the ICD10 code until end-of-follow-up (2015/MI/death)
  data2 <- subset(total.comorb, pred == 1)
  data2$tstart <- as.numeric(data2$icd10_date - data2$startfollowup)
  data2$tstop <- as.numeric(data2$endfollowup - data2$startfollowup)
  data2$outP <- data2$mi
  data2$trt <- 1
  
  # If you don't have the ICD10 code, but you have MI, you contribute as unexposed (trt=0) from when you enter until MI
  data3 <- subset(total.comorb, pred == 0 & mi == 1)
  data3$tstart <- 0
  data3$tstop <- as.numeric(data3$endfollowup - data3$startfollowup)
  data3$outP <- 1
  data3$trt <- 0
  
  dataF <- rbind(data1, data2, data3)
  
  # If you don't have the ICD10 code or MI, you contribute as unexposed (trt=0) from when you enter until end-of-follow-up
  otherdata <- subset(total.comorb, pred != 1 & mi != 1)
  
  otherdata$tstart <- 0
  otherdata$tstop <- as.numeric(otherdata$endfollowup - otherdata$startfollowup)
  otherdata$outP <- 0
  otherdata$trt <- 0
  
  findata <- rbind(dataF, otherdata)
  
  # remove one individual where tstart = tstop
  findata <- findata[findata$tstart < findata$tstop,]

  # base model
  mod <- tryCatch({
    coxph(Surv(tstart, tstop, outP) ~ trt, data = findata)
  }, error = function(w) {
    NA
  })
  cindex_icd_base <- survConcordance(Surv(tstart, tstop, outP) ~ predict(mod), data = findata)$concordance

  # model with age and sex
  mod1 <- tryCatch({
    coxph(Surv(tstart, tstop, outP) ~ age + trt + sex, data = findata)
  }, error = function(w) {
    NA
  })
  cindex_icd <- survConcordance(Surv(tstart, tstop, outP) ~ predict(mod1), data = findata)$concordance

  # model with risk factors
  mod2 <- tryCatch({
    coxph(Surv(tstart, tstop, outP) ~ age + trt + sex + sbp + smoke + bmi + diabetes + lipid_lowering, data = findata)
  }, error = function(w) {
    NA
  })
  cindex_icd_risk_factor <- survConcordance(Surv(tstart, tstop, outP) ~ predict(mod2), data = findata[!is.na(findata$sbp) & !is.na(findata$smoke) & !is.na(findata$bmi) & !is.na(findata$diabetes) & !is.na(findata$lipid_lowering),])$concordance


  # base model of interaction
  mod3 <- tryCatch({
    coxph(Surv(tstart, tstop, outP) ~ trt*PRS_0.5, data = findata)
  }, error = function(w) {
    NA
  })

  # base model of interaction with age and sex
  mod4 <- tryCatch({
    coxph(Surv(tstart, tstop, outP) ~ age + trt*PRS_0.5 + sex, data = findata)
  }, error = function(w) {
    NA
  })

  # base model of interaction with risk factors
  mod5 <- tryCatch({
    coxph(Surv(tstart, tstop, outP) ~ age + trt*PRS_0.5 + sex + sbp + smoke + bmi + diabetes + lipid_lowering, data = findata)
  }, error = function(w) {
    NA
  })

  
  res <- rbind(res, c(coef(summary(mod))[1, c(1, 4)], coef(summary(mod1))[2, c(1, 4)], coef(summary(mod2))[2, c(1, 4)], coef(summary(mod3))[3, c(1, 4)], coef(summary(mod4))[5, c(1, 4)], coef(summary(mod5))[10, c(1, 4)], mean_day_dist, sd_day_dist, cindex_icd_base, cindex_icd, cindex_icd_risk_factor))

  icd.final <- c(icd.final, icd.all[i])
  
  print(i)

}


#fix this stuff
resdf = as.data.frame(res)
resdf = cbind(icd10 = icd.all, resdf)
colnames(resdf) <- c("ICD10", "beta_base", "z_base", "beta_icd", "z_icd", "beta_risk_factor", "z_risk_factor", "beta_interaction_base", "z_interaction_base", "beta_interaction", "z_interaction", "beta_interaction_risk_factor", "z_interaction_risk_factor", "mean_day_dist", "sd_day_dist", "cindex_icd_base", "cindex_icd", "cindex_icd_risk_factor")
write.table(resdf, file = paste0(path, 'survivalanalysis_pre_sensitivity1.tsv'), sep = "\t", row.names = FALSE, col.names = TRUE, quote=F)

##################################################
# Goal: Find association between MI and PRS
##################################################

# create data set for survival analysis
total.newTEMP <- merge(total.newU, total.newMI, by = "eid", all.x = T)
total.newTEMP$mi <- ifelse(is.na(total.newTEMP$mi_date), 0, 1)
# create models and calculate C-index for different models of PRS
mod_base <- coxph(Surv(endfollowup - startfollowup,mi) ~ age + sex, data = total.newTEMP)
cindex_base <- survConcordance(Surv(endfollowup - startfollowup, mi) ~ predict(mod_base), data = total.newTEMP)$concordance
# C-index = 0.7175515
mod_risk_factor <- coxph(Surv(endfollowup - startfollowup,mi) ~ age + sex + sbp + smoke + bmi + diabetes + lipid_lowering, data = total.newTEMP)
cindex_risk_factor <- survConcordance(Surv(endfollowup - startfollowup, mi) ~ predict(mod_risk_factor), data = total.newTEMP[!is.na(total.newTEMP$sbp) & !is.na(total.newTEMP$smoke) & !is.na(total.newTEMP$bmi) & !is.na(total.newTEMP$diabetes) & !is.na(total.newTEMP$lipid_lowering),])$concordance
# C-index = 0.7440676
mod_prs_base <- coxph(Surv(endfollowup - startfollowup, mi) ~ scale(PRS_0.5), data = total.newTEMP)
cindex_prs_base <- survConcordance(Surv(endfollowup - startfollowup, mi) ~ predict(mod_prs_base), data = total.newTEMP[!is.na(total.newTEMP$PRS_0.5),])$concordance
# C-index = 0.5614979
mod_prs <- coxph(Surv(endfollowup - startfollowup,mi) ~ age + scale(PRS_0.5) + sex, data = total.newTEMP)
cindex_prs <- survConcordance(Surv(endfollowup - startfollowup, mi) ~ predict(mod_prs), data = total.newTEMP[!is.na(total.newTEMP$PRS_0.5),])$concordance
# C-index = 0.7252172
mod_prs_risk_factor <- coxph(Surv(endfollowup - startfollowup, mi) ~ age + scale(PRS_0.5) + sex + sbp + smoke + bmi + diabetes + lipid_lowering, data = total.newTEMP)
cindex_prs_risk_factor <- survConcordance(Surv(endfollowup - startfollowup, mi) ~ predict(mod_prs_risk_factor), data = total.newTEMP[!is.na(total.newTEMP$sbp) & !is.na(total.newTEMP$smoke) & !is.na(total.newTEMP$bmi) & !is.na(total.newTEMP$diabetes) & !is.na(total.newTEMP$lipid_lowering) & !is.na(total.newTEMP$PRS_0.5),])$concordance
# C-index = 0.7487259

# subset population with no previous comorbidities 
no.icd.ids <- unique(total.new$eid[is.na(total.new$icd10)])
icd.ids <- unique(total.new$eid[!is.na(total.new$icd10)])
no.icd.ids <- unique(no.icd.ids[!(no.icd.ids %in% icd.ids)])
total.newNA <- total.new[total.new$eid %in% no.icd.ids,]
total.newNA$mi <- ifelse(is.na(total.newNA$mi_date), 0, 1)
total.newNA <- total.newNA[!duplicated(total.newNA$eid), c("eid", "sex", "age", "death", "death_date", "mi", "mi_date", "icd10", "icd10_date", "PRS_0.5", "startfollowup", "endfollowup")]
total.newNA <- total.newNA[which(total.newNA$endfollowup > total.newNA$startfollowup),]
# check association between PRS and MI for individuals with no previous comorbidities
mod_prs_no_comorb <- coxph(Surv(endfollowup - startfollowup,mi) ~ scale(PRS_0.5), data = total.newNA)
cindex_prs_no_comorb <- survConcordance(Surv(endfollowup - startfollowup, mi) ~ predict(mod_prs_no_comorb), data = total.newNA[!is.na(total.newNA$PRS_0.5),])$concordance
err_prs_no_comorb <- as.numeric(survConcordance(Surv(endfollowup - startfollowup, mi) ~ predict(mod_prs_no_comorb), data = total.newNA[!is.na(total.newNA$PRS_0.5),])$std.err)
no_comorb <- cbind(group = "None", cindex = cindex_prs_no_comorb, error = err_prs_no_comorb)

# plot graphs based on statistics calculated across all ICD10 codes
resdf <- read.table(file = paste0(path, 'survivalanalysis_pre_sensitivity1.tsv'), sep = "\t", header = TRUE, stringsAsFactors = F)
resdf$logp_base <- -log10(2*pnorm(-abs(resdf$z_base)))
resdf$logp_base[resdf$logp_base == "Inf"] <- 300
resdf$logp_icd <- -log10(2*pnorm(-abs(resdf$z_icd)))
resdf$logp_icd[resdf$logp_icd == "Inf"] <- 300
resdf$logp_risk_factor <- -log10(2*pnorm(-abs(resdf$z_risk_factor)))
resdf$logp_risk_factor[resdf$logp_risk_factor == "Inf"] <- 300
resdf$hr_base <- exp(resdf$beta_base)
resdf$ICD10label <- resdf$ICD10
resdf$ICD10label[resdf$logp < -log10(0.05 / length(icd.all))] <- ""
resdf$ICD10 <- factor(resdf$ICD10, levels = resdf$ICD10[order(resdf$mean_day_dist)])
# What is the distribution of time to MI for each ICD10 code?
ggplot(resdf, aes(ICD10, mean_day_dist)) + 
  geom_bar(stat = "identity") + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1,size=4), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x = "ICD10 Code", y = "Mean Time Distance to MI (in days)")
# What is the association between each ICD10 code and MI? (p-value)
ggplot(data = resdf, mapping = aes(ICD10, logp_base)) + 
  geom_point() +
  geom_text_repel(mapping = aes(label = resdf$ICD10label), size = 3, segment.alpha = 0.5) + 
  geom_hline(yintercept = -log10(0.05 / length(icd.all)), size = 1.5, color = "red") + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = "ICD10 Code", y = "-log10 Transformed P-Value") 
# What is the association between each ICD10 code and MI? (p-value and hazard ratio)
ggplot(data = resdf, mapping = aes(ICD10, logp_base)) + 
  geom_point(mapping = aes(size = exp(beta_base))) +
  scale_size_continuous("Hazard Ratio", breaks = c(2, 5, 10, 15, 20, 25, 30), labels=c("2", "5", "10", "15", "20", "25", "30"))
  geom_text_repel(data = resdf, mapping = aes(label = resdf$ICD10label), size = 3, segment.alpha = 0.5) + 
  geom_hline(yintercept = -log10(0.05 / length(icd.all)), size = 1.5, color = "red") + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = "ICD10 Code", y = "-log10 Transformed P-Value") 
# What is the association between each ICD10 code and MI, adjusted for other risk factors? (adjusted and unadjusted p-value)
resdf2 <- melt(resdf, measure.vars = c("logp_risk_factor","logp_base"))
resdf2$ICD10label[resdf2$value < -log10(0.05 / length(icd.all)) | resdf2$variable == "logp_main"] <- ""
ggplot(data = resdf2, mapping = aes(ICD10, value)) + 
  geom_point(mapping = aes(size = exp(beta_main), col = variable)) + 
  geom_text_repel(mapping = aes(label = resdf2$ICD10label),size=3,segment.alpha=0.5) +
  geom_hline(yintercept = -log10(0.05 / length(icd.all)), size = 1.5, color = "red") + 
  scale_size_continuous("HR", breaks = c(2, 5, 10, 15, 20, 25, 30), labels=c("2","5","10","15","20","25","30")) +
  scale_colour_discrete("", breaks = c("logp_risk_factor", "logp_main"), labels=c("Adjusted", "Unadjusted")) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = "ICD Code", y = "-log10 Transformed P-Value") 
# What is the comparison of the hazard ratio for ICD10 codes and PRS?
ggplot(data = resdf, mapping = aes(ICD10, hr_base)) + 
  geom_point() + 
  geom_text_repel(mapping = aes(label = resdf$ICD10label), size = 3, segment.alpha = 0.5) +
  geom_hline(yintercept = coef(summary(mod_prs_dich))[2,2], size = 1.5, color = "red") + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = "ICD Code", y = "Hazard Ratio")
# What is the trend between predictive ability and hazard ratio with MI for each ICD10 code?
ggplot(resdf, aes(hr_base, cindex_base)) + 
  geom_point() + 
  geom_text_repel(aes(label = resdf$ICD10label), size = 3, segment.alpha = 0.5) +
  geom_vline(xintercept = coef(summary(mod_prs_base))[1, 2], size = 1.5, color = "red") + 
  geom_hline(yintercept = cindex_prs_base, size = 1.5, color = "red") + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = "Hazard Ratio", y = "C-index")
# Comparing C-index of ICD10 codes with PRS for different models
ggplot() + 
  geom_point(data = resdf, mapping = aes(ICD10, cindex_icd_base, color = "base line")) +
  geom_point(data = resdf, mapping = aes(ICD10, cindex_icd, color = "age and sex")) +
  geom_point(data = resdf, mapping = aes(ICD10, cindex_icd_risk_factor, color = "risk factor")) +
  geom_hline(yintercept = as.numeric(cindex_prs_base)) +
  geom_hline(yintercept = as.numeric(cindex_prs)) +
  geom_hline(yintercept = as.numeric(cindex_prs_risk_factor)) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = "ICD Code", y = "C-index") 

# calculate and plot C-index for population subsets (e.g., individuals with only I20) using model of PRS
sub.pop.table <- NULL
sub.pop.table <- data.frame(sub.pop.table)
for (i in 1:length(icd.all)) {
  ids <- unique(total.new$eid[total.new$icd10 == icd.all[i]])
  total.newI <- total.new[total.new$eid %in% ids,]
  total.newI$mi <- ifelse(is.na(total.newI$mi_date), 0, 1)
  total.newI <- total.newI[which(total.newI$icd10 == icd.all[i]), c("eid", "sex", "age", "death", "death_date", "mi", "mi_date", "icd10", "icd10_date", "PRS_0.5", "startfollowup", "endfollowup")]
  total.newI <- total.newI[!duplicated(total.newI$eid),]
  total.newI <- total.newI[which(total.newI$endfollowup > total.newI$startfollowup),]
  mod_prs_i <- coxph(Surv(endfollowup - startfollowup,mi) ~ scale(PRS_0.5), data = total.newI)
  cindex_prs_i <- survConcordance(Surv(endfollowup - startfollowup, mi) ~ predict(mod_prs_i), data = total.newI[!is.na(total.newI$PRS_0.5),])$concordance
  std.err <- as.numeric(survConcordance(Surv(endfollowup - startfollowup, mi) ~ predict(mod_prs_i), data = total.newI[!is.na(total.newI$PRS_0.5),])$std.err)
  sub.pop.table <- rbind(sub.pop.table, cbind(cindex_prs_i, std.err))
  print(i)
}
rownames(sub.pop.table) <- icd.all
colnames(sub.pop.table) <- c("cindex", "error")
sub.pop.table <- cbind(group = rownames(sub.pop.table), sub.pop.table)
sub.pop.table <- rbind(sub.pop.table, no_comorb)
# plot C-index
ggplot(sub.pop.table, aes(group, as.numeric(cindex))) + 
  geom_boxplot() + 
  geom_errorbar(aes(ymin = as.numeric(cindex)-as.numeric(error), ymax = as.numeric(cindex)+as.numeric(error)), width = 0.2, position = position_dodge(1)) + 
  xlab("Group") + ylab("C-index") + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

# calculate and plot C-index for entire population using model of PRS
total.newSS <- total.new[,c("eid","icd10","mi_date","PRS_0.5", "startfollowup", "endfollowup")]
total.newSS$mi <- ifelse(is.na(total.new$mi_date), 0, 1)
# loop across all ICD10 codes and create dichotomous variable (e.g., 1 = individual has I20 and 0 = individual does not have I20)
for (i in 1:length(icd.all)) {
  temp_ids <- total.newSS$eid[total.newSS$icd10 == icd.all[i]]
  total.newSS[,ncol(total.newSS)+1] <- ifelse(total.newSS$eid %in% temp_ids, 1, 0)
  colnames(total.newSS)[ncol(total.newSS)] <- icd.all[i]
  print(i)
}
total.newSS <- total.newSS[!duplicated(total.newSS$eid),]
# loop to create dichotomous variable for PRS at top 1%, 2%, 5%, and 10% of PRS
cutoffs <- c("top1", "top2", "top5", "top10")
total.newSS$top1 <- ifelse(total.newSS$PRS_0.5 > quantile(total.newSS$PRS_0.5, 0.99, na.rm = TRUE), 1, 0)
total.newSS$top2 <- ifelse(total.newSS$PRS_0.5 > quantile(total.newSS$PRS_0.5, 0.98, na.rm = TRUE), 1, 0)
total.newSS$top5 <- ifelse(total.newSS$PRS_0.5 > quantile(total.newSS$PRS_0.5, 0.95, na.rm = TRUE), 1, 0)
total.newSS$top10 <- ifelse(total.newSS$PRS_0.5 > quantile(total.newSS$PRS_0.5, 0.90, na.rm = TRUE), 1, 0)
# loop across all ICD10 codes and find sensitivity and specificity
pop.table <- NULL
for (i in 1:length(icd.all)) {
  sens <- confusionMatrix(data = factor(total.newSS[,icd.all[i]]), reference = factor(total.newSS$mi), positive = "1")$byClass["Sensitivity"]
  spec <- confusionMatrix(data = factor(total.newSS[,icd.all[i]]), reference = factor(total.newSS$mi), positive = "1")$byClass["Specificity"]
  ss <- cbind(sens, spec)
  pop.table <- rbind(pop.table, ss)
  print(i)
}
# loop across cutoffs for dichotomous PRS and find sensitivity and specificity
for (i in 1:length(cutoffs)) {
  sens <- confusionMatrix(data = factor(total.newSS[,cutoffs[i]]), reference = factor(total.newSS$mi), positive = "1")$byClass["Sensitivity"]
  spec <- confusionMatrix(data = factor(total.newSS[,cutoffs[i]]), reference = factor(total.newSS$mi), positive = "1")$byClass["Specificity"]
  ss <- cbind(sens, spec)
  pop.table <- rbind(pop.table, ss)
  print(i)
}
rownames(pop.table) <- c(icd.all, cutoffs)
colnames(pop.table) <- c("sens", "spec")
pop.table <- data.frame(pop.table)
pop.table$fpr <- 1-as.numeric(pop.table$spec)
pop.table$label <- rownames(pop.table)
pop.table$type <- ifelse(pop.table$label %in% icd.all, "icd10", "prs")

# find AUC of dichotomous PRS
auc.prs <- as.numeric(roc(total.newSS$mi, total.newSS$PRS_0.5)$auc)
# C-index = 0.562317
auc.prs.top1 <- as.numeric(roc(total.newSS$mi, total.newSS$top1)$auc)
# C-index = 0.5041281
auc.prs.top2 <- as.numeric(roc(total.newSS$mi, total.newSS$top2)$auc)
# C-index = 0.5078302
auc.prs.top5 <- as.numeric(roc(total.newSS$mi, total.newSS$top5)$auc)
# C-index = 0.5151975
auc.prs.top10 <- as.numeric(roc(total.newSS$mi, total.newSS$top10)$auc)
# C-index = 0.5228124
cindex_dich <- data.frame(cutoffs, cindex = c(auc.prs.top1, auc.prs.top2, auc.prs.top5, auc.prs.top10))
# add C-index
pop.table$cindex <- ifelse(pop.table$label %in% icd.all, resdf$cindex_icd_base, "")
pop.table[c("top1", "top2", "top5", "top10"), "cindex"] <- cindex_dich$cindex

write.table(pop.table, file = paste0(path, 'pred_power_icd10_base.tsv'), sep = "\t", row.names = FALSE, col.names = TRUE, quote=F)
pop.table <- read.table(file = paste0(path, 'pred_power_icd10_base.tsv'), sep = "\t", header = TRUE, stringsAsFactors = F)

# make labels
pop.table$glabel <- ""
for (i in 1:nrow(pop.table)) {
  icd10.label <- pop.table[i, "label"]
  cindex.label <- pop.table[i, "cindex"]
  glabel <- paste(icd10.label, "", substr(cindex.label, 0, 5))
  pop.table[i, "glabel"] <- glabel
  print(i)
}
pop.table$glabel <- ifelse(pop.table$cindex > 0.5041281, pop.table[, "glabel"], "")

# get sensitivities and specificities for PRS from ROC
roc.table <- cbind(roc(total.newSS$mi, total.newSS$PRS_0.5)$sensitivities, roc(total.newSS$mi, total.newSS$PRS_0.5)$specificities)
roc.table <- data.frame(roc.table)
colnames(roc.table) <- c("sens", "spec")
roc.table$fpr <- 1-roc.table$spec

# compare ROC for group with previous comorbidities and group with no previous comorbidities
prev.comorb.roc <- cbind(roc(total.newSS$mi, total.newSS$PRS_0.5)$sensitivities, roc(total.newSS$mi, total.newSS$PRS_0.5)$specificities)
prev.comorb.roc <- data.frame(prev.comorb.roc)
colnames(prev.comorb.roc) <- c("sens", "spec")
prev.comorb.roc$fpr <- 1-prev.comorb.roc$spec
no.prev.comorb.roc <- cbind(roc(total.newNA$mi, total.newNA$PRS_0.5)$sensitivities, roc(total.newNA$mi, total.newNA$PRS_0.5)$specificities)
no.prev.comorb.roc <- data.frame(no.prev.comorb.roc)
colnames(no.prev.comorb.roc) <- c("sens", "spec")
no.prev.comorb.roc$fpr <- 1-no.prev.comorb.roc$spec
ggplot() +
  geom_point(data = prev.comorb.roc, mapping = aes(fpr, sens, color = "ICD10"), size = 0.5) +
  geom_point(data = no.prev.comorb.roc, mapping = aes(fpr, sens, color = "No ICD10"), size = 0.5) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + 
  geom_abline(slope = 1, intercept = 0) + 
  annotate(geom = "text", x = 0.3, y = 0.75, label = paste("C-index (ICD10): ", substr(as.character(as.numeric(roc(total.newSS$mi, total.newSS$PRS_0.5)$auc)), 0, 5))) +
  annotate(geom = "text", x = 0.3, y = 0.85, label = paste("C-index (No ICD10): ", substr(as.character(as.numeric(roc(total.newNA$mi, total.newNA$PRS_0.5)$auc)), 0, 5)))
  
# plot several ROC with annotations
# plot ROC 
ggplot() + 
  geom_point(data = roc.table, mapping = aes(fpr, sens), size = 1) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + 
  geom_abline(slope = 1, intercept = 0) + 
  annotate(geom = "text", x = 0.5, y = 0.75, label = paste("C-index: ", substr(as.character(as.numeric(roc(total.newSS$mi, total.newSS$PRS_0.5)$auc)), 0, 5)))
# plot ROC with ICD10 codes and dichotomous PRS
ggplot() + 
  geom_point(data = roc.table, mapping = aes(fpr, sens), size = 1) +
  geom_point(data = pop.table, mapping = aes(fpr, sens, color = type)) +
  geom_label_repel(data = pop.table, aes(fpr, sens, label = glabel), size = 3, segment.alpha = 0.5) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + 
  geom_abline(slope = 1, intercept = 0) + 
  annotate(geom = "text", x = 0.5, y = 0.75, label = paste("C-index: ", substr(as.character(as.numeric(roc(total.newSS$mi, total.newSS$PRS_0.5)$auc)), 0, 5)))
# plot zoomed in ROC with ICD10 codes and dichotomous PRS
ggplot() + 
  geom_point(data = roc.table, mapping = aes(fpr, sens), size = 1) +
  geom_point(data = pop.table, mapping = aes(fpr, sens, color = type)) +
  geom_label_repel(data = pop.table, aes(fpr, sens, label = glabel), size = 3, segment.alpha = 0.5) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  coord_cartesian(xlim = c(0, 0.1), ylim = c(0, 0.15)) + 
  geom_abline(slope = 1, intercept = 0) + 
  annotate(geom = "text", x = 0.025, y = 0.1, label = paste("C-index: ", substr(as.character(as.numeric(roc(total.newSS$mi, total.newSS$PRS_0.5)$auc)), 0, 5)))