##################################################
# MYOCARDIAL INFARCTION SENSITIVITY ANALYSIS
##################################################

##################################################
# Install and load packages
##################################################

install.packages("survival")
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("reshape2")
install.packages("caret")
install.packages("pROC")
install.packages("plotROC")
install.packages("data.table")
library(survival)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(caret)
library(pROC)
library(plotROC)
library(data.table)

##################################################
# Load and format data
##################################################

total <- read.table(file = 'data.tsv')
# Columns in data should follow following format:
# eid            - character for ID
# age            - int/num for age
# sex            - boolean for sex
# birth_date     - date for day of birth in formt YYYY-MM-DD (e.g., 1999-10-16)
# death          - boolean for death
# death_date     - date for day of death in formt YYYY-MM-DD (e.g., 1999-10-16)
# prs            - int/num for polygenic risk score
# sbp            - int/num for systolic blood pressure
# smoke          - boolean for smoking status
# bmi            - int/num for body mass index
# diabetes       - boolean for diabetes status
# lipid_lowering - boolean for lipid-lowering medications
# tot_chol       - int/num for total cholesterol
# hdl_chol       - int/num for HDL cholesterol
# mi             - boolean for myocardial infarction
# mi_date        - date for day of FIRST myocardial infarction in formt YYYY-MM-DD (e.g., 1999-10-16)
# ac_date        - date for baseline in formt YYYY-MM-DD (e.g., 1999-10-16)
# icd10          - chr for ICD10 event in three character format (e.g., I20.1 and I20.2 will be I20)
# icd10_date     - date for day of ICD10 event in formt YYYY-MM-DD (e.g., 1999-10-16)
# Rows in data should be an ICD10 code for an individual (one individual can have multiple rows)
total$eid <- as.character(total$eid)
total$age <- as.numeric(total$age)
total$sex <- as.integer(total$sex)
total$birth_date <- as.Date(total$birth_date, "%Y-%m-%d")
total$death <- as.integer(total$death)
total$death_date <- as.Date(total$death_date, "%Y-%m-%d")
total$prs <- as.numeric(total$prs)
total$sbp <- as.numeric(total$sbp)
total$smoke <- as.integer(total$smoke)
total$bmi <- as.numeric(total$bmi)
total$diabetes <- as.integer(total$diabetes)
total$lipid_lowering <- as.integer(total$lipid_lowering)
total$tot_chol <- as.integer(total$tot_chol)
total$hdl_chol <- as.integer(total$hdl_chol)
total$mi <- as.integer(total$mi)
total$mi_date <- as.Date(total$mi_date, "%Y-%m-%d")
total$ac_date <- as.Date(total$ac_date, "%Y-%m-%d")
total$icd10 <- as.character(total$icd10)
total$icd10_date <- as.Date(total$icd10_date, "%Y-%m-%d")
nrow(total)
#_____ rows
length(unique(total$eid))
#_____ individuals
length(unique(total$icd10))
#_____ ICD codes
length(unique(total[which(!is.na(total$mi_date)),]$eid))
#_____ individuals with MI

##################################################
# Clean data based on exclusion criteria
##################################################

# create start-of-follow-up (entrance in the study) and end-of-follow-up (exit in the study, date of MI, or date of death death)
# CHANGE: date of beginning of study and date of end of study will be different
study.start <- as.Date("YYYY-MM-DD", "%Y-%m-%d")
study.end <- as.Date("YYYY-MM-DD", "%Y-%m-%d")
total$startfollowup <- total$ac_date
total$endfollowup <- ifelse(!is.na(total$death_date) & is.na(total$mi_date), as.Date(total$death_date, "%Y-%m-%d"), ifelse(!is.na(total$mi_date), as.Date(total$mi_date, "%Y-%m-%d"), study.end))
class(total$startfollowup) <- "Date"
class(total$endfollowup) <- "Date"

# exclude ICD10 codes that are not related to disease (e.g., car crash) 
total.new <- total[,c("eid", "age", "sex", "birth_date", "death", "death_date", "prs", "sbp", "smoke", "bmi", "diabetes", "lipid_lowering", "tot_chol", "hdl_chol", "mi", "mi_date", "ac_date", "icd10", "icd10_date", "startfollowup", "endfollowup")]
to_remove <- unique(total.new$icd10[grepl("^V|^W|^X|^Y|^Z", total.new$icd10)])
to_remove <- unique(to_remove[!is.na(to_remove)])
total.new$icd10[grepl("^V|^W|^X|^Y|^Z", total.new$icd10)] <- NA
total.new$icd10_date[grepl("^V|^W|^X|^Y|^Z", total.new$icd10)] <- NA
if (length(which(total.new$icd10 %in% to_remove == "TRUE")) != 0) {
  print("ERROR: Exclude ICD10 codes that are not related to disease (e.g., car crash).")
} else {
  print("NO ERROR.")
}
nrow(total.new)
#_____ rows
length(unique(total.new$eid))
#_____ individuals
length(unique(total.new$icd10))
#_____ ICD codes
length(unique(total.new[which(!is.na(total.new$mi_date)),]$eid))
#_____ individuals with MI

# exclude ICD10 codes that do not have dates
total.new$icd10 <- ifelse(!is.na(total.new$icd10) & is.na(total.new$icd10_date), NA, total.new$icd10)
if (length(which(!is.na(total.new$icd10) & is.na(total.new$icd10_date)) != 0)) {
  print("ERROR: Exclude ICD10 codes that do not have dates.") 
} else {
  print("NO ERROR.")
}
nrow(total.new)
#_____ rows
length(unique(total.new$eid))
#_____ individuals
length(unique(total.new$icd10))
#_____ ICD codes
length(unique(total.new[which(!is.na(total.new$mi_date)),]$eid))
#_____ individuals with MI

# remove individuals with MI outside of timeframe 
to_remove <- total.new$eid[total.new$mi_date > study.end | total.new$mi_date < total.new$startfollowup]
to_remove <- unique(to_remove[!is.na(to_remove)])
total.new = total.new[!total.new$eid %in% to_remove,]
if (length(which(total.new$eid %in% to_remove == "TRUE")) != 0) {
  print("ERROR: Remove individuals with MI outside of timeframe.")
} else {
  print("NO ERROR.")
}
nrow(total.new)
#_____ rows
length(unique(total.new$eid))
#_____ individuals
length(unique(total.new$icd10))
#_____ ICD codes
length(unique(total.new[which(!is.na(total.new$mi_date)),]$eid))
#_____ individuals with MI

# SENSITIVITY ANALYSIS: remove ICD10 codes related to MI (I20 through I25)
total.new2 <- total.new
to_remove <- total.new2$eid[total.new2$icd10 %in% c("I20", "I21", "I22", "I23", "I24", "I25")]
to_remove <- unique(to_remove[!is.na(to_remove)])
total.new2 = total.new2[!total.new2$eid %in% to_remove,]
if (length(which(total.new2$eid %in% to_remove == "TRUE")) != 0) {
  print("ERROR: Remove ICD10 codes related to MI.")
} else {
  print("NO ERROR.")
}
nrow(total.new2)
#_____ rows
length(unique(total.new2$eid))
#_____ individuals 
length(unique(total.new2$icd10))
#_____ ICD codes
length(unique(total.new2[which(!is.na(total.new2$mi_date)),]$eid))
#_____ individuals with MI 

# SENSITIVITY ANALYSIS: keep ICD10 codes only before baseline
total.new2$icd10[total.new2$icd10_date > total.new2$startfollowup | total.new2$icd10_date < study.start] <- NA
total.new2$icd10_date[total.new2$icd10_date > total.new2$startfollowup | total.new2$icd10_date < study.start] <- NA
total.new2$icd10[(total.new2$mi_date - total.new2$icd10_date) < 8] <- NA
total.new2$icd10_date[(total.new2$mi_date - total.new2$icd10_date) < 8] <- NA
# remove duplicate rows or second diagnoses (e.g., for I20 diagnosis on 1999-10-16 and 2005-11-14, keep I20 diagnosis on 1999-10-16 and remove I20 diagnosis on 2005-11-14)
total.new2 <- total.new2[order(total.new2$eid, total.new2$icd10, total.new2$icd10_date),]
total.new2 <- total.new2[!duplicated(total.new2[, c("eid", "icd10")], fromLast = FALSE),]
nrow(total.new2)
#_____ rows
length(unique(total.new2$eid))
#_____ individuals 
length(unique(total.new2$icd10))
#_____ ICD codes
length(unique(total.new2[which(!is.na(total.new2$mi_date)),]$eid))
#_____ individuals with MI 















##################################################
# Survival analysis for PRS
total.new <- total.new2 # for sensitivity analysis
##################################################

# create datasets for survival analysis
total.newMOD <- total.new[!duplicated(total.new$eid), c("eid", "age", "sex", "birth_date", "death", "death_date", "prs", "sbp", "smoke", "bmi", "diabetes", "lipid_lowering", "tot_chol", "hdl_chol", "mi", "mi_date", "ac_date", "startfollowup", "endfollowup")]
total.newMOD <- total.newMOD[!is.na(total.newMOD$prs) & !is.na(total.newMOD$sbp) & !is.na(total.newMOD$smoke) & !is.na(total.newMOD$bmi) & !is.na(total.newMOD$diabetes) & !is.na(total.newMOD$lipid_lowering) & !is.na(total.newMOD$tot_chol) & !is.na(total.newMOD$hdl_chol),]
# survival analysis models and c-index for PRS
mod_risk_factor <- coxph(Surv(endfollowup - startfollowup,mi) ~ age + sex + sbp + smoke + bmi + diabetes + lipid_lowering + tot_chol + hdl_chol, data = total.newMOD)
cindex_risk_factor <- survConcordance(Surv(endfollowup - startfollowup, mi) ~ predict(mod_risk_factor), data = total.newMOD)$concordance
# c-index = 
mod_prs <- coxph(Surv(endfollowup - startfollowup,mi) ~ scale(prs), data = total.newMOD)
cindex_prs <- survConcordance(Surv(endfollowup - startfollowup, mi) ~ predict(mod_prs), data = total.newMOD)$concordance
# c-index = 
mod_prs_base <- coxph(Surv(endfollowup - startfollowup,mi) ~ age + scale(prs) + sex, data = total.newMOD)
cindex_prs_base <- survConcordance(Surv(endfollowup - startfollowup, mi) ~ predict(mod_prs_base), data = total.newMOD)$concordance
# c-index = 
mod_prs_risk_factor <- coxph(Surv(endfollowup - startfollowup, mi) ~ age + scale(prs) + sex + sbp + smoke + bmi + diabetes + lipid_lowering + tot_chol + hdl_chol, data = total.newMOD)
cindex_prs_risk_factor <- survConcordance(Surv(endfollowup - startfollowup, mi) ~ predict(mod_prs_risk_factor), data = total.newMOD)$concordance
# c-index = 
numb <- rbind(cindex_risk_factor, coef(summary(mod_prs))["scale(prs)", "exp(coef)"], cindex_prs, 
  coef(summary(mod_prs_base))["scale(prs)", "exp(coef)"], cindex_prs_base,
  coef(summary(mod_prs_risk_factor))["scale(prs)", "exp(coef)"], cindex_prs_risk_factor)
numb <- data.frame(numb)
rownames(numb) <- c("cindex_risk_factor", "hr_prs", "cindex_prs",
  "hr_prs_base", "cindex_prs_base", 
  "hr_prs_risk_factor", "cindex_prs_risk_factor")
colnames(numb) <- c("numb")
write.table(numb, file = 'sensitivity_prs.tsv', sep = "\t", row.names = TRUE, col.names = TRUE)

##################################################
# Survival analysis for ICD10 codes
##################################################

# identify list of ICD10 codes for analysis
icd.mi <- unique(total.new[!is.na(total.new$mi_date), c("eid", "icd10")])
icd.mi <- icd.mi[!is.na(icd.mi$icd10), "icd10"]
icd.all <- names(table(icd.mi))[table(icd.mi) > 10]
res <- NULL

# create datasets for survival analysis
total.newT <- total.new[!is.na(total.new$mi_date),]
total.newMI <- total.newT[!duplicated(total.newT$eid), c("eid", "mi", "mi_date")]
total.newU <- total.new[!duplicated(total.new$eid), c("eid", "age", "sex", "birth_date", "death", "death_date", "prs", "sbp", "smoke", "bmi", "diabetes", "lipid_lowering", "tot_chol", "hdl_chol", "ac_date", "startfollowup", "endfollowup")]

for (i in 1:length(icd.all)) {
  
  # keep unique record for comorbity examined
  total.newT <- total.new[total.new$icd10 == icd.all[i] & !is.na(total.new$icd10),]
  # sort and keep first event for each individual
  total.new.sorted <- total.newT[order(total.newT$eid, as.numeric(total.newT$icd10_date)),]
  total.new.sortedU <- total.new.sorted[!duplicated(total.new.sorted$eid), c("eid", "icd10", "icd10_date")]
  # merge ICD10 dataset, main dataset, and MI-only dataset
  temp1 <- merge(total.newU, total.new.sortedU, all.x = T, by = "eid")
  total.comorb <- merge(temp1, total.newMI, all.x = T, by = "eid")
  # assign boolen for ICD10 codes (0/1) and MI (0/1)
  total.comorb$pred <- ifelse(is.na(total.comorb$icd10_date), 0, 1)
  total.comorb$mi <- ifelse(is.na(total.comorb$mi_date), 0, 1)
  
  # calculate time distance between ICD10 codes and MI
  total.comorbT <- total.comorb[total.comorb$pred == 1 & total.comorb$mi == 1,]
  mean_day_dist <- mean(total.comorbT$mi_date - total.comorbT$icd10_date)
  sd_day_dist <- sd(total.comorbT$mi_date - total.comorbT$icd10_date)
  
  # If you have the ICD10 code, you contribute as unexposed (trt=0) from when you enter to when you get the ICD10 code
  data1 <- subset(total.comorb, pred == 1)
  data1$tstart <- 0
  data1$tstop <- as.numeric(data1$icd10_date - data1$startfollowup)
  data1$outP <- 0
  data1$trt <- 0
  # If you have the ICD10 code, you contribute as exposed (trt=1) from when you get the ICD10 code until end-of-follow-up
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
  # create final dataset
  findata <- rbind(dataF, otherdata)
  findata <- findata[findata$tstart < findata$tstop,]
  findata <- findata[!is.na(findata$prs) & !is.na(findata$sbp) & !is.na(findata$smoke) & !is.na(findata$bmi) & !is.na(findata$diabetes) & !is.na(findata$lipid_lowering) & !is.na(findata$tot_chol) & !is.na(findata$hdl_chol),]

  # base model
  mod <- tryCatch({ coxph(Surv(tstart, tstop, outP) ~ trt, data = findata) }, error = function(w) { NA })
  cindex_icd_base <- survConcordance(Surv(tstart, tstop, outP) ~ predict(mod), data = findata)$concordance
  # model with age and sex
  mod1 <- tryCatch({ coxph(Surv(tstart, tstop, outP) ~ age + trt + sex, data = findata) }, error = function(w) { NA })
  cindex_icd <- survConcordance(Surv(tstart, tstop, outP) ~ predict(mod1), data = findata)$concordance
  # model with risk factors
  mod2 <- tryCatch({ coxph(Surv(tstart, tstop, outP) ~ age + trt + sex + sbp + smoke + bmi + diabetes + lipid_lowering + tot_chol + hdl_chol, data = findata) }, error = function(w) { NA })
  cindex_icd_risk_factor <- survConcordance(Surv(tstart, tstop, outP) ~ predict(mod2), data = findata)$concordance
  # base model of interaction
  mod3 <- tryCatch({ coxph(Surv(tstart, tstop, outP) ~ trt*prs, data = findata) }, error = function(w) { NA })
  cindex_interaction_base <- survConcordance(Surv(tstart, tstop, outP) ~ predict(mod3), data = findata)$concordance
  # base model of interaction with age and sex
  mod4 <- tryCatch({ coxph(Surv(tstart, tstop, outP) ~ age + trt*prs + sex, data = findata) }, error = function(w) { NA })
  cindex_interaction <- survConcordance(Surv(tstart, tstop, outP) ~ predict(mod4), data = findata)$concordance
  # base model of interaction with risk factors
  mod5 <- tryCatch({ coxph(Surv(tstart, tstop, outP) ~ age + trt*prs + sex + sbp + smoke + bmi + diabetes + lipid_lowering + tot_chol + hdl_chol, data = findata) }, error = function(w) { NA })
  cindex_interaction_risk_factor <- survConcordance(Surv(tstart, tstop, outP) ~ predict(mod5), data = findata)$concordance

  res <- rbind(res, c(coef(summary(mod))["trt", c("coef", "z")], coef(summary(mod1))["trt", c("coef", "z")], coef(summary(mod2))["trt", c("coef", "z")], coef(summary(mod3))["trt:prs", c("coef", "z")], coef(summary(mod4))["trt:prs", c("coef", "z")], coef(summary(mod5))["trt:prs", c("coef", "z")], 
    mean_day_dist, sd_day_dist, cindex_icd_base, cindex_icd, cindex_icd_risk_factor, cindex_interaction_base, cindex_interaction, cindex_interaction_risk_factor,
    length(unique(findata$eid[which(findata$pred == 1 & findata$mi == 1)])),
    length(unique(findata$eid[which(findata$pred == 0 & findata$mi == 1)])),
    length(unique(findata$eid[which(findata$pred == 1 & findata$mi == 0)])),
    length(unique(findata$eid[which(findata$pred == 0 & findata$mi == 0)]))
    ))

  print(i)

}

resdf = as.data.frame(res)
resdf = cbind(icd10 = icd.all, resdf)
rownames(resdf) <- 1:nrow(resdf)
colnames(resdf) <- c("ICD10", "beta_base", "z_base", "beta_icd", "z_icd", "beta_risk_factor", "z_risk_factor", "beta_interaction_base", "z_interaction_base", "beta_interaction", "z_interaction", "beta_interaction_risk_factor", "z_interaction_risk_factor", 
  "mean_day_dist", "sd_day_dist", "cindex_icd_base", "cindex_icd", "cindex_icd_risk_factor", "cindex_interaction_base", "cindex_interaction", "cindex_interaction_risk_factor",
  "icd_mi", "noicd_mi", "icd_nomi", "noicd_nomi")
write.table(resdf, file = 'sensitivity_analysis.tsv', sep = "\t", row.names = TRUE, col.names = TRUE)
resdf.new <- resdf[resdf$ICD10 %in% resdf$ICD10[resdf$logp_icd > -log10(0.05 / length(icd.all))],]
icd.new <- resdf.new$ICD10















##################################################
# Analysis of interaction for population subsets of ICD10 codes
##################################################

# calculate and plot c-index for population subsets (e.g., individuals with only I20) using model of PRS
sub.pop.table <- NULL
for (i in 1:length(icd.all)) {

  # create list of individuals with specific ICD10 code
  ids <- unique(total.new$eid[total.new$icd10 == icd.all[i]])

  # create dataset of individuals with specific ICD10 code
  total.newI <- total.new[total.new$eid %in% ids,]
  total.newI$mi <- ifelse(is.na(total.newI$mi_date), 0, 1)
  total.newI <- total.newI[which(total.newI$icd10 == icd.all[i]), c("eid", "sex", "age", "death", "death_date", "mi", "mi_date", "icd10", "icd10_date", "prs", "sbp", "smoke", "bmi", "diabetes", "lipid_lowering", "tot_chol", "hdl_chol", "startfollowup", "endfollowup")]
  total.newI <- total.newI[!duplicated(total.newI$eid),]
  # survival analysis of individuals with specific ICD10 code
  mod_prs_i <- coxph(Surv(endfollowup - startfollowup,mi) ~ age + scale(prs) + sex + sbp + smoke + bmi + diabetes + lipid_lowering + tot_chol + hdl_chol, data = total.newI)
  hr_i <- coef(summary(mod_prs_i))["scale(prs)", "exp(coef)"]
  hr_i_min <- exp(coef(summary(mod_prs_i))["scale(prs)", "coef"] - 1.96*coef(summary(mod_prs_i))["scale(prs)", "se(coef)"])
  hr_i_max <- exp(coef(summary(mod_prs_i))["scale(prs)", "coef"] + 1.96*coef(summary(mod_prs_i))["scale(prs)", "se(coef)"])

  # create datatset of individuals without specific ICD10 code
  total.newNI <- total.new[!(total.new$eid %in% ids),]
  total.newNI$mi <- ifelse(is.na(total.newNI$mi_date), 0, 1)
  total.newNI <- total.newNI[,c("eid", "sex", "age", "death", "death_date", "mi", "mi_date", "icd10", "icd10_date", "prs", "sbp", "smoke", "bmi", "diabetes", "lipid_lowering", "tot_chol", "hdl_chol", "startfollowup", "endfollowup")]
  total.newNI <- total.newNI[!duplicated(total.newNI$eid),]
  # survival analysis of individuals without specific ICD10 code
  mod_prs_ni <- coxph(Surv(endfollowup - startfollowup,mi) ~ age + scale(prs) + sex + sbp + smoke + bmi + diabetes + lipid_lowering + tot_chol + hdl_chol, data = total.newNI)
  hr_ni <- coef(summary(mod_prs_ni))["scale(prs)", "exp(coef)"]
  hr_ni_min <- exp(coef(summary(mod_prs_ni))["scale(prs)", "coef"] - 1.96*coef(summary(mod_prs_ni))["scale(prs)", "se(coef)"])
  hr_ni_max <- exp(coef(summary(mod_prs_ni))["scale(prs)", "coef"] + 1.96*coef(summary(mod_prs_ni))["scale(prs)", "se(coef)"])

  numb <- rbind(cbind(icd.all[i], hr_i, hr_i_min, hr_i_max), cbind(paste("No", icd.all[i]), hr_ni, hr_ni_min, hr_ni_max))
  sub.pop.table <- rbind(sub.pop.table, numb)
  print(i)
}
sub.pop.table <- data.frame(sub.pop.table)
rownames(sub.pop.table) <- 1:nrow(sub.pop.table)
colnames(sub.pop.table) <- c("group", "hr", "hrmin", "hrmax")
sub.pop.table$group <- factor(sub.pop.table$group, levels = sub.pop.table$group)
sub.pop.table$hr <- as.numeric(as.character(sub.pop.table$hr))
sub.pop.table$hrmin <- as.numeric(as.character(sub.pop.table$hrmin))
sub.pop.table$hrmax <- as.numeric(as.character(sub.pop.table$hrmax))
write.table(sub.pop.table, file = 'sensitivity_subset_interaction.tsv', sep = "\t", row.names = TRUE, col.names = TRUE)















##################################################
# Analysis of continuous PRS for population with and without previous revelant comorbidities
##################################################

# create list of people with no previous comorbidities and previous comorbidities
icd.ids <- unique(total.new$eid[total.new$icd10 %in% icd.new])
no.icd.ids <- unique(total.new$eid[!(total.new$eid %in% icd.ids)])
# subset population with no previous comorbidities 
total.newNA <- total.new[total.new$eid %in% no.icd.ids,]
total.newNA$mi <- ifelse(is.na(total.newNA$mi_date), 0, 1)
total.newNA <- total.newNA[!duplicated(total.newNA$eid), c("eid", "age", "sex", "birth_date", "death", "death_date", "mi", "mi_date", "icd10", "icd10_date", "prs", "sbp", "smoke", "bmi", "diabetes", "lipid_lowering", "tot_chol", "hdl_chol", "startfollowup", "endfollowup")]
total.newNA <- total.newNA[!is.na(total.newNA$prs) & !is.na(total.newNA$sbp) & !is.na(total.newNA$smoke) & !is.na(total.newNA$bmi) & !is.na(total.newNA$diabetes) & !is.na(total.newNA$lipid_lowering) & !is.na(total.newNA$tot_chol) & !is.na(total.newNA$hdl_chol),]
# check association between PRS and MI for individuals with no previous comorbidities
mod_prs_no_comorb <- coxph(Surv(endfollowup - startfollowup,mi) ~ age + scale(prs) + sex + sbp + smoke + bmi + diabetes + lipid_lowering + tot_chol + hdl_chol, data = total.newNA)
cindex_prs_no_comorb <- survConcordance(Surv(endfollowup - startfollowup, mi) ~ predict(mod_prs_no_comorb), data = total.newNA)$concordance
err_prs_no_comorb <- as.numeric(survConcordance(Surv(endfollowup - startfollowup, mi) ~ predict(mod_prs_no_comorb), data = total.newNA)$std.err)
mod_prs_no_comorb_base <- coxph(Surv(endfollowup - startfollowup,mi) ~ age + sex + sbp + smoke + bmi + diabetes + lipid_lowering + tot_chol + hdl_chol, data = total.newNA)
cindex_prs_no_comorb_base <- survConcordance(Surv(endfollowup - startfollowup, mi) ~ predict(mod_prs_no_comorb_base), data = total.newNA)$concordance
err_prs_no_comorb_base <- as.numeric(survConcordance(Surv(endfollowup - startfollowup, mi) ~ predict(mod_prs_no_comorb_base), data = total.newNA)$std.err)
no_comorb <- cbind(group = "None", cindex = cindex_prs_no_comorb, error = err_prs_no_comorb, cindex.base = cindex_prs_no_comorb_base, error.base = err_prs_no_comorb_base)
# subset population with previous comorbidities
total.newICD <- total.new[total.new$eid %in% icd.ids,]
total.newICD$mi <- ifelse(is.na(total.newICD$mi_date), 0, 1)
total.newICD <- total.newICD[!duplicated(total.newICD$eid), c("eid", "age", "sex", "birth_date", "death", "death_date", "mi", "mi_date", "icd10", "icd10_date", "prs", "sbp", "smoke", "bmi", "diabetes", "lipid_lowering", "tot_chol", "hdl_chol", "startfollowup", "endfollowup")]
total.newICD <- total.newICD[!is.na(total.newICD$prs) & !is.na(total.newICD$sbp) & !is.na(total.newICD$smoke) & !is.na(total.newICD$bmi) & !is.na(total.newICD$diabetes) & !is.na(total.newICD$lipid_lowering) & !is.na(total.newICD$tot_chol) & !is.na(total.newICD$hdl_chol),]
# check association between PRS and MI for individuals with previous comorbidities
mod_prs_comorb <- coxph(Surv(endfollowup - startfollowup,mi) ~ age + scale(prs) + sex + sbp + smoke + bmi + diabetes + lipid_lowering + tot_chol + hdl_chol, data = total.newICD)
cindex_prs_comorb <- survConcordance(Surv(endfollowup - startfollowup, mi) ~ predict(mod_prs_comorb), data = total.newICD)$concordance
err_prs_comorb <- as.numeric(survConcordance(Surv(endfollowup - startfollowup, mi) ~ predict(mod_prs_comorb), data = total.newICD)$std.err)
mod_prs_comorb_base <- coxph(Surv(endfollowup - startfollowup,mi) ~ age + sex + sbp + smoke + bmi + diabetes + lipid_lowering + tot_chol + hdl_chol, data = total.newICD)
cindex_prs_comorb_base <- survConcordance(Surv(endfollowup - startfollowup, mi) ~ predict(mod_prs_comorb_base), data = total.newICD)$concordance
err_prs_comorb_base <- as.numeric(survConcordance(Surv(endfollowup - startfollowup, mi) ~ predict(mod_prs_comorb_base), data = total.newICD)$std.err)
comorb <- cbind(group = "All", cindex = cindex_prs_comorb, error = err_prs_comorb, cindex.base = cindex_prs_comorb_base, error.base = err_prs_comorb_base)
# finalize data frame
cindex.comorb.nocomorb <- rbind(no_comorb, comorb)
cindex.comorb.nocomorb <- data.frame(cindex.comorb.nocomorb)
cindex.comorb.nocomorb$cindex <- as.numeric(as.character(cindex.comorb.nocomorb$cindex))
cindex.comorb.nocomorb$error <- as.numeric(as.character(cindex.comorb.nocomorb$error))
cindex.comorb.nocomorb$cindex.base <- as.numeric(as.character(cindex.comorb.nocomorb$cindex.base))
cindex.comorb.nocomorb$error.base <- as.numeric(as.character(cindex.comorb.nocomorb$error.base))
cindex.comorb.nocomorb$diff <- as.numeric(cindex.comorb.nocomorb$cindex - cindex.comorb.nocomorb$cindex.base)
cindex.comorb.nocomorb$diff.error <- cindex.comorb.nocomorb$error - cindex.comorb.nocomorb$error.base
write.table(cindex.comorb.nocomorb, file = 'sensitivity_cont_prs.tsv', sep = "\t", row.names = TRUE, col.names = TRUE)















##################################################
# Analysis of dichotomous PRS for population with and without previous relevant comorbidities
##################################################

# analysis for population without previous relevant comorbidities
total.newPRS <- total.newNA
cutoffs <- c("top1", "top2", "top5", "top10")
# create dichotomous variable for different percentiles of PRS (e.g., for top 1%, 1 = PRS higher than 99th percentile and 0 = PRS lower than 99th percentile)
total.newPRS$top1 <- ifelse(total.newPRS$prs >= quantile(total.newPRS$prs, 0.99, na.rm = TRUE), 1, 0)
total.newPRS$top2 <- ifelse(total.newPRS$prs >= quantile(total.newPRS$prs, 0.98, na.rm = TRUE), 1, 0)
total.newPRS$top5 <- ifelse(total.newPRS$prs >= quantile(total.newPRS$prs, 0.95, na.rm = TRUE), 1, 0)
total.newPRS$top10 <- ifelse(total.newPRS$prs >= quantile(total.newPRS$prs, 0.90, na.rm = TRUE), 1, 0)
# survival analysis for different cutoffs of dichotomous PRS
cindex.dich.prs <- NULL
for (i in 1:length(cutoffs)) {
  mod_dich_prs <- coxph(Surv(endfollowup - startfollowup, mi) ~ age + total.newPRS[, cutoffs[i]] + sex + sbp + smoke + bmi + diabetes + lipid_lowering + tot_chol + hdl_chol, data = total.newPRS)
  numb <- cbind(group = cutoffs[i], hr = coef(summary(mod_dich_prs))["total.newPRS[, cutoffs[i]]", "exp(coef)"], 
    hrmin = exp(coef(summary(mod_dich_prs))["total.newPRS[, cutoffs[i]]", "coef"] - 1.96*coef(summary(mod_dich_prs))["total.newPRS[, cutoffs[i]]", "se(coef)"]), 
    hrmax = exp(coef(summary(mod_dich_prs))["total.newPRS[, cutoffs[i]]", "coef"] + 1.96*coef(summary(mod_dich_prs))["total.newPRS[, cutoffs[i]]", "se(coef)"]))
  cindex.dich.prs <- rbind(cindex.dich.prs, numb)
  print(i)
}
rownames(cindex.dich.prs) <- 1:nrow(cindex.dich.prs)
cindex.dich.prs <- data.frame(cindex.dich.prs)
cindex.dich.prs$hr <- as.numeric(as.character(cindex.dich.prs$hr))
cindex.dich.prs$hrmin <- as.numeric(as.character(cindex.dich.prs$hrmin))
cindex.dich.prs$hrmax <- as.numeric(as.character(cindex.dich.prs$hrmax))
cindex.dich.prs1 <- cindex.dich.prs
cindex.dich.prs1$group <- paste("No Prev Comorb", cindex.dich.prs1$group)

# analysis for population with previous relevant comorbidities
total.newPRS <- total.newICD
cutoffs <- c("top1", "top2", "top5", "top10")
# create dichotomous variable for different percentiles of PRS (e.g., for top 1%, 1 = PRS higher than 99th percentile and 0 = PRS lower than 99th percentile)
total.newPRS$top1 <- ifelse(total.newPRS$prs >= quantile(total.newPRS$prs, 0.99, na.rm = TRUE), 1, 0)
total.newPRS$top2 <- ifelse(total.newPRS$prs >= quantile(total.newPRS$prs, 0.98, na.rm = TRUE), 1, 0)
total.newPRS$top5 <- ifelse(total.newPRS$prs >= quantile(total.newPRS$prs, 0.95, na.rm = TRUE), 1, 0)
total.newPRS$top10 <- ifelse(total.newPRS$prs >= quantile(total.newPRS$prs, 0.90, na.rm = TRUE), 1, 0)
# survival analysis for different cutoffs of dichotomous PRS
cindex.dich.prs <- NULL
for (i in 1:length(cutoffs)) {
  mod_dich_prs <- coxph(Surv(endfollowup - startfollowup, mi) ~ age + total.newPRS[, cutoffs[i]] + sex + sbp + smoke + bmi + diabetes + lipid_lowering + tot_chol + hdl_chol, data = total.newPRS)
  numb <- cbind(group = cutoffs[i], hr = coef(summary(mod_dich_prs))["total.newPRS[, cutoffs[i]]", "exp(coef)"], 
    hrmin = exp(coef(summary(mod_dich_prs))["total.newPRS[, cutoffs[i]]", "coef"] - 1.96*coef(summary(mod_dich_prs))["total.newPRS[, cutoffs[i]]", "se(coef)"]), 
    hrmax = exp(coef(summary(mod_dich_prs))["total.newPRS[, cutoffs[i]]", "coef"] + 1.96*coef(summary(mod_dich_prs))["total.newPRS[, cutoffs[i]]", "se(coef)"]))
  cindex.dich.prs <- rbind(cindex.dich.prs, numb)
  print(i)
}
rownames(cindex.dich.prs) <- 1:nrow(cindex.dich.prs)
cindex.dich.prs <- data.frame(cindex.dich.prs)
cindex.dich.prs$hr <- as.numeric(as.character(cindex.dich.prs$hr))
cindex.dich.prs$hrmin <- as.numeric(as.character(cindex.dich.prs$hrmin))
cindex.dich.prs$hrmax <- as.numeric(as.character(cindex.dich.prs$hrmax))
cindex.dich.prs2 <- cindex.dich.prs
cindex.dich.prs2$group <- paste("Prev Comorb", cindex.dich.prs2$group)

cindex.dich.prs <- rbind(cindex.dich.prs1, cindex.dich.prs2)
cindex.dich.prs <- cindex.dich.prs[c(1, 5, 2, 6, 3, 7, 4, 8),]
cindex.dich.prs$group <- factor(cindex.dich.prs$group, levels = cindex.dich.prs$group)
write.table(cindex.dich.prs, file = 'sensitivity_dich_prs.tsv', sep = "\t", row.names = TRUE, col.names = TRUE)















##################################################
# ROC analysis of PRS with ICD10 codes and dichotomous PRS
##################################################

# create dataset for ROC analysis
total.newSS <- total.new[which(total.new$icd10 %in% icd.new), c("eid", "icd10", "mi", "prs")]
# create dichotomous variable for ICD10 codes (e.g., 1 = individual has I20 and 0 = individual does not have I20)
for (i in 1:length(icd.new)) {
  temp_ids <- total.newSS$eid[total.newSS$icd10 == icd.new[i]]
  total.newSS[,ncol(total.newSS)+1] <- ifelse(total.newSS$eid %in% temp_ids, 1, 0)
  colnames(total.newSS)[ncol(total.newSS)] <- icd.new[i]
  print(i)
}
total.newSS <- total.newSS[!duplicated(total.newSS$eid),]
# create dichotomous variable for PRS at top 1%, 2%, 5%, and 10% of PRS
cutoffs <- c("top1", "top2", "top5", "top10")
total.newSS$top1 <- ifelse(total.newSS$prs > quantile(total.newSS$prs, 0.99, na.rm = TRUE), 1, 0)
total.newSS$top2 <- ifelse(total.newSS$prs > quantile(total.newSS$prs, 0.98, na.rm = TRUE), 1, 0)
total.newSS$top5 <- ifelse(total.newSS$prs > quantile(total.newSS$prs, 0.95, na.rm = TRUE), 1, 0)
total.newSS$top10 <- ifelse(total.newSS$prs > quantile(total.newSS$prs, 0.90, na.rm = TRUE), 1, 0)

# loop across all ICD10 codes and find sensitivity and specificity
pop.table <- NULL
for (i in 1:length(icd.new)) {
  sens <- confusionMatrix(data = factor(total.newSS[,icd.new[i]]), reference = factor(total.newSS$mi), positive = "1")$byClass["Sensitivity"]
  spec <- confusionMatrix(data = factor(total.newSS[,icd.new[i]]), reference = factor(total.newSS$mi), positive = "1")$byClass["Specificity"]
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
rownames(pop.table) <- c(icd.new, cutoffs)
colnames(pop.table) <- c("sens", "spec")
pop.table <- data.frame(pop.table)
pop.table$fpr <- 1-as.numeric(pop.table$spec)
pop.table$label <- rownames(pop.table)
pop.table$type <- ifelse(pop.table$label %in% icd.new, "icd10", "prs")
write.table(pop.table, file = 'sensitivity_icd_sens_spec.tsv', sep = "\t", row.names = TRUE, col.names = TRUE)

# find AUC of dichotomous PRS
auc.prs <- as.numeric(roc(total.newSS$mi, total.newSS$prs)$auc)
auc.prs.top1 <- as.numeric(roc(total.newSS$mi, total.newSS$top1)$auc)
auc.prs.top2 <- as.numeric(roc(total.newSS$mi, total.newSS$top2)$auc)
auc.prs.top5 <- as.numeric(roc(total.newSS$mi, total.newSS$top5)$auc)
auc.prs.top10 <- as.numeric(roc(total.newSS$mi, total.newSS$top10)$auc)
cindex_dich <- data.frame(cutoffs = c("all", cutoffs), cindex = c(auc.prs, auc.prs.top1, auc.prs.top2, auc.prs.top5, auc.prs.top10))
write.table(cindex_dich, file = 'sensitivity_prs_auc.tsv', sep = "\t", row.names = TRUE, col.names = TRUE)

# find sensitivities and specificities for PRS 
roc.table <- cbind(roc(total.newSS$mi, total.newSS$prs)$sensitivities, roc(total.newSS$mi, total.newSS$prs)$specificities)
roc.table <- data.frame(roc.table)
colnames(roc.table) <- c("sens", "spec")
roc.table$fpr <- 1-roc.table$spec
write.table(roc.table, file = 'sensitivity_prs_roc.tsv', sep = "\t", row.names = TRUE, col.names = TRUE)