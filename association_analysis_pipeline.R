##################################################
# MYOCARDIAL INFARCTION ASSOCIATION ANALYSIS
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

total <- read.table('data.tsv', sep = "\t", header = TRUE, stringsAsFactors = F)
# Columns in data should follow following format:
# eid            - boolean for ID
# age            - int/num for age
# sex            - boolean for sex
# birth_date     - date for day of birth
# death          - boolean for death
# death_date     - date for day of death
# prs            - int/num for polygenic risk score
# sbp            - int/num for systolic blood pressure
# smoke          - boolean for smoking status
# bmi            - int/num for body mass index
# diabetes       - boolean for diabetes status
# lipid_lowering - boolean for lipid-lowering medications
# tot_chol       - int/num for total cholesterol
# hdl_chol       - int/num for HDL cholesterol
# ldl_chol       - int/num for LDL cholesterol
# mi             - boolean for myocardial infarction
# mi_date        - date for day of myocardial infarction
# ac_date        - date for baseline
# icd10          - chr for ICD10 event
# icd10_date     - date for day of ICD10 event
# Rows in data should be unique ICD10 codes for different individuals
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
total$ldl_chol <- as.integer(total$ldl_chol)
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
study.start <- as.Date("1998-01-01", "%Y-%m-%d")
study.end <- as.Date("2015-03-01", "%Y-%m-%d")
total$startfollowup <- total$ac_date
total$endfollowup <- ifelse(!is.na(total$death_date) & is.na(total$mi_date), as.Date(total$death_date, "%Y-%m-%d"), ifelse(!is.na(total$mi_date), as.Date(total$mi_date, "%Y-%m-%d"), study.end))
class(total$startfollowup) <- "Date"
class(total$endfollowup) <- "Date"

# exclude ICD10 codes that are not related to disease (e.g., car crash) 
total.new <- total[,c("eid", "age", "sex", "birth_date", "death", "death_date", "prs", "sbp", "smoke", "bmi", "diabetes", "lipid_lowering", "tot_chol", "hdl_chol", "ldl_chol", "mi", "mi_date", "ac_date", "icd10", "icd10_date", "startfollowup", "endfollowup")]
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

# ASSOCIATION ANALYSIS: keep ICD10 codes both before and after baseline
total.new1 <- total.new
total.new1$icd10[total.new1$icd10_date > study.end | total.new1$icd10_date < study.start] <- NA
total.new1$icd10_date[total.new1$icd10_date > study.end | total.new1$icd10_date < study.start] <- NA
total.new1$icd10[(total.new1$mi_date - total.new1$icd10_date) < 8] <- NA
total.new1$icd10_date[(total.new1$mi_date - total.new1$icd10_date) < 8] <- NA
nrow(total.new1)
#_____ rows
length(unique(total.new1$eid))
#_____ individuals
length(unique(total.new1$icd10))
#_____ ICD codes
length(unique(total.new1[which(!is.na(total.new1$mi_date)),]$eid))
#_____ individuals with MI





















##################################################
# Survival analysis for PRS
total.new <- total.new1 # for association analysis
##################################################

# create datasets for survival analysis
total.newMOD <- total.new[!duplicated(total.new$eid), c("eid", "age", "sex", "birth_date", "death", "death_date", "prs", "sbp", "smoke", "bmi", "diabetes", "lipid_lowering", "tot_chol", "hdl_chol", "ldl_chol", "mi", "mi_date", "ac_date", "startfollowup", "endfollowup")]
total.newMOD <- total.newMOD[!is.na(total.newMOD$prs) & !is.na(total.newMOD$sbp) & !is.na(total.newMOD$smoke) & !is.na(total.newMOD$bmi) & !is.na(total.newMOD$diabetes) & !is.na(total.newMOD$lipid_lowering) & !is.na(total.newMOD$tot_chol) & !is.na(total.newMOD$hdl_chol) & !is.na(total.newMOD$ldl_chol),]
# survival analysis models and c-index for PRS
mod_risk_factor <- coxph(Surv(endfollowup - startfollowup,mi) ~ age + sex + sbp + smoke + bmi + diabetes + lipid_lowering + tot_chol + hdl_chol + ldl_chol, data = total.newMOD)
cindex_risk_factor <- survConcordance(Surv(endfollowup - startfollowup, mi) ~ predict(mod_risk_factor), data = total.newMOD)$concordance
mod_prs_risk_factor <- coxph(Surv(endfollowup - startfollowup, mi) ~ age + scale(prs) + sex + sbp + smoke + bmi + diabetes + lipid_lowering + tot_chol + hdl_chol + ldl_chol, data = total.newMOD)
cindex_prs_risk_factor <- survConcordance(Surv(endfollowup - startfollowup, mi) ~ predict(mod_prs_risk_factor), data = total.newMOD)$concordance

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
total.newU <- total.new[!duplicated(total.new$eid), c("eid", "age", "sex", "birth_date", "death", "death_date", "prs", "sbp", "smoke", "bmi", "diabetes", "lipid_lowering", "tot_chol", "hdl_chol", "ldl_chol", "ac_date", "startfollowup", "endfollowup")]

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
  findata <- findata[!is.na(findata$prs) & !is.na(findata$sbp) & !is.na(findata$smoke) & !is.na(findata$bmi) & !is.na(findata$diabetes) & !is.na(findata$lipid_lowering) & !is.na(findata$tot_chol) & !is.na(findata$hdl_chol) & !is.na(findata$ldl_chol),]

  # base model
  mod <- tryCatch({ coxph(Surv(tstart, tstop, outP) ~ trt, data = findata) }, error = function(w) { NA })
  cindex_icd_base <- survConcordance(Surv(tstart, tstop, outP) ~ predict(mod), data = findata)$concordance
  # model with age and sex
  mod1 <- tryCatch({ coxph(Surv(tstart, tstop, outP) ~ age + trt + sex, data = findata) }, error = function(w) { NA })
  cindex_icd <- survConcordance(Surv(tstart, tstop, outP) ~ predict(mod1), data = findata)$concordance
  # model with risk factors
  mod2 <- tryCatch({ coxph(Surv(tstart, tstop, outP) ~ age + trt + sex + sbp + smoke + bmi + diabetes + lipid_lowering + tot_chol + hdl_chol + ldl_chol, data = findata) }, error = function(w) { NA })
  cindex_icd_risk_factor <- survConcordance(Surv(tstart, tstop, outP) ~ predict(mod2), data = findata)$concordance
  # base model of interaction
  mod3 <- tryCatch({ coxph(Surv(tstart, tstop, outP) ~ trt*prs, data = findata) }, error = function(w) { NA })
  cindex_interaction_base <- survConcordance(Surv(tstart, tstop, outP) ~ predict(mod3), data = findata)$concordance
  # base model of interaction with age and sex
  mod4 <- tryCatch({ coxph(Surv(tstart, tstop, outP) ~ age + trt*prs + sex, data = findata) }, error = function(w) { NA })
  cindex_interaction <- survConcordance(Surv(tstart, tstop, outP) ~ predict(mod4), data = findata)$concordance
  # base model of interaction with risk factors
  mod5 <- tryCatch({ coxph(Surv(tstart, tstop, outP) ~ age + trt*prs + sex + sbp + smoke + bmi + diabetes + lipid_lowering + tot_chol + hdl_chol + ldl_chol, data = findata) }, error = function(w) { NA })
  cindex_interaction_risk_factor <- survConcordance(Surv(tstart, tstop, outP) ~ predict(mod5), data = findata)$concordance

  res <- rbind(res, c(coef(summary(mod))["trt", c("coef", "z")], coef(summary(mod1))["trt", c("coef", "z")], coef(summary(mod2))["trt", c("coef", "z")], coef(summary(mod3))["trt:prs", c("coef", "z")], coef(summary(mod4))["trt:prs", c("coef", "z")], coef(summary(mod5))["trt:prs", c("coef", "z")], 
    mean_day_dist, sd_day_dist, cindex_icd_base, cindex_icd, cindex_icd_risk_factor, cindex_interaction_base, cindex_interaction, cindex_interaction_risk_factor))

  print(i)

}

resdf = as.data.frame(res)
resdf = cbind(icd10 = icd.all, resdf)
rownames(resdf) <- 1:nrow(resdf)
colnames(resdf) <- c("ICD10", "beta_base", "z_base", "beta_icd", "z_icd", "beta_risk_factor", "z_risk_factor", "beta_interaction_base", "z_interaction_base", "beta_interaction", "z_interaction", "beta_interaction_risk_factor", "z_interaction_risk_factor", 
  "mean_day_dist", "sd_day_dist", "cindex_icd_base", "cindex_icd", "cindex_icd_risk_factor", "cindex_interaction_base", "cindex_interaction", "cindex_interaction_risk_factor")
write.table(resdf, file = paste0(path, 'association_analysis.tsv'), sep = "\t", row.names = FALSE, col.names = TRUE, quote=F)



















##################################################
# Plot c-index, hazard ratio, and p-value
##################################################

# load results of survival analysis
resdf <- read.table(file = paste0(path, 'association_analysis.tsv'), sep = "\t", header = TRUE, stringsAsFactors = F)
# calculate hazard ratio and p-value from results of survival analysis
resdf$logp_base <- -log10(2*pnorm(-abs(resdf$z_base)))
resdf$logp_base[resdf$logp_base == "Inf"] <- 300
resdf$logp_icd <- -log10(2*pnorm(-abs(resdf$z_icd)))
resdf$logp_icd[resdf$logp_icd == "Inf"] <- 300
resdf$logp_risk_factor <- -log10(2*pnorm(-abs(resdf$z_risk_factor)))
resdf$logp_risk_factor[resdf$logp_risk_factor == "Inf"] <- 300
resdf$logp_interaction_base <- -log10(2*pnorm(-abs(resdf$z_interaction_base)))
resdf$logp_interaction_base[resdf$logp_interaction_base == "Inf"] <- 300
resdf$logp_interaction <- -log10(2*pnorm(-abs(resdf$z_interaction)))
resdf$logp_interaction[resdf$logp_interaction == "Inf"] <- 300
resdf$logp_interaction_risk_factor <- -log10(2*pnorm(-abs(resdf$z_interaction_risk_factor)))
resdf$logp_interaction_risk_factor[resdf$logp_interaction_risk_factor == "Inf"] <- 300
resdf$hr_base <- exp(resdf$beta_base)
resdf$hr_icd <- exp(resdf$beta_icd)
resdf$hr_risk_factor <- exp(resdf$beta_risk_factor)
resdf$ICD10label <- resdf$ICD10
resdf$ICD10label[resdf$logp_icd < -log10(0.05 / length(icd.all))] <- ""
resdf.new <- resdf[resdf$ICD10 %in% resdf$ICD10[resdf$logp_icd > -log10(0.05 / length(icd.all)) | resdf$logp_risk_factor > -log10(0.05 / length(icd.all))],]
icd.new <- resdf.new$ICD10
# FIGURE 1A: plot association of MI and ICD10 codes (p-value)
ggplot() +
  geom_point(data = resdf, mapping = aes(x = ICD10, y = logp_risk_factor, size = hr_risk_factor, color = "Adjusted")) +
  geom_text_repel(data = resdf, mapping = aes(x = ICD10, y = logp_risk_factor, label = resdf$ICD10label), size = 3, segment.alpha = 0.5) +
  geom_point(data = resdf, mapping = aes(x = ICD10, y = logp_icd, size = hr_icd, color = "Unadjusted")) +
  geom_text_repel(data = resdf, mapping = aes(x = ICD10, y = logp_icd, label = resdf$ICD10label), size = 3, segment.alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05 / length(icd.all)), size = 1.5, color = "red") + 
  annotate(geom = "text", x = length(icd.all)/2, y = 3, label = "Significance Level", color = "red") +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = "ICD10 Code", y = "-log10 Transformed P-Value") 
# FIGURE 1B: plot predictive ability of ICD10 codes for MI (c-index)
ggplot() + 
  geom_point(data = resdf.new, mapping = aes(ICD10, as.numeric(cindex_icd_risk_factor - cindex_risk_factor))) +
  geom_text_repel(data = resdf.new, mapping = aes(ICD10, as.numeric(cindex_icd_risk_factor - cindex_risk_factor), label = resdf.new$ICD10label), size = 3, segment.alpha = 0.5) +
  geom_point(mapping = aes(x = length(icd.new)/2, y = as.numeric(cindex_prs_risk_factor - cindex_risk_factor)), color = "#F8766D", shape = 18, size = 10) + 
  annotate(geom = "text", x = length(icd.new)/2, y = as.numeric(cindex_prs_risk_factor - cindex_risk_factor) + 0.0005, label = "Polygenic Risk Score") +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = "ICD10 Code", y = "Difference in C-index Compared to Baseline") 
# FIGURE 2A: plot association of MI and interaction between ICD10 codes and PRS (p-value)
ggplot(data = resdf.new, mapping = aes(ICD10, logp_interaction_base)) +
  geom_point() +
  geom_text_repel(data = resdf.new, mapping = aes(label = resdf.new$ICD10label), size = 3, segment.alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05 / length(icd.all)), size = 1.5, color = "red") +
  annotate(geom = "text", x = length(icd.new)/2, y = 3.4, label = "Significance Level", color = "red") +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = "ICD10 Code", y = "-log10 Transformed P-Value") 




















##################################################
# Analysis of interaction for population subsets of ICD10 codes
##################################################

# calculate and plot c-index for population subsets (e.g., individuals with only I20) using model of PRS
sub.pop.table <- NULL
# CHANGE: replace main.icd with important ICD10 codes (e.g., ICD10 codes with large c-indices)
main.icd <- c("I20", "I24", "R07")
for (i in 1:length(main.icd)) {

  # create list of individuals with specific ICD10 code
  ids <- unique(total.new$eid[total.new$icd10 == main.icd[i]])

  # create dataset of individuals with specific ICD10 code
  total.newI <- total.new[total.new$eid %in% ids,]
  total.newI$mi <- ifelse(is.na(total.newI$mi_date), 0, 1)
  total.newI <- total.newI[which(total.newI$icd10 == main.icd[i]), c("eid", "sex", "age", "death", "death_date", "mi", "mi_date", "icd10", "icd10_date", "prs", "sbp", "smoke", "bmi", "diabetes", "lipid_lowering", "tot_chol", "hdl_chol", "ldl_chol", "startfollowup", "endfollowup")]
  total.newI <- total.newI[!duplicated(total.newI$eid),]
  # survival analysis of individuals with specific ICD10 code
  mod_prs_i <- coxph(Surv(endfollowup - startfollowup,mi) ~ age + scale(prs) + sex + sbp + smoke + bmi + diabetes + lipid_lowering + tot_chol + hdl_chol + ldl_chol, data = total.newI)
  hr_i <- coef(summary(mod_prs_i))["scale(prs)", "exp(coef)"]
  hr_i_min <- exp(coef(summary(mod_prs_i))["scale(prs)", "coef"] - 1.96*coef(summary(mod_prs_i))["scale(prs)", "se(coef)"])
  hr_i_max <- exp(coef(summary(mod_prs_i))["scale(prs)", "coef"] + 1.96*coef(summary(mod_prs_i))["scale(prs)", "se(coef)"])

  # create datatset of individuals without specific ICD10 code
  total.newNI <- total.new[!(total.new$eid %in% ids),]
  total.newNI$mi <- ifelse(is.na(total.newNI$mi_date), 0, 1)
  total.newNI <- total.newNI[,c("eid", "sex", "age", "death", "death_date", "mi", "mi_date", "icd10", "icd10_date", "prs", "sbp", "smoke", "bmi", "diabetes", "lipid_lowering", "tot_chol", "hdl_chol", "ldl_chol", "startfollowup", "endfollowup")]
  total.newNI <- total.newNI[!duplicated(total.newNI$eid),]
  # survival analysis of individuals without specific ICD10 code
  mod_prs_ni <- coxph(Surv(endfollowup - startfollowup,mi) ~ age + scale(prs) + sex + sbp + smoke + bmi + diabetes + lipid_lowering + tot_chol + hdl_chol + ldl_chol, data = total.newNI)
  hr_ni <- coef(summary(mod_prs_ni))["scale(prs)", "exp(coef)"]
  hr_ni_min <- exp(coef(summary(mod_prs_ni))["scale(prs)", "coef"] - 1.96*coef(summary(mod_prs_ni))["scale(prs)", "se(coef)"])
  hr_ni_max <- exp(coef(summary(mod_prs_ni))["scale(prs)", "coef"] + 1.96*coef(summary(mod_prs_ni))["scale(prs)", "se(coef)"])

  numb <- rbind(cbind(main.icd[i], hr_i, hr_i_min, hr_i_max), cbind(paste("No", main.icd[i]), hr_ni, hr_ni_min, hr_ni_max))
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

# FIGURE 2A: association of dichotomous PRS and MI in population with and without previous relevant comorbidites (hazard ratio)
ggplot(data = sub.pop.table, mapping = aes(group, hr)) +
  geom_point() +
  geom_errorbar(aes(ymin = hrmin, ymax = hrmax), width = 0.2, position = position_dodge(1)) + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x = "Group", y = "Hazard Ratio")




















##################################################
# Analysis of continuous PRS for population with and without previous revelant comorbidities
##################################################

# create list of people with no previous comorbidities and previous comorbidities
icd.ids <- unique(total.new$eid[total.new$icd10 %in% icd.new])
no.icd.ids <- unique(total.new$eid[!(total.new$eid %in% icd.ids)])
# subset population with no previous comorbidities 
total.newNA <- total.new[total.new$eid %in% no.icd.ids,]
total.newNA$mi <- ifelse(is.na(total.newNA$mi_date), 0, 1)
total.newNA <- total.newNA[!duplicated(total.newNA$eid), c("eid", "age", "sex", "birth_date", "death", "death_date", "mi", "mi_date", "icd10", "icd10_date", "prs", "sbp", "smoke", "bmi", "diabetes", "lipid_lowering", "tot_chol", "hdl_chol", "ldl_chol", "startfollowup", "endfollowup")]
total.newNA <- total.newNA[!is.na(total.newNA$prs) & !is.na(total.newNA$sbp) & !is.na(total.newNA$smoke) & !is.na(total.newNA$bmi) & !is.na(total.newNA$diabetes) & !is.na(total.newNA$lipid_lowering) & !is.na(total.newNA$tot_chol) & !is.na(total.newNA$hdl_chol) & !is.na(total.newNA$ldl_chol),]
# check association between PRS and MI for individuals with no previous comorbidities
mod_prs_no_comorb <- coxph(Surv(endfollowup - startfollowup,mi) ~ age + scale(prs) + sex + sbp + smoke + bmi + diabetes + lipid_lowering + tot_chol + hdl_chol + ldl_chol, data = total.newNA)
cindex_prs_no_comorb <- survConcordance(Surv(endfollowup - startfollowup, mi) ~ predict(mod_prs_no_comorb), data = total.newNA)$concordance
err_prs_no_comorb <- as.numeric(survConcordance(Surv(endfollowup - startfollowup, mi) ~ predict(mod_prs_no_comorb), data = total.newNA)$std.err)
mod_prs_no_comorb_base <- coxph(Surv(endfollowup - startfollowup,mi) ~ age + sex + sbp + smoke + bmi + diabetes + lipid_lowering + tot_chol + hdl_chol + ldl_chol, data = total.newNA)
cindex_prs_no_comorb_base <- survConcordance(Surv(endfollowup - startfollowup, mi) ~ predict(mod_prs_no_comorb_base), data = total.newNA)$concordance
err_prs_no_comorb_base <- as.numeric(survConcordance(Surv(endfollowup - startfollowup, mi) ~ predict(mod_prs_no_comorb_base), data = total.newNA)$std.err)
no_comorb <- cbind(group = "None", cindex = cindex_prs_no_comorb, error = err_prs_no_comorb, cindex.base = cindex_prs_no_comorb_base, error.base = err_prs_no_comorb_base)
# subset population with previous comorbidities
total.newICD <- total.new[total.new$eid %in% icd.ids,]
total.newICD$mi <- ifelse(is.na(total.newICD$mi_date), 0, 1)
total.newICD <- total.newICD[!duplicated(total.newICD$eid), c("eid", "age", "sex", "birth_date", "death", "death_date", "mi", "mi_date", "icd10", "icd10_date", "prs", "sbp", "smoke", "bmi", "diabetes", "lipid_lowering", "tot_chol", "hdl_chol", "ldl_chol", "startfollowup", "endfollowup")]
total.newICD <- total.newICD[!is.na(total.newICD$prs) & !is.na(total.newICD$sbp) & !is.na(total.newICD$smoke) & !is.na(total.newICD$bmi) & !is.na(total.newICD$diabetes) & !is.na(total.newICD$lipid_lowering) & !is.na(total.newICD$tot_chol) & !is.na(total.newICD$hdl_chol) & !is.na(total.newICD$ldl_chol),]
# check association between PRS and MI for individuals with previous comorbidities
mod_prs_comorb <- coxph(Surv(endfollowup - startfollowup,mi) ~ age + scale(prs) + sex + sbp + smoke + bmi + diabetes + lipid_lowering + tot_chol + hdl_chol + ldl_chol, data = total.newICD)
cindex_prs_comorb <- survConcordance(Surv(endfollowup - startfollowup, mi) ~ predict(mod_prs_comorb), data = total.newICD)$concordance
err_prs_comorb <- as.numeric(survConcordance(Surv(endfollowup - startfollowup, mi) ~ predict(mod_prs_comorb), data = total.newICD)$std.err)
mod_prs_comorb_base <- coxph(Surv(endfollowup - startfollowup,mi) ~ age + sex + sbp + smoke + bmi + diabetes + lipid_lowering + tot_chol + hdl_chol + ldl_chol, data = total.newICD)
cindex_prs_comorb_base <- survConcordance(Surv(endfollowup - startfollowup, mi) ~ predict(mod_prs_comorb_base), data = total.newICD)$concordance
err_prs_comorb_base <- as.numeric(survConcordance(Surv(endfollowup - startfollowup, mi) ~ predict(mod_prs_comorb_base), data = total.newICD)$std.err)
comorb <- cbind(group = "All", cindex = cindex_prs_comorb, error = err_prs_comorb, cindex.base = cindex_prs_comorb_base, error.base = err_prs_comorb_base)
# FIGURE 3: plot predictive ability of PRS in population with previous comorbidities and population with no previous comorbidities (c-index)
cindex.comorb.nocomorb <- rbind(no_comorb, comorb)
cindex.comorb.nocomorb <- data.frame(cindex.comorb.nocomorb)
cindex.comorb.nocomorb$cindex <- as.numeric(as.character(cindex.comorb.nocomorb$cindex))
cindex.comorb.nocomorb$error <- as.numeric(as.character(cindex.comorb.nocomorb$error))
cindex.comorb.nocomorb$cindex.base <- as.numeric(as.character(cindex.comorb.nocomorb$cindex.base))
cindex.comorb.nocomorb$error.base <- as.numeric(as.character(cindex.comorb.nocomorb$error.base))
cindex.comorb.nocomorb$diff <- as.numeric(cindex.comorb.nocomorb$cindex - cindex.comorb.nocomorb$cindex.base)
cindex.comorb.nocomorb$diff.error <- cindex.comorb.nocomorb$error - cindex.comorb.nocomorb$error.base
# FIGURE 3A: predictive ability of PRS in population with and without previous relevant comorbidities (c-index)
ggplot(data = cindex.comorb.nocomorb, mapping = aes(group, diff)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = diff - diff.error, ymax = diff + diff.error), width = 0.2, position = position_dodge(1)) + 
  xlab("Group") + ylab("Difference in C-index Compared to Baseline") + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))




















##################################################
# Analysis of dichotomous PRS for population with and without previous revelant comorbidities
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
  mod_dich_prs <- coxph(Surv(endfollowup - startfollowup, mi) ~ age + total.newPRS[, cutoffs[i]] + sex + sbp + smoke + bmi + diabetes + lipid_lowering + tot_chol + hdl_chol + ldl_chol, data = total.newPRS)
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
  mod_dich_prs <- coxph(Surv(endfollowup - startfollowup, mi) ~ age + total.newPRS[, cutoffs[i]] + sex + sbp + smoke + bmi + diabetes + lipid_lowering + tot_chol + hdl_chol + ldl_chol, data = total.newPRS)
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

# FIGURE 3B: association of dichotomous PRS and MI in population with and without previous relevant comorbidites (hazard ratio)
ggplot(data = cindex.dich.prs, mapping = aes(group, hr)) +
  geom_point() +
  geom_errorbar(aes(ymin = hrmin, ymax = hrmax), width = 0.2, position = position_dodge(1)) + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x = "Group", y = "Hazard Ratio")




















##################################################
# ROC analysis of PRS with ICD10 codes and dichotomous PRS
##################################################

# create dataset for ROC analysis
total.newSS <- total.new[which(total.new$icd10 %in% resdf.new$ICD10), c("eid", "icd10", "mi", "prs")]
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

# find AUC of dichotomous PRS
auc.prs <- as.numeric(roc(total.newSS$mi, total.newSS$prs)$auc)
auc.prs.top1 <- as.numeric(roc(total.newSS$mi, total.newSS$top1)$auc)
auc.prs.top2 <- as.numeric(roc(total.newSS$mi, total.newSS$top2)$auc)
auc.prs.top5 <- as.numeric(roc(total.newSS$mi, total.newSS$top5)$auc)
auc.prs.top10 <- as.numeric(roc(total.newSS$mi, total.newSS$top10)$auc)
cindex_dich <- data.frame(cutoffs, cindex = c(auc.prs.top1, auc.prs.top2, auc.prs.top5, auc.prs.top10))
pop.table$cindex <- ifelse(pop.table$label %in% icd.new, resdf$cindex_icd_base, "")
pop.table[c("top1", "top2", "top5", "top10"), "cindex"] <- cindex_dich$cindex

# make labels for figure
pop.table$sens1 <- signif(as.numeric(pop.table$sens), 3)
pop.table$spec1 <- signif(as.numeric(pop.table$spec), 3)
pop.table$cindex1 <- signif(as.numeric(pop.table$cindex), 3)
pop.table$glabel <- ""
for (i in 1:nrow(pop.table)) {
  icd10.label <- pop.table[i, "label"]
  sens.label <- paste("Sens: ", pop.table[i, "sens1"])
  spec.label <- paste("Spec: ", pop.table[i, "spec1"])
  glabel <- paste(icd10.label, sens.label, spec.label, sep="\n")
  pop.table[i, "glabel"] <- glabel
  print(i)
}
# CHANGE: replace list with important ICD10 codes (e.g., ICD10 codes with large AUCs)
pop.table$glabel <- ifelse(pop.table$label %in% c("I20", "I24", "R07", "top1", "top2", "top5", "top10"), pop.table[, "glabel"], "")

# find sensitivities and specificities for PRS 
roc.table <- cbind(roc(total.newSS$mi, total.newSS$prs)$sensitivities, roc(total.newSS$mi, total.newSS$prs)$specificities)
roc.table <- data.frame(roc.table)
colnames(roc.table) <- c("sens", "spec")
roc.table$fpr <- 1-roc.table$spec

# FIGURE 1C: plot ROC of PRS with ICD10 codes and dichotomous PRS
ggplot() + 
  geom_point(data = roc.table, mapping = aes(fpr, sens), size = 1) +
  geom_point(data = pop.table, mapping = aes(fpr, sens, color = type), size = 3) +
  geom_label_repel(data = pop.table, aes(fpr, sens, label = glabel, fill = type), size = 3, segment.alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0) +
  annotate(geom = "text", x = 0, y = 0, label = paste("C-index: ", substr(as.character(as.numeric(roc(total.newSS$mi, total.newSS$prs)$auc)), 0, 5))) +
  coord_cartesian(xlim = c(0, 0.25), ylim = c(0, 0.25)) + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  labs(x = "False Positive Rate", y = "True Positive Rate")