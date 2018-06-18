# load libraries 
# install.packages("survival")
library(survival)
library(ggplot2)
library(ggrepel)
library(reshape2)

# identify the data path
path <- "C:/Jiwoo Lee/Myocardial Infarction Research Project 2017/"
path <- "/Users/andreaganna/Documents/Work/Post_doc/jiwoo/"

source(paste0(path,"MI_project/utils.R"))

##########################################################################################################
# Filter data - Goal: Find association between ICD codes and MI - check how this relates to the PRS.
# Sensitivity analysis - keep only ICD codes after baseline
###########################################################################################################

# Load data from utils.R
total <- load_data(path)

# make mi_date as a date rather than a string
total$mi_date <- as.Date(total$mi_date, "%Y-%m-%d")

# make icd10_date as a date rather than a string
total$icd10_date <- as.Date(total$icd10_date, "%Y-%m-%d")

# create end-of-follow up analysis, which is first occurence of 2015-03-31 or MI or death
total$endfollowup <- ifelse(!is.na(total$death_date) & is.na(total$mi_date), as.Date(total$death_date, "%Y-%m-%d"),
                            ifelse(!is.na(total$mi_date), as.Date(total$mi_date, "%Y-%m-%d"), as.Date("2015-03-01", "%Y-%m-%d")))
class(total$endfollowup)<- "Date"

# create start-of-follow up, which is entrance in the study (THIS IS FOR SENSITIVITY ANALYSIS)
total$startfollowup <- as.Date(total$ac_date, "%Y-%m-%d")
nrow(total)
#2864652 rows
length(unique(total$eid))
#498828 individuals
length(unique(total$icd10))
#1472 ICD codes
length(unique(total[which(!is.na(total$mi_date)), ]$eid))
#17959 individuals with MI

# set to missing ICD10 codes/dates that are not related to disease (e.g., car crash) 
total.new = total
total.new$icd10[grepl("^V|^Z", total.new$icd10)] <- NA
total.new$icd10_date[grepl("^V|^Z", total.new$icd10)] <- NA
nrow(total.new)
#2864652 rows
length(unique(total.new$eid))
#498828 individuals
length(unique(total.new$icd10))
#1393 ICD codes
length(unique(total.new[which(!is.na(total.new$mi_date)), ]$eid))
#17959 individuals with MI

# remove individuals with myocardial infarction outside of timeframe
to_remove <- total.new$eid[total.new$mi_date > as.Date("2015-03-01", "%Y-%m-%d") | total.new$mi_date < total.new$startfollowup]
to_remove <- unique(to_remove[!is.na(to_remove)])
total.new = total.new[!total.new$eid %in% to_remove,]
nrow(total.new)
#2681645 rows
length(unique(total.new$eid))
#486626 individuals
length(unique(total.new$icd10))
#1389 ICD codes
length(unique(total.new[which(!is.na(total.new$mi_date)), ]$eid))
#5757 individuals with MI

# set to missing ICD10 codes/dates outside of timeframe or less than 7 days from MI date
total.new$icd10[total.new$icd10_date > as.Date("2015-03-01", "%Y-%m-%d") | total.new$icd10_date < total.new$startfollowup] <- NA
total.new$icd10_date[total.new$icd10_date > as.Date("2015-03-01", "%Y-%m-%d") | total.new$icd10_date < total.new$startfollowup] <- NA
# This is different between pre and post analysis
total.new$icd10[(total.new$mi_date - total.new$icd10_date) < 8] <- NA
total.new$icd10_date[(total.new$mi_date - total.new$icd10_date) < 8] <- NA
nrow(total.new)
#2681645 rows
length(unique(total.new$eid))
#486626 individuals
length(unique(total.new$icd10))
#1274 ICD codes
length(unique(total.new[which(!is.na(total.new$mi_date)), ]$eid))
#5757 individuals with MI

# remove individuals with I25 before MI
# LET'S CHECK WITH ERIK IF THIS EXCLUSION IS NEED IT
to_remove <- total.new$eid[total.new$icd10 == "I25"]
to_remove <- unique(to_remove[!is.na(to_remove)])
total.new = total.new[!total.new$eid %in% to_remove, ]
nrow(total.new)
#2590659 rows
length(unique(total.new$eid))
#479990 individuals
length(unique(total.new$icd10))
#1268 ICD codes
length(unique(total.new[which(!is.na(total.new$mi_date)), ]$eid))
#5245 individuals with MI

# keep only unique individuals with MI
total.newT <- total.new[!is.na(total.new$mi_date), ]
total.newMI <- total.newT[!duplicated(total.newT$eid), c("eid","mi_date")]

# keep only unique individuals (and a subset of variables from total.new)
total.newU <- total.new[!duplicated(total.new$eid), c("eid", "sex", "age", "death", "death_date", "endfollowup", "startfollowup","PRS_0.5","sbp","smoke","bmi")]


###############################################################################################################
# association analysis: Association between MI and ICD codes before MI. Including only events after baseline.
################################################################################################################

# Check association between the PRS and MI 
total.newTEMP <- merge(total.newU,total.newMI,by="eid",all.x=T)
total.newTEMP$mi <- ifelse(is.na(total.newTEMP$mi_date),0,1)
modprs <- coxph(Surv(endfollowup-startfollowup,mi) ~ age + scale(PRS_0.5) + sex, data = total.newTEMP)
base_cindex_prs <- survConcordance(Surv(endfollowup-startfollowup, mi)~predict(modprs), data = total.newTEMP[!is.na(total.newTEMP$PRS_0.5),])$concordance


# Calculate C-index for baseline model without inclusion of the ICD code
mod_base <- coxph(Surv(endfollowup-startfollowup,mi) ~ age + sex, data = total.newTEMP)
base_cindex_base <- survConcordance(Surv(endfollowup-startfollowup, mi)~predict(mod_base), data = total.newTEMP)$concordance
mod_base <- coxph(Surv(endfollowup-startfollowup,mi) ~ age  + sex + sbp + smoke + bmi, data = total.newTEMP)
risk_factor_base <- survConcordance(Surv(endfollowup-startfollowup, mi)~predict(mod_base), data = total.newTEMP[!is.na(total.newTEMP$sbp) & !is.na(total.newTEMP$smoke) & !is.na(total.newTEMP$bmi),])$concordance


# loop across ICD10 codes associated with MI and n > 10
icd.mi <- total.new[!is.na(total.new$mi_date),c("eid","icd10")]
icd.mi <- unique(icd.mi[,c('eid','icd10')])
icd.mi <- icd.mi[!is.na(icd.mi$icd10),"icd10"]
icd.all <- names(table(icd.mi))[table(icd.mi) > 10]
res <- NULL
icd.final <- NULL
for (i in 1:length(icd.all)) {
  
  # Keep unique record for the comorbity examined
  total.newT <- total.new[total.new$icd10 == icd.all[i] & !is.na(total.new$icd10), ]
  
  # Sort and keep the first event for each individual
  total.new.sorted <- total.newT[order(total.newT$eid, as.numeric(total.newT$icd10_date)), ]
  total.new.sortedU <- total.new.sorted[!duplicated(total.new.sorted$eid), c("eid", "icd10_date")]
  
  # Merge the ICD-only dataset with the main dataset
  temp1 <- merge(total.newU, total.new.sortedU, all.x = T, by = "eid")
  
  # Now add the MI-only dataset
  total.comorb <- merge(temp1, total.newMI, all.x = T, by = "eid")
  
  # Assign pred and MI
  total.comorb$pred <- ifelse(is.na(total.comorb$icd10_date), 0, 1)
  total.comorb$mi <- ifelse(is.na(total.comorb$mi_date), 0, 1)
  
  
  # Calculate mean time distance between ICD code and MI
  total.comorbT <- total.comorb[total.comorb$pred==1 & total.comorb$mi==1,]
  mean_day_dist <- mean(total.comorbT$mi_date-total.comorbT$icd10_date)
  sd_day_dist <- sd(total.comorbT$mi_date-total.comorbT$icd10_date)
  
  # If you have the ICD you contribute as unexposed (trt=0) from when you enter to when you get the ICD
  data1 <- subset(total.comorb, pred == 1)
  data1$tstart <- 0
  data1$tstop <- as.numeric(data1$icd10_date - data1$startfollowup)
  data1$outP <- 0
  data1$trt <- 0
  
  # If you have the ICD, you contribute as exposed (trt=1) from when you get the ICD untill end of follow-up (which can be death/MI/of 2015)
  data2 <- subset(total.comorb, pred == 1)
  data2$tstart <- as.numeric(data2$icd10_date - data2$startfollowup)
  data2$tstop <- as.numeric(data2$endfollowup - data2$startfollowup)
  data2$outP <- data2$mi
  data2$trt <- 1
  
  # If you don't have the ICD, but you have MI you contribute as unexposed (trt=0) from when you enter until the MI date
  data3 <- subset(total.comorb, pred == 0 & mi == 1)
  data3$tstart <- 0
  data3$tstop <- as.numeric(data3$endfollowup - data3$startfollowup)
  data3$outP <- 1
  data3$trt <- 0
  
  dataF <- rbind(data1, data2, data3)
  
  # If you don't have the ICD or MI you contribute as unexposed (trt=0) from when you enter untill you exit the study
  otherdata <- subset(total.comorb, pred != 1 & mi != 1)
  
  otherdata$tstart <- 0
  otherdata$tstop <- as.numeric(otherdata$endfollowup - otherdata$startfollowup)
  otherdata$outP <- 0
  otherdata$trt <- 0
  
  findata <- rbind(dataF, otherdata)
  
  # This is because there is one individuals where tstart=tstop
  findata <- findata[findata$tstart < findata$tstop, ]
  mod <- tryCatch({
    coxph(Surv(tstart, tstop, outP) ~ age + trt + sex, data = findata)
  }, error = function(w) {
    NA
  })
  
  # C-index / discrimination measure for base model: age, sex + icd code
  base_cindex <- survConcordance(Surv(tstart, tstop, outP)~predict(mod), data = findata)$concordance
  
  
  # Just a test to calculate the base C-index for the model with PRS
  modtt <- tryCatch({
    coxph(Surv(tstart, tstop, outP) ~ age + scale(PRS_0.5) + sex, data = findata[!is.na(findata$PRS_0.5),])
  }, error = function(w) {
    NA
  })
  
  base_cindex_prs_test <- survConcordance(Surv(tstart, tstop, outP)~predict(modtt), data = findata[!is.na(findata$PRS_0.5),])$concordance
  
  
  
  # Adjusting for risk factors
  mod3 <- tryCatch({
    coxph(Surv(tstart, tstop, outP) ~ age + trt + sex + sbp + smoke + bmi, data = findata)
  }, error = function(w) {
    NA
  })
  
  # C-index / discrimination measure for base model + risk factors
  risk_factor_cindex <- survConcordance(Surv(tstart, tstop, outP)~predict(mod3), data = findata[!is.na(findata$sbp) & !is.na(findata$smoke) & !is.na(findata$bmi),])$concordance
  
  
  # Testing the interaction between ICD code and PRS
  mod1 <- tryCatch({
    coxph(Surv(tstart, tstop, outP) ~ age + trt*PRS_0.5  + sex, data = findata)
  }, error = function(w) {
    NA
  })
  
  # Testing for association between previous ICD code and PRS in individuals with MI
  datax <- subset(total.comorb, mi == 1)
  datax$endfollowup[datax$pred==1] <- datax$icd10_date[datax$pred==1]
  datax$tstart <- 0
  datax$tstop <- as.numeric(datax$endfollowup - datax$startfollowup)
  datax$outP <- datax$pred
  
  # This is because there is one individuals where tstart=tstop
  datax <- datax[datax$tstart < datax$tstop, ]
  mod2 <- tryCatch({
    coxph(Surv(tstart, tstop, outP) ~ age + PRS_0.5 + sex, data = datax)
  }, error = function(w) {
    NA
  })
  
  
  
  res <- rbind(res, c(coef(summary(mod))[2, c(1, 4)],coef(summary(mod1))[5, c(1, 4)],coef(summary(mod2))[2, c(1, 4)],coef(summary(mod3))[2, c(1, 4)],mean_day_dist,sd_day_dist,base_cindex,risk_factor_cindex,base_cindex_prs_test))
  
  icd.final <- c(icd.final, icd.all[i])
  
  print(i)
}
rownames(res) <- icd.final
resdf = as.data.frame(res)
resdf = cbind(icd10 = rownames(resdf), resdf)
rownames(resdf) <- 1:nrow(resdf)
colnames(resdf) <- c("ICD10","beta_main","z_main","beta_interaction","z_interaction","beta_just_mi","z_just_mi","beta_main_adj","z_main_adj","mean_day_dist","sd_day_dist","cindex_base","cindex_risk_factor","base_cindex_prs_test")
write.table(resdf, file = paste0(path,'survivalanalysis_pre_sensitivity.tsv'), sep = "\t", row.names = FALSE, col.names = TRUE, quote=F)


#### DISTROBUTION TIME DISTANCE BETWEEN ICD CODE AND MI ####
resdf <- read.table(file = paste0(path,'survivalanalysis_pre_sensitivity.tsv'), sep = "\t", header = TRUE, stringsAsFactors = F)
resdf$ICD10 <- factor(resdf$ICD10, levels = resdf$ICD10[order(resdf$mean_day_dist)])
ggplot(resdf, aes(ICD10,mean_day_dist)) + 
  geom_bar(stat = "identity") + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1,size=4), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x = "ICD Code", y = "Mean distance to MI (in days)")



### PLOT ASSOCIATION BETWEEN ICD codes and MI ### 
resdf <- read.table(file = paste0(path,'survivalanalysis_pre_sensitivity.tsv'), sep = "\t", header = TRUE, stringsAsFactors = F)
resdf$p_main <- 2*pnorm(-abs(resdf$z_main))
resdf$logp <- -log10(resdf$p_main)
resdf$logp[resdf$logp=="Inf"] <- 300
resdf$ICD10label <- resdf$ICD10
resdf$ICD10label[resdf$logp < -log10(0.05 / length(icd.all))] <- ""
ggplot(resdf, aes(ICD10,logp)) + 
  geom_point(aes(size=exp(beta_main))) + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_hline(yintercept = -log10(0.05 / length(icd.all)), size = 1.5, color = "red") + 
  labs(x = "ICD Code", y = "-log10 Transformed P-Value") + 
  geom_text_repel(aes(label = resdf$ICD10label),size=3,segment.alpha=0.5) + expand_limits(y = 200) +
  scale_size_continuous("HR",breaks=c(2,5,10,15,20,25,30), labels=c("2","5","10","15","20","25","30"))



### PLOT ASSOCIATION BETWEEN ICD codes and MI - PLOT C-INDEX ### 
resdf <- read.table(file = paste0(path,'survivalanalysis_pre_sensitivity.tsv'), sep = "\t", header = TRUE, stringsAsFactors = F)
resdf$p_main <- 2*pnorm(-abs(resdf$z_main))
resdf$logp <- -log10(resdf$p_main)
resdf$logp[resdf$logp=="Inf"] <- 300
resdf$ICD10label <- resdf$ICD10
resdf$ICD10label[resdf$logp < -log10(0.05 / length(icd.all))] <- ""
ggplot(resdf, aes(ICD10,cindex_base)) + 
  geom_point(aes(size=exp(beta_main))) + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_hline(yintercept = mean(resdf$base_cindex_prs_test), size = 1.5, color = "red", linetype=2) +
  labs(x = "ICD Code", y = "C-index (including age and sex)") + 
  geom_text_repel(aes(label = resdf$ICD10label),size=3,segment.alpha=0.5)  +
  scale_size_continuous("HR",breaks=c(2,5,10,15,20,25,30), labels=c("2","5","10","15","20","25","30"))


### PLOT ASSOCIATION BETWEEN ICD codes and MI, ADJUSTED AND NON-ADJUSTED ### 
resdf <- read.table(file = paste0(path,'survivalanalysis_pre_sensitivity.tsv'), sep = "\t", header = TRUE, stringsAsFactors = F)
resdf$p_main <- 2*pnorm(-abs(resdf$z_main))
resdf$logp_main <- -log10(resdf$p_main)
resdf$logp_main[resdf$logp_main=="Inf"] <- 300
resdf$p_adj <- 2*pnorm(-abs(resdf$z_main_adj))
resdf$logp_adj <- -log10(resdf$p_adj)
resdf$logp_adj[resdf$logp_adj=="Inf"] <- 300
resdf$ICD10label <- resdf$ICD10

resdf2 <- melt(resdf,measure.vars = c("logp_adj","logp_main"))
resdf2$ICD10label[resdf2$value < -log10(0.05 / length(icd.all)) | resdf2$variable=="logp_main"] <- ""
ggplot(resdf2, aes(ICD10,value)) + 
  geom_point(aes(size=exp(beta_main),col=variable)) + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_hline(yintercept = -log10(0.05 / length(icd.all)), size = 1.5, color = "red") + 
  labs(x = "ICD Code", y = "-log10 Transformed P-Value") + 
  geom_text_repel(aes(label = resdf2$ICD10label),size=3,segment.alpha=0.5) + expand_limits(y = 200) +
  scale_size_continuous("HR",breaks=c(2,5,10,15,20,25,30), labels=c("2","5","10","15","20","25","30")) +
  scale_colour_discrete("",breaks=c("logp_adj","logp_main"), labels=c("Adjusted","Unadjusted"))


## PLOT ASSOCIATION BETWEEN PRS for MI AND ICD-CODES ###
resdf <- read.table(file = paste0(path,'survivalanalysis_pre_sensitivity.tsv'), sep = "\t", header = TRUE, stringsAsFactors = F)
resdf$p_icd_code <- 2*pnorm(-abs(resdf$z_icd_code))
resdf$logp <- -log10(resdf$p_icd_code)
resdf$ICD10label <- resdf$ICD10
resdf$ICD10label[resdf$logp < -log10(0.05 / length(icd.all))] <- ""
ggplot(resdf, aes(ICD10,logp)) + 
  geom_point() + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_hline(yintercept = -log10(0.05 / length(icd.all)), size = 1.5, color = "red") + 
  labs(x = "ICD Code", y = "-log10 Transformed P-Value") + 
  geom_text_repel(aes(label = resdf$ICD10label),size=3,segment.alpha=0.5)



## PLOT INTERACTION BETWEEN PRS and ICD for association with MI ##
resdf <- read.table(file = paste0(path,'survivalanalysis_pre_sensitivity.tsv'), sep = "\t", header = TRUE, stringsAsFactors = F)
icd.sel <- resdf$ICD10[2*pnorm(-abs(resdf$z_main)) < (0.05/nrow(resdf))]
resdfS <- resdf[resdf$ICD10 %in% icd.sel, ]
resdfS$p_interaction <- 2*pnorm(-abs(resdfS$z_interaction))
resdfS$logp <- -log10(resdfS$p_interaction)
resdfS$ICD10label <- resdfS$ICD10
resdfS$ICD10label[resdfS$logp < -log10(0.05 / length(icd.sel))] <- ""
ggplot(resdfS, aes(ICD10,logp)) + 
  geom_point() + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_hline(yintercept = -log10(0.05 / length(icd.all)), size = 1.5, color = "red") + 
  labs(x = "ICD Code", y = "-log10 Transformed P-Value") + 
  geom_text_repel(aes(label = resdfS$ICD10label),size=3,segment.alpha=0.5)



# top five ICD10 codes are I20, R07, I50, I24, and J22
resdf[which(resdf$logp>quantile(resdf$logp,0.95)),]


## PLOT AVERAGE PRS FOR MI, NO MI, ICD, NO ICD
plot_icd_mi_prs("R07",total.new,total.newMI)
plot_icd_mi_prs("I20",total.new,total.newMI)
plot_icd_mi_prs("I24",total.new,total.newMI)
plot_icd_mi_prs("I50",total.new,total.newMI)



#### PLOT ASSOCIATION BETWEEN I20 and MI as function of the PRS ####
plot_icd_mi_prs_int("I20",total.new,total.newMI)





### LOOP ACROSS ICD CODES SIGNIFICANTLY ASSOCIATED WITH MI 
# keep only unique individuals (and a subset of variables from total.new)
total.newU <- total.new[!duplicated(total.new$eid), c("eid", "sex", "age", "death", "death_date", "endfollowup", "startfollowup",colnames(total.new)[grep("^X",colnames(total.new))])]




# Check association between the SNPs and MI 
total.newTEMP <- merge(total.newU,total.newMI,by="eid",all.x=T)
total.newTEMP$mi <- ifelse(is.na(total.newTEMP$mi_date),0,1)
for (snp in colnames(total.newTEMP)[grep("^X",colnames(total.newTEMP))])
{
  modprs <- coxph(Surv(endfollowup-startfollowup,mi) ~ age + total.newTEMP[,snp] + sex, data = total.newTEMP)
  print(modprs)
}


icd.sel <- resdf$ICD10[2*pnorm(-abs(resdf$z_main)) < (0.05/nrow(resdf))]
res <- NULL
icd.final <- NULL
for (i in 1:length(icd.sel)) {
  
  # Keep unique record for the comorbity examined
  total.newT <- total.new[total.new$icd10 == icd.sel[i] & !is.na(total.new$icd10), ]
  
  # Sort and keep the first event for each individual
  total.new.sorted <- total.newT[order(total.newT$eid, as.numeric(total.newT$icd10_date)), ]
  total.new.sortedU <- total.new.sorted[!duplicated(total.new.sorted$eid), c("eid", "icd10_date")]
  
  # Merge the ICD-only dataset with the main dataset
  temp1 <- merge(total.newU, total.new.sortedU, all.x = T, by = "eid")
  
  # Now add the MI-only dataset
  total.comorb <- merge(temp1, total.newMI, all.x = T, by = "eid")
  
  # Assign pred and MI
  total.comorb$pred <- ifelse(is.na(total.comorb$icd10_date), 0, 1)
  total.comorb$mi <- ifelse(is.na(total.comorb$mi_date), 0, 1)
  
  # If you have the ICD you contribute as unexposed (trt=0) from when you enter to when you get the ICD
  data1 <- subset(total.comorb, pred == 1)
  data1$tstart <- 0
  data1$tstop <- as.numeric(data1$icd10_date - data1$startfollowup)
  data1$outP <- 0
  data1$trt <- 0
  
  # If you have the ICD, you contribute as exposed (trt=1) from when you get the ICD untill end of follow-up (which can be death/MI/of 2015)
  data2 <- subset(total.comorb, pred == 1)
  data2$tstart <- as.numeric(data2$icd10_date - data2$startfollowup)
  data2$tstop <- as.numeric(data2$endfollowup - data2$startfollowup)
  data2$outP <- data2$mi
  data2$trt <- 1
  
  # If you don't have the ICD, but you have MI you contribute as unexposed (trt=0) from when you enter until the MI date
  data3 <- subset(total.comorb, pred == 0 & mi == 1)
  data3$tstart <- 0
  data3$tstop <- as.numeric(data3$endfollowup - data3$startfollowup)
  data3$outP <- 1
  data3$trt <- 0
  
  dataF <- rbind(data1, data2, data3)
  
  # If you don't have the ICD or MI you contribute as unexposed (trt=0) from when you enter untill you exit the study
  otherdata <- subset(total.comorb, pred != 1 & mi != 1)
  
  otherdata$tstart <- 0
  otherdata$tstop <- as.numeric(otherdata$endfollowup - otherdata$startfollowup)
  otherdata$outP <- 0
  otherdata$trt <- 0
  
  findata <- rbind(dataF, otherdata)
  # This is because there is one individuals where tstart=tstop
  findata <- findata[findata$tstart < findata$tstop, ]
  
  for (snp in colnames(total.new)[grep("^X",colnames(total.new))])
  {
    mod  <- tryCatch({
      coxph(Surv(tstart, tstop, outP) ~ age + trt*findata[,snp]  + sex, data = findata)
    }, error = function(w) {
      NA
    })
    
    res <- rbind(res, c(as.character(icd.sel[i]),snp,coef(summary(mod))[5, c(1, 4)]))
  }
 
  icd.final <- c(icd.final, icd.sel[i])
  
  print(i)
}

rownames(res) <- icd.final
resdf_snp = as.data.frame(res)
resdf_snp = cbind(icd10 = rownames(resdf_snp), resdf_snp)
rownames(resdf_snp) <- 1:nrow(resdf_snp)
colnames(resdf_snp) <- c("ICD10","snp","beta_snp","z_snp")
write.table(resdf_snp, file = paste0(path,'survivalanalysis_pre_sensitivity_snp.tsv'), sep = "\t", row.names = FALSE, col.names = TRUE, quote=F)



