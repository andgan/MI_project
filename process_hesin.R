##################################################
# Process registry data and phenotype data 
##################################################

# load libraries 
# install.packages("survival")
# install.packages("ggplot2")
# install.packages("ggrepel")
library(survival)
library(ggplot2)
library(ggrepel)

# identify the data path
path <- "C:/Jiwoo Lee/Myocardial Infarction Research Project 2017/"
path <- "/Users/andreaganna/Documents/Work/Post_doc/jiwoo/"

# get ID, ICD10, and ICD10 date from registry data
total = read.table(file = paste0(path,'hesin_registry.tsv'), sep = '\t', header = TRUE, stringsAsFactors = FALSE, na.strings = "")
total = total[, c(1, 11, 23)]
icd10 = substr(total$diag_icd10, 0, 3)
total = cbind(total, icd10)
total = total[, c(1, 4, 3)]
write.table(total, file = paste0(path,'hesin_registry_new.tsv'), sep = '\t', row.names = F, quote = F)
rm(icd10)

# plot distribution of ICD10 codes to find timeframe for study
plot(table(total$epistart), xlab = "ICD Date", ylab = "Frequency")
# January 1998 to April 2015

# get and clean assessment center data
load(paste0(path,'out4.Rdata'))
# sex (f31), assessment center visit date (f53), assessment center visit age (f21003), assessment center location (f54), date of myocardial infarction (f42000), date of stroke (f42006), systolic blood pressure (f4080), diastolic blood pressure (f4079), body mass index (f21001), and smoking status (f20116), diabetes (f2443)
bdE4_new = bdE4[, c("f.eid", "f.31.0.0", "f.53.0.0", "f.21003.0.0", "f.54.0.0", "f.42000.0.0", "f.42006.0.0", "f.4080.0.0", "f.4079.0.0", "f.21001.0.0", "f.20116.0.0","f.2443.0.0", grep("20002", colnames(bdE4), value = TRUE), grep("20003", colnames(bdE4), value = TRUE), grep("6177", colnames(bdE4), value = TRUE), grep("6153", colnames(bdE4), value = TRUE))]

##################################################
# Diabetes
# Goal: Define type 2 diabetes from self-report information
##################################################

dm_t2dm_sr_ni <- NA
for (i in 0:28) { dm_t2dm_sr_ni[bdE4_new[[paste('f.20002.0.', i, sep="")]] == 1223] <- 1 }                     
dm_t2dm_sr_ni[is.na(dm_t2dm_sr_ni)] <- 0

# define medications for type two diabetes
meds_insulin_sr_ni <- NA
for (i in 0:47) { meds_insulin_sr_ni[bdE4_new[[paste('f.20003.0.', i, sep="")]] == 1140883066] <- 1 }
meds_insulin_sr_ni[is.na(meds_insulin_sr_ni)] <- 0

meds_metformin_sr_ni <- NA
v <- c(1140884600, 1140874686, 1141189090)
for (i in 0:47) { meds_metformin_sr_ni[bdE4_new[[paste('f.20003.0.', i, sep="")]]  %in% v] <- 1 }
meds_metformin_sr_ni[is.na(meds_metformin_sr_ni)] <- 0

meds_nonmet_oad_sr_ni <- NA
v <- c(1140874718, 1140874744, 1140874746, 1141152590, 1141156984, 1140874646, 1141157284, 1140874652, 1140874674, 1140874728, 1140868902, 1140868908, 1140857508, 1141173882, 1141173786, 1141168660, 1141171646, 1141171652, 1141153254, 1141177600, 1141177606)
for (i in 0:47) { meds_nonmet_oad_sr_ni[bdE4_new[[paste('f.20003.0.', i, sep="")]] %in% v] <- 1 }
meds_nonmet_oad_sr_ni[is.na(meds_nonmet_oad_sr_ni)] <- 0

dm_insulin_sr_ts <- NA
for (i in 0:2) { dm_insulin_sr_ts[bdE4_new[[paste('f.6177.0.', i, sep="")]] == "Insulin"] <- 1 }          
for (i in 0:3) { dm_insulin_sr_ts[bdE4_new[[paste('f.6153.0.', i, sep="")]] == "Insulin"] <- 1 }          
dm_insulin_sr_ts[is.na(dm_insulin_sr_ts)] <- 0

meds_any_sr_ni_ts <- NA 
# this single variable captures all 3 categories of DM medications and insulin TS
meds_any_sr_ni_ts[meds_insulin_sr_ni + meds_metformin_sr_ni + meds_nonmet_oad_sr_ni + dm_insulin_sr_ts > 0] <- 1
meds_any_sr_ni_ts[is.na(meds_any_sr_ni_ts)] <- 0

bdE4_new$diabetes <- ifelse(dm_t2dm_sr_ni == 1 | meds_any_sr_ni_ts == 1 | bdE4$f.2443.0.0 == "Yes", 1, 0)

##################################################
# Lipid-lowering medications
##################################################

lipid_medications <- NA
v <- c(1140861958, 1140888648, 1140888594, 1141146234, 1141192410, 1141192736, 140861868, 1140861954)
for (i in 0:47) { lipid_medications[bdE4_new[[paste('f.20003.0.', i, sep="")]] %in% v] <- 1 }
lipid_medications[is.na(lipid_medications)] <- 0

lipid_medications_self <- NA
for (i in 0:2) { lipid_medications_self[bdE4_new[[paste('f.6177.0.', i, sep="")]] == "Cholesterol lowering medication"] <- 1 }          
for (i in 0:3) { lipid_medications_self[bdE4_new[[paste('f.6153.0.', i, sep="")]] == "Cholesterol lowering medication"] <- 1 }          
lipid_medications_self[is.na(lipid_medications_self)] <- 0

bdE4_new$lipid_lowering <- ifelse(lipid_medications == 1 | lipid_medications_self == 1, 1, 0)

write.table(bdE4_new[,!colnames(bdE4_new) %in% c("f.2443.0.0", grep("20002", colnames(bdE4), value = TRUE), grep("20003", colnames(bdE4), value = TRUE), grep("6177", colnames(bdE4), value = TRUE), grep("6153", colnames(bdE4), value = TRUE))], file = 'assess_cent_all.tsv', sep = '\t', row.names = F, quote = F)
rm(bdE4, bdE4_new)

# merge registry data and assessment center data
reg = read.table(file = paste0(path,'hesin_registry_new.tsv'), sep = '\t', header = TRUE, stringsAsFactors = FALSE)
ac = read.table(file = paste0(path,'assess_cent_all.tsv'), sep = '\t', header = TRUE, stringsAsFactors = FALSE)
load(paste0(path,'dU.Rdata'))
ac.new = ac[, c("f.eid", "f.21003.0.0", "f.42000.0.0", "f.42006.0.0", "f.4080.0.0", "f.4079.0.0", "diabetes", "lipid_lowering")]
ac.all = merge(ac.new, dU, by.x = "f.eid", by.y = "eid", all.x = TRUE)
data = merge(reg, ac.all, by.x = "eid", by.y = "f.eid", all.y = TRUE)
colnames(data) = c("eid", "icd10", "icd10_date", "age", "mi_date", "stroke_date", "sbp", "dbp", "diabetes", "lipid_lowering", "birth_date", "sex", "death", "death_date", "ac_date", "ac_location", "bmi", "bp", "smoke", "age_fixed")
data = data[, c("eid", "sex", "birth_date", "death", "death_date", "bp", "sbp", "dbp", "bmi", "smoke", "mi_date", "stroke_date", "icd10", "icd10_date", "ac_date", "ac_location", "age", "diabetes", "lipid_lowering")]
write.table(data, file = paste0(path,'hesin_registry_assess_cent_all_v2.csv'), sep = '\t', row.names = F, quote = F)
rm(reg, ac, ac.new, ac.all, data)