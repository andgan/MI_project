
# load libraries 
# install.packages("survival")
library(survival)
library(ggplot2)
library(ggrepel)

# identify the data path
path <- "C:/Jiwoo Lee/Myocardial Infarction Research Project 2017/"
path <- "/Users/andreaganna/Documents/Work/Post_doc/jiwoo/"


##################################################
# Process registry data and phenotype data 
##################################################

# ID, ICD10, and ICD10 date from registry data
total = read.table(file = paste0(path,'hesin_registry.tsv'), sep = '\t', header = TRUE, stringsAsFactors = FALSE, na.strings = "")
total = total[, c(1, 11, 23)]
icd10 = substr(total$diag_icd10, 0, 3)
total = cbind(total, icd10)
total = total[, c(1, 4, 3)]
write.table(total, file = paste0(path,'hesin_registry_new.tsv'), sep = '\t', row.names=F, quote=F)
rm(icd10)

# plot distribution of ICD10 codes to find timeframe for study
plot(table(total$epistart), xlab = "ICD Date", ylab = "Frequency")
# January 1998 to April 2015

# get and clean assessment center data
load('out4.Rdata')
# sex (f31), assessment center visit date (f53), assessment center visit age (f21003), assessment center location (f54), date of myocardial infarction (f42000), date of stroke (f42006), systolic blood pressure (f4080), diastolic blood pressure (f4079), body mass index (f21001), and smoking status (f20116)
bdE4_new = bdE4[, c("f.eid", "f.31.0.0", "f.53.0.0", "f.21003.0.0", "f.54.0.0", "f.42000.0.0", "f.42006.0.0", "f.4080.0.0", "f.4079.0.0", "f.21001.0.0", "f.20116.0.0")]
write.table(bdE4_new, file='assess_cent_all.tsv',sep='\t', row.names=F, quote=F)
rm(bdE4, bdE4_new)

# merge registry data and assessment center data
reg = read.table(file = paste0(path,'hesin_registry_new.tsv'), sep = '\t', header = TRUE, stringsAsFactors = FALSE)
ac = read.table(file = paste0(path,'assess_cent_all.tsv'), sep = '\t', header = TRUE, stringsAsFactors = FALSE)
load(paste0(path,'dU.Rdata'))
ac.new = ac[, c("f.eid", "f.21003.0.0", "f.42000.0.0", "f.42006.0.0", "f.4080.0.0", "f.4079.0.0")]
ac.all = merge(ac.new, dU, by.x = "f.eid", by.y = "eid", all.x = TRUE)
data = merge(reg, ac.all, by.x = "eid", by.y = "f.eid", all.y = TRUE)
colnames(data) = c("eid", "icd10", "icd10_date", "ac_age", "mi_date", "stroke_date", "sbp", "dbp", "birth_date", "sex", "death", "death_date", "ac_date", "ac_location", "bmi", "bp", "smoke", "age")
data = data[, c("eid", "sex", "age", "birth_date", "death", "death_date", "bp", "sbp", "dbp", "bmi", "smoke", "mi_date", "stroke_date", "icd10", "icd10_date", "ac_date", "ac_location", "ac_age")]
write.table(data, file = paste0(path,'hesin_registry_assess_cent_all_v2.csv'), sep = '\t', row.names=F, quote=F)
rm(reg, ac, ac.new, ac.all, data)
