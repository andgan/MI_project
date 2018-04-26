# load libraries 
# install.packages("survival")
library(survival)
library(ggplot2)
library(ggrepel)

# identify the data path
path <- "C:/Jiwoo Lee/Myocardial Infarction Research Project 2017/"
path <- "/Users/andreaganna/Documents/Work/Post_doc/jiwoo/"

source(paste0(path,"MI_project/utils.R"))

########################################################################################################
#  Filter data - Goal: Find association between ICD codes and MI - check how this relates to the PRS
########################################################################################################

# Load data
total <- load_data(path)

# make mi_date as a date rather than a string
total$mi_date <- as.Date(total$mi_date, "%Y-%m-%d")

# make icd10_date as a date rather than a string
total$icd10_date <- as.Date(total$icd10_date, "%Y-%m-%d")

# create end-of-follow up analysis, which is first occurence of 2015-03-31 or MI or death
total$endfollowup <- ifelse(!is.na(total$death_date) & is.na(total$mi_date), as.Date(total$death_date, "%Y-%m-%d"),
					 ifelse(!is.na(total$mi_date), as.Date(total$mi_date, "%Y-%m-%d"), as.Date("2015-03-01", "%Y-%m-%d")))
class(total$endfollowup)<- "Date"

# create start-of-follow up, which is 1998-01-01
total$startfollowup <- as.Date("1998-01-01", "%Y-%m-%d")
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
#2816591 rows
length(unique(total.new$eid))
#495453 individuals
length(unique(total.new$icd10))
#1393 ICD codes
length(unique(total.new[which(!is.na(total.new$mi_date)), ]$eid))
#14584 individuals with MI

# set to missing ICD10 codes/dates outside of timeframe or less than 7 days from MI date
total.new$icd10[total.new$icd10_date > as.Date("2015-03-01", "%Y-%m-%d") | total.new$icd10_date < total.new$startfollowup] <- NA
total.new$icd10_date[total.new$icd10_date > as.Date("2015-03-01", "%Y-%m-%d") | total.new$icd10_date < total.new$startfollowup] <- NA
total.new$icd10[(total.new$mi_date - total.new$icd10_date) < 8] <- NA
total.new$icd10_date[(total.new$mi_date - total.new$icd10_date) < 8] <- NA
nrow(total.new)
#2816591 rows
length(unique(total.new$eid))
#495453 individuals
length(unique(total.new$icd10))
#1385 ICD codes
length(unique(total.new[which(!is.na(total.new$mi_date)), ]$eid))
#14584 individuals with MI

# remove individuals with I25 before MI
# LET'S CHECK WITH ERIK IF THIS EXCLUSION IS NEED IT
to_remove <- total.new$eid[total.new$icd10 == "I25"]
to_remove <- unique(to_remove[!is.na(to_remove)])
total.new = total.new[!total.new$eid %in% to_remove, ]
nrow(total.new)
#2644684 rows
length(unique(total.new$eid))
#483171 individuals
length(unique(total.new$icd10))
#1375 ICD codes
length(unique(total.new[which(!is.na(total.new$mi_date)), ]$eid))
#12765 individuals with MI

# keep only unique individuals with MI
total.newT <- total.new[!is.na(total.new$mi_date), ]
total.newMI <- total.newT[!duplicated(total.newT$eid), c("eid","mi_date")]

# keep only unique individuals (and a subset of variables from total.new)
total.newU <- total.new[!duplicated(total.new$eid), c("eid", "sex", "age", "death", "death_date", "endfollowup", "startfollowup","PRS_0.5")]


##########################################################################
# association analysis: Association between MI and ICD codes before MI
##########################################################################

# Check association between the PRS and MI 
total.newTEMP <- merge(total.newU,total.newMI,by="eid",all.x=T)
total.newTEMP$mi <- ifelse(is.na(total.newTEMP$mi_date),0,1)
modprs <- coxph(Surv(endfollowup-startfollowup,mi) ~ age + scale(PRS_0.5) + sex, data = total.newTEMP)

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
      
      res <- rbind(res, c(coef(summary(mod))[2, c(1, 4)],coef(summary(mod1))[5, c(1, 4)],coef(summary(mod2))[2, c(1, 4)]))
      
      icd.final <- c(icd.final, icd.all[i])
      
      print(i)
}
rownames(res) <- icd.final
resdf = as.data.frame(res)
resdf = cbind(icd10 = rownames(resdf), resdf)
rownames(resdf) <- 1:nrow(resdf)
colnames(resdf) <- c("ICD10","beta_main","z_main","beta_interaction","z_interaction","beta_just_mi","z_just_mi")
write.table(resdf, file = paste0(path,'survivalanalysis_pre.tsv'), sep = "\t", row.names = FALSE, col.names = TRUE, quote=F)



## PLOT INTERACTION BETWEEN PRS and ICD for association with MI ##
resdf <- read.table(file = paste0(path,'survivalanalysis_pre.tsv'), sep = "\t", header = TRUE, stringsAsFactors = F)
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





