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
#  Filter data - Goal: Find association between MI and ICD codes after MI and the relationship with PRS
########################################################################################################

# Load data
total <- load_data(path)

# make mi_date as a date rather than a string
total$mi_date <- as.Date(total$mi_date, "%Y-%m-%d")

# make icd10_date as a date rather than a string
total$icd10_date <- as.Date(total$icd10_date, "%Y-%m-%d")

# create end-of-follow up analysis, which is first occurence of 2015-03-31 or death
total$endfollowup <- ifelse(!is.na(total$death_date), as.Date(total$death_date, "%Y-%m-%d"), as.Date("2015-03-01", "%Y-%m-%d"))
class(total$endfollowup)<- "Date"

# create start-of-follow up, which is 1998-01-01
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


# set to missing ICD10 codes/dates outside of timeframe or < 30 days after MI
total.new$icd10[total.new$icd10_date > as.Date("2015-03-01", "%Y-%m-%d") | total.new$icd10_date < total.new$startfollowup] <- NA
total.new$icd10_date[total.new$icd10_date > as.Date("2015-03-01", "%Y-%m-%d") | total.new$icd10_date < total.new$startfollowup] <- NA
total.new$icd10[(total.new$mi_date - total.new$icd10_date) > -8] <- NA
total.new$icd10_date[(total.new$mi_date - total.new$icd10_date) > -8] <- NA

nrow(total.new)
#2681645 rows
length(unique(total.new$eid))
#486626 individuals
length(unique(total.new$icd10))
#1276 ICD codes
length(unique(total.new[which(!is.na(total.new$mi_date)), ]$eid))
#5757 individuals with MI

# If individuals have an I25 flag that, this is useful when we want to exclude those from the analysis because such a strong proxy for MI
I25ids <- unique(total.new$eid[total.new$icd10=="I25"])
total.new$I25 <- ifelse(total.new$eid %in% I25ids,1,0) 

# keep only unique individuals with MI
total.newT <- total.new[!is.na(total.new$mi_date), ]
total.newMI <- total.newT[!duplicated(total.newT$eid), c("eid","mi_date")]

# keep only unique individuals (and a subset of variables from total.new)
total.newU <- total.new[!duplicated(total.new$eid), c("eid", "sex", "age", "death", "death_date", "endfollowup", "startfollowup","PRS_0.5","sbp","smoke","bmi","I25")]


########################################################################
# association analysis: Association between MI and ICD codes after MI
#######################################################################

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
    
    # Recalculate end of follow-up, which is the ICD code date (or the existing definition if not ICD code date)
    total.comorb$endfollowup[total.comorb$pred==1] <- total.comorb$icd10_date[total.comorb$pred==1]
    
    # If you have the MI you contribute as unexposed (trt=0) from when you enter to when you get the MI
    data1 <- subset(total.comorb, mi == 1)
    data1$tstart <- 0
    data1$tstop <- as.numeric(data1$mi_date - data1$startfollowup)
    data1$outP <- 0
    data1$trt <- 0
    
    # If you have the MI, you contribute as exposed (trt=1) from when you get the MI until end of follow-up (which can be death/ICD/of 2015)
    data2 <- subset(total.comorb, mi == 1)
    data2$tstart <- as.numeric(data2$mi_date - data2$startfollowup)
    data2$tstop <- as.numeric(data2$endfollowup - data2$startfollowup)
    data2$outP <- data2$pred
    data2$trt <- 1
    
    # If you don't have the MI, but you have ICD you contribute as unexposed (trt=0) from when you enter until the ICD date
    # For some codes (e.g. I22) there is no one that didn't have an MI and so we skip this step. (For these individuals the model does not converge)
    if (nrow(subset(total.comorb, pred == 1 & mi == 0))>0)
    {
      data3 <- subset(total.comorb, pred == 1 & mi == 0)
      data3$tstart <- 0
      data3$tstop <- as.numeric(data3$endfollowup - data3$startfollowup)
      data3$outP <- 1
      data3$trt <- 0
      
      dataF <- rbind(data1, data2, data3)
    }else
    {
      dataF <- rbind(data1, data2)
    }
    
    
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
      coef(summary(coxph(Surv(tstart, tstop, outP) ~ age + trt + sex, data = findata)))[2, c(1, 4)]
    }, error = function(w) {
    c(NA,NA)
    })
    
    # Adjusting for risk factors
    mod3 <- tryCatch({
      coef(summary(coxph(Surv(tstart, tstop, outP) ~ age + trt + sex + sbp + smoke + bmi, data = findata)))[2, c(1, 4)]
    }, error = function(w) {
      c(NA,NA)
    })
    
    
    # Testing the interaction between ICD code and PRS
    mod1 <- tryCatch({
      coef(summary(coxph(Surv(tstart, tstop, outP) ~ age + trt*PRS_0.5  + sex, data = findata)))[5, c(1, 4)]
    }, error = function(w) {
      c(NA,NA)
    })
    
    # Testing for association between ICD code and PRS in individuals with MI
    # This is because there is one individuals where tstart=tstop
    data2x <- data2[data2$tstart < data2$tstop, ]
    mod2 <- tryCatch({
      coef(summary(coxph(Surv(tstart, tstop, outP) ~ age + PRS_0.5 + sex, data = data2x)))[2, c(1, 4)]
    }, error = function(w) {
      c(NA,NA)
    })
    
    
    # Association between PRS and ICD codes in individuals without MI and I25
    total.comorbM <-  total.comorb[total.comorb$mi==0 & total.comorb$I25==0,]
    total.comorbM$tstart <- 0
    total.comorbM$endfollowup  <- ifelse(!is.na(total.comorbM$death_date) & total.comorbM$pred==0, as.Date(total.comorbM$death_date, "%Y-%m-%d"),ifelse(total.comorbM$pred==1, as.Date(total.comorbM$icd10_date, "%Y-%m-%d"), as.Date("2015-03-01", "%Y-%m-%d")))
    total.comorbM$endfollowup <- as.Date(total.comorbM$endfollowup,origin="1970-01-01")
    total.comorbM$tstop <- as.numeric(total.comorbM$endfollowup -  total.comorbM$startfollowup)
    total.comorbM <- total.comorbM[total.comorbM$tstart < total.comorbM$tstop, ]
    mod4 <- tryCatch({
      coef(summary(coxph(Surv(tstart, tstop, pred) ~ age + PRS_0.5 + sex, data = total.comorbM)))[2, c(1, 4)]
    }, error = function(w) {
      c(NA,NA)
    })
    
    
    res <- rbind(res, c(mod,mod1,mod2,mod3,mod4))
    
    
    icd.final <- c(icd.final, icd.all[i])
    print(i)
}

rownames(res) <- icd.final
resdf = as.data.frame(res)
resdf = cbind(icd10 = rownames(resdf), resdf)
rownames(resdf) <- 1:nrow(resdf)
colnames(resdf) <- c("ICD10","beta_main","z_main","beta_interaction","z_interaction","beta_just_mi","z_just_mi","beta_main_adj","z_main_adj","beta_prs_icd","z_prs_icd")
write.table(resdf, file = paste0(path,'survivalanalysis_post_sensitivity.tsv'), sep = "\t", row.names = FALSE, col.names = TRUE, quote=F)

# rm(data1, data2, data3, dataF, findata, i, icd.all, icd.mi, mod, otherdata, res, temp1, to_remove, total.comorb, total.new.sorted, total.new.sortedU, total.newMI, total.newT, total.newU)


# PLOT ASSOCIATION BETWEEN MI AND ICD-CODES FOLLOWING MI
resdf <- read.table(file = paste0(path,'survivalanalysis_post_sensitivity.tsv'), sep = "\t", header = TRUE, stringsAsFactors = F)
resdf$p_main <- 2*pnorm(-abs(resdf$z_main))
resdf$logp <- -log10(resdf$p_main)
resdf$logp[resdf$logp=="Inf"] <- 300
resdf$ICD10label <- resdf$ICD10
resdf$ICD10label[resdf$logp < -log10(0.05 / length(icd.all))] <- ""
resdf <- resdf[!is.na(resdf$beta_main) & resdf$ICD10!="I22",] # This is because the model didn't converged for I22
ggplot(resdf, aes(ICD10,logp)) + 
  geom_point(aes(size=exp(beta_main))) + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_hline(yintercept = -log10(0.05 / length(icd.all)), size = 1.5, color = "red") + 
  labs(x = "ICD Code", y = "-log10 Transformed P-Value") + 
  geom_text_repel(aes(label = resdf$ICD10label),size=3,segment.alpha=0.5) + expand_limits(y = 200) +
  scale_size_continuous("HR",breaks=c(2,5,10,15,20,25,30), labels=c("2","5","10","15","20","25","30"))


### PLOT ASSOCIATION BETWEEN MI codes and ICD, ADJUSTED AND NON-ADJUSTED ### 
resdf <- read.table(file = paste0(path,'survivalanalysis_post_sensitivity.tsv'), sep = "\t", header = TRUE, stringsAsFactors = F)
resdf$p_main <- 2*pnorm(-abs(resdf$z_main))
resdf$logp_main <- -log10(resdf$p_main)
resdf$logp_main[resdf$logp_main=="Inf"] <- 300
resdf$p_adj <- 2*pnorm(-abs(resdf$z_main_adj))
resdf$logp_adj <- -log10(resdf$p_adj)
resdf$logp_adj[resdf$logp_adj=="Inf"] <- 300
resdf$ICD10label <- resdf$ICD10
resdf <- resdf[!is.na(resdf$beta_main) & resdf$ICD10!="I22",] # This is because the model didn't converged for I22
resdf2 <- melt(resdf,measure.vars = c("logp_adj","logp_main"))
resdf2$ICD10label[resdf2$value < -log10(0.05 / length(icd.all)) | resdf2$variable=="logp_main"] <- ""
ggplot(resdf2, aes(ICD10,value)) + 
  geom_point(aes(size=exp(beta_main),col=variable)) + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_hline(yintercept = -log10(0.05 / length(icd.all)), size = 1.5, color = "red") + 
  labs(x = "ICD Code", y = "-log10 Transformed P-Value") + 
  geom_text_repel(aes(label = resdf2$ICD10label),size=3,segment.alpha=0.5) + expand_limits(y = 200) +
  scale_size_continuous("HR",breaks=c(2,5,10,15,20,25,30), labels=c("2","5","10","15","20","25","30")) +
  scale_colour_discrete("",breaks=c("logp_adj","logp_main"), labels=c("Adjusted","Unadjusted"))



## PLOT ASSOCIATION BETWEEN PRS for MI AND ICD-CODES IN INDIVIDUALS WITHOUT MI ###
resdf <- read.table(file = paste0(path,'survivalanalysis_post_sensitivity.tsv'), sep = "\t", header = TRUE, stringsAsFactors = F)
icd.sel <- resdf$ICD10[2*pnorm(-abs(resdf$z_main)) < (0.05/nrow(resdf))]
resdfS <- resdf[resdf$ICD10 %in% icd.sel, ]
resdfS$p_icd_code <- 2*pnorm(-abs(resdfS$z_prs_icd))
resdfS$logp <- -log10(resdfS$p_icd_code)
resdfS$ICD10label <- resdfS$ICD10
resdfS$ICD10label[resdfS$logp < -log10(0.05 / length(icd.all))] <- ""
ggplot(resdfS, aes(ICD10,logp)) + 
  geom_point() + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_hline(yintercept = -log10(0.05 / length(icd.all)), size = 1.5, color = "red") + 
  labs(x = "ICD Code", y = "-log10 Transformed P-Value") + 
  geom_text_repel(aes(label = resdfS$ICD10label),size=3,segment.alpha=0.5)


# PLOT ASSOCIATION BETWEEN PRS for MI and ICD-CODES IN INDIVIDUALS WITH MI #### 
resdf <- read.table(file = paste0(path,'survivalanalysis_post_sensitivity.tsv'), sep = "\t", header = TRUE, stringsAsFactors = F)
icd.sel <- resdf$ICD10[2*pnorm(-abs(resdf$z_main)) < (0.05/nrow(resdf))]
resdfS <- resdf[resdf$ICD10 %in% icd.sel, ]
resdfS$p_just_mi <- 2*pnorm(-abs(resdfS$z_just_mi))
resdfS$logp <- -log10(resdfS$p_just_mi)
resdfS$ICD10label <- resdfS$ICD10
resdfS$ICD10label[resdfS$logp < -log10(0.05 / length(icd.all))] <- ""
ggplot(resdfS, aes(ICD10,logp)) + 
  geom_point() + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_hline(yintercept = -log10(0.05 / length(icd.all)), size = 1.5, color = "red") + 
  labs(x = "ICD Code", y = "-log10 Transformed P-Value") + 
  geom_text_repel(aes(label = resdfS$ICD10label),size=3,segment.alpha=0.5)




## PLOT INTERACTION BETWEEN PRS and MI for association with ICD ##
resdf <- read.table(file = paste0(path,'survivalanalysis_post_sensitivity.tsv'), sep = "\t", header = TRUE, stringsAsFactors = F)
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




## PLOT AVERAGE PRS FOR MI, NO MI, ICD, NO ICD
plot_icd_mi_prs("I20",total.new,total.newMI)


#### PLOT ASSOCIATION BETWEEN I20 and MI as function of the PRS ####
plot_icd_mi_prs_int("I20",total.new,total.newMI)


