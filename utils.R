### FILE CONTAINING FUNCTIONS TO USE IN THE ANALYSIS ###


## PLOT AVERAGE PRS FOR MI, NO MI, ICD, NO ICD
plot_icd_mi_prs <- function(icdcode, alldata, allmi)
{
  
  alldataU <- alldata[!duplicated( alldata$eid), c("eid", "sex", "age", "death", "death_date", "endfollowup", "startfollowup","PRS_0.5")]
  
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
  
  df <- data.frame(prs_mean=c(mean(icd.mi$PRS_0.5,na.rm=T),mean(icd.no.mi$PRS_0.5,na.rm=T),mean(no.icd.mi$PRS_0.5,na.rm=T),mean(no.icd.no.mi$PRS_0.5,na.rm=T)), prs_se=c(sd(icd.mi$PRS_0.5,na.rm=T)/sqrt(nrow(icd.mi)),sd(icd.no.mi$PRS_0.5,na.rm=T)/sqrt(nrow(icd.no.mi)),sd(no.icd.mi$PRS_0.5,na.rm=T)/sqrt(nrow(no.icd.mi)),sd(no.icd.no.mi$PRS_0.5,na.rm=T)/sqrt(nrow(no.icd.no.mi))),mi=c("MI","NO MI","MI","NO MI"),ICD=c(icdcode,icdcode,paste("NO",icdcode),paste("NO",icdcode)))
  df$CIU <- df$prs_mean + 1.96*df$prs_se
  df$CID <- df$prs_mean - 1.96*df$prs_se
  
  ggplot(aes(y=prs_mean,x=mi), data=df) + geom_bar(aes(fill=ICD), stat="identity", position="dodge") + geom_errorbar(width=0.2,aes(ymin=CID, ymax=CIU, fill=ICD),position=position_dodge(.9) ) + theme_bw()
  
}

## PLOT ASSOCIATION BETWEEN THE ICD code and MI as function of the PRS 

plot_icd_mi_prs_int <- function(icdcode, alldata, allmi)
{
  
  alldataU <- alldata[!duplicated( alldata$eid), c("eid", "sex", "age", "death", "death_date", "endfollowup", "startfollowup","PRS_0.5")]
  
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
  findata$PRS_quartiles <- as.numeric(cut(findata$PRS_0.5,quantile(findata$PRS_0.5,seq(0,1,0.25),na.rm=T)))
  
  mod1 <- coxph(Surv(tstart, tstop, outP) ~ age + trt + sex, data = findata[findata$PRS_quartiles==1,])
  mod2 <- coxph(Surv(tstart, tstop, outP) ~ age + trt + sex, data = findata[findata$PRS_quartiles==2,])
  mod3 <- coxph(Surv(tstart, tstop, outP) ~ age + trt + sex, data = findata[findata$PRS_quartiles==3,])
  mod4 <- coxph(Surv(tstart, tstop, outP) ~ age + trt + sex, data = findata[findata$PRS_quartiles==4,])
  
  df <- data.frame(BETA=c(coef(mod1)[2],coef(mod2)[2],coef(mod3)[2],coef(mod4)[2]),
                   SE=c(coef(summary(mod1))[2,3],coef(summary(mod2))[2,3],coef(summary(mod3))[2,3],coef(summary(mod4))[2,3]),
                   P=c(coef(summary(mod1))[2,5],coef(summary(mod2))[2,5],coef(summary(mod3))[2,5],coef(summary(mod4))[2,5]),
                   X=c("1st","2nd","3rd","4th"))
  
  df$CIU <- exp(df$BETA+1.96*df$SE)
  df$CID <- exp(df$BETA-1.96*df$SE)
  
  ggplot(aes(y=exp(BETA),x=X), data=df) + geom_errorbar(width=.1, aes(ymin=CID, ymax=CIU), colour="red") + geom_point(shape=21, size=3, fill="white") + theme_bw() + ylab(paste("HR for association between", icdcode, "and MI")) + xlab("Quartiles of PRS for MI")
  
}

