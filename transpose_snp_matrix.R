## TRANSPOSE THE SNPFILE IN THE ID by VARIANT MATRIX ##
library(data.table)
library(reshape2)

# identify the data path
path <- "C:/Jiwoo Lee/Myocardial Infarction Research Project 2017/"
path <- "/Users/andreaganna/Documents/Work/Post_doc/jiwoo/"

# Load information about the SNPs
info_snp <- read.csv(file=paste0(path,"MI_project/all_gwsign_snps_CHD.csv"))

# Load original matrix
d <- fread(paste0(path,"variant_gw_cad_all.tsv"), header=T, stringsAsFactor=F)

# Keep only known SNPs in snp_mat
ss <- strsplit(d$V,":")
pos <- sapply(ss,"[",2)
chr <- sapply(ss,"[",1)
pos_chr <- paste0(sub("^[0]+", "", chr),"_",pos)

# Select only variants that are known
d2 <- d[pos_chr %in% paste0(info_snp$chromosome,"_",info_snp$position)[info_snp$known==1],]

# transpose, as numeric and save
dT <- acast(d2, ID ~ V, value.var = 'GT', fill = '0')
dTN <- data.frame(apply(dT,2,as.numeric))
dTN$eid <- rownames(dT)

write.table(dTN,file=paste0(path,"variant_gw_cad_allT.tsv"),sep="\t",row.names=F, col.names=T, quote=F)
