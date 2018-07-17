##################################################
# File containing functions to use in analyses
##################################################

##################################################
# Function to load data, merge, and return a dataset for analyses
##################################################

load_data <- function(path) {
  
  # install.packages("data.table")
  library(data.table)

  # get data
  total1 = fread(file = paste0(path, 'hesin_registry_assess_cent_all_v2.csv'), sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  
  # polygenic risk score
  prs = fread(file = paste0(path, 'UKB_CardioCAD_PRS.txt'), sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  
  # merge phenotype data, polygenic risk score, and SNP matrix
  total = merge(total1, prs, by.x = "eid", by.y = "s", all.x = TRUE)
  
  return(data.frame(total))

}