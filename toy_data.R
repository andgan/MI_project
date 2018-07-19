# make toy data
n = 1000000
total <- NULL
eid <- c(1:n)
sex <- sample(0:1, n, replace = TRUE)
prs <- sample(-10:10, n, replace = TRUE)
sbp <- rnorm(n, 120, 5)
smoke <- sample(0:1, n, replace = TRUE)
bmi <- rnorm(n, 20, 5)
diabetes <- sample(0:1, n, replace = TRUE)
lipid_lowering <- sample(0:1, n, replace = TRUE)
tot_chol <- rnorm(n, 200, 10)
hdl_chol <- rnorm(n, 60, 10)
total <- data.frame(eid, sex, prs, sbp, smoke, bmi, diabetes, lipid_lowering, tot_chol, hdl_chol)

total$birth_date <- sample(seq(as.Date('1950/01/01'), as.Date('1999/10/16'), by = "day"), n, replace = TRUE)
total$ac_date <- sample(seq(as.Date('1995/10/16'), as.Date('2005/11/14'), by = "day"), n, replace = TRUE)
total$age <- round(as.numeric((total$ac_date - total$birth_date) / 365.25), 2)
total$death <- sample(0:1, n, replace = TRUE)
total$death_date <- as.Date(sample(seq(as.Date('1972/11/12'), as.Date('2000/01/01'), by = "day"), n, replace = TRUE))
total$death_date[which(total$death == 0)] <- NA
total$mi <- sample(0:1, n, replace = TRUE)
total$mi_date <- as.Date(sample(seq(as.Date('1972/11/12'), as.Date('1999/10/16'), by = "day"), n, replace = TRUE))
total$mi_date[which(total$mi == 0)] <- NA
total$freq <- sample(1:10, n, replace = TRUE)
total <- total[rep(row.names(total), total$freq), ]

total$icd10 <- sample(toupper(c(letters[1:26], NA)), nrow(total), replace = TRUE)
total$icd10_date <- as.Date(sample(seq(as.Date('1972/11/12'), as.Date('2010/11/12'), by = "day"), nrow(total), replace = TRUE))
total$icd10_date[which(is.na(total$icd10))] <- NA