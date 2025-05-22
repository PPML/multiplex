library(data.table)
library(foreach)
library(doParallel)

setwd("/scratch/users/mreitsma/ppml14/")

source("R/functions/utilities.R")
source("input/parameters.R")

crn_folder <- ("/scratch/users/mreitsma/clearance2_crn_draws/")

args <- commandArgs(trailingOnly=TRUE)

s_id <- as.numeric(args[1])

sens_list <- c("baseline", "test_coverage_low", "test_coverage_high", "init_high_25", "init_low_25")
i_df <- as.data.table(expand.grid(stochastic_num = 1:30, sens = sens_list))
stochastic_num <- as.numeric(i_df$stochastic_num[s_id])
sens <- as.character(i_df$sens[s_id])

p_ltfu_hiv_base <- rate_to_prob(prob_to_rate(0.3), 12)

ltfu_levels <- c(rate_to_prob(prob_to_rate(0.3), 12), rate_to_prob(prob_to_rate(.22), 12), rate_to_prob(prob_to_rate(0.14), 12), 
                 rate_to_prob(prob_to_rate(0.0625), 12), rate_to_prob(prob_to_rate(0.015), 12))

  for (param_level in 1:5) {
    
  p_ltfu_hiv <- ltfu_levels[param_level]                  
  
  set.seed(12345+stochastic_num)
                  
  load(paste0(crn_folder, sens, "/fixed/", stochastic_num, ".RData"))

## HIV Progression ----
set.seed(89101112+stochastic_num)
hiv_test_mat <- matrix(runif((max_time*nrow(df)), 0, 1), nrow = nrow(df))

set.seed(11121314+stochastic_num)
hiv_treat_mat <- matrix(runif((n_test_allowed*nrow(df)), 0, 1), nrow = nrow(df))
hiv_treat_mat <- ifelse(hiv_treat_mat<p_treat_hiv, 1, 0)

set.seed(13141516+stochastic_num)
hiv_ltfu_mat <- matrix(runif((max_time*nrow(df)), 0, 1), nrow = nrow(df))

set.seed(33333333+stochastic_num)
hiv_permltfu_mat <- matrix(runif((n_test_allowed*nrow(df)), 0, 1), nrow = nrow(df))

hiv_test_df <- copy(df)
hiv_test_df <- hiv_test_df[,.(pid, time_enter)]
for (i in 1:n_test_allowed) {
  hiv_test_df <- hiv_test_df[, paste0("hiv_test_time_", i):=apply(hiv_test_mat, 1, function(x) which(x < p_test_hiv)[i])]
}
hiv_test_df <- melt(hiv_test_df, id.vars = c("pid", "time_enter"), value.name = "time_hiv_test")
hiv_test_df <- hiv_test_df[, test_num:=seq_len(.N), by = "pid"]
hiv_test_df <- hiv_test_df[!is.na(time_hiv_test)]
hiv_test_df <- hiv_test_df[, variable:=NULL]

hiv_treat_df <- as.data.table(hiv_treat_mat)
colnames(hiv_treat_df) <- c(paste0(1:n_test_allowed))
hiv_treat_df <- hiv_treat_df[, pid:=1:nrow(df)]
hiv_treat_df <- melt(hiv_treat_df, id.vars = "pid", value.name = "hiv_treat", variable.name = "test_num")
hiv_treat_df <- hiv_treat_df[, test_num:=as.numeric(test_num)]
hiv_treat_df <- merge(hiv_test_df, hiv_treat_df, by = c("pid", "test_num"))
hiv_treat_df <- hiv_treat_df[hiv_treat==1, treat_num:=seq_len(.N), by = "pid"]

hiv_permltfu_df <- copy(df)
hiv_permltfu_df <- hiv_permltfu_df[,.(pid)]
hiv_permltfu_df <- hiv_permltfu_df[, permltfu_hiv:=apply(hiv_permltfu_mat, 1, function(x) which(x < (p_permltfu_hiv))[1])]

## LTFU INTERVENTION
hiv_ltfu_df <- copy(df)
hiv_ltfu_df <- hiv_ltfu_df[,.(pid)]
for (i in 1:(n_test_allowed)) {
  hiv_ltfu_df <- hiv_ltfu_df[, paste0("hiv_ltfu_time_", (i)):=apply(hiv_ltfu_mat, 1, function(x) which(x < (p_ltfu_hiv))[i])]
}

hiv_ltfu_df <- melt(hiv_ltfu_df, id.vars = c("pid"), value.name = "time_hiv_ltfu")
hiv_ltfu_df <- hiv_ltfu_df[, treat_num:=seq_len(.N), by = "pid"]
hiv_ltfu_df <- hiv_ltfu_df[order(pid, treat_num)]
hiv_ltfu_df <- hiv_ltfu_df[, l1:=shift(time_hiv_ltfu), by = "pid"]
hiv_ltfu_df <- hiv_ltfu_df[!is.na(l1), time_hiv_ltfu:=time_hiv_ltfu-l1]
hiv_ltfu_df <- hiv_ltfu_df[is.na(time_hiv_ltfu), time_hiv_ltfu:=censor_time]
hiv_ltfu_df <- hiv_ltfu_df[, c("variable", "l1"):=NULL]

## LTFU BASELINE
base_hiv_ltfu_df <- copy(df)
base_hiv_ltfu_df <- base_hiv_ltfu_df[,.(pid)]
for (i in 1:(n_test_allowed)) {
  base_hiv_ltfu_df <- base_hiv_ltfu_df[, paste0("hiv_ltfu_time_", (i)):=apply(hiv_ltfu_mat, 1, function(x) which(x < (p_ltfu_hiv_base))[i])]
}

base_hiv_ltfu_df <- melt(base_hiv_ltfu_df, id.vars = c("pid"), value.name = "time_hiv_ltfu_base")
base_hiv_ltfu_df <- base_hiv_ltfu_df[, treat_num:=seq_len(.N), by = "pid"]
base_hiv_ltfu_df <- base_hiv_ltfu_df[order(pid, treat_num)]
base_hiv_ltfu_df <- base_hiv_ltfu_df[, l1:=shift(time_hiv_ltfu_base), by = "pid"]
base_hiv_ltfu_df <- base_hiv_ltfu_df[!is.na(l1), time_hiv_ltfu_base:=time_hiv_ltfu_base-l1]
base_hiv_ltfu_df <- base_hiv_ltfu_df[is.na(time_hiv_ltfu_base), time_hiv_ltfu_base:=censor_time]
base_hiv_ltfu_df <- base_hiv_ltfu_df[, c("variable", "l1"):=NULL]

hiv_treat_df <- merge(hiv_treat_df, hiv_ltfu_df, by = c("pid", "treat_num"), all.x = T)
hiv_treat_df <- merge(hiv_treat_df, base_hiv_ltfu_df, by = c("pid", "treat_num"), all.x = T)
hiv_treat_df <- hiv_treat_df[order(pid, test_num)]
hiv_treat_df <- hiv_treat_df[, time_hiv_test_l1:=shift(time_hiv_test), by = "pid"]
hiv_treat_df <- hiv_treat_df[test_num>1, time_hiv_test:=time_hiv_test-time_hiv_test_l1]
hiv_treat_df <- hiv_treat_df[, c("time_hiv_test_l1", "treat_num"):=NULL]

hiv_treat_df <- merge(hiv_treat_df, hiv_permltfu_df, by = "pid", all.x=T)

hiv_treat_df <- hiv_treat_df[hiv_treat==1, hiv_treat_num:=seq_len(.N), by = "pid"]

baseline_df <- hiv_treat_df[,.(pid, test_num, time_enter, time_hiv_test, hiv_treat, time_hiv_ltfu_base)]
setnames(baseline_df, "time_hiv_ltfu_base", "time_hiv_ltfu")

baseline_df <- baseline_df[test_num==1, time_hiv_test_true:=time_hiv_test]
baseline_df <- baseline_df[, lag_time_ltfu:=shift(time_hiv_ltfu), by = "pid"]
baseline_df <- baseline_df[is.na(lag_time_ltfu), lag_time_ltfu:=0]

for (i in 1:30) {
  baseline_df <- baseline_df[, time_hiv_test_base:=shift(time_hiv_test_true), by = "pid"]
  baseline_df <- baseline_df[test_num==(i+1), time_hiv_test_true:=time_hiv_test_base + lag_time_ltfu + time_hiv_test]
} 

baseline_df <- baseline_df[, time_hiv_test:=time_hiv_test_true]
baseline_df <- baseline_df[, c("time_hiv_test_true", "lag_time_ltfu", "time_hiv_test_base"):=NULL]
baseline_df <- baseline_df[hiv_treat==1]
baseline_df <- baseline_df[, time_block_number:=seq_len(.N), by = "pid"]
baseline_df <- merge(baseline_df, hiv_permltfu_df, by = "pid", all.x=T)
baseline_df <- baseline_df[time_block_number>permltfu_hiv, time_hiv_test:=censor_time]
baseline_df <- baseline_df[time_hiv_test < censor_time]
setnames(baseline_df, "time_hiv_test", "time_hiv_treat")

## HIV Progression ----
set.seed(18192021+stochastic_num)
hiv_mort_mat <- matrix(runif((max_time*nrow(df)), 0, 1), nrow = nrow(df))

#### COMPLETE HISTORY FOR BASELINE LTFU
hiv_mort_df <- copy(df)
hiv_mort_df <- hiv_mort_df[, .(pid, ever_test_hiv)]
hiv_mort_df <- hiv_mort_df[, time_hiv_death := as.numeric(NA)]
hiv_mort_df <- hiv_mort_df[, c("cd4", "min_cd4"):=870]
hiv_mort_df <- hiv_mort_df[, time_art:=0]
hiv_mort_df <- hiv_mort_df[, time_no_art:=1]
hiv_mort_df <- hiv_mort_df[, test_num:=1]

cd4_art_df <- NULL
for (i in 1:max_time) {
  hiv_mort_df <- hiv_mort_df[, c("time_enter", "time_hiv_test", "hiv_treat", "time_hiv_ltfu", "time_hiv_ltfu_base",
                                 "permltfu_hiv", "hiv_treat_num", "art", "hiv_mort_prob", "hiv_mort_val", "hiv_mort", "time", "status"):=NULL]
  hiv_mort_df <- merge(hiv_mort_df, hiv_treat_df, by = c("pid", "test_num"), all.x=T)
  ## If test and treat
  hiv_mort_df <- hiv_mort_df[ever_test_hiv==1 & time_no_art==time_hiv_test & hiv_treat==1 & (hiv_treat_num<=permltfu_hiv | is.na(permltfu_hiv)), status:="Test and Treat"]
  hiv_mort_df <- hiv_mort_df[status=="Test and Treat", art:=1]
  hiv_mort_df <- hiv_mort_df[status=="Test and Treat", time_art:=1]
  hiv_mort_df <- hiv_mort_df[status=="Test and Treat", time_no_art:=0]
  ## If test and treat - but ltfu permanently
  hiv_mort_df <- hiv_mort_df[(is.na(status) & ever_test_hiv==1 & time_no_art>=time_hiv_test & hiv_treat==1 & (hiv_treat_num>permltfu_hiv & !is.na(permltfu_hiv))) | test_num>30, status:="Permanent LTFU"]
  hiv_mort_df <- hiv_mort_df[status=="Permanent LTFU", art:=0]
  hiv_mort_df <- hiv_mort_df[status=="Permanent LTFU", time_art:=0]
  hiv_mort_df <- hiv_mort_df[status=="Permanent LTFU", time_no_art:=time_no_art+1]
  ## If test but no treat
  hiv_mort_df <- hiv_mort_df[is.na(status) & ever_test_hiv==1 & time_no_art==time_hiv_test & hiv_treat==0 & is.na(art), status:="Test No Treat"]
  hiv_mort_df <- hiv_mort_df[status=="Test No Treat", art:=0]
  hiv_mort_df <- hiv_mort_df[status=="Test No Treat", test_num:=test_num+1]
  hiv_mort_df <- hiv_mort_df[status=="Test No Treat", time_art:=0]
  hiv_mort_df <- hiv_mort_df[status=="Test No Treat", time_no_art:=1]
  ## If no test
  hiv_mort_df <- hiv_mort_df[is.na(status) & ever_test_hiv==1 & time_no_art<time_hiv_test & time_no_art>0 & is.na(art), status:="No Test"]
  hiv_mort_df <- hiv_mort_df[status=="No Test", art:=0]
  hiv_mort_df <- hiv_mort_df[status=="No Test", time_art:=0]
  hiv_mort_df <- hiv_mort_df[status=="No Test", time_no_art:=time_no_art+1]
  ## No test ever
  hiv_mort_df <- hiv_mort_df[ever_test_hiv==0, status:="Never Test"]
  hiv_mort_df <- hiv_mort_df[status=="Never Test", art:=0]
  hiv_mort_df <- hiv_mort_df[status=="Never Test", time_art:=0]
  hiv_mort_df <- hiv_mort_df[status=="Never Test", time_no_art:=time_no_art+1]
  ## If on treatment
  ## Stop treatment
  hiv_mort_df <- hiv_mort_df[is.na(status) & ever_test_hiv==1 & time_art>0 & time_art>time_hiv_ltfu_base, status:="LTFU from Treatment"]
  hiv_mort_df <- hiv_mort_df[status == "LTFU from Treatment", art:=0]
  hiv_mort_df <- hiv_mort_df[status == "LTFU from Treatment", time_no_art:=1]
  hiv_mort_df <- hiv_mort_df[status == "LTFU from Treatment", test_num:=test_num+1]
  hiv_mort_df <- hiv_mort_df[status == "LTFU from Treatment", time_art:=0]
  ## Continue Treatment
  hiv_mort_df <- hiv_mort_df[is.na(status) & ever_test_hiv==1 & time_art>0 & time_art<=time_hiv_ltfu_base, status:="Continue Treatment"]
  hiv_mort_df <- hiv_mort_df[status=="Continue Treatment", art:=1]
  hiv_mort_df <- hiv_mort_df[status=="Continue Treatment", time_no_art:=0]
  hiv_mort_df <- hiv_mort_df[status=="Continue Treatment", time_art:=time_art+1]
  
  ## CD4 Decline if Untreated
  hiv_mort_df <- hiv_mort_df[art == 0, cd4:=pmin(cd4-(14.1/3), min_cd4)]
  hiv_mort_df <- hiv_mort_df[art == 1 & time_art %in% c(1:6), cd4:=cd4+(68/3)]
  hiv_mort_df <- hiv_mort_df[art == 1 & time_art %in% c(7:36), cd4:=cd4+(40/3)]
  hiv_mort_df <- hiv_mort_df[cd4<0, cd4:=0]
  hiv_mort_df <- hiv_mort_df[, min_cd4:=pmin(cd4, min_cd4)]
  hiv_mort_df <- hiv_mort_df[min_cd4<50, cd4:=pmin(cd4, 410)]
  hiv_mort_df <- hiv_mort_df[min_cd4<=200 & min_cd4>=50, cd4:=pmin(cd4, 548)]
  hiv_mort_df <- hiv_mort_df[min_cd4<=350 & min_cd4>200, cd4:=pmin(cd4, 660)]
  hiv_mort_df <- hiv_mort_df[min_cd4<=500 & min_cd4>350, cd4:=pmin(cd4, 780)]
  hiv_mort_df <- hiv_mort_df[min_cd4>500, cd4:=pmin(cd4, 870)]
  
  hiv_mort_df <- hiv_mort_df[, hiv_mort_prob:=rate_to_prob((.977^cd4)*0.59*(1+0.34*2)*2.03, 12)]
  hiv_mort_df <- hiv_mort_df[, hiv_mort_val:=hiv_mort_mat[, i]]
  
  hiv_mort_df <- hiv_mort_df[is.na(time_hiv_death), hiv_mort:=ifelse(hiv_mort_val < hiv_mort_prob, 1, 0)]
  hiv_mort_df <- hiv_mort_df[is.na(time_hiv_death) & hiv_mort==1, time_hiv_death:=i]
  
  temp <- hiv_mort_df[, time:=i]
  temp <- temp[is.na(time_hiv_death)]
  temp <- temp[,.(pid, cd4, art, time, status)]
  cd4_art_df <- rbind(cd4_art_df, temp, fill = T)
}

hiv_mort_df <- hiv_mort_df[, .(pid, time_hiv_death)]
hiv_mort_df <- hiv_mort_df[is.na(time_hiv_death), time_hiv_death:=censor_time]

hiv_mort_df_baseline <- copy(hiv_mort_df)
cd4_art_df_baseline <- copy(cd4_art_df)

## Seed HIV ----
set.seed(2001+stochastic_num)

seed_hiv_df <- copy(hiv_mort_df_baseline)
seed_hiv_df <- merge(seed_hiv_df, df[,.(pid, hiv, age, ever_test_hiv)], by = "pid")

## EVER TEST
if (nrow(seed_hiv_df[hiv==1 & age <= 24])>0) {
  seed_hiv_df <- seed_hiv_df[hiv==1 & age <= 24, duration_hiv:=sample(0:pmin((4*(1/timestep)), time_hiv_death-(2/timestep)), 1), by = "pid"]
}
if (nrow(seed_hiv_df[hiv==1 & age >=25 & age <= 29])>0) {
  seed_hiv_df <- seed_hiv_df[hiv==1 & age >=25 & age <= 29, duration_hiv:=sample(0:pmin((7*(1/timestep)), time_hiv_death-(2/timestep)), 1), by = "pid"]
}
if (nrow(seed_hiv_df[hiv==1 & age >=30 & age <= 39])>0) {
  seed_hiv_df <- seed_hiv_df[hiv==1 & age >=30 & age <= 39, duration_hiv:=sample(0:pmin((11*(1/timestep)), time_hiv_death-(2/timestep)), 1), by = "pid"]
}
if (nrow(seed_hiv_df[hiv==1 & age >=40])>0) {
  seed_hiv_df <- seed_hiv_df[hiv==1 & age >=40, duration_hiv:=sample(0:pmin((19*(1/timestep)), time_hiv_death-(2/timestep)), 1), by = "pid"]
}
seed_hiv_df <- seed_hiv_df[,.(pid, duration_hiv)]

#### COMPLETE HISTORY FOR INFECTED AT T-0
hiv_mort_df <- copy(df)
hiv_mort_df <- hiv_mort_df[, .(pid, ever_test_hiv)]
hiv_mort_df <- merge(hiv_mort_df, seed_hiv_df, by = "pid")
hiv_mort_df <- hiv_mort_df[!is.na(duration_hiv)]
hiv_mort_df <- hiv_mort_df[, time_hiv_death := as.numeric(NA)]
hiv_mort_df <- hiv_mort_df[, c("cd4", "min_cd4"):=870]
hiv_mort_df <- hiv_mort_df[, time_art:=0]
hiv_mort_df <- hiv_mort_df[, time_no_art:=1]
hiv_mort_df <- hiv_mort_df[, test_num:=1]

cd4_art_df <- NULL
# for (i in 1:66) {
for (i in 1:max_time) {
  hiv_mort_df <- hiv_mort_df[, c("time_enter", "time_hiv_test", "hiv_treat", "time_hiv_ltfu", "time_hiv_ltfu_base",
                                 "permltfu_hiv", "hiv_treat_num", "art", "hiv_mort_prob", "hiv_mort_val", "hiv_mort", "time", "status"):=NULL]
  hiv_mort_df <- merge(hiv_mort_df, hiv_treat_df, by = c("pid", "test_num"), all.x=T)
  ## If test and treat
  hiv_mort_df <- hiv_mort_df[ever_test_hiv==1 & time_no_art==time_hiv_test & hiv_treat==1 & (hiv_treat_num<=permltfu_hiv | is.na(permltfu_hiv)), status:="Test and Treat"]
  hiv_mort_df <- hiv_mort_df[status=="Test and Treat", art:=1]
  hiv_mort_df <- hiv_mort_df[status=="Test and Treat", time_art:=1]
  hiv_mort_df <- hiv_mort_df[status=="Test and Treat", time_no_art:=0]
  ## If test and treat - but ltfu permanently
  hiv_mort_df <- hiv_mort_df[(is.na(status) & ever_test_hiv==1 & time_no_art>=time_hiv_test & hiv_treat==1 & (hiv_treat_num>permltfu_hiv & !is.na(permltfu_hiv))) | test_num>30, status:="Permanent LTFU"]
  hiv_mort_df <- hiv_mort_df[status=="Permanent LTFU", art:=0]
  hiv_mort_df <- hiv_mort_df[status=="Permanent LTFU", time_art:=0]
  hiv_mort_df <- hiv_mort_df[status=="Permanent LTFU", time_no_art:=time_no_art+1]
  ## If test but no treat
  hiv_mort_df <- hiv_mort_df[is.na(status) & ever_test_hiv==1 & time_no_art==time_hiv_test & hiv_treat==0 & is.na(art), status:="Test No Treat"]
  hiv_mort_df <- hiv_mort_df[status=="Test No Treat", art:=0]
  hiv_mort_df <- hiv_mort_df[status=="Test No Treat", test_num:=test_num+1]
  hiv_mort_df <- hiv_mort_df[status=="Test No Treat", time_art:=0]
  hiv_mort_df <- hiv_mort_df[status=="Test No Treat", time_no_art:=1]
  ## If no test
  hiv_mort_df <- hiv_mort_df[is.na(status) & ever_test_hiv==1 & time_no_art<time_hiv_test & time_no_art>0 & is.na(art), status:="No Test"]
  hiv_mort_df <- hiv_mort_df[status=="No Test", art:=0]
  hiv_mort_df <- hiv_mort_df[status=="No Test", time_art:=0]
  hiv_mort_df <- hiv_mort_df[status=="No Test", time_no_art:=time_no_art+1]
  ## No test ever
  hiv_mort_df <- hiv_mort_df[ever_test_hiv==0, status:="Never Test"]
  hiv_mort_df <- hiv_mort_df[status=="Never Test", art:=0]
  hiv_mort_df <- hiv_mort_df[status=="Never Test", time_art:=0]
  hiv_mort_df <- hiv_mort_df[status=="Never Test", time_no_art:=time_no_art+1]
  ## If on treatment
    ## Stop treatment
    hiv_mort_df <- hiv_mort_df[is.na(status) & ever_test_hiv==1 & time_art>0 & time_art>time_hiv_ltfu_base & ((i < duration_hiv) | (i > (duration_hiv+120))), status:="LTFU from Treatment"]
    hiv_mort_df <- hiv_mort_df[is.na(status) & ever_test_hiv==1 & time_art>0 & time_art>time_hiv_ltfu & ((i >= duration_hiv ) & (i <= (duration_hiv+120))), status:="LTFU from Treatment"]
    hiv_mort_df <- hiv_mort_df[status == "LTFU from Treatment", art:=0]
    hiv_mort_df <- hiv_mort_df[status == "LTFU from Treatment", time_no_art:=1]
    hiv_mort_df <- hiv_mort_df[status == "LTFU from Treatment", test_num:=test_num+1]
    hiv_mort_df <- hiv_mort_df[status == "LTFU from Treatment", time_art:=0]
    ## Continue Treatment
    hiv_mort_df <- hiv_mort_df[(is.na(status) & ever_test_hiv==1 & time_art>0 & time_art<=time_hiv_ltfu_base & ((i < duration_hiv) | (i > (duration_hiv+120)) | (p_ltfu_hiv == p_ltfu_hiv_base))), status:="Continue Treatment"]
    hiv_mort_df <- hiv_mort_df[(is.na(status) & ever_test_hiv==1 & time_art>0 & time_art<=time_hiv_ltfu & ((i >= duration_hiv ) & (i <= (duration_hiv+120)) & (p_ltfu_hiv != p_ltfu_hiv_base))), status:="Continue Treatment"]
    hiv_mort_df <- hiv_mort_df[status=="Continue Treatment", art:=1]
    hiv_mort_df <- hiv_mort_df[status=="Continue Treatment", time_no_art:=0]
    hiv_mort_df <- hiv_mort_df[status=="Continue Treatment", time_art:=time_art+1]

  ## CD4 Decline if Untreated
  hiv_mort_df <- hiv_mort_df[art == 0, cd4:=pmin(cd4-(14.1/3), min_cd4)]
  hiv_mort_df <- hiv_mort_df[art == 1 & time_art %in% c(1:6), cd4:=cd4+(68/3)]
  hiv_mort_df <- hiv_mort_df[art == 1 & time_art %in% c(7:36), cd4:=cd4+(40/3)]
  hiv_mort_df <- hiv_mort_df[cd4<0, cd4:=0]
  hiv_mort_df <- hiv_mort_df[, min_cd4:=pmin(cd4, min_cd4)]
  hiv_mort_df <- hiv_mort_df[min_cd4<50, cd4:=pmin(cd4, 410)]
  hiv_mort_df <- hiv_mort_df[min_cd4<=200 & min_cd4>=50, cd4:=pmin(cd4, 548)]
  hiv_mort_df <- hiv_mort_df[min_cd4<=350 & min_cd4>200, cd4:=pmin(cd4, 660)]
  hiv_mort_df <- hiv_mort_df[min_cd4<=500 & min_cd4>350, cd4:=pmin(cd4, 780)]
  hiv_mort_df <- hiv_mort_df[min_cd4>500, cd4:=pmin(cd4, 870)]
  
  hiv_mort_df <- hiv_mort_df[, hiv_mort_prob:=rate_to_prob((.977^cd4)*0.59*(1+0.34*2)*2.03, 12)]
  hiv_mort_df <- hiv_mort_df[, hiv_mort_val:=hiv_mort_mat[unique(hiv_mort_df$pid), i]]
  
  hiv_mort_df <- hiv_mort_df[is.na(time_hiv_death), hiv_mort:=ifelse(hiv_mort_val < hiv_mort_prob, 1, 0)]
  hiv_mort_df <- hiv_mort_df[is.na(time_hiv_death) & hiv_mort==1, time_hiv_death:=i]
  
  temp <- hiv_mort_df[, time:=i]
  temp <- temp[is.na(time_hiv_death)]
  temp <- temp[,.(pid, cd4, art, time, status)]
  cd4_art_df <- rbind(cd4_art_df, temp, fill = T)
}

hiv_mort_df <- hiv_mort_df[, .(pid, time_hiv_death)]
hiv_mort_df <- hiv_mort_df[is.na(time_hiv_death), time_hiv_death:=censor_time]

hiv_mort_df_seed <- copy(hiv_mort_df)
cd4_art_df_seed <- copy(cd4_art_df)

  save(seed_hiv_df, hiv_mort_df_baseline, cd4_art_df_baseline,
       hiv_mort_df_seed, cd4_art_df_seed,
       hiv_treat_df, hiv_mort_mat, file=paste0(crn_folder, sens, "/hiv/", stochastic_num, "_ltfu_", round(p_ltfu_hiv, 4), ".RData"))
  }
