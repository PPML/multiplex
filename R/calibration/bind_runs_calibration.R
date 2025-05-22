library(data.table)
library(dplyr)
library(foreach)
library(doParallel)

setDTthreads(1)

setwd("/scratch/users/mreitsma/ppml14")
results_folder <- ("/scratch/users/mreitsma/clearance2_results/")

args <- commandArgs(trailingOnly=TRUE)

sherlock_table <- as.data.table(expand.grid(stochastic_num = 1:5,
                                            p_trans_hcv = seq(0.004, 0.01, .0005),
                                            rel_infect = seq(1, 2, .1),
                                            sx_rel = seq(2, 5, .5)))
sherlock_table <- sherlock_table[, sherlock_id_group:=ceiling(seq_len(.N)/50)]
sherlock_table <- sherlock_table[, sherlock_id:=seq_len(.N), by = "sherlock_id_group"]
sherlock_table <- sherlock_table[, edge_iter:=1]
sherlock_table <- sherlock_table[, scenario_id:=1]
sherlock_table <- sherlock_table[sherlock_id_group == args[1]]

temp_out <- NULL

for (s_id in 1:50) {                  
  if (s_id %in% sherlock_table$sherlock_id) {
    
  sherlock_table_run <- sherlock_table[sherlock_id==s_id]
  stochastic_num <- as.numeric(sherlock_table_run$stochastic_num)
  p_trans_hcv <- as.numeric(sherlock_table_run$p_trans_hcv)
  rel_infect <- as.numeric(sherlock_table_run$rel_infect)
  sx_rel <- as.numeric(sherlock_table_run$sx_rel)
                          
    temp <- fread(paste0(results_folder, "trans_calib/sim_", stochastic_num, "_phcv_", p_trans_hcv, "_relhcvhiv_", rel_infect, "_releqsx_", sx_rel, "_out_full.csv"))
    temp <- temp[analytic == 1]
    
    ## Incidence
    temp <- temp[hiv == 1, n_hiv:=.N, by = "time"]
    temp <- temp[hcv == 1, n_hcv:=.N, by = "time"]
    temp <- temp[hiv==1 & hcv==1, n_hiv_hcv:=.N, by = "time"]
    temp <- temp[art == 1 & ever_test_hiv == 1, n_treat_hiv:=.N, by = "time"]  
    temp <- temp[, hcv_diff:=hcv - shift(hcv), by = "pid"]
    temp <- temp[, hcv_times_diff:=times_hcv - shift(times_hcv), by = "pid"]
    temp <- temp[hcv_diff==1, new_hcv:=1]
    temp <- temp[hcv_diff==-1, hcv_cure:=1]
    temp <- temp[, ever_cure:=max(hcv_cure, na.rm=T), by = "pid"]
    temp <- temp[ever_cure==1, n_ever_cure:=.N, by = "time"]
    temp <- temp[, ever_hcv:=max(hcv, na.rm=T), by = "pid"]
    temp <- temp[ever_hcv==1, n_ever_hcv:=.N, by = "time"]
    temp <- temp[hcv_diff==1 & times_hcv>1, n_reinfect_hcv:=.N, by = "time"]
    temp <- temp[new_hcv==1, n_new_hcv:=.N, by = "time"]
    temp <- temp[hcv_cure==1, n_cure_hcv:=.N, by = "time"]
    temp <- temp[hcv_times_diff==1, new_hcv_incl_acute:=1]
    temp <- temp[new_hcv_incl_acute==1, new_hcv_incl_acute:=.N, by = "time"]
    temp <- temp[, hiv_diff:=hiv - shift(hiv), by = "pid"]
    temp <- temp[hiv_diff==1, new_hiv:=1]
    temp <- temp[new_hiv==1, n_new_hiv:=.N, by = "time"]
    
    ## Denominators
    temp <- temp[, n_pop_tot:=.N, by = "time"]
    temp <- temp[hcv==0 & current_pwid==1, n_hcv_sus_active:=.N, by = "time"]
    temp <- temp[current_pwid==1, n_active:=.N, by = "time"]
    temp <- temp[current_pwid==0 & permanent_quit == 0, n_former_temp:=.N, by = "time"]
    temp <- temp[permanent_quit == 1, n_former_permanent:=.N, by = "time"]
    temp <- temp[hiv==0, n_hiv_sus:=.N, by = "time"]
    temp <- temp[hcv==0, n_hcv_sus:=.N, by = "time"]
    
    ## QALYs
    ## HCV Utility
    temp <- temp[hcv == 1 & fibrosis_stage%in%c(0, 1, 2, 3), hcv_qol:= 0.96]
    temp <- temp[hcv == 0 & fibrosis_stage%in%c(0, 1, 2, 3), hcv_qol:= 1]
    temp <- temp[hcv == 1 & fibrosis_stage==4, hcv_qol:= .85]
    temp <- temp[hcv == 0 & fibrosis_stage==4, hcv_qol:= .96]
    temp <- temp[hcv == 1 & fibrosis_stage==5, hcv_qol:= .77]
    temp <- temp[hcv == 0 & fibrosis_stage==5, hcv_qol:= .93]
    
    ## HIV Utility
    temp <- temp[hiv==0 | is.na(cd4), hiv_qol:=1]
    temp <- temp[hiv==1 & cd4 >=500, hiv_qol:=1]
    temp <- temp[hiv==1 & cd4 < 500 & cd4>=200 & art == 1, hiv_qol:=.948]
    temp <- temp[hiv==1 & cd4 < 500 & cd4>=200 & (art == 0 | is.na(art)), hiv_qol:=0.92]
    temp <- temp[hiv==1 & cd4 < 200 & art == 1, hiv_qol:=.847]
    temp <- temp[hiv==1 & cd4 < 200 & (art == 0 | is.na(art)), hiv_qol:=0.764]
    
    ## PWID Utility
    temp <- temp[current_pwid == 1, pwid_qol:=0.68]
    temp <- temp[current_pwid == 0, pwid_qol:=.82]
    
    ## Background Utility
    temp <- temp[age<25, age_qol:=0.919]
    temp <- temp[age<35 & age>=25, age_qol:=.911]
    temp <- temp[age<45 & age>=35, age_qol:=.841]
    temp <- temp[age<55 & age>=45, age_qol:=.816]
    temp <- temp[age<65 & age>=55, age_qol:=.815]
    temp <- temp[age<75 & age>=65, age_qol:=.824]
    temp <- temp[age>=75, age_qol:=.811]
    
    temp <- temp[, utility:=age_qol*hiv_qol*pwid_qol*hcv_qol]
    temp <- temp[, sum_utility:=sum(utility), by = c("time")]
    
    temp <- temp[,.(scenario_id, stochastic_num, edge_iter_num, time, n_hiv, n_hcv, n_hiv_hcv, n_treat_hiv, n_new_hcv, n_reinfect_hcv, 
                    n_ever_cure, n_ever_hcv, new_hcv_incl_acute, n_cure_hcv, n_new_hiv, n_pop_tot, n_active, n_former_temp, n_former_permanent, 
                    n_hcv_sus_active, n_hiv_sus, n_hcv_sus, sum_utility)]
    temp <- unique(temp[, lapply(.SD, mean, na.rm=T), by = c("time")])
    temp <- temp[, p_trans_hcv:=p_trans_hcv]
    temp <- temp[, rel_infect:=rel_infect]
    temp <- temp[, sx_rel:=sx_rel]
    
    temp_out <- rbind(temp_out, temp, fill = T)
  }
}

fwrite(temp_out, paste0("/scratch/users/mreitsma/clearance2_summaries/trans_calib/", args[1], ".csv"), na = "", row.names = F)
