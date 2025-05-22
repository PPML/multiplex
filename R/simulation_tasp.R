library(EpiModel)
library(data.table)
library(dplyr)
library(foreach)
library(doParallel)

  setDTthreads(1)

  setwd("/scratch/users/mreitsma/ppml14")
  crn_folder <- ("/scratch/users/mreitsma/clearance2_crn_draws/")
  results_folder <- ("/scratch/users/mreitsma/clearance2_results/")
  
  source("R/functions/nw_construct.R")
  source("R/functions/transmission.R")
  source("R/functions/utilities.R")
  source("input/parameters.R")

  args <- commandArgs(trailingOnly=TRUE)
  sens <- as.character(args[2])
  sherlock_table_name <- as.character(args[3])
  n_loop <- as.numeric(args[4])
  years_run <- as.numeric(args[5])
  
  sherlock_table <- fread(paste0("input/sherlock/", sherlock_table_name))
  sherlock_table <- sherlock_table[sherlock_id_group == as.numeric(args[1])]
  
  if (sens == "hiv_tasp_low") {
    hiv_tasp_effect <- 0.25
  } else {
    hiv_tasp_effect <- 0
  }

  for (s_id in 1:n_loop) {                  
    if (s_id %in% sherlock_table$sherlock_id) {

  sherlock_table_run <- sherlock_table[sherlock_id==s_id]
  stochastic_num <- as.numeric(sherlock_table_run$stochastic_num)
  edge_iter <- as.numeric(sherlock_table_run$edge_iter)
  scenario_id <- as.numeric(sherlock_table_run$scenario_id)
  
  load(paste0(crn_folder, "est_net_eq.RData"))
  load(paste0(crn_folder, "est_net_sx.RData"))
  
  diss_coef_adj_eq <- diss_coef_ad(duration = (duration_days_eq/365)*(1/timestep), exit_rate = 0.004) #.002
  diss_coef_adj_sx <- diss_coef_ad(duration = (duration_days_sx/365)*(1/timestep), exit_rate = 0.004) #.004

  set.seed(12345+stochastic_num+(edge_iter*100))
  
  if (sens == "scaledrm_0.5") {
    load(paste0(crn_folder, "baseline/fixed/", stochastic_num, ".RData"))
    load(paste0(crn_folder, "baseline/compiled_main_scaledrm_0.5/variable_", stochastic_num, "_scenario_", scenario_id, ".RData"))
  } else if (sens == "scaledrm_0.75") {
    load(paste0(crn_folder, "baseline/fixed/", stochastic_num, ".RData"))
    load(paste0(crn_folder, "baseline/compiled_main_scaledrm_0.75/variable_", stochastic_num, "_scenario_", scenario_id, ".RData"))
  } else if (sens != "trans_calib") {
    load(paste0(crn_folder, sens, "/fixed/", stochastic_num, ".RData"))
    load(paste0(crn_folder, sens, "/compiled_main/variable_", stochastic_num, "_scenario_", scenario_id, ".RData"))
  } else {
    load(paste0(crn_folder, "baseline/fixed/", stochastic_num, ".RData"))
    load(paste0(crn_folder, "baseline/compiled_trans_calib/variable_", stochastic_num, "_scenario_", scenario_id, ".RData"))
    load(paste0(crn_folder, "baseline/trans_calib/", stochastic_num, "_phcv_", sherlock_table_run$p_trans_hcv, "_relhcvhiv_", sherlock_table_run$rel_infect, "_releqsx_", sherlock_table_run$sx_rel, ".RData"))
  }
  
  fibrosis_df[is.na(fibrosis_df)] <- censor_time
  
  hcv_trans_df <- hcv_trans_df[, times_hcv:=times_hcv-1]
  hcv_trans_df <- hcv_trans_df[pid%in%df$pid[df$times_hcv==1], never_tx_exposures:=shift(never_tx_exposures, type = "lag"), by = "pid"]
  setnames(hcv_trans_df, "times_hcv", "times_hcv_exposure_merge")
  
  df <- merge(df, seed_hcv_df, by = "pid")
  df <- merge(df, seed_hiv_df, by = "pid", all.x=T)
  df <- merge(df, hiv_trans_df, by = "pid")
  df <- merge(df, fibrosis_df, by = "pid")
  df <- merge(df, drm_df[,.(pid, time_drug_death)], by = "pid")
  df <- df[, time_background_death:=time_background_death+time_enter]  
  df <- merge(df, hiv_mort_df_seed[,.(pid, time_hiv_death)], by = "pid", all.x=T)
  
  cess_df <- cess_df[, time:=time-1]

  df <- df[, time_to_f2:=time_to_f1+time_to_f2]
  df <- df[, time_to_f3:=time_to_f2+time_to_f3]
  df <- df[, time_to_f4:=time_to_f3+time_to_f4]
  df <- df[, time_to_f5:=time_to_f4+time_to_f5]
  
  df <- df[hiv == 1, time_hiv_infection:=0]
  df <- df[hcv == 1, time_hcv_infection:=0]
  df <- df[is.na(duration_hcv), duration_hcv:=0]
  df <- df[is.na(duration_hiv), duration_hiv:=0]
  df <- df[, cod:=as.character(NA)]
  df <- df[, age:=as.numeric(age)]
  df <- df[, current_pwid:=NULL]
  df <- df[, hiv_sx_exposure_counter:=0]
  df <- df[, hiv_eq_exposure_counter:=0]
  df <- df[, hiv_gensx_exposure_counter:=0]
  df <- df[, hcv_exposure_counter:=0]
  df <- df[, hcv_test_num_intervention:=0]
  df <- df[, hcv_test_num_baseline:=0]
  df <- df[, hcv_times_treated:=0]
  df <- df[, times_hcv_exposure_merge:=times_hcv]
  
  full_df <- copy(df)
    
  df <- df[order(pid)]
  df <- df[, cid:=seq_len(.N)]
  df <- df[time_enter==0]

  sx_el <- return_edges_dt(est_net_sx$newnetwork, 0)
  eq_el <- return_edges_dt(est_net_eq$newnetwork, 0)
  
  cd4_art_df <- copy(cd4_art_df_seed)
  setnames(cd4_art_df, "time", "duration_hiv")
  setnames(hcv_test_df, "duration_hcv", "duration_hcv_floor")

  out_full <- NULL
  deaths <- NULL

  interventions_stop <- 120
  
  for (t in 1:(12*years_run)) {
    if (nrow(df) > 0) {
    # message(t)
    
    ## Clear old info
    df <- df[, c("ssp", "condom", "ssp_year", "eq_iso", "permanent_quit", "current_pwid", "ever_tx_exposures", "never_tx_exposures", "hcv_exposures",
                 "acute_clear", "hcv_treat", "time", "cd4", "art", "status", "hcv_treat_intervention",
                 "hcv_treat_baseline", "hcv_test_intervention", "hcv_test_baseline", "time_hiv_death_new"):=NULL]
    
    ## Age
    df <- df[, age:=age+timestep]
    
    ## New Entrants
    df <- rbind(df, full_df[time_enter==t], fill = T)
    
    df <- df[order(pid)]
    df <- df[, new_cid:=seq_len(.N)]
    cid_update <- df[,.(cid, new_cid)]
    sx_el <- update_cid(el = sx_el, cid_update = cid_update)
    eq_el <- update_cid(el = eq_el, cid_update = cid_update)
    df <- df[, cid:=new_cid]
    df <- df[, new_cid:=NULL]
    df <- df[, ssp_year:=ceiling(t/12)]
    df <- merge(df, ssp_df, by = c("pid", "ssp_year"), all.x=T)

    df <- merge(df, cess_df[time == t], by = "pid")
    df <- merge(df, hcv_trans_df[,.(pid, ever_tx_exposures, never_tx_exposures, times_hcv_exposure_merge)], by = c("pid", "times_hcv_exposure_merge"), all.x=T)
    df <- df[, hcv_exposures:=ifelse(hcv_times_treated>=1, ever_tx_exposures, never_tx_exposures)]
    
    ## Add to cumulative duration of infection
    df <- df[hiv==1, duration_hiv:=duration_hiv+1]
    df <- merge(df, cd4_art_df, by = c("pid", "duration_hiv"), all.x=T)
    df <- df[, duration_hcv:=as.numeric(duration_hcv)]
    df <- df[hcv==1 & hiv == 0, duration_hcv:=duration_hcv+1]
    df <- df[hcv==1 & hiv == 1 & art == 0, duration_hcv:=duration_hcv+2.5]
    df <- df[hcv==1 & hiv == 1 & art == 1, duration_hcv:=duration_hcv+1.7]
    df <- df[, duration_hcv_floor:=floor(duration_hcv)]
    
    if (t <= interventions_stop) {
      new_el <- attach_edges(data = df[permanent_quit==0], sx = sx_el, eq = eq_el, time = t, n_sx = 1, n_eq = 1)
      eq_el <- rbind(eq_el, new_el[[1]])
      sx_el <- rbind(sx_el, new_el[[2]])
    }
    
    df <- df[order(pid)]
    df <- df[permanent_quit==0, new_cid:=seq_len(.N)]
    cid_update <- df[,.(cid, new_cid)]
    sx_el <- update_cid(el = sx_el, cid_update = cid_update)
    eq_el <- update_cid(el = eq_el, cid_update = cid_update)
    df <- df[, cid:=new_cid]
    df <- df[, new_cid:=NULL]
    
    if (nrow(df[current_pwid==1]) >= 25) {
    ## Update nodal attributes
    sx_nwd <- nw_from_el(data = df[permanent_quit==0], el = sx_el)
    eq_nwd <- nw_from_el(data = df[permanent_quit==0], el = eq_el)
    
    ## Adjust edges coefficient to account for changing population size
    edge_coef_eq <- edges_ad(orig_coef = est_net_eq$coef.form[1], orig_n = 1000, new_n = nrow(df[permanent_quit==0]))
    edge_coef_sx <- edges_ad(orig_coef = est_net_sx$coef.form[1], orig_n = 1000, new_n = nrow(df[permanent_quit==0]))
    
    edgecov_coef_eq <- edgecov_ad(orig_coef = est_net_eq$coef.form[2], orig_n = 1000, new_n = nrow(df[permanent_quit==0]))
    edgecov_coef_sx <- edgecov_ad(orig_coef = est_net_sx$coef.form[2], orig_n = 1000, new_n = nrow(df[permanent_quit==0]))
    
    ## Simulate next time step of network configuration
    eq_nwd %n% "ec" <- as.matrix(sx_nwd)
    eq_nwd <- nw_sim(eq_nwd, est_net = est_net_eq, adj_edges = edge_coef_eq, adj_edgecov = edgecov_coef_eq, adj_diss = diss_coef_adj_eq, time.start = t, time.interval = 1, time.offset = 1)[,2:4]
    if (is.null(dim(eq_nwd))) {
      eq_nwd <- as.data.table(t(eq_nwd))
    } else {
      eq_nwd <- as.data.table(eq_nwd)
    }
    n_lost_eq <- nrow(eq_nwd[to==0])
    n_gained_eq <- nrow(eq_nwd[to==1])
    eq_el <- eq_el[!eq_nwd[to==0], on = .(tail, head)] # Remove edges
    eq_el <- rbind(eq_el, eq_nwd[to==1, c("tail", "head"), with = F]) # Add edges

    eq_nwd <- nw_from_el(data = df[permanent_quit==0], el = eq_el)
    sx_nwd %n% "ec" <- as.matrix(eq_nwd)
    sx_nwd <- nw_sim(sx_nwd, est_net = est_net_sx, adj_edges = edge_coef_sx, adj_edgecov = edgecov_coef_sx, adj_diss = diss_coef_adj_sx, time.start = t, time.interval = 1, time.offset = 1)[,2:4]
    if (is.null(dim(sx_nwd))) {
      sx_nwd <- as.data.table(t(sx_nwd))
    } else {
      sx_nwd <- as.data.table(sx_nwd)
    }
    n_lost_sx <- nrow(sx_nwd[to==0])
    n_gained_sx <- nrow(sx_nwd[to==1])
    sx_el <- sx_el[!sx_nwd[to==0], on = .(tail, head)] # Remove edges
    sx_el <- rbind(sx_el, sx_nwd[to==1, c("tail", "head"), with = F]) # Add edges
    }
    
    ## Infection transmission
    exposures <- count_exposures(df = df, sx_el = sx_el, eq_el = eq_el, t=t, ssp_effect = ssp_effect, condom_effect = condom_effect, hiv_equip_tasp = hiv_tasp_effect)
    df <- update_exposures(df, exposures)
    df <- update_infections(df, t)
    df <- merge(df, hcv_acute_df, by = c("times_hcv", "pid"), all.x=T)
    df <- df[time_hcv_infection==t & acute_clear==1 & hcv==1, hcv:=0]
    
    ## Test and Treat (HCV)
    df <- merge(df, hcv_test_df, by = c("pid", "duration_hcv_floor"), all.x=T)
    if (t <= interventions_stop) {
      df <- df[!is.na(hcv_test_intervention) & hcv==1 & hcv_test_intervention == 1, hcv_test_num_intervention:=hcv_test_num_intervention+hcv_test_intervention]
      df <- merge(df, hcv_treat_df, by = c("pid", "hcv_test_num_intervention"), all.x=T)
      df <- df[hcv==0 | ever_test_hcv==0 | hcv_test_intervention==0 | hcv_treat_intervention==0, hcv_treat:=0]
      df <- df[hcv==1 & ever_test_hcv==1 & hcv_test_intervention==1 & hcv_treat_intervention==1, hcv_treat:=1]
    } else {
      df <- df[!is.na(hcv_test_baseline) & hcv==1 & hcv_test_baseline == 1, hcv_test_num_baseline:=hcv_test_num_baseline+hcv_test_baseline]
      df <- merge(df, base_hcv_treat_df, by = c("pid", "hcv_test_num_baseline"), all.x=T)
      df <- df[hcv==0 | ever_test_hcv==0 | hcv_test_baseline==0 | hcv_treat_baseline==0, hcv_treat:=0]
      df <- df[hcv==1 & ever_test_hcv==1 & hcv_test_baseline==1 & hcv_treat_baseline==1, hcv_treat:=1]
    }
    df <- df[ever_test_hcv == 1 & hcv_treat==1 & hcv==1, hcv_times_treated:=hcv_times_treated+1]
    df <- df[ever_test_hcv == 1 & hcv_treat==1 & hcv==1 & hcv_times_treated==1, times_hcv_exposure_merge:=0]
    df <- df[ever_test_hcv == 1 & hcv_treat==1 & hcv==1 & hcv_times_treated==1, hcv_exposure_counter:=0]
    df <- df[ever_test_hcv == 1 & hcv_treat==1 & hcv==1, hcv:=0]
    
    ## HCV Progression
    df <- df[duration_hcv < time_to_f1, fibrosis_stage:=0]
    df <- df[duration_hcv >= time_to_f1 & duration_hcv < time_to_f2, fibrosis_stage:=1]
    df <- df[duration_hcv >= time_to_f2 & duration_hcv < time_to_f3, fibrosis_stage:=2]
    df <- df[duration_hcv >= time_to_f3 & duration_hcv < time_to_f4, fibrosis_stage:=3]
    df <- df[duration_hcv >= time_to_f4 & duration_hcv < time_to_f5, fibrosis_stage:=4]
    df <- df[duration_hcv >= time_to_f5 , fibrosis_stage:=5]
    
    if (nrow(df[time_hiv_infection==t]) > 0) {
      if (t <= interventions_stop & scenario_id > 9) {
        ## HIV Progression - New Infections
        hiv_mort_df <- df[time_hiv_infection==t]
        hiv_mort_df <- hiv_mort_df[, .(pid, ever_test_hiv, duration_hiv)]
        hiv_mort_df <- hiv_mort_df[, time_hiv_death := as.numeric(NA)]
        hiv_mort_df <- hiv_mort_df[, c("cd4", "min_cd4"):=870]
        hiv_mort_df <- hiv_mort_df[, time_art:=0]
        hiv_mort_df <- hiv_mort_df[, time_no_art:=1]
        hiv_mort_df <- hiv_mort_df[, test_num:=1]
        
        cd4_art_new <- NULL
        for (i in 1:(max_time-t)) {
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
          if (((t+i)>interventions_stop)) {
            hiv_mort_df <- hiv_mort_df[is.na(status) & ever_test_hiv==1 & time_art>0 & time_art>time_hiv_ltfu_base, status:="LTFU from Treatment"]
          } else {
            hiv_mort_df <- hiv_mort_df[is.na(status) & ever_test_hiv==1 & time_art>0 & time_art>time_hiv_ltfu, status:="LTFU from Treatment"]
          }
          hiv_mort_df <- hiv_mort_df[status == "LTFU from Treatment", art:=0]
          hiv_mort_df <- hiv_mort_df[status == "LTFU from Treatment", time_no_art:=1]
          hiv_mort_df <- hiv_mort_df[status == "LTFU from Treatment", test_num:=test_num+1]
          hiv_mort_df <- hiv_mort_df[status == "LTFU from Treatment", time_art:=0]
          ## Continue Treatment
          if (((t+i)>interventions_stop)) {
            hiv_mort_df <- hiv_mort_df[is.na(status) & ever_test_hiv==1 & time_art>0 & time_art<=time_hiv_ltfu_base, status:="Continue Treatment"]
          } else {
            hiv_mort_df <- hiv_mort_df[is.na(status) & ever_test_hiv==1 & time_art>0 & time_art<=time_hiv_ltfu, status:="Continue Treatment"]
          }
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
          cd4_art_new <- rbind(cd4_art_new, temp, fill = T)
        }
      }
      else {
        cd4_art_new <- cd4_art_df_baseline[pid%in%df$pid[df$hiv==1 & df$time_hiv_infection==t]]
        hiv_mort_df <- hiv_mort_df_baseline[pid%in%df$pid[df$hiv==1 & df$time_hiv_infection==t]]
      }
      setnames(cd4_art_new, "time", "duration_hiv")
      cd4_art_df <- rbind(cd4_art_df, cd4_art_new)
      
      hiv_mort_df <- hiv_mort_df[,.(pid, time_hiv_death)]
      setnames(hiv_mort_df, "time_hiv_death", "time_hiv_death_new")
      df <- merge(df, hiv_mort_df, by = "pid", all.x=T)
      df <- df[is.na(time_hiv_death) & !is.na(time_hiv_death_new), time_hiv_death:=time_hiv_death_new]
    }

    ## HCV Mortality
    df <- df[, hcv_mort_prob:=mapply(FUN = hcv_mort, hcv = hcv, fibrosis_stage = fibrosis_stage, return_prob = T), by = "pid"]
    df <- df[(hcv_mort_df[pid%in%df$pid][[t]]) < hcv_mort_prob, time_hcv_death:=t]    
    
    ## Cause of Death
    df <- df[t == time_background_death, cod:="Background"]
    df <- df[t == time_drug_death, cod:="Drug"]
    df <- df[duration_hiv == time_hiv_death & hiv==1, cod:="HIV"]
    df <- df[t == time_hcv_death, cod:="HCV"]
    
    # if (t %in% seq(3, 800, 3)) {
    #   temp_eq <- eq_el[head%in%df$cid[df$current_pwid==1 & df$eq_iso==0] & tail%in%df$cid[df$current_pwid==1 & df$eq_iso==0]]
    #   joint_el <- merge(temp_eq, sx_el)
    #   temp <- data.table(n_pop = nrow(df[is.na(cod)]), n_pop_active = nrow(df[is.na(cod) & permanent_quit==0]),
    #                      n_new_hiv = nrow(df[time_hiv_infection == t]),
    #                      n_new_hcv = nrow(df[time_hcv_infection == t]),
    #                      hcv_pct = nrow(df[permanent_quit==0 & is.na(cod) & hcv==1])/nrow(df[is.na(cod) & permanent_quit==0]),
    #                      hiv_pct = nrow(df[permanent_quit==0 & is.na(cod) & hiv==1])/nrow(df[permanent_quit==0 & is.na(cod)]),
    #                      n_eq = nrow(temp_eq)*2/nrow(df[current_pwid==1]), n_sx = nrow(sx_el)*2/nrow(df[permanent_quit==0]),
    #                      joint_eq = nrow(joint_el)/nrow(temp_eq),
    #                      joint_sx = nrow(joint_el)/nrow(sx_el), n_lost_sx = n_lost_sx, n_gained_sx = n_gained_sx,
    #                      n_lost_eq = n_lost_eq, n_gained_eq = n_gained_eq)
    #   print(temp)
    # }
    
    ## Summaries
    df <- df[, time:=t]
    df <- df[(cod%in%c("Background", "HIV", "Drug", "HCV")), time_death:=t]
    deaths <- rbind(deaths, df[(cod%in%c("Background", "HIV", "Drug", "HCV"))], fill = T)
    out_full <- rbind(out_full, df[!(cod%in%c("Background", "HIV", "Drug", "HCV"))], fill = T)

    ## UPDATE CID AND CORRESPONDING EDGELISTS
    df <- df[order(pid)]
    df <- df[!(cod%in%c("Background", "HIV", "Drug", "HCV")), new_cid:=seq_len(.N)]
    cid_update <- df[,.(cid, new_cid)]
    sx_el <- update_cid(el = sx_el, cid_update = cid_update)
    eq_el <- update_cid(el = eq_el, cid_update = cid_update)
    df <- df[, cid:=new_cid]
    df <- df[, new_cid:=NULL]
    df <- df[!(cod%in%c("Background", "HIV", "Drug", "HCV"))]
    }
  }
  # })
  
  out_full <- out_full[, stochastic_num:=stochastic_num]
  out_full <- out_full[, scenario_id:=scenario_id]
  out_full <- out_full[, edge_iter_num:=edge_iter]
  out_full <- out_full[, scale_drm:=1]
  if (sens == "scaledrm_0.5") {
    out_full <- out_full[, scale_drm:=.5]
  } else if (sens == "scaledrm_0.75") {
    out_full <- out_full[, scale_drm:=.75]
  } else {
    out_full <- out_full[, scale_drm:=1]
  }
  
  deaths <- deaths[, stochastic_num:=stochastic_num]
  deaths <- deaths[, scenario_id:=scenario_id]
  deaths <- deaths[, edge_iter_num:=edge_iter]
  if (sens == "scaledrm_0.5") {
    deaths <- deaths[, scale_drm:=.5]
  } else if (sens == "scaledrm_0.75") {
    deaths <- deaths[, scale_drm:=.75]
  } else {
    deaths <- deaths[, scale_drm:=1]
  }
  
  if (sens != "trans_calib") {
    fwrite(out_full, paste0(results_folder, sens, "/sim_", stochastic_num, "_scenario_", scenario_id, "_edge_", edge_iter, "_out_full.csv"), na = "", row.names = F)
    fwrite(deaths, paste0(results_folder, sens, "/sim_", stochastic_num, "_scenario_", scenario_id, "_edge_", edge_iter, "_deaths.csv"), na = "", row.names = F)
  } else {
    out_full <- out_full[, p_trans_hcv:=sherlock_table_run$p_trans_hcv]
    out_full <- out_full[, rel_infect:=sherlock_table_run$rel_infect]
    out_full <- out_full[, sx_rel:=sherlock_table_run$p_trans_hcv]
    
    fwrite(out_full, paste0(results_folder, sens, "/sim_", stochastic_num, "_phcv_", sherlock_table_run$p_trans_hcv, "_relhcvhiv_", sherlock_table_run$rel_infect, "_releqsx_", sherlock_table_run$sx_rel, "_out_full.csv"), na = "", row.names = F)
    fwrite(deaths, paste0(results_folder, sens, "/sim_", stochastic_num, "_phcv_", sherlock_table_run$p_trans_hcv, "_relhcvhiv_", sherlock_table_run$rel_infect, "_releqsx_", sherlock_table_run$sx_rel, "_deaths.csv"), na = "", row.names = F)
  }
  
  rm(out_full, deaths)
}
}
