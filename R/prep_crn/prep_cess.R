message("Running Cessation CRN for Scenario: ", sens)

n.cores <- 5

my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)

out <- foreach (s_id = 1:n.cores, .combine = 'rbind', .packages = c("data.table")) %dopar% {
  
  cessation_rate_base <- 0.013
  cess_levels <- c(cessation_rate_base, ((cessation_rate_base*12)+.15)/12, ((cessation_rate_base*12)+.3)/12)

    for (stochastic_num in (s_id*6-5):(s_id*6)) {
      
      for (cessation_rate in cess_levels) {

        load(paste0(crn_folder, sens, "/fixed/", stochastic_num, ".RData"))
        
        ## Cessation Status ----
        message("Running Cessation and Drug-Related Mortality CRN")
        set.seed(34567+stochastic_num)
        cess_mat <- matrix(runif((10*max_time*nrow(df)), 0, 1), nrow = nrow(df))
        set.seed((34567*2)+stochastic_num)
        cess_mat_baseline <- matrix(runif((10*max_time*nrow(df)), 0, 1), nrow = nrow(df))
        set.seed(45678+stochastic_num)
        perm_cess_mat <- matrix(runif((n_cess_allowed*nrow(df)), 0, 1), nrow = nrow(df))
        
        cess_df <- copy(df)
        cess_df <- cess_df[, .(pid, current_pwid, time_enter)]
        cess_df <- cess_df[, perm_quit_num:=apply(perm_cess_mat, 1, function(x) which(x < permanent_quit_prob)[1])]
        
        for (i in 1:n_cess_allowed) {
          cess_df <- cess_df[, paste0("cess_time_", i):=apply(cess_mat, 1, function(x) which(x < cessation_rate)[i])]
        }
        ## Interval between cessation events (time counts for active use)
        for (i in n_cess_allowed:2) {
          cess_df <- cess_df[, paste0("cess_time_", i):=get(paste0("cess_time_", i))-get(paste0("cess_time_", (i-1)))]
        }
        
        ## BASELINE CESSATION ##
        for (i in 1:n_cess_allowed) {
          cess_df <- cess_df[, paste0("base_cess_time_", i):=apply(cess_mat_baseline, 1, function(x) which(x < cessation_rate_base)[i])]
        }
        ## Interval between cessation events (time counts for active use)
        for (i in n_cess_allowed:2) {
          cess_df <- cess_df[, paste0("base_cess_time_", i):=get(paste0("base_cess_time_", i))-get(paste0("base_cess_time_", (i-1)))]
        }
        
        set.seed(56789+stochastic_num)
        ## Duration of cessation before relapse
        for (i in 0:n_cess_allowed) {
          cess_df <-cess_df[, paste0("quit_time_", i):=round(rexp(.N, rate = relapse_rate))]
        }
        
        cess_df <- cess_df[, cess_time_true_0:=time_enter]
        cess_df <- cess_df[current_pwid==1, quit_time_0:=0]
        cess_df <- cess_df[, post_intervention_counter:=1]
        ## Real model time
        for (i in 1:n_cess_allowed) {
          ## Intervention Period
          cess_df <- cess_df[get(paste0("cess_time_true_", (i-1))) <= (intervention_stop), 
                             paste0("cess_time_true_", i):=get(paste0("cess_time_true_", (i-1)))+get(paste0("quit_time_", (i-1)))+get(paste0("cess_time_", (i)))]
          ## Follow-Up (No-Intervention) Period
          for (pic in 1:n_cess_allowed) {
            cess_df <- cess_df[get(paste0("cess_time_true_", (i-1))) > (intervention_stop) & post_intervention_counter==pic, 
                               paste0("cess_time_true_", i):=get(paste0("cess_time_true_", (i-1)))+get(paste0("quit_time_", (i-1)))+get(paste0("base_cess_time_", (pic)))]
          }
          cess_df <- cess_df[get(paste0("cess_time_true_", (i-1))) > (intervention_stop), post_intervention_counter:=post_intervention_counter+1]
        }
        
        cess_df <- cess_df[, c("pid", "perm_quit_num", names(cess_df)[names(cess_df)%like%"true"], names(cess_df)[names(cess_df)%like%"quit_time"]), with = FALSE]
        cess_df <- melt(cess_df, id.vars = c("pid", "perm_quit_num"))
        cess_df <- cess_df[, variable:=gsub("_true", "", variable)]
        cess_df <- cess_df[, c("var1", "var2", "time_block_number"):=tstrsplit(variable, "_")]
        cess_df <- cess_df[, variable:=paste0(var1, "_", var2)]
        cess_df <- cess_df[, c("var1", "var2"):=NULL]
        cess_df <- dcast(cess_df, pid + perm_quit_num + time_block_number ~ variable, value.var = "value")
        cess_df <- cess_df[, time_block_number:=as.numeric(time_block_number)]
        
        rm(list = c("cess_mat", "perm_cess_mat", "cess_mat_baseline"))
        
        ## Expand for each timestep ----
        cess_df_full <- as.data.table(expand.grid(pid = unique(df$pid), time = 1:max_time)) 
        cess_df_full <- merge(cess_df_full, cess_df, by = "pid", allow.cartesian = TRUE)
        cess_df_full <- cess_df_full[time>=cess_time & time<cess_time+quit_time, current_pwid:=0]
        cess_df_full <- cess_df_full[, current_pwid:=mean(current_pwid, na.rm=T), by = c("pid", "time")]
        cess_df_full <- cess_df_full[is.na(current_pwid), current_pwid:=1]
        cess_df_full <- cess_df_full[time >= cess_time & perm_quit_num == time_block_number, permanent_quit:=1]
        cess_df_full <- cess_df_full[, permanent_quit:=mean(permanent_quit, na.rm=T), by = c("pid", "time")]
        cess_df_full <- cess_df_full[is.na(permanent_quit), permanent_quit:=0]
        cess_df_full <- cess_df_full[permanent_quit==1, current_pwid:=0]
        cess_df_full <- unique(cess_df_full[,.(pid, time, current_pwid, permanent_quit)])
        
        cess_df <- copy(cess_df_full)
        rm(cess_df_full)
        
        ## Drug-Related Mortality ----
        set.seed(678910+stochastic_num)
        drug_mort_mat <- matrix(runif((max_time*nrow(df)), 0, 1), nrow = nrow(df))
        
        drm_df <- copy(df)
        drm_df <- drm_df[, time_drug_death:=as.numeric(NA)]
        drm_df <- drm_df[, age:=as.numeric(age)]
        for (t in 1:max_time) {
          drm_df <- drm_df[t >= time_enter, age:=age+timestep]
          drm_df <- drm_df[t >= time_enter, age_floor:=pmin(floor(age), 100), by = "pid"]
          drm_df <- merge(drm_df[, c("current_pwid", "permanent_quit", "time", "current_excess_prob", "former_excess_prob"):=NULL], cess_df[time == t], by = "pid")
          drm_df <- merge(drm_df, lt[, .(age, current_excess_prob, former_excess_prob)], by.x = "age_floor", by.y = "age", all.x=T)
          if (t <= (intervention_stop)) {
            drm_df <- drm_df[, mort_prob:=ifelse(current_pwid == 1, current_excess_prob*scale_drug_mort, former_excess_prob*scale_drug_mort)]
          } else {
            drm_df <- drm_df[, mort_prob:=ifelse(current_pwid == 1, current_excess_prob*1, former_excess_prob*1)] ## scale_drug_mort == 1 in post-intervention period
          }
          drm_df <- drm_df[, drug_mort:=ifelse(drug_mort_mat[, t] <= mort_prob, 1, 0)]
          drm_df <- drm_df[is.na(time_drug_death) & drug_mort==1 & t >= time_enter, time_drug_death:=t]
        }
        
        drm_df <- drm_df[, .(pid, time_drug_death)]
        
        save(drm_df, cess_df, 
             file=paste0(crn_folder, sens, "/cess/", stochastic_num, "_cess_", round(cessation_rate, 4), "_scaledrm_", scale_drug_mort, ".RData"))
      }
    }
}

parallel::stopCluster(cl = my.cluster)
