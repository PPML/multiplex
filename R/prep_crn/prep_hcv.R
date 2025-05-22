message("Running HCV CRN for Scenario: ", sens)

n.cores <- 30

my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)

out <- foreach (stochastic_num = 1:+n.cores, .combine = 'rbind', .packages = c("data.table")) %dopar% {
    for (parameter_id in 1:5) {
      
      if (sens == "hcv_tx_complete") {
        load(paste0(crn_folder, "baseline/fixed/", stochastic_num, ".RData"))
      } else {
        load(paste0(crn_folder, sens, "/fixed/", stochastic_num, ".RData"))
      }
      
      fibrosis_df[is.na(fibrosis_df)] <- censor_time
      
      if (sens == "hcv_tx_complete") {
        p_tx_complete <- 0.75
      } else {
        p_tx_complete <- 0.96
      }
      
      p_test_hcv_base <- rate_to_prob(prob_to_rate(0.077/p_ever_test_hcv), time = 12)
      p_treat_hcv_base <- rate_to_prob(prob_to_rate(.23*p_tx_complete*.94))
      
      treat_hcv_levels <- c(p_treat_hcv_base, 
                            rate_to_prob(prob_to_rate((.23+.075)*p_tx_complete*.94)), 
                            rate_to_prob(prob_to_rate((.23+.15)*p_tx_complete*.94)),
                            rate_to_prob(prob_to_rate((.25+.225)*p_tx_complete*.94)),
                            rate_to_prob(prob_to_rate((.23+.3)*p_tx_complete*.94)))
      
      treat_filenames <- c(rate_to_prob(prob_to_rate(.23*.96*.94)), 
                            rate_to_prob(prob_to_rate((.23+.075)*.96*.94)), 
                            rate_to_prob(prob_to_rate((.23+.15)*.96*.94)),
                            rate_to_prob(prob_to_rate((.25+.225)*.96*.94)),
                            rate_to_prob(prob_to_rate((.23+.3)*.96*.94)))
      
      test_hcv_levels <- c(p_test_hcv_base, 
                           rate_to_prob(prob_to_rate(0.077/p_ever_test_hcv)+0.25, time = 12),
                           rate_to_prob(prob_to_rate(0.077/p_ever_test_hcv)+0.75, time = 12),
                           rate_to_prob(prob_to_rate(0.077/p_ever_test_hcv)+1.5, time = 12),
                           rate_to_prob(prob_to_rate(0.077/p_ever_test_hcv)+2.25, time = 12))
            
      p_treat_hcv <- treat_hcv_levels[parameter_id]
      p_treat_filename <- treat_filenames[parameter_id]
      p_test_hcv <- test_hcv_levels[parameter_id]
      
      message("Running HCV Clinical CRN")
      
      ## Seed fibrosis distribution
      set.seed(22222 + stochastic_num)
      hcv_seed_df <- copy(df)
      hcv_seed_df <- hcv_seed_df[times_hcv == 1 & age >= 45 & hcv==1, fibrosis_stage_temp := sample(0:5, 1, 
                                                                                                   prob = c(0.1464913, 0.36409872, 0.21760741, 0.13590128, 0.10872103, 0.02718026)), by = "pid"]
      hcv_seed_df <- hcv_seed_df[times_hcv == 1 & age >= 45 & hcv==0, fibrosis_stage_temp := sample(0:5, 1, 
                                                                                                   prob = c(0.33033041, 0.46931939, 0.13898899, 0.03068061, 0.02454449, 0.00613612)), by = "pid"]
      hcv_seed_df <- hcv_seed_df[times_hcv == 1 & age < 45 & hcv==1, fibrosis_stage_temp := sample(0:5, 1, 
                                                                                                   prob = c(0.39812065, 0.47403048, 0.07590984, 0.02596952, 0.02077562, 0.0051939)), by = "pid"]
      hcv_seed_df <- hcv_seed_df[times_hcv == 1 & age < 45 & hcv==0, fibrosis_stage_temp := sample(0:5, 1, 
                                                                                                   prob = c(0.47131268, 0.49457965, 0.02326696, 0.00542035, 0.00433628, 0.00108407)), by = "pid"]

      hcv_seed_df <- merge(hcv_seed_df, fibrosis_df, by = "pid")
      hcv_seed_df <- hcv_seed_df[fibrosis_stage_temp==0, duration_hcv:=sample(1:time_to_f1, 1), by = "pid"]
      hcv_seed_df <- hcv_seed_df[fibrosis_stage_temp==1, duration_hcv:=sample(time_to_f1:(time_to_f1+time_to_f2), 1), by = "pid"]
      hcv_seed_df <- hcv_seed_df[fibrosis_stage_temp==2, duration_hcv:=sample((time_to_f1+time_to_f2):(time_to_f1+time_to_f2+time_to_f3), 1), by = "pid"]
      hcv_seed_df <- hcv_seed_df[fibrosis_stage_temp==3, duration_hcv:=sample((time_to_f1+time_to_f2+time_to_f3):(time_to_f1+time_to_f2+time_to_f3+time_to_f4), 1), by = "pid"]
      hcv_seed_df <- hcv_seed_df[fibrosis_stage_temp==4, duration_hcv:=sample((time_to_f1+time_to_f2+time_to_f3+time_to_f4):(time_to_f1+time_to_f2+time_to_f3+time_to_f4+time_to_f5), 1), by = "pid"]
      hcv_seed_df <- hcv_seed_df[fibrosis_stage_temp==5, duration_hcv:=(time_to_f1+time_to_f2+time_to_f3+time_to_f4+time_to_f5), by = "pid"]
      hcv_seed_df <- hcv_seed_df[, fibrosis_stage_temp:=NULL]
      
      seed_hcv_df <- hcv_seed_df[,.(pid, duration_hcv)]
      
      ## Clinical: HCV ----
      set.seed(7891011+stochastic_num)
      hcv_test_mat <- matrix(runif((max_time*nrow(df)), 0, 1), nrow = nrow(df))
      hcv_test_df <- melt(data.table(pid = 1:nrow(df), ifelse(hcv_test_mat<=p_test_hcv, 1, 0)), id.vars = "pid", variable.name = "duration_hcv", value.name = "hcv_test_intervention")
      hcv_test_df <- hcv_test_df[, duration_hcv:=as.numeric(gsub("V", "", duration_hcv))]
      
      base_hcv_test_df <- melt(data.table(pid = 1:nrow(df), ifelse(hcv_test_mat<=p_test_hcv_base, 1, 0)), id.vars = "pid", variable.name = "duration_hcv", value.name = "hcv_test_baseline")
      base_hcv_test_df <- base_hcv_test_df[, duration_hcv:=as.numeric(gsub("V", "", duration_hcv))]
      hcv_test_df <- merge(hcv_test_df, base_hcv_test_df, by = c("pid", "duration_hcv"))
      
      set.seed(10111213+stochastic_num)
      hcv_treat_mat <- matrix(runif((n_test_allowed*nrow(df)), 0, 1), nrow = nrow(df))
      hcv_treat_mat <- ifelse(hcv_treat_mat<=p_treat_hcv, 1, 0)
      hcv_treat_df <- melt(data.table(pid = 1:nrow(df), hcv_treat_mat), id.vars = "pid", variable.name = "hcv_test_num_intervention", value.name = "hcv_treat_intervention")
      hcv_treat_df <- hcv_treat_df[, hcv_test_num_intervention:=as.numeric(gsub("V", "", hcv_test_num_intervention))]
      
      base_hcv_treat_mat <- matrix(runif((n_test_allowed*nrow(df)), 0, 1), nrow = nrow(df))
      base_hcv_treat_mat <- ifelse(base_hcv_treat_mat<=p_treat_hcv_base, 1, 0)
      base_hcv_treat_df <- melt(data.table(pid = 1:nrow(df), base_hcv_treat_mat), id.vars = "pid", variable.name = "hcv_test_num_baseline", value.name = "hcv_treat_baseline")
      base_hcv_treat_df <- base_hcv_treat_df[, hcv_test_num_baseline:=as.numeric(gsub("V", "", hcv_test_num_baseline))]
      
      save(seed_hcv_df, hcv_test_df, hcv_treat_df, base_hcv_treat_df, 
             file=paste0(crn_folder, sens, "/hcv/", stochastic_num, "_test_", round(p_test_hcv, 4), "_treat_", round(p_treat_filename, 4), ".RData"))
    }
}

parallel::stopCluster(cl = my.cluster)
