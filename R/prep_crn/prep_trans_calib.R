n.cores <- 5

my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)

out <- foreach (stochastic_num = 1:n.cores, .combine = 'rbind', .packages = c("data.table")) %dopar% {

      for (p_trans_hcv in seq(0.004, 0.01, .0005)) {
        for (rel_infect in seq(1, 2.0, .1)) {
          for (sx_rel in seq(2, 5, .5)) {
            
            p_trans_hiv_eq <- p_trans_hcv*((1-.67)/(1-.35))/(rel_infect)
            p_trans_hiv_sx <- p_trans_hiv_eq/sx_rel
            p_trans_hiv_genpop <- 1-((1-(p_trans_hiv_sx*(1.2/258.3)*.87*.65))^(mean_deg_sx*(1-p_sx_pwid))) 
          
            reinfect_scale <- 1
            
            load(paste0(crn_folder, sens, "/fixed/", stochastic_num, ".RData"))
            
            max_time_calib <- 12*10+1
            
            ## Transmission HIV----
            set.seed(14151617+stochastic_num)
            hiv_sx_trans_mat <- matrix(runif((max_time_calib*20*nrow(df)), 0, 1), nrow = nrow(df))
            hiv_eq_trans_mat <- matrix(runif((max_time_calib*20*nrow(df)), 0, 1), nrow = nrow(df))
            hiv_genpop_trans_mat <- matrix(runif((max_time_calib*nrow(df)), 0, 1), nrow = nrow(df))
            
            hiv_trans_df <- copy(df)
            hiv_trans_df <- hiv_trans_df[,.(pid)]
            hiv_trans_df <- hiv_trans_df[, hiv_sx_exposures:=apply(hiv_sx_trans_mat, 1, function(x) which(x < (p_trans_hiv_sx))[1])]
            hiv_trans_df <- hiv_trans_df[, hiv_eq_exposures:=apply(hiv_eq_trans_mat, 1, function(x) which(x < (p_trans_hiv_eq))[1])]
            hiv_trans_df <- hiv_trans_df[, hiv_genpop_exposures:=apply(hiv_genpop_trans_mat, 1, function(x) which(x < (p_trans_hiv_genpop))[1])]
            rm(hiv_sx_trans_mat, hiv_eq_trans_mat, hiv_genpop_trans_mat)
            
            ## Transmission HCV----
            set.seed(15161718+stochastic_num)
            hcv_eq_trans_mat <- matrix(runif((max_time_calib*10*nrow(df)), 0, 1), nrow = nrow(df))
            set.seed((15161718*2)+stochastic_num)
            hcv_eq_trans_mat_tx <- matrix(runif((max_time_calib*10*nrow(df)), 0, 1), nrow = nrow(df))
            
            hcv_trans_df <- copy(df)
            hcv_trans_df <- hcv_trans_df[,.(pid)]
            for (i in 1:20) {
              hcv_trans_df <- hcv_trans_df[, paste0("hcv_trans_", (i)):=apply(hcv_eq_trans_mat, 1, function(x) which(x < (p_trans_hcv))[i])]
              hcv_trans_df <- hcv_trans_df[, paste0("hcv_transriskred_", (i)):=apply(hcv_eq_trans_mat_tx, 1, function(x) which(x < (p_trans_hcv*reinfect_scale))[i])]
            }
            hcv_trans_df <- melt(hcv_trans_df, id.vars = "pid")
            hcv_trans_df <- hcv_trans_df[, ever_svr:=ifelse(variable%like%"riskred", "ever_tx", "never_tx")]
            hcv_trans_df <- hcv_trans_df[, c("v1", "v2", "times_hcv"):=tstrsplit(variable, "_")]
            hcv_trans_df <- hcv_trans_df[, c("v1", "v2"):=NULL]
            hcv_trans_df <- hcv_trans_df[, times_hcv:=as.numeric(times_hcv)]
            setnames(hcv_trans_df, "value", "exposures")
            hcv_trans_df <- hcv_trans_df[, "variable":=NULL]
            hcv_trans_df <- dcast(hcv_trans_df, pid + times_hcv ~ ever_svr, value.var = "exposures")
            rm("hcv_eq_trans_mat", "hcv_eq_trans_mat_tx")
            setnames(hcv_trans_df, c("ever_tx", "never_tx"), c("ever_tx_exposures", "never_tx_exposures"))
            
            save(hcv_trans_df, hiv_trans_df, file=paste0(crn_folder, sens, "/trans_calib/", stochastic_num, "_phcv_", p_trans_hcv, "_relhcvhiv_", rel_infect, "_releqsx_", sx_rel, ".RData"))
            }
        }
      }
    }

parallel::stopCluster(cl = my.cluster)

sherlock_table <- as.data.table(expand.grid(stochastic_num = 1:5,
                                            p_trans_hcv = seq(0.004, 0.01, .0005),
                                            rel_infect = seq(1, 2, .1),
                                            sx_rel = seq(2, 5, .5)))
sherlock_table <- sherlock_table[, sherlock_id_group:=ceiling(seq_len(.N)/10)]
sherlock_table <- sherlock_table[, sherlock_id:=seq_len(.N), by = "sherlock_id_group"]
sherlock_table <- sherlock_table[, edge_iter:=1]
sherlock_table <- sherlock_table[, scenario_id:=1]
max(sherlock_table$sherlock_id_group)
write.csv(sherlock_table, "input/sherlock/trans_calibration.csv", row.names = F, na = "")
