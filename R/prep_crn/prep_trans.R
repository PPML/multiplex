message("Running Transmission CRN for Scenario: ", sens)

n.cores <- 30

my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)

out <- foreach (stochastic_num = 1:n.cores, .combine = 'rbind', .packages = c("data.table")) %dopar% {

  p_trans_hcv <- 0.006
  p_trans_hiv_eq <- p_trans_hcv*((1-.67)/(1-.35))/(1.3)
  p_trans_hiv_sx <- p_trans_hiv_eq/2.5
  p_trans_hiv_genpop <- 1-((1-(p_trans_hiv_sx*(1.2/258.3)*.87*.65))^(mean_deg_sx*(1-p_sx_pwid)))
  
    if (sens == "hcv_rr_tx_on") {
      reinfect_scale <- 0.34
    } else {
      reinfect_scale <- 1
    }
  
    if (sens == "trans_prob_high") {
      p_trans_hiv_eq <- p_trans_hiv_eq * 1.25
      p_trans_hiv_sx <- p_trans_hiv_sx * 1.25
      p_trans_hiv_genpop <- p_trans_hiv_genpop * 1.25
      p_trans_hcv <- p_trans_hcv * 1.25
    }
    if (sens == "trans_prob_low") {
      p_trans_hiv_eq <- p_trans_hiv_eq * 0.75
      p_trans_hiv_sx <- p_trans_hiv_sx * 0.75
      p_trans_hiv_genpop <- p_trans_hiv_genpop * 0.75
      p_trans_hcv <- p_trans_hcv * 0.75
    }
      
  if (sens %in% c("init_low_25", "init_high_25")) {
    load(paste0(crn_folder, sens, "/fixed/", stochastic_num, ".RData"))
  } else {
    load(paste0(crn_folder, "baseline/fixed/", stochastic_num, ".RData"))
  }     
  
      ## Transmission HIV----
      set.seed(14151617+stochastic_num)
      hiv_sx_trans_mat <- matrix(runif((max_time*20*nrow(df)), 0, 1), nrow = nrow(df))
      hiv_eq_trans_mat <- matrix(runif((max_time*20*nrow(df)), 0, 1), nrow = nrow(df))
      hiv_genpop_trans_mat <- matrix(runif((max_time*nrow(df)), 0, 1), nrow = nrow(df))
      
      hiv_trans_df <- copy(df)
      hiv_trans_df <- hiv_trans_df[,.(pid)]
      hiv_trans_df <- hiv_trans_df[, hiv_sx_exposures:=apply(hiv_sx_trans_mat, 1, function(x) which(x < (p_trans_hiv_sx))[1])]
      hiv_trans_df <- hiv_trans_df[, hiv_eq_exposures:=apply(hiv_eq_trans_mat, 1, function(x) which(x < (p_trans_hiv_eq))[1])]
      hiv_trans_df <- hiv_trans_df[, hiv_genpop_exposures:=apply(hiv_genpop_trans_mat, 1, function(x) which(x < (p_trans_hiv_genpop))[1])]
      rm(hiv_sx_trans_mat, hiv_eq_trans_mat, hiv_genpop_trans_mat)
      
      ## Transmission HCV----
      set.seed(15161718+stochastic_num)
      hcv_eq_trans_mat <- matrix(runif((max_time*10*nrow(df)), 0, 1), nrow = nrow(df))
      set.seed((15161718*2)+stochastic_num)
      hcv_eq_trans_mat_tx <- matrix(runif((max_time*10*nrow(df)), 0, 1), nrow = nrow(df))
      
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
      
      save(hcv_trans_df, hiv_trans_df, file=paste0(crn_folder, sens, "/trans/", stochastic_num, ".RData"))
    }

parallel::stopCluster(cl = my.cluster)
