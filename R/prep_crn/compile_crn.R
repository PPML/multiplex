message("Running Compile CRN for Scenario: ", sens)

n.cores <- 30

my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)

out <- foreach (stochastic_num = 1:n.cores, .combine = 'rbind',
                .packages = c("data.table")) %dopar% {
    
    if (sens == "baseline") {
      scenario_df <- fread("input/final_scenario_table_full.csv")
    } else {
      scenario_df <- fread("input/final_scenario_table_main.csv")
    }
        
      for (s_id in 1:max(scenario_df$scenario_id)) {
        
      if (sens != "trans_calib" | (sens == "trans_calib" & s_id == 1)) {
          
        scenario_df_filtered <- scenario_df[scenario_id==s_id]  
        for (col_name in names(scenario_df_filtered)) {
          assign(col_name, scenario_df_filtered[[col_name]])
        }
        rm(col_name)
      
      if (sens %in% c("trans_prob_low", "trans_prob_high", "hcv_rr_tx_on")) {
        load(paste0(crn_folder, "baseline/fixed/", stochastic_num, ".RData"))
        load(paste0(crn_folder, sens, "/trans/", stochastic_num, ".RData"))
        load(paste0(crn_folder, "baseline/ssp/", stochastic_num, "_ssp_", p_ssp, ".RData"))
        load(paste0(crn_folder, "baseline/hcv/", stochastic_num, "_test_", round(p_test_hcv, 4), "_treat_", round(p_treat_hcv, 4), ".RData"))
        load(paste0(crn_folder, "baseline/hiv/", stochastic_num, "_ltfu_", round(p_ltfu_hiv, 4), ".RData"))
        load(paste0(crn_folder, "baseline/cess/", stochastic_num, "_cess_", round(cessation_rate, 4), "_scaledrm_", scale_drug_mort, ".RData"))
      } else if (sens %in% c("test_coverage_low", "test_coverage_high")) {
        load(paste0(crn_folder, sens, "/fixed/", stochastic_num, ".RData"))
        load(paste0(crn_folder, "baseline/trans/", stochastic_num, ".RData"))
        load(paste0(crn_folder, "baseline/ssp/", stochastic_num, "_ssp_", p_ssp, ".RData"))
        load(paste0(crn_folder, sens, "/hcv/", stochastic_num, "_test_", round(p_test_hcv, 4), "_treat_", round(p_treat_hcv, 4), ".RData"))
        load(paste0(crn_folder, sens, "/hiv/", stochastic_num, "_ltfu_", round(p_ltfu_hiv, 4), ".RData"))
        load(paste0(crn_folder, "baseline/cess/", stochastic_num, "_cess_", round(cessation_rate, 4), "_scaledrm_", scale_drug_mort, ".RData"))
      } else if (sens %in% c("baseline", "init_high_25", "init_low_25")) {
        load(paste0(crn_folder, sens, "/fixed/", stochastic_num, ".RData"))
        load(paste0(crn_folder, sens, "/trans/", stochastic_num, ".RData"))
        load(paste0(crn_folder, sens, "/ssp/", stochastic_num, "_ssp_", p_ssp, ".RData"))
        load(paste0(crn_folder, sens, "/hcv/", stochastic_num, "_test_", round(p_test_hcv, 4), "_treat_", round(p_treat_hcv, 4), ".RData"))
        load(paste0(crn_folder, sens, "/hiv/", stochastic_num, "_ltfu_", round(p_ltfu_hiv, 4), ".RData"))
        load(paste0(crn_folder, sens, "/cess/", stochastic_num, "_cess_", round(cessation_rate, 4), "_scaledrm_", scale_drug_mort, ".RData"))
      } else if (sens == "hiv_tasp_low") {
        load(paste0(crn_folder,"/baseline/fixed/", stochastic_num, ".RData"))
        load(paste0(crn_folder,"/baseline/trans/", stochastic_num, ".RData"))
        load(paste0(crn_folder,"/baseline/ssp/", stochastic_num, "_ssp_", p_ssp, ".RData"))
        load(paste0(crn_folder,"/baseline/hcv/", stochastic_num, "_test_", round(p_test_hcv, 4), "_treat_", round(p_treat_hcv, 4), ".RData"))
        load(paste0(crn_folder,"/baseline/hiv/", stochastic_num, "_ltfu_", round(p_ltfu_hiv, 4), ".RData"))
        load(paste0(crn_folder,"/baseline/cess/", stochastic_num, "_cess_", round(cessation_rate, 4), "_scaledrm_", scale_drug_mort, ".RData"))
      } else if (sens == "hcv_tx_complete") {
        load(paste0(crn_folder,"/baseline/fixed/", stochastic_num, ".RData"))
        load(paste0(crn_folder, "baseline/trans/", stochastic_num, ".RData"))
        load(paste0(crn_folder,"/baseline/ssp/", stochastic_num, "_ssp_", p_ssp, ".RData"))
        load(paste0(crn_folder, sens, "/hcv/", stochastic_num, "_test_", round(p_test_hcv, 4), "_treat_", round(p_treat_hcv, 4), ".RData"))
        load(paste0(crn_folder,"/baseline/hiv/", stochastic_num, "_ltfu_", round(p_ltfu_hiv, 4), ".RData"))
        load(paste0(crn_folder,"/baseline/cess/", stochastic_num, "_cess_", round(cessation_rate, 4), "_scaledrm_", scale_drug_mort, ".RData"))
      } else if (sens == "trans_calib") {
        load(paste0(crn_folder,"/baseline/fixed/", stochastic_num, ".RData"))
        load(paste0(crn_folder,"/baseline/ssp/", stochastic_num, "_ssp_", p_ssp, ".RData"))
        load(paste0(crn_folder,"/baseline/hcv/", stochastic_num, "_test_", round(p_test_hcv, 4), "_treat_", round(p_treat_hcv, 4), ".RData"))
        load(paste0(crn_folder,"/baseline/hiv/", stochastic_num, "_ltfu_", round(p_ltfu_hiv, 4), ".RData"))
        load(paste0(crn_folder,"/baseline/cess/", stochastic_num, "_cess_", round(cessation_rate, 4), "_scaledrm_", scale_drug_mort, ".RData"))
      }
      
      if (scale_drug_mort == 1 & sens != "trans_calib") {
        save(cess_df, drm_df, ssp_df, 
             seed_hcv_df, hcv_test_df, hcv_treat_df, base_hcv_treat_df,
             hcv_trans_df, hiv_trans_df,
             seed_hiv_df, hiv_mort_df_baseline, cd4_art_df_baseline,
             hiv_mort_df_seed, cd4_art_df_seed, hiv_treat_df, hiv_mort_mat, 
             file=paste0(crn_folder, sens, "/compiled_main/variable_", stochastic_num, "_scenario_", s_id, ".RData"))
      } else if (sens == "trans_calib") {
        save(cess_df, drm_df, ssp_df, 
             seed_hcv_df, hcv_test_df, hcv_treat_df, base_hcv_treat_df,
             seed_hiv_df, hiv_mort_df_baseline, cd4_art_df_baseline,
             hiv_mort_df_seed, cd4_art_df_seed, hiv_treat_df, hiv_mort_mat, 
             file=paste0(crn_folder, "baseline/compiled_trans_calib/variable_", stochastic_num, "_scenario_", s_id, ".RData"))
      }  else {
        save(cess_df, drm_df, ssp_df, 
             seed_hcv_df, hcv_test_df, hcv_treat_df, base_hcv_treat_df,
             hcv_trans_df, hiv_trans_df,
             seed_hiv_df, hiv_mort_df_baseline, cd4_art_df_baseline,
             hiv_mort_df_seed, cd4_art_df_seed, hiv_treat_df, hiv_mort_mat, 
             file=paste0(crn_folder, sens, "/compiled_main_scaledrm_", scale_drm, "/variable_", stochastic_num, "_scenario_", s_id, ".RData"))
      }
    }
  }
}

parallel::stopCluster(cl = my.cluster)

