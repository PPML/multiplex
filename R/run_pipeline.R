library(EpiModel)
library(data.table)
library(readxl)
library(foreach)
library(doParallel)
library(ggplot2)

setwd("/scratch/users/mreitsma/ppml14")

source("R/functions/utilities.R")
source("R/functions/nw_construct.R")
source("R/functions/transmission.R")
source("input/parameters.R")

make_folders <- FALSE

if (isTRUE(make_folders)) {
  base <- ("/scratch/users/mreitsma/clearance2_crn_draws")
  for (i in c("baseline", "test_coverage_low", "test_coverage_high", "trans_prob_low", "trans_prob_high",
              "hiv_tasp_low", "hcv_rr_tx_on", "init_high_25", "init_low_25", "hcv_tx_complete")) {
    for (j in c("ssp", "cess", "trans", "fixed", "hcv", "hiv", "compiled_main"))
      dir.create(file.path(base, i, j), recursive = T)
      if (i == "baseline" & j == "compiled_main") {
      dir.create(file.path(base, i, paste0(j, "_scaledrm_0.5")), recursive = T)
      dir.create(file.path(base, i, paste0(j, "_scaledrm_0.75")), recursive = T)
    }
  }

  base <- ("/scratch/users/mreitsma/clearance2_results")
  for (i in c("baseline", "test_coverage_low", "test_coverage_high", "trans_prob_low", "trans_prob_high",
              "hiv_tasp_low", "hcv_rr_tx_on", "init_high_25", "init_low_25", "hcv_tx_complete")) {
    dir.create(file.path(base, i))
    if (i == "baseline") {
      dir.create(file.path(base, paste0(i, "_scaledrm_0.5")), recursive = T)
      dir.create(file.path(base, paste0(i, "_scaledrm_0.75")), recursive = T)
    }
  }

  base <- ("/scratch/users/mreitsma/clearance2_summaries")
  for (i in c("baseline", "trans_calib", "test_coverage_low", "test_coverage_high", "trans_prob_low", "trans_prob_high",
              "hiv_tasp_low", "hcv_rr_tx_on", "init_high_25", "init_low_25", "hcv_tx_complete")) {
    dir.create(file.path(base, i))
    if (i == "baseline") {
      dir.create(file.path(base, paste0(i, "_scaledrm_0.5")), recursive = T)
      dir.create(file.path(base, paste0(i, "_scaledrm_0.75")), recursive = T)
    }
  }
}

lt <- load_life_tables(timestep = 1) # Lifetable is in monthly timesteps as default

crn_folder <- ("/scratch/users/mreitsma/clearance2_crn_draws/")
results_folder <- ("/scratch/users/mreitsma/clearance2_results/")
summaries_folder <- ("/scratch/users/mreitsma/clearance2_summaries/")

## Prepare CRN Draws for Analysis (assumes running on a machine with >= 30 cores)

for (sens in c("baseline", "test_coverage_high", "test_coverage_low", "init_high_25", "init_low_25", 
               "trans_prob_low", "trans_prob_high", "hcv_rr_tx_on", "hcv_tx_complete")) {
  # FIXED
  source("R/prep_crn/prep_fixed.R")

  # CESSATION
  if (sens %in% c("baseline", "init_high_25", "init_low_25")) {
    for (scale_drug_mort in c(1, .75, .5)) {
  if (sens == "baseline" | scale_drug_mort == 1) {
      source("R/prep_crn/prep_cess.R")
  }
    }
  }

  # SSP
  if (sens %in% c("baseline", "init_high_25", "init_low_25")) {
    source("R/prep_crn/prep_ssp.R")
  }

  # HCV
  if (sens %in% c("baseline", "test_coverage_low", "test_coverage_high", "init_high_25", "init_low_25", "hcv_tx_complete")) {
    source("R/prep_crn/prep_hcv.R")
  }

  # Transmission Calibration
  if (sens == "baseline") {
    source("R/prep_crn/prep_trans_calib.R")
  }
}

## RUN HIV ON CLUSTER
# sbatch R/sbatch/crn_hiv.sbatch

## FIT STERGM 
source("R/prep_crn/fit_network_coefficients.R")

## RUN CALIBRATION
sens <- "trans_calib"
scale_drug_mort <- 1
source("R/prep_crn/compile_crn.R")

# sbatch R/sbatch/calibrate_trans.sbatch

## SELECT TRANSMISSION PROBABILITY
source("R/calibration/select_calibrated_params.R")

## PREP CRN TRANS
for (sens in c("baseline", "hcv_rr_tx_on", "trans_prob_low", "trans_prob_high", "init_high_25", "init_low_25")) {
  source("R/prep_crn/prep_trans.R")
}

## COMPILE CRN DRAWS
for (sens in c("baseline",  "test_coverage_low", "test_coverage_high",
               "trans_prob_low", "trans_prob_high", "hiv_tasp_low", "hcv_rr_tx_on", 
               "init_high_25", "init_low_25", "hcv_tx_complete")) {
  scale_drug_mort <- 1
  source("R/prep_crn/compile_crn.R")
}

for (scale_drm in c(0.5, .75)) {
  sens <- "baseline"
  scale_drug_mort <- scale_drm
  source("R/prep_crn/compile_crn.R")
}

for (sens in c("trans_prob_low", "trans_prob_high", "hiv_tasp_low", "hcv_rr_tx_on", "hcv_tx_complete")) {
  for (i in c(1:30)) {
    file.copy(paste0(crn_folder, "baseline/fixed/", i, ".RData"), 
              paste0(crn_folder, sens, "/fixed/", i, ".RData"))
  }
}

## PREP RUN TABLE
sherlock_df <- fread("input/final_scenario_table_full.csv")
sherlock_df <- sherlock_df[, temp:=1]
sherlock_table <- merge(as.data.table(expand.grid(stochastic_num = 1:30,
                                            edge_iter = 1:5, temp  = 1)), sherlock_df, by = "temp", allow.cartesian = TRUE)[, temp:=NULL]
sherlock_table <- sherlock_table[, sherlock_id_group:=ceiling(seq_len(.N)/15)]
sherlock_table <- sherlock_table[, sherlock_id:=seq_len(.N), by = "sherlock_id_group"]
max(sherlock_table$sherlock_id_group)
write.csv(sherlock_table, "input/sherlock/run_full.csv", row.names = F, na = "")

## PREP RUN TABLE SENS
sherlock_df <- fread("input/final_scenario_table_main.csv")
sherlock_df <- sherlock_df[, temp:=1]
sherlock_table <- merge(as.data.table(expand.grid(stochastic_num = 1:30,
                                                  edge_iter = 1:5, temp  = 1)), sherlock_df, by = "temp", allow.cartesian = TRUE)[, temp:=NULL]
sherlock_table <- sherlock_table[, sherlock_id_group:=ceiling(seq_len(.N)/10)]
sherlock_table <- sherlock_table[, sherlock_id:=seq_len(.N), by = "sherlock_id_group"]
max(sherlock_table$sherlock_id_group)
write.csv(sherlock_table, "input/sherlock/run_main.csv", row.names = F, na = "")

## RUN
# sbatch run_full.sbatch

## RERUN (some will time out due to variation in shared cluster)
# files <- list.files(path = paste0(results_folder, "init_high_25"), full.names = F, pattern = "deaths")
# 
# files_df <- data.table(paths = files)
# files_df <- files_df[, paths_ids:=gsub("sim_", "", paths)]
# files_df <- files_df[, paths_ids:=gsub("scenario_", "", paths_ids)]
# files_df <- files_df[, paths_ids:=gsub("edge_", "", paths_ids)]
# files_df <- files_df[, c("stochastic_num", "scenario_id", "edge_iter", "gar"):=tstrsplit(paths_ids, "_", fixed = TRUE)]
# files_df <- files_df[, run:=1]
# files_df <- files_df[, lapply(.SD, as.numeric), by = c("paths", "paths_ids", "gar")]
# 
# scenario_df <- fread("input/sherlock/run_main.csv")
# 
# sherlock_table <- as.data.table(expand.grid(scenario_id = unique(scenario_df$scenario_id), stochastic_num = 1:30, edge_iter = 1:5))
# sherlock_table <- merge(sherlock_table, files_df[,.(stochastic_num, scenario_id, edge_iter, run)], by = c("stochastic_num", "scenario_id", "edge_iter"), all = T)
# sherlock_table <- sherlock_table[is.na(run)]
# sherlock_table <- sherlock_table[, run:=NULL]
# sherlock_table <- sherlock_table[order(stochastic_num, edge_iter)]
# sherlock_table <- sherlock_table[, sherlock_id_group:=ceiling(seq_len(.N)/1)]
# sherlock_table <- sherlock_table[, sherlock_id:=seq_len(.N), by = "sherlock_id_group"]
# max(sherlock_table$sherlock_id_group)

# write.csv(sherlock_table, "input/sherlock/run_main_scale5.csv", row.names = F, na = "")
# write.csv(sherlock_table, "input/sherlock/run_main_scale75.csv", row.names = F, na = "")
# write.csv(sherlock_table, "input/sherlock/run_main_test_low.csv", row.names = F, na = "")
# write.csv(sherlock_table, "input/sherlock/run_main_hcv_rr.csv", row.names = F, na = "")
# write.csv(sherlock_table, "input/sherlock/run_main_init_high.csv", row.names = F, na = "")

