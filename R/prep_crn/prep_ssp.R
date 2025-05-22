message("Running SSP CRN for Scenario: ", sens)

n.cores <- 30

my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)

base_ssp <- p_ssp_base
base_condom <- p_condom_base
intervention_stop_yr <- intervention_stop/12
ssp_levels <- c(base_ssp, base_ssp + 0.15, base_ssp + 0.3)
condom_levels <- c(base_condom, base_condom + 0.15, base_condom + 0.3)  

out <- foreach (stochastic_num = 1:n.cores, .combine = 'rbind', .packages = c("data.table")) %dopar% {
    
  for (int_level in 1:3) {
    p_ssp <- ssp_levels[int_level]
    p_condom <- condom_levels[int_level]
      
    load(paste0(crn_folder, sens, "/fixed/", stochastic_num, ".RData"))
    
      ## SSP ----
      message("Running SSP CRN")
      
      set.seed(98989898+stochastic_num)
      condom_mat <- matrix(runif(((max_time/12)*nrow(df)), 0, 1), nrow = nrow(df))
      condom_mat <- cbind(ifelse(condom_mat[,1:intervention_stop_yr]<=p_condom, 1, 0),
                       ifelse(condom_mat[,(intervention_stop_yr+1):(max_time/12)]<=base_condom, 1, 0))
      set.seed(1617181920+stochastic_num)
      ssp_mat <- matrix(runif(((max_time/12)*nrow(df)), 0, 1), nrow = nrow(df))
      ssp_mat <- cbind(ifelse(ssp_mat[,1:intervention_stop_yr]<=p_ssp, 1, 0),
                          ifelse(ssp_mat[,(intervention_stop_yr+1):(max_time/12)]<=base_ssp, 1, 0))
      set.seed(1718192021+stochastic_num)
      eq_iso_mat <- matrix(runif(((max_time/12)*nrow(df)), 0, 1), nrow = nrow(df))
      eq_iso_mat <- cbind(ifelse(eq_iso_mat[,1:intervention_stop_yr]<=p_iso_eq/p_ssp, 1, 0),
                          ifelse(eq_iso_mat[,(intervention_stop_yr+1):(max_time/12)]<=p_iso_eq/(base_ssp), 1, 0))
      
      condom_df <- melt(data.table(pid=1:nrow(df), condom_mat), id.vars = "pid", variable.name = "ssp_year", value.name = "condom")
      ssp_df <- melt(data.table(pid=1:nrow(df), ssp_mat), id.vars = "pid", variable.name = "ssp_year", value.name = "ssp")
      eq_iso_df <- melt(data.table(pid=1:nrow(df), eq_iso_mat), id.vars = "pid", variable.name = "ssp_year", value.name = "eq_iso")

      ssp_df <- merge(condom_df, ssp_df, by = c("pid", "ssp_year"))
      ssp_df <- merge(ssp_df, eq_iso_df, by = c("pid", "ssp_year"))
      ssp_df <- ssp_df[, ssp_year:=as.numeric(gsub("V", "", ssp_year))]
      ssp_df <- ssp_df[ssp==0, eq_iso:=0]
      
        save(ssp_df, file=paste0(crn_folder, sens, "/ssp/", stochastic_num, "_ssp_", p_ssp, ".RData"))
    }
  }

parallel::stopCluster(cl = my.cluster)
