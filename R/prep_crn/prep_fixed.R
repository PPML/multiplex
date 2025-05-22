message("Running Fixed CRN for Scenario: ", sens)

parallel::detectCores()
n.cores <- 30
my.cluster <- parallel::makeCluster(
 n.cores,
 type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)

out <- foreach (stochastic_num = 1:n.cores, .combine = 'rbind',
                .packages = c("data.table", "EpiModel", "readxl")) %dopar% {
    
    set.seed(12345+stochastic_num)
    
    if (sens == "test_coverage_high") {
      p_ever_test_hiv <- 1
      p_ever_test_hcv <- 0.9
    } 
    if (sens == "test_coverage_low") {
      p_ever_test_hiv <- .8
      p_ever_test_hcv <- 0.7
    }       
    if (sens == "init_high_25") {
      num_enter <- 5
      n_pop_enter <- num_enter*(intervention_stop)
      n_pop <- n_pop_nw + n_pop_enter
    }
    if (sens == "init_low_25") {
      num_enter <- 3
      n_pop_enter <- num_enter*(intervention_stop)
      n_pop <- n_pop_nw + n_pop_enter
    }
    
    df <- data.table(time_enter = 0,
                     pid=1:n_pop_nw,
                     current_pwid = rbinom(n = n_pop_nw, size = 1, prob = p_current),
                     male = rbinom(n = n_pop_nw, size = 1, prob = pct_male),
                     ever_test_hiv = rbinom(n = n_pop_nw, size = 1, prob = p_ever_test_hiv),
                     analytic = 1)
    df <- df[, age:=sample(18:65, n_pop_nw, replace = TRUE,
                           prob = c(rep(pct_age_18_24/length(18:24), length(18:24)),
                                    rep(pct_age_25_29/length(25:29), length(25:29)),
                                    rep(pct_age_30_39/length(30:39), length(30:39)),
                                    rep(pct_age_40_49/length(40:49), length(40:49)),
                                    rep(pct_age_50plus/length(50:65), length(50:65))))]
    df <- df[ever_test_hiv==1, ever_test_hcv := rbinom(n = .N, size = 1, prob = p_ever_test_hcv/p_ever_test_hiv)]
    df <- df[ever_test_hiv==0, ever_test_hcv := 0]
    init_df <- data.table(time_enter = rep(1:(n_pop_enter/num_enter), times = num_enter),
                          current_pwid = 1,
                          male = rbinom(n = n_pop_enter, size = 1, prob = pct_male),
                          ever_test_hiv = rbinom(n = n_pop_enter, size = 1, prob = p_ever_test_hiv),
                          hcv = 0, hiv = 0, times_hcv = 0)
    init_df <- init_df[order(time_enter)]
    init_df <- init_df[, pid:=(n_pop_nw+1):(n_pop_nw+n_pop_enter)]
    init_df <- init_df[, age:=sample(18:24, n_pop_enter, replace = TRUE)]
    init_df <- init_df[, analytic:=ifelse(time_enter<=120, 1, 0)]
    init_df <- init_df[ever_test_hiv==1, ever_test_hcv := rbinom(n = .N, size = 1, prob = p_ever_test_hcv/p_ever_test_hiv)]
    init_df <- init_df[ever_test_hiv==0, ever_test_hcv := 0]
    
    ## SEED HCV + HIV
    df <- df[age < 25, hcv:=rbinom(.N, 1, .243)]
    df <- df[age >=25 & age <30, hcv:=rbinom(.N, 1, .341)]
    df <- df[age >=30 & age <40, hcv:=rbinom(.N, 1, .420)]
    df <- df[age >=40 & age <50, hcv:=rbinom(.N, 1, .433)]
    df <- df[age >=50, hcv:=rbinom(.N, 1, .495)]
    df <- df[hcv==1, times_hcv:= 1]
    df <- df[hcv==0 & age < 30, times_hcv:=rbinom(.N, 1, 0.127/(0.127+0.446))] ## Based on antibody +
    df <- df[hcv==0 & age >=30, times_hcv:=rbinom(.N, 1, 0.21/(0.21+0.345))] ## Based on antibody +
    df <- df[hcv==1, hiv:=rbinom(.N, 1, (0.04)/(0.401))]
    df <- df[hcv==0, hiv:=rbinom(.N, 1, (0.086-0.04)/(1-0.401))]
    
    df <- rbind(df, init_df)
    
    ## Time to background death ----
    df <- df[, age_temp:=age]
    df <- df[, death:=0]
    
    set.seed(23456 + stochastic_num)
    for (i in 1:(max_time)) {
      df <- df[, age_temp:=age_temp+(timestep)]
      df <- df[death==0, death:=background_mort(age_temp), by = "pid"]
      df <- df[death==1, time_background_death:=i]
      df <- df[death==1, death:=2]
    }
    
    df <- df[, c("age_temp", "death"):=NULL]
    
    ## HCV Progression ----
    set.seed(19202122 + stochastic_num)
    f1_mat <- matrix(runif((max_time*nrow(df)), 0, 1), nrow = nrow(df))
    f2_mat <- matrix(runif((max_time*nrow(df)), 0, 1), nrow = nrow(df))
    f3_mat <- matrix(runif((max_time*nrow(df)), 0, 1), nrow = nrow(df))
    f4_mat <- matrix(runif((max_time*nrow(df)), 0, 1), nrow = nrow(df))
    f5_mat <- matrix(runif((max_time*nrow(df)), 0, 1), nrow = nrow(df))
    
    fibrosis_df <- copy(df)
    fibrosis_df <- fibrosis_df[,.(pid)]
    fibrosis_df <- fibrosis_df[, time_to_f1:=apply(f1_mat, 1, function(x) which(x < (f0_f1))[1])]
    fibrosis_df <- fibrosis_df[, time_to_f2:=apply(f2_mat, 1, function(x) which(x < (f1_f2))[1])]
    fibrosis_df <- fibrosis_df[, time_to_f3:=apply(f3_mat, 1, function(x) which(x < (f2_f3))[1])]
    fibrosis_df <- fibrosis_df[, time_to_f4:=apply(f4_mat, 1, function(x) which(x < (f3_f4))[1])]
    fibrosis_df <- fibrosis_df[, time_to_f5:=apply(f5_mat, 1, function(x) which(x < (f4_decomp))[1])]
    
    rm(f1_mat, f2_mat, f3_mat, f4_mat, f5_mat)
    
    ## HCV Mortality ----
    set.seed(20212223 + stochastic_num)
    hcv_mort_mat <- matrix(runif((max_time*nrow(df)), 0, 1), nrow = nrow(df))
    hcv_mort_df <- as.data.table(hcv_mort_mat)
    hcv_mort_df <- hcv_mort_df[, pid:=1:nrow(df)]
    
    ## HCV Acute ----
    set.seed(21222324 + stochastic_num)
    hcv_acute_mat <- matrix(runif((n_test_allowed*nrow(df)), 0, 1), nrow = nrow(df))
    hcv_acute_mat <- ifelse(hcv_acute_mat<p_acute_clear, 1, 0)
    
    hcv_acute_df <- as.data.table(hcv_acute_mat)
    colnames(hcv_acute_df) <- c(paste0(1:n_test_allowed))
    hcv_acute_df <- hcv_acute_df[, pid:=1:nrow(df)]
    hcv_acute_df <- melt(hcv_acute_df, id.vars = "pid")
    setnames(hcv_acute_df, c("variable", "value"), c("times_hcv", "acute_clear"))
    hcv_acute_df <- hcv_acute_df[, times_hcv:=as.numeric(times_hcv)]
    
    save(df, fibrosis_df, hcv_acute_df, hcv_mort_df, lt, 
         file=paste0(crn_folder, sens, "/fixed/", stochastic_num, ".RData"))
  }

parallel::stopCluster(cl = my.cluster)

