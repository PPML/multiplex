df <- NULL

for (f in list.files("/scratch/users/mreitsma/clearance2_summaries/baseline/", full.names = T)) {
  temp <- fread(paste(f))
  df <- rbind(df, temp, fill = T)
}

df <- df[time <= 120]
df <- df[, lapply(.SD, sum, na.rm=T), by = c("scenario_id", "stochastic_num", "edge_iter_num")]
df <- df[scenario_id == 1]
df <- df[, hcv_inc:=new_hcv_incl_acute]
df <- df[, hiv_inc:=n_new_hiv]
df <- df[order(stochastic_num, edge_iter_num)]

df <- df[,.(stochastic_num, edge_iter_num, hcv_inc, hiv_inc)]
df <- df[, s_id:=seq_len(.N)]

convergence <- NULL
for (j in 1:30) {
  message(j)
  s_order <- sample(1:150, replace = T)
  for (i in 1:150) {
    test <- df[s_id%in%s_order[0:i]]
    test <- test[, lapply(.SD, mean)]
    test <- test[, i := i]
    test <- test[, j := j]
    convergence <- rbind(convergence, test, fill = T)
  }
}

ggplot(data = convergence, aes(x = i, y = hiv_inc, group = j)) + geom_line(alpha = .6) +
  theme_bw() +
  labs(x = "Number of Iterations", y = "Baseline New HIV Infections")

ggplot(data = convergence, aes(x = i, y = hcv_inc, group = j)) + geom_line(alpha = .6) +
  theme_bw() +
  labs(x = "Number of Iterations", y = "Baseline New HCV Infections")

df <- NULL

for (f in list.files("/scratch/users/mreitsma/clearance2_summaries/baseline/", full.names = T)) {
  temp <- fread(paste(f))
  df <- rbind(df, temp, fill = T)
}

df <- df[time <= 120]
df <- df[, lapply(.SD, sum, na.rm=T), by = c("scenario_id", "stochastic_num", "edge_iter_num")]
df <- df[, hcv_inc:=new_hcv_incl_acute]
df <- df[, hiv_inc:=n_new_hiv]
df <- df[scenario_id == 1, base_hcv:=hcv_inc]
df <- df[scenario_id == 1, base_hiv:=hiv_inc]
df <- df[, base_hcv:=mean(base_hcv, na.rm=T), by = c("stochastic_num", "edge_iter_num")]
df <- df[, base_hiv:=mean(base_hiv, na.rm=T), by = c("stochastic_num", "edge_iter_num")]
df <- df[, pct_hcv:=(hcv_inc - base_hcv)/base_hcv]
df <- df[, pct_hiv:=(hiv_inc - base_hiv)/base_hiv]

df <- df[order(scenario_id, stochastic_num, edge_iter_num)]

df <- df[,.(stochastic_num, scenario_id, pct_hcv, pct_hiv)]
df <- df[, s_id:=seq_len(.N), by = "scenario_id"]

convergence <- NULL
for (j in 1:150) {
  message(j)
  s_order <- sample(1:150, replace = T)
  for (i in 1:150) {
    test <- df[s_id%in%s_order[0:i]]
    test <- test[, lapply(.SD, mean), by = "scenario_id"]
    test <- test[, i := i]
    test <- test[, j := j]
    convergence <- rbind(convergence, test, fill = T)
  }
}

scenario_df <- fread("input/final_scenario_table_main.csv")

convergence <- merge(convergence, scenario_df, by = "scenario_id")

convergence <- convergence[ cessation_rate == 0.013, cess_level:="Base"]
convergence <- convergence[ cessation_rate > .014 & cessation_rate < .03, cess_level:="Mid"]
convergence <- convergence[ cessation_rate > .035, cess_level:="High"]
convergence <- convergence[ p_ssp == .53, ssp_level:="Base"]
convergence <- convergence[ p_ssp == .68, ssp_level:="Mid"]
convergence <- convergence[ p_ssp == 0.83, ssp_level:="High"]
convergence <- convergence[ p_treat_hcv < 0.21, test_level:="Base"]
convergence <- convergence[ p_treat_hcv > 0.3 & p_treat_hcv < 0.4, test_level:="Mid"]
convergence <- convergence[ p_treat_hcv > 0.45 & p_treat_hcv < 0.5, test_level:="High"]

convergence <- convergence[test_level != ""]

convergence <- convergence[, fac_lab:=paste0(cess_level, " / ", ssp_level, " / ", test_level)]

fac_lab_levels <- unique(convergence[order(cessation_rate, p_ssp, p_treat_hcv)]$fac_lab)
convergence <- convergence[, fac_lab:=factor(fac_lab, levels = fac_lab_levels, ordered = TRUE)]

ggplot(data = convergence, aes(x = i, y = pct_hiv, group = j)) + geom_line(alpha = .4) +
  facet_wrap(~fac_lab) +
  theme_bw() +
  labs(x = "Number of Iterations", y = "Percent Change in HIV Incidence") +
  scale_y_continuous(labels = scales::percent)
ggsave("R/figures_tables/convergence_hiv.pdf", width = 8, height = 7)

ggplot(data = convergence, aes(x = i, y = pct_hcv, group = j)) + geom_line(alpha = .4) +
  facet_wrap(~fac_lab) +
  theme_bw() +
  labs(x = "Number of Iterations", y = "Percent Change in HCV Incidence")  +
  scale_y_continuous(labels = scales::percent)
ggsave("R/figures_tables/convergence_hcv.pdf", width = 8, height = 7)
