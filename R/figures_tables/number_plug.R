df <- NULL

for (f in list.files("/scratch/users/mreitsma/clearance2_summaries/baseline/", full.names = T)) {
  temp <- fread(paste(f))
  df <- rbind(df, temp, fill = T)
}

df <- df[time<=120]

df <- df[, lapply(.SD, sum, na.rm=T), by = c("scenario_id", "stochastic_num", "edge_iter_num")]
df <- df[, hcv_inc:=new_hcv_incl_acute]
df <- df[, hiv_inc:=n_new_hiv]

mean(df$hcv_inc[df$scenario_id==1])
quantile(df$hcv_inc[df$scenario_id==1], c(.025, .975))
mean(df$hiv_inc[df$scenario_id==1])
quantile(df$hiv_inc[df$scenario_id==1], c(.025, .975))

mean(df$n_hcv[df$scenario_id==1]/df$n_pop_tot[df$scenario_id==1])
quantile(df$n_hcv[df$scenario_id==1]/df$n_pop_tot[df$scenario_id==1], c(.025, .975))
mean(df$n_hiv[df$scenario_id==1]/df$n_pop_tot[df$scenario_id==1])
quantile(df$n_hiv[df$scenario_id==1]/df$n_pop_tot[df$scenario_id==1], c(.025, .975))
mean(df$n_hiv_hcv[df$scenario_id==1]/df$n_pop_tot[df$scenario_id==1])
quantile(df$n_hiv_hcv[df$scenario_id==1]/df$n_pop_tot[df$scenario_id==1], c(.025, .975))

## QALYS
df <- NULL

for (f in list.files("/scratch/users/mreitsma/clearance2_summaries/baseline/", full.names = T)) {
  temp <- fread(paste(f))
  df <- rbind(df, temp, fill = T)
}

disc <- (1+.03)^(1/12)-1

df <- df[, disc_utility:=sum_utility*(1/((1+disc)^(time-1)))/12/1240]
df <- df[, disc_ly:=n_pop_tot*(1/((1+disc)^(time-1)))/12/1240]
df <- df[, raw_utility:=sum_utility/12/1240]
df <- df[, raw_ly:=n_pop_tot/12/1240]
df <- df[, lapply(.SD, sum, na.rm=T), by = c("scenario_id", "stochastic_num", "edge_iter_num")]

mean(df$disc_ly[df$scenario_id==1])
quantile(df$disc_ly[df$scenario_id==1], c(.025, .975))

mean(df$disc_utility[df$scenario_id==1])
quantile(df$disc_utility[df$scenario_id==1], c(.025, .975))

mean(df$raw_ly[df$scenario_id==1])
quantile(df$raw_ly[df$scenario_id==1], c(.025, .975))

mean(df$raw_utility[df$scenario_id==1])
quantile(df$raw_utility[df$scenario_id==1], c(.025, .975))

## Cases averted
df <- NULL

for (f in list.files("/scratch/users/mreitsma/clearance2_summaries/baseline/", full.names = T)) {
  temp <- fread(paste(f))
  df <- rbind(df, temp, fill = T)
}

df <- df[ time <= 120]
df <- df[, lapply(.SD, sum, na.rm=T), by = c("scenario_id", "stochastic_num", "edge_iter_num")]

mean(df$n_new_hiv[df$scenario_id==1]) / 1240 * 3694500 / 1000000
mean(df$new_hcv_incl_acute[df$scenario_id==1]) / 1240 * 3694500 / 1000000

files <- list.files(path = "/scratch/users/mreitsma/clearance2_results/baseline/", pattern = "deaths", full.names = T)
files <- files[files%like%"scenario_1_"]

out <- NULL
for (f in files) {
  temp <- fread(paste0(f))
  out <- rbind(out, temp, fill = T)
}

out <- out[analytic==1]

summary <- out %>%
  group_by(stochastic_num, edge_iter_num) %>%
  summarize(age = mean(age))
mean(summary$age)
quantile(summary$age, c(.025, .975))

summary <- out[, drm:=ifelse(cod=="Drug", 1, 0)] %>%
  group_by(stochastic_num, edge_iter_num) %>%
  summarize(drm = mean(drm))
mean(summary$drm)
quantile(summary$drm, c(.025, .975))

## Intervention Effects

df <- NULL

for (f in list.files("/scratch/users/mreitsma/clearance2_summaries/baseline/", full.names = T)) {
  temp <- fread(paste(f))
  df <- rbind(df, temp, fill = T)
}

df <- df[time<=120]

df <- df[, lapply(.SD, sum, na.rm=T), by = c("scenario_id", "stochastic_num", "edge_iter_num")]
df <- df[, hcv_inc:=new_hcv_incl_acute]
df <- df[, hiv_inc:=n_new_hiv]

df <- df[scenario_id == 1, base_hcv:=hcv_inc]
df <- df[scenario_id == 1, base_hiv:=hiv_inc]
df <- df[, base_hcv:=mean(base_hcv, na.rm=T), by = c("stochastic_num", "edge_iter_num")]
df <- df[, base_hiv:=mean(base_hiv, na.rm=T), by = c("stochastic_num", "edge_iter_num")]
df <- df[base_hcv!=0, pct_hcv:=(hcv_inc - base_hcv)/base_hcv]
df <- df[base_hiv!=0, pct_hiv:=(hiv_inc - base_hiv)/base_hiv]
df <- df[base_hiv==0, pct_hiv:=0]
df <- df[base_hcv==0, pct_hcv:=0]

df <- unique(df[,.(scenario_id, stochastic_num, edge_iter_num, pct_hcv, pct_hiv)])

scenario_df <- fread("input/final_scenario_table_full.csv")

out_summary <- merge(df, scenario_df, by = "scenario_id")

out_summary <- out_summary[ cessation_rate == 0.013, cess_level:="Baseline Cessation"]
out_summary <- out_summary[ cessation_rate > .014 & cessation_rate < .03, cess_level:="Moderate Cessation"]
out_summary <- out_summary[ cessation_rate > .035, cess_level:="High Cessation"]
out_summary <- out_summary[ p_ssp == .53, ssp_level:="Baseline Transmission"]
out_summary <- out_summary[ p_ssp == .68, ssp_level:="Moderate Transmission Reduction"]
out_summary <- out_summary[ p_ssp == 0.83, ssp_level:="High Transmission Reduction"]
out_summary <- out_summary[ p_treat_hcv < 0.21, test_level:="Baseline Test and Treat"]
out_summary <- out_summary[ p_treat_hcv > 0.25 & p_treat_hcv < 0.3, test_level:="*Low Test and Treat"]
out_summary <- out_summary[ p_treat_hcv > 0.3 & p_treat_hcv < 0.4, test_level:="Moderate Test and Treat"]
out_summary <- out_summary[ p_treat_hcv > 0.4 & p_treat_hcv < 0.45, test_level:="*Mod-High Test and Treat"]
out_summary <- out_summary[ p_treat_hcv > 0.45 & p_treat_hcv < 0.5, test_level:="High Test and Treat"]

unique(out_summary$scenario_id[out_summary$test_level=="High Test and Treat" & out_summary$ssp_level=="High Transmission Reduction" &
                          out_summary$cess_level=="High Cessation"])
unique(out_summary$scenario_id[out_summary$test_level=="Moderate Test and Treat" & out_summary$ssp_level=="Moderate Transmission Reduction" &
                                 out_summary$cess_level=="Moderate Cessation"])
temp <- copy(out_summary)
temp <- temp[, lower_hcv:=quantile(pct_hcv, .025), by = "scenario_id"]
temp <- temp[, lower_hiv:=quantile(pct_hiv, .025), by = "scenario_id"]
temp <- temp[, upper_hcv:=quantile(pct_hcv, .975), by = "scenario_id"]
temp <- temp[, upper_hiv:=quantile(pct_hiv, .975), by = "scenario_id"]
temp <- temp[, pct_hcv:=mean(pct_hcv), by = "scenario_id"]
temp <- temp[, pct_hiv:=mean(pct_hiv), by = "scenario_id"]
temp <- unique(temp[,.(pct_hiv, pct_hcv, lower_hcv, lower_hiv, upper_hcv, upper_hiv, cess_level, ssp_level, test_level)])

temp[cess_level == "High Cessation" & ssp_level=="High Transmission Reduction" & test_level == "High Test and Treat"][, lapply(.SD, FUN = function(x) round(x*100)), by = c("cess_level", "ssp_level", "test_level")]
temp[cess_level == "Moderate Cessation" & ssp_level=="Moderate Transmission Reduction" & test_level == "Moderate Test and Treat"][, lapply(.SD, FUN = function(x) round(x*100)), by = c("cess_level", "ssp_level", "test_level")]
