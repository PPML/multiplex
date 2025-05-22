df <- NULL

for (f in list.files("/scratch/users/mreitsma/clearance2_summaries/baseline/", full.names = T)) {
  temp <- fread(paste(f))
  df <- rbind(df, temp, fill = T)
}

df <- df[time<=120]

for (i in c(1, 10, 19, 28, 37)) {
  print(rate_to_prob(prob_to_rate(mean(df$n_cure_hcv[df$scenario_id==i]/(df$n_hcv[df$scenario_id==i] + df$n_cure_hcv[df$scenario_id==i]), na.rm=T)), time = 1/12))
}

for (i in c(1, 10, 19, 28, 37)) {
  print(mean(df$n_treat_hiv[df$scenario_id==i]/(df$n_hiv[df$scenario_id==i]), na.rm=T))
}

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
out_summary <- out_summary[ p_ssp == .53, ssp_level:="Baseline Prevention Interventions"]
out_summary <- out_summary[ p_ssp == .68, ssp_level:="Moderate Prevention Interventions"]
out_summary <- out_summary[ p_ssp == 0.83, ssp_level:="High Prevention Interventions"]
out_summary <- out_summary[ p_treat_hcv < 0.21, test_level:="Baseline Test and Treat"]
out_summary <- out_summary[ p_treat_hcv > 0.25 & p_treat_hcv < 0.3, test_level:="*Low Test and Treat"]
out_summary <- out_summary[ p_treat_hcv > 0.3 & p_treat_hcv < 0.4, test_level:="Moderate Test and Treat"]
out_summary <- out_summary[ p_treat_hcv > 0.4 & p_treat_hcv < 0.45, test_level:="*Mod-High Test and Treat"]
out_summary <- out_summary[ p_treat_hcv > 0.45 & p_treat_hcv < 0.5, test_level:="High Test and Treat"]

temp <- copy(out_summary)
temp <- temp[, pct_hcv:=mean(pct_hcv), by = "scenario_id"]
temp <- temp[, pct_hiv:=mean(pct_hiv), by = "scenario_id"]
temp <- temp[, lower_hcv:=quantile(pct_hcv, .025), by = "scenario_id"]
temp <- temp[, lower_hiv:=quantile(pct_hiv, .025), by = "scenario_id"]
temp <- temp[, upper_hcv:=quantile(pct_hcv, .975), by = "scenario_id"]
temp <- temp[, upper_hiv:=quantile(pct_hiv, .975), by = "scenario_id"]

temp <- unique(temp[,.(pct_hiv, pct_hcv, lower_hcv, lower_hiv, upper_hcv, upper_hiv, cess_level, ssp_level, test_level)])
temp <- temp[, y_var := paste0(cess_level, " + ", ssp_level)]
temp <- temp[, test_level:=factor(test_level, levels = c("Baseline Test and Treat", "*Low Test and Treat",
                                                         "Moderate Test and Treat", "*Mod-High Test and Treat", "High Test and Treat"))]
temp <- unique(temp[,.(y_var, test_level, pct_hiv, pct_hcv)])
temp <- melt(temp, id.vars = c("y_var", "test_level"))
temp <- temp[, y_var:=factor(y_var, levels = c("Baseline Cessation + Baseline Prevention Interventions",
                                               "Baseline Cessation + Moderate Prevention Interventions",
                                               "Baseline Cessation + High Prevention Interventions",
                                               "Moderate Cessation + Baseline Prevention Interventions",
                                               "Moderate Cessation + Moderate Prevention Interventions",
                                               "Moderate Cessation + High Prevention Interventions",
                                               "High Cessation + Baseline Prevention Interventions",
                                               "High Cessation + Moderate Prevention Interventions",
                                               "High Cessation + High Prevention Interventions"))]

temp <- temp[variable == "pct_hiv", variable:="HIV"]
temp <- temp[variable == "pct_hcv", variable:="HCV"]

ggplot(data = temp[!is.na(test_level)], aes(x = test_level, y = y_var, fill = value)) +
  geom_tile() +
  geom_text(aes(label = scales::percent(round(value, digits = 2), accuracy = 1)), color = "black", fontface = "bold") +
  scale_fill_distiller(palette = "Spectral", labels = scales::percent) +
  facet_wrap(~variable) +
  theme_bw() +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  labs(x = "", y = "", fill = "Percent Reduction in\nIncidence vs. Baseline") +
  theme(legend.position = "bottom", legend.key.width = unit(2, "cm"),
        legend.key.height = unit(1, "cm"),
        text = element_text(size = 14),
        axis.text.x = element_text(angle = 60, hjust = 1),
        strip.background = element_rect(fill = "NA"),
        strip.text = element_text(size = 18))
ggsave("R/figures_tables/Figure_Heatmap_Incidence_Full_Bold.pdf", width = 12, height = 7)

ggplot(data = temp[!is.na(test_level) & test_level!="*Low Test and Treat" & test_level!="*Mod-High Test and Treat"], aes(x = test_level, y = y_var, fill = value)) +
  geom_tile() +
  geom_text(aes(label = scales::percent(round(value, digits = 2), accuracy = 1)), color = "black", fontface = "bold") +
  # scale_fill_viridis(labels = scales::percent) +
  scale_fill_distiller(palette = "Spectral", labels = scales::percent) +
  facet_wrap(~variable) +
  theme_bw() +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  labs(x = "", y = "", fill = "Percent Reduction in\nIncidence vs. Baseline") +
  theme(legend.position = "bottom", legend.key.width = unit(2, "cm"),
        legend.key.height = unit(1, "cm"),
        text = element_text(size = 14),
        axis.text.x = element_text(angle = 60, hjust = 1),
        strip.background = element_rect(fill = "NA"),
        strip.text = element_text(size = 18))
ggsave("R/figures_tables/Figure_Heatmap_Incidence_ThreeLvl.pdf", width = 12, height = 7)