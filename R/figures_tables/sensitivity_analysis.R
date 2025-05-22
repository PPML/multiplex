df <- NULL
s_id <- "output"
  message(s_id)
  for (f in list.files(paste0("/scratch/users/mreitsma/clearance2_summaries/baseline/"), full.names = T)) {
    temp <- fread(paste(f))
    temp <- temp[, sens:=paste0(s_id)]
    df <- rbind(df, temp, fill = T)
  }

sens_ids <- c("test_coverage_low", "test_coverage_high", "trans_prob_low", "trans_prob_high", "hiv_tasp_low", 
              "hcv_rr_tx_on", "init_high_25", "init_low_25", "hcv_tx_complete")

for (s_id in sens_ids) {
  message(s_id)
  for (f in list.files(paste0("/scratch/users/mreitsma/clearance2_summaries/", s_id, "/"), full.names = T)) {
    temp <- fread(paste(f))
    temp <- temp[, sens:=paste0(s_id)]
    df <- rbind(df, temp, fill = T)
  }
}

scenario_df <- fread("input/final_scenario_table_main.csv")

backup <- copy(df)
df <- copy(backup)

df <- df[scenario_id%in%c(1, 2, 3, 4, 7, 10, 14, 19, 27)]

df <- df[time <= 120]

df <- df[, lapply(.SD, sum, na.rm=T), by = c("scenario_id", "stochastic_num", "edge_iter_num", "sens")]
df <- df[, hcv_inc:=new_hcv_incl_acute]
df <- df[, hiv_inc:=n_new_hiv]

df <- df[scenario_id == 1 & sens=="output", base_hcv_main:=hcv_inc]
df <- df[scenario_id == 1 & sens=="output" , base_hiv_main:=hiv_inc]
df <- df[scenario_id == 1 , base_hcv:=hcv_inc]
df <- df[scenario_id == 1 , base_hiv:=hiv_inc]
df <- df[, base_hcv:=mean(base_hcv, na.rm=T), by = c("stochastic_num", "edge_iter_num", "sens")]
df <- df[, base_hiv:=mean(base_hiv, na.rm=T), by = c("stochastic_num", "edge_iter_num", "sens")]
df <- df[, base_hcv_main:=mean(base_hcv_main, na.rm=T), by = c("stochastic_num", "edge_iter_num")]
df <- df[, base_hiv_main:=mean(base_hiv_main, na.rm=T), by = c("stochastic_num", "edge_iter_num")]
df <- df[, pct_hcv:=(hcv_inc - base_hcv)/base_hcv]
df <- df[, pct_hiv:=(hiv_inc - base_hiv)/base_hiv]
df <- df[scenario_id==1, pct_hcv:=(hcv_inc - base_hcv_main)/base_hcv_main]
df <- df[scenario_id==1, pct_hiv:=(hiv_inc - base_hiv_main)/base_hiv_main]

df <- unique(df[,.(scenario_id, stochastic_num, edge_iter_num, sens, pct_hcv, pct_hiv)])
df <- melt(df, id.vars = c("scenario_id", "stochastic_num", "edge_iter_num", "sens"))

df <- df[, mean:=mean(value, na.rm=T), by = c("scenario_id", "sens", "variable")]
df <- df[, samples:=.N, by = c("scenario_id", "sens", "variable")]
df <- df[, lower:=mean-(1.96*(sd(value)/sqrt(samples))), by = c("scenario_id", "sens", "variable")]
df <- df[, upper:=mean+(1.96*(sd(value)/sqrt(samples))), by = c("scenario_id", "sens", "variable")]

df <- unique(df[,.(scenario_id, sens, variable, mean, lower, upper)])

scenario_df <- fread("input/final_scenario_table_main.csv")[,.(p_ssp, cessation_rate, p_treat_hcv, scenario_id)]

df <- merge(df, scenario_df, by = "scenario_id")
df <- df[p_ssp==.53 & cessation_rate==.013 & p_treat_hcv < 0.21, sens_id:="Baseline"]
df <- df[p_ssp==.68 & cessation_rate==.013 & p_treat_hcv < 0.21, sens_id:="Moderate Prevention Interventions"]
df <- df[p_ssp==.83 & cessation_rate==.013 & p_treat_hcv < 0.21, sens_id:="High Prevention Interventions"]
df <- df[p_ssp==.53 & cessation_rate>.014 & cessation_rate<.03 & p_treat_hcv < 0.21, sens_id:="Moderate Cessation"]
df <- df[p_ssp==.53 & cessation_rate>.035 & p_treat_hcv < 0.21, sens_id:="High Cessation"]
df <- df[p_ssp==.68 & cessation_rate>.014 & cessation_rate<.03 & p_treat_hcv > 0.3 & p_treat_hcv < 0.4, sens_id:="Moderate All"]
df <- df[p_ssp==.83 & cessation_rate>.035 & p_treat_hcv > 0.45 & p_treat_hcv < 0.5, sens_id:="High All"]
df <- df[p_ssp==.53 & cessation_rate==.013 & p_treat_hcv > 0.3 & p_treat_hcv < 0.4, sens_id:="Moderate Test and Treat"]
df <- df[p_ssp==.53 & cessation_rate==.013 & p_treat_hcv > 0.45 & p_treat_hcv < 0.5, sens_id:="High Test and Treat"]

df <- df[, sens_id:=factor(sens_id, levels = c("Baseline", "Moderate Prevention Interventions", "Moderate Cessation", "Moderate Test and Treat",
                                               "Moderate All", "High Prevention Interventions", "High Cessation", "High Test and Treat", "High All"))]
df <- df[variable == "pct_hiv", variable:="HIV"]
df <- df[variable == "pct_hcv", variable:="HCV"]

df <- df[sens == "test_coverage_high", sens_lab:="Ever Test HIV 100%; Ever Test HCV 90%"]
df <- df[sens == "test_coverage_low", sens_lab:="Ever Test HIV 80%; Ever Test HCV 70%"]
df <- df[sens == "hiv_tasp_low", sens_lab:="HIV TasP 75% Effective for Injection Exposure"]
df <- df[sens == "hcv_rr_tx_on", sens_lab:="Includes Risk Reduction following HCV Tx"]
df <- df[sens == "trans_prob_high", sens_lab:="Transmission Probability 25% Higher"]
df <- df[sens == "trans_prob_low", sens_lab:="Transmission Probability 25% Lower"]
df <- df[sens == "init_high_25", sens_lab:="Initiation Rate 25% Higher"]
df <- df[sens == "init_low_25", sens_lab:="Initiation Rate 25% Lower"]
df <- df[sens == "hcv_tx_complete", sens_lab:="HCV Treatment Completion Probability 0.75"]
df <- df[sens == "output", sens_lab:="Baseline"]

df <- df[, sens_lab:=factor(sens_lab, levels = c("Baseline",
                                                 "Ever Test HIV 100%; Ever Test HCV 90%", "Ever Test HIV 80%; Ever Test HCV 70%",
                                                 "Transmission Probability 25% Higher", "Transmission Probability 25% Lower",
                                                 "Initiation Rate 25% Higher", "Initiation Rate 25% Lower",
                                                 "Includes Risk Reduction following HCV Tx", "HCV Treatment Completion Probability 0.75",
                                                 "HIV TasP 75% Effective for Injection Exposure"))]

pdf("R/figures_tables/Figure_Sensitivity_Points_AddInit_Taller.pdf", width = 10, height = 9)
ggplot(data = df[sens_id != "Baseline" & variable!="QALYs"], aes(x = sens_id)) + 
  geom_segment(data = df[sens_id != "Baseline" & sens=="output" & variable !="QALYs"], 
               aes(x = as.numeric(as.factor(sens_id))-.4, xend = as.numeric(as.factor(sens_id))+.4, y = mean), 
               color = "black", size = 1.2) +
  geom_point(data = df[sens_id != "Baseline" & sens!="output" & variable !="QALYs"], 
             aes(x = as.numeric(as.factor(sens_id)), y = mean, color = sens_lab), 
             position = position_dodge(width = .8), alpha = .9, size = 1.5) +
  geom_errorbar(data = df[sens_id != "Baseline" & sens!="output" & variable !="QALYs"], 
             aes(x = as.numeric(as.factor(sens_id)), ymin = lower, ymax = upper, color = sens_lab), 
             position = position_dodge(width = .8), alpha = .9, width = 0, linewidth = 0.8) +
  geom_rect(data = df[sens_id != "Baseline" & sens=="output" & variable !="QALYs"], 
            aes(x = as.numeric(as.factor(sens_id)),
                xmin = as.numeric(as.factor(sens_id))-.4, 
                xmax = as.numeric(as.factor(sens_id))+.4,
                ymin = lower,
                ymax = upper), 
            fill = "black", size = 1.2, alpha = .2) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  facet_wrap(~variable) +
  scale_x_continuous(breaks = c(1:9), labels = levels(df$sens_id)) +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  theme(text = element_text(size = 14)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "none")  +
  labs(x = "", y = "Percent Change in Incidence", color = "Sensitivity Analysis", 
       caption = "Black line is the main analysis estimate. Uncertainty displayed is the standard error of the mean estimate.") +
  scale_color_brewer(palette = "Paired") +
  guides(color=guide_legend(ncol=1,byrow=TRUE))
dev.off()

pdf("R/figures_tables/Figure_Sensitivity_Bars_AddInit.pdf", width = 10, height = 4)
ggplot() + 
  geom_bar(data = df[sens_id == "Baseline" & sens!="output"], 
               aes(y = sens_lab, x = mean, fill = sens_lab), stat = "identity") +
  theme_bw() +
  theme(text = element_text(size = 14)) +
  scale_fill_brewer(palette = "Paired") +
  scale_x_continuous(labels = scales::percent) +
  facet_wrap(~variable) +
  theme(legend.position = "none") +
  guides(fill=guide_legend(nrow=4,byrow=TRUE)) +
  labs(y = "", x = "",
       title = "Percent Change in Baseline Incidence vs. Main Analysis", fill = "Sensitivity Analysis")
dev.off()  

pdf("R/figures_tables/Figure_Sensitivity_Legend_AddInit.pdf", width = 12, height = 6)
ggplot() + 
  geom_bar(data = df[sens_id == "Baseline" & sens!="output"], 
           aes(y = sens_lab, x = mean, fill = sens_lab), stat = "identity") +
  theme_bw() +
  theme(text = element_text(size = 14)) +
  scale_fill_brewer(palette = "Paired") +
  scale_x_continuous(labels = scales::percent) +
  facet_wrap(~variable) +
  theme(legend.position = "bottom") +
  guides(fill=guide_legend(nrow=5,byrow=TRUE)) +
  labs(y = "", x = "",
       title = "Percent Change in Baseline Incidence vs. Main Analysis", fill = "")
dev.off()  


