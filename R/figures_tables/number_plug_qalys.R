library(data.table)
library(ggplot2)
library(forcats)
library(ggupset)

df <- NULL

for (f in list.files("/scratch/users/mreitsma/clearance2_summaries/baseline/", full.names = T)) {
  temp <- fread(paste(f))
  temp <- temp[, scale_drm:="Baseline Excess Mortality"]
  df <- rbind(df, temp, fill = T)
}

for (f in list.files("/scratch/users/mreitsma/clearance2_summaries/scaledrm_0.75/", full.names = T)) {
  temp <- fread(paste(f))
  temp <- temp[, scale_drm:="25% Reduction in Excess Mortality"]
  df <- rbind(df, temp, fill = T)
}

for (f in list.files("/scratch/users/mreitsma/clearance2_summaries/scaledrm_0.5/", full.names = T)) {
  temp <- fread(paste(f))
  temp <- temp[, scale_drm:="50% Reduction in Excess Mortality"]
  df <- rbind(df, temp, fill = T)
}

disc <- (1+.03)^(1/12)-1

df <- df[, disc_utility:=sum_utility*(1/((1+disc)^(time-1)))/12/1240]
df <- df[, disc_ly:=n_pop_tot*(1/((1+disc)^(time-1)))/12/1240]
df <- df[, raw_utility:=sum_utility/12/1240]
df <- df[, raw_ly:=n_pop_tot/12/1240]

df <- df[, lapply(.SD, sum, na.rm=T), by = c("scenario_id", "stochastic_num", "edge_iter_num", "scale_drm")]

df <- unique(df[,.(scenario_id, stochastic_num, edge_iter_num, scale_drm, disc_utility, disc_ly, raw_utility, raw_ly)])

scenario_df <- fread("input/final_scenario_table_main.csv")

out_summary <- merge(df, scenario_df, by = "scenario_id")

out_summary <- out_summary[ cessation_rate == 0.013, cess_level:="Baseline Cessation"]
out_summary <- out_summary[ cessation_rate > .014 & cessation_rate < .03, cess_level:="Moderate Cessation"]
out_summary <- out_summary[ cessation_rate > .035, cess_level:="High Cessation"]
out_summary <- out_summary[ p_ssp == .53, ssp_level:="Baseline Prevention Interventions"]
out_summary <- out_summary[ p_ssp == .68, ssp_level:="Moderate Prevention Interventions"]
out_summary <- out_summary[ p_ssp == 0.83, ssp_level:="High Prevention Interventions"]
out_summary <- out_summary[ p_treat_hcv < 0.21, test_level:="Baseline Test and Treat"]
out_summary <- out_summary[ p_treat_hcv > 0.3 & p_treat_hcv < 0.4, test_level:="Moderate Test and Treat"]
out_summary <- out_summary[ p_treat_hcv > 0.45 & p_treat_hcv < 0.5, test_level:="High Test and Treat"]

out_summary <- out_summary[!is.na(test_level)]

temp <- copy(out_summary)
temp <- temp[scenario_id==1 & scale_drm == "Baseline Excess Mortality", base:=disc_utility]
temp <- temp[, base:=mean(base, na.rm=T), by = c("stochastic_num", "edge_iter_num")]
temp <- temp[, mean:=mean(disc_utility-base), by = c("scenario_id", "scale_drm")]
temp <- temp[, lower:=quantile(disc_utility-base, .025), by = c("scenario_id", "scale_drm")]
temp <- temp[, upper:=quantile(disc_utility-base, .975), by = c("scenario_id", "scale_drm")]

temp <- unique(temp[,.(scale_drm, mean, lower, upper, cess_level, ssp_level, test_level)])

temp <- temp[, scenario_level:= paste(cess_level, ssp_level, test_level, sep = "-")]
temp <- temp[, scenario_level:=fct_reorder(scenario_level, mean)]

temp <- temp[, scale_drm:=factor(scale_drm, levels = c("Baseline Excess Mortality", "25% Reduction in Excess Mortality", "50% Reduction in Excess Mortality"))]

plot <- ggplot(data = temp[test_level%in%c("Baseline Test and Treat", "Moderate Test and Treat", "High Test and Treat")], aes(x = scenario_level, y = mean, color = scale_drm)) +
  geom_point(size = 3, alpha = 1, aes(shape = scale_drm)) +
  # geom_point(size = 3, alpha = 1, position = position_dodge(width = .1, preserve = "total")) +
  # geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, size = 1.2, alpha = 1,
  #               position = position_dodge(width = .1, preserve = "total"))+
  theme_bw() +
  axis_combmatrix(sep = "-", levels = c("Baseline Cessation", "Baseline Prevention Interventions", "Baseline Test and Treat",
                                        "Moderate Cessation", "Moderate Prevention Interventions", "Moderate Test and Treat",
                                        "High Cessation", "High Prevention Interventions", "High Test and Treat")) +
  theme_combmatrix(combmatrix.panel.line.size = .2) +
  labs(x = "", y = "Difference in Discounted Quality-Adjusted\nLife-Years (Per Person)", color = "", shape = "") +
  theme(text = element_text(size = 14), plot.margin = unit(c(1,1,1,1), "cm"),
        legend.position = "bottom") +
  scale_color_manual(values = c("#009AB4", "#8F993E", "#E98300"))
ggsave(plot = plot, filename = "R/figures_tables/Figure_Upset_QALYs.pdf", width = 12, height = 7)

## High All

temp[scale_drm == "Baseline Excess Mortality" & test_level == "High Test and Treat" & cess_level == "High Cessation" &
       ssp_level == "High Prevention Interventions"]

## Only moderate cessation

temp[scale_drm == "Baseline Excess Mortality" & test_level == "Baseline Test and Treat" & cess_level == "Moderate Cessation" &
       ssp_level == "Baseline Prevention Interventions"]

## High test and treat and transmission reduction 

temp[scale_drm == "Baseline Excess Mortality" & test_level == "High Test and Treat" & cess_level == "Baseline Cessation" &
       ssp_level == "High Prevention Interventions"]

## Moderate all

temp[scale_drm == "Baseline Excess Mortality" & test_level == "Moderate Test and Treat" & cess_level == "Moderate Cessation" &
       ssp_level == "Moderate Prevention Interventions"]


