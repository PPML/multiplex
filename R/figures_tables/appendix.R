baseline <- fread("../clearance2_summaries/baseline/out_full_scenario_1.csv")

ggplot(data = baseline, aes(x = time/12, y = n_pop_tot/1000, group = interaction(stochastic_num, edge_iter_num))) +
  geom_line(alpha = .1) +
  theme_classic() +
  labs(x = "Years", y = "PWID Population Size", caption = "Note: Interventions stop and population is closed after 10 years (dashed line).") +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(labels = scales::percent, breaks = c(0, 0.25, 0.5, 0.75, 1, 1.25)) +
  theme(text = element_text(size = 16)) +
  geom_vline(aes(xintercept = 10), linetype = "dashed") +
  theme(plot.caption = element_text(hjust = 0), plot.caption.position = "plot", plot.margin = unit(c(1,1,1,1), "cm"))
ggsave("R/figures_tables/pwid_size_appendix.pdf", width = 8, height = 4)

baseline <- baseline[, current:=n_active/n_pop_tot]
baseline <- baseline[, former:=(n_former_temp)/n_pop_tot]
baseline <- baseline[, permanent_former:=(n_former_permanent)/n_pop_tot]

plot_data <- baseline[,.(time, current, former, permanent_former, stochastic_num, edge_iter_num)]
plot_data <- melt(plot_data, id.vars = c("time", "stochastic_num", "edge_iter_num"))
plot_data <- plot_data[, variable:=fcase(variable == "current", "Current Injection Drug Use",
                                         variable == "former", "Temporary Cessation",
                                         variable == "permanent_former", "Permanent Cessation")]

ggplot(data = plot_data, aes(x = time/12, y = value, group = interaction(stochastic_num, edge_iter_num, variable))) +
  geom_vline(aes(xintercept = 10), linetype = "dashed") +
  geom_line(alpha = .1, aes(color = variable)) +
  theme_classic() +
  labs(x = "Years", y = "Prevalence", color = "", caption = "Note: Interventions stop and population is closed after 10 years (dashed line).") +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(labels = scales::percent, breaks = c(0, 0.25, 0.5, 0.75, 1, 1.25), limits = c(0, 1)) +
  theme(text = element_text(size = 16), legend.position = "bottom") +
  scale_color_manual(values = c("#007C92", "#E98300", "#7F2D48")) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 3))) +
  theme(plot.caption = element_text(hjust = 0), plot.caption.position = "plot",plot.margin = unit(c(1,1,1,1), "cm"))
ggsave("R/figures_tables/pwid_prevalence_appendix.pdf", width = 8, height = 4)

baseline <- baseline[, current:=n_active]
baseline <- baseline[, former:=(n_former_temp)]
baseline <- baseline[, permanent_former:=(n_former_permanent)]

plot_data <- baseline[,.(time, current, former, permanent_former, stochastic_num, edge_iter_num)]
plot_data <- melt(plot_data, id.vars = c("time", "stochastic_num", "edge_iter_num"))
plot_data <- plot_data[, variable:=fcase(variable == "current", "Current Injection Drug Use",
                                         variable == "former", "Temporary Cessation",
                                         variable == "permanent_former", "Permanent Cessation")]

ggplot(data = plot_data, aes(x = time/12, y = value, group = interaction(stochastic_num, edge_iter_num, variable))) +
  geom_vline(aes(xintercept = 10), linetype = "dashed") +
  geom_line(alpha = .1, aes(color = variable)) +
  theme_classic() +
  labs(x = "Years", y = "Population Size", color = "", caption = "Note: Interventions stop and population is closed after 10 years (dashed line).") +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  # scale_y_continuous(labels = scales::percent, breaks = c(0, 0.25, 0.5, 0.75, 1, 1.25), limits = c(0, 1)) +
  theme(text = element_text(size = 16), legend.position = "bottom") +
  scale_color_manual(values = c("#007C92", "#E98300", "#7F2D48")) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 3))) +
  theme(plot.caption = element_text(hjust = 0), plot.caption.position = "plot", plot.margin = unit(c(1,1,1,1), "cm"))
ggsave("R/figures_tables/pwid_number_appendix.pdf", width = 8, height = 4)

trajectories <- fread("/scratch/users/mreitsma/clearance2_results/baseline/sim_1_scenario_1_edge_1_out_full.csv")
trajectories <- trajectories[current_pwid==1, pwid_status:="Current"]
trajectories <- trajectories[current_pwid==0 & permanent_quit==0, pwid_status:="Temporary Cessation"]
trajectories <- trajectories[current_pwid==0 & permanent_quit==1, pwid_status:="Permanent Cessation"]
trajectories <- trajectories[, pwid_status:=factor(pwid_status, levels = c("Current", "Temporary Cessation",
                                                                              "Permanent Cessation"))]
trajectories <- trajectories[, pid_label:=paste0("ID ", pid)]

ggplot(data = trajectories[pid<=20], aes(x = time/12, y = as.numeric(as.factor(pwid_status)))) + 
  geom_line() +
  facet_wrap(~reorder(pid_label, pid)) +
  theme_classic() +
  scale_y_continuous(labels = c("Current", "Temporary Cessation", "Permanent Cessation"), breaks = c(1, 2, 3)) +
  labs(x = "Year", y = "") +
  theme(text = element_text(size = 16))
ggsave("R/figures_tables/pwid_trajectory_appendix.pdf", width = 10, height = 8)

cd4_mort_df <- data.table(cd4 = 1:1000)
cd4_mort_df <- cd4_mort_df[, mort:=(.977^cd4)*0.59*(1+0.34*2)*2.03]

ggplot(data = cd4_mort_df, aes(x = cd4, y = mort)) + geom_line(linewidth = 1.5) +
  coord_cartesian(xlim = c(10, 800), ylim = c(0, 1.5)) +
  theme_classic() +
  labs(x = "CD4 Count", y = "Annual HIV-Related Mortality Rate") +
  geom_vline(aes(xintercept = 200), linetype = "dashed") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02)), breaks = seq(100, 800, 100)) +
  theme(text = element_text(size = 16))
ggsave("../../crn_ppml14/code/clean/figures_tables/cd4_mortality.pdf", width = 8, height = 4)

library(ggridges)

trajectories <- fread("../clearance2_results/baseline/sim_1_scenario_1_edge_1_out_full.csv")
trajectories <- trajectories[, stoch_id:=1]
temp <- fread("../clearance2_results/baseline/sim_2_scenario_1_edge_1_out_full.csv")
temp <- temp[, stoch_id:=2]
trajectories <- rbind(trajectories, temp)
temp <- fread("../clearance2_results/baseline/sim_3_scenario_1_edge_1_out_full.csv")
temp <- temp[, stoch_id:=3]
trajectories <- rbind(trajectories, temp)

trajectories <- trajectories[, pid:=paste0(pid, "_", stoch_id)]
trajectories <- trajectories[, pid_label:=paste0("ID ", pid, "_", stoch_id)]

trajectories <- trajectories[, fib_lag:=shift(fibrosis_stage), by = "pid"]
trajectories <- trajectories[fibrosis_stage != fib_lag & !is.na(fib_lag), change_time:=duration_hcv, by = "pid"]
dens_data <- unique(trajectories[!is.na(change_time),.(pid, change_time, fibrosis_stage)])

ggplot(data = dens_data[change_time/12<=70]) + 
  geom_line(data = trajectories[duration_hcv/12<=70], aes(x = duration_hcv/12, 
                                                          y = as.numeric(as.factor(fibrosis_stage)), 
                                                          group = as.factor(pid)), alpha = .1) +
  geom_density_ridges(scale = .8, stat = "density", aes(x = change_time/12, 
                                                        y = as.numeric(as.factor(fibrosis_stage)),
                                                        group = fibrosis_stage,
                                                        height = stat(density)),
                      alpha = .8, fill = "#7F2D48") +
  theme_classic() +
  scale_y_continuous(labels = c("F0", "F1", "F2", "F3", "Compensated Cirrhosis", "Decompensated Cirrhosis"),
                     breaks = c(1:6)) +
  labs(x = "Cumulative Duration of HCV Infection (Years)", y = "") +
  theme(text = element_text(size = 16))
ggsave("R/figures_tables/fibrosis_trajectory_appendix.pdf", width = 10, height = 8)


trajectories <- fread("../clearance2_results/baseline/sim_1_scenario_1_edge_1_out_full.csv")
trajectories <- trajectories[hiv==1]
trajectories <- trajectories[, pid_label:=paste0("ID ", pid)]

set.seed(1234)
plot_data <- trajectories[pid%in%sample(unique(trajectories$pid[time_hiv_infection>0]), size = 20)]
plot_data <- plot_data[, id:=.GRP, by = "pid"]
plot_data <- plot_data[, pid_label:=paste0("ID ", id)]

ggplot(data = plot_data, aes(x = duration_hiv/12, y = cd4)) + 
  geom_line() +
  facet_wrap(~reorder(pid_label, id)) +
  theme_classic() +
  labs(x = "Duration of HIV Infection (Years)", y = "CD4 Count") +
  theme(text = element_text(size = 16)) +
  geom_hline(yintercept = 200, linetype = "dashed")
ggsave("R/figures_tables/hiv_trajectory_appendix.pdf", width = 10, height = 8)


trajectories <- fread("../../output_main_final/sim_1_scenario_1_edge_1_out_full.csv")
trajectories <- trajectories[hiv==1]
trajectories <- trajectories[, pid_label:=paste0("ID ", pid)]
trajectories <- trajectories[, mean_art:=mean(art, na.rm=T), by = "pid"]

ggplot(data = trajectories[mean_art > .05 & mean_art < .25], aes(x = duration_hiv/12, y = cd4, group = pid)) + 
  geom_line(alpha = .1) +
  theme_classic() +
  labs(x = "Duration of HIV Infection", y = "CD4 Count") +
  theme(text = element_text(size = 16)) +
  geom_hline(yintercept = 200, linetype = "dashed")
ggplot(data = trajectories[mean_art < .5 & mean_art > .25], aes(x = duration_hiv/12, y = cd4, group = pid)) + 
  geom_line(alpha = .1) +
  theme_classic() +
  labs(x = "Duration of HIV Infection", y = "CD4 Count") +
  theme(text = element_text(size = 16)) +
  geom_hline(yintercept = 200, linetype = "dashed")

deaths <- fread("../../output_main_final/sim_1_scenario_1_edge_1_deaths.csv")

trajectories <- fread("../../output_main_final/sim_1_scenario_1_edge_1_out_full.csv")
trajectories <- trajectories[hiv==1]
trajectories <- trajectories[, max_duration:=max(duration_hiv), by = "pid"]
trajectories <- trajectories[, max_art:=max(art, na.rm=T), by = "pid"]
trajectories <- trajectories[is.na(max_art), max_art:=0]
trajectories <- trajectories[, mean_art:=mean(art, na.rm=T), by = "pid"]
trajectories <- unique(trajectories[,.(pid, max_duration, max_art, mean_art, time_hiv_infection)])

trajectories <- merge(trajectories, deaths[,.(pid, cod)])

ggplot(data = trajectories[time_hiv_infection>0 & cod=="HIV"], aes(x = max_duration/12)) + 
  geom_histogram() +
  facet_wrap(~max_art) +
  theme_classic() +
  labs(x = "Year", y = "CD4 Count") +
  theme(text = element_text(size = 16))

mean(trajectories$max_duration[trajectories$time_hiv_infection>0 & trajectories$max_art==0 & trajectories$cod=="HIV"])/12
mean(trajectories$max_duration[trajectories$time_hiv_infection>0 & trajectories$max_art==1 & trajectories$cod=="HIV"])/12

mean(trajectories$max_duration[trajectories$time_hiv_infection>0 & trajectories$max_art==0 & trajectories$cod=="HIV"])/12
mean(trajectories$max_duration[trajectories$time_hiv_infection>0 & trajectories$max_art==1 & trajectories$cod=="HIV"])/12
mean(trajectories$max_duration[trajectories$time_hiv_infection>0 & trajectories$mean_art>.5 & trajectories$cod=="HIV"])/12
mean(trajectories$max_duration[trajectories$time_hiv_infection>0 & trajectories$mean_art>0 & trajectories$mean_art<.5 & trajectories$cod=="HIV"])/12
mean(trajectories$max_duration[trajectories$time_hiv_infection>0 & trajectories$mean_art==0 & trajectories$cod=="HIV"])/12

## Calibration Plots
baseline <- fread("/scratch/users/mreitsma/clearance2_summaries/baseline/out_full_scenario_1.csv")

baseline <- baseline[, hcv_prev:=n_hcv/n_pop_tot]
baseline <- baseline[, hiv_prev:=n_hiv/n_pop_tot]
baseline <- baseline[, hiv_hcv_prev:=n_hiv_hcv/n_pop_tot]

baseline <- baseline[, hcv_inc:=(new_hcv_incl_acute/n_hcv_sus_active)]
baseline <- baseline[, hiv_inc:=(n_new_hiv/n_hiv_sus)]
baseline <- baseline[is.na(hcv_inc), hcv_inc:=0]
baseline <- baseline[is.na(hiv_inc), hiv_inc:=0]

baseline <- baseline[time <= 120]

baseline <- baseline[,.(hcv_prev, hiv_prev, hiv_hcv_prev, hcv_inc, hiv_inc, stochastic_num, edge_iter_num)]
baseline <- baseline[, lapply(.SD, mean), by = c("stochastic_num", "edge_iter_num")]

long <- melt(baseline, id.vars = c("stochastic_num", "edge_iter_num"))

long <- long[variable == "hiv_hcv_prev", target:=0.04]
long <- long[variable == "hcv_prev", target:=0.441]
long <- long[variable == "hiv_prev", target:=0.086]
long <- long[variable == "hcv_inc", target:=0.12]
long <- long[variable == "hiv_inc", target:=0.009]

long <- long[, variable:=factor(variable, levels = c("hcv_prev", "hiv_prev", "hiv_hcv_prev", "hcv_inc", "hiv_inc"))]

long <- long[variable == "hcv_prev", fac_lab:="HCV\nPrevalence"]
long <- long[variable == "hiv_prev", fac_lab:="HIV\nPrevalence"]
long <- long[variable == "hiv_hcv_prev", fac_lab:="HIV/HCV\nCoinfection"]
long <- long[variable == "hcv_inc", fac_lab:="HCV\nIncidence"]
long <- long[variable == "hiv_inc", fac_lab:="HIV\nIncidence"]

long <- long[variable%like%"_inc", value:=value*12]
long <- long[variable == "hcv_inc", lower:=0.084]
long <- long[variable == "hcv_inc", upper:=0.175]
long <- long[variable == "hiv_inc", lower:=0.006]
long <- long[variable == "hiv_inc", upper:=0.012]

ggplot(data = long[variable%in%c("hcv_prev", "hiv_prev", "hiv_hcv_prev")], aes(x = variable, y = value)) +
  geom_boxplot(alpha = 1, color = "darkgray", width = .4, size = 1.1) +
  geom_boxplot(aes(y = target), color = "maroon", size = 1.2) +
  theme_classic() +
  facet_wrap(~fac_lab, scales = "free_x", ncol = 4) +
  scale_y_continuous(limits = c(0, .55), expand = c(0,0), labels = scales::percent) +
  theme(text = element_text(size = 20), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(x = "", y = "Prevalence")

ggplot(data = long[variable%in%c("hcv_inc", "hiv_inc")], aes(x = variable, y = value*100)) +
  geom_boxplot(alpha = 1, color = "darkgray", width = .4, size = 1.1) +
  geom_boxplot(aes(y = target*100), color = "maroon", size = 1.2) +
  geom_ribbon(aes(xmin = variable,
                xmax = variable,
                ymin = lower*100, ymax = upper*100), alpha = .2, fill = "maroon") +
  theme_classic() +
  facet_wrap(~fac_lab, scales = "free_x", ncol = 4) +
  scale_y_continuous(limits = c(0, 15), expand = c(0,0), breaks = c(0, 2, 4, 6, 8, 10, 12, 14)) +
  theme(text = element_text(size = 20), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(x = "", y = "Incidence Per 100 Person-Years")
