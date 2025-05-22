files <- list.files(paste0(summaries_folder, "trans_calib/"), full.names = T)

baseline <- NULL

for (f in files) {
  temp <- fread(paste0(f))
  baseline <- rbind(baseline, temp, fill = T)
}

baseline <- baseline[, mean_hiv_hcv:=n_hiv_hcv/n_pop_tot]
baseline <- baseline[, mean_hcv:=n_hcv/n_pop_tot]
baseline <- baseline[, hcv_inc:=(new_hcv_incl_acute/n_hcv_sus_active)*12]
baseline <- baseline[, mean_hiv:=n_hiv/n_pop_tot]
baseline <- baseline[, hiv_inc:=(n_new_hiv/n_hiv_sus)*12]
baseline <- baseline[is.na(hiv_inc), hiv_inc:=0]
baseline <- baseline[is.na(hcv_inc), hcv_inc:=0]

baseline <- baseline[,.(time, p_trans_hcv, rel_infect, sx_rel, mean_hiv_hcv, mean_hcv, hcv_inc, mean_hiv, hiv_inc)]
baseline <- unique(baseline[, lapply(.SD, mean), by = c("time", "p_trans_hcv", "rel_infect", "sx_rel")])
baseline <- melt(baseline, id.vars = c("p_trans_hcv", "rel_infect", "sx_rel", "time"))

baseline <- baseline[variable == "mean_hiv_hcv", target:=0.04]
baseline <- baseline[variable == "mean_hcv", target:=0.441]
baseline <- baseline[variable == "mean_hiv", target:=0.086]
baseline <- baseline[variable == "hcv_inc", target:=0.12]
baseline <- baseline[variable == "hiv_inc", target:=0.009]

baseline <- baseline[variable%in%c("mean_hiv", "mean_hcv")]

baseline <- baseline[, variable:=factor(variable, levels = c("mean_hcv", "mean_hiv", "mean_hiv_hcv", "hcv_inc", "hiv_inc"))]

baseline <- baseline[, mape:=mean(abs((target - value)/target)), by = c("p_trans_hcv", "rel_infect", "sx_rel")]
baseline <- baseline[, sum_sq:=sum((value - target)^2), by = c("p_trans_hcv", "rel_infect", "sx_rel")]
baseline <- unique(baseline[,.(p_trans_hcv, rel_infect, sx_rel, sum_sq, mape)])
baseline <- baseline[, mape_rank:=frank(mape)]
baseline <- baseline[, sse_rank:=frank(sum_sq)]
baseline[sse_rank==min(sse_rank)]
baseline[mape_rank==min(mape_rank)]

baseline <- baseline[, weight:=1/mape]

weights <- copy(baseline)
weights <- weights[, weight_norm:=weight/sum(weight)]

baseline <- NULL

for (f in files) {
  temp <- fread(paste0(f))
  baseline <- rbind(baseline, temp, fill = T)
}

baseline <- baseline[is.na(n_hcv), n_hcv:=0]
baseline <- baseline[is.na(n_hiv), n_hcv:=0]
baseline <- baseline[, mean_hcv:=n_hcv/n_pop_tot]
baseline <- baseline[, mean_hiv:=n_hiv/n_pop_tot]

baseline <- unique(merge(baseline, weights, by = c("p_trans_hcv", "rel_infect", "sx_rel")))
baseline <- baseline[, plot_group:=.GRP, by = c("p_trans_hcv", "rel_infect", "sx_rel", "stochastic_num")]

ggplot(data = baseline[time<=120], aes(x = (time-1)/12, y = mean_hcv, group = plot_group)) +
  geom_line(alpha = .025) +
  geom_line(data = baseline[weight_norm==max(weight_norm) & time<=120], color = "maroon", alpha = 1, linewidth = 2) +
  coord_cartesian(xlim = c(0, 10)) +
  theme_bw() +
  geom_hline(yintercept = 0.441, linewidth = 1.2, color = "#007C92") +
  scale_y_continuous(labels = scales::percent, n.breaks = 6) +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10), labels = c("0", "2", "4", "6", "8", "10")) +
  labs(x = "Year", y = "HCV Infection Prevalence", caption = "Note: Horizontal blue line is the calibration target. Maroon lines are from the selected parameter set.") +
  theme(text = element_text(size = 14), plot.caption = element_text(hjust = 0), plot.caption.position = "plot", plot.margin = unit(c(1,1,1,1), "cm"))
ggsave("R/figures_tables/calibration_hcv.pdf", width = 8, height = 4)

ggplot(data = baseline[time<=120], aes(x = (time-1)/12, y = mean_hiv, group = plot_group)) +
  geom_line(alpha = .025) +
  geom_line(data = baseline[weight_norm==max(weight_norm) & time<=120], color = "maroon", alpha = 1, linewidth = 2) +
  coord_cartesian(xlim = c(0, 10)) +
  theme_bw() +
  geom_hline(yintercept = 0.086, linewidth = 1.2, color = "#007C92") +
  scale_y_continuous(labels = scales::percent, n.breaks = 6) +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10), labels = c("0", "2", "4", "6", "8", "10")) +
  labs(x = "Year", y = "HIV Infection Prevalence", caption = "Note: Horizontal blue line is the calibration target. Maroon lines are from the selected parameter set.") +
  theme(text = element_text(size = 14), plot.caption = element_text(hjust = 0), plot.caption.position = "plot", plot.margin = unit(c(1,1,1,1), "cm"))
ggsave("R/figures_tables/calibration_hiv.pdf", width = 8, height = 4)

