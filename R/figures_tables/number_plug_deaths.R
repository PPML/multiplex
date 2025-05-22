
files <- list.files(path = "/scratch/users/mreitsma/clearance2_results/baseline/", pattern = "deaths")
files <- files[files%like%"scenario_1_"]

out <- NULL
for (f in files) {
  temp <- fread(paste0("/scratch/users/mreitsma/clearance2_results/baseline/", f))
  out <- rbind(out, temp, fill = T)
}

out <- out[analytic==1]

prop.table(table(out$cod))

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

summary <- out[, drm:=ifelse(hiv==1 & hcv==1, 1, 0)] %>%
  group_by(stochastic_num, edge_iter_num) %>%
  summarize(drm = mean(drm))
mean(summary$drm)
quantile(summary$drm, c(.025, .975))

summary <- out[, drm:=ifelse(hiv==1 & cod=="HCV", 1, 0)] %>%
  filter(cod == "HCV") %>%
  group_by(stochastic_num, edge_iter_num) %>%
  summarize(drm = mean(drm))
mean(summary$drm)
quantile(summary$drm, c(.025, .975))

summary <- out[, drm:=ifelse(hcv==1 & cod=="HIV", 1, 0)] %>%
  filter(cod == "HIV") %>%
  group_by(stochastic_num, edge_iter_num) %>%
  summarize(drm = mean(drm))
mean(summary$drm)
quantile(summary$drm, c(.025, .975))

## 125 & 63
files <- list.files(path = "/scratch/users/mreitsma/clearance2_results/baseline/", pattern = "deaths")
files <- files[files%like%"scenario_1_"]

out <- NULL
for (f in files) {
  temp <- fread(paste0("/scratch/users/mreitsma/clearance2_results/baseline/", f))
  out <- rbind(out, temp, fill = T)
}

out <- out[analytic==1]

out <- out[cod == "HIV", n_hiv_deaths:=.N, by = c("edge_iter_num", "stochastic_num")]
out <- out[cod == "HCV", n_hcv_deaths:=.N, by = c("edge_iter_num", "stochastic_num")]
out <- out[,.(edge_iter_num, stochastic_num, n_hiv_deaths, n_hcv_deaths)]
out <- out[, lapply(.SD, mean, na.rm=T), by = c("edge_iter_num", "stochastic_num")]

baseline <- out

## 125 & 63
files <- list.files(path = "/scratch/users/mreitsma/clearance2_results/baseline/", pattern = "deaths")
files <- files[files%like%"scenario_27_"]

out <- NULL
for (f in files) {
  temp <- fread(paste0("/scratch/users/mreitsma/clearance2_results/baseline/", f))
  out <- rbind(out, temp, fill = T)
}

out <- out[analytic==1]

out <- out[cod == "HIV", n_hiv_deaths:=.N, by = c("edge_iter_num", "stochastic_num")]
out <- out[cod == "HCV", n_hcv_deaths:=.N, by = c("edge_iter_num", "stochastic_num")]
out <- out[,.(edge_iter_num, stochastic_num, n_hiv_deaths, n_hcv_deaths)]
out <- out[, lapply(.SD, mean, na.rm=T), by = c("edge_iter_num", "stochastic_num")]

high <- copy(out)

files <- list.files(path = "/scratch/users/mreitsma/clearance2_results/baseline/", pattern = "deaths")
files <- files[files%like%"scenario_14_"]

out <- NULL
for (f in files) {
  temp <- fread(paste0("/scratch/users/mreitsma/clearance2_results/baseline/", f))
  out <- rbind(out, temp, fill = T)
}

out <- out[analytic==1]

out <- out[cod == "HIV", n_hiv_deaths:=.N, by = c("edge_iter_num", "stochastic_num")]
out <- out[cod == "HCV", n_hcv_deaths:=.N, by = c("edge_iter_num", "stochastic_num")]
out <- out[,.(edge_iter_num, stochastic_num, n_hiv_deaths, n_hcv_deaths)]
out <- out[, lapply(.SD, mean, na.rm=T), by = c("edge_iter_num", "stochastic_num")]

moderate <- copy(out)

setnames(high, c("n_hiv_deaths", "n_hcv_deaths"), c("hiv_high", "hcv_high"))
setnames(moderate, c("n_hiv_deaths", "n_hcv_deaths"), c("hiv_mod", "hcv_mod"))

baseline <- merge(baseline, high, by = c("edge_iter_num", "stochastic_num"))
baseline <- merge(baseline, moderate, by = c("edge_iter_num", "stochastic_num"))

mean((baseline$hiv_high-baseline$n_hiv_deaths)/baseline$n_hiv_deaths)
quantile((baseline$hiv_high-baseline$n_hiv_deaths)/baseline$n_hiv_deaths, .025)
quantile((baseline$hiv_high-baseline$n_hiv_deaths)/baseline$n_hiv_deaths, .975)

mean((baseline$hcv_high-baseline$n_hcv_deaths)/baseline$n_hcv_deaths)
quantile((baseline$hcv_high-baseline$n_hcv_deaths)/baseline$n_hcv_deaths, .025)
quantile((baseline$hcv_high-baseline$n_hcv_deaths)/baseline$n_hcv_deaths, .975)

mean(out$age)
prop.table(table(out$cod))

mean(out$age[out$cod=="HCV"])
mean(out$age[out$cod=="HIV"])
mean(out$age[out$cod=="Drug"])
mean(out$age[out$cod=="Background"])
mean(out$age[out$hiv==1 & out$hcv==1])
nrow(out[hiv==1 & hcv==1])/nrow(out)
nrow(out[hiv==1 & hcv==1 & ever_test==0])/nrow(out)

nrow(out[cod=="HCV" & hiv==1])/nrow(out[cod=="HCV"])


