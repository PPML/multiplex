################
## Parameters
################

integrated_test <- FALSE

max_time <- 85*12
n_pop_nw <- 1000
num_enter <- 4
censor_time <- 2000
intervention_stop <- 120
n_pop_enter <- num_enter*(intervention_stop)
n_pop <- n_pop_nw + n_pop_enter

p_current <- 0.68
timestep <- 1/12

## Network Targets
mean_deg_eq <- 4.54
p_iso_eq <- 0.35
p_same_gender_eq <- 0.31
p_eq_and_sx <- 0.21
mean_deg_sx <- 6.46
p_iso_sx <- 0.15
p_same_gender_sx <- 0.06
duration_days_sx <- 1624
duration_days_eq <- 365*3
p_sx_pwid <- 0.54

## Interventions
p_ssp_base <- 0.53
p_condom_base <- 0.2
ssp_effect <- 0.5
condom_effect <- 0.2

## DEMOGRAPHICS
pct_male <- 0.69
pct_age_18_24 <- 0.036
pct_age_25_29 <- 0.106
pct_age_30_39 <- 0.264
pct_age_40_49 <- 0.232
pct_age_50plus <- 0.362

## INFECTIONS
p_hcv <- .437 # https://www.natap.org/2021/HCV/052821_01.htm
p_hcv_antibody <- .64 # https://www.natap.org/2021/HCV/052821_01.htm
p_acute_hcv <- .04 # https://www.natap.org/2021/HCV/052821_01.htm   
p_hiv <- 0.09 # https://www.natap.org/2021/HCV/052821_01.htm
p_hiv_hcv <- .04 # https://www.natap.org/2021/HCV/052821_01.htm

## HCV NATURAL HISTORY AND CARE CASCADE
p_ever_test_hcv <- 0.8
p_acute_clear <- 0.26

# Fibrosis progression probabilities for each cycle
f0_f1 <- rate_to_prob(prob_to_rate(0.008877), time = 1)
f1_f2 <- rate_to_prob(prob_to_rate(0.00681), time = 1)
f2_f3 <- rate_to_prob(prob_to_rate(0.0097026), time = 1)
f3_f4 <- rate_to_prob(prob_to_rate(0.0096201), time = 1)
f4_decomp <- rate_to_prob(prob_to_rate(0.005584349), time = 1)

## HIV NATURAL HISTORY AND CARE CASCADE
p_ever_test_hiv <- 0.9 
n_test_allowed <- 30
p_test_hiv <- 0.55/p_ever_test_hiv # This is annual probability from NHBS
p_test_hiv <- prob_to_rate(p_test_hiv)
p_test_hiv <- rate_to_prob(p_test_hiv, time = 1/timestep)
p_treat_hiv <- 0.56/p_ever_test_hiv # https://www.cdc.gov/hiv/pdf/library/reports/surveillance/cdc-hiv-surveillance-report-vol-26-no-2.pdf

## CESSATION
n_cess_allowed <- 30
relapse_rate <- .033
permanent_quit_prob <- 0.13
p_permltfu_hiv <- 0.12

# Create scenario table
p_test_hcv_base <- rate_to_prob(prob_to_rate(0.077/p_ever_test_hcv), time = 12)
p_treat_hcv_base <- rate_to_prob(prob_to_rate(.23*.96*.94))
p_ltfu_hiv_base <- rate_to_prob(prob_to_rate(0.3), 12)
cessation_rate_base <- 0.013
base_ssp <- 0.53

treat_hcv_levels <- c(p_treat_hcv_base, 
                      rate_to_prob(prob_to_rate((.23+.15)*.96*.94)),
                      rate_to_prob(prob_to_rate((.23+.3)*.96*.94)),
                      rate_to_prob(prob_to_rate((.23+.075)*.96*.94)), 
                      rate_to_prob(prob_to_rate((.25+.225)*.96*.94)))
test_hcv_levels <- c(p_test_hcv_base, 
                     rate_to_prob(prob_to_rate(0.077/p_ever_test_hcv)+0.75, time = 12),
                     rate_to_prob(prob_to_rate(0.077/p_ever_test_hcv)+2.25, time = 12),
                     rate_to_prob(prob_to_rate(0.077/p_ever_test_hcv)+0.25, time = 12),
                     rate_to_prob(prob_to_rate(0.077/p_ever_test_hcv)+1.5, time = 12))
ltfu_levels <- c(rate_to_prob(prob_to_rate(0.3), 12), rate_to_prob(prob_to_rate(0.125), 12), rate_to_prob(prob_to_rate(0.02), 12),
                 rate_to_prob(prob_to_rate(0.2125), 12), 
                 rate_to_prob(prob_to_rate(0.0725), 12)) 
ssp_levels <- c(base_ssp, base_ssp + 0.15, base_ssp + 0.3)
cess_levels <- c(cessation_rate_base, ((cessation_rate_base*12)+.15)/12, ((cessation_rate_base*12)+.3)/12)

combined_table <- as.data.table(expand.grid(p_ssp = ssp_levels, cessation_rate = cess_levels))

out <- NULL
for (i in 1:5) {
  temp <- copy(combined_table)
  temp <- temp[, p_treat_hcv:=treat_hcv_levels[i]]
  temp <- temp[, p_test_hcv:=test_hcv_levels[i]]
  temp <- temp[, p_ltfu_hiv:=ltfu_levels[i]]
  out <- rbind(out, temp, fill = T)
}

out <- out[, scenario_id:=seq_len(.N)]
write.csv(out, "input/final_scenario_table_full.csv", na = "", row.names = F)

out <- NULL
for (i in c(1:3)) {
  temp <- copy(combined_table)
  temp <- temp[, p_treat_hcv:=treat_hcv_levels[i]]
  temp <- temp[, p_test_hcv:=test_hcv_levels[i]]
  temp <- temp[, p_ltfu_hiv:=ltfu_levels[i]]
  out <- rbind(out, temp, fill = T)
}

out <- out[, scenario_id:=seq_len(.N)]
write.csv(out, "input/final_scenario_table_main.csv", na = "", row.names = F)