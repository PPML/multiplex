set.seed(12345)

load(paste0(crn_folder, "baseline/fixed/1.RData"))

df <- df[time_enter==0]
df <- df[, cid:=seq_len(.N)]

## Generate Network Targets
eq_deg_target_raw <- mean_deg_eq/(1-p_iso_eq) # Mean degree among non-isolates
eq_deg_target <- (eq_deg_target_raw*(nrow(df))/2)*(1/((nrow(df)*(1-p_iso_eq)*(p_current))/nrow(df))) # Some partnerships are not active due to isolate and former
eq_gender_target <- eq_deg_target*(p_same_gender_eq)

sx_deg_target <- (mean_deg_sx*nrow(df)/2)*p_sx_pwid
sx_isolates_target <- nrow(df)*p_iso_sx
sx_gender_target <- sx_deg_target*(p_same_gender_sx)
sx_and_eq_target <- eq_deg_target*p_eq_and_sx

equipment_nw <- nw_construct(data = df)
sexual_nw <- nw_construct(data = df)

## FIT EQUIPMENT NETWORK
est_net_eq <- nw_fit(nw = equipment_nw,
                     formation = as.formula("~ edges"),
                     target.stats = c(eq_deg_target),
                     coef.diss = dissolution_coefs(dissolution = ~offset(edges),
                                                   duration = (duration_days_eq/365)*(1/timestep)),
                     edapprox = TRUE)

## EXTRACT EDGE COVARIATE
sexual_nw %n% "ec" <- as.matrix(est_net_eq$newnetwork)

est_net_sx <- nw_fit(nw = sexual_nw,
                     formation = as.formula("~ edges + edgecov('ec')"),
                     target.stats = c(sx_deg_target, sx_and_eq_target),
                     coef.diss = dissolution_coefs(dissolution = ~offset(edges), duration = (duration_days_sx/365)*(1/timestep)), edapprox = TRUE)

message("Saving SX Coefficients")
save(est_net_sx, file = paste0(crn_folder, "est_net_sx.RData"))

## EXTRACT EDGE COVARIATE
equipment_nw %n% "ec" <- as.matrix(est_net_sx$newnetwork)

est_net_eq <- nw_fit(nw = equipment_nw,
                     formation = as.formula("~ edges + edgecov('ec')"),
                     target.stats = c(eq_deg_target, sx_and_eq_target),
                     coef.diss = dissolution_coefs(dissolution = ~offset(edges), duration = (duration_days_eq/365)*(1/timestep)), edapprox = TRUE)

message("Saving EQ Coefficients")
save(est_net_eq, file = paste0(crn_folder, "est_net_eq.RData"))
