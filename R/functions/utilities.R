rate_to_prob <- function(rate, time = 1) {
  prob <- 1-exp(-rate/time)
  return(prob)
}

prob_to_rate <- function(prob) {
  rate <- -log(1-prob)
  return(rate)
}

load_life_tables <- function(file = "input/age_specific_smr_estimates_0_100_splitexcess.csv", timestep = 12, bothsex = FALSE) {
    lt <- fread(paste0(file))
    lt <- lt[,.(age, qx, current_excess_prob, former_excess_prob)]
    lt <- lt[, qx:=rate_to_prob(prob_to_rate(qx), time = timestep)]
    lt <- lt[, current_excess_prob:=rate_to_prob(prob_to_rate(current_excess_prob), time = timestep)]
    lt <- lt[, former_excess_prob:=rate_to_prob(prob_to_rate(former_excess_prob), time = timestep)]
    return(lt)
}

get_hiv_tx_status <- function(hiv_treat_df, df, t) {
  hiv_treat <- copy(hiv_treat_df)
  hiv_treat <- merge(hiv_treat, df[hiv==1,.(pid, ever_test, hiv, time_hiv_infection)], by = "pid")
  hiv_treat <- hiv_treat[, art:=ifelse(t >= (time_hiv_treat + time_hiv_infection) & t <= (time_hiv_ltfu + time_hiv_treat + time_hiv_infection), 1, 0)]
  hiv_treat <- hiv_treat[, art:=max(art, na.rm=T), by = "pid"]
  hiv_treat <- unique(hiv_treat[,.(pid, art)])
  return(hiv_treat)
}

get_pwid_status <- function(data, t, return_permanent = FALSE) {
  cess_df <- copy(data)
  cess_df <- cess_df[t >= cess_time & t <= cess_time+quit_time, quit_num:=time_block_number]
  cess_df <- cess_df[, quit_num:=mean(quit_num, na.rm=T), by = "pid"]
  cess_df <- cess_df[t >= cess_time & perm_quit_num==time_block_number, permanent_quit:=1]
  cess_df <- cess_df[t >= cess_time & perm_quit_num==time_block_number, current_pwid:=0]
  cess_df <- cess_df[, current_pwid:=mean(current_pwid, na.rm=T), by = "pid"]
  cess_df <- cess_df[, permanent_quit:=mean(permanent_quit, na.rm=T), by = "pid"]
  cess_df <- cess_df[!is.na(quit_num), current_pwid:=0]
  cess_df <- cess_df[is.na(current_pwid), current_pwid:=1]
  if (return_permanent==TRUE) {
    cess_df <- unique(cess_df[,.(pid, current_pwid, permanent_quit)])
    return(cess_df)
  } else {
    cess_df <- unique(cess_df[,.(pid, current_pwid)])
    return(cess_df)
  }
}

background_mort <- function(age) {
  m_base <- lt$qx[lt$age == floor(age)]
  m_base <- rbinom(1, 1, prob = m_base)
  return(m_base)
} 

drug_mort <- function(age, current_pwid, scalar = 1, return_prob = FALSE) {
  if (current_pwid == 1) {
    m_drug <- lt$current_excess_prob[lt$age == floor(age)]*scalar
  } else if (current_pwid == 0) {
    m_drug <- lt$former_excess_prob[lt$age == floor(age)]*scalar
  }
  if (return_prob==TRUE) {
    return(m_drug)
  } else {
    m_drug <- rbinom(1, 1, prob = m_drug)
    return(m_drug)
  }
} 

hiv_mort <- function(age, cd4, phi1 = .977, phi2 = 0.59, phi3 = 0.34, phi4 = 2.03, return_prob = TRUE, timestep = 1/12) {
  age_cat <- hiv_natural_history_age_cat(age)
  m_hiv <- (phi1^cd4)*phi2*(1+phi3*age)*phi4
  m_hiv <- rate_to_prob(m_hiv, time = 1/timestep)
  if (return_prob==TRUE) {
    return(m_hiv)
  } else {
    m_hiv <- rbinom(1, 1, prob = m_hiv)
    return(m_hiv)
  }
}

hcv_mort <- function(hcv, fibrosis_stage, m_f4 = prob_to_rate(0.027), m_decomp = prob_to_rate(0.316), m_svr_f4_decomp = 0.29, return_prob = TRUE) {
  if (hcv == 1 & fibrosis_stage == 4) {
    m_hcv <- rate_to_prob(m_f4, time = 60)
  } else if (hcv == 0 & fibrosis_stage == 4) {
    m_hcv <- rate_to_prob(m_f4*m_svr_f4_decomp, time = 60)
  } else if (hcv == 1 & fibrosis_stage == 5) {
    m_hcv <- rate_to_prob(m_decomp, time = 60)
  } else if (hcv == 0 & fibrosis_stage == 5) {
    m_hcv <- rate_to_prob(m_decomp * m_svr_f4_decomp, time = 60)
  } else {
    m_hcv <- 0
  }
  if (return_prob == TRUE) {
    return(m_hcv)
  } else {
    m_hcv <- rbinom(1, 1, prob = m_hcv)
    return(m_hcv)
  }
}