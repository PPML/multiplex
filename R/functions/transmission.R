## TRANSMISSION
count_exposures <- function(df, sx_el, eq_el, t, ssp_effect = 0.5, condom_effect = 0.2, hiv_equip_tasp = 0) {
  data <- copy(df[,.(cid, condom, hcv, hiv, ssp, eq_iso, art, current_pwid)])
  
  sx_ties <- sx_el
  sx_ties <- merge(sx_ties, data[,.(cid, hiv, condom, art)], by.x = "tail", by.y = "cid")
  sx_ties <- merge(sx_ties, data[,.(cid, hiv, condom, art)], by.x = "head", by.y = "cid", suffixes = c("V1", "V2"))
  sx_ties <- sx_ties[hivV1!=hivV2 & (artV1!=1 | is.na(artV1)) & (artV2!=1 | is.na(artV2))]
  sx_ties <- sx_ties[condomV1==0 & condomV2==0, transmit:=1]
  sx_ties <- sx_ties[(condomV1==1 & condomV2==0) | (condomV1==0 & condomV2==1), transmit:=condom_effect]
  sx_ties <- sx_ties[condomV1==1 & condomV2==1, transmit:=condom_effect^2]
  sx_ties <- sx_ties[hivV1==0 & transmit > 0, cid:=tail]
  sx_ties <- sx_ties[hivV2==0 & transmit > 0, cid:=head]
  sx_ties <- sx_ties[!is.na(cid)]
  sx_ties <- sx_ties[, sx_hiv_exposures:=sum(transmit), by = "cid"]
  sx_ties <- unique(sx_ties[,.(cid, sx_hiv_exposures)])
  
  eq_ties <- eq_el
  eq_ties <- merge(eq_ties, data[,.(cid, hiv, ssp, eq_iso, art, current_pwid)], by.x = "tail", by.y = "cid")
  eq_ties <- merge(eq_ties, data[,.(cid, hiv, ssp, eq_iso, art, current_pwid)], by.x = "head", by.y = "cid", suffixes = c("V1", "V2"))
  eq_ties <- eq_ties[hivV1!=hivV2 & current_pwidV1==1 & current_pwidV2==1 & eq_isoV1==0 & eq_isoV2==0]
  eq_ties <- eq_ties[sspV1==0 & sspV2==0 & (artV1==0 | artV2==0), transmit:=1]
  eq_ties <- eq_ties[sspV1==0 & sspV2==0 & (artV1==1 | artV2==1), transmit:=hiv_equip_tasp]
  
  eq_ties <- eq_ties[((sspV1==0 & sspV2==1) | (sspV1==1 & sspV2==0)) & (artV1==0 | artV2==0),  transmit:=ssp_effect]
  eq_ties <- eq_ties[((sspV1==0 & sspV2==1) | (sspV1==1 & sspV2==0)) & (artV1==1 | artV2==1),  transmit:=ssp_effect*hiv_equip_tasp]
  
  eq_ties <- eq_ties[(sspV1==1 & sspV2==1) & (artV1==0 | artV2==0), transmit:=(ssp_effect^2)]
  eq_ties <- eq_ties[(sspV1==1 & sspV2==1) & (artV1==1 | artV2==1), transmit:=hiv_equip_tasp*(ssp_effect^2)]
  
  eq_ties <- eq_ties[hivV1==0 & transmit > 0, cid:=tail]
  eq_ties <- eq_ties[hivV2==0 & transmit > 0, cid:=head]
  eq_ties <- eq_ties[!is.na(cid)]
  eq_ties <- eq_ties[, eq_hiv_exposures:=sum(transmit), by = "cid"]
  eq_ties <- unique(eq_ties[,.(cid, eq_hiv_exposures)])
  
  hcv_ties <- eq_el
  hcv_ties <- merge(hcv_ties, data[,.(cid, hcv, ssp, eq_iso, current_pwid)], by.x = "tail", by.y = "cid")
  hcv_ties <- merge(hcv_ties, data[,.(cid, hcv, ssp, eq_iso, current_pwid)], by.x = "head", by.y = "cid", suffixes = c("V1", "V2"))
  hcv_ties <- hcv_ties[hcvV1!=hcvV2 & current_pwidV1==1 & current_pwidV2==1 & eq_isoV1==0 & eq_isoV2==0]
  hcv_ties <- hcv_ties[(sspV1==0 & sspV2==0), transmit:=1]
  hcv_ties <- hcv_ties[((sspV1==0 & sspV2==1) | (sspV1==1 & sspV2==0)), transmit:=ssp_effect]
  hcv_ties <- hcv_ties[(sspV1==1 & sspV2==1), transmit:=ssp_effect^2]
  hcv_ties <- hcv_ties[hcvV1==0 & transmit > 0, cid:=tail]
  hcv_ties <- hcv_ties[hcvV2==0 & transmit > 0, cid:=head]
  hcv_ties <- hcv_ties[!is.na(cid)]
  hcv_ties <- hcv_ties[, eq_hcv_exposures:=sum(transmit), by = "cid"]
  hcv_ties <- unique(hcv_ties[,.(cid, eq_hcv_exposures)])
  
  return(list(sx_ties, eq_ties, hcv_ties))
}

update_exposures <- function(df, exposures) {
  data <- copy(df)
  data <- merge(data, exposures[1], by = "cid", all.x=T)
  data <- merge(data, exposures[2], by = "cid", all.x=T)
  data <- merge(data, exposures[3], by = "cid", all.x=T)
  data <- data[!is.na(sx_hiv_exposures), hiv_sx_exposure_counter:=hiv_sx_exposure_counter+sx_hiv_exposures]
  data <- data[!is.na(eq_hiv_exposures), hiv_eq_exposure_counter:=hiv_eq_exposure_counter+eq_hiv_exposures]
  data <- data[!is.na(eq_hcv_exposures), hcv_exposure_counter:=hcv_exposure_counter+eq_hcv_exposures]
  data <- data[hiv==0 & condom == 0, hiv_gensx_exposure_counter:=hiv_gensx_exposure_counter+1]
  data <- data[hiv==0 & condom == 1, hiv_gensx_exposure_counter:=hiv_gensx_exposure_counter+condom_effect]
  data <- data[, c("sx_hiv_exposures", "eq_hiv_exposures", "eq_hcv_exposures"):=NULL]
  return(data)
}

update_infections <- function(df, t) {
  data <- copy(df)
  ## Sexual - HIV
  data <- data[hiv == 0 & hiv_sx_exposure_counter>=hiv_sx_exposures, hiv_method:="Sexual"]
  data <- data[hiv == 0 & hiv_sx_exposure_counter>=hiv_sx_exposures, time_hiv_infection:=t]
  data <- data[hiv == 0 & hiv_sx_exposure_counter>=hiv_sx_exposures, hiv:=1]
  ## Equipment - HIV
  data <- data[hiv == 0 & hiv_eq_exposure_counter>=hiv_eq_exposures, hiv_method:="Equipment"]
  data <- data[hiv == 0 & hiv_eq_exposure_counter>=hiv_eq_exposures, time_hiv_infection:=t]
  data <- data[hiv == 0 & hiv_eq_exposure_counter>=hiv_eq_exposures, hiv:=1]
  ## GenPop - HIV
  data <- data[hiv == 0 & hiv_gensx_exposure_counter>=hiv_genpop_exposures, hiv_method:="GenPop"]
  data <- data[hiv == 0 & hiv_gensx_exposure_counter>=hiv_genpop_exposures, time_hiv_infection:=t]
  data <- data[hiv == 0 & hiv_gensx_exposure_counter>=hiv_genpop_exposures, hiv:=1]
  ## HCV
  data <- data[hcv == 0 & hcv_exposure_counter>=hcv_exposures, times_hcv:=times_hcv+1]
  data <- data[hcv == 0 & hcv_exposure_counter>=hcv_exposures, times_hcv_exposure_merge:=times_hcv_exposure_merge+1]
  data <- data[hcv == 0 & hcv_exposure_counter>=hcv_exposures, time_hcv_infection:=t]
  data <- data[hcv == 0 & hcv_exposure_counter>=hcv_exposures, hcv:=1]
}
