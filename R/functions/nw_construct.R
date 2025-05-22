nw_construct <- function(data) {
  nw <- network.initialize(n = nrow(data), directed = FALSE)
  
  network.vertex.names(nw) <- data$cid
  nw %v% "male" <- data$male

  return(nw)
}

nw_fit <- function(nw, formation, target.stats, coef.diss, edapprox = TRUE) {
  formation <- formation
  target.stats = target.stats
  coef.diss <- coef.diss
  est_net <- netest(nw, formation, target.stats, coef.diss, edapprox = edapprox, 
                    set.control.ergm = control.ergm(MCMC.burnin = 1e6))
  
  return(est_net)
}

init_sim <- function(nw, est_net, time.burnin = 1000) {
  nwd <- simulate(nw,
                  formation = est_net$formation,
                  dissolution = est_net$coef.diss$dissolution,
                  coef.form = est_net$coef.form,
                  coef.diss = est_net$coef.diss$coef.adj,
                  time.burnin = time.burnin,
                  output = "final")
  return(nwd)
}

nw_sim <- function(nw, est_net, adj_edges, adj_edgecov, adj_diss, time.start, time.interval = 1, time.offset = 0) {
  nwd <- simulate(nw,
                  formation = est_net$formation,
                  dissolution = est_net$coef.diss$dissolution,
                  coef.form = c(edges = adj_edges, edgecov.ec = adj_edgecov),
                  coef.diss = adj_diss,
                  time.start = time.start, time.interval = time.interval, time.offset = time.offset, time.slices = 1, output = "changes",
                  control = list(MCMC.burnin.min = 1e5))
  return(nwd)
}

diss_coef_ad <- function(duration, exit_rate){
  #P(Et) overall dissolution probability
  e <- 1/duration
  #P(Nt) either node die/leave
  n <- 1-(1-exit_rate)^2
  #P(N_t) both survive
  s <- (1-exit_rate)^2
  #logit(1-(e-n)/s)=log((s-e+n)/(e-n))=log((1-e)/(e-n))
  adj_diss <- log((1-e)/(e-n))
  return(adj_diss)
}

edges_ad <- function(orig_coef, orig_n, new_n) {
  #https://statnet.org/nme/d4-s2-DynamicNets.pdf
  new_coef <- orig_coef + log(orig_n) - log(new_n)
  return(new_coef)
}

edgecov_ad <- function(orig_coef, orig_n, new_n) {
  new_coef <- (log(exp(orig_coef)*(new_n/orig_n)))
  return(new_coef)
}

return_edges_dt <- function(nw, time) {
  dt <- as.data.table(get.dyads.active(nw, at = time))
  colnames(dt) <- c("tail", "head")
  return(dt)
}

attach_edges <- function(data, sx, eq, time, n_sx, n_eq) {
  # New Individuals
  new_eq_cid <- data$cid[data$time_enter == time & data$eq_iso == 0]
  new_eq_el <- NULL  
  if (length(new_eq_cid > 0)) {
    new_eq_el <- data.table(tail = sample(x = data$cid[!(data$cid%in%new_eq_cid) & data$eq_iso==0], size = length(new_eq_cid)*n_eq, replace = F),
                            head = rep(new_eq_cid, n_eq))
    
    # Joint SX & EQ  
    joint_el <- new_eq_el[, .SD[sample(.N, 1)], by = "head"] 
  } 
  
  new_sx_cid <- data$cid[data$time_enter==time & data$eq_iso==1]
  new_sx_el <- NULL
  if (length(new_sx_cid) > 0) {
    new_sx_el <- data.table(tail = sample(x = data$cid[!(data$cid%in%new_sx_cid)], size = length(new_sx_cid), replace = F),
                            head = rep(new_sx_cid, n_sx))
  }
  if (length(new_eq_cid > 0)) {
    new_sx_el <- rbind(new_sx_el, joint_el)
  }
  return(list(new_eq_el, new_sx_el))
}

nw_from_el <- function(data, el) {
  nw <- nw_construct(data)
  nw <- network::add.edges(nw, tail = el$tail, head = el$head)
  
  return(nw)
}

update_cid <- function(el, cid_update) {
  el <- merge(el, cid_update, by.x = "tail", by.y = "cid", all.x = T)
  el <- el[!is.na(new_cid)]
  el <- el[, tail:=new_cid]
  el <- el[, new_cid:=NULL]
  el <- merge(el, cid_update, by.x = "head", by.y = "cid", all.x = T)
  el <- el[!is.na(new_cid)]
  el <- el[, head:=new_cid]
  el <- el[, new_cid:=NULL]
  return(el)
}




