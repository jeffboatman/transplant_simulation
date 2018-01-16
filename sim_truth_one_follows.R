# ------------------------------------------------------------
# simulation code for estimating true survival distribution
# if one patient follows regime.
# R version 3.1.3
# ------------------------------------------------------------

# Simulation parameters ---------------------------------------------------
n.mc <- 1e4             # number of Monte Carlo Iterations
n.cores <- 4            # number of cores to use for parallel processing
# -------------------------------------------------------------------------

options(digits = 3)
library(parallel)
library(survival)
library(dplyr)


# designed to run on 5 servers
host <- system2("hostname", stdout = TRUE)
hosts <- paste0(c("carbon", "cesium", "chromium", 
  "potassium", "silicon"), 
  ".ccbr.umn.edu")
n.s <- length(hosts)
j <- match(host, hosts)



# source("parameters.R")


expit <- function(x) {
  expit <- exp(x)/(1+exp(x))
  return(expit)
}

outfunc <- function(x) {
  outfunc <- c(outer(x,x))
  return(outfunc)
}

mlast <- function(x) {
  mlast <- x[-length(x)]
  return(mlast)
}

last <- function(x) { 
  return(x[length(x)]) 
}

offered_func <- function(x) {
  offered_func <- cumprod(c(1, x))
  offered_func <- mlast(offered_func)
}

sim <- function(sim_num, days_in_study) {
	
	n <- 20000 # base number of patient and organ arrivals
	tau <- 730 # length of follow-up after listing
	knot_pts <- c(0, 180, 365, 540, 730)
	burn_in <- 10000 # length of burn-in period
	max_cal_time <- burn_in + days_in_study
	
	# blood types and their associated probabilities
	b_types <- c("O", "A", "B", "AB")
	b_probs <- c(0.45, 0.4, 0.11, 0.04)
	
	# set seed
	set.seed(sim_num + 1101985)
	
	
	# Simulate data -----------------------------------------------
	
	# Generate organ arrival process. Organs continue to arrive until the last 
    # patient arrives on WL 
	# o_arv_times <- round(cumsum(rexp(n, 0.35)))
	o_arv_times <- round(cumsum(rexp(n, 0.32)))
	o_arv_times <- o_arv_times[o_arv_times < max_cal_time]
	n_o <- length(o_arv_times)
	o_blood_type <- sample(b_types, n_o, replace = T, prob = b_probs)
	good_organ <- rbinom(n_o, 1, 1/2)
	
	#  Generate patient arrivals, death time in absence of transplant, and 
    # patient characteristics 
	# (blood type and time-varying covariate)
	p_arv_times <- round(cumsum(rexp(n, 0.5)))
	p_arv_times <- p_arv_times[p_arv_times < max_cal_time]
	n_p <- length(p_arv_times)
	p_blood_type <- sample(b_types, n_p, replace = TRUE, prob = b_probs)
	blood_type <- sample(b_types, n_p, replace = T, prob = b_probs)

    # which subject to follow
    lucky_id <- sample(which(p_arv_times > burn_in), 1)
	
	# max times we need to consider
	max_p_arv_time <- max(p_arv_times)
	months_flup <- ceiling(max_p_arv_time / 30)
	
	# Generate time-varying covariate process and death in absence of transplant
	cov_int <- rnorm(n_p, -1, 1)
	cov_slope <- rnorm(n_p, 1/365, (1/365)/4)
	
	mk_death_time <- function(i) {
	  cov_month <- (0:months_flup)*30* cov_slope[i] + cov_int[i]
	  death_time_month <- rexp(months_flup + 1, exp(cov_month)/365)
	  death_time_interval <- which(death_time_month < 30)
	  if(length(death_time_interval) > 0) 
	  {
	    death_time_interval <- death_time_interval[1] 
	    death_time_i <- (death_time_interval - 1)*30 + death_time_month[death_time_interval]
	  } else {
	    death_time_i <- months_flup*30
	  }
	  death_time_i	
	}
	
	death_time <- sapply(1:n_p, mk_death_time)
	
	# create vectors for death time (patient and calendar time) and 
    # removal from waitlist
	rm_time <- p_arv_times + death_time
	rm_pt_time <- death_time
	rm_type <- rep(1, n_p)
	ev_time <- rm_time
	ev_pt_time <- death_time
	ev_ind <- rep(1, n_p)
	tx_ed <- rep(0, n_p)
	organ_received <- rep(0, n_p)
	
	# create vectors for first organ eligable and ineligable
	first_organ <- sapply(1:n_p, function(x) {
	  pt_arv <- p_arv_times[x]
	  pt_rm <- rm_time[x]
	  organ_avail <- which(o_arv_times >= pt_arv & o_arv_times < pt_rm)
	  if(length(organ_avail > 0)) {
	    organ_avail[1]
	  } else {
	    NA
	  }
	})
	
	last_organ <- sapply(1:n_p, function(x) {
	  pt_arv <- p_arv_times[x]
	  pt_rm <- rm_time[x]
	  organ_avail <- which(o_arv_times >= pt_arv & o_arv_times < pt_rm)
	  if(length(organ_avail > 0)) {
	    last(organ_avail)
	  } else {
	    NA
	  }
	})
	
	# matrix for whether or not a person would have accepted organ.
    # leave as NA if not observed
	accept.matrix <- matrix(NA, nrow = n_p, ncol = n_o)
	
	# create dataset for organ allocation process
	data_set <- list()
	
	for(i in seq_along(o_arv_times)) {
	  compatible <- rep(0,n_p)
	  exact <- rep(0,n_p)
	  o_blood_type_i <- o_blood_type[i]
	  good_organ_i <- good_organ[i]
	  time_i <- o_arv_times[i]
	  
	  if(o_blood_type_i == "O") {
	    compatible <- rep(1,n_p)
	    exact[p_blood_type == "O"] <- 1
	  }
	  
	  if(o_blood_type_i == "A") {
	    compatible[(p_blood_type == "A")|(p_blood_type == "AB")] <- 1
	    exact[p_blood_type == "A"] <- 1
	  }
	  if(o_blood_type_i == "B") {
	    compatible[(p_blood_type == "B")|(p_blood_type == "AB")] <- 1
	    exact[p_blood_type == "B"] <- 1
	  }
	  if(o_blood_type_i == "AB") {
	    compatible[(p_blood_type == "AB")] <- 1
	    exact[p_blood_type == "AB"] <- 1
	  }
	  
	  cov_value <- floor((time_i - p_arv_times)/30) * 30 * cov_slope + cov_int
	  prob_accept <- expit(-2.5 + 0.25 * cov_value)
	  accept <- rbinom(n_p, 1, prob_accept)

      # force lucky_obs to decline bad organs.
      if(!good_organ_i) accept[lucky_id] <- 0

	  elig <- which((p_arv_times <= time_i) & 
        (p_arv_times + death_time >= time_i) & 
	      tx_ed == 0 & compatible == 1 & accept == 1)	
	  if(length(elig) > 0) 
	  {
	    who_got_it <- elig[order(-exact[elig], -cov_value[elig])[1]]
	    cov_value_match <- cov_value[who_got_it]
	    exact_match <- exact[who_got_it]
	    
	    offered <- which((p_arv_times <= time_i) & 
          (p_arv_times + death_time >= time_i) & tx_ed == 0
	      & compatible  == 1 & ((exact > exact_match) | 
	          (exact == exact_match & cov_value >= cov_value_match)))
	    offered_vec <- rep(0,n_p)
	    offered_vec[offered] <- 1
	    elig_to_receive <- which((p_arv_times <= time_i) & 
          (p_arv_times + death_time >= time_i) & 
	        tx_ed == 0 & compatible == 1)
	    # Figure out the probability of being offered organ 
	    prob_offered <- rep(0,n_p)
	    order_use <- elig_to_receive[order(-exact[elig_to_receive], 
          -cov_value[elig_to_receive])]
	    prob_offered[order_use] <- mlast(cumprod(c(1, 1-prob_accept[order_use])))
	    
	    tx_ed_i <- who_got_it
	    tx_ed[tx_ed_i] <- 1
	    rm_time[tx_ed_i] <- time_i
	    rm_pt_time[tx_ed_i] <- time_i - p_arv_times[tx_ed_i]
	    rm_type[tx_ed_i] <- 2
	    last_organ[tx_ed_i] <- i
	    organ_received[tx_ed_i] <- i
	    
	    # generate post-transplant survival
	      # in this scenario good and bad organs differ
	    post_tx_surv <- rexp(1, ifelse(good_organ_i == 1, exp(-2.5)/365, 
        exp(-0.5)/365)) 
	      # in this scenario good and bad organs do not differ
	    #post_tx_surv <- rexp(1,exp(-2.5)/365) 
	    ev_time[tx_ed_i] <- time_i + post_tx_surv
	    ev_pt_time[tx_ed_i] <- time_i + post_tx_surv - p_arv_times[tx_ed_i]
	    
	    # determine who we know was willing to accept organ
	    accept.matrix[offered_vec == 1, i] <- accept[offered_vec == 1] 
	    
	    data_set_i <- cbind(data.frame(cov_value, prob_accept, prob_offered, 
          offered_vec, tx_ed, exact,
	      i, time_i)[elig_to_receive, ], elig_to_receive)
	    data_set[[i]] <- data_set_i
	  }
	  elig_to_receive <-  which((p_arv_times <= time_i) & 
          (p_arv_times + death_time >= time_i) & 
	      tx_ed == 0 & compatible == 1 )
	  if(length(elig) == 0 & length(elig_to_receive) > 0) {
	    offered_vec <- rep(0, n_p)
	    offered_vec[elig_to_receive] <- 1
	    
	    # Figure out the probability of being offered organ 
	    prob_offered <- rep(0,n_p)
	    order_use <- elig_to_receive[order(-exact[elig_to_receive], 
          -cov_value[elig_to_receive])]
	    prob_offered[order_use] <- mlast(cumprod(c(1, 1-prob_accept[order_use])))
	    
	    # determine who we know was willing to accept organ
	    accept.matrix[offered_vec == 1, i] <- accept[offered_vec == 1] 
	    
	    data_set_i <- cbind(data.frame(cov_value, prob_accept, prob_offered, 
          offered_vec, tx_ed, exact,
	      i, time_i)[elig_to_receive, ], elig_to_receive)	
	    data_set[[i]] <- data_set_i
	  }
	}
	data_set <- bind_rows(data_set)
	
	colnames(data_set) <- c("cov_value", "prob_accept", "prob_offered", 
      "offered", "tx_ed",
	  "exact_blood", "organ_num", "cal_time","id")
	data_set <- data_set[order(data_set$organ_num, -data_set$exact_blood, 
      -data_set$cov_value),]
	
	#  We assume that follow-up continues until the max_cal_time and then 
    # obs are censored
	rm_type <- ifelse(rm_time > max_cal_time, 0, rm_type)
	rm_pt_time <- ifelse(rm_time > max_cal_time, max_cal_time - p_arv_times, 
      rm_pt_time)
	rm_time <- ifelse(rm_time > max_cal_time, max_cal_time, rm_time)
	
	ev_ind <- ifelse(ev_time > max_cal_time, 0, ev_ind)
	ev_pt_time <- ifelse(ev_time > max_cal_time, max_cal_time - p_arv_times, 
      ev_pt_time)
	ev_time <- ifelse(ev_time > max_cal_time, max_cal_time, ev_time)
	
	
	# Combine the patient, organ, and time-varying information together
	pt_data_set <- data.frame(
	  id = 1:n_p, 
	  p_arv_times, 
	  p_blood_type, 
	  dth_time_notx = death_time, 
	  rm_time,
	  rm_type,
	  rm_pt_time, 
	  ev_time, 
	  ev_pt_time, 
	  ev_ind, 
	  ever_tx_ed = tx_ed,
	  cov_int,
	  cov_slope,
	  first_organ,
	  last_organ,
	  organ_received)
	
	org_data_set <- data.frame(
	  organ_num = 1:n_o, 
	  o_blood_type, 
	  o_arv_times, 
	  good_organ) 
	
	# data_set_final <- left_join(data_set, pt_data_set, by = "id")
    data_set_final <- left_join(pt_data_set, data_set, by = "id")
	data_set_final <- left_join(data_set_final, org_data_set, by = "organ_num")	
    data_set_final <- as.data.frame(data_set_final)
    data_set_final <- data_set_final[with(data_set_final, order(id, o_arv_times)), ]
	
    lucky_data <- subset(data_set_final, id == lucky_id, 
      select = c(id, p_arv_times, ev_pt_time, ev_ind, ever_tx_ed, tx_ed, good_organ))
    # lucky_data <- data.frame(sim.num = sim_num, lucky_data)
    
    lucky_data <- lucky_data[nrow(lucky_data), ]
    lucky_data <- data.frame(sim.num = sim_num, lucky_data)
    out <- lucky_data
    write.table(out, outfile, quote = FALSE, row.names = FALSE,
      col.names = FALSE, append = TRUE)
}

outfile <- paste0("./outfiles/sim_truth_one_follows_out_", j, ".txt")
cols <- c("sim.num", "id", "p_arv_times", "ev_pt_time", "ev_ind", "ever_tx_ed", 
  "tx_ed", "good_organ")



write.table(t(cols), 
  file = outfile, 
  row.names = FALSE,   
  col.names = FALSE,
  append = FALSE, 
  quote = FALSE)

sims <- ((j * n.mc / n.s) - (n.mc / n.s - 1)):(j * n.mc / n.s)
out.list <- mclapply(sims,
  sim,
  mc.cores = n.cores,
  days_in_study = 20* 365)
