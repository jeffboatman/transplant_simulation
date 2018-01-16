# Simulation parameters ---------------------------------------------------
n.mc <- 1000             # number of Monte Carlo Iterations
n.cores <- 4            # number of cores to use for parallel processing
# -------------------------------------------------------------------------

options(digits = 3)
library(parallel)

# designed to run on 5 servers
host <- system2("hostname", stdout = TRUE)
hosts <- paste0(c("carbon", "cesium", "chromium", 
  "potassium", "silicon"), 
  ".ccbr.umn.edu")
n.s <- length(hosts)
j <- match(host, hosts)

library(survival)
# library(rms)
library(dplyr)
print.tbl_df <- print.data.frame

quiet <- function(f) {
  function(...) {
    suppressMessages(f(...))
  }
}

l_join <- quiet(left_join)
r_join <- quiet(right_join)

# miscellaneous functions
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
  return( x[length(x)] ) 
}

offered_func <- function(x) {
  offered_func <- cumprod(c(1, x))
  offered_func <- mlast(offered_func)
}

sample_df <- function(df) {
  which_rows <- sort(sample(1:nrow(df), replace = TRUE))
  df[which_rows, ]
}

`%!in%` <- Negate(`%in%`)

first.tx.fun <- function(x) {
  cumsum(x) == 1
}


sim <- function(sim_num, days_in_study, n.boot = 100, include.se = FALSE) {
	### define parameters for simulation --------------------------------------------------------------
	n <- 20000 # base number of patient and organ arrivals
	tau <- 730 # length of follow-up after listing
	knot_pts <- c(0,180,365,540,730)
	burn_in <- 10000 # length of burn-in period
	max_cal_time <- burn_in + days_in_study
	
	# blood types and their associated probabilities
	b_types <- c("O", "A", "B", "AB")
	b_probs <- c(0.45, 0.4, 0.11, 0.04)
	
	# set seed
	set.seed(sim_num + 1101985)
	
	
	# Simulate data -----------------------------------------------------------------------------------
	
	# Generate organ arrival process. Organs continue to arrive until the last patient arrives on WL 
	# o_arv_times <- round(cumsum(rexp(n, 0.35)))
	o_arv_times <- round(cumsum(rexp(n, 0.32)))
	o_arv_times <- o_arv_times[o_arv_times < max_cal_time]
	n_o <- length(o_arv_times)
	o_blood_type <- sample(b_types, n_o, replace = T, prob = b_probs)
	good_organ <- rbinom(n_o, 1, 1/2)
	
	#  Generate patient arrivals, death time in absence of transplant, and patient characteristics 
	# (blood type and time-varying covariate)
	p_arv_times <- round(cumsum(rexp(n, 0.5)))
	p_arv_times <- p_arv_times[p_arv_times < max_cal_time]
	n_p <- length(p_arv_times)
	p_blood_type <- sample(b_types, n_p, replace = TRUE, prob = b_probs)
	blood_type <- sample(b_types, n_p, replace = T, prob = b_probs)
	
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
	
	# create vectors for death time (patient and calendar time) and removal from waitlist
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
	
	# matrix for whether or not a person would have accepted organ leave as NA if not observed
	accept.matrix <- matrix(NA, nrow = n_p, ncol = n_o)
	
	# create dataset for organ allocation process
	start <- proc.time()
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
	  elig <- which((p_arv_times <= time_i) & (p_arv_times + death_time >= time_i) & 
	      tx_ed == 0 & compatible == 1 & accept == 1)	
	  if(length(elig) > 0) 
	  {
	    who_got_it <- elig[order(-exact[elig], -cov_value[elig])[1]]
	    cov_value_match <- cov_value[who_got_it]
	    exact_match <- exact[who_got_it]
	    
	    offered <- which((p_arv_times <= time_i) & (p_arv_times + death_time >= time_i) & tx_ed == 0
	      & compatible  == 1 & ((exact > exact_match) | 
	          (exact == exact_match & cov_value >= cov_value_match)))
	    offered_vec <- rep(0,n_p)
	    offered_vec[offered] <- 1
	    elig_to_receive <- which((p_arv_times <= time_i) & (p_arv_times + death_time >= time_i) & 
	        tx_ed == 0 & compatible == 1)
	    # Figure out the probability of being offered organ 
	    prob_offered <- rep(0,n_p)
	    order_use <- elig_to_receive[order(-exact[elig_to_receive], -cov_value[elig_to_receive])]
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
	    post_tx_surv <- rexp(1, ifelse(good_organ_i == 1, exp(-2.5)/365, exp(-0.5)/365)) 
	      # in this scenario good and bad organs do not differ
	    #post_tx_surv <- rexp(1,exp(-2.5)/365) 
	    ev_time[tx_ed_i] <- time_i + post_tx_surv
	    ev_pt_time[tx_ed_i] <- time_i + post_tx_surv - p_arv_times[tx_ed_i]
	    
	    # determine who we know was willing to accept organ
	    accept.matrix[offered_vec == 1, i] <- accept[offered_vec == 1] 
	    
	    data_set_i <- cbind(data.frame(cov_value, prob_accept, prob_offered, offered_vec, tx_ed, exact,
	      i, time_i)[elig_to_receive, ], elig_to_receive)
	    data_set[[i]] <- data_set_i
	  }
	  elig_to_receive <-  which((p_arv_times <= time_i) & (p_arv_times + death_time >= time_i) & 
	      tx_ed == 0 & compatible == 1 )
	  if(length(elig) == 0 & length(elig_to_receive) > 0) {
	    offered_vec <- rep(0, n_p)
	    offered_vec[elig_to_receive] <- 1
	    
	    # Figure out the probability of being offered organ 
	    prob_offered <- rep(0,n_p)
	    order_use <- elig_to_receive[order(-exact[elig_to_receive], -cov_value[elig_to_receive])]
	    prob_offered[order_use] <- mlast(cumprod(c(1, 1-prob_accept[order_use])))
	    
	    # determine who we know was willing to accept organ
	    accept.matrix[offered_vec == 1, i] <- accept[offered_vec == 1] 
	    
	    data_set_i <- cbind(data.frame(cov_value, prob_accept, prob_offered, offered_vec, tx_ed, exact,
	      i, time_i)[elig_to_receive, ], elig_to_receive)	
	    data_set[[i]] <- data_set_i
	  }
	}
	data_set <- bind_rows(data_set)
	
	colnames(data_set) <- c("cov_value", "prob_accept", "prob_offered", 
      "offered", "tx_ed", "exact_blood", "organ_num", "cal_time","id")
	data_set <- data_set[order(data_set$organ_num, -data_set$exact_blood, -data_set$cov_value),]
	
	#  We assume that follow-up continues until the max_cal_time and then obs are censored
	rm_type <- ifelse(rm_time > max_cal_time, 0, rm_type)
	rm_pt_time <- ifelse(rm_time > max_cal_time, max_cal_time - p_arv_times, rm_pt_time)
	rm_time <- ifelse(rm_time > max_cal_time, max_cal_time, rm_time)
	
	ev_ind <- ifelse(ev_time > max_cal_time, 0, ev_ind)
	ev_pt_time <- ifelse(ev_time > max_cal_time, max_cal_time - p_arv_times, ev_pt_time)
	ev_time <- ifelse(ev_time > max_cal_time, max_cal_time, ev_time)
	
	# No longer round these times to make life easier
	#ev_time <- round(ev_time)
	#ev_pt_time <- round(ev_pt_time)
	
	# Combine the patient, organ, and time-varying information together
	pt_data_set <- data.frame(
	  id = 1:n_p, # n -> n_p 
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
	
	data_set_final <- l_join(data_set, pt_data_set, by = "id")
	data_set_final <- l_join(data_set_final, org_data_set, by = "organ_num")	

    dsf <- data_set_final
	
	#  Obtain estimates of parameters in the logistic model for accpeting organ -----------------------
	data4model <- filter(data_set_final, cal_time >= burn_in & offered == 1) 
	m1 <- glm(tx_ed ~ cov_value, family = "binomial", data = data4model)
	logit_param <- m1$coef
	
	# Get estimated probability of accepting organ, being offered organ, receiving tx,
	# and receiving particular treatment of each subject
	data_set_final$prob_accept_est <- as.numeric(expit(cbind(1, data_set_final$cov_value) %*% 
	    logit_param))
	data_set_final <- data_set_final[order(data_set_final$organ_num, -data_set_final$exact_blood,
	  -data_set_final$cov_value),]
	data_set_final$prob_offered_est <- unlist(by(1-data_set_final$prob_accept_est, 
	  data_set_final$organ_num, offered_func, simplify = TRUE))
	data_set_final <- mutate(data_set_final, 
	  prob_tx_est = prob_accept_est*prob_offered_est,
	  prob_tx_received_est = tx_ed*prob_tx_est + (1 - tx_ed)*(1 - prob_tx_est))
	data_set_final <- data_set_final[order(data_set_final$id, data_set_final$cal_time), ]
	
	# Create a dataset with time-varying weights ------------------------------------------------------
	
	data_set_final <- data_set_final[order(data_set_final$id, data_set_final$cal_time),]
	data_set_final <- mutate(data_set_final, pt_time = cal_time - p_arv_times)
	pt_first_obs <- data.frame(id = 1:n_p, p_arv_times, p_blood_type, dth_time_notx = death_time, 
	  rm_time, rm_pt_time, ev_time, ev_pt_time, ev_ind, ever_tx_ed = tx_ed, pt_time = 0)
	data_set_final <- merge(data_set_final, pt_first_obs, all = TRUE)
	data_set_final<- data_set_final[order(data_set_final$id, data_set_final$pt_time),]
	data_set_final <- mutate(data_set_final, 
	  prob_tx_received_est = ifelse(is.na(prob_tx_received_est), 1, 
        prob_tx_received_est),
      prob_accept_est = ifelse(is.na(prob_accept_est), 0, prob_accept_est))
	
	# Get cumulative probability of following treatment trajectory
	data_set_final <- group_by(data_set_final, id)
	data_set_final <- mutate(data_set_final, pt_time_l = lead(pt_time), 
	  # den_weight = cumprod(prob_tx_received_est),
	  compliant = 1 - tx_ed * (1 - good_organ),
	  compliant = ifelse(is.na(compliant), 1, compliant)) 
	data_set_final <- as.data.frame(data_set_final)
	
	# create dataset in counting process format
	# not quite sure why we had so many -1's here
	data_set_final <- mutate(data_set_final, start = pt_time, 
	  stop = ifelse(is.na(pt_time_l)==FALSE, pt_time_l, ev_pt_time),
	  event_window =  ifelse(is.na(pt_time_l)==FALSE, 0, ev_ind),
	  cal_start = start + p_arv_times, cal_stop = stop + p_arv_times,
      tx_ed = ifelse(is.na(tx_ed), 0, tx_ed),
      c_weight = ((1-tx_ed )/(1-prob_accept_est))^((1-good_organ)*offered),
      c_weight = ifelse(is.na(c_weight), 1, c_weight),
      c_weight2 = ((1-tx_ed )/(1-prob_tx_est)) ^ ((1-good_organ)),
      c_weight2 = ifelse(is.na(c_weight2), 1, c_weight2))

    # deal with multiple tx per day
    data_set_final$sameday <- with(data_set_final, start == stop)
    data_set_final <- data_set_final[with(data_set_final, order(id, start, -sameday)), ]
    data_set_final <- group_by(data_set_final, id)
    data_set_final <- mutate(data_set_final, 
      den_weight = cumprod(prob_tx_received_est),
      ipw  = cumprod(c_weight), 
      ipw2 = cumprod(c_weight2))
    data_set_final <- as.data.frame(data_set_final)
    data_set_final$sameday <- NULL
    data_set_final <- subset(data_set_final, start < stop)
	
	# limit information to only within tau years post-listing survival
	# and information after the burn-in period of burn_in days
	data_set_final <- filter(data_set_final, start < tau)
	data_set_final <- mutate(data_set_final, 
	  event_window = ifelse(stop > tau, 0, event_window), 
	  cal_stop = ifelse(stop > tau, tau + p_arv_times, cal_stop),
	  stop = ifelse(stop > tau, tau, stop))
	
	data_set_final <- subset(data_set_final, cal_stop >= burn_in)
	data_set_final <- mutate(data_set_final, 
	  start = ifelse(cal_start < burn_in, start + (burn_in - cal_start), start))
	data_set_final <- subset(data_set_final, start < stop)

    #exclude participants transplanted prior to burn-in
    fo <- which(o_arv_times > burn_in)[1]
    data_set_final <- subset(data_set_final, !(ever_tx_ed & organ_received < fo))

	
	# Estimate the probability of following treatment trajectory if all subjects followed the
	# regime of interest ------------------------------------------------------------------------------
	
	# subjects observed after burn-in period
	unique_ids <- unique(data_set_final$id)
	n.subj <- length(unique_ids)
	# patient data for those patients that are actually observed in study
	pt_obs <- filter(pt_data_set, id %in% unique_ids)
	
	# Get a matrix with covariate value on each day after the burn-in period --------------------------
	# If the subject is not observed on that particular calendar day replace with NA ------------------
	
	# add 1 to months because you could overlap 2 months
	months_in_study <- ceiling((days_in_study+1)/30)+1
	months_start_study <- floor((burn_in - pt_obs$p_arv_times + 1)/30)
	months_subject <- unlist(sapply(months_start_study, 
	  function(x) {(x:(x + months_in_study - 1))*30}))
	
	# covariate by month
	pt_cov_by_month <- matrix(pt_obs$cov_int, nrow = n.subj , ncol = months_in_study) + 
	  matrix(pt_obs$cov_slope, nrow = n.subj , ncol = months_in_study)*
	  matrix(months_subject, nrow = n.subj , ncol = months_in_study, byrow = TRUE)
	# number of days in month starting on burn-in day
	days_subject <- sapply(1:n.subj, function(x) {
	  days_start <- 30
	  if(pt_obs$p_arv_times[x] <= burn_in) {
	    days_start <- 30 - ((burn_in - pt_obs$p_arv_times[x])%%30)
	  }
	  days_subject <- c(rep(pt_cov_by_month[x, 1], days_start), 
	    rep(pt_cov_by_month[x, 2:months_in_study], each = 30))
	  days_subject <- days_subject[1:(days_in_study+1)]
	  return(days_subject)
	})
	pt_cov_by_day <- t(days_subject)
	
	## Replace non-active calendar days with NA
	for (i in 1: n.subj){
	  if((floor(pt_obs$rm_time[i]) - burn_in + 1) < (days_in_study+1)) {
	    days_active <- max(floor(pt_obs$rm_time[i]) - burn_in + 1, 1)
	    pt_cov_by_day[i, (days_active + 1) : (days_in_study + 1)] <- NA
	  }
	  if(pt_obs$p_arv_times[i] > burn_in) {
	    pt_cov_by_day[i, 1:(pt_obs$p_arv_times[i] - burn_in) ] <- NA
	  }
	  #if((i %% 100)==0) {print(i)}	
	}
	
	# Impute covariate trajectory if never transplanted using hot deck imputation --------------------- 
	
	pt_cov_by_day_impute <- pt_cov_by_day
	start <- proc.time()   
	
	for (i in 1:n.subj) {
	  impute.need <- (pt_obs$rm_type[i] == 2 &  pt_obs$rm_time[i] >= burn_in)
	  first.time <- 1
	  if(pt_obs$p_arv_times[i] > burn_in) {
	    first.time <- pt_obs$p_arv_times[i] - burn_in
	  }
	  for (q in 1:100) {
	    if (impute.need==T) {
	      # minimum time with NAs
     	  min.time <- which(is.na(pt_cov_by_day_impute[i,]) == T & 
	          1:(days_in_study + 1) > first.time)[1]
	      study.day <- min.time + burn_in - 1 - pt_obs$p_arv_times[i] + 1
	      # column in dataset corresponding to pt day = study.day
	      # pt_obs$p_arv_times + study.day - 1 = calendar day corresponding to study.day in study 
	      # - burn_in + 1 = finding column in dataset 
	      cal.day.compare <- pt_obs$p_arv_times + study.day -1 - burn_in + 1 
	      possible.compare.group <- (1:n.subj)[(cal.day.compare - 1) > 0 & 
	          cal.day.compare <= (days_in_study +1) ]
	      t1 <- cbind(possible.compare.group, cal.day.compare[possible.compare.group])
	      compare.group <- possible.compare.group[which(is.na(pt_cov_by_day[t1]) == FALSE)]
	      if(length(compare.group) > 0 ){
	        data.compare.group <- pt_cov_by_day[ cbind(compare.group, 
	          (cal.day.compare - 1)[compare.group])]
	        data.i <-  pt_cov_by_day_impute[i, c(min.time - 1)]
	        difference <- abs(data.compare.group-data.i)
	        candidate <- (compare.group)[which(difference == min(difference))]
	        candidate.data <- pt_cov_by_day[candidate, ]
	        candidate.data.start <- cal.day.compare[candidate]
	        length.impute <- (days_in_study + 1) - max(min.time, candidate.data.start)
	        pt_cov_by_day_impute[i, min.time:(min.time + length.impute)] <- 
	          candidate.data[candidate.data.start:(candidate.data.start + length.impute)]
	        impute.need <- (pt_obs$rm_type[candidate] == 2)
	      }	    
	    }
	  }
	  # if(i %% 100 == 0) {print(i)}
	}
	proc.time() - start
	
	# Simulate organ allocation if everyone followed regime of interest -------------------------------
	tx_ed_s1 <- rep(0, n.subj)
	o_use <- which(o_arv_times >= burn_in & good_organ == 1)
	accept.obs <- accept.matrix[which(pt_data_set$id %in% unique_ids), ]
	num.WL <- NULL
	data_set_sim <- list()
	
	for(i in o_use) {
	  compatible <- rep(0,n.subj)
	  exact <- rep(0,n.subj)
	  o_blood_type_i <- o_blood_type[i]
	  good_organ_i <- good_organ[i]
	  time_i <- o_arv_times[i]
	  
	  if(o_blood_type_i == "O") {
	    compatible <- rep(1,n.subj)
	    exact[pt_obs$p_blood_type == "O"] <- 1
	  }
	  
	  if(o_blood_type_i == "A") {
	    compatible[(pt_obs$p_blood_type == "A")|
	        (pt_obs$p_blood_type == "AB")] <- 1
	    exact[pt_obs$p_blood_type == "A"] <- 1
	  }
	  if(o_blood_type_i == "B") {
	    compatible[(pt_obs$p_blood_type == "B")|
	        (pt_obs$p_blood_type == "AB")] <- 1
	    exact[pt_obs$p_blood_type == "B"] <- 1
	  }
	  if(o_blood_type_i == "AB") {
	    compatible[(pt_obs$p_blood_type == "AB")] <- 1
	    exact[pt_obs$p_blood_type == "AB"] <- 1
	  }
	  
	  cov_value <- pt_cov_by_day_impute[, (time_i - burn_in  + 1)]
	  prob_accept <- expit(logit_param[1] + logit_param[2] * cov_value)
	  accept <- ifelse(is.na(accept.obs[, i]) == FALSE, accept.obs[, i], 
	    rbinom(n.subj, 1, prob_accept))
	  elig <- which(is.na(prob_accept) == FALSE & 
	      tx_ed_s1 == 0 & compatible == 1 & accept == 1)	
	  if(length(elig) > 0) {
	    who_got_it <- elig[order(-exact[elig], -cov_value[elig])[1]]
	    cov_value_match <- cov_value[who_got_it]
	    exact_match <- exact[who_got_it]
	    
	    offered <- which(is.na(cov_value) == FALSE & 
	        tx_ed_s1 == 0 & compatible  == 1 &
	        ((exact > exact_match) | (exact == exact_match & cov_value >= cov_value_match)))
	    offered_vec <- rep(0,n.subj)
	    offered_vec[offered] <- 1
	    elig_to_receive <- which(is.na(cov_value) == FALSE  &
	        tx_ed_s1 == 0 & compatible == 1)
	    # Figure out the probability of being offered organ 
	    prob_offered <- rep(0, n.subj)
	    order_use <- elig_to_receive[order(-exact[elig_to_receive], -cov_value[elig_to_receive])]
	    prob_offered[order_use] <- mlast(cumprod(c(1, 1-prob_accept[order_use])))
	    
	    tx_ed_s1_i <- who_got_it
	    tx_ed_s1[tx_ed_s1_i] <- 1
	    
	    data_set_i <- cbind(data.frame(cov_value, prob_accept, prob_offered, offered_vec, tx_ed_s1,
	      exact , i, time_i)[elig_to_receive,],elig_to_receive)
	    data_set_sim[[i]] <- data_set_i
	    # if(i %% 100 ==0) {print(i)}
	  }
	  elig_to_receive <-  which(is.na(cov_value) == FALSE & tx_ed_s1 == 0 & compatible == 1 )
	  num.WL <- c(num.WL, length(which(is.na(cov_value) == FALSE & tx_ed_s1 == 0)))
	  
	  if(length(elig) == 0 & length(elig_to_receive) > 0) {
	    offered_vec <- rep(0, n.subj)
	    offered_vec[elig_to_receive] <- 1
	    # Figure out the probability of being offered organ 
	    prob_offered <- rep(0,n.subj)
	    order_use <- elig_to_receive[order(-exact[elig_to_receive], -cov_value[elig_to_receive])]
	    prob_offered[order_use] <- mlast(cumprod(c(1, 1-prob_accept[order_use])))
	    
	    data_set_i <- cbind(data.frame(cov_value, prob_accept, prob_offered, offered_vec, tx_ed_s1,
	      exact, i, time_i)[elig_to_receive,],elig_to_receive)
	    data_set_sim[[i]] <- data_set_i
	    # if(i %% 100 ==0) {print(i)}
	  }
	}
	data_set_sim <- bind_rows(data_set_sim)
	
	colnames(data_set_sim) <- c("cov_value", "prob_accept", "prob_offered", "offered", "tx_ed",
	  "exact_blood", "organ_num", "cal_time","id")
	data_set_sim <- data_set_sim[order(data_set_sim$organ_num, -data_set_sim$exact_blood, 
	  -data_set_sim$cov_value),]
	
	
	# Calculate probability of following treatment trajectory if all subjects followed the regime of 
	# interest by estimating the probability of being offered organ given simulation above. 
	# -------------------------------------------------------------------------------------------------
	
	data_set_final_1 <- filter(data_set_final, p_arv_times >= burn_in)
	unique_ids_1 <- unique(data_set_final_1$id)
	n.subj_1 <- length(unique_ids_1)
	pt_obs_1 <- filter(pt_data_set, id %in% unique_ids_1)
	pt_cov_by_day_impute_s <- pt_cov_by_day_impute[unique_ids %in% unique_ids_1, ]
	
	tolong <- function(id, first_organ, last_organ) {
	  if(is.na(first_organ)) return(NULL)
	  cbind(id, first_organ:last_organ)
	}
	
	pt_obs_long <- mapply(tolong, pt_obs_1$id, pt_obs_1$first_organ, pt_obs_1$last_organ)
	pt_obs_long <- data.frame(do.call(rbind, pt_obs_long))
	names(pt_obs_long) <- c("id", "organ_num")
	pt_obs_long <- merge(pt_obs_long, pt_obs_1, by = "id")
	pt_obs_long$o_arv_times <- o_arv_times[pt_obs_long$organ_num] 
	covlookup <- cbind(match(pt_obs_long$id, pt_obs_1$id), 
	  pt_obs_long$o_arv_times - burn_in  + 1)
	pt_obs_long$cov_value <- pt_cov_by_day_impute_s[covlookup]
	pt_obs_long$prob_accept <- predict(m1, newdata = pt_obs_long, type = "response")
	pt_obs_long$orig <- TRUE
	pt_obs_long$o_blood_type <- o_blood_type[pt_obs_long$organ_num]
	pt_obs_long <- mutate(pt_obs_long, 
	  compatible =  (o_blood_type == "O") + 
	  (o_blood_type == "A")*((p_blood_type == "A")|(p_blood_type == "AB")) + 
	  (o_blood_type == "B")*((p_blood_type == "B")|(p_blood_type == "AB")) + 
	  (o_blood_type == "AB")*(p_blood_type == "AB"),
	  exact_blood = as.numeric(o_blood_type == p_blood_type))
	
	pt_obs_long <- subset(pt_obs_long, compatible == 1, 
	  select = c(id, cov_value, prob_accept, exact_blood, organ_num))
	pt_obs_long$orig <- TRUE 
	
	dsm <- subset(data_set_sim, 
	  select = c(id, cov_value, prob_accept, exact_blood, organ_num))
	dsm$orig <- FALSE
	
	ptll <- split(pt_obs_long, pt_obs_long$id)
	
	calc_prob_offered <- function(df) {
	  # odf <- subset(dsm, organ_num %in% df$organ_num)
      odf <- dsm[dsm$organ_num %in% df$organ_num, ]    
	  df <- rbind(df, odf)
	  df <- df[order(df$organ_num, -df$exact_blood, -df$cov_value, -df$orig), ]
	  df$prob_offered <- unlist(tapply(1 - df$prob_accept, df$organ_num, offered_func))
	  subset(df, orig)
	}
	
	dsw <- bind_rows(lapply(ptll, calc_prob_offered))
	
	dsw$good_organ <- good_organ[dsw$organ_num]
	dsw$cal_time <- o_arv_times[dsw$organ_num]
	dsw <- merge(dsw, subset(pt_obs_1, select = c(id, organ_received)),
	  by = "id")
	dsw <- mutate(dsw, prob_accept = ifelse(!good_organ, 0, prob_accept),
	                   prob_offered = ifelse(!good_organ, 1, prob_offered),
	                   tx_ed = as.numeric(organ_num == organ_received))
	dsw <- subset(dsw, select = c(cov_value, prob_accept, prob_offered, tx_ed, 
	  organ_num, cal_time, id, good_organ))

	dsw <- mutate(dsw, 
	  prob_tx_est = prob_accept*prob_offered,
	  prob_tx_received_est = tx_ed*prob_tx_est + (1 - tx_ed)*(1 - prob_tx_est))
	dsw <- dsw[order(dsw$id, dsw$cal_time), ]
	
	# Create a dataset with time-varying weights ------------------------------------------------------
	# Here the weights are based on probability of following treatement trajectory if everyone followed
	# regime of interest ------------------------------------------------------------------------------
	
    dsw <- dsw[with(dsw, order(id, cal_time, -organ_num)), ]
    dsw <- group_by(dsw, id, cal_time)
    dsw <- mutate(dsw, 
      prob_tx_received_est = rev(cumprod(prob_tx_received_est)))
    dsw <- as.data.frame(dsw)
    dsw <- dsw[!duplicated(dsw[c("id", "cal_time")]), ]

	
	# create dataset in counting process format
	pt_first_obs <- data.frame(id = unique_ids_1,  cal_time = pt_obs_1$p_arv_times)
	dsw <- merge(dsw, pt_first_obs, all = TRUE)
	dsw <- dsw [order(dsw$id, dsw$cal_time),]
	dsw <- mutate(dsw, 
	  prob_tx_received_est = ifelse(is.na(prob_tx_received_est), 1, prob_tx_received_est))
	
	# Get cumulative probability of following treatment trajectory
	dsw <- group_by(dsw, id)
	dsw <- mutate(dsw, num_weight = cumprod(prob_tx_received_est)) 
	dsw <- as.data.frame(dsw)
	
	# Merge with dataset with the denominator weights, create final weights
	dsw <- select(dsw, cal_time, id, num_weight)
	data_set_final_1 <- mutate(data_set_final_1 , 
      cal_time = ifelse(is.na(cal_time) == TRUE, p_arv_times, cal_time))
	data_set_final_1 <- r_join(dsw, data_set_final_1 )
	data_set_final_1 <- mutate(data_set_final_1, total_weight = num_weight/den_weight)
	
	# Use IPCW Kaplan-Meier to estimate marginal survival ---------------------------------------------
	surv1 <- survfit(Surv(start, stop, event_window) ~ 1, 
      data = data_set_final_1, 
	  weight = total_weight)
	t1main <- summary(surv1, times = c(180, 360, 540, 720))

	surv3 <- survfit(Surv(start, stop, event_window) ~ 1, 
      data = data_set_final_1, 
	  weights = ipw2)
	t3main <- summary(surv3, times = c(180, 360, 540, 720))

	surv4 <- survfit(Surv(start, stop, event_window) ~ 1, 
      data = data_set_final_1, 
	  weights = compliant)
	t4main <- summary(surv4, times = c(180, 360, 540, 720))

    ### begin bootstrap ----------------------------------
	# boots <- matrix(0, nrow = n.boot, ncol = 12)
    
    boots <- matrix(nrow = n.boot, ncol = 14)
    if(include.se) {
      for (b in seq_len(n.boot)){
        dsf.boot <- dsf[order(dsf$id, dsf$cal_time), ]
        dsf.boot <- group_by(dsf.boot, id)
        dsf.boot <- mutate(dsf.boot, 
          pt_time = cal_time - p_arv_times,
          start = pt_time,
          pt_time_l = lead(pt_time), 
          stop = ifelse(!is.na(pt_time_l), pt_time_l, 
            ev_pt_time),
          event_window =  ifelse(!is.na(pt_time_l), 0, 
            ev_ind))
        dsf.boot <- as.data.frame(dsf.boot)
        dsf.boot <- sample_df(dsf.boot)
        
        logit.data.boot <- filter(dsf.boot, cal_time >= burn_in & offered == 1) 
        m1.boot <- glm(tx_ed ~ cov_value, family = "binomial", data = logit.data.boot)
        logit.coef.boot <- m1.boot$coef
        
        dsf.boot$prob_accept_est <- predict(m1.boot, newdata = dsf.boot, type = "response")
        dsf.boot <- dsf.boot[order(dsf.boot$organ_num, -dsf.boot$exact_blood, -dsf.boot$cov_value),]
        dsf.boot$prob_offered_est <- unlist(by(1-dsf.boot$prob_accept_est, dsf.boot$organ_num, offered_func, simplify = TRUE))
        dsf.boot <- mutate(dsf.boot, 
          prob_tx_est = prob_accept_est * 
            prob_offered_est,
          prob_tx_received_est = tx_ed * prob_tx_est + 
           (1 - tx_ed) * (1 - prob_tx_est))
        
        dsf.boot <- dsf.boot[order(dsf.boot$id, dsf.boot$cal_time),]
        pt.first.obs.boot <- data.frame(
          id = 1:n_p, 
          p_arv_times, 
          p_blood_type, 
          dth_time_notx = death_time, 
          rm_time, 
          rm_pt_time, 
          ev_time, 
          ev_pt_time, 
          ev_ind, 
          ever_tx_ed = tx_ed, 
          pt_time = 0,
          start = 0)
        # pt.first.obs.boot <- subset(pt.first.obs.boot, id %in% dsf.boot$id)
        dsf.boot <- merge(dsf.boot, pt.first.obs.boot, all = TRUE)
        dsf.boot <- dsf.boot[order(dsf.boot$id, dsf.boot$start), ]
                
        dsf.boot$ran <- runif(nrow(dsf.boot))
        dsf.boot.g <- group_by(dsf.boot, id)
        dsf.boot.g <- mutate(dsf.boot.g, 
          pt_time_l = lead(pt_time), 
          event_window = ifelse(is.na(stop), 0, 
            event_window),
          stop = ifelse(is.na(stop), pt_time_l, stop),
          cal_start = start + p_arv_times, 
          cal_stop = stop + p_arv_times,
          tx_ed = ifelse(is.na(tx_ed), 0, tx_ed),
          first.tx = first.tx.fun(tx_ed),
          event_window = ifelse(tx_ed == 1 & !first.tx, 
            0, event_window),
          stop = ifelse(tx_ed == 1 & !first.tx, start, 
            stop),
          drop = ifelse(tx_ed == 1 & !first.tx & ran < 
            1 / 2, TRUE, FALSE),
          tx_ed = ifelse(tx_ed == 1 & first.tx, 1, 0),
          compliant = 1 - tx_ed * (1 - good_organ),
          compliant = ifelse(is.na(compliant), 1, 
            compliant),
          prob_tx_received_est = tx_ed * prob_tx_est + (1 - tx_ed) * (1 - prob_tx_est),
          prob_tx_received_est =  
            ifelse(is.na(prob_tx_received_est), 1, 
            prob_tx_received_est),
          c_weight = ((1 - tx_ed ) / (1 - 
            prob_accept_est))^((1 - good_organ) * 
            offered),
          c_weight = ifelse(is.na(c_weight), 1, 
            c_weight),
          c_weight2 = ((1 - tx_ed )/(1 - prob_tx_est)) ^ 
            ((1 - good_organ)),
          c_weight2 = ifelse(is.na(c_weight2), 1, 
            c_weight2)
        )
        dsf.boot <- as.data.frame(dsf.boot.g)
        dsf.boot <- subset(dsf.boot, !drop)
        dsf.boot <- dsf.boot[order(dsf.boot$id, dsf.boot$start, dsf.boot$stop), ]
        dsf.boot$ran <- dsf.boot$drop <- dsf.boot$first.tx <- NULL
        dsf.boot$sameday <- with(dsf.boot, start == stop)
        dsf.boot <- dsf.boot[with(dsf.boot, order(id, start, -sameday)), ]
        dsf.boot <- group_by(dsf.boot, id)
        dsf.boot <- mutate(dsf.boot, 
                           den_weight = cumprod(prob_tx_received_est),
                           ipw  = cumprod(c_weight), 
                           ipw2 = cumprod(c_weight2))
        dsf.boot <- as.data.frame(dsf.boot)
        dsf.boot$sameday <- NULL
        dsf.boot <- subset(dsf.boot, start < stop)
        
        dsf.boot <- filter(dsf.boot, start < tau)
        dsf.boot <- mutate(dsf.boot, 
                           event_window = ifelse(stop > tau, 0, event_window), 
                           cal_stop = ifelse(stop > tau, tau + p_arv_times, cal_stop),
                           stop = ifelse(stop > tau, tau, stop))
        
        dsf.boot <- subset(dsf.boot, cal_stop >= burn_in)
        dsf.boot <- mutate(dsf.boot, 
                           start = ifelse(cal_start < burn_in, start + (burn_in - cal_start), start))
        dsf.boot <- subset(dsf.boot, start < stop)
        
        fo <- which(o_arv_times > burn_in)[1]
        dsf.boot <- subset(dsf.boot, !(ever_tx_ed & organ_received < fo))
        
        unique.ids.boot <- unique(dsf.boot$id)
        n.subj <- length(unique.ids.boot)
        pt.obs.boot <- filter(pt_data_set, id %in% unique.ids.boot)
        
        months_start_study <- floor((burn_in - pt.obs.boot$p_arv_times + 1) / 30)
        months_subject <- unlist(sapply(months_start_study, 
                                        function(x) {(x:(x + months_in_study - 1)) * 30}))
        
        cov.by.mont.boot <- matrix(pt.obs.boot$cov_int, nrow = n.subj, ncol = months_in_study) + 
          matrix(pt.obs.boot$cov_slope, nrow = n.subj , ncol = months_in_study) *
          matrix(months_subject, nrow = n.subj , ncol = months_in_study, byrow = TRUE)
        
        days_subject <- sapply(1:n.subj, function(x) {
          days_start <- 30
          if(pt.obs.boot$p_arv_times[x] <= burn_in) {
            days_start <- 30 - ((burn_in - pt.obs.boot$p_arv_times[x]) %% 30)
          }
          days_subject <- c(rep(cov.by.mont.boot[x, 1], days_start), 
                            rep(cov.by.mont.boot[x, 2:months_in_study], each = 30))
          days_subject <- days_subject[1:(days_in_study+1)]
          return(days_subject)
        })
        cov.by.day.boot <- t(days_subject)
        
        for (i in 1: n.subj){
          if((floor(pt.obs.boot$rm_time[i]) - burn_in + 1) < (days_in_study+1)) {
            days_active <- max(floor(pt.obs.boot$rm_time[i]) - burn_in + 1, 1)
            cov.by.day.boot[i, (days_active + 1) : (days_in_study + 1)] <- NA
          }
          if(pt.obs.boot$p_arv_times[i] > burn_in) {
            cov.by.day.boot[i, 1:(pt.obs.boot$p_arv_times[i] - burn_in) ] <- NA
          }
        }
        
        cov.by.day.impute.boot <- cov.by.day.boot
        
        for (i in 1:n.subj) {
          impute.need <- (pt.obs.boot$rm_type[i] == 2 &  pt.obs.boot$rm_time[i] >= burn_in)
          first.time <- 1
          if(pt.obs.boot$p_arv_times[i] > burn_in) {
            first.time <- pt.obs.boot$p_arv_times[i] - burn_in
          }
          for (q in 1:100) {
            if (impute.need==T) {
              min.time <- which(is.na(cov.by.day.impute.boot[i,]) == T & 
                                  1:(days_in_study + 1) > first.time)[1]
              study.day <- min.time + burn_in - 1 - pt.obs.boot$p_arv_times[i] + 1
              cal.day.compare <- pt.obs.boot$p_arv_times + study.day -1 - burn_in + 1 
              possible.compare.group <- (1:n.subj)[(cal.day.compare - 1) > 0 & 
                                                     cal.day.compare <= (days_in_study +1) ]
              t1 <- cbind(possible.compare.group, cal.day.compare[possible.compare.group])
              compare.group <- possible.compare.group[which(is.na(cov.by.day.boot[t1]) == FALSE)]
              if(length(compare.group) > 0 ){
                data.compare.group <- cov.by.day.boot[ cbind(compare.group, 
                                                             (cal.day.compare - 1)[compare.group])]
                data.i <-  cov.by.day.impute.boot[i, c(min.time - 1)]
                difference <- abs(data.compare.group-data.i)
                candidate <- (compare.group)[which(difference == min(difference))]
                candidate.data <- cov.by.day.boot[candidate, ]
                candidate.data.start <- cal.day.compare[candidate]
                length.impute <- (days_in_study + 1) - max(min.time, candidate.data.start)
                cov.by.day.impute.boot[i, min.time:(min.time + length.impute)] <- 
                  candidate.data[candidate.data.start:(candidate.data.start + length.impute)]
                impute.need <- (pt.obs.boot$rm_type[candidate] == 2)
              }	    
            }
          }
        }
        
        
        tx_ed_s1 <- rep(0, n.subj)
        o_use <- which(o_arv_times >= burn_in & good_organ == 1)
        accept.obs.boot <- accept.matrix[which(pt_data_set$id %in% unique.ids.boot), ]
        dss.boot <- list()
        
        for(i in o_use) {
          compatible <- rep(0,n.subj)
          exact <- rep(0,n.subj)
          o_blood_type_i <- o_blood_type[i]
          good_organ_i <- good_organ[i]
          time_i <- o_arv_times[i]
          
          if(o_blood_type_i == "O") {
            compatible <- rep(1,n.subj)
            exact[pt.obs.boot$p_blood_type == "O"] <- 1
          }
          
          if(o_blood_type_i == "A") {
            compatible[(pt.obs.boot$p_blood_type == "A")|
                         (pt.obs.boot$p_blood_type == "AB")] <- 1
            exact[pt.obs.boot$p_blood_type == "A"] <- 1
          }
          if(o_blood_type_i == "B") {
            compatible[(pt.obs.boot$p_blood_type == "B")|
                         (pt.obs.boot$p_blood_type == "AB")] <- 1
            exact[pt.obs.boot$p_blood_type == "B"] <- 1
          }
          if(o_blood_type_i == "AB") {
            compatible[(pt.obs.boot$p_blood_type == "AB")] <- 1
            exact[pt.obs.boot$p_blood_type == "AB"] <- 1
          }
          
          cov_value <- cov.by.day.impute.boot[, (time_i - burn_in  + 1)]
          prob_accept <- expit(logit.coef.boot[1] + logit.coef.boot[2] * cov_value)
          accept <- ifelse(is.na(accept.obs.boot[, i]) == FALSE, accept.obs.boot[, i], 
                           rbinom(n.subj, 1, prob_accept))
          elig <- which(is.na(prob_accept) == FALSE & 
                          tx_ed_s1 == 0 & compatible == 1 & accept == 1)	
          if(length(elig) > 0) {
            who_got_it <- elig[order(-exact[elig], -cov_value[elig])[1]]
            cov_value_match <- cov_value[who_got_it]
            exact_match <- exact[who_got_it]
            
            offered <- which(is.na(cov_value) == FALSE & 
                               tx_ed_s1 == 0 & compatible  == 1 &
                               ((exact > exact_match) | (exact == exact_match & cov_value >= cov_value_match)))
            offered_vec <- rep(0,n.subj)
            offered_vec[offered] <- 1
            elig_to_receive <- which(is.na(cov_value) == FALSE  &
                                       tx_ed_s1 == 0 & compatible == 1)
            prob_offered <- rep(0, n.subj)
            order_use <- elig_to_receive[order(-exact[elig_to_receive], -cov_value[elig_to_receive])]
            prob_offered[order_use] <- mlast(cumprod(c(1, 1-prob_accept[order_use])))
            
            tx_ed_s1_i <- who_got_it
            tx_ed_s1[tx_ed_s1_i] <- 1
            
            data_set_i <- cbind(data.frame(cov_value, prob_accept, prob_offered, offered_vec, tx_ed_s1,
                                           exact , i, time_i)[elig_to_receive,],elig_to_receive)
            dss.boot[[i]] <- data_set_i
          }
          elig_to_receive <-  which(is.na(cov_value) == FALSE & tx_ed_s1 == 0 & compatible == 1 )
          
          if(length(elig) == 0 & length(elig_to_receive) > 0) {
            offered_vec <- rep(0, n.subj)
            offered_vec[elig_to_receive] <- 1
            # Figure out the probability of being offered organ 
            prob_offered <- rep(0,n.subj)
            order_use <- elig_to_receive[order(-exact[elig_to_receive], -cov_value[elig_to_receive])]
            prob_offered[order_use] <- mlast(cumprod(c(1, 1-prob_accept[order_use])))
            
            data_set_i <- cbind(data.frame(cov_value, prob_accept, prob_offered, offered_vec, tx_ed_s1,
                                           exact, i, time_i)[elig_to_receive,],elig_to_receive)
            dss.boot[[i]] <- data_set_i
          }
        }
        dss.boot <- bind_rows(dss.boot)
        
        colnames(dss.boot) <- c("cov_value", "prob_accept", "prob_offered", "offered", "tx_ed",
                                "exact_blood", "organ_num", "cal_time","id")
        dss.boot <- dss.boot[order(dss.boot$organ_num, -dss.boot$exact_blood, 
                                   -dss.boot$cov_value),]
        
        dsf.final.boot <- filter(dsf.boot, p_arv_times >= burn_in)
        unique_ids_1 <- unique(dsf.final.boot$id)
        n.subj_1 <- length(unique_ids_1)
        pt.obs.1.boot <- filter(pt_data_set, id %in% unique_ids_1)
        cov.by.day.impute.s.boot <- cov.by.day.impute.boot[unique.ids.boot %in% unique_ids_1, ]
        
        tolong <- function(id, first_organ, last_organ) {
          if(is.na(first_organ)) return(NULL)
          cbind(id, first_organ:last_organ)
        }
        
        pt.obs.long.boot <- mapply(tolong, pt.obs.1.boot$id, pt.obs.1.boot$first_organ, pt.obs.1.boot$last_organ)
        pt.obs.long.boot <- data.frame(do.call(rbind, pt.obs.long.boot))
        names(pt.obs.long.boot) <- c("id", "organ_num")
        pt.obs.long.boot <- merge(pt.obs.long.boot, pt.obs.1.boot, by = "id")
        pt.obs.long.boot$o_arv_times <- o_arv_times[pt.obs.long.boot$organ_num] 
        covlookup <- cbind(match(pt.obs.long.boot$id, pt.obs.1.boot$id), 
                           pt.obs.long.boot$o_arv_times - burn_in  + 1)
        pt.obs.long.boot$cov_value <- cov.by.day.impute.s.boot[covlookup]
        pt.obs.long.boot$prob_accept <- predict(m1.boot, newdata = pt.obs.long.boot, type = "response")
        pt.obs.long.boot$orig <- TRUE
        pt.obs.long.boot$o_blood_type <- o_blood_type[pt.obs.long.boot$organ_num]
        pt.obs.long.boot <- mutate(pt.obs.long.boot, 
                                   compatible =  (o_blood_type == "O") + 
                                     (o_blood_type == "A")*((p_blood_type == "A")|(p_blood_type == "AB")) + 
                                     (o_blood_type == "B")*((p_blood_type == "B")|(p_blood_type == "AB")) + 
                                     (o_blood_type == "AB")*(p_blood_type == "AB"),
                                   exact_blood = as.numeric(o_blood_type == p_blood_type))
        
        pt.obs.long.boot <- subset(pt.obs.long.boot, compatible == 1, 
                                   select = c(id, cov_value, prob_accept, exact_blood, organ_num))
        pt.obs.long.boot$orig <- TRUE 
        
        dsm.boot <- subset(dss.boot, 
                           select = c(id, cov_value, prob_accept, exact_blood, organ_num))
        dsm.boot$orig <- FALSE
        
        ptll.boot <- split(pt.obs.long.boot, pt.obs.long.boot$id)
        
        calc_prob_offered <- function(df) {
          odf <- dsm.boot[dsm.boot$organ_num %in% df$organ_num, ]    
          df <- rbind(df, odf)
          df <- df[order(df$organ_num, -df$exact_blood, -df$cov_value, -df$orig), ]
          df$prob_offered <- unlist(tapply(1 - df$prob_accept, df$organ_num, offered_func))
          subset(df, orig)
        }
        
        dsw.boot <- bind_rows(lapply(ptll.boot, calc_prob_offered))
        
        dsw.boot$good_organ <- good_organ[dsw.boot$organ_num]
        dsw.boot$cal_time <- o_arv_times[dsw.boot$organ_num]
        dsw.boot <- merge(dsw.boot, subset(pt.obs.1.boot, select = c(id, organ_received)),
                          by = "id")
        dsw.boot <- mutate(dsw.boot, prob_accept = ifelse(!good_organ, 0, prob_accept),
                           prob_offered = ifelse(!good_organ, 1, prob_offered),
                           tx_ed = as.numeric(organ_num == organ_received))
        dsw.boot <- subset(dsw.boot, select = c(cov_value, prob_accept, prob_offered, tx_ed, 
                                                organ_num, cal_time, id, good_organ))
        
        dsw.boot <- mutate(dsw.boot, 
                           prob_tx_est = prob_accept*prob_offered,
                           prob_tx_received_est = tx_ed*prob_tx_est + (1 - tx_ed)*(1 - prob_tx_est))
        dsw.boot <- dsw.boot[order(dsw.boot$id, dsw.boot$cal_time), ]
        
        dsw.boot <- dsw.boot[with(dsw.boot, order(id, cal_time, -organ_num)), ]
        dsw.boot <- group_by(dsw.boot, id, cal_time)
        dsw.boot <- mutate(dsw.boot, 
                           prob_tx_received_est = rev(cumprod(prob_tx_received_est)))
        dsw.boot <- as.data.frame(dsw.boot)
        dsw.boot <- dsw.boot[!duplicated(dsw.boot[c("id", "cal_time")]), ]
        
        pt.first.obs.boot <- data.frame(id = unique_ids_1,  cal_time = pt.obs.1.boot$p_arv_times)
        dsw.boot <- merge(dsw.boot, pt.first.obs.boot, all = TRUE)
        dsw.boot <- dsw.boot[order(dsw.boot$id, dsw.boot$cal_time),]
        dsw.boot <- mutate(dsw.boot, 
                           prob_tx_received_est = ifelse(is.na(prob_tx_received_est), 1, prob_tx_received_est))
        
        dsw.boot <- group_by(dsw.boot, id)
        dsw.boot <- mutate(dsw.boot, num_weight = cumprod(prob_tx_received_est)) 
        dsw.boot <- as.data.frame(dsw.boot)
        
        dsw.boot <- select(dsw.boot, cal_time, id, num_weight)
        dsf.final.boot <- mutate(dsf.final.boot , 
                                 cal_time = ifelse(is.na(cal_time) == TRUE, p_arv_times, cal_time))
        dsf.final.boot <- r_join(dsw.boot, dsf.final.boot )
        dsf.final.boot <- mutate(dsf.final.boot, total_weight = num_weight / den_weight)
        
        surv1 <- survfit(Surv(start, stop, event_window) ~ 1, 
                         weights = total_weight,
                         data = dsf.final.boot)
        t1boot <- summary(surv1, times = c(180, 360, 540, 720))
        
        surv3 <- survfit(Surv(start, stop, event_window) ~ 1, 
                         weights = ipw2,
                         data = dsf.final.boot)
        t3boot <- summary(surv3, times = c(180, 360, 540, 720))
        
        surv4 <- survfit(Surv(start, stop, event_window) ~ 1,
                         weights = compliant,
                         data = dsf.final.boot)
        t4boot <- summary(surv4, times = 1:4 * 180)
        
        out <- c(t1boot$surv, t3boot$surv, t4boot$surv, coef(m1.boot))
        boots[b, ] <- out
      }      
    }

    if(include.se) {
      mean.boot <- colMeans(boots, na.rm = TRUE)
      se_est <- apply(boots, 2, sd, na.rm = TRUE)
    } else {
      mean.boot <- rep(NA, ncol(boots))
      se_est <- rep(NA, ncol(boots))
    }

    # end boot ------------------------------------------
    out <- c(days_in_study, sim_num, t1main$surv,
      t3main$surv, t4main$surv, unname(coef(m1)), mean.boot, se_est)
	write.table(t(out), 
      file = outfile, 
      quote = FALSE, 
      row.names = FALSE,
      col.names = FALSE, 
      append = TRUE)
}

cols <- c("days_in_study", "sim_num",
  paste0("allsurv", 1:4 * 180), 
  # paste0("onesurv", 1:4 * 180),
  paste0("one2surv", 1:4 * 180),
  paste0("compsurv", 1:4 * 180),
  c("b0", "b1"),
  paste0("allsurvboot", 1:4 * 180),
  # paste0("onesurvboot", 1:4 * 180),
  paste0("one2survboot", 1:4 * 180),
  paste0("compsurvboot", 1:4 * 180),
  c("b0boot", "b1boot"),
  paste0("allse", 1:4 * 180),
  # paste0("onese", 1:4 * 180),
  paste0("one2se", 1:4 * 180),
  paste0("compse", 1:4 * 180),
  c("b0se", "b1se")
)

outfile <- paste0("./outfiles/sim_out_", j, ".txt")

write.table(t(cols), 
  file = outfile, 
  row.names = FALSE,   
  col.names = FALSE,
  append = FALSE, 
  quote = FALSE)


# run simulation in parallel
sims <- ((j * n.mc / n.s) - (n.mc / n.s - 1)):(j * n.mc / n.s)
out.list <- mclapply(sims,
  failwith(NA, sim),
  mc.cores = n.cores,
  days_in_study = 10* 365,
  include.se = TRUE)
