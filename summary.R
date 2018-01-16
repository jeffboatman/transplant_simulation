library(survival)
library(dplyr)
library(xtable)

options(digits = 3)

# survival times
times <- 1:4 * 180

# true survival if everyone follows the regime.
all.files <- paste0("./outfiles/sim_truth_all_follow_out_", 1:5, ".txt")
tt_all <- do.call(rbind, lapply(all.files, read.table, header = TRUE))
(tmeans_all <- colMeans(tt_all[, paste0("surv", times)]))
sprintf("%.13f", tmeans_all)

one.files <- paste0("./outfiles/sim_truth_one_follows_out_", 1:5, ".txt")
tt_one <- do.call(rbind, lapply(one.files, read.table, header = TRUE))
event.time <- tt_one$ev_pt_time
event.indicator <- tt_one$ev_ind
surv <- survfit(Surv(event.time, event.indicator) ~ 1)
(tmeans_one <- summary(surv, times)$surv)
sprintf("%.13f", tmeans_one)

outfiles <- paste0("./outfiles/sim_out_", 1:5, ".txt")
dfs <- lapply(outfiles, function(file) read.table(file, header = TRUE))
est <- do.call(rbind, dfs)

### bias --------------------------------------------------------------------
(emeans_all <- colMeans(est[, paste0("allsurv", times)]))
(bias_all <- emeans_all - tmeans_all)

emeans_one2 <- colMeans(est[, paste0("one2surv", times)])
(bias_one2 <- emeans_one2 - tmeans_one)


nc_means <- colMeans(est[, paste0("compsurv", times)])

### coverage probability ----------------------------------------------------
z <- qnorm(0.975)
tdf <- data.frame(matrix(c(tmeans_all, tmeans_one), nrow = 1))
names(tdf) <- c(paste0("ta", times), paste0("to", times))
cpdf <- data.frame(est, tdf)

cpdf <- mutate(cpdf,
  cp.all.180 = allsurv180 - z * allse180 < ta180 & ta180 <  allsurv180 + z * allse180,
  cp.all.360 = allsurv360 - z * allse360 < ta360 & ta360 <  allsurv360 + z * allse360,
  cp.all.540 = allsurv540 - z * allse540 < ta540 & ta540 <  allsurv540 + z * allse540,
  cp.all.720 = allsurv720 - z * allse720 < ta720 & ta720 <  allsurv720 + z * allse720,       
  xcp.all.180 = allsurv180 - z * allse180 < to180 & to180 <  allsurv180 + z * allse180,
  xcp.all.360 = allsurv360 - z * allse360 < to360 & to360 <  allsurv360 + z * allse360,
  xcp.all.540 = allsurv540 - z * allse540 < to540 & to540 <  allsurv540 + z * allse540,
  xcp.all.720 = allsurv720 - z * allse720 < to720 & to720 <  allsurv720 + z * allse720
)

cpdf <- mutate(cpdf,
  cp.one.180 = one2surv180 - z * one2se180 < to180 & to180 <  one2surv180 + z * one2se180,
  cp.one.360 = one2surv360 - z * one2se360 < to360 & to360 <  one2surv360 + z * one2se360,
  cp.one.540 = one2surv540 - z * one2se540 < to540 & to540 <  one2surv540 + z * one2se540,
  cp.one.720 = one2surv720 - z * one2se720 < to720 & to720 <  one2surv720 + z * one2se720,
  xcp.one.180 = one2surv180 - z * one2se180 < ta180 & ta180 <  one2surv180 + z * one2se180,
  xcp.one.360 = one2surv360 - z * one2se360 < ta360 & ta360 <  one2surv360 + z * one2se360,
  xcp.one.540 = one2surv540 - z * one2se540 < ta540 & ta540 <  one2surv540 + z * one2se540,
  xcp.one.720 = one2surv720 - z * one2se720 < ta720 & ta720 <  one2surv720 + z * one2se720
)

(cpse_all <- colMeans(cpdf[, paste0("cp.all.", times)]))
(cpse_one <- colMeans(cpdf[, paste0("cp.one.", times)]))

xcpse_all <- colMeans(cpdf[, paste0("xcp.all.", times)])
xcpse_one <- colMeans(cpdf[, paste0("xcp.one.", times)])

# Mean SE / Monte Carlo SE
colSd <- function(x) apply(x, 2, sd)

colMeans(est[, paste0("allse", times)]) / colSd(est[, paste0("allsurv", times)])
colMeans(est[, paste0("one2se", times)]) / colSd(est[, paste0("one2surv", times)])


### Table showing Bias and coverage probabilties
(bias_cp <- cbind(rep(1:4 * 180, 2), 
  c(tmeans_one, tmeans_all),
  c(nc_means - tmeans_one, nc_means - tmeans_all),
  c(emeans_one2 - tmeans_one, emeans_one2 - tmeans_all),
  c(emeans_all - tmeans_one, emeans_all - tmeans_all),
  c(cpse_one, xcpse_one),
  c(xcpse_all, cpse_all))
)

print(xtable(bias_cp, digits = c(0, 0, 3, 3, 3, 3, 3, 3)), include.rownames = FALSE,
  include.colnames = FALSE)
