library(dplyr);library(dlnm);library(gnm);library(splines);library(mixmeta);library(trend)

setwd("E:/OneDrive - University of Maryland/Research/Climate whiplash/")

data <- readRDS("G:/Whiplash US/data/data.final.rds")
# data <- readRDS("data/data.final.2ndmo.rds")

# main analysis 
# define the strata
# data <- data %>%
#   mutate(stratum = factor(paste(GEOID, year, sep=":")))
data.main <- data %>%
  group_by(GEOID) %>%
  filter(sum(total) > 0) %>%
  ungroup()

argvar = list(fun="lin")
arglag = list(fun="ns", df=4)

# 2-stage time-series analysis
whip_vars <- c(
  "whip03_05y_all","whip03_10y_all","whip03_20y_all","whip03_05y_d2w","whip03_10y_d2w","whip03_20y_d2w","whip03_05y_w2d","whip03_10y_w2d","whip03_20y_w2d",
  "whip06_05y_all","whip06_10y_all","whip06_20y_all","whip06_05y_d2w","whip06_10y_d2w","whip06_20y_d2w","whip06_05y_w2d","whip06_10y_w2d","whip06_20y_w2d"
)

main <- list()

for (wv in whip_vars) {
  message("Processing ", wv)
  
  # -------------------------
  # Containers for stage 1
  # -------------------------
  coef_full <- vcov_full <- list()
  
  for(cty in unique(data.main$GEOID)) {
    data.cty <- subset(data.main, GEOID == cty)
    
    # skip if exposure never varies
    if (length(unique(data.cty[[wv]])) == 1) next
    
    # define splines for month of the year
    spltavg <- onebasis(data.cty$tavg, "ns", df=3)
    splprcp <- onebasis(data.cty$prcp, "ns", df=3)
    splmoy <- onebasis(data.cty$month, "ns", df=3)
    
    cb <- crossbasis(data.cty[[wv]], lag=5, argvar=argvar, arglag=arglag) 
    
    mod <- try(glm(total ~ cb + spltavg + splprcp + pm25 + splmoy + factor(year),
                   family=quasipoisson(link="log"), data=data.cty), silent=TRUE)
    if(inherits(mod,"try-error") || any(is.na(coef(mod)))) next
    
    coef_full[[cty]] <- coef(mod)[grep("cb", names(coef(mod)))]
    vcov_full[[cty]] <- vcov(mod)[grep("cb", names(coef(mod))), grep("cb", names(coef(mod)))]
  }
  
  # -------------------------
  # Stack results
  # -------------------------
  coef_mat <- do.call(rbind, coef_full)
  
  meta_df <- cbind(coef_mat, data.main[match(rownames(coef_mat),data.main$GEOID),c("pct_65plus","income","svi","urban")])
  
  meta_df$pct_65plus_c <- scale(meta_df$pct_65plus, scale = FALSE)
  meta_df$income_c     <- scale(meta_df$income, scale = FALSE)
  meta_df$svi_c        <- scale(meta_df$svi, scale = FALSE)
  meta_df$urban <- factor(meta_df$urban, levels = c(0, 1),
                          labels = c("rural", "urban"))
  
  mv_full <- mixmeta(coef_mat ~ pct_65plus_c + income_c + svi_c + urban, data=meta_df, S=vcov_full, method="reml") # convergency issue

  # -------------------------
  # Prediction bases
  # -------------------------
  p_urban <- mean(meta_df$urban == "urban", na.rm = TRUE)
  
  idx_int   <- grep("\\(Intercept\\)", names(coef(mv_full)))
  idx_urban <- grep("urbanurban", names(coef(mv_full)))
  
  beta_avg <- coef(mv_full)[idx_int] + p_urban * coef(mv_full)[idx_urban]
  
  sigma_avg <- vcov(mv_full)[idx_int, idx_int] +
    p_urban^2 * vcov(mv_full)[idx_urban, idx_urban] +
    p_urban * vcov(mv_full)[idx_int, idx_urban] +
    p_urban * vcov(mv_full)[idx_urban, idx_int]
  
  cp <- crosspred(cb, coef=beta_avg, vcov=sigma_avg, model.link="log", 
                  at=0:1, cen=0, lag=5, bylag=0.1)
  
  # -------------------------
  # Extract Cumulative RR for Lag 0–3
  # -------------------------
  lag_basis_b <- do.call("onebasis", c(list(x = 0:5), attr(cb, "arglag")))
  w_b <- colSums(lag_basis_b)
  logRR_b <- as.numeric(w_b %*% beta_avg)
  se_b <- sqrt(as.numeric(t(w_b) %*% sigma_avg %*% w_b))
  
  # -------------------------
  # Save all outputs
  # -------------------------
  main[[wv]] <- list(mvmeta_full = mv_full, 
                     cp_final    = cp,
                     cum_b       = data.frame(estimate = exp(logRR_b), 
                                              se       = se_b,
                                              low      = exp(logRR_b - 1.96 * se_b), 
                                              high     = exp(logRR_b + 1.96 * se_b)))
}
saveRDS(main,"study2/results/main.rds")


# attributable number
for (c in c("03","06")) {
  for (cc in c("all","d2w","w2d")) {
    data.main[, paste0("whip", c, "_", cc)] <- ifelse(data.main[,paste0("whip", c, "_05y_", cc)]==1 & data.main[,paste0("whip", c, "_10y_", cc)]==0 & data.main[,paste0("whip", c, "_20y_", cc)]==0, 1,
                                                      ifelse(data.main[,paste0("whip", c, "_05y_", cc)]==1 & data.main[,paste0("whip", c, "_10y_", cc)]==1 & data.main[,paste0("whip", c, "_20y_", cc)]==0, 2,
                                                             ifelse(data.main[,paste0("whip", c, "_05y_", cc)]==1 & data.main[,paste0("whip", c, "_10y_", cc)]==1 & data.main[,paste0("whip", c, "_20y_", cc)]==1, 3, 0)))
  }
}
data.main$county_year <- paste(data.main$GEOID, data.main$year, sep = "_")
data.main$region <- ifelse(data.main$state %in% c("CT","DL","ME","MD","MA","NH","NJ","NY","PA","RI","VT"),"Northeast",
                           ifelse(data.main$state %in% c("AL","FL","GA","NC","SC","VA"),"Southeast",
                                  ifelse(data.main$state %in% c("IL","IN","KY","MO","OH","TN","WV"),"Ohio Valley",
                                         ifelse(data.main$state %in% c("IA","MI","MN","WI"),"Upper Midwest",
                                                ifelse(data.main$state %in% c("MT","NE","ND","SD","WY"),"Northern Rockies and Plains",
                                                       ifelse(data.main$state %in% c("AR","KS","LA","MS","OK","TX"),"South",
                                                              ifelse(data.main$state %in% c("AZ","CO","NM","UT"),"Southwest",
                                                                     ifelse(data.main$state %in% c("CA","NV"),"West","Northwest"))))))))


set.seed(123)
nsim <- 5000
types <- c("all","d2w","w2d")

an_sim         <- list()
an_county      <- list()
an_year        <- list()
an_county_year <- list()
an_state       <- list()
an_region      <- list()

n_county      <- length(unique(data.main$GEOID))
n_year        <- length(unique(data.main$year))
n_county_year <- length(unique(data.main$county_year))
n_state       <- length(unique(data.main$state))
n_region      <- length(unique(data.main$region))

for (tp in types) {
  an_county[[tp]]      <- matrix(0, nrow=n_county,      ncol=nsim)
  an_year[[tp]]        <- matrix(0, nrow=n_year,        ncol=nsim)
  an_county_year[[tp]] <- matrix(0, nrow=n_county_year, ncol=nsim)
  an_state[[tp]]       <- matrix(0, nrow=n_state,       ncol=nsim)
  an_region[[tp]]      <- matrix(0, nrow=n_region,      ncol=nsim)
}

whip_vars <- c("whip06_05y_all", "whip06_10y_all", "whip06_20y_all", 
               "whip06_05y_d2w", "whip06_10y_d2w", "whip06_20y_d2w", 
               "whip06_05y_w2d", "whip06_10y_w2d", "whip06_20y_w2d")

get_type <- function(wv) {
  if (grepl("all", wv)) return("all")
  if (grepl("d2w", wv)) return("d2w")
  if (grepl("w2d", wv)) return("w2d")
}
get_level <- function(wv) {
  if (grepl("05y", wv)) return(1)
  if (grepl("10y", wv)) return(2)
  if (grepl("20y", wv)) return(3)
}

for (wv in whip_vars) {
  message("Processing ", wv)
  
  level <- get_level(wv)
  type  <- get_type(wv)
  
  # choose multi-level exposure (ONCE)
  var <- switch(type,
                "all" = data.main$whip03_all,
                "d2w" = data.main$whip03_d2w,
                "w2d" = data.main$whip03_w2d
  )
  
  # mutually exclusive indicator (ONCE)
  exposed <- as.numeric(var == level)
  
  coef_full <- vcov_full <- list()
  # an_sim[[wv]] <- matrix(NA, nrow=nrow(data.main), ncol=nsim)
  
  for(cty in unique(data.main$GEOID)) {
    data.cty <- subset(data.main, GEOID == cty)
    
    # skip if exposure never varies
    if (length(unique(data.cty[[wv]])) == 1) next
    
    # define splines for month of the year
    spltavg <- onebasis(data.cty$tavg, "ns", df=3)
    splprcp <- onebasis(data.cty$prcp, "ns", df=3)
    splmoy <- onebasis(data.cty$month, "ns", df=3)
    
    cb <- crossbasis(data.cty[[wv]], lag=5, argvar=argvar, arglag=arglag) 
    
    mod <- try(glm(total ~ cb + spltavg + splprcp + pm25 + splmoy + factor(year),
                   family=quasipoisson(link="log"), data=data.cty), silent=TRUE)
    if(inherits(mod,"try-error") || any(is.na(coef(mod)))) next
    
    coef_full[[cty]] <- coef(mod)[grep("cb", names(coef(mod)))]
    vcov_full[[cty]] <- vcov(mod)[grep("cb", names(coef(mod))), grep("cb", names(coef(mod)))]
  }
  
  coef_mat <- do.call(rbind, coef_full)
  
  meta_df <- cbind(coef_mat, data.main[match(rownames(coef_mat),data.main$GEOID),c("pct_65plus","income","svi","urban")])
  
  meta_df$pct_65plus_c <- scale(meta_df$pct_65plus, scale = FALSE)
  meta_df$income_c     <- scale(meta_df$income, scale = FALSE)
  meta_df$svi_c        <- scale(meta_df$svi, scale = FALSE)
  meta_df$urban <- factor(meta_df$urban, levels = c(0, 1),
                          labels = c("rural", "urban"))
  
  mv_full <- mixmeta(coef_mat~pct_65plus_c+income_c+svi_c+urban, data=meta_df, S=vcov_full, method="reml")
  
  p_urban <- mean(meta_df$urban == "urban", na.rm = TRUE)
  idx_int   <- grep("\\(Intercept\\)", names(coef(mv_full)))
  idx_urban <- grep("urbanurban", names(coef(mv_full)))
  
  beta_avg <- coef(mv_full)[idx_int] + p_urban * coef(mv_full)[idx_urban]
  
  sigma_avg <- vcov(mv_full)[idx_int, idx_int] +
    p_urban^2 * vcov(mv_full)[idx_urban, idx_urban] +
    p_urban * vcov(mv_full)[idx_int, idx_urban] +
    p_urban * vcov(mv_full)[idx_urban, idx_int]
  
  beta_sim <- MASS::mvrnorm(n=nsim, mu=beta_avg, Sigma=sigma_avg)
  lag_basis_03 <- do.call("onebasis", c(list(x = 0:5), attr(cb, "arglag")))
  w_03 <- colSums(lag_basis_03)
  
  for (i in 1:nsim) {
    logRR <- as.numeric(w_03 %*% beta_sim[i, ])
    rr <- exp(logRR)
    af <- (rr - 1) / rr # attributable fraction

    contrib <- (af * exposed) * data.main$total
    
    # Level-specific results (5y vs 10y vs 20y)
    # an_sim[[wv]][ ,i] <- contrib 
    
    # Total burden by type
    an_county[[type]][, i] <-
      an_county[[type]][, i] + rowsum(contrib, data.main$GEOID)[,1] / 21
    
    an_year[[type]][, i] <-
      an_year[[type]][, i] + rowsum(contrib, data.main$year)[,1]
    
    an_county_year[[type]][, i] <-
      an_county_year[[type]][, i] + rowsum(contrib, data.main$county_year)[,1]
    
    an_state[[type]][, i] <-
      an_state[[type]][, i] + rowsum(contrib, data.main$state)[,1] / 21
    
    an_region[[type]][, i] <-
      an_region[[type]][, i] + rowsum(contrib, data.main$region)[,1] / 21
  }
  gc()
}

total.state <- rowsum(data.main$total, data.main$state)/21  # annual average total deaths
total.region <- rowsum(data.main$total, data.main$region)/21

af_state  <- list()
af_region <- list()

an_county_df <- data.frame(matrix(nrow=length(unique(data.main$GEOID)), ncol=20))
an_year_df <- data.frame(matrix(nrow=length(2003:2023), ncol=10))
an_county_year_df <- data.frame(matrix(nrow=length(unique(data.main$GEOID)), ncol=65))
an_state_year_df <- data.frame(matrix(nrow=length(unique(data.main$state))*21, ncol=11))
an_state_df <- af_state_df <- data.frame(matrix(nrow=49, ncol=10))
an_region_df <- af_region_df <- data.frame(matrix(nrow=9, ncol=10))

names(an_county_df) <- c("GEOID","state",
                         paste0(rep(c("mean_","lower_","upper_",
                                      "abs_change_","change_lower_","change_upper_"),times=3), rep(c("all","d2w","w2d"),each=6)))
names(an_year_df) <- c("year",paste0(rep(c("mean_","lower_","upper_"),times=3), rep(c("all","d2w","w2d"),each=3)))
names(an_county_year_df) <- c("GEOID","state",paste0("all_",2003:2023),paste0("d2w_",2003:2023),paste0("w2d_",2003:2023))
names(an_state_year_df) <- c("state","year",paste0(rep(c("all_","d2w_","w2d_"),each=3),c("mean","lower","upper")))
names(an_region_df) <- names(af_region_df) <- 
  c("region",paste0(rep(c("mean_","lower_","upper_"),times=3), rep(c("all","d2w","w2d"),each=3)))
names(an_state_df) <- names(af_state_df) <- 
  c("state",paste0(rep(c("mean_","lower_","upper_"),times=3), rep(c("all","d2w","w2d"),each=3)))

an_county_df$GEOID <- an_county_year_df$GEOID <- unique(data.main$GEOID)
an_county_df$state <- an_county_year_df$state <- data.main[data.main$yearmon=="2003-01",]$state
an_state_year_df$state <- rep(sort(unique(data.main$state)),each=21)
an_year_df$year <- 2003:2023
an_state_year_df$year <- rep(2003:2023,times=49)
an_state_df$state <- af_state_df$state <- sort(unique(data.main$state))
an_region_df$region <- af_region_df$region <- sort(unique(data.main$region))

for (tp in c("all","d2w","w2d")) {
  data.unique <- unique(data.main[, c("GEOID", "state", "year")])
  county_year_df <- data.frame(
    GEOID = data.unique$GEOID,
    state = data.unique$state,
    year  = data.unique$year,
    an_county_year[[tp]]
  )
  
  p1_rows <- county_year_df$year >= 2003 & county_year_df$year <= 2012
  p2_rows <- county_year_df$year >= 2013 & county_year_df$year <= 2023
  
  an_county_0312 <- rowsum(an_county_year[[tp]][p1_rows, ], county_year_df$GEOID[p1_rows])
  an_county_1323 <- rowsum(an_county_year[[tp]][p2_rows, ], county_year_df$GEOID[p2_rows])
  
  an_county_0312_annual <- an_county_0312 / 10
  an_county_1323_annual <- an_county_1323 / 11
  
  abs_change <- an_county_1323_annual - an_county_0312_annual
  
  county_year_long <- data.frame(
    GEOID   = county_year_df$GEOID,
    year    = county_year_df$year,
    an      = rowMeans(county_year_df[,4:5003]) # equal to rowMeans(an_sim_county_year1[[wv]])
  )
  county_year_wide <- reshape(county_year_long, idvar = "GEOID", timevar = "year", direction = "wide")

  tmp <- county_year_df[,2:5003] %>%
    group_by(state, year) %>%
    summarise(across(everything(), sum), .groups = "drop")
  mat <- as.matrix(tmp[, -(1:2)])
  state_year_long <- tmp %>%
    mutate(
      mean  = rowMeans(mat),
      lower = apply(mat, 1, quantile, 0.025),
      upper = apply(mat, 1, quantile, 0.975)
    ) %>%
    dplyr::select(state, year, mean, lower, upper)
  
  tp.idx <- ifelse(tp=="all", 1, ifelse(tp=="d2w", 2, 3))
  
  an_county_df[,tp.idx*6-3] <- rowMeans(an_county[[tp]]) # average annual attributable deaths
  an_county_df[,tp.idx*6-2] <- apply(an_county[[tp]], 1, quantile, 0.025)
  an_county_df[,tp.idx*6-1] <- apply(an_county[[tp]], 1, quantile, 0.975)
  an_county_df[,tp.idx*6]   <- rowMeans(abs_change)
  an_county_df[,tp.idx*6+1] <- apply(abs_change, 1, quantile, 0.025)
  an_county_df[,tp.idx*6+2] <- apply(abs_change, 1, quantile, 0.975)
  
  an_year_df[,tp.idx*3-1] <- rowMeans(an_year[[tp]])
  an_year_df[,tp.idx*3]   <- apply(an_year[[tp]], 1, quantile, 0.025)
  an_year_df[,tp.idx*3+1] <- apply(an_year[[tp]], 1, quantile, 0.975)
  
  an_county_year_df[,(tp.idx*21-18):(tp.idx*21+2)] <- county_year_wide[,2:22]
  an_state_year_df[,(tp.idx*3):(tp.idx*3+2)] <- state_year_long[,3:5]
  
  af_state[[tp]]  <- an_state[[tp]]/as.vector(total.state)*100
  af_region[[tp]] <- an_region[[tp]]/as.vector(total.region)*100
  
  an_state_df[,tp.idx*3-1] <- rowMeans(an_state[[tp]])
  an_state_df[,tp.idx*3]   <- apply(an_state[[tp]], 1, quantile, 0.025)
  an_state_df[,tp.idx*3+1] <- apply(an_state[[tp]], 1, quantile, 0.975)
  
  af_state_df[,tp.idx*3-1] <- rowMeans(af_state[[tp]])
  af_state_df[,tp.idx*3]   <- apply(af_state[[tp]], 1, quantile, 0.025)
  af_state_df[,tp.idx*3+1] <- apply(af_state[[tp]], 1, quantile, 0.975)
  
  an_region_df[,tp.idx*3-1] <- rowMeans(an_region[[tp]])
  an_region_df[,tp.idx*3]   <- apply(an_region[[tp]], 1, quantile, 0.025)
  an_region_df[,tp.idx*3+1] <- apply(an_region[[tp]], 1, quantile, 0.975)
  
  af_region_df[,tp.idx*3-1] <- rowMeans(af_region[[tp]])
  af_region_df[,tp.idx*3]   <- apply(af_region[[tp]], 1, quantile, 0.025)
  af_region_df[,tp.idx*3+1] <- apply(af_region[[tp]], 1, quantile, 0.975)
}

trend_fun <- function(x) {
  if (all(is.na(x))) return(c(NA, NA))
  res <- try(mk.test(x), silent = TRUE)
  c(res$estimates, res$p.value)
}

trend_all <- as.data.frame(t(apply(an_county_year_df[,3:23],1,trend_fun)))
trend_d2w <- as.data.frame(t(apply(an_county_year_df[,24:44],1,trend_fun)))
trend_w2d <- as.data.frame(t(apply(an_county_year_df[,45:65],1,trend_fun)))
names(trend_all) <- names(trend_d2w) <- names(trend_w2d) <- c("S","varS","tau", "pval")
an_county_year_df$trend_all <- trend_all$tau
an_county_year_df$sig_all <- trend_all$pval
an_county_year_df$trend_d2w <- trend_d2w$tau
an_county_year_df$sig_d2w <- trend_d2w$pval
an_county_year_df$trend_w2d <- trend_w2d$tau
an_county_year_df$sig_w2d <- trend_w2d$pval

an_county_df <- cbind(an_county_df, 
                      an_county_year_df$trend_all, an_county_year_df$sig_all, 
                      an_county_year_df$trend_d2w, an_county_year_df$sig_d2w, 
                      an_county_year_df$trend_w2d, an_county_year_df$sig_w2d)
names(an_county_df)[21:26] <- c("trend_all","sig_all","trend_d2w","sig_d2w","trend_w2d","sig_w2d")

write.csv(an_county_df,"study2/results/burden_new2/an_county.csv", row.names=FALSE)
write.csv(an_year_df,"study2/results/burden_new2/an_year.csv", row.names=FALSE)
write.csv(an_county_year_df,"study2/results/burden_new2/an_county_year.csv", row.names=FALSE)
write.csv(an_state_year_df,"study2/results/burden_new2/an_state_year.csv", row.names=FALSE)
write.csv(an_state_df,"study2/results/burden_new2/an_state.csv", row.names=FALSE)
write.csv(an_region_df,"study2/results/burden_new2/an_region.csv", row.names=FALSE)
write.csv(af_state_df,"study2/results/burden_new2/af_state.csv", row.names=FALSE)
write.csv(af_region_df,"study2/results/burden_new2/af_region.csv", row.names=FALSE)



