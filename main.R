## set up

library(kernlab)
library(scoringRules)
library(Rcpp)
library(ggplot2)
library(WeightedForecastVerification)
library(crch)
source("utility_funcs.R")
#sourceCpp("crps_div.cpp")
sourceCpp("kernels.cpp")

## similarity metric between the ensembles (could also add additional covariates)
## discretise the similarity metric and treat the groups separately
## parametrise a weight function that decreases with similarity (needed?)
## may need a cut off so that not everything is considered
## treat everything at all locations together?
## extend to continuous distributions somehow - make the link to existing post-processing methods
## as a comparison, do the discretisation in space-time (or using topological covariates)
## this extends easily to multivariate case, and to other spaces (e.g. functions or graphs)
## analog-based post-processing - read the literature on existing analog-based methods


## initialise variables

# # univariate kernels
# gaussk <- function(x1, x2, sigma = 1) exp(-((x1 - x2)^2)/sigma)
# laplk <- function(x1, x2, sigma = 1) exp(-abs(x1 - x2)/sigma)
# enerk <- function(x1, x2) abs(x1) + abs(x2) - abs(x1 - x2)
# 
# # multivariate kernels
# eucabs <- function(x) scoringRules:::euclnormC(x)
# gaussk_mv <- function(x1, x2, sigma = 1) exp(-(eucabs(x1 - x2)^2)/sigma)
# laplk_mv <- function(x1, x2, sigma = 1) exp(-eucabs(x1 - x2)/sigma)
# enerk <- function(x1, x2) eucabs(x1) + eucabs(x2) - eucabs(x1 - x2)

# lead times
lt_vec <- 1:33

# functions
get_H <- function(x, kernel = "Gaussian") {
  M <- ncol(x)
  if (kernel == "Energy") {
    kern <- enerk_mat
  } else if (kernel == "Laplace") {
    kern <- laplk_mat
  } else {
    kern <- gaussk_mat
  }
  kmat <- kern(x, x)
  # kmat <- matrix(NA, M, M)
  # for (i in 1:M) {
  #   for (j in i:M) {
  #     kmat[i, j] <- kmat[j, i] <- mean(kern(x[, i], x[, j]))
  #   }
  # }
  return(kmat)
}

get_c <- function(x, y, kernel = "Gaussian") {
  M <- ncol(x)
  if (kernel == "Energy") {
    kern <- enerk_mat
  } else if (kernel == "Laplace") {
    kern <- laplk_mat
  } else {
    kern <- gaussk_mat
  }
  -sapply(1:M, function(i) mean(kern(as.matrix(x[, i]), as.matrix(y))))
}

get_mv_c <- function(x, y, kernel = "Gaussian") {
  n <- dim(x)[1]
  M <- dim(x)[3]
  if (kernel == "Energy") {
    kern <- enerk_mv_mat
  } else if (kernel == "Laplace") {
    kern <- laplk_mv_mat
  } else {
    kern <- gaussk_mv_mat
  }
  -sapply(1:M, function(i) mean(sapply(1:n, function(t) kern(as.matrix(x[t, , i]), as.matrix(y[t, ])))))
}

get_mv_H <- function(x, kernel = "Gaussian") {
  n <- dim(x)[1]
  M <- dim(x)[3]
  if (kernel == "Energy") {
    kern <- enerk_mv_mat
  } else if (kernel == "Laplace") {
    kern <- laplk_mv_mat
  } else {
    kern <- gaussk_mv_mat
  }
  kmat <- lapply(1:n, function(i) kern(x[i, , ], x[i, , ]))
  kmat <- Reduce("+", kmat) / n
  # kmat <- matrix(NA, M, M)
  # for (i in 1:M) {
  #   for (j in i:M) {
  #     kmat[i, j] <- kmat[j, i] <- mean(sapply(1:n, function(t) kern(x[t, , i], x[t, , j])))
  #   }
  # }
  return(kmat)
}

get_weights <- function(x, y, kernel = "Gaussian", ind = NULL) {
  
  H <- get_H(x, kernel)
  c <- get_c(x, y, kernel)
  
  if (!is.null(ind)) {
    M <- length(ind)
    H <- sapply(1:M, function(i) sapply(1:M, function(j) mean(H[ind[[i]], ind[[j]]])))
    c <- sapply(1:M, function(i) mean(c[ind[[i]]]))
  } else {
    M <- ncol(x)
  }
  
  r <- 0
  l <- rep(0, M)
  u <- rep(1, M)
  b <- 1
  A <- rep(1, M)
  
  primal(ipop(c, H, A, b, l, u, r))
}

get_mv_weights <- function(x, y, kernel = "Gaussian") {
  
  n <- dim(x)[1]
  d <- dim(x)[2]
  M <- dim(x)[3]
  
  r <- 0
  l <- rep(0, M)
  u <- rep(1, M)
  b <- 1
  A <- rep(1, M)
  
  H <- get_mv_H(x, kernel)
  c <- get_mv_c(x, y, kernel)
  
  primal(ipop(c, H, A, b, l, u, r))
}

get_cal <- function(y, x, w = rep(1/ncol(x), ncol(x)), type = NULL) {
  n <- length(y)
  F_y <- sapply(1:n, function(i) sum(w[x[i, ] <= y[i]]))
  Fm_y <- sapply(1:n, function(i) sum(w[x[i, ] < y[i]]))
  Z_F <- Fm_y + runif(n)*(F_y - Fm_y)
  if (!is.null(type)) {
    if (type == "diV") {
      crps_div(Z_F) 
    } else if (type == "var") {
      var(Z_F)
    } else if (type == "mean") {
      mean(Z_F)
    }
  } else {
    Z_F
  }
}


w_list <- mse_list <- cal_list <- crps_list <- es_list <- vector(mode = "list", length = 3)


################################################################################
# Experiment 0: exploratory analysis

experiment0 <- function(lt = 18) {

  dat <- load_data(lt, stat_id = NULL, path = "C:/Users/sa20i493/Documents/Data/MeteoSwiss/") # load dat
  dat <- dat %>% # restrict attention to stations with complete set of observations
    group_by(nat_abbr) %>%
    filter(n() == 1030) %>% 
    ungroup()
  
  stat_ids <- unique(dat$nat_abbr)
  
  crps_mat <- array(NA, c(1030, length(stat_ids), 4))
  
  for (i in seq_along(stat_ids)) {
    
    loc_dat <- dat %>% filter(nat_abbr == stat_ids[i])
    y <- loc_dat$obs
    
    print(c(lt, i))
    
    # COSMO-1E
    x <- loc_dat %>% select(`COSMO-1E`) %>% as.matrix()
    crps_mat[, i, 1] <- crps_sample(y, x)
    
    # COSMO-2E
    x <- loc_dat %>% select(`COSMO-2E`) %>% as.matrix()
    crps_mat[, i, 2] <- crps_sample(y, x)
    
    # IFS
    x <- loc_dat %>% select(ECMWF_IFS) %>% as.matrix()
    crps_mat[, i, 3] <- crps_sample(y, x)
    
    # Multi-model
    x <- loc_dat %>% select(`COSMO-1E`, `COSMO-2E`, ECMWF_IFS) %>% as.matrix()
    crps_mat[, i, 4] <- crps_sample(y, x)
    
  }
  
  df <- data.frame(s = crps_mat[, , 1:3] %>% apply(., c(1, 3), mean) %>% as.vector(),
                   mth = rep(c("COSMO-1E", "COSMO-2E", "ECMWF_IFS"), each = 1030),
                   date = unique(dat$reftime))
  ggplot(df %>% filter(date > "2022-06-01")) + geom_line(aes(x = date, y = s, col = mth)) + 
    scale_x_datetime(name = "Date", expand = c(0, 0)) +
    scale_y_continuous(name = "CRPS (mm)") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.position = "bottom")
  ggsave("Figures/exp0_crps_vs_date.png", width = 6, height = 4)
  
  crps_mat %>% apply(., 3, mean)
    
  crps_mat_rol <- sapply(1:1030, function(i) {
    if (i < 30) {
      numeric(4)*NA
    } else {
      crps_mat[(i - 29):i, , ] %>% apply(., 3, mean)
    }
  })
  df <- data.frame(s = as.vector(crps_mat_rol), 
                   mth = rep(c("COSMO-1E", "COSMO-2E", "ECMWF_IFS", "Multimodel"), each = 1030),
                   date = unique(dat$reftime))
  ggplot(df) + geom_line(aes(x = date, y = s, col = mth)) + 
    scale_x_datetime(name = "Date", expand = c(0, 0)) +
    scale_y_continuous(name = "CRPS (mm)") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.position = "bottom")
  
}
  
 
################################################################################
# Experiment 1: all stations

# notes: can we also calculate the energy score for the univariate methods? Global yes, local no?

experiment1 <- function(kernel, lt_vec = 1:33) {
  w_mat <- mse_mat <- matrix(NA, nrow = 33, ncol = 83)
  crps_mat <- cal_mat <- matrix(NA, nrow = 33, ncol = 5)
  colnames(crps_mat) <- colnames(cal_mat) <- c("COSMO-1E", "COSMO-2E", "ECMWF_IFS", "Multi", "Weighted")
  
  for (lt in lt_vec) {
    
    print(lt)
    if (lt == 1) start <- Sys.time()
    
    dat <- load_data(lt, stat_id = NULL) # load dat
    dat <- dat %>% # restrict attention to stations with complete set of observations
      group_by(nat_abbr) %>%
      filter(n() == 1030) %>% 
      ungroup()
    tr_dat <- dat %>% filter(reftime < "2022-06-01")
    ts_dat <- dat %>% filter(reftime >= "2022-06-01")
    rm(dat)
    
    y_tr <- tr_dat$obs
    y_ts <- ts_dat$obs
    
    
    # COSMO-1E
    x <- ts_dat %>% select(`COSMO-1E`) %>% as.matrix()
    crps_mat[lt, 1] <- mean(crps_sample(y_ts, x))
    #cal_mat[lt, 1] <- get_cal(y_ts, x)
    
    # COSMO-2E
    x <- ts_dat %>% select(`COSMO-2E`) %>% as.matrix()
    crps_mat[lt, 2] <- mean(crps_sample(y_ts, x))
    #cal_mat[lt, 2] <- get_cal(y_ts, x)
    
    # IFS
    x <- ts_dat %>% select(ECMWF_IFS) %>% as.matrix()
    crps_mat[lt, 3] <- mean(crps_sample(y_ts, x))
    #cal_mat[lt, 3] <- get_cal(y_ts, x)
    
    # Multi-model
    x <- ts_dat %>% select(`COSMO-1E`, `COSMO-2E`, ECMWF_IFS) %>% as.matrix()
    crps_mat[lt, 4] <- mean(crps_sample(y_ts, x))
    #cal_mat[lt, 4] <- get_cal(y_ts, x)
    
    # Weighted
    x_tr <- tr_dat %>% select(`COSMO-1E`, `COSMO-2E`, ECMWF_IFS) %>% as.matrix()
    w <- get_weights(x_tr, y_tr, kernel)
    crps_mat[lt, 5] <- mean(crps_sample(y_ts, x, w = t(replicate(length(y_ts), w))))
    #cal_mat[lt, 5] <- get_cal(y_ts, x, w = t(replicate(length(y_ts), w)))
    
    w_mat[lt, ] <- w
    mse_mat[lt, ] <- sapply(1:ncol(x), function(i) mean(abs(x[, i] - y_ts)^2))
    
    if (lt == 1) print(Sys.time() - start)
  }
  
  return(list(crps = crps_mat, mse = mse_mat, w = w_mat))
}

out <- experiment1(kernel = "Energy")

# save results in a list
crps_list[[1]] <- out$crps
w_list[[1]] <- out$w
mse_list[[1]] <- out$mse

# get pit values
get_pit <- function(lt, w_mat) {
  dat <- load_data(lt, stat_id = NULL) # load dat
  dat <- dat %>% # restrict attention to stations with complete set of observations
    group_by(nat_abbr) %>%
    filter(n() == 1030)
  y <- dat$obs
  
  pit_df <- data.frame(obs = y)
  
  # COSMO-1E
  x <- dat %>% select(`COSMO-1E`) %>% as.matrix()
  pit_df$C1 <- get_cal(y, x)
  
  # COSMO-2E
  x <- dat %>% select(`COSMO-2E`) %>% as.matrix()
  pit_df$C2 <- get_cal(y, x)
  
  # ECMWF_IFS
  x <- dat %>% select(ECMWF_IFS) %>% as.matrix()
  pit_df$IFS <- get_cal(y, x)
  
  # Multimodel
  x <- dat %>% select(c(`COSMO-1E`, `COSMO-2E`, ECMWF_IFS)) %>% as.matrix()
  pit_df$Multi <- get_cal(y, x)
  
  # Weighted
  pit_df$Weighted <- get_cal(y, x, w = w_mat[lt, ])
  
  return(pit_df)
}
get_pit(18, w_list[[1]])





################################################################################
# Experiment 2: each station individually, with climatology added

experiment2 <- function(kernel, lt_vec = 1:33) {
  w_mat <- mse_mat <- matrix(0, nrow = 33, ncol = 83)
  w_mat_cl <- matrix(0, nrow = 33, ncol = 752)
  crps_mat <- cal_mat <- matrix(0, nrow = 33, ncol = 7)
  colnames(crps_mat) <- colnames(cal_mat) <- c("COSMO-1E", "COSMO-2E", "ECMWF_IFS", "Climatology", "Multimodel", "Weighted", "Weighted+Clim")

  for (lt in lt_vec) {
    
    if (lt == 1) start <- Sys.time()
    
    dat <- load_data(lt, stat_id = NULL) # load dat
    dat <- dat %>% # restrict attention to stations with complete set of observations
      group_by(nat_abbr) %>%
      filter(n() == 1030) %>% 
      ungroup()
    
    stat_ids <- unique(dat$nat_abbr)
    for (i in seq_along(stat_ids)) {
      
      tr_dat <- dat %>% filter(nat_abbr == stat_ids[i], reftime < "2022-06-01")
      ts_dat <- dat %>% filter(nat_abbr == stat_ids[i], reftime >= "2022-06-01")
      y_tr <- tr_dat$obs
      y_ts <- ts_dat$obs
      
      print(c(lt, i))
      
      # COSMO-1E
      x <- ts_dat %>% select(`COSMO-1E`) %>% as.matrix()
      crps_mat[lt, 1] <- crps_mat[lt, 1] + sum(crps_sample(y_ts, x))
      #cal_mat[lt, 1] <- get_cal(y, x)
      
      # COSMO-2E
      x <- ts_dat %>% select(`COSMO-2E`) %>% as.matrix()
      crps_mat[lt, 2] <- crps_mat[lt, 2] + sum(crps_sample(y_ts, x))
      #cal_mat[lt, 2] <- get_cal(y, x)
      
      # IFS
      x <- ts_dat %>% select(ECMWF_IFS) %>% as.matrix()
      crps_mat[lt, 3] <- crps_mat[lt, 3] + sum(crps_sample(y_ts, x))
      #cal_mat[lt, 3] <- get_cal(y, x)
      
      # Climatology
      x <- t(replicate(length(y_ts), y_tr))
      crps_mat[lt, 4] <- crps_mat[lt, 4] + sum(crps_sample(y_ts, x))
      #cal_mat[lt, 4] <- get_cal(y, x)
      
      # Multi-model
      x <- ts_dat %>% select(`COSMO-1E`, `COSMO-2E`, ECMWF_IFS) %>% as.matrix()
      crps_mat[lt, 5] <- crps_mat[lt, 5] + sum(crps_sample(y_ts, x))
      #cal_mat[lt, 5] <- get_cal(y, x)
      
      mse_mat[lt, ] <- mse_mat[lt, ] + sapply(1:ncol(x), function(j) sum((y_ts - x[, j])^2))
      
      # Weighted
      x_tr <- tr_dat %>% select(`COSMO-1E`, `COSMO-2E`, ECMWF_IFS) %>% as.matrix()
      w <- get_weights(x_tr, y_tr, kernel)
      crps_mat[lt, 6] <- crps_mat[lt, 6] + sum(crps_sample(y_ts, x, w = t(replicate(length(y_ts), w))))
      #cal_mat[lt, 6] <- get_cal(y_ts, x, w = t(replicate(length(y_ts), w)))
      w_mat[lt, ] <- w_mat[lt, ]  + w
      
      # Weighted+Clim
      x_tr <- cbind(x_tr, t(replicate(length(y_tr), y_tr)))
      w <- get_weights(x_tr, y_tr, kernel)
      x <- cbind(x, t(replicate(length(y_ts), y_tr)))
      crps_mat[lt, 7] <- crps_mat[lt, 7] + sum(crps_sample(y_ts, x, w = t(replicate(length(y_ts), w))))
      #cal_mat[lt, 7] <- get_cal(y, x, w = t(replicate(length(y), w)))
      
      #w <- c(w[1:83], sum(w[84:length(w)]))
      w_mat_cl[lt, ] <- w_mat_cl[lt, ]  + w
      
    }
    
    crps_mat[lt, ] <- crps_mat[lt, ]/sum(dat$reftime >= "2022-06-01")
    mse_mat[lt, ] <- mse_mat[lt, ]/sum(dat$reftime >= "2022-06-01")
    w_mat[lt, ] <- w_mat[lt, ]/length(stat_ids)
    w_mat_cl[lt, ] <- w_mat_cl[lt, ]/length(stat_ids)
    
    if (lt == 1) print(Sys.time() - start)
  }
  
  return(list(crps = crps_mat, mse = mse_list, w = list(w = w_mat, w_cl = w_mat_cl)))
}

out <- experiment2(kernel = "Energy")

# save results in a list
crps_list[[2]] <- out$crps
w_list[[2]] <- out$w
mse_list[[2]] <- out$mse

# get pit values
get_pit <- function(lt, stat_id, kernel = "Energy") {
  
  dat <- load_data(lt, stat_id = stat_id) # load dat
  dat <- dat %>% # restrict attention to stations with complete set of observations
    group_by(nat_abbr) %>%
    filter(n() == 1030) %>%
    ungroup()
  
  tr_dat <- dat %>% filter(reftime < "2022-06-01")
  ts_dat <- dat %>% filter(reftime >= "2022-06-01")
  
  y_tr <- tr_dat$obs
  y <- ts_dat$obs

  pit_df <- data.frame(obs = y)
  
  # COSMO-1E
  x <- ts_dat %>% select(`COSMO-1E`) %>% as.matrix()
  pit_df$C1 <- get_cal(y, x)
    
  # COSMO-2E
  x <- ts_dat %>% select(`COSMO-2E`) %>% as.matrix()
  pit_df$C2 <- get_cal(y, x)
  
  # ECMWF_IFS
  x <- ts_dat %>% select(ECMWF_IFS) %>% as.matrix()
  pit_df$IFS <- get_cal(y, x)
  
  # Climatology
  x <- t(replicate(length(y), y_tr))
  pit_df$Clim <- get_cal(y, x)
  
  # Multimodel
  x <- ts_dat %>% select(c(`COSMO-1E`, `COSMO-2E`, ECMWF_IFS)) %>% as.matrix()
  pit_df$Multi <- get_cal(y, x)
  
  # Weighted
  x_tr <- tr_dat %>% select(c(`COSMO-1E`, `COSMO-2E`, ECMWF_IFS)) %>% as.matrix()
  w <- get_weights(x_tr, y_tr, kernel)
  pit_df$Weighted <- get_cal(y, x, w = w)
    
  # Weighted+Clim
  x_tr <- cbind(x_tr, t(replicate(length(y_tr), y_tr)))
  w <- get_weights(x_tr, y_tr, kernel)
  x <- cbind(x, t(replicate(length(y), y_tr)))
  pit_df$WeightedClim <- get_cal(y, x, w = w)
  
  return(pit_df)
}
pit_df <- get_pit(18, "PAY")


################################################################################
# Experiment 2: get local weights

experiment2_getlocalw <- function(kernel, lt = 18) {
  start <- Sys.time()
    
  dat <- load_data(lt, stat_id = NULL, path = paste0(getwd(), "/")) # load dat
  dat <- dat %>% # restrict attention to stations with complete set of observations
    group_by(nat_abbr) %>%
    filter(n() == 1030) %>% 
    ungroup()
  
  SMN_stations <- readRDS("C:/Users/sa20i493/Documents/R/WeightedLoss/WindGusts/SMN_stations.rds")
  
  stat_ids <- unique(dat$nat_abbr)
  w_df <- data.frame(COSMO.1E = numeric(length(stat_ids)),
                     COSMO.2E = numeric(length(stat_ids)),
                     ECMWF_IFS = numeric(length(stat_ids)),
                     longitude = numeric(length(stat_ids)),
                     latitude = numeric(length(stat_ids)),
                     height = numeric(length(stat_ids)))
  
  for (i in seq_along(stat_ids)) {
    
    print(i)
    
    if (stat_ids[i] %in% SMN_stations$nat_abbr) {
      w_df[c("longitude", "latitude", "height")][i, ] <- 
        SMN_stations %>% filter(nat_abbr == stat_ids[i]) %>% select(longitude, latitude, station_height)
    } else {
      w_df[c("longitude", "latitude")][i, ] <- NA
    }
    
    
    tr_dat <- dat %>% filter(nat_abbr == stat_ids[i], reftime < "2022-06-01")
    ts_dat <- dat %>% filter(nat_abbr == stat_ids[i], reftime >= "2022-06-01")
    y_tr <- tr_dat$obs
    y_ts <- ts_dat$obs

    x <- ts_dat %>% select(`COSMO-1E`, `COSMO-2E`, ECMWF_IFS) %>% as.matrix()
    x_tr <- tr_dat %>% select(`COSMO-1E`, `COSMO-2E`, ECMWF_IFS) %>% as.matrix()
    w <- get_weights(x_tr, y_tr, kernel)
    w_df$COSMO.1E[i] <- sum(w[1:11])
    w_df$COSMO.2E[i] <- sum(w[12:32])
    w_df$ECMWF_IFS[i] <- sum(w[33:83])
  }

  print(Sys.time() - start)
  
  return(w_df)
}

localw_df <- experiment2_getlocalw(kernel = "Energy")

plot_map <- function(z, la, lo){
  SW_map <- map_data(map = "world", region = c("Switzerland"))
  ggplot() + geom_polygon(data = SW_map, aes(x = long, y = lat, group = group),
                          col = "black", fill = "white") +
    geom_point(data = data.frame(Lat = la, Lon = lo, z = z),
               aes(x = Lon, y = Lat, fill = z), shape = 21, size = 2) +
    scale_fill_gradient(low = "yellow", high = "red", name = "") +
    xlab("Longitude") + ylab("Latitude") + theme_void()
}
plot_map(localw_df$COSMO.1E, localw_df$latitude, localw_df$longitude)
plot_map(localw_df$COSMO.2E, localw_df$latitude, localw_df$longitude)
plot_map(localw_df$ECMWF_IFS, localw_df$latitude, localw_df$longitude)

plot(localw_df)


################################################################################
# Experiment 3: multivariate

experiment3 <- function(kernel, lt_vec = 1:33) {
  w_mat <- mse_mat <- matrix(0, nrow = 33, ncol = 83)
  w_mat_cl <- matrix(0, nrow = 33, ncol = 752)
  es_mat <- crps_mat <- matrix(0, nrow = 33, ncol = 7)
  colnames(es_mat) <- colnames(crps_mat) <- c("COSMO-1E", "COSMO-2E", "ECMWF_IFS", "Climatology", "Multimodel", "Weighted", "Weighted+Clim")
  
  for (lt in lt_vec) {

    print(lt)
    
    if (lt == 1) start <- Sys.time()
    
    dat <- load_data(lt, stat_id = NULL) # load dat
    dat <- dat %>% # restrict attention to stations with complete set of observations
      group_by(nat_abbr) %>%
      filter(n() == 1030) %>%
      ungroup()
    
    tr_dat <- dat %>% filter(reftime < "2022-06-01")
    ts_dat <- dat %>% filter(reftime >= "2022-06-01")
    rm(dat)

    y_tr <- tr_dat %>% select(nat_abbr, obs) %>% 
      group_by(nat_abbr) %>%
      mutate(row = row_number()) %>%
      pivot_wider(names_from = nat_abbr, values_from = obs) %>%
      select(-row) %>% as.matrix()
    y_ts <- ts_dat %>% select(nat_abbr, obs) %>% 
      group_by(nat_abbr) %>%
      mutate(row = row_number()) %>%
      pivot_wider(names_from = nat_abbr, values_from = obs) %>%
      select(-row) %>% as.matrix()
    n_tr <- nrow(y_tr)
    n_ts <- nrow(y_ts)
    d <- ncol(y_tr)
    
    stat_ids <- unique(ts_dat$nat_abbr)
    
    # COSMO-1E
    x <- lapply(stat_ids, function(z) {
      ts_dat %>% filter(nat_abbr == z) %>% 
        select(`COSMO-1E`) %>%
        cbind() %>% as.matrix()
    }) %>% simplify2array() %>% aperm(c(1, 3, 2))
    es_mat[lt, 1] <- mean(sapply(1:n_ts, function(i) es_sample(y_ts[i, ], x[i, , ])))
    crps_mat[lt, 1] <- mean(sapply(1:d, function(j) crps_sample(y_ts[, j], x[, j, ])))
    
    # COSMO-2E
    x <- lapply(stat_ids, function(z) {
      ts_dat %>% filter(nat_abbr == z) %>% 
        select(`COSMO-2E`) %>%
        cbind() %>% as.matrix()
    }) %>% simplify2array() %>% aperm(c(1, 3, 2))
    es_mat[lt, 2] <- mean(sapply(1:n_ts, function(i) es_sample(y_ts[i, ], x[i, , ])))
    crps_mat[lt, 2] <- mean(sapply(1:d, function(j) crps_sample(y_ts[, j], x[, j, ])))
    
    # IFS
    x <- lapply(stat_ids, function(z) {
      ts_dat %>% filter(nat_abbr == z) %>% 
        select(ECMWF_IFS) %>%
        cbind() %>% as.matrix()
    }) %>% simplify2array() %>% aperm(c(1, 3, 2))
    es_mat[lt, 3] <- mean(sapply(1:n_ts, function(i) es_sample(y_ts[i, ], x[i, , ])))
    crps_mat[lt, 3] <- mean(sapply(1:d, function(j) crps_sample(y_ts[, j], x[, j, ])))
    
    # Climatology
    y_mat <- replicate(n_ts, t(y_tr)) %>% aperm(c(3, 1, 2))
    es_mat[lt, 4] <- mean(sapply(1:n_ts, function(i) es_sample(y_ts[i, ], t(y_tr))))
    crps_mat[lt, 4] <- mean(sapply(1:d, function(j) crps_sample(y_ts[, j], y_mat[, j, ])))
    
    # Multi-model
    x <- lapply(stat_ids, function(z) {
      ts_dat %>% filter(nat_abbr == z) %>% 
        select(`COSMO-1E`, `COSMO-2E`, ECMWF_IFS) %>%
        cbind() %>% as.matrix()
    }) %>% simplify2array() %>% aperm(c(1, 3, 2))
    es_mat[lt, 5] <- mean(sapply(1:n_ts, function(i) es_sample(y_ts[i, ], x[i, , ])))
    crps_mat[lt, 5] <- mean(sapply(1:d, function(j) crps_sample(y_ts[, j], x[, j, ])))
    
    mse_mat[lt, ] <- sapply(1:dim(x)[3], function(j) mean(sapply(1:n_ts, function(i) euclnormC(x[i, , j] - y_ts[i, ])^2)))
    
    # Weighted
    x_tr <- lapply(stat_ids, function(z) {
      tr_dat %>% filter(nat_abbr == z) %>% 
        select(`COSMO-1E`, `COSMO-2E`, ECMWF_IFS) %>%
        cbind() %>% as.matrix()
    }) %>% simplify2array() %>% aperm(c(1, 3, 2))
    w <- get_mv_weights(x_tr, y_tr, kernel)
    es_mat[lt, 6] <- mean(sapply(1:n_ts, function(i) es_sample(y_ts[i, ], x[i, , ], w = w)))
    crps_mat[lt, 6] <- mean(sapply(1:d, function(j) crps_sample(y_ts[, j], x[, j, ], w = t(replicate(n_ts, w)))))
    w_mat[lt, ] <- w
    
    # Weighted+Clim
    y_mat_tr <- replicate(n_tr, t(y_tr)) %>% aperm(c(3, 1, 2))
    x_tr <- abind::abind(x_tr, y_mat_tr, along = 3)
    w <- get_mv_weights(x_tr, y_tr, kernel)
    x <- abind::abind(x, y_mat, along = 3)
    es_mat[lt, 7] <- mean(sapply(1:n_ts, function(i) es_sample(y_ts[i, ], x[i, , ], w = w)))
    crps_mat[lt, 7] <- mean(sapply(1:d, function(j) crps_sample(y_ts[, j], x[, j, ], w = t(replicate(n_ts, w)))))
    w_mat_cl[lt, ] <- w
    
    if (lt == 1) print(Sys.time() - start)
    
  }
  
  # library(foreach)
  # library(doParallel)
  # library(doSNOW)
  # cl <- makeCluster(6)
  # registerDoSNOW(cl)
  # pb <- txtProgressBar(min = 1, max = max(lt_vec), style = 3)
  # progress <- function(n) setTxtProgressBar(pb, n)
  # opts <- list(progress = progress)
  # output <- foreach(lt = 1:max(lt_vec), .packages = c("Rcpp", "scoringRules", "kernlab"), .options.snow = opts) %dopar% {
  #   source("example_usage.R")
  #   sourceCpp("kernels.cpp")
  #   get_out_lt(lt)
  # }
  # close(pb)
  # stopCluster(cl)
  
  return(list(es = es_mat, crps = crps_mat, w = list(w = w_mat, w_cl = w_mat_cl), mse = mse_mat))
}

out <- experiment3(kernel = "Energy")

# save results in a list
crps_list[[3]] <- out$crps
es_list[[3]] <- out$es
w_list[[3]] <- out$w
mse_list[[3]] <- out$mse


################################################################################
# Experiment 3: add baseline multivariate methods

# function to apply ensemble copula coupling
apply_template <- function(ref_ens, rank_template = NULL) {
  ens <- array(NA, dim(ref_ens))
  n_test <- dim(ens)[1]
  n_locs <- dim(ens)[2]
  n_ens <- dim(ens)[3]
  if (!is.null(rank_template)) {
    for(i in 1:n_test){
      for(j in 1:n_locs){
        rank_ord <- rank_template[i, j, ]
        ens[i, j, ] <- ref_ens[i, j, rank_ord]
      }
    }
  } else {
    for(i in 1:n_test){
      for(j in 1:n_locs){
        ens[i, j, ] <- ref_ens[i, j, sample(1:n_ens)]
      }
    }
  }
  return(ens)
}

experiment_baseline <- function(lt_vec = 1:33) {
  es_mat <- crps_mat <- matrix(0, nrow = 33, ncol = 2)
  colnames(es_mat) <- colnames(crps_mat) <- c("ECC", "SS")
  
  for (lt in lt_vec) {
    
    print(lt)
    
    start <- Sys.time()
    
    dat <- load_data(lt, stat_id = NULL) # load dat
    dat <- dat %>% # restrict attention to stations with complete set of observations
      group_by(nat_abbr) %>%
      filter(n() == 1030) %>%
      ungroup()
    
    stat_ids <- unique(dat$nat_abbr)
    
    mvensemble <- array(NA, c(361, length(stat_ids), 83))
    
    for (i in seq_along(stat_ids)) {
      
      tr_dat <- dat %>% filter(nat_abbr == stat_ids[i], reftime < "2022-06-01")
      ts_dat <- dat %>% filter(nat_abbr == stat_ids[i], reftime >= "2022-06-01")
      y_tr <- tr_dat$obs
      y_ts <- ts_dat$obs
      
      # Univariate EMOS
      x_tr <- tr_dat %>% select(`COSMO-1E`, `COSMO-2E`, ECMWF_IFS) %>% as.matrix()
      ensmn <- rowMeans(x_tr)
      enssd <- apply(x_tr, 1, sd)
      model <- crch(y_tr ~ ensmn|enssd, dist = "logistic", left = 0, truncated = TRUE, type = "crps")
      x_ts <- ts_dat %>% select(`COSMO-1E`, `COSMO-2E`, ECMWF_IFS) %>% as.matrix()
      ensmn <- rowMeans(x_ts)
      enssd <- apply(x_ts, 1, sd)
      mu <- predict(model, newdata = data.frame(ensmn = ensmn, enssd = enssd), type = "location")
      sig <- predict(model, newdata = data.frame(ensmn = ensmn, enssd = enssd), type = "scale")
      crps_mat[lt, ] <- crps_mat[lt, ] + sum(crps_tlogis(y_ts, mu, sig, lower = 0))
      
      # Sample ensemble
      mvensemble[, i, ] <- sapply(1:length(y_ts), function(t) qtlogis((1:83)/84, mu[t], sig[t], left = 0)) %>% t()
    }
    
    crps_mat[lt, ] <- crps_mat[lt, ] / sum(dat$reftime >= "2022-06-01")

    
    # mv data
    y_mv <- dat %>% filter(reftime >= "2022-06-01") %>%
      select(nat_abbr, obs) %>% 
      group_by(nat_abbr) %>%
      mutate(row = row_number()) %>%
      pivot_wider(names_from = nat_abbr, values_from = obs) %>%
      select(-row) %>% as.matrix()
    
    y_mv_tr <- dat %>% filter(reftime < "2022-06-01") %>%
      select(nat_abbr, obs) %>% 
      group_by(nat_abbr) %>%
      mutate(row = row_number()) %>%
      pivot_wider(names_from = nat_abbr, values_from = obs) %>%
      select(-row) %>% as.matrix()
    
    # ECC
    x <- lapply(unique(dat$nat_abbr), function(z) {
      dat %>% filter(nat_abbr == z, reftime >= "2022-06-01") %>% 
        select(`COSMO-1E`, `COSMO-2E`, ECMWF_IFS) %>%
        cbind() %>% as.matrix()
    }) %>% simplify2array() %>% aperm(c(1, 3, 2))
    
    rank_mat <- apply(x, c(1, 2), rank, ties.method = "random") %>% 
      aperm(c(2, 3, 1)) %>% 
      unname() # ensemble dependence template
    x_ecc <- apply_template(mvensemble, rank_mat) # apply ECC
    es_mat[lt, 1] <- mean(sapply(1:nrow(y_mv), function(t) es_sample(y_mv[t, ], x_ecc[t, , ])))
    
    # SS
    x <- replicate(nrow(y_mv), y_mv_tr[sample(1:nrow(y_mv_tr), size = 83), ]) %>% aperm(c(3, 2, 1))
    rank_mat <- apply(x, c(1, 2), rank, ties.method = "random") %>% 
      aperm(c(2, 3, 1)) %>% 
      unname() # (random) Schaake shuffle dependence template
    x_ss <- apply_template(mvensemble, rank_mat) # apply ECC
    es_mat[lt, 2] <- mean(sapply(1:nrow(y_mv), function(t) es_sample(y_mv[t, ], x_ss[t, , ])))
  
    print(Sys.time() - start)  
    
  }
  
  return(list(es = es_mat, crps = crps_mat))
}

out <- experiment_baseline()

crps_list[[2]] <- cbind(crps_list[[2]], out$crps)
es_list[[3]] <- cbind(es_list[[3]], out$es)

#save(crps_list, es_list, w_list, mse_list, file = "all_experiments.RData")


################################################################################
# Experiment 4: shared weights for ensemble systems

experiment4 <- function(kernel, lt = 18) {
  start <- Sys.time()
  
  dat <- load_data(lt, stat_id = NULL, path = "C:/Users/sa20i493/Documents/Data/MeteoSwiss/") # load dat
  dat <- dat %>% # restrict attention to stations with complete set of observations
    group_by(nat_abbr) %>%
    filter(n() == 1030) %>% 
    ungroup()
  
  SMN_stations <- readRDS("C:/Users/sa20i493/Documents/R/WeightedLoss/WindGusts/SMN_stations.rds")
  
  stat_ids <- unique(dat$nat_abbr)
  w_df <- data.frame(COSMO.1E = numeric(length(stat_ids)),
                     COSMO.2E = numeric(length(stat_ids)),
                     ECMWF_IFS = numeric(length(stat_ids)),
                     longitude = numeric(length(stat_ids)),
                     latitude = numeric(length(stat_ids)),
                     height = numeric(length(stat_ids)))
  
  for (i in seq_along(stat_ids)) {
    
    print(i)
    
    if (stat_ids[i] %in% SMN_stations$nat_abbr) {
      w_df[c("longitude", "latitude", "height")][i, ] <- 
        SMN_stations %>% filter(nat_abbr == stat_ids[i]) %>% select(longitude, latitude, station_height)
    } else {
      w_df[c("longitude", "latitude")][i, ] <- NA
    }
    
    
    tr_dat <- dat %>% filter(nat_abbr == stat_ids[i], reftime < "2022-06-01")
    ts_dat <- dat %>% filter(nat_abbr == stat_ids[i], reftime >= "2022-06-01")
    y_tr <- tr_dat$obs
    y_ts <- ts_dat$obs
    
    x <- ts_dat %>% select(`COSMO-1E`, `COSMO-2E`, ECMWF_IFS) %>% as.matrix()
    x_tr <- tr_dat %>% select(`COSMO-1E`, `COSMO-2E`, ECMWF_IFS) %>% as.matrix()
    w <- get_weights(x_tr, y_tr, kernel, ind = list(1:11, 12:32, 33:83))
    w_df$COSMO.1E[i] <- w[1]
    w_df$COSMO.2E[i] <- w[2]
    w_df$ECMWF_IFS[i] <- w[3]
  }
  
  print(Sys.time() - start)
  
  return(w_df)
}

shared_w <- experiment4(kernel = "Energy")

plot_map(shared_w$COSMO.1E, shared_w$latitude, shared_w$longitude)
plot_map(shared_w$COSMO.2E, shared_w$latitude, shared_w$longitude)
plot_map(shared_w$ECMWF_IFS, shared_w$latitude, shared_w$longitude)

plot(shared_w)


################################################################################
# Experiment 5: shared weights with rolling window

experiment4_roll <- function(kernel, lt = 18, win_len = 30) {
  start <- Sys.time()
  
  dat <- load_data(lt, stat_id = NULL, path = "C:/Users/sa20i493/Documents/Data/MeteoSwiss/") # load dat
  dat <- dat %>% # restrict attention to stations with complete set of observations
    group_by(nat_abbr) %>%
    filter(n() == 1030) %>% 
    ungroup()
  
  SMN_stations <- readRDS("C:/Users/sa20i493/Documents/R/WeightedLoss/WindGusts/SMN_stations.rds")
  
  stat_ids <- unique(dat$nat_abbr)
  w_list <- vector("list", length(stat_ids))

  date_list <<- unique(dat$reftime)

  for (i in seq_along(stat_ids)) {
    
    w_df <- data.frame(COSMO.1E = numeric(length(date_list)),
                       COSMO.2E = numeric(length(date_list)),
                       ECMWF_IFS = numeric(length(date_list)))
    
    if (stat_ids[i] %in% SMN_stations$nat_abbr) {
      w_list[[i]]$lon <- SMN_stations %>% filter(nat_abbr == stat_ids[i]) %>% select(longitude)
      w_list[[i]]$lat <- SMN_stations %>% filter(nat_abbr == stat_ids[i]) %>% select(latitude)
      w_list[[i]]$hei <- SMN_stations %>% filter(nat_abbr == stat_ids[i]) %>% select(station_height)
    } else {
      w_list[[i]]$lon <- NA
      w_list[[i]]$lat <- NA
      w_list[[i]]$hei <- NA
    }
    
    for (j in seq_along(date_list)) {
      
      print(c(i, j))
      
      ele <- date_list[j]
      
      tr_dates <- tail(sort(date_list[date_list < ele]), win_len)
      if (sum(tr_dates %in% date_list) < win_len) {
        next
      }

      tr_dat <- dat %>% filter(nat_abbr == stat_ids[i], reftime %in% tr_dates)
      ts_dat <- dat %>% filter(nat_abbr == stat_ids[i], reftime == ele)
      y_tr <- tr_dat$obs
      y_ts <- ts_dat$obs
      
      x <- ts_dat %>% select(`COSMO-1E`, `COSMO-2E`, ECMWF_IFS) %>% as.matrix()
      x_tr <- tr_dat %>% select(`COSMO-1E`, `COSMO-2E`, ECMWF_IFS) %>% as.matrix()
      w <- get_weights(x_tr, y_tr, kernel, ind = list(1:11, 12:32, 33:83))

      w_df[j, ] <- w
      
    }

    w_list[[i]]$w <- w_df
    
  }
  
  print(Sys.time() - start)
  
  return(w_list)
}

shared_w_roll <- experiment4_roll(kernel = "Energy")

average_w_roll <- lapply(shared_w_roll, function(x) x$w) %>% Reduce("+", .)/length(shared_w_roll)
average_w_roll[rowSums(average_w_roll) == 0, ] <- NA
average_w_roll$date <- date_list
df <- pivot_longer(average_w_roll, cols = c("COSMO.1E", "COSMO.2E", "ECMWF_IFS"))
ggplot(df) + geom_line(aes(x = date, y = value, col = name)) + 
  scale_x_datetime(name = "Date") +
  scale_y_continuous(name = "Weight", limits = c(0, 0.7), expand = c(0, 0)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom")
ggsave("Figures/exp5_weight_vs_date.png", width = 6, height = 4)


################################################################################
# Experiment 6: shared weights with seasonal window

experiment6_seas <- function(kernel, lt_vec = 1:33) {
 
  SMN_stations <- readRDS("C:/Users/sa20i493/Documents/R/WeightedLoss/WindGusts/SMN_stations.rds")
  
  month_ind <- list(c(12, 1, 2), c(3, 4, 5), c(6, 7, 8), c(9, 10, 11))

  w_list <- vector("list", length(lt_vec))
  
  for (lt in lt_vec) {
    
    if (lt == 1) start <- Sys.time()
    
    dat <- load_data(lt, stat_id = NULL, path = "C:/Users/sa20i493/Documents/Data/MeteoSwiss/") # load dat
    dat <- dat %>% # restrict attention to stations with complete set of observations
      group_by(nat_abbr) %>%
      filter(n() == 1030) %>% 
      ungroup()
    
    stat_ids <- unique(dat$nat_abbr)
    
    w_list[[lt]] <- vector("list", length(stat_ids))
    
    for (i in seq_along(stat_ids)) {

      w_df <- data.frame(COSMO.1E = numeric(4),
                         COSMO.2E = numeric(4),
                         ECMWF_IFS = numeric(4))
      
      if (stat_ids[i] %in% SMN_stations$nat_abbr) {
        w_list[[lt]][[i]]$lon <- SMN_stations %>% filter(nat_abbr == stat_ids[i]) %>% select(longitude)
        w_list[[lt]][[i]]$lat <- SMN_stations %>% filter(nat_abbr == stat_ids[i]) %>% select(latitude)
        w_list[[lt]][[i]]$hei <- SMN_stations %>% filter(nat_abbr == stat_ids[i]) %>% select(station_height)
      } else {
        w_list[[lt]][[i]]$lon <- NA
        w_list[[lt]][[i]]$lat <- NA
        w_list[[lt]][[i]]$hei <- NA
      }
    
      for (j in 1:4) {
        print(c(lt, i, j))
        
        tr_dat <- dat %>% filter(nat_abbr == stat_ids[i], reftime < "2022-06-01", lubridate::month(reftime) %in% month_ind[[j]])
        ts_dat <- dat %>% filter(nat_abbr == stat_ids[i], reftime >= "2022-06-01", lubridate::month(reftime) %in% month_ind[[j]])
        y_tr <- tr_dat$obs
        y_ts <- ts_dat$obs
        
        x <- ts_dat %>% select(`COSMO-1E`, `COSMO-2E`, ECMWF_IFS) %>% as.matrix()
        x_tr <- tr_dat %>% select(`COSMO-1E`, `COSMO-2E`, ECMWF_IFS) %>% as.matrix()
        w <- get_weights(x_tr, y_tr, kernel, ind = list(1:11, 12:32, 33:83))
        
        w_df[j, ] <- w
        
      }
    
      w_list[[lt]][[i]]$w <- w_df
    }
  
    if (lt == 1) print(Sys.time() - start)
  
  }
  
  return(w_list)
}

shared_w_seas <- experiment6_seas(kernel = "Energy")

average_w_roll <- lapply(shared_w_roll, function(x) x$w) %>% Reduce("+", .)/length(shared_w_roll)
average_w_roll[rowSums(average_w_roll) == 0, ] <- NA
average_w_roll$date <- date_list
df_roll <- pivot_longer(average_w_roll, cols = c("COSMO.1E", "COSMO.2E", "ECMWF_IFS"))
df_roll <- df_roll %>% subset(date >= "2022-06-01")

lt <- 18
average_w_seas <- lapply(shared_w_seas[[lt]], function(x) x$w) %>% Reduce("+", .)/length(shared_w_seas[[lt]])
df_seas <- data.frame(date = date_list[date_list >= "2022-06-01"])
df_seas$seas <- 1
df_seas$seas[lubridate::month(df_seas$date) %in% c(3, 4, 5)] <- 2
df_seas$seas[lubridate::month(df_seas$date) %in% c(6, 7, 8)] <- 3
df_seas$seas[lubridate::month(df_seas$date) %in% c(9, 10, 11)] <- 4
df_seas$COSMO.1E <- df_seas$COSMO.2E <- df_seas$ECMWF_IFS <- NA
for (i in 1:4) {
  df_seas$COSMO.1E[df_seas$seas == i] <- average_w_seas$COSMO.1E[i]
  df_seas$COSMO.2E[df_seas$seas == i] <- average_w_seas$COSMO.2E[i]
  df_seas$ECMWF_IFS[df_seas$seas == i] <- average_w_seas$ECMWF_IFS[i]
}
df_seas <- pivot_longer(df_seas, cols = c("COSMO.1E", "COSMO.2E", "ECMWF_IFS"))

ggplot() + 
  geom_line(data = df_roll, aes(x = date, y = value, col = name)) + 
  geom_line(data = df_seas, aes(x = date, y = value, col = name)) +
  scale_x_datetime(name = "Date") +
  scale_y_continuous(name = "Weight", limits = c(0, 0.7), expand = c(0, 0)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom")
ggsave("Figures/exp6_weight_vs_date.png", width = 6, height = 4)


################################################################################
# Experiment 7: shared weights with monthly window

experiment7_mon <- function(kernel, lt_vec = 1:33) {
  
  SMN_stations <- readRDS("C:/Users/sa20i493/Documents/R/WeightedLoss/WindGusts/SMN_stations.rds")
  
  w_list <- vector("list", length(lt_vec))
  
  for (lt in lt_vec) {
    
    if (lt == 1) start <- Sys.time()
    
    dat <- load_data(lt, stat_id = NULL, path = "C:/Users/sa20i493/Documents/Data/MeteoSwiss/") # load dat
    dat <- dat %>% # restrict attention to stations with complete set of observations
      group_by(nat_abbr) %>%
      filter(n() == 1030) %>% 
      ungroup()
    
    stat_ids <- unique(dat$nat_abbr)
    
    w_list[[lt]] <- vector("list", length(stat_ids))
    
    for (i in seq_along(stat_ids)) {
      
      w_df <- data.frame(COSMO.1E = numeric(12),
                         COSMO.2E = numeric(12),
                         ECMWF_IFS = numeric(12))
      
      if (stat_ids[i] %in% SMN_stations$nat_abbr) {
        w_list[[lt]][[i]]$lon <- SMN_stations %>% filter(nat_abbr == stat_ids[i]) %>% select(longitude)
        w_list[[lt]][[i]]$lat <- SMN_stations %>% filter(nat_abbr == stat_ids[i]) %>% select(latitude)
        w_list[[lt]][[i]]$hei <- SMN_stations %>% filter(nat_abbr == stat_ids[i]) %>% select(station_height)
      } else {
        w_list[[lt]][[i]]$lon <- NA
        w_list[[lt]][[i]]$lat <- NA
        w_list[[lt]][[i]]$hei <- NA
      }
      
      for (j in 1:12) {
        print(c(lt, i, j))
        
        tr_dat <- dat %>% filter(nat_abbr == stat_ids[i], reftime < "2022-06-01", lubridate::month(reftime) == j)
        ts_dat <- dat %>% filter(nat_abbr == stat_ids[i], reftime >= "2022-06-01", lubridate::month(reftime) == j)
        y_tr <- tr_dat$obs
        y_ts <- ts_dat$obs
        
        x <- ts_dat %>% select(`COSMO-1E`, `COSMO-2E`, ECMWF_IFS) %>% as.matrix()
        x_tr <- tr_dat %>% select(`COSMO-1E`, `COSMO-2E`, ECMWF_IFS) %>% as.matrix()
        w <- get_weights(x_tr, y_tr, kernel, ind = list(1:11, 12:32, 33:83))
        
        w_df[j, ] <- w
        
      }
      
      w_list[[lt]][[i]]$w <- w_df
    }
    
    if (lt == 1) print(Sys.time() - start)
    
  }
  
  return(w_list)
}

shared_w_mon <- experiment7_mon(kernel = "Energy")


average_w_roll <- lapply(shared_w_roll, function(x) x$w) %>% Reduce("+", .)/length(shared_w_roll)
average_w_roll[rowSums(average_w_roll) == 0, ] <- NA
average_w_roll$date <- date_list
df_roll <- pivot_longer(average_w_roll, cols = c("COSMO.1E", "COSMO.2E", "ECMWF_IFS"))
df_roll <- df_roll %>% subset(date >= "2022-06-01")

lt <- 18
average_w_mon <- lapply(shared_w_mon[[lt]], function(x) x$w) %>% Reduce("+", .)/length(shared_w_mon[[lt]])
df_mon <- data.frame(date = date_list[date_list >= "2022-06-01"])
df_mon$mon <- lubridate::month(df_mon$date)
df_mon$COSMO.1E <- df_mon$COSMO.2E <- df_mon$ECMWF_IFS <- NA
for (i in 1:12) {
  df_mon$COSMO.1E[df_mon$mon == i] <- average_w_mon$COSMO.1E[i]
  df_mon$COSMO.2E[df_mon$mon == i] <- average_w_mon$COSMO.2E[i]
  df_mon$ECMWF_IFS[df_mon$mon == i] <- average_w_mon$ECMWF_IFS[i]
}
df_mon <- pivot_longer(df_mon, cols = c("COSMO.1E", "COSMO.2E", "ECMWF_IFS"))

ggplot() + 
  geom_line(data = df_roll, aes(x = date, y = value, col = name)) + 
  geom_line(data = df_mon, aes(x = date, y = value, col = name)) +
  scale_x_datetime(name = "Date") +
  scale_y_continuous(name = "Weight", limits = c(0, 0.7), expand = c(0, 0)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom")
ggsave("Figures/exp7_weight_vs_date.png", width = 6, height = 4)


################################################################################
# plot results

load("all_experiments.RData")


##### Experiment 1

s_mat <- crps_list[[1]]
w_mat <- w_list[[1]]
m_mat <- mse_list[[1]]

# plot crps vs lead time
df <- data.frame(s = as.vector(s_mat), 
                 lt = lt_vec, 
                 mth = rep(colnames(s_mat), each = length(lt_vec)))
ggplot(df) + geom_line(aes(x = lt, y = s, col = mth)) +
  scale_x_continuous(name = "Lead time (hours)", limits = c(0, 35), expand = c(0, 0)) +
  scale_y_continuous(name = "CRPS", limits = c(0.5, 1.5), expand = c(0, 0)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom")
ggsave("Figures/exp1_crps_vs_lead.png", width = 6, height = 4)

# plot w vs lead time
df <- data.frame(s = as.vector(w_mat),
                 lt = lt_vec, 
                 mth = rep(c(rep("COSMO-1E", 11), rep("COSMO-2E", 21), rep("ECMWF_IFS", 51)), each = length(lt_vec)),
                 mem = rep(1:83, each = length(lt_vec)), 
                 control = rep(c(1, rep(0, 10), 1, rep(0, 20), 1, rep(0, 50)), each = length(lt_vec)))

ggplot(df) + geom_line(aes(x = lt, y = s, group = interaction(as.factor(mem), as.factor(control)), lty = as.factor(control), col = mth)) +
  scale_x_continuous(name = "Lead time (hours)", limits = c(0, 35), expand = c(0, 0)) +
  scale_y_continuous(name = "Weight", limits = c(0, 0.25), expand = c(0, 0)) +
  theme_bw() +
  theme(legend.position = "none")


# plot aggregated w vs lead time
df <- data.frame(s = c(rowSums(w_mat[, 1:11]), rowSums(w_mat[, 12:32]), rowSums(w_mat[, 33:83])),
                 lt = lt_vec, 
                 mth = rep(c("COSMO-1E", "COSMO-2E", "ECMWF_IFS"), each = length(lt_vec)))

ggplot(df) + geom_line(aes(x = lt, y = s, col = mth)) +
  scale_x_continuous(name = "Lead time (hours)", limits = c(0, 35), expand = c(0, 0)) +
  scale_y_continuous(name = "Weight", limits = c(0, 1), expand = c(0, 0)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom")
ggsave("Figures/exp1_weight_vs_lead.png", width = 6, height = 4)

# plot mse vs weight
lt <- 18
df <- data.frame(s = w_mat[lt, ], mse = m_mat[lt, ], 
                 mth = c(rep("COSMO-1E", 11), rep("COSMO-2E", 21), rep("ECMWF_IFS", 51)),
                 control = c(1, rep(0, 10), 1, rep(0, 20), 1, rep(0, 50)))
ggplot(df) + geom_point(aes(x = mse, y = s, col = mth, shape = as.factor(control))) +
  scale_x_continuous(name = "MSE") +
  scale_y_continuous(name = "Weight") +
  scale_shape_manual(values = c(19, 4)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom") +
  guides(shape = "none")
ggsave("Figures/exp1_weight_vs_mse_18.png", width = 6, height = 4)



##### Experiment 2

s_mat <- crps_list[[2]]
w_mat <- w_list[[2]]$w
w_mat_cl <- w_list[[2]]$w_cl

# plot crps vs lead time
df <- data.frame(s = as.vector(s_mat[, -9]), 
                 lt = lt_vec, 
                 mth = rep(c(colnames(s_mat[, -c(8,9)]), "EMOS"), each = length(lt_vec)))
ggplot(df) + geom_line(aes(x = lt, y = s, col = mth)) +
  scale_x_continuous(name = "Lead time (hours)", limits = c(0, 35), expand = c(0, 0)) +
  scale_y_continuous(name = "CRPS", limits = c(0.4, 1.5), expand = c(0, 0)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom")
ggsave("Figures/exp2_crps_vs_lead.png", width = 6, height = 4)

# plot w vs lead time
df <- data.frame(s = as.vector(w_mat),
                 lt = lt_vec, 
                 mth = rep(c(rep("COSMO-1E", 11), rep("COSMO-2E", 21), rep("ECMWF_IFS", 51)), each = length(lt_vec)),
                 mem = rep(1:83, each = length(lt_vec)), 
                 control = rep(c(1, rep(0, 10), 1, rep(0, 20), 1, rep(0, 50)), each = length(lt_vec)))

ggplot(df) + geom_line(aes(x = lt, y = s, group = interaction(as.factor(mem), as.factor(control)), lty = as.factor(control), col = mth)) +
  scale_x_continuous(name = "Lead time (hours)", limits = c(0, 35), expand = c(0, 0)) +
  scale_y_continuous(name = "Weight", limits = c(0, 0.25), expand = c(0, 0)) +
  theme_bw() +
  theme(legend.position = "none")

# plot aggregated w vs lead time
df <- data.frame(s = c(rowSums(w_mat[, 1:11]), rowSums(w_mat[, 12:32]), rowSums(w_mat[, 33:83])),
                 lt = lt_vec, 
                 mth = rep(c("COSMO-1E", "COSMO-2E", "ECMWF_IFS"), each = length(lt_vec)))

ggplot(df) + geom_line(aes(x = lt, y = s, col = mth)) +
  scale_x_continuous(name = "Lead time (hours)", limits = c(0, 35), expand = c(0, 0)) +
  scale_y_continuous(name = "Weight", limits = c(0, 1), expand = c(0, 0)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom")
ggsave("Figures/exp2_weight_vs_lead.png", width = 6, height = 4)

# plot aggregated w_cl vs lead time
df <- data.frame(s = c(rowSums(w_mat_cl[, 1:11]), rowSums(w_mat_cl[, 12:32]), rowSums(w_mat_cl[, 33:83]), rowSums(w_mat_cl[, 84:ncol(w_mat_cl)])),
                 lt = lt_vec, 
                 mth = rep(c("COSMO-1E", "COSMO-2E", "ECMWF_IFS", "Past obs"), each = length(lt_vec)))

ggplot(df) + geom_line(aes(x = lt, y = s, col = mth)) +
  scale_x_continuous(name = "Lead time (hours)", limits = c(0, 35), expand = c(0, 0)) +
  scale_y_continuous(name = "Weight", limits = c(0, 1), expand = c(0, 0)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom")
ggsave("Figures/exp2_weight_cl_vs_lead.png", width = 6, height = 4)


# pit histograms

#df <- pit_df %>% pivot_longer(-obs)
pit_hist(pit_df$C1, bins = 11, ranks = F, ymax = 0.4, title = "COSMO-1E")
ggsave("Figures/exp2_pit_C1_PAY_18.png", width = 3, height = 3)

pit_hist(pit_df$C2, bins = 11, ranks = F, ymax = 0.4, title = "COSMO-2E")
ggsave("Figures/exp2_pit_C2_PAY_18.png", width = 3, height = 3)

pit_hist(pit_df$IFS, bins = 11, ranks = F, ymax = 0.4, title = "IFS")
ggsave("Figures/exp2_pit_IFS_PAY_18.png", width = 3, height = 3)

pit_hist(pit_df$Clim, bins = 11, ranks = F, ymax = 0.4, title = "Clim.")
ggsave("Figures/exp2_pit_cl_PAY_18.png", width = 3, height = 3)

pit_hist(pit_df$Multi, bins = 11, ranks = F, ymax = 0.4, title = "Multimodel")
ggsave("Figures/exp2_pit_mm_PAY_18.png", width = 3, height = 3)

pit_hist(pit_df$Weighted, bins = 11, ranks = F, ymax = 0.4, title = "Weighted")
ggsave("Figures/exp2_pit_w_PAY_18.png", width = 3, height = 3)

pit_hist(pit_df$WeightedClim, bins = 11, ranks = F, ymax = 0.4, title = "Weighted + Clim.")
ggsave("Figures/exp2_pit_wcl_PAY_18.png", width = 3, height = 3)


##### Experiment 3

s_mat <- es_list[[3]]
w_mat <- w_list[[3]]$w
w_mat_cl <- w_list[[3]]$w_cl
m_mat <- mse_list[[3]]

# plot es vs lead time
df <- data.frame(s = as.vector(s_mat), 
                 lt = lt_vec, 
                 mth = rep(colnames(s_mat), each = length(lt_vec)))
ggplot(df) + geom_line(aes(x = lt, y = s, col = mth)) +
  scale_x_continuous(name = "Lead time (hours)", limits = c(0, 35), expand = c(0, 0)) +
  scale_y_continuous(name = "ES", limits = c(6, 22), expand = c(0, 0)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom")
ggsave("Figures/exp3_es_vs_lead.png", width = 6, height = 4)

# plot w vs lead time
df <- data.frame(s = as.vector(w_mat),
                 lt = lt_vec, 
                 mth = rep(c(rep("COSMO-1E", 11), rep("COSMO-2E", 21), rep("ECMWF_IFS", 51)), each = length(lt_vec)),
                 mem = rep(1:83, each = length(lt_vec)), 
                 control = rep(c(1, rep(0, 10), 1, rep(0, 20), 1, rep(0, 50)), each = length(lt_vec)))

ggplot(df) + geom_line(aes(x = lt, y = s, group = interaction(as.factor(mem), as.factor(control)), lty = as.factor(control), col = mth)) +
  scale_x_continuous(name = "Lead time (hours)", limits = c(0, 35), expand = c(0, 0)) +
  scale_y_continuous(name = "Weight", limits = c(0, 0.25), expand = c(0, 0)) +
  theme_bw() +
  theme(legend.position = "none")

# plot aggregated w vs lead time
df <- data.frame(s = c(rowSums(w_mat[, 1:11]), rowSums(w_mat[, 12:32]), rowSums(w_mat[, 33:83])),
                 lt = lt_vec, 
                 mth = rep(c("COSMO-1E", "COSMO-2E", "ECMWF_IFS"), each = length(lt_vec)))

ggplot(df) + geom_line(aes(x = lt, y = s, col = mth)) +
  scale_x_continuous(name = "Lead time (hours)", limits = c(0, 35), expand = c(0, 0)) +
  scale_y_continuous(name = "Weight", limits = c(0, 1), expand = c(0, 0)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom")
ggsave("Figures/exp3_weight_vs_lead.png", width = 6, height = 4)

# plot aggregated w_cl vs lead time
df <- data.frame(s = c(rowSums(w_mat_cl[, 1:11]), rowSums(w_mat_cl[, 12:32]), rowSums(w_mat_cl[, 33:83]), rowSums(w_mat_cl[, 84:ncol(w_mat_cl)])),
                 lt = lt_vec, 
                 mth = rep(c("COSMO-1E", "COSMO-2E", "ECMWF_IFS", "Past obs"), each = length(lt_vec)))

ggplot(df) + geom_line(aes(x = lt, y = s, col = mth)) +
  scale_x_continuous(name = "Lead time (hours)", limits = c(0, 35), expand = c(0, 0)) +
  scale_y_continuous(name = "Weight", limits = c(0, 1), expand = c(0, 0)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom")
ggsave("Figures/exp3_weight_cl_vs_lead.png", width = 6, height = 4)

# plot mse vs weight
lt <- 18
df <- data.frame(s = w_mat[lt, ], mse = m_mat[lt, ], 
                 mth = c(rep("COSMO-1E", 11), rep("COSMO-2E", 21), rep("ECMWF_IFS", 51)),
                 control = c(1, rep(0, 10), 1, rep(0, 20), 1, rep(0, 50)))
ggplot(df) + geom_point(aes(x = mse, y = s, col = mth, shape = as.factor(control))) +
  scale_x_continuous(name = "MSE") +
  scale_y_continuous(name = "Weight") +
  scale_shape_manual(values = c(19, 4)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom") +
  guides(shape = "none")
ggsave("Figures/exp3_weight_vs_mse_18.png", width = 6, height = 4)








################################################################################
### analog-based post-processing (local)

library(provenance)

experimenta0 <- function(kernel = "energy", lt = 18, M = 21, K = 10) {
  
  dat <- load_data(lt, stat_id = NULL, path = "C:/Users/sa20i493/Documents/Data/MeteoSwiss/") # load dat
  dat <- dat %>% # restrict attention to stations with complete set of observations
    group_by(nat_abbr) %>%
    filter(n() == 1030) %>% 
    ungroup()
  
  stat_ids <- unique(dat$nat_abbr)

  # restrict to cosmo-2e
  dat <- dat %>% select(-c("COSMO-1E", "ECMWF_IFS"))
  
  crps_mat <- matrix(NA, nrow = length(stat_ids), ncol = 5)
  colnames(crps_mat) <- c("COSMO-2E", "Analog", "Weight", "Weight_tap", "RKHS")
  
  for (i in seq_along(stat_ids)) {
    print(c(lt, i))
    
    tr_dat <- dat %>% filter(nat_abbr == stat_ids[i], reftime < "2022-06-01")
    n_tr <- tr_dat %>% nrow()
    ts_dat <- dat %>% filter(nat_abbr == stat_ids[i], reftime >= "2022-06-01")
    n_ts <- ts_dat %>% nrow()
    y_tr <- tr_dat$obs
    y_ts <- ts_dat$obs
    
    ksmat <- sapply(1:n_ts, function(j) {
      print(j)
      sapply(1:n_tr, function(k) {
        x_tr <- tr_dat$`COSMO-2E`[k, ] %>% unlist() %>% unname()
        x_ts <- ts_dat$`COSMO-2E`[j, ] %>% unlist() %>% unname()
        ks.test(x_tr, x_ts)$statistic %>% unname()
    })})
    
    ksmat_tr <- sapply(1:n_tr, function(j) { # can be simplified by only calculating half since it's diagonal
      print(j)
      sapply(1:n_tr, function(k) {
        x_1 <- tr_dat$`COSMO-2E`[k, ] %>% unlist() %>% unname()
        x_2 <- tr_dat$`COSMO-2E`[j, ] %>% unlist() %>% unname()
        ks.test(x_1, x_2)$statistic %>% unname()
      })})
    
    ## cosmo-1e
    x <- ts_dat %>% select(`COSMO-2E`) %>% as.matrix()
    crps_mat[i, 1] <- crps_sample(y_ts, x) %>% mean()
    
    ## analogs
    x <- sapply(1:n_ts, function(j) {
      ind <- which(rank(ksmat[ ,j], ties.method = "random") <= M)
      tr_dat$obs[ind]
    }) %>% t()
    crps_mat[i, 2] <- crps_sample(y_ts, x) %>% mean()
    
    ## inverse weighted analog
    w <- t(1/ksmat)
    w <- w/rowSums(w)
    x <- replicate(n_ts, tr_dat$obs) %>% t()
    crps_mat[i, 3] <- crps_sample(y_ts, x, w = w) %>% mean()
    
    ## tapered weighted analog
    for (j in 1:n_ts) {
      ind <- which(rank(-w[j, ], ties.method = "random") > M)
      w[j, ind] <- 0
    }
    w <- w/rowSums(w)
    crps_mat[i, 4] <- crps_sample(y_ts, x, w = w) %>% mean()
    
    ## rkhs weights
    x_tr <- sapply(1:n_tr, function(j) x[1, order(ksmat_tr[, j])]) %>% t()
    ind <- split(1:n_tr, cut(1:n_tr, breaks = K, labels = F))
    w <- get_weights(x_tr, y_tr, kernel, ind = ind)
    x_ts <- sapply(1:n_ts, function(j) x[j, order(ksmat[, j])]) %>% t()
    wmat <- replicate(n_ts, rep(w, sapply(ind, length))) %>% t()
    crps_mat[i, 5] <- crps_sample(y_ts, x_ts, w = wmat) %>% mean()
  }
  
  return(crps_mat)
}

result <- experimenta0()




