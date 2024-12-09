################################################################################
##### set up

library(kernlab)
library(scoringRules)
library(Rcpp)
library(ggplot2)
source("utility_funcs.R")
sourceCpp("kernels.cpp")



# functions
get_stations <- function(lt_vec = 1:33, new = FALSE) {
  
  if (new) {
    stat_ids <- vector("list", length(lt_vec))
    
    for (lt in lt_vec) {
      print(lt)
      dat <- load_data(lt, stat_id = NULL, path = "C:/Users/sa20i493/Documents/Data/MeteoSwiss/") # load dat
      dat <- dat %>% # restrict attention to stations with complete set of observations
        group_by(nat_abbr) %>%
        filter(n() == 1030) %>% 
        ungroup()
      
      stat_ids[[lt]] <- unique(dat$nat_abbr)
    }
    
    stat_list <- Reduce(intersect, stat_ids)
    save(stat_list, file = "Data/stat_list.RData")
  } else {
    load("Data/stat_list.RData")
    return(stat_list)
  }
  
}

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

get_mv_weights <- function(x, y, kernel = "Gaussian", ind = NULL) {
  
  H <- get_mv_H(x, kernel)
  c <- get_mv_c(x, y, kernel)
  
  if (!is.null(ind)) {
    M <- length(ind)
    H <- sapply(1:M, function(i) sapply(1:M, function(j) mean(H[ind[[i]], ind[[j]]])))
    c <- sapply(1:M, function(i) mean(c[ind[[i]]]))
  } else {
    M <- dim(x)[3]
  }
  
  r <- 0
  l <- rep(0, M)
  u <- rep(1, M)
  b <- 1
  A <- rep(1, M)
  
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

get_scores <- function(tr_dat, ts_dat, kernel = "Gaussian", clim = FALSE) {
  
  if (clim) {
    crps_vec <- rep(NA, 11)
    names(crps_vec) <- c("Clim.", "C1", "C2", "IFS", "MM", "LP", "Wtd", "Wtd-Ord", "LP+Clim", "Wtd+Clim", "Wtd-Ord+Clim")
  } else {
    crps_vec <- rep(NA, 8)
    names(crps_vec) <- c("Clim.", "C1", "C2", "IFS", "MM", "LP", "Wtd", "Wtd-Ord")
  }
  
  y_tr <- tr_dat$obs
  y_ts <- ts_dat$obs
  
  
  # Climatology
  x <- t(replicate(length(y_ts), y_tr))
  crps_vec[1] <- mean(crps_sample(y_ts, x))
  
  
  # COSMO-1E
  x <- ts_dat %>% dplyr::select(`COSMO-1E`) %>% as.matrix()
  crps_vec[2] <- mean(crps_sample(y_ts, x))
  
  
  # COSMO-2E
  x <- ts_dat %>% dplyr::select(`COSMO-2E`) %>% as.matrix()
  crps_vec[3] <- mean(crps_sample(y_ts, x))
  
  
  # IFS
  x <- ts_dat %>% dplyr::select(ECMWF_IFS) %>% as.matrix()
  crps_vec[4] <- mean(crps_sample(y_ts, x))
  
  
  # Multi-model
  x <- ts_dat %>% dplyr::select(`COSMO-1E`, `COSMO-2E`, ECMWF_IFS) %>% as.matrix()
  crps_vec[5] <- mean(crps_sample(y_ts, x))
  
  
  # Linear pool
  x_tr <- tr_dat %>% dplyr::select(`COSMO-1E`, `COSMO-2E`, ECMWF_IFS) %>% as.matrix()
  w_lp <- get_weights(x_tr, y_tr, kernel, ind = list(1:11, 12:32, 33:83))
  w <- t(replicate(length(y_ts), rep(w_lp, c(11, 21, 51))))
  
  crps_vec[6] <- mean(crps_sample(y_ts, x, w = w))
  
  
  # Weighted
  w_wtd <- get_weights(x_tr, y_tr, kernel)
  w <- t(replicate(length(y_ts), w_wtd))
  
  crps_vec[7] <- mean(crps_sample(y_ts, x, w = w))
  
  
  # Weighted order statistics
  x_c1 <- tr_dat %>% dplyr::select(`COSMO-1E`) %>% as.matrix() %>% apply(., 1, sort) %>% t()
  x_c2 <- tr_dat %>% dplyr::select(`COSMO-2E`) %>% as.matrix() %>% apply(., 1, sort) %>% t()
  x_ifs <- tr_dat %>% dplyr::select(ECMWF_IFS) %>% as.matrix() %>% apply(., 1, sort) %>% t()
  x_ord_tr <- cbind(x_c1, x_c2, x_ifs)
  w_ord <- get_weights(x_ord_tr, y_tr, kernel)
  w <- t(replicate(length(y_ts), w_ord))
  
  x_c1 <- ts_dat %>% dplyr::select(`COSMO-1E`) %>% as.matrix() %>% apply(., 1, sort) %>% t()
  x_c2 <- ts_dat %>% dplyr::select(`COSMO-2E`) %>% as.matrix() %>% apply(., 1, sort) %>% t()
  x_ifs <- ts_dat %>% dplyr::select(ECMWF_IFS) %>% as.matrix() %>% apply(., 1, sort) %>% t()
  x_ord <- cbind(x_c1, x_c2, x_ifs)
  crps_vec[8] <- mean(crps_sample(y_ts, x_ord, w = w))
  
  
  
  if (clim) {
    # Linear pool + Clim
    x_tr <- tr_dat %>% dplyr::select(`COSMO-1E`, `COSMO-2E`, ECMWF_IFS) %>% as.matrix()
    x_tr <- cbind(x_tr, t(replicate(length(y_tr), y_tr)))
    w_lp_cl <- get_weights(x_tr, y_tr, kernel, ind = list(1:11, 12:32, 33:83, 84:ncol(x_tr)))
    w <- t(replicate(length(y_ts), rep(w_lp_cl, c(11, 21, 51, (ncol(x_tr) - 83)))))
    
    x <- ts_dat %>% dplyr::select(`COSMO-1E`, `COSMO-2E`, ECMWF_IFS) %>% as.matrix()
    x <- cbind(x, t(replicate(length(y_ts), y_tr)))
    crps_vec[9] <- mean(crps_sample(y_ts, x, w = w))
    
    
    # Weighted + Clim
    w_wtd_cl <- get_weights(x_tr, y_tr, kernel)
    w <- t(replicate(length(y_ts), w_wtd_cl))
    
    x <- cbind(x, t(replicate(length(y_ts), y_tr)))
    crps_vec[10] <- mean(crps_sample(y_ts, x, w = w))
    
    
    # Weighted order + Clim
    x_tr <- cbind(x_ord_tr, t(replicate(length(y_tr), y_tr)))
    w_ord_cl <- get_weights(x_tr, y_tr, kernel)
    w <- t(replicate(length(y_ts), w_ord_cl))
    
    x <- cbind(x_ord, t(replicate(length(y_ts), y_tr)))
    crps_vec[11] <- mean(crps_sample(y_ts, x, w = w))
    
    
    w <- list(lp = w_lp, lp_cl = w_lp_cl, wtd = w_wtd, 
              wtd_cl = w_wtd_cl, ord = w_ord, ord_cl = w_ord_cl)
  } else {
    w <- list(lp = w_lp, wtd = w_wtd, ord = w_ord)
  }
  
  return(list(crps = crps_vec, w = w))
}

get_mbm_scores <- function(tr_dat, ts_dat, kernel = "Gaussian", clim = FALSE) {
  
  if (clim) {
    crps_vec <- rep(NA, 11)
    names(crps_vec) <- c("Clim.", "C1", "C2", "IFS", "MM", "LP", "Wtd", "Wtd-Ord", "LP+Clim", "Wtd+Clim", "Wtd-Ord+Clim")
  } else {
    crps_vec <- rep(NA, 8)
    names(crps_vec) <- c("Clim.", "C1", "C2", "IFS", "MM", "LP", "Wtd", "Wtd-Ord")
  }
  
  y_tr <- tr_dat$obs
  y_ts <- ts_dat$obs
  
  
  # Climatology
  x <- t(replicate(length(y_ts), y_tr))
  crps_vec[1] <- mean(crps_sample(y_ts, x))
  
  
  # COSMO-1E
  x_tr_c1 <- tr_dat %>% dplyr::select(`COSMO-1E`) %>% as.matrix()
  x_c1 <- ts_dat %>% dplyr::select(`COSMO-1E`) %>% as.matrix()
  x_c1 <- mbm_mom_est(y_tr, x_tr_c1, x_c1)
  crps_vec[2] <- mean(crps_sample(y_ts, x_c1))
  
  
  # COSMO-2E
  x_tr_c2 <- tr_dat %>% dplyr::select(`COSMO-2E`) %>% as.matrix()
  x_c2 <- ts_dat %>% dplyr::select(`COSMO-2E`) %>% as.matrix()
  x_c2 <- mbm_mom_est(y_tr, x_tr_c2, x_c2)
  crps_vec[3] <- mean(crps_sample(y_ts, x_c2))
  
  
  # IFS
  x_tr_ifs <- tr_dat %>% dplyr::select(ECMWF_IFS) %>% as.matrix()
  x_ifs <- ts_dat %>% dplyr::select(ECMWF_IFS) %>% as.matrix()
  x_ifs <- mbm_mom_est(y_tr, x_tr_ifs, x_ifs)
  crps_vec[4] <- mean(crps_sample(y_ts, x_ifs))
  
  
  # Multi-model
  x <- cbind(x_c1, x_c2, x_ifs)
  crps_vec[5] <- mean(crps_sample(y_ts, x))
  
  
  # Linear pool
  x_tr <- cbind(x_tr_c1, x_tr_c2, x_tr_ifs)
  x_tr <- mbm_mom_est(y_tr, x_tr, x_tr)
  w_lp <- get_weights(x_tr, y_tr, kernel, ind = list(1:11, 12:32, 33:83))
  w <- t(replicate(length(y_ts), rep(w_lp, c(11, 21, 51))))
  
  crps_vec[6] <- mean(crps_sample(y_ts, x, w = w))
  
  
  # Weighted
  w_wtd <- get_weights(x_tr, y_tr, kernel)
  w <- t(replicate(length(y_ts), w_wtd))
  
  crps_vec[7] <- mean(crps_sample(y_ts, x, w = w))
  
  
  # Weighted order statistics
  x_c1_tr <- tr_dat %>% dplyr::select(`COSMO-1E`) %>% as.matrix() 
  x_c1 <- x_c1_tr %>% mbm_mom_est(y_tr, ., .) %>% apply(., 1, sort) %>% t()
  x_c2_tr <- tr_dat %>% dplyr::select(`COSMO-2E`) %>% as.matrix() 
  x_c2 <- x_c2_tr %>% mbm_mom_est(y_tr, ., .) %>% apply(., 1, sort) %>% t()
  x_ifs_tr <- tr_dat %>% dplyr::select(ECMWF_IFS) %>% as.matrix() 
  x_ifs <- x_ifs_tr %>% mbm_mom_est(y_tr, ., .) %>% apply(., 1, sort) %>% t()
  x_ord_tr <- cbind(x_c1, x_c2, x_ifs)
  w_ord <- get_weights(x_ord_tr, y_tr, kernel)
  w <- t(replicate(length(y_ts), w_ord))
  
  x_c1 <- ts_dat %>% dplyr::select(`COSMO-1E`) %>% as.matrix() %>% 
    mbm_mom_est(y_tr, x_c1_tr, .) %>% apply(., 1, sort) %>% t()
  x_c2 <- ts_dat %>% dplyr::select(`COSMO-2E`) %>% as.matrix() %>%
    mbm_mom_est(y_tr, x_c2_tr, .) %>% apply(., 1, sort) %>% t()
  x_ifs <- ts_dat %>% dplyr::select(ECMWF_IFS) %>% as.matrix() %>%
    mbm_mom_est(y_tr, x_ifs_tr, .) %>% apply(., 1, sort) %>% t()
  x_ord <- cbind(x_c1, x_c2, x_ifs)
  crps_vec[8] <- mean(crps_sample(y_ts, x_ord, w = w))
  
  
  if (clim) {
    # Linear pool + Clim
    x_tr <- tr_dat %>% dplyr::select(`COSMO-1E`, `COSMO-2E`, ECMWF_IFS) %>% as.matrix()
    x_tr <- mbm_mom_est(y_tr, x_tr, x_tr)
    x_tr <- cbind(x_tr, t(replicate(length(y_tr), y_tr)))
    w_lp_cl <- get_weights(x_tr, y_tr, kernel, ind = list(1:11, 12:32, 33:83, 84:ncol(x_tr)))
    w <- t(replicate(length(y_ts), rep(w_lp_cl, c(11, 21, 51, (ncol(x_tr) - 83)))))
    
    x <- ts_dat %>% dplyr::select(`COSMO-1E`, `COSMO-2E`, ECMWF_IFS) %>% as.matrix()
    x <- mbm_mom_est(y_tr, x_tr, x)
    x <- cbind(x, t(replicate(length(y_ts), y_tr)))
    crps_vec[9] <- mean(crps_sample(y_ts, x, w = w))
    
    
    # Weighted + Clim
    w_wtd_cl <- get_weights(x_tr, y_tr, kernel)
    w <- t(replicate(length(y_ts), w_wtd_cl))
    
    x <- cbind(x, t(replicate(length(y_ts), y_tr)))
    crps_vec[10] <- mean(crps_sample(y_ts, x, w = w))
    
    
    # Weighted order + Clim
    x_tr <- cbind(x_ord_tr, t(replicate(length(y_tr), y_tr)))
    w_ord_cl <- get_weights(x_tr, y_tr, kernel)
    w <- t(replicate(length(y_ts), w_ord_cl))
    
    x <- cbind(x_ord, t(replicate(length(y_ts), y_tr)))
    crps_vec[11] <- mean(crps_sample(y_ts, x, w = w))
    
    
    w <- list(lp = w_lp, lp_cl = w_lp_cl, wtd = w_wtd, 
              wtd_cl = w_wtd_cl, ord = w_ord, ord_cl = w_ord_cl)
  } else {
    w <- list(lp = w_lp, wtd = w_wtd, ord = w_ord)
  }
  
  return(list(crps = crps_vec, w = w))
}

get_mv_scores  <- function(tr_dat, ts_dat, stat_ids, kernel = "Gaussian", clim = FALSE) {
  
  if (clim) {
    crps_mat <- matrix(NA, length(stat_ids), 9)
    es_vec <- rep(NA, 9)
    colnames(crps_mat) <- names(es_vec) <- c("Clim.", "C1", "C2", "IFS", "MM", "LP", "Wtd", "LP+Clim", "Wtd+Clim")
  } else {
    crps_mat <- matrix(NA, length(stat_ids), 7)
    es_vec <- rep(NA, 7)
    colnames(crps_mat) <- names(es_vec) <- c("Clim.", "C1", "C2", "IFS", "MM", "LP", "Wtd")
  }
  
  y_tr <- tr_dat %>% dplyr::select(nat_abbr, obs) %>% 
    group_by(nat_abbr) %>%
    mutate(row = row_number()) %>%
    pivot_wider(names_from = nat_abbr, values_from = obs) %>%
    dplyr::select(-row) %>% as.matrix()
  y_ts <- ts_dat %>% dplyr::select(nat_abbr, obs) %>% 
    group_by(nat_abbr) %>%
    mutate(row = row_number()) %>%
    pivot_wider(names_from = nat_abbr, values_from = obs) %>%
    dplyr::select(-row) %>% as.matrix()
  n_tr <- nrow(y_tr)
  n_ts <- nrow(y_ts)
  d <- ncol(y_tr)
  
  # Climatology
  y_mat <- replicate(n_ts, t(y_tr)) %>% aperm(c(3, 1, 2))
  es_vec[1] <- mean(sapply(1:n_ts, function(i) es_sample(y_ts[i, ], t(y_tr))))
  crps_mat[, 1] <- sapply(1:d, function(j) mean(crps_sample(y_ts[, j], y_mat[, j, ])))
  
  
  # COSMO-1E
  x <- lapply(stat_ids, function(z) {
    ts_dat %>% filter(nat_abbr == z) %>% 
      dplyr::select(`COSMO-1E`) %>%
      cbind() %>% as.matrix()
  }) %>% simplify2array() %>% aperm(c(1, 3, 2))
  es_vec[2] <- mean(sapply(1:n_ts, function(i) es_sample(y_ts[i, ], x[i, , ])))
  crps_mat[, 2] <- sapply(1:d, function(j) mean(crps_sample(y_ts[, j], x[, j, ])))
  
  
  # COSMO-2E
  x <- lapply(stat_ids, function(z) {
    ts_dat %>% filter(nat_abbr == z) %>% 
      dplyr::select(`COSMO-2E`) %>%
      cbind() %>% as.matrix()
  }) %>% simplify2array() %>% aperm(c(1, 3, 2))
  es_vec[3] <- mean(sapply(1:n_ts, function(i) es_sample(y_ts[i, ], x[i, , ])))
  crps_mat[, 3] <- sapply(1:d, function(j) mean(crps_sample(y_ts[, j], x[, j, ])))
  
  
  # IFS
  x <- lapply(stat_ids, function(z) {
    ts_dat %>% filter(nat_abbr == z) %>% 
      dplyr::select(ECMWF_IFS) %>%
      cbind() %>% as.matrix()
  }) %>% simplify2array() %>% aperm(c(1, 3, 2))
  es_vec[4] <- mean(sapply(1:n_ts, function(i) es_sample(y_ts[i, ], x[i, , ])))
  crps_mat[, 4] <- sapply(1:d, function(j) mean(crps_sample(y_ts[, j], x[, j, ])))
  
  
  # Multi-model
  x <- lapply(stat_ids, function(z) {
    ts_dat %>% filter(nat_abbr == z) %>% 
      dplyr::select(`COSMO-1E`, `COSMO-2E`, ECMWF_IFS) %>%
      cbind() %>% as.matrix()
  }) %>% simplify2array() %>% aperm(c(1, 3, 2))
  es_vec[5] <- mean(sapply(1:n_ts, function(i) es_sample(y_ts[i, ], x[i, , ])))
  crps_mat[, 5] <- sapply(1:d, function(j) mean(crps_sample(y_ts[, j], x[, j, ])))
  
  
  # Linear pool
  x_tr <- lapply(stat_ids, function(z) {
    tr_dat %>% filter(nat_abbr == z) %>% 
      dplyr::select(`COSMO-1E`, `COSMO-2E`, ECMWF_IFS) %>%
      cbind() %>% as.matrix()
  }) %>% simplify2array() %>% aperm(c(1, 3, 2))
  w_lp <- get_mv_weights(x_tr, y_tr, kernel, ind = list(1:11, 12:32, 33:83))
  
  w <- rep(w_lp, c(11, 21, 51))
  es_vec[6] <- mean(sapply(1:n_ts, function(i) es_sample(y_ts[i, ], x[i, , ], w = w / sum(w))))
  w <- t(replicate(n_ts, w))
  if (nrow(y_ts) == 1) w <- as.vector(w)
  crps_mat[, 6] <- sapply(1:d, function(j) mean(crps_sample(y_ts[, j], x[, j, ], w = w)))
  
  
  # Weighted
  w_wtd <- get_mv_weights(x_tr, y_tr, kernel)
  w <- t(replicate(n_ts, w_wtd))
  if (nrow(y_ts) == 1) w <- as.vector(w)
  
  es_vec[7] <- mean(sapply(1:n_ts, function(i) es_sample(y_ts[i, ], x[i, , ], w = w_wtd)))
  crps_mat[, 7] <- sapply(1:d, function(j) mean(crps_sample(y_ts[, j], x[, j, ], w = w)))
  
  
  if (clim) {
    # Linear pool + Clim
    y_mat_tr <- replicate(n_tr, t(y_tr)) %>% aperm(c(3, 1, 2))
    x_tr <- abind::abind(x_tr, y_mat_tr, along = 3)
    w_lp_cl <- get_mv_weights(x_tr, y_tr, kernel, ind = list(1:11, 12:32, 33:83, 84:dim(x)[3]))
    
    x <- abind::abind(x, y_mat, along = 3)
    w <- rep(w_lp_cl, c(11, 21, 51, dim(x)[3] - 83))
    es_vec[8] <- mean(sapply(1:n_ts, function(i) es_sample(y_ts[i, ], x[i, , ], w = w)))
    w <- t(replicate(n_ts, w))
    if (nrow(y_ts) == 1) w <- as.vector(w)
    crps_mat[, 8] <- sapply(1:d, function(j) mean(crps_sample(y_ts[, j], x[, j, ], w = w)))
    
    
    # Weighted + Clim
    w_wtd_cl <- get_mv_weights(x_tr, y_tr, kernel)
    w <- t(replicate(n_ts, w_wtd_cl))
    if (nrow(y_ts) == 1) w <- as.vector(w)
    
    es_vec[9] <- mean(sapply(1:n_ts, function(i) es_sample(y_ts[i, ], x[i, , ], w = w_wtd_cl)))
    crps_mat[, 9] <- sapply(1:d, function(j) mean(crps_sample(y_ts[, j], x[, j, ], w = w)))
    
    w <- list(lp = w_lp, lp_cl = w_lp_cl, wtd = w_wtd, wtd_cl = w_wtd_cl)
    
  } else {
    w <- list(lp = w_lp, wtd = w_wtd)
  }
  
  return(list(crps = crps_mat, es = es_vec, w = w))
}

get_mv_mbm_scores  <- function(tr_dat, ts_dat, stat_ids, kernel = "Gaussian", clim = FALSE) {
  
  if (clim) {
    crps_mat <- matrix(NA, length(stat_ids), 9)
    es_vec <- rep(NA, 9)
    colnames(crps_mat) <- names(es_vec) <- c("Clim.", "C1", "C2", "IFS", "MM", "LP", "Wtd", "LP+Clim", "Wtd+Clim")
  } else {
    crps_mat <- matrix(NA, length(stat_ids), 7)
    es_vec <- rep(NA, 7)
    colnames(crps_mat) <- names(es_vec) <- c("Clim.", "C1", "C2", "IFS", "MM", "LP", "Wtd")
  }
  
  y_tr <- tr_dat %>% dplyr::select(nat_abbr, obs) %>% 
    group_by(nat_abbr) %>%
    mutate(row = row_number()) %>%
    pivot_wider(names_from = nat_abbr, values_from = obs) %>%
    dplyr::select(-row) %>% as.matrix()
  y_ts <- ts_dat %>% dplyr::select(nat_abbr, obs) %>% 
    group_by(nat_abbr) %>%
    mutate(row = row_number()) %>%
    pivot_wider(names_from = nat_abbr, values_from = obs) %>%
    dplyr::select(-row) %>% as.matrix()
  n_tr <- nrow(y_tr)
  n_ts <- nrow(y_ts)
  d <- ncol(y_tr)
  
  # Climatology
  y_mat <- replicate(n_ts, t(y_tr)) %>% aperm(c(3, 1, 2))
  es_vec[1] <- mean(sapply(1:n_ts, function(i) es_sample(y_ts[i, ], t(y_tr))))
  crps_mat[, 1] <- sapply(1:d, function(j) mean(crps_sample(y_ts[, j], y_mat[, j, ])))
  
  
  # COSMO-1E
  x_tr <- lapply(stat_ids, function(z) {
    tr_dat %>% filter(nat_abbr == z) %>% 
      dplyr::select(`COSMO-1E`) %>%
      cbind() %>% as.matrix()
  }) %>% simplify2array() %>% aperm(c(1, 3, 2))
  x <- lapply(stat_ids, function(z) {
    ts_dat %>% filter(nat_abbr == z) %>% 
      dplyr::select(`COSMO-1E`) %>%
      cbind() %>% as.matrix()
  }) %>% simplify2array() %>% aperm(c(1, 3, 2))
  x <- lapply(1:d, function(j) mbm_mom_est(y_tr[, j], x_tr[, j, ], x[, j, ]))
  x <- x %>% simplify2array() %>% aperm(c(1, 3, 2))
  es_vec[2] <- mean(sapply(1:n_ts, function(i) es_sample(y_ts[i, ], x[i, , ])))
  crps_mat[, 2] <- sapply(1:d, function(j) mean(crps_sample(y_ts[, j], x[, j, ])))
  
  
  # COSMO-2E
  x_tr <- lapply(stat_ids, function(z) {
    tr_dat %>% filter(nat_abbr == z) %>% 
      dplyr::select(`COSMO-2E`) %>%
      cbind() %>% as.matrix()
  }) %>% simplify2array() %>% aperm(c(1, 3, 2))
  x <- lapply(stat_ids, function(z) {
    ts_dat %>% filter(nat_abbr == z) %>% 
      dplyr::select(`COSMO-2E`) %>%
      cbind() %>% as.matrix()
  }) %>% simplify2array() %>% aperm(c(1, 3, 2))
  x <- lapply(1:d, function(j) mbm_mom_est(y_tr[, j], x_tr[, j, ], x[, j, ]))
  x <- x %>% simplify2array() %>% aperm(c(1, 3, 2))
  es_vec[3] <- mean(sapply(1:n_ts, function(i) es_sample(y_ts[i, ], x[i, , ])))
  crps_mat[, 3] <- sapply(1:d, function(j) mean(crps_sample(y_ts[, j], x[, j, ])))
  
  
  # IFS
  x_tr <- lapply(stat_ids, function(z) {
    tr_dat %>% filter(nat_abbr == z) %>% 
      dplyr::select(ECMWF_IFS) %>%
      cbind() %>% as.matrix()
  }) %>% simplify2array() %>% aperm(c(1, 3, 2))
  x <- lapply(stat_ids, function(z) {
    ts_dat %>% filter(nat_abbr == z) %>% 
      dplyr::select(ECMWF_IFS) %>%
      cbind() %>% as.matrix()
  }) %>% simplify2array() %>% aperm(c(1, 3, 2))
  x <- lapply(1:d, function(j) mbm_mom_est(y_tr[, j], x_tr[, j, ], x[, j, ]))
  x <- x %>% simplify2array() %>% aperm(c(1, 3, 2))
  es_vec[4] <- mean(sapply(1:n_ts, function(i) es_sample(y_ts[i, ], x[i, , ])))
  crps_mat[, 4] <- sapply(1:d, function(j) mean(crps_sample(y_ts[, j], x[, j, ])))

  
  # Multi-model
  x_tr <- lapply(stat_ids, function(z) {
    tr_dat %>% filter(nat_abbr == z) %>% 
      dplyr::select(`COSMO-1E`, `COSMO-2E`, ECMWF_IFS) %>%
      cbind() %>% as.matrix()
  }) %>% simplify2array() %>% aperm(c(1, 3, 2))
  x <- lapply(stat_ids, function(z) {
    ts_dat %>% filter(nat_abbr == z) %>% 
      dplyr::select(`COSMO-1E`, `COSMO-2E`, ECMWF_IFS) %>%
      cbind() %>% as.matrix()
  }) %>% simplify2array() %>% aperm(c(1, 3, 2))
  x <- lapply(1:d, function(j) mbm_mom_est(y_tr[, j], x_tr[, j, ], x[, j, ]))
  x <- x %>% simplify2array() %>% aperm(c(1, 3, 2))
  es_vec[5] <- mean(sapply(1:n_ts, function(i) es_sample(y_ts[i, ], x[i, , ])))
  crps_mat[, 5] <- sapply(1:d, function(j) mean(crps_sample(y_ts[, j], x[, j, ])))
  
  
  # Linear pool
  x_tr <- lapply(1:d, function(j) mbm_mom_est(y_tr[, j], x_tr[, j, ], x_tr[, j, ]))
  x_tr <- x_tr %>% simplify2array() %>% aperm(c(1, 3, 2))
  w_lp <- get_mv_weights(x_tr, y_tr, kernel, ind = list(1:11, 12:32, 33:83))
  
  w <- rep(w_lp, c(11, 21, 51))
  es_vec[6] <- mean(sapply(1:n_ts, function(i) es_sample(y_ts[i, ], x[i, , ], w = w / sum(w))))
  w <- t(replicate(n_ts, w))
  if (nrow(y_ts) == 1) w <- as.vector(w)
  crps_mat[, 6] <- sapply(1:d, function(j) mean(crps_sample(y_ts[, j], x[, j, ], w = w)))
  
  
  # Weighted
  w_wtd <- get_mv_weights(x_tr, y_tr, kernel)
  w <- t(replicate(n_ts, w_wtd))
  if (nrow(y_ts) == 1) w <- as.vector(w)
  
  es_vec[7] <- mean(sapply(1:n_ts, function(i) es_sample(y_ts[i, ], x[i, , ], w = w_wtd)))
  crps_mat[, 7] <- sapply(1:d, function(j) mean(crps_sample(y_ts[, j], x[, j, ], w = w)))
  
  
  if (clim) {
    # Linear pool + Clim
    y_mat_tr <- replicate(n_tr, t(y_tr)) %>% aperm(c(3, 1, 2))
    x_tr <- abind::abind(x_tr, y_mat_tr, along = 3)
    w_lp_cl <- get_mv_weights(x_tr, y_tr, kernel, ind = list(1:11, 12:32, 33:83, 84:dim(x)[3]))
    
    x <- abind::abind(x, y_mat, along = 3)
    w <- rep(w_lp_cl, c(11, 21, 51, dim(x)[3] - 83))
    es_vec[8] <- mean(sapply(1:n_ts, function(i) es_sample(y_ts[i, ], x[i, , ], w = w)))
    w <- t(replicate(n_ts, w))
    if (nrow(y_ts) == 1) w <- as.vector(w)
    crps_mat[, 8] <- sapply(1:d, function(j) mean(crps_sample(y_ts[, j], x[, j, ], w = w)))
    
    
    # Weighted + Clim
    w_wtd_cl <- get_mv_weights(x_tr, y_tr, kernel)
    w <- t(replicate(n_ts, w_wtd_cl))
    if (nrow(y_ts) == 1) w <- as.vector(w)
    
    es_vec[9] <- mean(sapply(1:n_ts, function(i) es_sample(y_ts[i, ], x[i, , ], w = w_wtd_cl)))
    crps_mat[, 9] <- sapply(1:d, function(j) mean(crps_sample(y_ts[, j], x[, j, ], w = w)))
  
    w <- list(lp = w_lp, lp_cl = w_lp_cl, wtd = w_wtd, wtd_cl = w_wtd_cl)
    
  } else {
    w <- list(lp = w_lp, wtd = w_wtd)
  }

  return(list(crps = crps_mat, es = es_vec, w = w))
}


# member-by-member post-processing

mbm_mom_est <- function(y, dat_tr, dat_ts, sc = "sqrt") {
  
  if (sc == "log") {
    dat_tr <- log(dat_tr)
    dat_ts <- log(dat_ts)
    y <- log(y)
  } else if (sc == "sqrt") {
    dat_tr <- sqrt(dat_tr)
    dat_ts <- sqrt(dat_ts)
    y <- sqrt(y)
  }
  
  xbar <- rowMeans(dat_tr)
  ybar <- mean(y)
  s2 <- apply(dat_tr, 1, var)
  
  b <- cov(xbar, y) / var(xbar)
  a <- ybar - b*mean(xbar)
  c <- (cov(s2, y^2) - 2*a*b*cov(s2, xbar) - (b^2)*cov(s2, xbar^2)) / var(s2)
  d <- var(y) - c*mean(s2) - (b^2)*var(xbar)
  
  if (d < 0) {
    d <- 0
    c <- (var(y) - (b^2)*var(xbar)) / mean(s2)
  }
  c <- max(c, 0)
  
  xbar <- rowMeans(dat_ts)
  s2 <- apply(dat_ts, 1, var)
  gamma <- sqrt(c + d/s2)
  newdat <- (a + b*xbar) + gamma*(dat_ts - xbar)
  
  if (sc == "log") {
    newdat <- exp(newdat)
  } else if (sc == "sqrt") {
    newdat <- newdat^2
  }
  
  return(newdat)
}


# lead times
lt_vec <- 1:33

# stations
stat_list <- get_stations(lt_vec)


################################################################################
##### univariate weight estimation (fixed window)

experiment1 <- function(kernel, lt_vec = 1:33, stat_ids = stat_list, clim = FALSE, mbm = FALSE) {
  w_lp <- array(NA, c(length(lt_vec), length(stat_ids), 3))
  w_wtd <- w_ord <- array(NA, c(length(lt_vec), length(stat_ids), 83))
  if (clim) {
    w_lp_cl <- array(NA, c(length(lt_vec), length(stat_ids), 4))
    w_wtd_cl <- w_ord_cl <- array(NA, c(length(lt_vec), length(stat_ids), 752))
    crps_mat <- array(NA, c(length(lt_vec), length(stat_ids), 11))
    dimnames(crps_mat)[[3]] <- c("Clim.", "C1", "C2", "IFS", "MM", "LP", "Wtd", 
                                 "Wtd-Ord", "LP+Clim", "Wtd+Clim", "Wtd-Ord+Clim")
  } else {
    crps_mat <- array(NA, c(length(lt_vec), length(stat_ids), 8))
    dimnames(crps_mat)[[3]] <- c("Clim.", "C1", "C2", "IFS", "MM", "LP", "Wtd", "Wtd-Ord")
  }

  
  for (lt in seq_along(lt_vec)) {
    lead <- lt_vec[lt]
    
    if (lt == 1) start <- Sys.time()
    
    dat <- load_data(lead, stat_id = NULL) # load dat
    dat <- dat %>% # restrict attention to stations with complete set of observations
      group_by(nat_abbr) %>%
      filter(n() == 1030) %>% 
      ungroup()
    
    for (i in seq_along(stat_ids)) {
      
      print(c(lead, i))
      
      tr_dat <- dat %>% filter(nat_abbr == stat_ids[i], reftime < "2022-06-01")
      ts_dat <- dat %>% filter(nat_abbr == stat_ids[i], reftime >= "2022-06-01")

      if (mbm) {
        scores <- get_mbm_scores(tr_dat, ts_dat, kernel, clim)
      } else {
        scores <- get_scores(tr_dat, ts_dat, kernel, clim)
      }
      
      crps_mat[lt, i, ] <- scores$crps
      w_lp[lt, i, ] <- scores$w$lp
      w_wtd[lt, i, ] <- scores$w$wtd
      w_ord[lt, i, ] <- scores$w$ord
      if (clim) {
        w_lp_cl[lt, i, ] <- scores$w$lp_cl
        w_wtd_cl[lt, i, ] <- scores$w$wtd_cl
        w_ord_cl[lt, i, ] <- scores$w$ord_cl
      }
      
    }

    if (lt == 1) print(Sys.time() - start)
  }
  
  if (clim) {
    return(list(crps = crps_mat, 
                w = list(lp = w_lp, lp_cl = w_lp_cl, wtd = w_wtd, 
                         wtd_cl = w_wtd_cl, ord = w_ord, ord_cl = w_ord_cl)))
  } else {
    return(list(crps = crps_mat, w = list(lp = w_lp, wtd = w_wtd, ord = w_ord)))
  }

}

results1 <- experiment1(kernel = "Energy", lt_vec = 1:33)
results1_mbm <- experiment1(kernel = "Energy", lt_vec = 1:33, mbm = TRUE)


################################################################################
##### univariate weight estimation (rolling window)

experiment2 <- function(kernel, lt_vec = 1:33, win_len = 45, stat_ids = stat_list, clim = FALSE, mbm = FALSE) {
  w_lp <- array(NA, c(length(lt_vec), length(stat_ids), 3, 730))
  w_wtd <- w_ord <- array(NA, c(length(lt_vec), length(stat_ids), 83, 730))
  if (clim) {
    w_lp_cl <- array(NA, c(length(lt_vec), length(stat_ids), 4, 730))
    w_wtd_cl <- w_ord_cl <- array(NA, c(length(lt_vec), length(stat_ids), 83 + win_len, 730))
    crps_mat <- array(NA, c(length(lt_vec), length(stat_ids), 11, 730))
    dimnames(crps_mat)[[3]] <- c("Clim.", "C1", "C2", "IFS", "MM", "LP", "Wtd", 
                                 "Wtd-Ord", "LP+Clim", "Wtd+Clim", "Wtd-Ord+Clim")
  } else {
    crps_mat <- array(NA, c(length(lt_vec), length(stat_ids), 8, 730))
    dimnames(crps_mat)[[3]] <- c("Clim.", "C1", "C2", "IFS", "MM", "LP", "Wtd", "Wtd-Ord")
  }
  

  
  for (lt in seq_along(lt_vec)) {
    lead <- lt_vec[lt]
  
    if (lt == 1) start <- Sys.time()
    
    dat <- load_data(lead, stat_id = NULL) # load dat
    dat <- dat %>% # restrict attention to stations with complete set of observations
      group_by(nat_abbr) %>%
      filter(n() == 1030) %>% 
      ungroup()

    date_list <- unique(dat$reftime)  # some dates are missing - these are the same for all lead times
    date_list <- sort(date_list) # dates are not chronological
    date_list_ts <- seq(as.Date("2021-06-01"), by = "days", length.out = 730)
    
    for (i in seq_along(stat_ids)) {
      
      for (j in seq_along(date_list_ts)) {
        
        print(c(lead, i, j))
        
        ele <- date_list_ts[j]
        
        if (!(ele %in% date_list)) next
        
        tr_dates <- tail(sort(date_list[date_list < ele]), win_len)
        if (sum(tr_dates %in% date_list) < win_len) {
          next
        }
        
        tr_dat <- dat %>% filter(nat_abbr == stat_ids[i], reftime %in% tr_dates)
        ts_dat <- dat %>% filter(nat_abbr == stat_ids[i], reftime == ele)
        
        if (mbm) {
          scores <- get_mbm_scores(tr_dat, ts_dat, kernel, clim)
        } else {
          scores <- get_scores(tr_dat, ts_dat, kernel, clim)
        }
        
        crps_mat[lt, i, , j] <- scores$crps
        w_lp[lt, i, , j] <- scores$w$lp
        w_wtd[lt, i, , j] <- scores$w$wtd
        w_ord[lt, i, , j] <- scores$w$ord
        if (clim) {
          w_lp_cl[lt, i, , j] <- scores$w$lp_cl
          w_wtd_cl[lt, i, , j] <- scores$w$wtd_cl
          w_ord_cl[lt, i, , j] <- scores$w$ord_cl
        }
        
      }
      
    }
    
    print(Sys.time() - start)
  
  }
  
  if (clim) {
    return(list(crps = crps_mat, 
                w = list(lp = w_lp, lp_cl = w_lp_cl, wtd = w_wtd, 
                         wtd_cl = w_wtd_cl, ord = w_ord, ord_cl = w_ord_cl)))
  } else {
    return(list(crps = crps_mat, w = list(lp = w_lp, wtd = w_wtd, ord = w_ord)))
  }
  
}

results2 <- experiment2(kernel = "Energy", lt_vec = 18)
results2_mbm <- experiment2(kernel = "Energy", lt_vec = 18, mbm = TRUE)


################################################################################
##### multivariate weight estimation (fixed window)

experiment3 <- function(kernel, lt_vec = 1:33, stat_ids = stat_list, clim = FALSE, mbm = FALSE) {
  w_lp <- array(NA, c(length(lt_vec), 3))
  w_wtd <- array(NA, c(length(lt_vec), 83))
  if (clim) {
    w_lp_cl <- array(NA, c(length(lt_vec), 4))
    w_wtd_cl <- array(NA, c(length(lt_vec), 752))
    crps_mat <- array(NA, c(length(lt_vec), length(stat_ids), 9))
    es_mat <- array(NA, c(length(lt_vec), 9))
    dimnames(crps_mat)[[3]] <- colnames(es_mat) <- 
      c("Clim.", "C1", "C2", "IFS", "MM", "LP", "Wtd", "LP+Clim", "Wtd+Clim")
  } else {
    crps_mat <- array(NA, c(length(lt_vec), length(stat_ids), 7))
    es_mat <- array(NA, c(length(lt_vec), 7))
    dimnames(crps_mat)[[3]] <- colnames(es_mat) <-
      c("Clim.", "C1", "C2", "IFS", "MM", "LP", "Wtd")
  }
  
  for (lt in seq_along(lt_vec)) {
    lead <- lt_vec[lt]
    
    print(lt)
    
    if (lt == 1) start <- Sys.time()
    
    dat <- load_data(lead, stat_id = NULL) # load dat
    dat <- dat %>% # restrict attention to stations with complete set of observations
      group_by(nat_abbr) %>%
      filter(n() == 1030) %>%
      ungroup()
    
    tr_dat <- dat %>% filter(reftime < "2022-06-01", nat_abbr %in% stat_ids)
    ts_dat <- dat %>% filter(reftime >= "2022-06-01", nat_abbr %in% stat_ids)
    rm(dat)
    
    if (mbm) {
      scores <- get_mv_mbm_scores(tr_dat, ts_dat, stat_ids, kernel, clim)
    } else {
      scores <- get_mv_scores(tr_dat, ts_dat, stat_ids, kernel, clim)
    }
    
    crps_mat[lt, , ] <- scores$crps
    es_mat[lt, ] <- scores$es
    w_lp[lt, ] <- scores$w$lp
    w_wtd[lt, ] <- scores$w$wtd
    if (clim) {
      w_lp_cl[lt, ] <- scores$w$lp_cl
      w_wtd_cl[lt, ] <- scores$w$wtd_cl
    }
    
    if (lt == 1) print(Sys.time() - start)
    
  }

  
  if (clim) {
    return(list(crps = crps_mat, es = es_mat,
                w = list(lp = w_lp, lp_cl = w_lp_cl, wtd = w_wtd, 
                         wtd_cl = w_wtd_cl, ord_cl = w_ord_cl)))
  } else {
    return(list(crps = crps_mat, es = es_mat, w = list(lp = w_lp, wtd = w_wtd)))
  }
  
}

results3 <- experiment3(kernel = "Energy", lt_vec = 1:33)
results3_mbm <- experiment3(kernel = "Energy", lt_vec = 1:33, mbm = TRUE)


################################################################################
##### multivariate weight estimation (rolling window)

experiment4 <- function(kernel, lt_vec = 1:33, win_len = 45, stat_ids = stat_list, clim = FALSE, mbm = FALSE) {
  w_lp <- array(NA, c(length(lt_vec), 3, 730))
  w_wtd <- array(NA, c(length(lt_vec), 83, 730))
  if (clim) {
    w_lp_cl <- array(NA, c(length(lt_vec), 4, 730))
    w_wtd_cl <- array(NA, c(length(lt_vec), 83 + win_len, 730))
    crps_mat <- array(NA, c(length(lt_vec), length(stat_ids), 9, 730))
    es_mat <- array(NA, c(length(lt_vec), 9, 730))
    dimnames(crps_mat)[[3]] <- colnames(es_mat) <- 
      c("Clim.", "C1", "C2", "IFS", "MM", "LP", "Wtd", "LP+Clim", "Wtd+Clim")
  } else {
    crps_mat <- array(NA, c(length(lt_vec), length(stat_ids), 7, 730))
    es_mat <- array(NA, c(length(lt_vec), 7, 730))
    dimnames(crps_mat)[[3]] <- colnames(es_mat) <- 
      c("Clim.", "C1", "C2", "IFS", "MM", "LP", "Wtd")
  }
  
  for (lt in seq_along(lt_vec)) {
    lead <- lt_vec[lt]
    
    if (lt == 1) start <- Sys.time()
    
    dat <- load_data(lead, stat_id = NULL) # load dat
    dat <- dat %>% # restrict attention to stations with complete set of observations
      group_by(nat_abbr) %>%
      filter(n() == 1030) %>% 
      ungroup()
    
    date_list <- unique(dat$reftime)
    date_list <- sort(date_list) # dates are not chronological
    date_list_ts <- seq(as.Date("2021-06-01"), by = "days", length.out = 730)
    
    for (j in seq_along(date_list_ts)) {
      
      print(c(lead, j))
      
      ele <- date_list_ts[j]
      
      if (!(ele %in% date_list)) next
      
      tr_dates <- tail(sort(date_list[date_list < ele]), win_len)
      if (sum(tr_dates %in% date_list) < win_len) {
        next
      }
      
      tr_dat <- dat %>% filter(reftime %in% tr_dates, nat_abbr %in% stat_ids)
      ts_dat <- dat %>% filter(reftime == ele, nat_abbr %in% stat_ids)
      
      if (mbm) {
        scores <- get_mv_mbm_scores(tr_dat, ts_dat, stat_ids, kernel, clim)
      } else {
        scores <- get_mv_scores(tr_dat, ts_dat, stat_ids, kernel, clim)
      }
      
      crps_mat[lt, , , j] <- scores$crps
      es_mat[lt, , j] <- scores$es
      w_lp[lt, , j] <- scores$w$lp
      w_wtd[lt, , j] <- scores$w$wtd
      if (clim) {
        w_lp_cl[lt, , j] <- scores$w$lp_cl
        w_wtd_cl[lt, , j] <- scores$w$wtd_cl
      }
      
    }
    
    print(Sys.time() - start)
    
  }
  
  if (clim) {
    return(list(crps = crps_mat, 
                w = list(lp = w_lp, lp_cl = w_lp_cl, wtd = w_wtd, wtd_cl = w_wtd_cl)))
  } else {
    return(list(crps = crps_mat, w = list(lp = w_lp, wtd = w_wtd)))
  }
  
}

results4 <- experiment4(kernel = "Energy", lt_vec = 18)
results4_mbm <- experiment4(kernel = "Energy", lt_vec = 18, mbm = TRUE)


################################################################################
##### evaluate

save(results1, results1_mbm, results3, results3_mbm, 
     file = "C:/Users/sa20i493/Documents/R/KernelEmbeddings/Data/results_data.RData")


### plot stations

library(terra)
library(viridis)
library(rnaturalearth)
library(sf)  # For handling vector data

SMN_stations <- readRDS("C:/Users/sa20i493/Documents/R/WeightedLoss/WindGusts/SMN_stations.rds")
SMN_stations <- SMN_stations %>% filter(nat_abbr %in% stat_list)

plot_map <- function(lons, lats, z, title = NULL){
  if (is.matrix(z)) z <- as.vector(z)
  
  ## elevation data
  #dem <- elevation_30s(country = "CHE", path = tempdir())
  dem <- elevation_global(res = 0.5, path = tempdir())
  
  xmin <- 5.8  # Western edge
  xmax <- 10.5 # Eastern edge
  ymin <- 45.8 # Southern edge
  ymax <- 47.9 # Northern edge
  
  region_extent <- terra::ext(xmin, xmax, ymin, ymax)
  cropped_dem <- terra::crop(dem, region_extent)
  dem_df <- as.data.frame(cropped_dem, xy = TRUE)
  
  colnames(dem_df) <- c("lon", "lat", "elevation")
  
  ## boundary data
  switzerland <- ne_countries(scale = "large", country = "Switzerland", returnclass = "sf")

  ## plot data
  df <- data.frame(lat = lats, lon = lons, z = z)
  
  # Plot the elevation data
  plot_obj <- 
    ggplot() +
    geom_raster(data = dem_df, aes(x = lon, y = lat, fill = elevation)) +
    geom_sf(data = switzerland, fill = NA, color = "black") +
    geom_point(data = df, aes(lon, lat, fill = z), shape = 21) +
    scale_x_continuous(name = "Longitude", expand = c(0, 0)) +
    scale_y_continuous(name = "Latitude", expand = c(0, 0)) +
    scale_fill_gradient(name = "Elevation (m)", low = "white", high = "grey40") +
    theme_minimal() +
    theme(panel.border = element_rect(color = "black", fill = NA))
  
  
  #world <- map_data("world")
  #ind <- (world$long >= xmin & world$long <= xmax) & (world$lat >= ymin & world$lat <= ymax)
  #world <- world[ind, ]
  
  #plot_obj <- ggplot() + borders("world") +
    #coord_fixed(ylim = range(lats), xlim = range(lons)) +
    #geom_point(data = df, aes(lon, lat, fill = z), shape = 21) +
    #scale_fill_gradient(low = "white", high = "red", limits = c(ymin, ymax),
    #                   guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
    #theme_void() + theme(legend.title = element_blank(),
    #                     panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    #                     legend.key.width = unit(0.3, "in")) +
    #ggtitle(title)
  
  return(plot_obj)
  
}

plot_map(SMN_stations$longitude, SMN_stations$latitude, SMN_stations$station_height)
ggsave("Figures/cs_alt_map.png", height = 3, width = 5)


### plot weights vs lead time

w <- apply(results1$w$wtd, c(1, 3), mean)
w <- cbind(rowSums(w[, 1:11]), rowSums(w[, 12:32]), rowSums(w[, 33:83]))
df <- data.frame(lt = 1:33, 
                 w = as.vector(w),
                 mth = rep(c("COSMO-1E", "COSMO-2E", "ECMWF IFS"), each = 33))
ggplot(df) + geom_line(aes(x = lt, y = w, col = mth, linetype = mth)) +
  scale_x_continuous(name = "Lead time (hours)", expand = c(0, 0)) +
  scale_y_continuous(name = "Weight", limits = c(0, 1), expand = c(0, 0)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.justification = c(0.5, 1),
        legend.position = c(0.5, 0.99),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        legend.direction="horizontal")
ggsave("Figures/cs_weight_vs_lead.png", height = 2.5, width = 4)


### plot weights vs altitude

lt <- 18

w <- rowSums(results1$w$wtd[lt, , 1:11])
w <- w[stat_list %in% SMN_stations$nat_abbr]
plot_map(SMN_stations$longitude, SMN_stations$latitude, w, ymax = 1, title = "COSMO-1E")
ggsave("Figures/cs_weight_vs_stat_c1.png", height = 1.7, width = 3.4)

w <- rowSums(results1$w$wtd[lt, , 22:32])
w <- w[stat_list %in% SMN_stations$nat_abbr]
plot_map(SMN_stations$longitude, SMN_stations$latitude, w, ymax = 1, title = "COSMO-2E")
ggsave("Figures/cs_weight_vs_stat_c2.png", height = 1.7, width = 3.4)

w <- rowSums(results1$w$wtd[lt, , 33:83])
w <- w[stat_list %in% SMN_stations$nat_abbr]
plot_map(SMN_stations$longitude, SMN_stations$latitude, w, ymax = 1, title = "ECMWF IFS")
ggsave("Figures/cs_weight_vs_stat_ifs.png", height = 1.7, width = 3.4)


# plot station with highest weight at each location

w <- cbind(rowSums(results1$w$wtd[lt, , 1:11]), 
           rowSums(results1$w$wtd[lt, , 22:32]), 
           rowSums(results1$w$wtd[lt, , 33:83]))
w <- w[stat_list %in% SMN_stations$nat_abbr, ]
plot_map2 <- function(lons, lats, z, title = NULL){
  
  ## elevation data
  #dem <- elevation_30s(country = "CHE", path = tempdir())
  dem <- elevation_global(res = 0.5, path = tempdir())
  
  xmin <- 5.8  # Western edge
  xmax <- 10.5 # Eastern edge
  ymin <- 45.8 # Southern edge
  ymax <- 47.9 # Northern edge
  
  region_extent <- terra::ext(xmin, xmax, ymin, ymax)
  cropped_dem <- terra::crop(dem, region_extent)
  dem_df <- as.data.frame(cropped_dem, xy = TRUE)
  
  colnames(dem_df) <- c("lon", "lat", "elevation")
  
  ## boundary data
  switzerland <- ne_countries(scale = "large", country = "Switzerland", returnclass = "sf")
  
  ## plot data
  z <- apply(z, 1, which.max)
  z[z == 1] <- c("COSMO-1E")
  z[z == 2] <- c("COSMO-2E")
  z[z == 3] <- c("ECMWF IFS")
  df <- data.frame(lat = lats, lon = lons, Model = z)
  
  # Plot the elevation data
  plot_obj <- 
    ggplot() +
    geom_raster(data = dem_df, aes(x = lon, y = lat, fill = elevation)) +
    geom_sf(data = switzerland, fill = NA, color = "black") +
    geom_point(data = df, aes(lon, lat, shape = Model, color = Model), size = 2) +
    scale_x_continuous(name = "Longitude", expand = c(0, 0)) +
    scale_y_continuous(name = "Latitude", expand = c(0, 0)) +
    scale_fill_gradient(name = "Elevation (m)", low = "white", high = "grey40") +
    theme_minimal() +
    theme(panel.border = element_rect(color = "black", fill = NA))
  
  
  #world <- map_data("world")
  #ind <- (world$long >= xmin & world$long <= xmax) & (world$lat >= ymin & world$lat <= ymax)
  #world <- world[ind, ]
  
  #plot_obj <- ggplot() + borders("world") +
  #coord_fixed(ylim = range(lats), xlim = range(lons)) +
  #geom_point(data = df, aes(lon, lat, fill = z), shape = 21) +
  #scale_fill_gradient(low = "white", high = "red", limits = c(ymin, ymax),
  #                   guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  #theme_void() + theme(legend.title = element_blank(),
  #                     panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
  #                     legend.key.width = unit(0.3, "in")) +
  #ggtitle(title)
  
  return(plot_obj)
  
}
plot_map2(SMN_stations$longitude, SMN_stations$latitude, w)
ggsave("Figures/cs_weight_vs_stat.png", height = 3*1.2, width = 5*1.2)


### plot weights vs mse

lt <- 18
dat <- load_data(lt, stat_id = NULL) # load dat
dat <- dat %>% # restrict attention to stations with complete set of observations
  group_by(nat_abbr) %>%
  filter(n() == 1030) %>% 
  ungroup() %>% filter(reftime >= "2022-06-01", nat_abbr %in% stat_list)

mse <- c(colMeans((dat$`COSMO-1E` - dat$obs)^2),
         colMeans((dat$`COSMO-2E` - dat$obs)^2),
         colMeans((dat$ECMWF_IFS - dat$obs)^2))

w <- colMeans(results1$w$wtd[lt, , ])
df <- data.frame(s = mse, w = w, 
                 mth = c(rep("COSMO-1E", 11), rep("COSMO-2E", 21), rep("ECMWF IFS", 51)),
                 control = c(1, rep(0, 10), 1, rep(0, 20), 1, rep(0, 50)))
ggplot(df) + geom_point(aes(x = s, y = w, col = mth, shape = as.factor(control))) +
  scale_x_continuous(name = "MSE") +
  scale_y_continuous(name = "Weight", limits = c(0, 0.07)) +
  scale_shape_manual(values = c(19, 4)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.justification = c(0.5, 1),
        legend.position = c(0.5, 0.99),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        legend.direction="horizontal") +
  guides(shape = "none")
ggsave("Figures/cs_weight_vs_mse.png", height = 2.5, width = 4)


### role of the control member

# COSMO-1E 
x <- as.matrix(dat$`COSMO-1E`)
M <- ncol(x)
crps_dsum <- rowSums(sapply(1:M, function(i) rowSums(abs(x - x[, i]))))
crps_wc <- mean(rowMeans(abs(x - dat$obs)) -  crps_dsum/ (2 * M * (M - 1)))
x <- x[, -1]
M <- ncol(x)
crps_dsum <- rowSums(sapply(1:M, function(i) rowSums(abs(x - x[, i]))))
crps_woc <- mean(rowMeans(abs(x - dat$obs)) -  crps_dsum/ (2 * M * (M - 1)))
1 - crps_woc/crps_wc

# COSMO-2E
x <- as.matrix(dat$`COSMO-2E`)
M <- ncol(x)
crps_dsum <- rowSums(sapply(1:M, function(i) rowSums(abs(x - x[, i]))))
crps_wc <- mean(rowMeans(abs(x - dat$obs)) -  crps_dsum/ (2 * M * (M - 1)))
x <- x[, -1]
M <- ncol(x)
crps_dsum <- rowSums(sapply(1:M, function(i) rowSums(abs(x - x[, i]))))
crps_woc <- mean(rowMeans(abs(x - dat$obs)) -  crps_dsum/ (2 * M * (M - 1)))
1 - crps_woc/crps_wc

# ECMWF IFS
x <- as.matrix(dat$ECMWF_IFS)
M <- ncol(x)
crps_dsum <- rowSums(sapply(1:M, function(i) rowSums(abs(x - x[, i]))))
crps_wc <- mean(rowMeans(abs(x - dat$obs)) -  crps_dsum/ (2 * M * (M - 1)))
x <- x[, -1]
M <- ncol(x)
crps_dsum <- rowSums(sapply(1:M, function(i) rowSums(abs(x - x[, i]))))
crps_woc <- mean(rowMeans(abs(x - dat$obs)) -  crps_dsum/ (2 * M * (M - 1)))
1 - crps_woc/crps_wc

# Multi-model
x <- cbind(as.matrix(dat$`COSMO-1E`), as.matrix(dat$`COSMO-2E`), as.matrix(dat$ECMWF_IFS))
M <- ncol(x)
crps_dsum <- rowSums(sapply(1:M, function(i) rowSums(abs(x - x[, i]))))
crps_wc <- mean(rowMeans(abs(x - dat$obs)) -  crps_dsum/ (2 * M * (M - 1)))
x <- x[, -c(1, 12, 33)]
M <- ncol(x)
crps_dsum <- rowSums(sapply(1:M, function(i) rowSums(abs(x - x[, i]))))
crps_woc <- mean(rowMeans(abs(x - dat$obs)) -  crps_dsum/ (2 * M * (M - 1)))
1 - crps_woc/crps_wc


### plot mv weights vs lead time

w <- results3$w$wtd
w <- cbind(rowSums(w[, 1:11]), rowSums(w[, 12:32]), rowSums(w[, 33:83]))
df <- data.frame(lt = 1:33, 
                 w = as.vector(w),
                 mth = rep(c("COSMO-1E", "COSMO-2E", "ECMWF IFS"), each = 33))
ggplot(df) + geom_line(aes(x = lt, y = w, col = mth, linetype = mth)) +
  scale_x_continuous(name = "Lead time (hours)", expand = c(0, 0)) +
  scale_y_continuous(name = "Weight", limits = c(0, 1), expand = c(0, 0)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.justification = c(0.5, 1),
        legend.position = c(0.5, 0.99),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        legend.direction="horizontal")
ggsave("Figures/cs_mvweight_vs_lead.png", height = 2.5, width = 4)


### plot mv weights vs mse

lt <- 18
dat <- load_data(lt, stat_id = NULL) # load dat
dat <- dat %>% # restrict attention to stations with complete set of observations
  group_by(nat_abbr) %>%
  filter(n() == 1030) %>% 
  ungroup() %>% filter(reftime >= "2022-06-01", nat_abbr %in% stat_list)

y <- dat %>% dplyr::select(nat_abbr, obs) %>% 
  group_by(nat_abbr) %>%
  mutate(row = row_number()) %>%
  pivot_wider(names_from = nat_abbr, values_from = obs) %>%
  dplyr::select(-row) %>% as.matrix()

x <- lapply(stat_list, function(z) {
  dat %>% filter(nat_abbr == z) %>% 
    dplyr::select(`COSMO-1E`, `COSMO-2E`, ECMWF_IFS) %>%
    cbind() %>% as.matrix()
}) %>% simplify2array() %>% aperm(c(1, 3, 2))

mse <- colMeans(sqrt(apply((x - replicate(83, y))^2, c(1, 3), sum)))


w <- results3$w$wtd[lt, ]
df <- data.frame(s = mse, w = w, 
                 mth = c(rep("COSMO-1E", 11), rep("COSMO-2E", 21), rep("ECMWF IFS", 51)),
                 control = c(1, rep(0, 10), 1, rep(0, 20), 1, rep(0, 50)))
ggplot(df) + geom_point(aes(x = s, y = w, col = mth, shape = as.factor(control))) +
  scale_x_continuous(name = "MSE") +
  scale_y_continuous(name = "Weight", limits = c(0, 0.1)) +
  scale_shape_manual(values = c(19, 4)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.justification = c(0.5, 1),
        legend.position = c(0.5, 0.99),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        legend.direction="horizontal") +
  guides(shape = "none")
ggsave("Figures/cs_mvweight_vs_mse.png", height = 5*0.51, width = 8*0.51)


### plot crps vs lead time

s <- apply(results1$crps, c(1, 3), mean)
s <- s[, -c(1, 8)]
colnames(s) <- c("COSMO-1E", "COSMO-2E", "ECMWF IFS", "LP  Equal", "LP Discrete", "LP Point")

df <- data.frame(lt = 1:33, s = as.vector(s), mth = rep(colnames(s), each = 33))
ggplot(df) + geom_line(aes(x = lt, y = s, col = mth, linetype = mth)) +
  scale_x_continuous(name = "Lead time (hours)", expand = c(0, 0)) +
  scale_y_continuous(name = "CRPS (m/s)", limits = c(0.5, 2)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.justification = c(0.5, 1),
        legend.position = c(0.5, 0.99),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA)) +
  guides(col = guide_legend(nrow = 2, byrow = TRUE))
ggsave("Figures/cs_crps_vs_lead.png", height = 5*0.5, width = 8*0.5)
range(1 - s[, 5]/s[, 1])
range(1 - s[, 6]/s[, 5])


### plot crps of order-statistics approach

s <- apply(results1$crps, c(1, 3), mean)
s <- s[, c(5, 6, 7, 8)]
colnames(s) <- c("LP  Equal", "LP Discrete", "LP Point", "LP Ordered")

df <- data.frame(lt = 1:33, s = as.vector(s), mth = rep(c("Linear pool", "Weighted", "Weighted (OS)"), each = 33))
ggplot(df) + geom_line(aes(x = lt, y = s, col = mth)) +
  geom_hline(aes(yintercept = 0), lty = "dotted") +
  scale_x_continuous(name = "Lead time (hours)", expand = c(0, 0)) +
  scale_y_continuous(name = "CRPSS", limits = c(-0.05, 0.5)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.justification = c(0.5, 1),
        legend.position = c(0.5, 0.99),
        legend.direction="horizontal")

df <- data.frame(lt = 1:33, s = as.vector(s), mth = rep(colnames(s), each = 33))
ggplot(df) + geom_line(aes(x = lt, y = s, col = mth, linetype = mth)) +
  scale_x_continuous(name = "Lead time (hours)", expand = c(0, 0)) +
  scale_y_continuous(name = "CRPS (m/s)", limits = c(0.5, 1.3)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.justification = c(0.5, 1),
        legend.position = c(0.5, 0.99),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA)) +
  guides(col = guide_legend(nrow = 2, byrow = TRUE))
ggsave("Figures/cs_wcrpss_vs_lead.png", height = 5*0.5, width = 8*0.5)


### plot es vs lead time

s <- results3$es
s <- s[, -1]
colnames(s) <- c("COSMO-1E", "COSMO-2E", "ECMWF IFS", "LP  Equal", "LP Discrete", "LP Point")

df <- data.frame(lt = 1:33, s = as.vector(s), mth = rep(colnames(s), each = 33))
ggplot(df) + geom_line(aes(x = lt, y = s, col = mth, linetype = mth)) +
  scale_x_continuous(name = "Lead time (hours)", expand = c(0, 0)) +
  scale_y_continuous(name = "Energy score", limits = c(8, 22)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.justification = c(0.5, 1),
        legend.position = c(0.5, 0.99),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA)) +
  guides(col = guide_legend(nrow = 2, byrow = TRUE))
ggsave("Figures/cs_es_vs_lead.png", height = 5*0.5, width = 8*0.5)
1 - s[, 5]/s[, 1]
1 - s[, 6]/s[, 5]



### plot weights of order statistics

lt <- 18
#w <- colMeans(results1$w$ord[lt, , 1:11])
#w <- colMeans(results1$w$ord[lt, , 12:32])
w <- colMeans(results1$w$ord[lt, , 33:83])
M <- length(w)
df <- data.frame(ord = 1:M, w = w)
ggplot(df) + geom_bar(aes(x = ord, y = w), stat = "identity") + 
  #scale_x_continuous(name = "Order statistic", breaks = 1:M, labels = 1:M, expand = c(0, 0)) +
  #scale_x_continuous(name = "Order statistic", breaks = c(1, 5, 10, 15, 21), 
  #                   labels = c(1, 5, 10, 15, 21), expand = c(0, 0)) +
  scale_x_continuous(name = "Order statistic", breaks = c(1, seq(10, 40, 10), 51), 
                     labels = c(1, seq(10, 40, 10), 51), expand = c(0, 0)) +
  scale_y_continuous(name = "Weight", limits = c(0, 0.1), expand = c(0, 0)) +
  ggtitle("ECMWF IFS") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 12))
ggsave("Figures/cs_weight_vs_ord_ifs.png", height = 5*0.4 + 0.3, width = 8*0.4)


### plot average cdf before and after weighting

x <- cbind(as.matrix(dat$`COSMO-1E`), as.matrix(dat$`COSMO-2E`), as.matrix(dat$ECMWF_IFS))
xx <- seq(0, 11, 0.01)  

F_mm <- sapply(xx, function(z) rowMeans(x <= z))
F_mm <- colMeans(F_mm)

w <- results1$w$lp[lt, , ]
w <- cbind(replicate(11, w[, 1]/11), replicate(21, w[, 2]/21), replicate(51, w[, 3]/51))
w <- lapply(1:nrow(w), function(i) replicate(361, w[i, ]) %>% t()) %>% do.call(rbind, .)
F_lp <- sapply(xx, function(z) rowSums(w*(x <= z)))
F_lp <- colMeans(F_lp)

x1 <- as.matrix(dat$`COSMO-1E`)%>% apply(1, sort) %>% t()
x2 <- as.matrix(dat$`COSMO-2E`)%>% apply(1, sort) %>% t()
x3 <- as.matrix(dat$ECMWF_IFS)%>% apply(1, sort) %>% t()
x <- cbind(x1, x2, x3)
w <- results1$w$ord[lt, , ]
w <- lapply(1:nrow(w), function(i) replicate(361, w[i, ]) %>% t()) %>% do.call(rbind, .)
F_ord <- sapply(xx, function(z) rowSums(w*(x <= z)))
F_ord <- colMeans(F_ord)


df <- data.frame(x = xx, y = c(F_mm, F_lp, F_ord), 
                 mth = rep(c("Multi-model", "Linear pool", "Weighted (OS)"), each = length(xx)))
ggplot(df) + geom_line(aes(x = x, y = y, col = mth)) +
  scale_x_continuous(name = "x", limits = c(0, 11), breaks = seq(0, 10, 2), expand = c(0, 0)) +
  scale_y_continuous(name = "F(x)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.justification = c(1, 0),
        legend.position = c(0.99, 0.01))
ggsave("Figures/cs_average_dist.png", height = 5*0.5, width = 8*0.5)



### PIT histograms

lt <- 18
dat <- load_data(lt, stat_id = NULL) # load dat
dat <- dat %>% # restrict attention to stations with complete set of observations
  group_by(nat_abbr) %>%
  filter(n() == 1030) %>% 
  ungroup() %>% filter(reftime >= "2022-06-01", nat_abbr %in% stat_list)

y <- dat$obs

# COSMO-1E
x <- as.matrix(dat$`COSMO-1E`)
pit <- rowMeans(x < y) + runif(length(y))*(rowMeans(x <= y) - rowMeans(x < y))
pit_hist(pit, ranks = FALSE, ymax = 0.6, xticks = FALSE, xlab = NULL, yticks = FALSE, ylab = NULL, title = "COSMO-1E") +
  theme(plot.title = element_text(size = 11))
ggsave("Figures/cs_pit_c1.png", height = 7*0.25, width = 8*0.25)

# COSMO-2E
x <- as.matrix(dat$`COSMO-2E`)
pit <- rowMeans(x < y) + runif(length(y))*(rowMeans(x <= y) - rowMeans(x < y))
pit_hist(pit, ranks = FALSE, ymax = 0.6, xticks = FALSE, xlab = NULL, yticks = FALSE, ylab = NULL, title = "COSMO-2E") +
  theme(plot.title = element_text(size = 11))
ggsave("Figures/cs_pit_c2.png", height = 7*0.25, width = 8*0.25)

# IFS
x <- as.matrix(dat$ECMWF_IFS)
pit <- rowMeans(x < y) + runif(length(y))*(rowMeans(x <= y) - rowMeans(x < y))
pit_hist(pit, ranks = FALSE, ymax = 0.6, xticks = FALSE, xlab = NULL, yticks = FALSE, ylab = NULL, title = "ECMWF IFS") +
  theme(plot.title = element_text(size = 11))
ggsave("Figures/cs_pit_ifs.png", height = 7*0.25, width = 8*0.25)

# Multi-model
x <- cbind(as.matrix(dat$`COSMO-1E`), as.matrix(dat$`COSMO-2E`), as.matrix(dat$ECMWF_IFS))
pit <- rowMeans(x < y) + runif(length(y))*(rowMeans(x <= y) - rowMeans(x < y))
pit_hist(pit, ranks = FALSE, ymax = 0.6, xticks = FALSE, xlab = NULL, yticks = FALSE, ylab = NULL, title = "LP Equal") +
  theme(plot.title = element_text(size = 11))
ggsave("Figures/cs_pit_mm.png", height = 7*0.25, width = 8*0.25)

# Linear pool
w <- results1$w$lp[lt, , ]
w <- cbind(replicate(11, w[, 1]/11), replicate(21, w[, 2]/21), replicate(51, w[, 3]/51))
w <- lapply(1:nrow(w), function(i) replicate(361, w[i, ]) %>% t()) %>% do.call(rbind, .)
pit <- rowSums(w*(x < y)) + runif(length(y))*(rowMeans(w*(x <= y)) - rowMeans(w*(x < y)))
pit_hist(pit, ranks = FALSE, ymax = 0.6, xticks = FALSE, xlab = NULL, yticks = FALSE, ylab = NULL, title = "LP Discrete") +
  theme(plot.title = element_text(size = 11))
ggsave("Figures/cs_pit_lp.png", height = 7*0.25, width = 8*0.25)

# Weighted
w <- results1$w$wtd[lt, , ]
w <- lapply(1:nrow(w), function(i) replicate(361, w[i, ]) %>% t()) %>% do.call(rbind, .)
pit <- rowSums(w*(x < y)) + runif(length(y))*(rowMeans(w*(x <= y)) - rowMeans(w*(x < y)))
pit_hist(pit, ranks = FALSE, ymax = 0.6, xticks = FALSE, xlab = NULL, yticks = FALSE, ylab = NULL, title = "LP Point") +
  theme(plot.title = element_text(size = 11))
ggsave("Figures/cs_pit_wtd.png", height = 7*0.25, width = 8*0.25)

# Weighted
x1 <- as.matrix(dat$`COSMO-1E`)%>% apply(1, sort) %>% t()
x2 <- as.matrix(dat$`COSMO-2E`)%>% apply(1, sort) %>% t()
x3 <- as.matrix(dat$ECMWF_IFS)%>% apply(1, sort) %>% t()
x <- cbind(x1, x2, x3)
w <- results1$w$ord[lt, , ]
w <- lapply(1:nrow(w), function(i) replicate(361, w[i, ]) %>% t()) %>% do.call(rbind, .)
pit <- rowSums(w*(x < y)) + runif(length(y))*(rowMeans(w*(x <= y)) - rowMeans(w*(x < y)))
pit_hist(pit, ranks = FALSE, ymax = 0.6, xticks = FALSE, xlab = NULL, yticks = FALSE, ylab = NULL, title = "LP Ordered") +
  theme(plot.title = element_text(size = 11))
ggsave("Figures/cs_pit_ord.png", height = 7*0.25, width = 8*0.25)




################################################################################
##### evaluate post-processed



### plot weights vs lead time

w <- apply(results1_mbm$w$wtd, c(1, 3), mean)
w <- cbind(rowSums(w[, 1:11]), rowSums(w[, 12:32]), rowSums(w[, 33:83]))
df <- data.frame(lt = 1:33, 
                 w = as.vector(w),
                 mth = rep(c("COSMO-1E", "COSMO-2E", "ECMWF IFS"), each = 33))
ggplot(df) + geom_line(aes(x = lt, y = w, col = mth, linetype = mth)) +
  scale_x_continuous(name = "Lead time (hours)", expand = c(0, 0)) +
  scale_y_continuous(name = "Weight", limits = c(0, 1), expand = c(0, 0)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.justification = c(0.5, 1),
        legend.position = c(0.5, 0.99),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        legend.direction="horizontal")
ggsave("Figures/cs_weight_vs_lead_mbm.png", height = 2.5, width = 4)


### plot weights vs altitude

lt <- 18

w <- rowSums(results1_mbm$w$wtd[lt, , 1:11])
w <- w[stat_list %in% SMN_stations$nat_abbr]
plot_map(SMN_stations$longitude, SMN_stations$latitude, w, ymax = 1, title = "COSMO-1E")
ggsave("Figures/cs_weight_vs_stat_c1_mbm.png", height = 1.7, width = 3.4)

w <- rowSums(results1_mbm$w$wtd[lt, , 22:32])
w <- w[stat_list %in% SMN_stations$nat_abbr]
plot_map(SMN_stations$longitude, SMN_stations$latitude, w, ymax = 1, title = "COSMO-2E")
ggsave("Figures/cs_weight_vs_stat_c2_mbm.png", height = 1.7, width = 3.4)

w <- rowSums(results1_mbm$w$wtd[lt, , 33:83])
w <- w[stat_list %in% SMN_stations$nat_abbr]
plot_map(SMN_stations$longitude, SMN_stations$latitude, w, ymax = 1, title = "ECMWF IFS")
ggsave("Figures/cs_weight_vs_stat_ifs_mbm.png", height = 1.7, width = 3.4)


# plot station with highest weight at each location

w <- cbind(rowSums(results1_mbm$w$wtd[lt, , 1:11]), 
           rowSums(results1_mbm$w$wtd[lt, , 22:32]), 
           rowSums(results1_mbm$w$wtd[lt, , 33:83]))
w <- w[stat_list %in% SMN_stations$nat_abbr, ]
plot_map2(SMN_stations$longitude, SMN_stations$latitude, w)
ggsave("Figures/cs_weight_vs_stat_mbm.png", height = 3*1.2, width = 5*1.2)


### plot mv weights vs lead time

w <- results3_mbm$w$wtd
w <- cbind(rowSums(w[, 1:11]), rowSums(w[, 12:32]), rowSums(w[, 33:83]))
df <- data.frame(lt = 1:33, 
                 w = as.vector(w),
                 mth = rep(c("COSMO-1E", "COSMO-2E", "ECMWF IFS"), each = 33))
ggplot(df) + geom_line(aes(x = lt, y = w, col = mth, linetype = mth)) +
  scale_x_continuous(name = "Lead time (hours)", expand = c(0, 0)) +
  scale_y_continuous(name = "Weight", limits = c(0, 1), expand = c(0, 0)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.justification = c(0.5, 1),
        legend.position = c(0.5, 0.99),        
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        legend.direction="horizontal")
ggsave("Figures/cs_mvweight_vs_lead_mbm.png", height = 2.5, width = 4)



### plot crps vs lead time

s <- apply(results1_mbm$crps, c(1, 3), mean)
s <- s[, -c(1, 8)]
colnames(s) <- c("COSMO-1E", "COSMO-2E", "ECMWF IFS", "LP  Equal", "LP Discrete", "LP Point")

df <- data.frame(lt = 1:33, s = as.vector(s), mth = rep(colnames(s), each = 33))
ggplot(df) + geom_line(aes(x = lt, y = s, col = mth, linetype = mth)) +
  scale_x_continuous(name = "Lead time (hours)", expand = c(0, 0)) +
  scale_y_continuous(name = "CRPS (m/s)", limits = c(0.5, 0.95)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.justification = c(0.5, 1),
        legend.position = c(0.5, 0.99),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA)) +
  guides(col = guide_legend(nrow = 2, byrow = TRUE))
ggsave("Figures/cs_crps_vs_lead_mbm.png", height = 5*0.5, width = 8*0.5)
range(1 - s[, 5]/s[, 1])
range(1 - s[, 6]/s[, 1])


### plot crps of order-statistics approach

s <- apply(results1_mbm$crps, c(1, 3), mean)
s <- s[, c(5, 6, 7, 8)]
colnames(s) <- c("LP  Equal", "LP Discrete", "LP Point", "LP Ordered")

df <- data.frame(lt = 1:33, s = as.vector(s), mth = rep(c("Linear pool", "Weighted", "Weighted (OS)"), each = 33))
ggplot(df) + geom_line(aes(x = lt, y = s, col = mth)) +
  geom_hline(aes(yintercept = 0), lty = "dotted") +
  scale_x_continuous(name = "Lead time (hours)", expand = c(0, 0)) +
  scale_y_continuous(name = "CRPSS", limits = c(-0.05, 0.5)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.justification = c(0.5, 1),
        legend.position = c(0.5, 0.99),
        legend.direction="horizontal")

df <- data.frame(lt = 1:33, s = as.vector(s), mth = rep(colnames(s), each = 33))
ggplot(df) + geom_line(aes(x = lt, y = s, col = mth, linetype = mth)) +
  scale_x_continuous(name = "Lead time (hours)", expand = c(0, 0)) +
  scale_y_continuous(name = "CRPS (m/s)", limits = c(0.45, 0.85)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.justification = c(0.5, 1),
        legend.position = c(0.5, 0.99),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA)) +
  guides(col = guide_legend(nrow = 2, byrow = TRUE))
ggsave("Figures/cs_wcrpss_vs_lead_mbm.png", height = 5*0.5, width = 8*0.5)


### plot es vs lead time

s <- results3_mbm$es
s <- s[, -1]
colnames(s) <- c("COSMO-1E", "COSMO-2E", "ECMWF IFS", "LP  Equal", "LP Discrete", "LP Point")

df <- data.frame(lt = 1:33, s = as.vector(s), mth = rep(colnames(s), each = 33))
ggplot(df) + geom_line(aes(x = lt, y = s, col = mth, linetype = mth)) +
  scale_x_continuous(name = "Lead time (hours)", expand = c(0, 0)) +
  scale_y_continuous(name = "Energy score", limits = c(6.5, 11.3)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.justification = c(0.5, 1),
        legend.position = c(0.5, 0.99),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA)) +
  guides(col = guide_legend(nrow = 2, byrow = TRUE))
ggsave("Figures/cs_es_vs_lead_mbm.png", height = 5*0.5, width = 8*0.5)
1 - s[, 5]/s[, 1]
1 - s[, 6]/s[, 5]



### plot weights of order statistics

lt <- 18
w <- colMeans(results1_mbm$w$ord[lt, , 1:11])
w <- colMeans(results1_mbm$w$ord[lt, , 12:32])
w <- colMeans(results1_mbm$w$ord[lt, , 33:83])
M <- length(w)
df <- data.frame(ord = 1:M, w = w)
ggplot(df) + geom_bar(aes(x = ord, y = w), stat = "identity") + 
  #scale_x_continuous(name = "Order statistic", breaks = 1:M, labels = 1:M, expand = c(0, 0)) +
  scale_x_continuous(name = "Order statistic", breaks = c(1, 5, 10, 15, 21), 
                     labels = c(1, 5, 10, 15, 21), expand = c(0, 0)) +
  #scale_x_continuous(name = "Order statistic", breaks = c(1, seq(10, 40, 10), 51), 
  #                   labels = c(1, seq(10, 40, 10), 51), expand = c(0, 0)) +
  scale_y_continuous(name = "Weight", limits = c(0, 0.1), expand = c(0, 0)) +
  ggtitle("COSMO-2E") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 12))
ggsave("Figures/cs_weight_vs_ord_c2_mbm.png", height = 5*0.4 + 0.3, width = 8*0.4)



### PIT histograms

lt <- 18
dat <- load_data(lt, stat_id = NULL) # load dat
dat <- dat %>% # restrict attention to stations with complete set of observations
  group_by(nat_abbr) %>%
  filter(n() == 1030) %>% 
  ungroup() %>% filter(nat_abbr %in% stat_list)

tr_dat <- dat %>% filter(reftime < "2022-06-01")
ts_dat <- dat %>% filter(reftime >= "2022-06-01")
rm(dat)

y_tr <- tr_dat$obs
y <- ts_dat$obs

get_mbm_dat <- function(tr_dat, ts_dat, stat_ids) {
  
  x_c1 <- x_c2 <- x_ifs <- y <- c()
  for (i in seq_along(stat_ids)) {
    print(i)
    
    tr_dat_i <- tr_dat %>% filter(nat_abbr == stat_ids[i])
    ts_dat_i <- ts_dat %>% filter(nat_abbr == stat_ids[i])
    
    y_tr <- tr_dat_i$obs
    y <- c(y, ts_dat_i$obs)
    
    # COSMO-1E 
    x_tr <- tr_dat_i %>% dplyr::select(`COSMO-1E`) %>% as.matrix()
    x <- ts_dat_i %>% dplyr::select(`COSMO-1E`) %>% as.matrix()
    x <- mbm_mom_est(y_tr, x_tr, x)
    x_c1 <- rbind(x_c1, x)
    
    # COSMO-2E 
    x_tr <- tr_dat_i %>% dplyr::select(`COSMO-2E`) %>% as.matrix()
    x <- ts_dat_i %>% dplyr::select(`COSMO-2E`) %>% as.matrix()
    x <- mbm_mom_est(y_tr, x_tr, x)
    x_c2 <- rbind(x_c2, x)
    
    # IFS
    x_tr <- tr_dat_i %>% dplyr::select(ECMWF_IFS) %>% as.matrix()
    x <- ts_dat_i %>% dplyr::select(ECMWF_IFS) %>% as.matrix()
    x <- mbm_mom_est(y_tr, x_tr, x)
    x_ifs <- rbind(x_ifs, x)
  }
  
  return(list(c1 = x_c1, c2 = x_c2, ifs = x_ifs, y = y))
}

x_mbm <- get_mbm_dat(tr_dat, ts_dat, stat_list)

y <- x_mbm$y

# COSMO-1E
x <- x_mbm$c1
pit <- rowMeans(x < y) + runif(length(y))*(rowMeans(x <= y) - rowMeans(x < y))
pit_hist(pit, ranks = FALSE, ymax = 0.6, xticks = FALSE, xlab = NULL, yticks = FALSE, ylab = NULL, title = "COSMO-1E") +
  theme(plot.title = element_text(size = 11))
ggsave("Figures/cs_pit_c1_mbm.png", height = 7*0.25, width = 8*0.25)

# COSMO-2E
x <- x_mbm$c2
pit <- rowMeans(x < y) + runif(length(y))*(rowMeans(x <= y) - rowMeans(x < y))
pit_hist(pit, ranks = FALSE, ymax = 0.6, xticks = FALSE, xlab = NULL, yticks = FALSE, ylab = NULL, title = "COSMO-2E") +
  theme(plot.title = element_text(size = 11))
ggsave("Figures/cs_pit_c2_mbm.png", height = 7*0.25, width = 8*0.25)

# IFS
x <- x_mbm$ifs
pit <- rowMeans(x < y) + runif(length(y))*(rowMeans(x <= y) - rowMeans(x < y))
pit_hist(pit, ranks = FALSE, ymax = 0.6, xticks = FALSE, xlab = NULL, yticks = FALSE, ylab = NULL, title = "ECMWF IFS") +
  theme(plot.title = element_text(size = 11))
ggsave("Figures/cs_pit_ifs_mbm.png", height = 7*0.25, width = 8*0.25)

# Multi-model
x <- cbind(x_mbm$c1, x_mbm$c2, x_mbm$ifs)
pit <- rowMeans(x < y) + runif(length(y))*(rowMeans(x <= y) - rowMeans(x < y))
pit_hist(pit, ranks = FALSE, ymax = 0.6, xticks = FALSE, xlab = NULL, yticks = FALSE, ylab = NULL, title = "LP Equal") +
  theme(plot.title = element_text(size = 11))
ggsave("Figures/cs_pit_mm_mbm.png", height = 7*0.25, width = 8*0.25)

# Linear pool
w <- results1_mbm$w$lp[lt, , ]
w <- cbind(replicate(11, w[, 1]/11), replicate(21, w[, 2]/21), replicate(51, w[, 3]/51))
w <- lapply(1:nrow(w), function(i) replicate(361, w[i, ]) %>% t()) %>% do.call(rbind, .)
pit <- rowSums(w*(x < y)) + runif(length(y))*(rowMeans(w*(x <= y)) - rowMeans(w*(x < y)))
pit_hist(pit, ranks = FALSE, ymax = 0.6, xticks = FALSE, xlab = NULL, yticks = FALSE, ylab = NULL, title = "LP Discrete") +
  theme(plot.title = element_text(size = 11))
ggsave("Figures/cs_pit_lp_mbm.png", height = 7*0.25, width = 8*0.25)

# Weighted
w <- results1_mbm$w$wtd[lt, , ]
w <- lapply(1:nrow(w), function(i) replicate(361, w[i, ]) %>% t()) %>% do.call(rbind, .)
pit <- rowSums(w*(x < y)) + runif(length(y))*(rowMeans(w*(x <= y)) - rowMeans(w*(x < y)))
pit_hist(pit, ranks = FALSE, ymax = 0.6, xticks = FALSE, xlab = NULL, yticks = FALSE, ylab = NULL, title = "LP Point") +
  theme(plot.title = element_text(size = 11))
ggsave("Figures/cs_pit_wtd_mbm.png", height = 7*0.25, width = 8*0.25)

# Weighted
x1 <- x_mbm$c1 %>% apply(1, sort) %>% t()
x2 <- x_mbm$c2 %>% apply(1, sort) %>% t()
x3 <- x_mbm$ifs %>% apply(1, sort) %>% t()
x <- cbind(x1, x2, x3)
w <- results1_mbm$w$ord[lt, , ]
w <- lapply(1:nrow(w), function(i) replicate(361, w[i, ]) %>% t()) %>% do.call(rbind, .)
pit <- rowSums(w*(x < y)) + runif(length(y))*(rowMeans(w*(x <= y)) - rowMeans(w*(x < y)))
pit_hist(pit, ranks = FALSE, ymax = 0.6, xticks = FALSE, xlab = NULL, yticks = FALSE, ylab = NULL, title = "LP Ordered") +
  theme(plot.title = element_text(size = 11))
ggsave("Figures/cs_pit_ord_mbm.png", height = 7*0.25, width = 8*0.25)



### plot average cdf before and after weighting

x <- cbind(x_mbm$c1, x_mbm$c2, x_mbm$ifs)
xx <- seq(0, 11, 0.01)  

F_mm <- sapply(xx, function(z) rowMeans(x <= z))
F_mm <- colMeans(F_mm)

w <- results1_mbm$w$lp[lt, , ]
w <- cbind(replicate(11, w[, 1]/11), replicate(21, w[, 2]/21), replicate(51, w[, 3]/51))
w <- lapply(1:nrow(w), function(i) replicate(361, w[i, ]) %>% t()) %>% do.call(rbind, .)
F_lp <- sapply(xx, function(z) rowSums(w*(x <= z)))
F_lp <- colMeans(F_lp)

x1 <- x_mbm$c1 %>% apply(1, sort) %>% t()
x2 <- x_mbm$c2 %>% apply(1, sort) %>% t()
x3 <- x_mbm$ifs %>% apply(1, sort) %>% t()
x <- cbind(x1, x2, x3)
w <- results1_mbm$w$ord[lt, , ]
w <- lapply(1:nrow(w), function(i) replicate(361, w[i, ]) %>% t()) %>% do.call(rbind, .)
F_ord <- sapply(xx, function(z) rowSums(w*(x <= z)))
F_ord <- colMeans(F_ord)


df <- data.frame(x = xx, y = c(F_mm, F_lp, F_ord), 
                 mth = rep(c("Multi-model", "Linear pool", "Weighted (OS)"), each = length(xx)))
ggplot(df) + geom_line(aes(x = x, y = y, col = mth)) +
  scale_x_continuous(name = "x (m/s)", limits = c(0, 11), breaks = seq(0, 10, 2), expand = c(0, 0)) +
  scale_y_continuous(name = "F(x)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.justification = c(1, 0),
        legend.position = c(0.99, 0.01))
ggsave("Figures/cs_average_dist_mbm.png", height = 5*0.5, width = 8*0.5)
