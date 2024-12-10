
load_data <- function(lt = NULL, stat_id = NULL, path = "C:/Users/sa20i493/Documents/Data/MeteoSwiss/") {

  library(arrow)
  library(dplyr)
  library(tidyr)
  library(ggplot2)

  # unfortunately, I stored the fcsts as data.frame columns
  # therefore you cannot concatenate the three forecast sources
  # into a single dataset (as the number of members differ)
  fcst <- sapply(
    paste0(path, c("COSMO-1E", "COSMO-2E", "ECMWF_IFS")),
    open_dataset,
    simplify = FALSE
  )
  obs <- open_dataset(paste0(path, "obs"))


  ## read in single station
  ## (for performance, but reading in all for one lead time might work)
  ## you may want to consider repartitioning the dataset by lead time -
  ## now partitioned by reference time - to optimize for analysis access
  ## patterns
  if (!is.null(stat_id) & !is.null(lt)) {
    ff <- lapply(
      names(fcst),
      function(nn) {
        print(paste0("load ", nn))
        fcst[[nn]] %>%
          dplyr::filter(lead == lt, nat_abbr == stat_id) %>%
          dplyr::collect() %>%
          dplyr::rename(!!nn := fcst) # change name of forecast column so that datasets can be merged
      }
    )
  } else if (is.null(stat_id) & !is.null(lt)) {
    ff <- lapply(
      names(fcst),
      function(nn) {
        print(paste0("load ", nn))
        fcst[[nn]] %>%
          dplyr::filter(lead == lt) %>%
          dplyr::collect() %>%
          dplyr::rename(!!nn := fcst) # change name of forecast column so that datasets can be merged
      }
    )
  } else if (!is.null(stat_id) & is.null(lt)) {
    ff <- lapply(
      names(fcst),
      function(nn) {
        print(paste0("load ", nn))
        fcst[[nn]] %>%
          dplyr::filter(nat_abbr == stat_id) %>%
          dplyr::collect() %>%
          dplyr::rename(!!nn := fcst) # change name of forecast column so that datasets can be merged
      }
    )
  } else {
    ff <- lapply(
      names(fcst),
      function(nn) {
        print(paste0("load ", nn))
        fcst[[nn]] %>%
          dplyr::collect() %>%
          dplyr::rename(!!nn := fcst) # change name of forecast column so that datasets can be merged
      }
    )
  }


  ## join the datasets, remove source column to allow for join
  ff <- lapply(ff, function(x) dplyr::select(x, -source)) %>%
    Reduce(dplyr::inner_join, .)

  time_of_day <- unique(lubridate::hour(ff$time))

  if (!is.null(stat_id)) {
    oo <- obs %>%
      dplyr::collect() %>%
      dplyr::filter(
        nat_abbr == stat_id,
        lubridate::hour(time) == time_of_day
      )
  } else {
    oo <- obs %>%
      dplyr::collect() %>%
      dplyr::filter(
        lubridate::hour(time) == time_of_day
      )
      
  }

  data <- ff %>% dplyr::inner_join(oo, by = c("time", "nat_abbr")) %>%
    dplyr::rename("COSMO-1E" := paste0(path, "COSMO-1E"),
                  "COSMO-2E" := paste0(path, "COSMO-2E"),
                  "ECMWF_IFS" := paste0(path, "ECMWF_IFS"))

  return(data)
}


plot_crps <- function(dat) {
  dat %>%
   dplyr::mutate(
     across(
       c("COSMO-1E", "COSMO-2E", "ECMWF_IFS"),
       function(x) {
         SpecsVerification::EnsCrps(as.matrix(x), obs)
       }
     )
   ) %>%
   dplyr::summarize(
     across(c("COSMO-1E", "COSMO-2E", "ECMWF_IFS"), mean),
     .by = "lead"
   ) %>%
   tidyr::pivot_longer(c("COSMO-1E", "COSMO-2E", "ECMWF_IFS"), names_to = "source") %>%
   ggplot(aes(x = source, y = value, fill = source)) +
   geom_bar(stat = "identity") +
   labs(y = "CRPS (m/s)") +
   theme(legend.position = "none")
}



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



plot_map <- function(lons, lats, z, title = NULL){
  if (is.matrix(z)) z <- as.vector(z)
  
  library(terra)
  library(viridis)
  library(rnaturalearth)
  library(sf)  # For handling vector data
  
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
  
  
  return(plot_obj)
  
}
