################################################################################
####################### utility functions used in main.R #######################
################################################################################


################################################################################
##### set up functions

load_data <- function(lt = NULL, stat_id = NULL, path = "C:/Users/sa20i493/Documents/Data/MeteoSwiss/") {
  
  # the fcsts are stored as data.frame columns
  # therefore, the three forecast sources cannot be concatenated
  # into a single dataset (as the number of members differ)
  fcst <- sapply(
    paste0(path, c("COSMO-1E", "COSMO-2E", "ECMWF_IFS")),
    open_dataset,
    simplify = FALSE
  )
  obs <- open_dataset(paste0(path, "obs"))
  
  
  ## read in single station
  ## (for performance, but reading in all for one lead time might work)
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


get_stations <- function(lt_vec = 1:33, new = FALSE) {
  
  if (new) {
    stat_ids <- vector("list", length(lt_vec))
    
    for (i in lt_vec) {
      print(i)
      dat <- load_data(i, path = "C:/Users/sa20i493/Documents/Data/MeteoSwiss/") # load dat
      # restrict attention to stations with complete set of observations
      dat <- dat %>% group_by(nat_abbr) %>% filter(n() == 1030) %>% ungroup()
      stat_ids[[i]] <- unique(dat$nat_abbr)
    }
    
    stat_list <- Reduce(intersect, stat_ids)
    save(stat_list, file = "Data/stat_list.RData")
  } else {
    load("Data/stat_list.RData")
    return(stat_list)
  }
  
}


################################################################################
##### wrapper to get results functions

## univariate
get_results_uv <- function(kernel, lt_vec = 1:33, stat_ids = stat_list, mbm = FALSE) {
  w_dsc <- array(NA, c(length(lt_vec), length(stat_ids), 3))
  w_poi <- w_ord <- array(NA, c(length(lt_vec), length(stat_ids), 83))
  
  crps_mat <- array(NA, c(length(lt_vec), length(stat_ids), 7))
  dimnames(crps_mat)[[3]] <- c("C1", "C2", "IFS", "LP-Eq", "LP-Ds", "LP-Po", "LP-Or")
  
  for (i in seq_along(lt_vec)) {
    lead <- lt_vec[i]
    
    dat <- load_data(lead, stat_id = NULL) # load dat
    # restrict attention to stations with complete set of observations
    dat <- dat %>% group_by(nat_abbr) %>% filter(n() == 1030) %>% ungroup()
    
    for (j in seq_along(stat_ids)) {
      print(c(lead, j))
      
      tr_dat <- dat %>% filter(nat_abbr == stat_ids[j], reftime < "2022-06-01")
      ts_dat <- dat %>% filter(nat_abbr == stat_ids[j], reftime >= "2022-06-01")
      
      if (mbm) {
        scores <- get_mbm_scores(tr_dat, ts_dat, kernel)
      } else {
        scores <- get_scores(tr_dat, ts_dat, kernel)
      }
      
      crps_mat[i, j, ] <- scores$crps
      w_dsc[i, j, ] <- scores$w$dsc
      w_poi[i, j, ] <- scores$w$poi
      w_ord[i, j, ] <- scores$w$ord
    }
  }
  
  return(list(crps = crps_mat, w = list(lp = w_lp, wtd = w_wtd, ord = w_ord)))
  
}

## multivariate
get_results_mv <- function(kernel, lt_vec = 1:33, stat_ids = stat_list, mbm = FALSE) {
  w_lp <- array(NA, c(length(lt_vec), 3))
  w_wtd <- array(NA, c(length(lt_vec), 83))
  crps_mat <- array(NA, c(length(lt_vec), length(stat_ids), 6))
  es_mat <- array(NA, c(length(lt_vec), 6))
  dimnames(crps_mat)[[3]] <- colnames(es_mat) <- c("C1", "C2", "IFS", "LP-Eq", "LP-Ds", "LP-Po")
  
  for (i in seq_along(lt_vec)) {
    lead <- lt_vec[i]
    print(lead)
    
    dat <- load_data(lead) # load dat
    # restrict attention to stations with complete set of observations
    dat <- dat %>% filter(nat_abbr %in% stat_ids) %>% group_by(nat_abbr) %>% filter(n() == 1030) %>% ungroup() 
    
    tr_dat <- dat %>% filter(reftime < "2022-06-01")
    ts_dat <- dat %>% filter(reftime >= "2022-06-01")
    rm(dat)
    
    if (mbm) {
      scores <- get_mv_mbm_scores(tr_dat, ts_dat, stat_ids, kernel)
    } else {
      scores <- get_mv_scores(tr_dat, ts_dat, stat_ids, kernel)
    }
    
    crps_mat[i, , ] <- scores$crps
    es_mat[i, ] <- scores$es
    w_lp[i, ] <- scores$w$lp
    w_wtd[i, ] <- scores$w$wtd
  }
  
  return(list(crps = crps_mat, es = es_mat, w = list(lp = w_lp, wtd = w_wtd)))
}


################################################################################
##### kernel functions

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


################################################################################
##### post-processing functions

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


################################################################################
##### evaluate functions

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


################################################################################
##### plot functions

plot_map <- function(lons, lats, z, title = NULL, filename = NULL){

  ## elevation data
  #dem <- geodata::elevation_30s(country = "CHE", path = tempdir())
  dem <- geodata::elevation_global(res = 0.5, path = tempdir())
  
  xmin <- 5.8  # Western edge
  xmax <- 10.5 # Eastern edge
  ymin <- 45.8 # Southern edge
  ymax <- 47.9 # Northern edge
  
  region_extent <- terra::ext(xmin, xmax, ymin, ymax)
  cropped_dem <- terra::crop(dem, region_extent)
  dem_df <- as.data.frame(cropped_dem, xy = TRUE)
  
  colnames(dem_df) <- c("lon", "lat", "elevation")
  
  ## boundary data
  switzerland <- rnaturalearth::ne_countries(scale = "large", country = "Switzerland", returnclass = "sf")
  
  maxplot <- is.matrix(z)
  
  if (maxplot) {

    # plot data
    z <- apply(z, 1, which.max)
    z[z == 1] <- c("COSMO-1E")
    z[z == 2] <- c("COSMO-2E")
    z[z == 3] <- c("ECMWF IFS")
    df <- data.frame(lat = lats, lon = lons, Model = z)
    
    # plot
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
    
  } else {
    
    # plot data
    df <- data.frame(lat = lats, lon = lons, z = z)
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
  }
  

  if (!is.null(filename)) {
    if (maxplot) {
      ggsave(plot = plot_obj, filename, height = 3.6, width = 6)
    } else {
      ggsave(plot = plot_obj, filename, height = 3, width = 5)
    }
  } else {
    return(plot_obj)
  }
  
}

# plot weight vs lead time
plot_w_vs_lt <- function(w, filename = NULL) {
  
  if (ncol(w) == 83) w <- cbind(rowSums(w[, 1:11]), rowSums(w[, 12:32]), rowSums(w[, 33:83]))
  
  df <- data.frame(lt = 1:33, 
                   w = as.vector(w),
                   mth = rep(c("COSMO-1E", "COSMO-2E", "ECMWF IFS"), each = 33))
  plot_obj <- ggplot(df) + 
    geom_line(aes(x = lt, y = w, col = mth, linetype = mth)) +
    scale_x_continuous(name = "Lead time (hours)", expand = c(0, 0)) +
    scale_y_continuous(name = "Weight", limits = c(0, 1), expand = c(0, 0)) +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.justification = c(0.5, 1),
          legend.position = c(0.5, 0.99),
          legend.background = element_rect(fill = "transparent", color = NA),
          legend.key = element_rect(fill = "transparent", color = NA),
          legend.direction="horizontal")
  
  if (!is.null(filename)) {
    ggsave(plot = plot_obj, filename, height = 2.5, width = 4)
  } else {
    return(plot_obj)
  }
}

# plot weight vs mean squared error
plot_w_vs_mse <- function(lt, uv = TRUE, ylims = c(0, 0.1), filename = NULL) {
  
  dat <- load_data(lt) %>% group_by(nat_abbr) %>% filter(n() == 1030) %>% 
    ungroup() %>% filter(reftime >= "2022-06-01", nat_abbr %in% stat_list)
  
  if (uv) {
    mse <- c(colMeans((dat$`COSMO-1E` - dat$obs)^2),
             colMeans((dat$`COSMO-2E` - dat$obs)^2),
             colMeans((dat$ECMWF_IFS - dat$obs)^2))
    
    w <- colMeans(results_uv$w$wtd[lt, , ])
    
    height <- 2.5
    width <- 4
  } else {
    y <- dat %>% select(nat_abbr, obs) %>% 
      group_by(nat_abbr) %>%
      mutate(row = row_number()) %>%
      pivot_wider(names_from = nat_abbr, values_from = obs) %>%
      select(-row) %>% as.matrix()
    
    x <- lapply(stat_list, function(z) {
      dat %>% filter(nat_abbr == z) %>% 
        select(`COSMO-1E`, `COSMO-2E`, ECMWF_IFS) %>%
        cbind() %>% as.matrix()
    }) %>% simplify2array() %>% aperm(c(1, 3, 2))
    
    mse <- colMeans(sqrt(apply((x - replicate(83, y))^2, c(1, 3), sum)))
    
    w <- results_mv$w$wtd[lt, ]
    
    height <- 2.55
    width <- 4.08
  }
  
  df <- data.frame(s = mse, w = w, 
                   mth = c(rep("COSMO-1E", 11), rep("COSMO-2E", 21), rep("ECMWF IFS", 51)),
                   control = c(1, rep(0, 10), 1, rep(0, 20), 1, rep(0, 50)))
  plot_obj <- ggplot(df) + geom_point(aes(x = s, y = w, col = mth, shape = as.factor(control))) +
    scale_x_continuous(name = "MSE") +
    scale_y_continuous(name = "Weight", limits = ylims) +
    scale_shape_manual(values = c(19, 4)) +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.justification = c(0.5, 1),
          legend.position = c(0.5, 0.99),
          legend.background = element_rect(fill = "transparent", color = NA),
          legend.key = element_rect(fill = "transparent", color = NA),
          legend.direction="horizontal") +
    guides(shape = "none")
  
  if (!is.null(filename)) {
    ggsave(plot = plot_obj, filename, height = height, width = width)
  } else {
    return(plot_obj)
  }
}

# plot score vs lead time
plot_score_vs_lt <- function(s, uv = TRUE, ylims = c(0.5, 2), filename = NULL) {
  
  if (ncol(s) == 6) {
    colnames(s) <- c("COSMO-1E", "COSMO-2E", "ECMWF IFS", "LP  Equal", "LP Discrete", "LP Point")
  } else if (ncol(s) == 4) {
    colnames(s) <- c("LP  Equal", "LP Discrete", "LP Point", "LP Ordered")
  }
  
  if (uv) {
    ylab <- "CRPS (m/s)"
  } else {
    ylab <- "Energy score"
  }
  
  df <- data.frame(lt = 1:33, s = as.vector(s), mth = rep(colnames(s), each = 33))
  plot_obj <- ggplot(df) + geom_line(aes(x = lt, y = s, col = mth, linetype = mth)) +
    scale_x_continuous(name = "Lead time (hours)", expand = c(0, 0)) +
    scale_y_continuous(name = "CRPS (m/s)", limits = c(0.5, 2)) +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.justification = c(0.5, 1),
          legend.position = c(0.5, 0.99),
          legend.background = element_rect(fill = "transparent", color = NA),
          legend.key = element_rect(fill = "transparent", color = NA)) +
    guides(col = guide_legend(nrow = 2, byrow = TRUE))
  
  if (!is.null(filename)) {
    ggsave(plot = plot_obj, filename, height = 2.5, width = 4)
  } else {
    return(plot_obj)
  }
}

# plot weights of order statistics
plot_order_w <- function(lt, w_arr, method = c("c1", "c2", "ifs"), filename = NULL) {
  
  method <- match.arg(method)
  
  if (method == "c1") {
    title <- "COSMO-1E"
    w <- colMeans(w_arr[lt, , 1:11])
    sc_x <- scale_x_continuous(name = "Order statistic", breaks = 1:M, labels = 1:M, expand = c(0, 0))
  } else if (method == "c2") {
    title <- "COSMO-2E"
    sc_x <- scale_x_continuous(name = "Order statistic", breaks = c(1, 5, 10, 15, 21),
                               labels = c(1, 5, 10, 15, 21), expand = c(0, 0))
    w <- colMeans(w_arr[lt, , 12:32])
  } else if (method == "ifs") {
    title <- "ECMWF IFS"
    w <- colMeans(w_arr[lt, , 33:83])
    sc_x <- scale_x_continuous(name = "Order statistic", breaks = c(1, seq(10, 40, 10), 51), 
                               labels = c(1, seq(10, 40, 10), 51), expand = c(0, 0))
  }
  
  M <- length(w)
  df <- data.frame(ord = 1:M, w = w)
  plot_obj <- ggplot(df) + geom_bar(aes(x = ord, y = w), stat = "identity") + sc_x +
    scale_y_continuous(name = "Weight", limits = c(0, 0.1), expand = c(0, 0)) +
    ggtitle("ECMWF IFS") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(size = 12))
  
  if (!is.null(filename)) {
    ggsave(plot = plot_obj, filename, height = 2.3, width = 3.2)
  } else {
    return(plot_obj)
  }
  
}

# plot PIT histograms
plot_pit(lt, method = c("c1", "c2", "ifs", "lp-eq", "lp-ds", "lp-po", "lp-or"), w_arr = NULL, filename = NULL) {
  
  method <- match.arg(method)
  
  # load test data
  dat <- load_data(lt) 
  dat <- dat %>% group_by(nat_abbr) %>% filter(n() == 1030) %>% ungroup() %>% 
    filter(reftime >= "2022-06-01", nat_abbr %in% stat_list)
  y <- dat$obs
  
  if (method == "c1") {
    x <- as.matrix(dat$`COSMO-1E`)
    title <- "COSMO-1E"
  } else if (method == "c2") {
    x <- as.matrix(dat$`COSMO-2E`)
    title <- "COSMO-2E"
  } else if (method == "ifs") {
    x <- as.matrix(dat$ECMWF_IFS)
    title <- "ECMWF IFS"
  } else if (method == "lp-eq") {
    x <- cbind(as.matrix(dat$`COSMO-1E`), as.matrix(dat$`COSMO-2E`), as.matrix(dat$ECMWF_IFS))
    title <- "LP Equal"
  } else if (method == "lp-ds") {
    x <- cbind(as.matrix(dat$`COSMO-1E`), as.matrix(dat$`COSMO-2E`), as.matrix(dat$ECMWF_IFS))
    w <- w_arr[lt, , ]
    w <- cbind(replicate(11, w[, 1]/11), replicate(21, w[, 2]/21), replicate(51, w[, 3]/51))
    w <- lapply(1:nrow(w), function(i) replicate(361, w[i, ]) %>% t()) %>% do.call(rbind, .)
    title <- "LP Discrete"
  } else if (method == "lp-po") {
    x <- cbind(as.matrix(dat$`COSMO-1E`), as.matrix(dat$`COSMO-2E`), as.matrix(dat$ECMWF_IFS))
    w <- w_arr[lt, , ]
    w <- lapply(1:nrow(w), function(i) replicate(361, w[i, ]) %>% t()) %>% do.call(rbind, .)
    title <- "LP Point"
  } else if (method == "lp-or") {
    x1 <- as.matrix(dat$`COSMO-1E`)%>% apply(1, sort) %>% t()
    x2 <- as.matrix(dat$`COSMO-2E`)%>% apply(1, sort) %>% t()
    x3 <- as.matrix(dat$ECMWF_IFS)%>% apply(1, sort) %>% t()
    x <- cbind(x1, x2, x3)
    w <- w_arr[lt, , ]
    w <- lapply(1:nrow(w), function(i) replicate(361, w[i, ]) %>% t()) %>% do.call(rbind, .)
    title <- "LP Ordered"
  }
  
  if (method %in% c("c1", "c2", "ifs", "lp-eq")) w <- array(1, dim(x))
  
  pit <- rowSums(w*(x < y)) + runif(length(y))*(rowMeans(w*(x <= y)) - rowMeans(w*(x < y)))
  
  plot_obj <- pit_hist(pit, ranks = FALSE, ymax = 0.6, xticks = FALSE, 
                       xlab = NULL, yticks = FALSE, ylab = NULL, title = title) +
    theme(plot.title = element_text(size = 11))
  
  if (!is.null(filename)) {
    ggsave(plot = plot_obj, filename, height = 7*0.25, width = 8*0.25)
  } else {
    return(plot_obj)
  }
  
}


