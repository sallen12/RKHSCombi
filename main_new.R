################################################################################
##### set up

library(kernlab)
library(scoringRules)
library(Rcpp)
library(ggplot2)
source("utility_funcs.R")
sourceCpp("kernels.cpp")


# lead times
lt_vec <- 1:33

# stations
stat_list <- get_stations(lt_vec)

# apply post-processing?
mbm <- FALSE


################################################################################
##### univariate weight estimation

experiment_uv <- function(kernel, lt_vec = 1:33, stat_ids = stat_list, mbm = FALSE) {
  w_lp <- array(NA, c(length(lt_vec), length(stat_ids), 3))
  w_wtd <- w_ord <- array(NA, c(length(lt_vec), length(stat_ids), 83))

  crps_mat <- array(NA, c(length(lt_vec), length(stat_ids), 8))
  dimnames(crps_mat)[[3]] <- c("Clim.", "C1", "C2", "IFS", "MM", "LP", "Wtd", "Wtd-Ord")

  for (lt in seq_along(lt_vec)) {
    lead <- lt_vec[lt]

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
        scores <- get_mbm_scores(tr_dat, ts_dat, kernel)
      } else {
        scores <- get_scores(tr_dat, ts_dat, kernel)
      }
      
      crps_mat[lt, i, ] <- scores$crps
      w_lp[lt, i, ] <- scores$w$lp
      w_wtd[lt, i, ] <- scores$w$wtd
      w_ord[lt, i, ] <- scores$w$ord
    }
  }

  return(list(crps = crps_mat, w = list(lp = w_lp, wtd = w_wtd, ord = w_ord)))

}

results_uv <- experiment_uv(kernel = "Energy", mbm = mbm)


################################################################################
##### multivariate weight estimation

experiment_mv <- function(kernel, lt_vec = 1:33, stat_ids = stat_list, mbm = FALSE) {
  w_lp <- array(NA, c(length(lt_vec), 3))
  w_wtd <- array(NA, c(length(lt_vec), 83))
  crps_mat <- array(NA, c(length(lt_vec), length(stat_ids), 7))
  es_mat <- array(NA, c(length(lt_vec), 7))
  dimnames(crps_mat)[[3]] <- colnames(es_mat) <-
    c("Clim.", "C1", "C2", "IFS", "MM", "LP", "Wtd")

  for (lt in seq_along(lt_vec)) {
    lead <- lt_vec[lt]
    print(lt)
    
    dat <- load_data(lead, stat_id = NULL) # load dat
    dat <- dat %>% # restrict attention to stations with complete set of observations
      group_by(nat_abbr) %>%
      filter(n() == 1030) %>%
      ungroup()
    
    tr_dat <- dat %>% filter(reftime < "2022-06-01", nat_abbr %in% stat_ids)
    ts_dat <- dat %>% filter(reftime >= "2022-06-01", nat_abbr %in% stat_ids)
    rm(dat)
    
    if (mbm) {
      scores <- get_mv_mbm_scores(tr_dat, ts_dat, stat_ids, kernel)
    } else {
      scores <- get_mv_scores(tr_dat, ts_dat, stat_ids, kernel)
    }
    
    crps_mat[lt, , ] <- scores$crps
    es_mat[lt, ] <- scores$es
    w_lp[lt, ] <- scores$w$lp
    w_wtd[lt, ] <- scores$w$wtd
  }

  return(list(crps = crps_mat, es = es_mat, w = list(lp = w_lp, wtd = w_wtd)))
}

results_mv <- experiment_mv(kernel = "Energy", mbm = mbm)


################################################################################
##### evaluate

filename <- if (mbm) {"Data/results_data_mbm.RData"} else {"Data/results_data.RData"}
save(results_uv, results_mv, file = filename)


### plot stations

SMN_stations <- readRDS("C:/Users/sa20i493/Documents/R/WeightedLoss/WindGusts/SMN_stations.rds")
SMN_stations <- SMN_stations %>% filter(nat_abbr %in% stat_list)

plot_map(SMN_stations$longitude, SMN_stations$latitude, SMN_stations$station_height)
ggsave("Figures/cs_alt_map.png", height = 3, width = 5)


### plot weights vs lead time

w <- apply(results_uv$w$wtd, c(1, 3), mean)
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

w <- rowSums(results_uv$w$wtd[lt, , 1:11])
w <- w[stat_list %in% SMN_stations$nat_abbr]
plot_map(SMN_stations$longitude, SMN_stations$latitude, w, ymax = 1, title = "COSMO-1E")
ggsave("Figures/cs_weight_vs_stat_c1.png", height = 1.7, width = 3.4)

w <- rowSums(results_uv$w$wtd[lt, , 22:32])
w <- w[stat_list %in% SMN_stations$nat_abbr]
plot_map(SMN_stations$longitude, SMN_stations$latitude, w, ymax = 1, title = "COSMO-2E")
ggsave("Figures/cs_weight_vs_stat_c2.png", height = 1.7, width = 3.4)

w <- rowSums(results_uv$w$wtd[lt, , 33:83])
w <- w[stat_list %in% SMN_stations$nat_abbr]
plot_map(SMN_stations$longitude, SMN_stations$latitude, w, ymax = 1, title = "ECMWF IFS")
ggsave("Figures/cs_weight_vs_stat_ifs.png", height = 1.7, width = 3.4)


# plot station with highest weight at each location

w <- cbind(rowSums(results_uv$w$wtd[lt, , 1:11]), 
           rowSums(results_uv$w$wtd[lt, , 22:32]), 
           rowSums(results_uv$w$wtd[lt, , 33:83]))
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

w <- colMeans(results_uv$w$wtd[lt, , ])
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

w <- results_mv$w$wtd
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


w <- results_mv$w$wtd[lt, ]
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

s <- apply(results_uv$crps, c(1, 3), mean)
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

s <- apply(results_uv$crps, c(1, 3), mean)
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

s <- results_mv$es
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
#w <- colMeans(results_uv$w$ord[lt, , 1:11])
#w <- colMeans(results_uv$w$ord[lt, , 12:32])
w <- colMeans(results_uv$w$ord[lt, , 33:83])
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

w <- results_uv$w$lp[lt, , ]
w <- cbind(replicate(11, w[, 1]/11), replicate(21, w[, 2]/21), replicate(51, w[, 3]/51))
w <- lapply(1:nrow(w), function(i) replicate(361, w[i, ]) %>% t()) %>% do.call(rbind, .)
F_lp <- sapply(xx, function(z) rowSums(w*(x <= z)))
F_lp <- colMeans(F_lp)

x1 <- as.matrix(dat$`COSMO-1E`)%>% apply(1, sort) %>% t()
x2 <- as.matrix(dat$`COSMO-2E`)%>% apply(1, sort) %>% t()
x3 <- as.matrix(dat$ECMWF_IFS)%>% apply(1, sort) %>% t()
x <- cbind(x1, x2, x3)
w <- results_uv$w$ord[lt, , ]
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
w <- results_uv$w$lp[lt, , ]
w <- cbind(replicate(11, w[, 1]/11), replicate(21, w[, 2]/21), replicate(51, w[, 3]/51))
w <- lapply(1:nrow(w), function(i) replicate(361, w[i, ]) %>% t()) %>% do.call(rbind, .)
pit <- rowSums(w*(x < y)) + runif(length(y))*(rowMeans(w*(x <= y)) - rowMeans(w*(x < y)))
pit_hist(pit, ranks = FALSE, ymax = 0.6, xticks = FALSE, xlab = NULL, yticks = FALSE, ylab = NULL, title = "LP Discrete") +
  theme(plot.title = element_text(size = 11))
ggsave("Figures/cs_pit_lp.png", height = 7*0.25, width = 8*0.25)

# Weighted
w <- results_uv$w$wtd[lt, , ]
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
w <- results_uv$w$ord[lt, , ]
w <- lapply(1:nrow(w), function(i) replicate(361, w[i, ]) %>% t()) %>% do.call(rbind, .)
pit <- rowSums(w*(x < y)) + runif(length(y))*(rowMeans(w*(x <= y)) - rowMeans(w*(x < y)))
pit_hist(pit, ranks = FALSE, ymax = 0.6, xticks = FALSE, xlab = NULL, yticks = FALSE, ylab = NULL, title = "LP Ordered") +
  theme(plot.title = element_text(size = 11))
ggsave("Figures/cs_pit_ord.png", height = 7*0.25, width = 8*0.25)



