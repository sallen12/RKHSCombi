################################################################################
##### set up

library(dplyr)
library(tidyr)
library(kernlab)
library(scoringRules)
library(Rcpp)
library(ggplot2)
# devtools::install_github("sallen12/WeightedForecastVerification")
library(WeightedForecastVerification)

source("utility_funcs.R")
sourceCpp("kernels.cpp")


# lead times
lt_vec <- 1:33

# stations
stat_list <- readRDS("Data/stat_list.RDS")
stat_info <- readRDS("Data/stat_info.RDS")

# apply member-by-member post-processing?
mbm <- FALSE


################################################################################
##### estimate weights

## univariate

results_uv <- get_results_uv(kernel = "Energy", mbm = mbm)


## multivariate
results_mv <- get_results_mv(kernel = "Energy", mbm = mbm)


## save data
filename <- if (mbm) {"Results/results_data_mbm.RData"} else {"Results/results_data.RData"}
save(results_uv, results_mv, file = filename)


################################################################################
##### plot results

## load data
load(filename)


### plot stations

## station height
plot_map(stat_info$longitude, stat_info$latitude, stat_info$station_height,
         filename = "Figures/cs_alt_map.png")

### model with highest weight
lt <- 18
w <- cbind(rowSums(results_uv$w$poi[lt, , 1:11]), 
           rowSums(results_uv$w$poi[lt, , 22:32]), 
           rowSums(results_uv$w$poi[lt, , 33:83]))
plot_map(stat_info$longitude, stat_info$latitude, w, filename = "Figures/cs_weight_vs_stat.png")


### plot weights vs lead time

## univariate
plot_w_vs_lt(w = apply(results_uv$w$poi, c(1, 3), mean), lt_vec = lt_vec, filename = "Figures/cs_weight_vs_lead.png")

## multivariate
plot_w_vs_lt(w = results_mv$w$poi, lt_vec = lt_vec, filename = "Figures/cs_mvweight_vs_lead.png")


### plot weights vs mse

lt <- 18

## univariate
plot_w_vs_mse(lt, w_arr = results_uv$w$poi, ylims = c(0, 0.07), filename = "Figures/cs_weight_vs_mse.png")

## multivariate
plot_w_vs_mse(lt, w_arr = results_mv$w$poi, ylims = c(0, 0.1), filename = "Figures/cs_mvweight_vs_mse.png")


### plot scores vs lead time

## univariate (crps)
plot_sc_vs_lt(s = apply(results_uv$crps[, , -7], c(1, 3), mean), lt_vec = lt_vec,
              filename = "Figures/cs_crps_vs_lead.png")

## multivariate (energy score)
plot_sc_vs_lt(s = results_mv$es, lt_vec = lt_vec, uv = F, ylims = c(8, 22),
             filename = "Figures/cs_es_vs_lead.png")

## univariate linear pools
plot_sc_vs_lt(s = apply(results_uv$crps[, , c(4, 5, 6, 7)], c(1, 3), mean), lt_vec = lt_vec, ylims = c(0.5, 1.3),
             filename = "Figures/cs_wcrpss_vs_lead.png")


### plot weights of order statistics

lt <- 18
path <- "Figures/cs_weight_vs_ord_"

plot_order_w(lt, results_uv$w$ord, method = "c1", filename = paste0(path, "c1.png"))
plot_order_w(lt, results_uv$w$ord, method = "c2", filename = paste0(path, "c2.png"))
plot_order_w(lt, results_uv$w$ord, method = "ifs", filename = paste0(path, "ifs.png"))


### PIT histograms

lt <- 18

plot_pit(lt, method = "c1", filename = "Figures/cs_pit_c1.png")
plot_pit(lt, method = "c2", filename = "Figures/cs_pit_c2.png")
plot_pit(lt, method = "ifs", filename = "Figures/cs_pit_ifs.png")
plot_pit(lt, method = "lp-eq", filename = "Figures/cs_pit_mm.png")
plot_pit(lt, method = "lp-ds", w_arr = results_uv$w$dsc, filename = "Figures/cs_pit_lp.png")
plot_pit(lt, method = "lp-po", w_arr = results_uv$w$poi, filename = "Figures/cs_pit_wtd.png")
plot_pit(lt, method = "lp-or", w_arr = results_uv$w$ord, filename = "Figures/cs_pit_ord.png")

