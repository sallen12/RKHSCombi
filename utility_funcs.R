
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


