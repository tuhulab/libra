calibrateBatch.intra <- function(data = ...,
                           Feature = ...){
  require(tidyverse)

  intensity <- data %>%
    filter(feature == Feature) %>%
    pull(intensity)
  injection_sequence <- data %>%
    filter(feature == Feature) %>%
    pull(injection_sequence)

  rlm <- MASS::rlm(intensity ~ injection_sequence)
  rlm_summary <- rlm %>% summary

  slope <- rlm_summary[["coefficients"]][2,1]
  intercept <- rlm_summary[["coefficients"]][1,1]

  #calibration center
  center_injec_seq <- length(injection_sequence)/2
  center_intensity <- center_injec_seq * slope + intercept

  #calibrted data
  intensity_calibration <- rlm_summary[["residuals"]] + center_intensity

  calibrated_data <-
    data.frame(intensity = intensity_calibration,
               injection_sequence,
               feature = Feature)
  return(calibrated_data)
  }
