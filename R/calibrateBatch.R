calibrateBatch.intra.rlm <- function(data = ...) {
    require(tidyverse)
    intensity <- data %>% pull(intensity)
    injection_sequence <- data %>% pull(injection_sequence)
    
    rlm <- MASS::rlm(intensity ~ injection_sequence)
    rlm_summary <- rlm %>% summary
    
    slope <- rlm_summary[["coefficients"]][2, 1]
    intercept <- rlm_summary[["coefficients"]][1, 1]
    
    # calibration center
    center_injec_seq <- length(injection_sequence)/2
    center_intensity <- center_injec_seq * slope + intercept
    
    # calibrted data
    intensity_calibration <- rlm_summary[["residuals"]] + center_intensity
    
    calibrated_data <- data.frame(intensity_calibrated = intensity_calibration, injection_sequence)
    return(calibrated_data)
}

calibrateBatch.intra.ar1 <- function(data = ...) {
    require(tidyverse)
    intensity <- data %>% pull(intensity)
    injection_sequence <- data %>% pull(injection_sequence)
    rlm <- MASS::rlm(intensity ~ injection_sequence)
    rlm_summary <- rlm %>% summary
    slope <- rlm_summary[["coefficients"]][2, 1]
    intercept <- rlm_summary[["coefficients"]][1, 1]
    # calibration center
    center_injec_seq <- length(injection_sequence)/2
    center_intensity <- center_injec_seq * slope + intercept
    # calibrted data
    intensity_calibration <- rlm_summary[["residuals"]] + center_intensity
    calibrated_data <- data.frame(intensity_calibrated = intensity_calibration, injection_sequence)
    return(calibrated_data)
}
