#' Deprecated calibration method (without the support of custom column name)
calibrateBatch.intra.rlm.old <- function(data = ...,
                                     intensity = intensity,
                                     injection_sequence = injection_sequence) {
    intensity <- data %>% dplyr::pull(intensity)
    injection_sequence <- data %>% dplyr::pull(injection_sequence)
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

#' Intra batch calibration by robust linear modelling (rlm)
#'
#' @param data Metabolomics data in long-format
#' @param intensity The column name of intensity (by default intensity)
#' @param injection_sequence The column name of injection sequence (by default injection_sequence)
#' @example
#' data.frame(intensity = runif(50, min=10, max=50), injection_sequence = 1:50) %>% libra::calibrateBatch.intra.rlm()
calibrateBatch.intra.rlm <- function(data=...,
                                     intensity = intensity,
                                     injection_sequence = injection_sequence){
    intensity <- rlang::enexpr(intensity)
    injection_sequence <- rlang::enexpr(injection_sequence)
    intensity. <- data %>% dplyr::pull(!!intensity)
    injection_sequence. <- data %>% dplyr::pull(!!injection_sequence)

    rlm <- MASS::rlm(intensity. ~ injection_sequence.)
    rlm_summary <- rlm %>% summary
    slope <- rlm_summary[["coefficients"]][2, 1]
    intercept <- rlm_summary[["coefficients"]][1, 1]

    # calibration center
    center_injec_seq <- length(injection_sequence)/2
    center_intensity <- center_injec_seq * slope + intercept

    # calibrted data
    intensity_calibration <- rlm_summary[["residuals"]] + center_intensity
    calibrated_data <- data %>% mutate(intensity_calibrated = intensity_calibration, center_intensity)
    return(calibrated_data)
}

calibrateBatch.intra.ar1 <- function(data = ...) {
    intensity <- data %>% dplyr::pull(intensity)
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
