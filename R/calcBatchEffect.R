calcBatchEffect.p <- function(data = ...) {
    require(tidyverse)
    feature_intensity <- data %>% pull(intensity)
    injection_sequence <- data %>% pull(injection_sequence)
    lm <- stats::lm(feature_intensity ~ injection_sequence) %>% summary()
    corr_coefficient <- lm$coefficients[2, 4]  #p-value
    return(corr_coefficient)
}

calcBatchEffect.vcf <- function(data = ...) {
    
}


calcBatchEffect.CVQC <- function(data = ..., QC_index = ...) {
    data_n <- data %>% dplyr::group_by(feature) %>% tidyr::nest()
    
    calcBatchEffect.CVQC_one <- function(data = ...) {
        intensity_QC_index <- data[QC_index, ] %>% pull(intensity)
        CVQC <- sd(intensity_QC_index)/mean(intensity_QC_index)
        return(CVQC)
    }
    data_n_CVQC <- data_n %>% dplyr::mutate(CVQC = map_dbl(data, calcBatchEffect.CVQC_one)) %>% tidyr::unnest(cols = c(data))
    return(data_n_CVQC)
}
