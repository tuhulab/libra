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
