calcBatchEffect <- function(data = ...,
                            Feature = ...){
  require(tidyverse)
  feature_intensity <- data %>% filter(feature == Feature) %>%
    pull(intensity)
  injection_sequence <- data %>% filter(feature == Feature) %>%
    pull(injection_sequence)
  lm <- stats::lm(feature_intensity ~ injection_sequence) %>% summary()
  corr_coefficient <- lm$coefficients[2,4] #p-value
  return(corr_coefficient)
}
