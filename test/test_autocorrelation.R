# install.packages(c("AER",
#                    "dynlm",
#                    "forecast",
#                    "readxl",
#                    "stargazer",
#                    "scales",
#                    "quantmod",
#                    "urca"))

library(AER)
library(dynlm)
library(forecast)
library(readxl)
library(stargazer)
library(scales)
library(quantmod)
library(urca)
library(dplyr)

load("data/multiomics_ad_serum_pos")

Feature <- "F2"
plotFeature(multiomics_ad_serum_pos,
            Feature)

data_rlm <-
  multiomics_ad_serum_pos %>%
  filter(feature == Feature) %>%
  calibrateBatch.intra()

data_rlm %>%
  mutate(feature = Feature,
         intensity = intensity_calibrated) %>%
  plotFeature.method(Feature, calibration_method = "rlm")

data_ar1 <-
  multiomics_ad_serum_pos %>%
  filter(feature == Feature) %>%
  pull(intensity) %>% ar.ols(order.max = 1)
data.frame(intensity = data_ar1$resid,
           injection_sequence = 1: length(data_ar1$resid),
           feature = Feature) %>%
  plotFeature.method(Feature = Feature, calibration_method = "Ar(1)")

data_ar1_fill1 <- data_ar1$resid
data_ar1_fill1[1] <-
  multiomics_ad_serum_pos %>%
  filter(feature == Feature,
         injection_sequence ==1) %>%
  pull(intensity)

data_ar1_rlm <-
  data.frame(intensity = data_ar1_fill1,
             injection_sequence = 1: length(data_ar1_fill1)) %>%
  calibrateBatch.intra()
data.frame(intensity = data_ar1_fill1,
           injection_sequence = 1: length(data_ar1$resid),
           feature = Feature) %>%
  plotFeature.method(Feature = Feature, calibration_method = "Ar(1)+rlm")


data_rlm_ar1 <-
data_rlm %>%
  pull(intensity_calibrated) %>%
  ar.ols(order.max = 1)
data.frame(intensity = data_rlm_ar1$resid,
           injection_sequence = 1: length(data_rlm_ar1$resid),
           feature = Feature) %>%
  plotFeature.method(Feature = Feature, calibration_method = "rlm+Ar(1)")



#### Calculate effect_matrix of raw data ####
feature_n <-
  multiomics_ad_serum_pos %>%
  pull(feature) %>%
  unique() %>%
  length()

acf_ar1_compute <- function(Feature = ...,
                            data = ...,
                            intensity_colname = ...){
acf <-
  data %>%
  filter(feature == Feature) %>% pull(intensity_colname) %>% acf
acf[["acf"]][2]}

p_value_compute <- function(Feature = ...,
                            data = ...,
                            intensity_colname = ...){
  feature_data <-
    multiomics_ad_serum_pos %>%
    filter(feature == Feature)
  lm <- lm(intensity_colname ~ injection_sequence, feature_data) %>% summary()
  p <- lm$coefficients[2,4]
}

i <- 1 : feature_n
# acf_ar1 <- sapply(i, function(i){acf_ar1_compute(Feature = paste0("F",i))})
# p_value <- sapply(i, function(i){p_value_compute(Feature = paste0("F",i))})

acf_ar1 %>% qplot
p_value %>% qplot

effect_matrix_raw <-
  data.frame(acf_ar1,
             p_value,
             Feature = paste0("F",i))

acf_ar1_feature %>% top_n(10, acf_ar1)

#### 1 - Calibrate batch effects by rlm ####

multiomics_ad_serum_pos_rlm <-
  multiomics_ad_serum_pos %>%
  mutate(intensity = ifelse(is.na(intensity), .01, intensity)) %>%
  split(.$feature) %>%
  map_dfr(function(x){
    data.frame(
      intensity = x$intensity,
      injection_sequence = x$injection_sequence
    ) %>% calibrateBatch.intra()
  }) %>%
  bind_cols(multiomics_ad_serum_pos) %>%
  rename(intensity_rlm = intensity_calibrated)


#### 1 - Calculate effect matrix of rlm data ####
acf_ar1_rlm <-
  sapply(i,
         function(i){
           acf_ar1_compute(Feature = paste0("F",i),
                           data = multiomics_ad_serum_pos_rlm,
                           intensity_colname = "intensity_rlm")
           })

p_value_rlm <-
  sapply(i,
         function(i){
           p_value_compute(Feature = paste0("F",i),
                           data = multiomics_ad_serum_pos_rlm,
                           intensity_colname = "intensity_rlm")})


effect_matrix_rlm <-
  data.frame(acf_ar1_rlm,
             p_value_rlm,
             Feature = paste0("F",i))



#### 2 - Calibrate batch effects by AR(1) ####




#### 2 - Calculate effect matrix of AR(1) data ####



































