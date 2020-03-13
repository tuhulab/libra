library(libra)
library(tidyverse)
####generate ID####
param_set <- "param1"
plate <- "1"
project_name <- "Life" #upper-case
analysis_name <- "Serum" #lower-case
mode <- "neg" #"pos"or"neg"
id <- paste(project_name,analysis_name,mode,param_set,sep = "_")
id_1 <- paste(project_name,analysis_name,param_set,sep = "_")

path_peaklist <-
  file.path("xcms_result/",
            paste0(id, "peaklist.csv"))

####read XCMS result####
data_table <- readr::read_csv("/Users/tuhu/Projects/multiomics-ad-metabolome-life/xcms_result/Life_Serum_pos_param1_peaklist.csv") %>%
  mutate(feature = paste0("F",rownames(.)))
sample_table <- readr::read_csv("/Users/tuhu/Projects/multiomics-ad-metabolome-life/xcms_result/Life_Serum_pos_param1_samplelist.csv")
data <-
  data_table %>%
  pivot_longer(colnames(.) %>%
                 str_match("X\\d{1,}") %>%
                 discard(is.na),
               names_to = "code",
               values_to = "intensity") %>%
  left_join(sample_table) %>%
  filter(sample_name != "Blank",
         sample_name != "Blank_IS",
         sample_name != "Assay_blank",
         sample_name != "MetStd") %>%
  group_by(feature) %>%
  mutate(injection_sequence = 1: sum(.$feature == "F1")) %>%
  select(feature, intensity, sample_name, injection_sequence)
data_calibrate_rlm <-
  data %>%
  mutate(intensity = ifelse(is.na(intensity), .01, intensity)) %>%
  split(.$feature) %>%
  map_dfr(function(x){
    data.frame(
      intensity = x$intensity,
      injection_sequence = x$injection_sequence
    ) %>% calibrateBatch.intra()
  }) %>% bind_cols(data) %>%
  select(-injection_sequence1,
         -intensity) %>%
  rename(intensity = intensity_calibrated)
data_calibrate_rlm_wide <-
  data_calibrate_rlm %>%
  pivot_wider(names_from = feature,
              values_from = intensity)
data_calibrate_rlm <- data_calibrate_wide

#####Raw data#####

#######calculate raw data p value########
p_raw_value <-
  data %>%
  mutate(intensity = ifelse(is.na(intensity), .01, intensity)) %>%
  split(.$feature) %>%
  map_dfr(function(x){
    lm <- lm(intensity~injection_sequence,
             data.frame(intensity = x$intensity,injection_sequence = x$injection_sequence)) %>%
    summary()
    lm$coefficients[2,4]
  }) %>% t()
p_raw <-
  data.frame(feature = rownames(p_raw_value),
                    p_raw_value) %>%
  remove_rownames()

######## calculate raw data vcf value #######
acf_raw_value <-
  data %>%
  mutate(intensity = ifelse(is.na(intensity), .01, intensity)) %>%
  split(.$feature) %>%
  map_dfr(function(x){
    acf <- acf(x$intensity)
    acf[["acf"]][2]
      }) %>% t()

acf_raw <-
  data.frame(feature = rownames(acf_raw_value),
             acf_raw_value) %>%
  remove_rownames()

### calculate rlm-calibrated data p value ###
p_rlm_value <-
  data_calibrate_rlm %>%
  mutate(intensity = ifelse(is.na(intensity), .01, intensity)) %>%
  split(.$feature) %>%
  map_dfr(function(x){
    lm <- lm(intensity~injection_sequence,
             data.frame(intensity = x$intensity,
                        injection_sequence = x$injection_sequence)) %>%
      summary()
    lm$coefficients[2,4]
  }) %>% t()
p_rlm <-
  data.frame(feature = rownames(p_rlm_value),
             p_rlm_value) %>%
  remove_rownames()

### calculate rlm-calibrated data vcf value ###
acf_rlm_value <-
  data_calibrate_rlm %>%
  mutate(intensity = ifelse(is.na(intensity), .01, intensity)) %>%
  split(.$feature) %>%
  map_dfr(function(x){
    acf <- acf(x$intensity)
    acf[["acf"]][2]
  }) %>% t()

acf_rlm <-
  data.frame(feature = rownames(acf_rlm_value),
             acf_rlm_value) %>%
  remove_rownames()

### calibrate data with AR(1)
data_calibrate_ar1 <-
  data %>%
  mutate(intensity = ifelse(is.na(intensity), .01, intensity)) %>%
  split(.$feature) %>%
  map(function(x){
    ar1 <- x$intensity %>% ar.ols(order.max = 1)
  })

data_calibrate_ar1_value <-
  data_calibrate_ar1 %>%
  map_dfc(function(x){
    x$resid
  })

# data_calibrate_ar1 <- data.frame(data_calibrate_ar1_model$resid)

# inverse prediction
feature <- "F1"
ar_inverse_predict <- function(feature=...){
x <- feature %>%
  str_match("\\d{1,}") %>% as.character()
Fx.1 <-
  data %>%
  filter(feature == paste0("F",x),
         injection_sequence == 1) %>%
  pull(intensity)
Fx.2 <-
  data %>%
  filter(feature == paste0("F",x),
         injection_sequence == 2) %>%
  pull(intensity)
Ar.1 <-
  data_calibrate_ar1[[paste0("F",x)]][["ar"]] %>% as.numeric()
m <-
  data_calibrate_ar1[[paste0("F",x)]][["x.mean"]]
ex.1 <-
  Fx.2-m-Ar.1*Fx.1+Ar.1*m
ex.1 <- ifelse(length(ex.1)==0, 0, ex.1)
return(ex.1)
}

data_calibrate_ar1_value_Na <-
data_calibrate_ar1_value[1,] %>%
  pivot_longer(colnames(data_calibrate_ar1_value),
               names_to = "Feature",
               values_to = "Intensity") %>%
  split(.$Feature) %>%
  map_dfr(function(x){
      x %>% mutate(e1 = .$Feature %>% ar_inverse_predict)
      })

data_calibrate_ar1_value_e1 <-
  data_calibrate_ar1_value_Na %>%
  mutate(Intensity = ifelse(is.na(Intensity),
                            e1,
                            Intensity)) %>%
  select(-e1) %>%
  pivot_wider(names_from = Feature, values_from = Intensity)

col_max <-
  colnames(data_calibrate_ar1_value) %>%
  str_match("\\d{1,4}") %>%
  as.numeric() %>%
  max
data_calibrate_ar1_value[1,] <- data_calibrate_ar1_value_e1 %>% as.numeric()

data_calibrate_ar1_value_long <-
  data_calibrate_ar1_value %>%
  select(paste0("F",1:col_max)) %>%
  mutate(injection_sequence = data_calibrate_rlm_wide$injection_sequence,
         sample_name = data_calibrate_rlm_wide$sample_name) %>%
  pivot_longer(paste0("F",1:col_max),
               names_to = "feature",
               values_to = "intensity") %>%
  arrange(feature, injection_sequence)

### calculate AR(1)-calibrated data p value ###
p_ar1_value <-
  data_calibrate_ar1_value_long %>%
  mutate(intensity = ifelse(is.na(intensity), .01, intensity)) %>%
  split(.$feature) %>%
  map_dfr(function(x){
    lm <- lm(intensity~injection_sequence,
             data.frame(intensity = x$intensity,
                        injection_sequence = x$injection_sequence)) %>%
      summary()
    lm$coefficients[2,4]
  }) %>% t()
p_ar1 <-
  data.frame(feature = rownames(p_ar1_value),
             p_ar1_value) %>%
  remove_rownames()

### calculate AR(1)-calibrated data vcf value ###
acf_ar1_value <-
  data_calibrate_ar1_value_long %>%
  mutate(intensity = ifelse(is.na(intensity), .01, intensity)) %>%
  split(.$feature) %>%
  map_dfr(function(x){
    acf <- acf(x$intensity)
    acf[["acf"]][2]
  }) %>% t()

acf_ar1 <-
  data.frame(feature = rownames(acf_ar1_value),
             acf_ar1_value) %>%
  remove_rownames()



#####calibrate with I2-AR(1)#####
data_calibrate_i2ar1 <-
  data %>%
  mutate(intensity = ifelse(is.na(intensity), .01, intensity)) %>%
  mutate(intensity_square = intensity**2) %>%
  split(.$feature) %>%
  map(function(x){
    ar1 <- x$intensity %>% ar.ols(order.max = 1)
  })
data_calibrate_i2ar1_value <-
  data_calibrate_ar1 %>%
  map_dfc(function(x){
    x$resid
  })

