# Test SysDiet data -----------------
library(tidyr)
library(dplyr)
library(purrr)
library(stringr)


sys.d <- readr::read_csv("data/Sysdiet_plasma_wo_bc.csv")

# Extract rt and mz
rt <- sys.d[1,] %>% as.numeric() %>% discard(is.na)
mz <- sys.d[2,] %>% as.numeric() %>% discard(is.na)
polarity <- sys.d[6,-1] %>% discard(is.na) %>% as.character()
group_camera <- sys.d[23,-1] %>% discard(is.na) %>% as.character() %>% str_remove("Class ")

sys.d_1 <- sys.d[-1:-24,-1]
column_names_raw <- sys.d_1 %>% colnames()
column_names_new <- c(column_names_raw[1:9] %>% str_replace_all(" ","_"),paste0("X",1:length(rt)))
colnames(sys.d_1) <- column_names_new

sys.d_2 <- sys.d_1 %>% arrange(File_name) %>%
    mutate(Samples = ifelse(Samples == "ps 1", "ps1", Samples))# a minor change, change "ps 1" to "ps1"

index_std <- sys.d_2$Samples %>% str_detect("std") %>% which()

sys.d_2_n <- sys.d_2[-index_std,] %>% group_by(Plate) %>% nest() ; mutate.injection_sequence <- function(df=...){mutate(df, injection_sequence = 1:nrow(df))}

sys.d_3 <- sys.d_2_n %>% mutate(data_n = map(data, mutate.injection_sequence)) %>% select(-data) %>% unnest(cols = data_n) %>%
    select(column_names_new[1:9], injection_sequence, paste0("X",1:length(rt)))

feature_info <- data.frame(rt,
                           mz,
                           group_camera,
                           polarity,
                           Feature=paste0("X",1:length(rt)))

data_tidy <- sys.d_3 %>%
  pivot_longer(cols = paste0("X",1:length(rt)), names_to = "Feature", values_to = "Intensity") %>%
  left_join(feature_info) %>%
  mutate(Intensity = ifelse(Intensity == 0, NA, Intensity %>% as.numeric())) #change 0 to NA

data_tidy_pos <- data_tidy %>% filter(polarity == "pos")
data_tidy_neg <- data_tidy %>% filter(polarity == "neg")

# impute NA (0) ----------------------
data_t_p_na <- data_tidy_pos %>% libra::imputeNA(feature = Feature, intensity = Intensity)
data_t_n_na <- data_tidy_neg %>% libra::imputeNA(feature = Feature, intensity = Intensity)

# calibrate inter batch --------------
data_t_p_na_n <- data_t_p_na %>% group_by(Feature, Plate) %>% nest()

# start from intra batch effect ------

# reshape data
data_reshape_n <- data_t_p_na %>% mutate(intensity = Intensity) %>% group_by(Feature, Plate) %>% nest()

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

data_reshape_n_intra <- data_reshape_n %>% mutate(data_intra_calibrate = map(data, calibrateBatch.intra.rlm)) %>% select(-data)
data_reshape_n_intra_inter <- data_reshape_n_intra %>% mutate(intra_batch_center=map(data_intra_calibrate, function(x)x[["center_intensity"]] %>% unique)) %>% ungroup()

multi_batch_center <- data_reshape_n_intra_inter %>% mutate(intra_batch_center = as.numeric(intra_batch_center)) %>% group_by(Feature) %>% summarize(multi_batch_center = mean(intra_batch_center))

calibration_matrix <- data_reshape_n_intra_inter %>% mutate(intra_batch_center = as.numeric(intra_batch_center)) %>% left_join(multi_batch_center) %>% mutate(facotr = intra_batch_center/multi_batch_center)

data_reshape_n_intra_inter_indi <- calibration_matrix %>% unnest()

indi_intensity <- data_reshape_n_intra_inter_indi %>% mutate(intensity_final = intensity_calibrated/facotr)
