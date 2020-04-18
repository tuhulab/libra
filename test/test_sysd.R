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
data_t_p_na_interc <- calibrateBatch.inter.rlm(data = data_t_p_na, feature = Feature, batch = Plate, intensity = Intensity)
data_t_n_na_interc <- calibrateBatch.inter.rlm(data = data_t_n_na, feature = Feature, batch = Plate, intensity = Intensity)

readr::write_csv(data_t_p_na_interc, "test/SysDiet_pos_interc.csv")
readr::write_csv(data_t_n_na_interc, "test/SysDiet_neg_interc.csv")


# output ----------------------------
data_l_p_c <- data_t_p_na_interc %>% select(-Intensity, -intensity_intra_calibrated, -intra_batch_center, -center_intensity, -multi_batch_center, -factor) %>%
    rename(Intensity_calibrated = intensity_intra_inter_calibrated)
data_l_n_c <- data_t_n_na_interc %>% select(-Intensity, -intensity_intra_calibrated, -intra_batch_center, -center_intensity, -multi_batch_center, -factor) %>%
    rename(Intensity_calibrated = intensity_intra_inter_calibrated)

data_w_p_c <- data_l_p_c %>% select(Feature, Intensity_calibrated,
                      File_name, Subject, Samples, Visit, Diet, Center, S_and_PS) %>%
               pivot_wider(names_from = Feature, values_from = Intensity_calibrated)

data_w_n_c <- data_l_n_c %>% select(Feature, Intensity_calibrated,
                                    File_name, Subject, Samples, Visit, Diet, Center, S_and_PS) %>%
  pivot_wider(names_from = Feature, values_from = Intensity_calibrated)

readr::write_csv(data_l_p_c, "test/SysDiet/SysDiet_pos_interc_long.csv")
readr::write_csv(data_l_n_c, "test/SysDiet/SysDiet_neg_interc_long.csv")

readr::write_csv(data_w_p_c, "test/SysDiet/SysDiet_pos_interc_wide.csv")
readr::write_csv(data_w_n_c, "test/SysDiet/SysDiet_neg_interc_wide.csv")

readr::write_csv(data_l_p_c %>% select(Feature, mz, rt, group_camera) %>% distinct(), "test/SysDiet/SysDiet_pos_featureinfo.csv")
readr::write_csv(data_l_n_c %>% select(Feature, mz, rt, group_camera) %>% distinct(), "test/SysDiet/SysDiet_neg_featureinfo.csv")

### test -----------------------------
feature <- rlang::expr(Feature)
batch <- rlang::expr(Plate)
intensity <- rlang::expr(Intensity)
data <- data_t_p_na

data_n <- data %>%
  group_by(!! feature, !! batch) %>% nest()


data_n_1 <- data_n %>% mutate(data_intra_calibrate = map(data, libra::calibrateBatch.intra.rlm, intensity = !!intensity))

