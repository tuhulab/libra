library(tidyverse)
data_table <- readr::read_csv("/Users/tuhu/Projects/multiomics-ad-metabolome-life/xcms_result/") %>%
              mutate(feature = paste0("F",rownames(.)))
sample_table <- readr::read_csv("/Users/tuhu/Projects/multiomics-ad-metabolome-life/xcms_result/Life_Serum_pos_param1_samplelist.csv")

multiomics_ad_serum_pos <- data_table %>% pivot_longer(colnames(.) %>%
                                           str_match("X\\d{1,}") %>%
                                           discard(is.na),
                                           names_to = "code",
                                           values_to = "intensity") %>% left_join(sample_table) %>%
             filter(sample_name != "Blank",
                    sample_name != "Blank_IS",
                    sample_name != "Assay_blank",
                    sample_name != "MetStd") %>%
             group_by(feature) %>%
             mutate(injection_sequence = 1:90) %>%
             select(feature, intensity, sample_name, injection_sequence)
save(multiomics_ad_serum_pos, file = "data/multiomics_ad_serum_pos")

load("data/multiomics_ad_serum_pos")
Feature <- "F2"

####plotFeature####
intensity <- multiomics_ad_serum_pos %>%
                filter(feature == Feature) %>%
                pull(intensity)
injection_sequence <- multiomics_ad_serum_pos %>%
                         filter(feature == Feature) %>%
                         pull(injection_sequence)
lm <- lm(intensity ~ injection_sequence) %>%
  summary()

lm$coefficients[2,4] #p-value

multiomics_ad_serum_pos %>% filter(feature == Feature) %>%
  ggplot(aes(injection_sequence, intensity)) +
  geom_point() +
  geom_abline(slope = lm$coefficients[2,1],
              intercept = lm$coefficients[1,1]) +
  ggtitle(paste("Feature -", Feature, "     ",
                "p-value = ", lm$coefficients[2,4]))


libra::plotFeature(multiomics_ad_serum_pos, "F2")

###fit rlm####
intensity <- multiomics_ad_serum_pos %>%
                filter(feature == Feature) %>%
                pull(intensity)
injection_sequence <- multiomics_ad_serum_pos %>%
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

calibrated_data <-
calibrateBatch(multiomics_ad_serum_pos,
               Feature)
plotFeature(multiomics_ad_serum_pos,
            Feature)
plotFeature(calibrated_data,
            Feature)


multiomics_ad_serum_pos %>% view()




calcBatchEffect <- function(data = ...,
                            Feature = ...){
  require(tidyverse)
  feature_intensity <- data %>% filter(feature == Feature) %>% pull(intensity)
  injection_sequence <- data %>% filter(feature == Feature) %>% pull(injection_sequence)
  rlm <- MASS::rlm(feature_intensity ~ injection_sequence) %>% summary(correlation = TRUE)
  r_cor_cof <- rlm[["correlation"]][1,2]
  return(r_cor_cof)
}

i <- data %>% pull(feature) %>% str_match("\\d{1,}") %>% as.numeric() %>% max()
cor_cof <-
sapply(1:i, function(x){
  calcBatchEffect(data = multiomics_ad_serum_pos, paste0("F", x))
})
