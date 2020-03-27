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

####read XCMS result####
data_table <- readr::read_csv("/Users/tuhu/Projects/multiomics-ad-metabolome-life/xcms_result/Life_Serum_neg_param1_peaklist.csv") %>%
  mutate(feature = paste0("F",rownames(.)))
sample_table <- readr::read_csv("/Users/tuhu/Projects/multiomics-ad-metabolome-life/xcms_result/Life_Serum_neg_param1_samplelist.csv")
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

data_imputeNA <- data %>% libra::imputeNA()

#### calculate CVQC ####
QC_index <- data_imputeNA %>% dplyr::filter(feature=="F1") %>% dplyr::pull(sample_name) %>% stringr::str_detect("Global_pool_IS") %>% which()
data_imputeNA_n <- data_imputeNA %>% dplyr::group_by(feature) %>% tidyr::nest()

calcBatchEffect.CVQC_one <- function(data=..., QC_index=...){
  intensity_QC_index <- data_imputeNA_n[[2]][[1]][QC_index,] %>% pull(intensity)
  CVQC <- sd(intensity_QC_index)/mean(intensity_QC_index)
  return(CVQC)
}

calcBatchEffect.CVQC_one(data_imputeNA_n[[2]][[1]],
                         QC_index)

data <- data_imputeNA

calcBatchEffect.CVQC <- function(data = ...,
                                 QC_index = ...) {
  data_n <- data %>% dplyr::group_by(feature) %>% tidyr::nest()

  calcBatchEffect.CVQC_one <- function(data=...){
    intensity_QC_index <- data[QC_index,] %>% pull(intensity)
    CVQC <- sd(intensity_QC_index)/mean(intensity_QC_index)
    return(CVQC)
  }
  data_n_CVQC <- data_n %>% dplyr::mutate(CVQC = map_dbl(data, calcBatchEffect.CVQC_one)) %>% tidyr::unnest(cols = c(data))
  return(data_n_CVQC)
}
