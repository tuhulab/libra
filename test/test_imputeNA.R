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

data_n <- data %>% dplyr::group_by(feature) %>% tidyr::nest()
impute_na_one <- function(data=...){
  test_na <- data$intensity %>% is.na() %>% sum()
  ifelse(test_na == 0,
         return(data),
         {na_index <- data$intensity %>% is.na()
         intensity_min <- data$intensity %>% min(na.rm = TRUE)
         data$intensity[na_index] <- runif(na_index %>% sum, 0, intensity_min)
         return(data)}
  )
}
data_n_impute <- data_n %>% dplyr::mutate(data_impute_na = map(data, impute_na_one)) %>% dplyr::select(-data) %>%
    tidyr::unnest(cols = c(data_impute_na))


data_imputeNA <- data %>% libra::imputeNA()
