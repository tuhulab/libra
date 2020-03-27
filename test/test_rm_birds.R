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


rt <- data_table %>% pull(rt)

detect.birds <- function(rt=...,
                         limit_low = .3,
                         limit_high = 6.5){
    rt.df <- data.frame(rt) %>% mutate(feature=paste0("F",1:length(rt))) %>% mutate(bird=ifelse( (rt>limit_low)&(rt<limit_high), TRUE, FALSE))
    birds_id <- rt.df %>% pull(bird) %>% as.logical()
    return(rt.df)
}

rt_bird <- rt %>% detect.birds()

