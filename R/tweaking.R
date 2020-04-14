# Tweaking functions -----------------------------


rm_sample <- function(data.df = ...,
                      sample.df = ...,
                      sample_rm = ...){
  data <-
    data.df %>% dplyr::mutate(feature = paste0("F", 1:nrow(.))) %>%
    pivot_longer(colnames(.) %>% stringr::str_match("X\\d{1,}") %>% purrr::discard(is.na), names_to = "code", values_to = "intensity") %>%
    left_join(sample.df, by = "code") %>%
    filter(! sample_name %in% sample_rm) %>%
    group_by(feature) %>% mutate(injection_sequence = 1:sum(.$feature == "F1")) %>%
    select(feature, code, mz, rt, pcgroup, adduct, id_predret, id_kudb, intensity, sample_name, injection_sequence)
  return(data)
}


rm_feature.bird <- function(data.df = ...,
                            col_rt = rt,
                            bird.limit.low = 0.3,
                            bird.limit.high = 6.5){
  col_rt <- rlang::enexpr(col_rt)
  data <- data.df %>% mutate(bird = ifelse(!! col_rt > bird.limit.high | col_rt < bird.limit.low, # Very important to !!unquote!!
                                           TRUE,
                                           FALSE)) %>%
    filter(bird == FALSE) %>% select(-bird)
  return(data)
}

#' The old method to detect early and late eluting compounds (early- and late- birds)
#'
#' @param rt retention time in the numeric format
#' @param limit_low lower limit of features to be removed
#' @param limit_high higher limit of features to be removed
#' @return \code{rt} between \code{limit_low} and \code{limit-high}
detect.birds.old <- function(rt = ...,
                             limit_low = 0.3,
                             limit_high = 6.5) {
  rt.df <- data.frame(rt) %>% mutate(feature = paste0("F", 1:length(rt))) %>% mutate(bird = ifelse(rt > limit_low|rt < limit_high, TRUE, FALSE))
  birds_id <- rt.df %>% pull(bird) %>% as.logical()
  return(rt.df)
}

#' Clean features based on blank samples
#'
#' @param data.df Data table
#' @param pattern.blank The pattern(regex) of blank samples for libra::pull.data
#' @param pattern.pool The pattern(regex) of pooled samples for libra::pull.data
#' @param threshold.prevalence Features have the prevalence higher than this threshold will be marked as noise and removed
#' @param threshold.intensity.fraction Features have the intensity higher than the fraction (blank mean intensity >  fraction * mean intensity of pooled samples) will be marked as noise and removed
rm_feature.blank <- function(data.df = data_table,
                             pattern.blank = "Blank$|Assay blank",
                             pattern.pool = "Global pool IS$",
                             threshold.prevalence = 0.6,
                             threshold.intensity.fraction = 0.67){
    data.blank <- libra::pull.data(pattern.blank)
    data.pool <- libra::pull.data(pattern.pool)
    noise_index <- (data.blank$prevalence >= threshold.prevalence) & (data.blank$intensity_mean > (threshold.intensity.fraction * data.pool$intensity_mean))
    message(paste("Fraction of Remove Feature (based on blank):", round(mean(noise_index),2)))
    return(data.df[!noise_index,])
}
