imputeNA <- function(data = ...) {
    data_n <- data %>% dplyr::group_by(feature) %>% tidyr::nest()
    impute_na_one <- function(data = ...) {
        test_na <- data$intensity %>% is.na() %>% sum()
        ifelse(test_na == 0, return(data), {
            na_index <- data$intensity %>% is.na()
            intensity_min <- data$intensity %>% min(na.rm = TRUE)
            data$intensity[na_index] <- runif(na_index %>% sum, 0, intensity_min * 0.67)
            return(data)
        })
    }
    data_n_impute <- data_n %>% dplyr::mutate(data_impute_na = map(data, impute_na_one)) %>% dplyr::select(-data) %>% tidyr::unnest(cols = c(data_impute_na))
    return(data_n_impute)
}

#' Detect early and late eluting compounds (early- and late- birds)
#'
#' @param rt retention time in the numeric format
#' @param limit_low lower limit of features to be removed
#' @param limit_high higher limit of features to be removed
#' @return \code{rt} between \code{limit_low} and \code{limit-high}
detect.birds <- function(rt = ..., limit_low = 0.3, limit_high = 6.5) {
    rt.df <- data.frame(rt) %>% mutate(feature = paste0("F", 1:length(rt))) %>% mutate(bird = ifelse((rt > limit_low) & (rt < limit_high), TRUE, FALSE))
    birds_id <- rt.df %>% pull(bird) %>% as.logical()
    return(rt.df)
}
