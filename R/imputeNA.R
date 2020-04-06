#' Impute missing data of metabolomics experiment
#'
#' @param data Tidy format dataframe consisting cols intensity and feature
#'
imputeNA.old <- function(data = ...) {
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

#' Impute missing data of metabolomics experiment
#'
#' @param data Tidy format dataframe consisting cols intensity and feature
#' @param feature The column name of feature (by default feature)
#' @param intensity The column name of intensity (by default intensity)
imputeNA <- function(data=...,
                     feature=feature,
                     intensity=intensity) {
    feature <- rlang::enexpr(feature)
    intensity <- rlang::enexpr(intensity)
    data_n <- data %>% dplyr::group_by(!!feature) %>% tidyr::nest()

    impute_na_one <- function(data. = ...) {
        test_na <- data. %>% select(!!intensity) %>% is.na() %>% sum()
        ifelse(test_na == 0, return(data.), {
            na_index <- data. %>% select(!!intensity) %>% is.na()
            intensity_min <- data. %>% pull(!!intensity) %>% min(na.rm = TRUE)

            intensity_n <- data. %>% pull(!!intensity)
            intensity_n[na_index] <- runif(na_index %>% sum, 0, intensity_min * 0.67)
            data. <- data. %>% mutate(!!intensity := intensity_n)
            return(data.)
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
