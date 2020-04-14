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
    data_n_impute <- data_n %>% dplyr::mutate(data_impute_na = purrr::map(data, impute_na_one)) %>% dplyr::select(-data) %>% tidyr::unnest(cols = c(data_impute_na))
    return(data_n_impute)
}




