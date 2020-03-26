imputeNA <- function(data = ...) {
    data_n <- data %>% dplyr::group_by(feature) %>% tidyr::nest()
    impute_na_one <- function(data = ...) {
        test_na <- data$intensity %>% is.na() %>% sum()
        ifelse(test_na == 0, return(data), {
            na_index <- data$intensity %>% is.na()
            intensity_min <- data$intensity %>% min(na.rm = TRUE)
            data$intensity[na_index] <- runif(na_index %>% sum, 0, intensity_min*0.67)
            return(data)
        })
    }
    data_n_impute <- data_n %>% dplyr::mutate(data_impute_na = map(data, impute_na_one)) %>% dplyr::select(-data) %>% tidyr::unnest(cols = c(data_impute_na))
    return(data_n_impute)
}
