#' Pull pre-processed data and summary statistics based on the regular expression (naming pattern)
#'
#' @param pattern regular expression in sample_table
#' @param data_table. data table
#' @param sample_table. sample information table
pull.data <- function(pattern = ...,
                      data_table. = data_table,
                      sample_table. = sample_table){
    datalist <- list()
    datalist$pattern <- pattern
    datalist$index <- sample_table.$filename %>% stringr::str_detect(pattern) %>% which()
    datalist$data <- data_table. %>% dplyr::select(1:6, paste0("X", datalist$index))
    datalist$intensity <- data_table. %>% dplyr::select(paste0("X", datalist$index))
    datalist$intensity_mean <- datalist$intensity %>% rowMeans(na.rm = TRUE)
    datalist$prevalence <- 1 - (datalist$intensity %>% is.na() %>% rowMeans())
    datalist$intensity[is.na(datalist$intensity)] <- 0
    return(datalist)
}
