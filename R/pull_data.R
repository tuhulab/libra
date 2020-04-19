#' Pull pre-processed data and summary statistics based on the regular expression (naming pattern)
#'
#' @param pattern regular expression in sample_table
#' @param data_table. data table
#' @param sample_table. sample information table
pull.data <- function(pattern = ...,
                      var.pattern = filename,
                      data_table. = data_table,
                      sample_table. = sample_table){
    var.pattern = enexpr(var.pattern)
    datalist <- list()
    datalist$pattern <- pattern
    datalist$index <- sample_table. %>% pull(!! var.pattern) %>% stringr::str_detect(pattern) %>% which()
    xn_index <- data_table. %>% colnames %>% stringr::str_detect("X\\d{1,}")
    non_xn_index <- which(!xn_index)
    datalist$data <- data_table. %>% dplyr::select(non_xn_index, paste0("X", datalist$index))
    datalist$intensity <- data_table. %>% dplyr::select(paste0("X", datalist$index))
    datalist$intensity_mean <- datalist$intensity %>% rowMeans(na.rm = TRUE)
    datalist$prevalence <- 1 - (datalist$intensity %>% is.na() %>% rowMeans())
    datalist$intensity[is.na(datalist$intensity)] <- 0
    return(datalist)
}
