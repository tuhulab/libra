#' Inter batch calibration by robust linear modelling (rlm)
#'
#' @param data metabolomics in tidy long-format
#' @param feature column name of feature (by default feature)
#' @param batch column name of batch (by default batch)
#' @param intensity column name of intensity (by default intensity)

calibrateBatch.inter.rlm <- function(data = ...,
                                     feature = feature,
                                     batch = batch,
                                     intensity = intensity){
  feature <- rlang::enexpr(feature)
  batch <- rlang::enexpr(batch)
  intensity <- rlang::enexpr(intensity)

  data_n <- data %>%
    group_by(!! feature, !! batch) %>% nest() %>%
    mutate(data_intra_calibrate = map(data, libra::calibrateBatch.intra.rlm, intensity= !! intensity)) %>%
    select(-data)

  data_n_intra <- data_n %>%
    mutate(intra_batch_center=map(data_intra_calibrate, function(x)x[["center_intensity"]] %>% unique)) %>%
    ungroup() %>% mutate(intra_batch_center = as.numeric(intra_batch_center))

  data_n_inter <- data_n_intra %>%
    group_by(!! feature) %>% summarize(multi_batch_center = mean(intra_batch_center))

  data_n_intra_inter <- data_n_intra %>%
    left_join(data_n_inter, by = rlang::as_string(feature)) %>% mutate(factor = intra_batch_center / multi_batch_center) %>%
    unnest(cols = c(data_intra_calibrate)) %>% mutate(intensity_intra_inter_calibrated = intensity_intra_calibrated / factor)
}
