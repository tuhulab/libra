#' Inter batch calibration by robust linear modelling (rlm)
#'
#' @param data metabolomics data in tidy long-format
#' @param feature column name of feature (by default feature)
#' @param batch column name of batch (by default batch)
#' @param intensity column name of intensity (by default intensity)

calibrateBatch.inter.rlm <- function (data = ...,
                                      feature = feature,
                                      batch = batch,
                                      intensity = intensity)
{
  feature <- rlang::enexpr(feature)
  batch <- rlang::enexpr(batch)
  intensity <- rlang::enexpr(intensity)
  data_n <- data %>% dplyr::group_by(!!batch, !!feature) %>%
    tidyr::nest() %>% dplyr::mutate(data_intra_calibrate = purrr::map(data,
                                                                      libra::calibrateBatch.intra.rlm.single.feature,
                                                                      intensity = !!intensity))
  data_n_intra <- data_n %>% dplyr::mutate(intra_batch_center = data_intra_calibrate %>%
                                             unlist %>% mean) %>% dplyr::ungroup()
  data_n_inter <- data_n_intra %>% dplyr::group_by(!!feature) %>%
    dplyr::summarize(multi_batch_center = mean(intra_batch_center))
  data_n_intra_inter <- data_n_intra %>% dplyr::left_join(data_n_inter,
                                                          by = c(rlang::as_string(feature))) %>%
    dplyr::mutate(factor = intra_batch_center/multi_batch_center) %>%
    tidyr::unnest(cols = c(data, data_intra_calibrate)) %>%
    dplyr::mutate(intensity_intra_inter_calibrated = intensity_intra_calibrated/factor)
  data_n_intra_inter_output <- data_n_intra_inter %>% dplyr::select(-intensity_intra_calibrated,
                                                                    -intra_batch_center, -multi_batch_center, -factor) %>%
    dplyr::arrange(!! feature, !! batch, injection_sequence)
  return(data_n_intra_inter_output)
}

calibrateBatch.inter.rlm.group <- function(data = ...,
                                           feature = feature,
                                           batch = batch,
                                           intensity = intensity,
                                           group = ...){
  feature <- rlang::enexpr(feature)
  batch <- rlang::enexpr(batch)
  intensity <- rlang::enexpr(intensity)
  group_factor <- rlang::enexpr(group)
  data_n <- data %>% dplyr::group_by(!!batch, !!feature, !!group_factor) %>%
    tidyr::nest() %>% dplyr::mutate(data_intra_calibrate = purrr::map(data,
                                                                      libra::calibrateBatch.intra.rlm.single.feature, intensity = !!intensity))
  data_n_intra <- data_n %>% dplyr::mutate(intra_batch_center = data_intra_calibrate %>%
                                             unlist %>% mean) %>% dplyr::ungroup()
  data_n_inter <- data_n_intra %>% dplyr::group_by(!!feature, !! group_factor) %>%
    dplyr::summarize(multi_batch_center = mean(intra_batch_center))
  data_n_intra_inter <- data_n_intra %>% dplyr::left_join(data_n_inter,
                                                          by = c(rlang::as_string(feature), rlang::as_string(group_factor))) %>% dplyr::mutate(factor = intra_batch_center/multi_batch_center) %>%
    tidyr::unnest(cols = c(data, data_intra_calibrate)) %>%
    dplyr::mutate(intensity_intra_inter_calibrated = intensity_intra_calibrated/factor)
  data_n_intra_inter_output <- data_n_intra_inter %>% dplyr::select(-intensity_intra_calibrated,
                                                                    -intra_batch_center, -multi_batch_center, -factor) %>%
    dplyr::arrange(!! feature, !! batch, injection_sequence)
  return(data_n_intra_inter_output)
}
