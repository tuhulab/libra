#' Plot features by injection sequence

#' @param data a dataframe consisting intensity, feature and injection sequence
#' @param Feature

plotFeature <- function(data = ..., Feature = ..., intensity = data$intensity) {
    intensity <- data %>% dplyr::filter(feature == Feature) %>% dplyr::pull(intensity)
    injection_sequence <- data %>% dplyr::filter(feature == Feature) %>% dplyr::pull(injection_sequence)
    lm <- stats::lm(intensity ~ injection_sequence) %>% summary()
    data %>% dplyr::filter(feature == Feature) %>% ggplot2::ggplot(aes(injection_sequence, intensity)) + ggplot2::geom_point() + geom_abline(slope = lm$coefficients[2, 
        1], intercept = lm$coefficients[1, 1]) + ggtitle(paste("Feature -", Feature, "     ", "p-value = ", lm$coefficients[2, 4]))
}


plotFeature.method <- function(data = ..., Feature = ..., intensity = data$intensity, calibration_method = ...) {
    intensity <- data %>% filter(feature == Feature) %>% pull(intensity)
    injection_sequence <- data %>% filter(feature == Feature) %>% pull(injection_sequence)
    lm <- lm(intensity ~ injection_sequence) %>% summary()
    data %>% dplyr::filter(feature == Feature) %>% ggplot(aes(injection_sequence, intensity)) + geom_point() + geom_abline(slope = lm$coefficients[2, 
        1], intercept = lm$coefficients[1, 1]) + ggtitle(paste("Feature -", Feature, "     ", "p-value = ", lm$coefficients[2, 4])) + annotate(geom = "text", 
        x = max(injection_sequence) * 0.9, y = max(intensity), label = paste0("calibrate=", calibration_method))
}
