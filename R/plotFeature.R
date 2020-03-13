plotFeature <- function(data = ..., Feature = ..., intensity = data$intensity) {
    require(tidyverse)
    require(grid)
    intensity <- data %>% filter(feature == Feature) %>% pull(intensity)
    injection_sequence <- data %>% filter(feature == Feature) %>% pull(injection_sequence)
    lm <- lm(intensity ~ injection_sequence) %>% summary()
    data %>% filter(feature == Feature) %>% ggplot(aes(injection_sequence, intensity)) + geom_point() + geom_abline(slope = lm$coefficients[2, 
        1], intercept = lm$coefficients[1, 1]) + ggtitle(paste("Feature -", Feature, "     ", "p-value = ", lm$coefficients[2, 
        4]))
}

plotFeature.method <- function(data = ..., Feature = ..., intensity = data$intensity, calibration_method = ...) {
    require(tidyverse)
    intensity <- data %>% filter(feature == Feature) %>% pull(intensity)
    injection_sequence <- data %>% filter(feature == Feature) %>% pull(injection_sequence)
    lm <- lm(intensity ~ injection_sequence) %>% summary()
    data %>% filter(feature == Feature) %>% ggplot(aes(injection_sequence, intensity)) + geom_point() + geom_abline(slope = lm$coefficients[2, 
        1], intercept = lm$coefficients[1, 1]) + ggtitle(paste("Feature -", Feature, "     ", "p-value = ", lm$coefficients[2, 
        4])) + annotate(geom = "text", x = max(injection_sequence) * 0.9, y = max(intensity), label = paste0("calibrate=", 
        calibration_method))
}
