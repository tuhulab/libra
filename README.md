# libra
Repo for libra - the data toolbox for working with LC-MS based metabolomics data (cleaning, calibration and visualization)

## Install
devtools::install_github("tuhulab/libra")

## Before start - data transformation
Metabolomics data needs to be transformed as data_table and sample_sample to start libra calibration. These two tables can be easily generated from previous [XCMS/CAMERA workflow](https://nexs-metabolomics.gitlab.io/INM/coffee-walkthrough/06_CAMERA.html).

data_table should compromize columns: mz, rt, X1, X2, X3, ... Xn. "Xi" denotes each sample, which each row stores intensity of each feature.

sample_table should compromize column: sample_name, code, batch (if multiple batch is involved). 

### Two tables


## Workflow
	# Load data ----------------------------------
    library(tidyr)
    library(dplyr)
    project <- "my_project"
    biospecimen <- "urine"
    mode <- "pos"

    data_table <- readr::read_csv(file.path("data", project, biospecimen,
                                            paste0(paste(biospecimen, mode, "raw" ,sep = "_"),".csv")))
    sample_table <- readr::read_csv(file.path("data", project, biospecimen,
                                            paste0(paste(biospecimen, mode, "sampleinfo" ,sep = "_"),".csv"))) %>%
    mutate(sample_name = ifelse(sample_name == "AS21-1", "AD21-1", sample_name))

    # Define parameter ---------------------------
    sample_rm <- c("Blank", "Blank_IS", "Assay_blank", "MetStd", "Global_pool_x3_IS", "Global_pool_x10_IS") # for removing samples
    pattern.blank <- "Blank$|Assay blank"
    pattern.pool <- "Global pool IS$"

    # Remove features (Clean Blanks, early and late birds)  --------
    d_rmF <- data_table %>% libra::rm_feature.blank() %>% libra::rm_feature.bird()

    # Remove samples(rmS) ------------------------
    d_rmFS <- d_rmF %>% libra::rm_sample(sample.df = sample_table, sample_rm) %>% left_join(sample_table %>% select(code, batch))

    # Impute NAs ---------------------------------
    d_rmSFNa <- d_rmFS %>% libra::imputeNA()

    # Calibrate: batch effect (intra)-------------
    d_rmSFNaC <- d_rmSFNa %>% libra::calibrateBatch.intra.rlm()

    # Calibrate: batch effect (inter)-------------
    d_rmSFNaC <- d_rmSFNa %>% libra::calibrateBatch.inter.rlm()

    # Output data wide format---------------------
    d_output_w <- d_rmSFNaC %>% select(feature, mz, rt, pcgroup, adduct, id_predret, id_kudb, code, intensity_intra_inter_calibrated) %>%
    pivot_wider(names_from = code, values_from = intensity_intra_inter_calibrated)
    readr::write_csv(d_output_w, file.path("data", project, biospecimen, paste0(paste(biospecimen, mode, "calibrated", "wide", sep = "_"), ".csv")))
    readr::write_csv(sample_table %>% select(-mzMLs) %>% filter(code %in% (d_output_w %>% colnames() %>% stringr::str_match("X\\d{1,}") %>% purrr::discard(is.na))),
                    file.path("data", project, biospecimen, paste0(paste(biospecimen, mode, "sampleinfo", "calibrated",sep = "_"), ".csv")))

    # Output data long data ----------------------
    readr::write_csv(d_rmSFNaC %>% select(-intensity, -injection_sequence, -intensity_intra_calibrated, -center_intensity, -intra_batch_center, -multi_batch_center, -factor),
                    file.path("data", project, biospecimen,
                            paste0(paste(biospecimen, mode, "calibrated", "long",sep = "_"),".csv")))
