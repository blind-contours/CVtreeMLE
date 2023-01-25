library(stringr)
library(zoo)
library(imputeTS)
library(here)

nhanes_data <- readRDS(here("sandbox/NHANES/output/NHANES_data.RDS"))

nhanes_data_telomere <- nhanes_data[!is.na(nhanes_data$mean_telomere),]
lead_exposure_data <- nhanes_data_telomere[!is.na(nhanes_data_telomere$lead),]

nhanes_data_telomere_na_thresh <- lead_exposure_data[ , colSums(is.na(lead_exposure_data)) < nrow(lead_exposure_data) ]
nhanes_data_telomere_na_thresh <- as.data.frame(unclass(nhanes_data_telomere_na_thresh),stringsAsFactors=TRUE)


metals <- c("barium", "cadmium", "cobalt", "cesium", "molybdenum", "lead",
            "antimony", "thallium", "tungsten")

nhanes_data_telomere_na_thresh[, metals] <- na_mean(nhanes_data_telomere_na_thresh[, metals])                                # Replace NA in all columns


outcome <- "mean_telomere"

covariates <- c("age_screen_years", "gender", "race", "educ_level",
                "marital_status", "alc_gm_1999", "cotinine_ng_ml", "bmi_kg_m2",
                "fam_poverty_ratio", "systolic_1", "diastolic_1", "birth_country",
                "avg_daily_physical_act", "muscle_training",
                "vigorous_intense_30_days", "fasting_glucose_mg_dl",
                "two_year_exam_weight", "two_year_interview_weight")


nhanes_results_neg <- CVtreeMLE(data = as.data.frame(nhanes_data_telomere_na_thresh),
                           w = covariates,
                           a = metals,
                           y = "mean_telomere",
                           fit_marginals = FALSE,
                           n_folds = 10,
                           direction = "negative",
                           seed = 143421,
                           parallel_cv = TRUE,
                           parallel = TRUE,
                           family = "continuous",
                           num_cores = 7,
                           max_iter = 10)

nhanes_results_neg$`Pooled TMLE Mixture Results` %>%
  dplyr::filter(Proportion_Folds >= 0.7)

mixture_plots <- plot_mixture_results(
  v_intxn_results = nhanes_results_neg$`V-Specific Mix Results`,
  hjust = 0.8)
mixture_plots$`cadmium-thallium`


qcomp <- qgcomp(Y~X1*X2*X3*X4*X5*X6*X7+Z+Z2+Z3, expnms=c(paste("X", seq(1,7), sep = "")),
              data = niehs_data,q=4, degree = 2, B=10)
plot(qcomp)


