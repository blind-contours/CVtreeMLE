library(ggplot2)
library(rgl)
library(viridis)

simulate_2d_data <- function(n_obs = 3000,
                          mu = c(2.5, 3),
                          sigma_mod = matrix(c(1, 0.65, 0.65, 1),
                                             nrow = 2,
                                             ncol = 2))
  {

  mixtures <- MASS::mvrnorm(n = 10000,
                            mu = mu,
                            Sigma = sigma_mod)

  mixtures[mixtures < 0] <- 0.001
  colnames(mixtures) <- c("M1", "M2")
  mixtures <- as.data.frame(mixtures)

  age <- rnorm(n_obs, 37, 3)
  bmi <- rnorm(n_obs, 20, 1)
  sex <- rbinom(n_obs, size = 1, prob = 0.5)

  indicator <- rbinom(n_obs,
                      size = 1,
                      prob = 1 / (1+ exp(0.07 + (0.03 *
                                                   bmi) +
                                               (0.08 *sex))))

  m1_thresh_stds <- abs(0.02 * age/5)
  m2_thresh_stds <- abs(0.05 * bmi/5)

  m1_thresholds <- rnorm(n = length(m1_thresh_stds),
                         mean = rep(2.5, length(m1_thresh_stds)),
                         sd = m1_thresh_stds)

  m2_thresholds <- rnorm(n = length(m2_thresh_stds),
                         mean = rep(3.5, length(m2_thresh_stds)),
                         sd = m2_thresh_stds)

  data <- cbind.data.frame(age,
                           bmi,
                           sex,
                           indicator,
                           m1_thresholds,
                           m2_thresholds)

  ms <- as.data.frame(matrix(data = NA, ncol = 2, nrow = n_obs))
  colnames(ms) <- c("M1", "M2")

  for (i in 1:nrow(data)) {
    row <- data[i,]
    if (row$indicator == 1) {
      subspace <- subset(mixtures, M1 > row$m1_thresholds & M2 > row$m2_thresholds)
      sample <- subspace[sample(nrow(subspace), 1, replace = FALSE), ]
    }else{
      subspace <- subset(mixtures, !(M1 > row$m1_thresholds & M2 > row$m2_thresholds))
      sample <- subspace[sample(nrow(subspace), 1, replace = FALSE), ]
    }

    ms[i,] <- sample
  }

  data_w_mixtures <- cbind.data.frame(data, ms)

  y <- (data_w_mixtures$indicator * (2/(0.3+exp(-data_w_mixtures$M1-2.5))) *
          (3/(0.3+exp(-data_w_mixtures$M2-2.5)))) +
    (abs(data_w_mixtures$indicator-1)* (2/(0.1+exp(-data_w_mixtures$M1-3)))) +
    (abs(data_w_mixtures$indicator-1)* (3/(0.1+exp(-data_w_mixtures$M2-4)))) +
    2 * data_w_mixtures$sex +
    rnorm(n = nrow(data_w_mixtures), mean = 0, sd = 4)

  data_w_mixtures$y <- y
  ind <- ifelse(data_w_mixtures$M1 > 2.361 &data_w_mixtures$M1 < 6.086 &
                  data_w_mixtures$M2 > 3.502 & data_w_mixtures$M2 < 6.593, 1, 0)

  region_y <- (ind * (2/(0.3+exp(-data_w_mixtures$M1-2.5))) *
                 (3/(0.3+exp(-data_w_mixtures$M2-2.5))))

  non_region_y <- (abs(ind-1)* (2/(0.1+exp(-data_w_mixtures$M1-3)))) +
    (abs(data_w_mixtures$indicator-1)* (3/(0.1+exp(-data_w_mixtures$M2-4))))

  mean_region_y <- mean(region_y[region_y!=0])
  mean_nonregion_y <- mean(non_region_y[non_region_y!=0])

  mean_region_diff <- mean_region_y - mean_nonregion_y

  return(list("data" = data_w_mixtures, "truth" = mean_region_diff
              ))
}

sim_data <- simulate_2d_data(n_obs = 5000)

ggplot(sim_data$data, aes(M1, M2))+
  geom_point(aes(color = y)) +
  scale_color_viridis(discrete = FALSE, option = "D")+
  scale_fill_viridis(discrete = FALSE) +
  theme_minimal() +
  theme(legend.position = "bottom")

scatter3d(x = sim_data$data$M1,
          y = sim_data$data$y,
          z = sim_data$data$M2,
          groups = as.factor(sim_data$data$indicator),
          grid = FALSE,
          fit = "smooth",
          xlab = "Exposure 1",
          ylab = "Outcome",
          zlab = "Exposure 2",
          residuals = FALSE)


sim_results <- CVtreeMLE(data = sim_data$data,
                         w = c("age", "bmi", "sex"),
                         a = c(paste("M", seq(2), sep = "")),
                         y = "y",
                         n_folds = 2,
                         parallel_cv = TRUE,
                         seed = 2333,
                         parallel_type = "multi_session",
                         family = "gaussian",
                         num_cores = 6)
