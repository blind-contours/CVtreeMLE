library(tidyverse)
library(here)
library(ggpubr)
library(purrr)
library(tidyr)

n_obs <- c(200, 350, 500, 750, 1000, 1500, 2000, 3000)

sim_results_1 <- readRDS(
  here("sandbox/data/CVtreeMLE_run_1.rds")
)

sim_results_2 <- readRDS(
  here("sandbox/data/CVtreeMLE_run_2.rds")
)

sim_results_3 <- readRDS(
  here("sandbox/data/CVtreeMLE_run_3.rds")
)

sim_results_4 <- readRDS(
  here("sandbox/data/CVtreeMLE_run_4.rds")
)

sim_results <- rbind(sim_results_1, sim_results_2, sim_results_3)


sim_statistics <- sim_results %>%
  group_by(n_obs) %>%
  summarise(
    tmle_pooled_da_bias = abs(mean(tmle_pooled_da_bias)),
    tmle_pooled_gt_bias = abs(mean(tmle_pooled_gt_bias)),
    v_spec_da_mean_bias = abs(mean(v_spec_da_mean_bias)),
    v_spec_gt_mean_bias = abs(mean(v_spec_gt_mean_bias)),
    v_pooled_da_bias = abs(mean(v_pooled_da_bias)),
    v_pooled_gt_bias = abs(mean(v_pooled_gt_bias)),
    tmle_pooled_sd = sd(tmle_pooled_ate),
    v_spec_mean_sd = sd(v_spec_mean_ate),
    v_pooled_mean_sd = sd(v_pooled_ate),
    tmle_pooled_da_mse = tmle_pooled_da_bias^2 + tmle_pooled_sd^2,
    tmle_pooled_gt_mse = tmle_pooled_gt_bias^2 + tmle_pooled_sd^2,
    v_spec_da_mean_mse = v_spec_da_mean_bias^2 + v_spec_mean_sd^2,
    v_spec_gt_mean_mse = v_spec_gt_mean_bias^2 + v_spec_mean_sd^2,
    v_pooled_da_mse = v_pooled_da_bias^2 + v_pooled_mean_sd^2,
    v_pooled_gt_mse = v_pooled_gt_bias^2 + v_pooled_mean_sd^2,
    tmle_pooled_gt_coverage = mean(tmle_pooled_gt_coverage),
    tmle_pooled_da_coverage = mean(tmle_pooled_da_coverage),
    v_spec_mean_da_cov = mean(v_spec_mean_da_cov),
    v_spec_mean_gt_cov = mean(v_spec_mean_gt_cov),
    pooled_da_cov = mean(pooled_da_cov),
    pooled_gt_cov = mean(pooled_gt_cov),
    True_pos = mean(`True Pos`),
    True_neg = mean(`True Neg`),
    False_pos = mean(`False Pos`),
    False_neg = mean(`False Neg`),

  ) %>%
  mutate(
    sqrt_n_abs_tmle_pooled_da_bias = sqrt(n_obs) * tmle_pooled_da_bias,
    sqrt_n_abs_tmle_pooled_gt_bias = sqrt(n_obs) * tmle_pooled_gt_bias,
    sqrt_n_abs_v_spec_da_mean_bias = sqrt(n_obs) * v_spec_da_mean_bias,
    sqrt_n_abs_v_spec_gt_mean_bias = sqrt(n_obs) * v_spec_gt_mean_bias,
    sqrt_n_abs_v_pooled_da_bias = sqrt(n_obs) * v_pooled_da_bias,
    sqrt_n_abs_v_pooled_gt_bias = sqrt(n_obs) * v_pooled_gt_bias,
    n_tmle_pooled_da_mse = n_obs*tmle_pooled_da_mse,
    n_tmle_pooled_gt_mse = n_obs*tmle_pooled_gt_mse,
    n_v_spec_da_mean_mse = n_obs*v_spec_da_mean_mse,
    n_v_spec_gt_mean_mse = n_obs*v_spec_gt_mean_mse,
    n_v_pooled_da_mse = n_obs*v_pooled_da_mse,
    n_v_spec_mean_gt_cov = n_obs*v_spec_mean_gt_cov,
  ) %>%
  ungroup

sim_statistics_long <- sim_statistics %>%
  tidyr::gather(statistic, value, -c(n_obs))


make_sim_statistics_plot <- function(sim_statistics_long,
                                     stats,
                                     labels,
                                     color,
                                     title) {
  filtered_sim_stats <- sim_statistics_long %>%
    filter(
      statistic %in% stats
    )

  ggplot(filtered_sim_stats, aes(n_obs, value)) +
    geom_point(color = color) +
    geom_line(color = color, size = 1) +
    geom_point(color= color) +
    facet_wrap("statistic",
               scales = "free_y",
               labeller = as_labeller(labels)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle=45, vjust=1, hjust=1)
    ) +
    scale_x_sqrt(breaks=n_obs)  + labs(x = "Number Observations",
                                       title = title)

}

plot_labels <- c(
  "tmle_pooled_da_bias" = "Pooled TMLE Bias Compared to Data-Adaptive Rule ATE",
  "tmle_pooled_gt_bias" = "Pooled TMLE Bias Compared to Ground-Truth Rule ATE",
  "v_spec_da_mean_bias" = "Mean V-Specific TMLE Bias Compared to Data-Adaptive Rule ATE",
  "v_spec_gt_mean_bias" = "Mean V-Specific TMLE Bias Compared to Ground-Truth Rule ATE",
  "v_pooled_da_bias" = "Harmonic Mean Bias Compared to Data-Adaptive Rule ATE",
  "v_pooled_gt_bias" = "Harmonic Mean Bias Compared to Ground-Truth Rule ATE"
)

CVtreeMLE_bias_plot <- make_sim_statistics_plot(
  sim_statistics_long,
  stats = names(plot_labels),
  labels = plot_labels,
  color = "steelblue",
  title = "Bias Measures"
)

plot_labels <- c(
  "True_pos" = "Rule Indicator True Positive Compared to Ground-Truth",
  "True_neg" = "Rule Indicator True Negative Compared to Ground-Truth",
  "False_pos" = "Rule Indicator False Positive Compared to Ground-Truth",
  "False_neg" = "Rule Indicator False Negative Compared to Ground-Truth"
)

CVtreeMLE_rule_plot <- make_sim_statistics_plot(
  sim_statistics_long,
  stats = names(plot_labels),
  labels = plot_labels,
  color = "orange",
  title = "Rule Coverage Measures"
)

plot_labels <- c(
  "tmle_pooled_gt_coverage" = "Pooled TMLE Ground-Truth ATE Coverage",
  "tmle_pooled_da_coverage" = "Pooled TMLE Data-Adaptive ATE Coverage",
  "v_spec_mean_da_cov" = "Average Fold Specific Data-Adaptive ATE Coverage",
  "v_spec_mean_gt_cov" = "Average Fold Specific Ground-Truth ATE Coverage",
  "pooled_da_cov" = "Harmonic Pooled Data-Adaptive Coverage",
  "pooled_gt_cov" =  "Harmonic Pooled Ground-Truth Coverage"
)

CVtreeMLE_cov_plot <- make_sim_statistics_plot(
  sim_statistics_long,
  stats = names(plot_labels),
  labels = plot_labels,
  color = "darkgreen",
  title = "Coverage Measures"
)

plot_labels <- c(
  "tmle_pooled_da_mse" = "Pooled TMLE MSE for Data-Adaptive Rule ATE",
  "tmle_pooled_gt_mse" = "Pooled TMLE MSE for Ground-Truth ATE",
  "v_spec_da_mean_mse" = "V-Specific TMLE MSE for Data-Adaptive Rule ATE",
  "v_spec_gt_mean_mse" = "V-Specific TMLE MSE for Ground-Truth Rule ATE",
  "v_pooled_da_mse" = "Harmonic Pooled MSE for Data-Adaptive Rule ATE",
  "v_pooled_gt_mse" = "Harmonic Pooled MSE for Ground-Truth Rule ATE"
)

CVtreeMLE_mse_plot <- make_sim_statistics_plot(
  sim_statistics_long,
  stats = names(plot_labels),
  labels = plot_labels,
  color = "darkgreen",
  title = "MSE Measures"
)

plot_labels <- c(
  "sqrt_n_abs_tmle_pooled_da_bias" = "Square Root Scaled Pooled TMLE Bias Compared to Data-Adaptive Rule ATE",
  "sqrt_n_abs_tmle_pooled_gt_bias" = "Square Root Scaled Pooled TMLE Bias Compared to Ground-Truth Rule ATE",
  "sqrt_n_abs_v_spec_da_mean_bias" = "Square Root Scaled Mean V-Specific TMLE Bias Compared to Data-Adaptive Rule ATE",
  "sqrt_n_abs_v_spec_gt_mean_bias" = "Square Root Scaled Mean V-Specific TMLE Bias Compared to Ground-Truth Rule ATE",
  "sqrt_n_abs_v_pooled_da_bias" = "Square Root Scaled Harmonic Mean Bias Compared to Data-Adaptive Rule ATE",
  "sqrt_n_abs_v_pooled_gt_bias" = "Square Root ScaledHarmonic Mean Bias Compared to Ground-Truth Rule ATE"
)

CVtreeMLE_bias_scaled_plot <- make_sim_statistics_plot(
  sim_statistics_long,
  stats = names(plot_labels),
  labels = plot_labels,
  color = "orange",
  title = "MSE Measures"
)

ggsave(
  here('sandbox/plots/CVtreeMLE_stats.png'),
  CVtreeMLE_stats_plot,
  width = 13,
  height = 7,
  device = png()
)

ggsave(
  here('sandbox/plots/CVtreeMLE_rule_stats.png'),
  CVtreeMLE_rule_stats_plot,
  width = 13,
  height = 7,
  device = png()
)

