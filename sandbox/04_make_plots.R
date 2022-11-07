library(tidyverse)
library(here)
library(ggpubr)
library(purrr)
library(tidyr)
library(hrbrthemes)

n_obs <- c(200, 350, 500, 750, 1000, 1500, 2000, 3000, 5000)

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

sim_results_5 <- readRDS(
  here("sandbox/data/CVtreeMLE_run_5.rds")
)

sim_results <- rbind(sim_results_2, sim_results_3)


sim_statistics <- sim_results %>%
  group_by(n_obs) %>%
  summarise(
    abs_tmle_pooled_da_bias = abs(mean(tmle_pooled_da_bias)),
    abs_tmle_pooled_gt_bias = abs(mean(tmle_pooled_gt_bias)),
    abs_v_spec_da_mean_bias = abs(mean(v_spec_da_mean_bias)),
    abs_v_spec_gt_mean_bias = abs(mean(v_spec_gt_mean_bias)),
    abs_v_pooled_da_bias = abs(mean(v_pooled_da_bias)),
    abs_v_pooled_gt_bias = abs(mean(v_pooled_gt_bias)),
    tmle_pooled_sd = sd(tmle_pooled_ate),
    v_spec_mean_sd = sd(v_spec_mean_ate),
    v_pooled_mean_sd = sd(v_pooled_ate),
    tmle_pooled_da_mse = abs_tmle_pooled_da_bias^2 + tmle_pooled_sd^2,
    tmle_pooled_gt_mse = abs_tmle_pooled_gt_bias^2 + tmle_pooled_sd^2,
    v_spec_da_mean_mse = abs_v_spec_da_mean_bias^2 + v_spec_mean_sd^2,
    v_spec_gt_mean_mse = abs_v_spec_gt_mean_bias^2 + v_spec_mean_sd^2,
    v_pooled_da_mse = abs_v_pooled_da_bias^2 + v_pooled_mean_sd^2,
    v_pooled_gt_mse = abs_v_pooled_gt_bias^2 + v_pooled_mean_sd^2,
    tmle_pooled_da_z_bias = mean(tmle_pooled_da_bias)/tmle_pooled_sd,
    tmle_pooled_gt_z_bias = mean(tmle_pooled_gt_bias)/tmle_pooled_sd,
    v_spec_da_z_bias = mean(v_spec_da_mean_bias)/v_spec_mean_sd,
    v_spec_gt_z_bias = mean(v_spec_gt_mean_bias)/v_spec_mean_sd,
    v_pooled_da_z_bias = mean(v_pooled_da_bias) /v_pooled_mean_sd,
    v_pooled_gt_z_bias = mean(v_pooled_gt_bias) /v_pooled_mean_sd,
    tmle_pooled_gt_coverage = mean(tmle_pooled_gt_coverage),
    tmle_pooled_da_coverage = mean(tmle_pooled_da_coverage),
    v_spec_mean_da_cov = mean(v_spec_mean_da_cov),
    v_spec_mean_gt_cov = mean(v_spec_mean_gt_cov),
    v_pooled_da_cov = mean(pooled_da_cov),
    v_pooled_gt_cov = mean(pooled_gt_cov),
    True_pos = mean(`True Pos`),
    True_neg = mean(`True Neg`),
    False_pos = mean(`False Pos`),
    False_neg = mean(`False Neg`),
    tmle_pooled_sd = sd(tmle_pooled_ate),
    v_spec_mean_sd = sd(v_spec_mean_ate),
    v_pooled_mean_sd = sd(v_pooled_ate),
    tmle_pooled_da_z_bias = tmle_pooled_da_bias/tmle_pooled_sd,
    tmle_pooled_gt_z_bias = tmle_pooled_gt_bias/tmle_pooled_sd,
    v_spec_da_z_bias = v_spec_da_mean_bias/v_spec_mean_sd,
    v_spec_gt_z_bias = v_spec_gt_mean_bias/v_spec_mean_sd,
    v_pooled_da_z_bias = v_pooled_da_bias /v_pooled_mean_sd,
    v_pooled_gt_z_bias = v_pooled_gt_bias /v_pooled_mean_sd,
    sqrt_n_abs_tmle_pooled_da_bias = sqrt(n_obs) * abs_tmle_pooled_da_bias,
    sqrt_n_abs_tmle_pooled_gt_bias = sqrt(n_obs) * abs_tmle_pooled_gt_bias,
    sqrt_n_abs_v_spec_da_mean_bias = sqrt(n_obs) * abs_v_spec_da_mean_bias,
    sqrt_n_abs_v_spec_gt_mean_bias = sqrt(n_obs) * abs_v_spec_gt_mean_bias,
    sqrt_n_abs_v_pooled_da_bias = sqrt(n_obs) * abs_v_pooled_da_bias,
    sqrt_n_abs_v_pooled_gt_bias = sqrt(n_obs) * abs_v_pooled_gt_bias,
    n_tmle_pooled_da_mse = n_obs*tmle_pooled_da_mse,
    n_tmle_pooled_gt_mse = n_obs*tmle_pooled_gt_mse,
    n_v_spec_da_mean_mse = n_obs*v_spec_da_mean_mse,
    n_v_spec_gt_mean_mse = n_obs*v_spec_gt_mean_mse,
    n_v_pooled_da_mse = n_obs*v_pooled_da_mse,
    n_v_spec_mean_gt_cov = n_obs*v_spec_mean_gt_cov)


sim_statistics_long <- sim_statistics %>%
  tidyr::gather(statistic, value, -c(n_obs))


make_sim_statistics_plot <- function(sim_statistics_long,
                                     stats,
                                     labels,
                                     title) {
  filtered_sim_stats <- sim_statistics_long %>%
    filter(
      statistic %in% stats
    )

  ggplot(filtered_sim_stats, aes(x=n_obs,
                                 y=value,
                                 group=statistic,
                                 color=statistic)) +
    geom_line(position = position_jitter(w=0.02, h=0.01)) +
    scale_color_viridis(discrete = TRUE) +
    ggtitle("Popularity of American names in the previous 30 years") +
    theme_ipsum() +
    ylab("Bias") +
    xlab("Number Observations") +
    ggtitle(title)

}

make_density_plot <- function(sim_statistics,
                                     z_stat,
                                     label) {


  sim_statistics$n_obs <- as.factor(sim_statistics$n_obs)

  ggplot(sim_statistics, aes(x=eval(parse(text = z_stat)), color=n_obs, fill=n_obs)) +
    geom_density(alpha=0.4) +
    scale_fill_viridis(discrete=TRUE) +
    scale_color_viridis(discrete=TRUE) +
    theme_ipsum() +
    theme(
      legend.position="left",
      panel.spacing = unit(0.1, "lines"),
      strip.text.x = element_text(size = 8)
    ) +
    ggtitle(label) +
    xlab("Z-Score Standardized Bias (Bias/SD)") +
    ylab("Probability")

}


plot_labels <- c(
  "abs_tmle_pooled_da_bias" = "Pooled TMLE (DA)",
  "abs_v_spec_da_mean_bias" = "Mean Fold TMLE (DA)",
  "abs_v_pooled_da_bias" = "Harmonic Pooled (DA)"
)

CVtreeMLE_bias_plot <- make_sim_statistics_plot(
  sim_statistics_long,
  stats = names(plot_labels),
  labels = plot_labels,
  title = "DA Bias Measures"
)

plot_labels <- c("abs_tmle_pooled_gt_bias" = "Pooled TMLE (GT)",
                 "abs_v_spec_gt_mean_bias" = "Mean Fold TMLE (GT)",
                 "abs_v_pooled_gt_bias" = "Harmonic Pooled (GT)")

CVtreeMLE_bias_plot <- make_sim_statistics_plot(
  sim_statistics_long,
  stats = names(plot_labels),
  labels = plot_labels,
  title = "GT Bias Measures"
)

plot_labels <- c(
  "True_pos" = "True Positive",
  "True_neg" = "True Negative",
  "False_pos" = "False Positive",
  "False_neg" = "False Negative"
)

CVtreeMLE_rule_plot <- make_sim_statistics_plot(
  sim_statistics_long,
  stats = names(plot_labels),
  labels = plot_labels,
  title = "Rule Coverage Measures"
)


plot_labels <- c(
  "tmle_pooled_da_coverage" = "Pooled TMLE (DA)",
  "v_pooled_da_cov" = "Harmonic Pooled (DA)",
  "v_spec_mean_da_cov" = "Mean Fold (DA)"
)

CVtreeMLE_cov_plot <- make_sim_statistics_plot(
  sim_statistics_long,
  stats = names(plot_labels),
  labels = plot_labels,
  title = "DA Coverage Measures"
)

plot_labels <- c("tmle_pooled_gt_coverage" = "Pooled TMLE (GT)",
                 "v_spec_mean_gt_cov" = "Mean Fold (GT)",
                 "v_pooled_gt_cov" =  "Harmonic Pooled (GT)")

CVtreeMLE_cov_plot <- make_sim_statistics_plot(
  sim_statistics_long,
  stats = names(plot_labels),
  labels = plot_labels,
  title = "GT Coverage Measures"
)


plot_labels <- c(
  "tmle_pooled_da_mse" = "Pooled TMLE (DA)",
  "v_spec_da_mean_mse" = "Mean Fold (DA)",
  "v_pooled_da_mse" = "Harmonic Pooled (DA)"
)

CVtreeMLE_mse_plot <- make_sim_statistics_plot(
  sim_statistics_long,
  stats = names(plot_labels),
  labels = plot_labels,
  title = "DA MSE Measures"
)

plot_labels <- c(
  "tmle_pooled_gt_mse" = "Pooled TMLE (GT)",
  "v_spec_gt_mean_mse" = "Mean Fold (GT)",
  "v_pooled_gt_mse" = "Harmonic Pooled (GT)"
)

CVtreeMLE_mse_plot <- make_sim_statistics_plot(
  sim_statistics_long,
  stats = names(plot_labels),
  labels = plot_labels,
  title = "GT MSE Measures"
)

plot_labels <- c("sqrt_n_abs_tmle_pooled_gt_bias" = "Pooled TMLE (GT)",
  "sqrt_n_abs_v_spec_gt_mean_bias" = "Mean Fold (GT)",
  "sqrt_n_abs_v_pooled_gt_bias" = "Harmonic Pooled (GT)")

plot_labels <- c(
  "sqrt_n_abs_tmle_pooled_da_bias" = "Pooled TMLE (DA)",
  "sqrt_n_abs_v_spec_da_mean_bias" = "Mean Fold (DA)",
  "sqrt_n_abs_v_pooled_da_bias" = "Harmonic Pooled (DA)"
)

CVtreeMLE_bias_scaled_plot <- make_sim_statistics_plot(
  sim_statistics_long,
  stats = names(plot_labels),
  labels = plot_labels,
  color = "#70DF00",
  title = "Square Root Scaled Bias"
)

plot_labels <- c(
  "tmle_pooled_da_z_bias" = "Pooled TMLE (DA)",
  "tmle_pooled_gt_z_bias" = "Pooled TMLE (GT)",
  # "v_spec_da_z_bias" = "Mean Fold (DA)",
  # "v_spec_gt_z_bias" = "Mean Fold (GT)",
  # "v_pooled_da_z_bias" = "Harmonic Pooled (DA)",
  # "v_pooled_gt_z_bias" = "Harmonic Pooled (GT)",
  "Normal" = "Normal Mean 0, SD 1"
)


make_density_plot(sim_statistics, z_stat = "tmle_pooled_da_z_bias",
                  label = "Pooled TMLE (DA)")


make_density_plot(sim_statistics_long_w_norm, stats = c("v_spec_da_z_bias",
                                                        "v_spec_gt_z_bias",
                                                        "Normal"),
                  label = c("Normal Mean 0, SD 1","Mean Fold TMLE (DA)", "Mean Fold TMLE (GT)" ))

make_density_plot(sim_statistics_long_w_norm, stats = c("v_pooled_da_z_bias",
                                                        "v_pooled_gt_z_bias",
                                                        "Normal"),
                  label = c("Normal Mean 0, SD 1","Harmonic Pooled TMLE (DA)", "Harmonic Pooled TMLE (GT)" ))


ggsave(
  here('sandbox/plots/CVtreeMLE_bias.png'),
  CVtreeMLE_bias_plot,
  width = 15,
  height = 6,
  device = png()
)
dev.off()

ggsave(
  here('sandbox/plots/CVtreeMLE_rule_stats.png'),
  CVtreeMLE_rule_plot,
  width = 15,
  height = 6,
  device = png()
)

ggsave(
  here('sandbox/plots/CVtreeMLE_cov.png'),
  CVtreeMLE_cov_plot,
  width = 15,
  height = 6,
  device = png()
)

ggsave(
  here('sandbox/plots/CVtreeMLE_mse.png'),
  CVtreeMLE_mse_plot,
  width = 15,
  height = 6,
  device = png()
)

ggsave(
  here('sandbox/plots/CVtreeMLE_bias_scaled.png'),
  CVtreeMLE_bias_scaled_plot,
  width = 15,
  height = 6,
  device = png()
)




