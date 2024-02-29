library(tidyverse)
library(here)
library(ggpubr)
library(purrr)
library(tidyr)
library(hrbrthemes)
library(viridis)

n_obs <- c(300, 500, 750, 1000, 1500, 2000, 5000)

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

sim_results_6 <- readRDS(
  here("sandbox/data/CVtreeMLE_run_6.rds")
)

sim_results_7 <- readRDS(
  here("sandbox/data/CVtreeMLE_run_7.rds")
)

sim_results_8 <- readRDS(
  here("sandbox/data/CVtreeMLE_run_8.rds")
)

sim_results_9 <- readRDS(
  here("sandbox/data/CVtreeMLE_run_9.rds")
)

sim_results_10 <- readRDS(
  here("sandbox/data/CVtreeMLE_run_10.rds")
)

sim_results_11 <- readRDS(
  here("sandbox/data/CVtreeMLE_run_11.rds")
)

sim_results <- rbind(sim_results_1, sim_results_2)

sim_statistics <- sim_results %>%
  group_by(n_obs) %>%
  summarise(
    abs_tmle_pooled_da_bias = abs(mean(tmle_pooled_da_bias)),
    abs_tmle_pooled_gt_bias = abs(mean(tmle_pooled_gt_bias)),
    abs_v_spec_da_mean_bias = abs(mean(v_spec_da_mean_bias)),
    abs_v_spec_gt_mean_bias = abs(mean(v_spec_gt_mean_bias)),
    pooled_tmle_CI_range = mean(pooled_tmle_CI_range),
    v_spec_tmle_CI_range = mean(v_spec_tmle_CI_range),
    abs_tree_mean_bias = abs(mean(rule_min_bias)),
    tmle_pooled_sd = sd(tmle_pooled_ate),
    v_spec_mean_sd = sd(v_spec_mean_ate),
    tmle_pooled_da_mse = abs_tmle_pooled_da_bias^2 + tmle_pooled_sd^2,
    tmle_pooled_gt_mse = abs_tmle_pooled_gt_bias^2 + tmle_pooled_sd^2,
    v_spec_da_mean_mse = abs_v_spec_da_mean_bias^2 + v_spec_mean_sd^2,
    v_spec_gt_mean_mse = abs_v_spec_gt_mean_bias^2 + v_spec_mean_sd^2,
    tmle_pooled_da_z_bias = mean(tmle_pooled_da_bias)/tmle_pooled_sd,
    tmle_pooled_gt_z_bias = mean(tmle_pooled_gt_bias)/tmle_pooled_sd,
    v_spec_da_z_bias = mean(v_spec_da_mean_bias)/v_spec_mean_sd,
    v_spec_gt_z_bias = mean(v_spec_gt_mean_bias)/v_spec_mean_sd,
    tmle_pooled_gt_coverage = mean(tmle_pooled_gt_coverage),
    tmle_pooled_da_coverage = mean(tmle_pooled_da_coverage),
    v_spec_mean_da_cov = mean(v_spec_mean_da_cov),
    v_spec_mean_gt_cov = mean(v_spec_mean_gt_cov),
    True_pos = mean(`True Pos`),
    True_neg = mean(`True Neg`),
    False_pos = mean(`False Pos`),
    False_neg = mean(`False Neg`),
    tmle_pooled_sd = sd(tmle_pooled_ate),
    v_spec_mean_sd = sd(v_spec_mean_ate)
    )


sim_statistics <- sim_statistics[sim_statistics$n_obs < 10000 ,]

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
    geom_smooth(method = "loess", se = FALSE, position = position_jitter(w=0.0, h=0.0)) +
    scale_color_viridis(discrete = TRUE) +
    ggtitle("Popularity of American names in the previous 30 years") +
    theme_ipsum() +
    ylab("Bias") +
    xlab("Number Observations") +
    ggtitle(title)

}

make_bias_plot <- function(sim_statistics_long, stats, title, legend_labels) {

  filtered_sim_stats <- sim_statistics_long %>%
    filter(
      statistic %in% stats
    )

  filtered_sim_stats <- filtered_sim_stats %>%
    group_by(statistic) %>%
    mutate(
      theoretical_bias = ifelse(n_obs == min(n_obs), value, NA)
    ) %>%
    fill(theoretical_bias) %>%
    ungroup() %>%
    mutate(
      theoretical_bias = theoretical_bias / sqrt(n_obs / min(n_obs))
    )


  filtered_sim_stats %>%
    ggplot(aes(x = n_obs, y = value, group = statistic, color = statistic)) +
    geom_line() +
    geom_line(aes(y = theoretical_bias), linetype = "dashed") +
    scale_color_manual(values = viridis(3, option = "plasma"), labels = legend_labels) +
    labs(
      title = title,
      subtitle = "Dashed lines represent the expected bias decrease if the estimators were sqrt(n)-consistent.",
      x = "Number of Observations",
      y = "Absolute MSE",
      color = "Effect"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom",
      strip.text = element_text(face = "bold")
    )
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
  "abs_tree_mean_bias" = "Average Minimum Region Bias From Decision Tree"
)

CVtreeMLE_DA_bias_plot <- make_bias_plot(
  sim_statistics_long,
  stats = names(plot_labels),
  title = "Decision Tree Bias",
  legend_labels = plot_labels
)


plot_labels <- c(
  "abs_tree_mean_bias" = "Average Minimum Region Bias From Decision Tree"
)

CVtreeMLE_DA_bias_plot <- make_bias_plot(
  sim_statistics_long,
  stats = names(plot_labels),
  title = "Decision Tree Bias",
  legend_labels = plot_labels
)



plot_labels <- c(
  "tmle_pooled_da_mse" = "Pooled TMLE (DA)",
  "v_spec_da_mean_mse" = "Mean Fold TMLE (DA)"
)

CVtreeMLE_DA_bias_plot <- make_bias_plot(
  sim_statistics_long,
  stats = names(plot_labels),
  title = "DA MSE Measures",
  legend_labels = plot_labels
)

plot_labels <- c("abs_tmle_pooled_gt_bias" = "Pooled TMLE (GT)",
                 "abs_v_spec_gt_mean_bias" = "Mean Fold TMLE (GT)")

CVtreeMLE_GT_bias_plot <- make_bias_plot(
  sim_statistics_long,
  stats = names(plot_labels),
  legend_labels = plot_labels,
  title = "GT MSE Measures"
)

plot_labels <- c("tmle_pooled_gt_mse" = "Pooled TMLE (GT)",
                 "v_spec_gt_mean_mse" = "Mean Fold TMLE (GT)")

CVtreeMLE_GT_bias_plot <- make_bias_plot(
  sim_statistics_long,
  stats = names(plot_labels),
  legend_labels = plot_labels,
  title = "GT MSE Measures"
)


plot_labels <- c("tmle_pooled_gt_mse" = "Pooled TMLE (DA)",
                 "v_spec_da_mean_mse" = "Mean Fold TMLE (DA)")

CVtreeMLE_GT_bias_plot <- make_bias_plot(
  sim_statistics_long,
  stats = names(plot_labels),
  legend_labels = plot_labels,
  title = "DA MSE Measures"
)




plot_labels <- c(
  "True_neg" = "True Negative",
  "False_pos" = "False Positive",
  "True_pos" = "True Positive",
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
  "v_spec_mean_da_cov" = "Mean Fold (DA)"
)

CVtreeMLE_DA_cov_plot <- make_sim_statistics_plot(
  sim_statistics_long,
  stats = names(plot_labels),
  labels = plot_labels,
  title = "DA Coverage Measures"
)

plot_labels <- c("tmle_pooled_gt_coverage" = "Pooled TMLE (GT)",
                 "v_spec_mean_gt_cov" = "Mean Fold (GT)")

CVtreeMLE_GT_cov_plot <- make_sim_statistics_plot(
  sim_statistics_long,
  stats = names(plot_labels),
  labels = plot_labels,
  title = "GT Coverage Measures"
)


plot_labels <- c(
  "tmle_pooled_da_mse" = "Pooled TMLE (DA)",
  "v_spec_da_mean_mse" = "Mean Fold (DA)")

CVtreeMLE_DA_mse_plot <- make_sim_statistics_plot(
  sim_statistics_long,
  stats = names(plot_labels),
  labels = plot_labels,
  title = "DA MSE Measures"
)

plot_labels <- c(
  "tmle_pooled_gt_mse" = "Pooled TMLE (GT)",
  "v_spec_gt_mean_mse" = "Mean Fold (GT)")

CVtreeMLE_GT_mse_plot <- make_sim_statistics_plot(
  sim_statistics_long,
  stats = names(plot_labels),
  labels = plot_labels,
  title = "GT MSE Measures"
)

plot_labels <- c("sqrt_n_abs_tmle_pooled_gt_bias" = "Pooled TMLE (GT)",
  "sqrt_n_abs_v_spec_gt_mean_bias" = "Mean Fold (GT)")

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


make_density_plot(sim_statistics, z_stat = "v_spec_da_z_bias",
                  label = c("Mean Fold TMLE (DA)"))

make_density_plot(sim_statistics, z_stat = "v_pooled_da_z_bias",
                  label = c("Harmonic Pooled (DA)"))



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




