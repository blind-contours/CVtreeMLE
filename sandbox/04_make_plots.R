library(tidyverse)
library(here)
library(ggpubr)
library(purrr)
library(tidry)

psi_true <- 18

sim_results_1 <- readRDS(
  here("sandbox/data/CVtreeMLE_run_1.rds")
)

sim_results_2 <- readRDS(
  here("sandbox/data/CVtreeMLE_run_2.rds")
)

sim_results_3 <- readRDS(
  here("sandbox/data/CVtreeMLE_run_3.rds")
)


sim_statistics <- CVtreeMLE_run_1 %>%
  group_by(n_obs) %>%
  summarise(
    est_bias = mean(Bias),
    est_da_bias = mean(`DA Rule Bias`),
    est_sd = sd(ATE),
    est_true_MSE = est_bias^2 + est_sd^2,
    est_da_MSE = est_da_bias^2 + est_sd^2,
    DA_CI_coverage = mean(`DA Cov`),
    Truth_CI_coverage = mean(`True Cov`),
    Mix_found = mean(`Mix ind`),
    True_pos = mean(`True Pos`),
    True_neg = mean(`True Neg`),
    False_pos = mean(`False Pos`),
    False_neg = mean(`False Neg`),

  ) %>%
  mutate(
    abs_true_bias = abs(est_bias),
    abs_da_bias = abs(est_da_bias),
    sqrt_n_abs_true_bias = sqrt(n_obs) * abs(est_bias),
    sqrt_n_abs_da_bias = sqrt(n_obs) * abs(est_da_bias),
    n_true_MSE = n_obs*est_true_MSE,
    n_da_MSE = n_obs*est_da_MSE
  ) %>%
  ungroup

sim_statistics_long <- sim_statistics %>%
  tidyr::gather(statistic, value, -c(n_obs))


make_sim_statistics_plot <- function(sim_statistics_long, stats, labels,) {
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
    scale_x_sqrt(breaks=n_obs)  + labs(x = "Number Observations", title = "Bias Measures")

}

plot_labels <- c(
  "abs_true_bias" = "Absolute Bias Compared to Truth",
  "abs_da_bias" = "Absolute Bias Compared to Data-Adaptive Truth",
  "sqrt_n_abs_true_bias" = "Absolute Bias Compared to Truth Scaled by Root N",
  "sqrt_n_abs_da_bias" = "Absolute Bias Compared to Data-Adaptive Truth Scaled by Root N"
)

CVtreeMLE_stats_plot <- make_sim_statistics_plot(
  sim_statistics_long,
  stats = c("abs_true_bias", "abs_da_bias",
            "sqrt_n_abs_true_bias", "sqrt_n_abs_da_bias"),
  labels = plot_labels
)

plot_labels <- c(
  "Mix_found" = "Proportion Mixture Rule Found",
  "True_pos" = "Rule Indicator True Positive Compared to Ground-Truth",
  "True_neg" = "Rule Indicator True Negative Compared to Ground-Truth",
  "False_pos" = "Rule Indicator False Positive Compared to Ground-Truth",
  "False_neg" = "Rule Indicator False Negative Compared to Ground-Truth"
)

CVtreeMLE_rule_stats_plot <- make_sim_statistics_plot(
  sim_statistics_long,
  stats = c("Mix_found", "True_pos", "True_neg", "False_pos","False_neg"),
  color = "orange"
)

CVtreeMLE_rule_stats_plot <- make_sim_statistics_plot(
  sim_statistics_long,
  stats = c("est_sd", "est_true_MSE", "est_da_MSE", "DA_CI_coverage","Truth_CI_coverage")
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

