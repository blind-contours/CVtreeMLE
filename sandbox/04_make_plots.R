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


sim_statistics <- sim_results_1 %>%
  group_by(n_obs) %>%
  summarise(
    est_bias = mean(Bias),
    est_sd = sd(ATE),
    est_MSE = est_bias^2 + est_sd^2,
    DA_CI_coverage = mean(`DA Cov`),
    Truth_CI_coverage = mean(`True Cov`),
    Mix_found = mean(`Mix ind`),
    True_pos = mean(`True Pos`),
    True_neg = mean(`True Neg`),
    False_pos = mean(`False Pos`),
    False_neg = mean(`False Neg`),

  ) %>%
  mutate(
    abs_bias = abs(est_bias),
    sqrt_n_abs_bias = sqrt(n_obs) * abs(est_bias),
    n_MSE = n_obs*est_MSE,
  ) %>%
  ungroup

sim_statistics_long <- sim_statistics %>%
  tidyr::gather(statistic, value, -c(n_obs))


make_sim_statistics_plot <- function(sim_statistics_long, stats) {
  filtered_sim_stats <- sim_statistics_long %>%
    filter(
      statistic %in% stats
    )

  ggplot(filtered_sim_stats, aes(n_obs, value)) +
    geom_point() + geom_line(linetype = "dashed") +
    facet_wrap("statistic", scales = "free_y") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle=45, vjust=1, hjust=1)
    ) +
    scale_x_sqrt(breaks=n_obs)
}


CVtreeMLE_stats_plot <- make_sim_statistics_plot(
  sim_statistics_long,
  stats = c("abs_bias", "est_sd",
            "sqrt_n_abs_bias", "est_MSE", "DA_CI_coverage", "n_MSE")
)

CVtreeMLE_rule_stats_plot <- make_sim_statistics_plot(
  sim_statistics_long,
  stats = c("Mix_found", "True_pos", "True_neg", "False_pos","False_neg")
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

