#' Create dot-whisker plots for the marginal results found
#'
#' @param v_marginal_results Table of the marginal results for mixture
#' components
#' @param mix_comps Character list of mixture variables
#' @param hjust Degree to horizontally adjust the placement of rules on plot
#'
#' @export

plot_marginal_results <- function(v_marginal_results,
                                  mix_comps,
                                  hjust) {
  plot_list <- list()

  for (var in mix_comps) {
    marg_data <- v_marginal_results[v_marginal_results$var == var, ]
    if (length(marg_data) != 0) {
      marg_data$`Marginal ATE` <- round(as.numeric(marg_data$`Marginal ATE`), 3)
      marg_data$`Lower CI` <- round(as.numeric(marg_data$`Lower CI`), 3)
      marg_data$`Upper CI` <- round(as.numeric(marg_data$`Upper CI`), 3)

      title <- "Marginal Results by Fold"
      text_size <- 10
      line_size <- 2
      point_size <- 2
      text_theme <- ggplot2::element_text(size = text_size,
                                          color = "black")
      axis_text_theme <- ggplot2::element_text(size = text_size,
                                               color = "black")

      marg_data$Type <- factor(marg_data$Levels,
        levels = unique(marg_data$Levels)
      )

      plot <- ggplot2::ggplot(
        marg_data,
        ggplot2::aes_string(
          x = "`Marginal ATE`", y = "Type", color = "Type",
          xmin = "`Lower CI`",
          xmax = "`Upper CI`",
          label = "`Comparison_Rule`"
        )
      ) +
        ggplot2::facet_wrap(~fold) +
        ggplot2::geom_errorbarh(size = line_size) +
        ggplot2::geom_point(size = point_size) +
        ggplot2::geom_vline(xintercept = 0, alpha = .25, linetype = "dotted",
                            size = line_size) +
        ggplot2::labs(x = "ATE", y = "Threshold", color = "") +
        ggplot2::ggtitle(title) +
        ggplot2::theme_classic() +
        ggplot2::theme(text = text_theme, axis.text = axis_text_theme,
                       legend.position = "none") +
        ggplot2::geom_text(size = 3,
                           hjust = hjust,
                           vjust = 0,
                           colour = "#3C3C3C",
                           nudge_x = 0,
                           nudge_y = 0.0)


      plot_list[[var]] <- plot
    }
  }

  return(plot_list)
}
