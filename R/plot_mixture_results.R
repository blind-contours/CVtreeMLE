#' Create dot-whisker plots for the mixture results found
#'
#' @param v_intxn_results Table of the interaction results for a mixture
#' @param hjust Horizontal adjustment of rule placement on plots relative to
#' point estimate
#'
#' @import viridis
#' @import hrbrthemes
#'
#' @export

plot_mixture_results <- function(v_intxn_results, hjust) {
  plot_list <- list()
  for (i in seq(v_intxn_results)) {
    intxn_results <- v_intxn_results[[i]]
    intxn_results$ate <- round(as.numeric(intxn_results$ate), 3)
    intxn_results$lower_ci <- round(as.numeric(intxn_results$lower_ci), 3)
    intxn_results$upper_ci <- round(as.numeric(intxn_results$upper_ci), 3)

    title <- paste("Interaction for:", unique(intxn_results$variables))
    text_size <- 12
    line_size <- 2
    point_size <- 3
    text_theme <- ggplot2::element_text(size = text_size, color = "black")
    axis_text_theme <- ggplot2::element_text(size = text_size, color = "black")

    plot <- ggplot2::ggplot(
      intxn_results,
      ggplot2::aes_string(
        x = "ate",
        y = "fold",
        xmin = "lower_ci",
        xmax = "upper_ci",
        color = "fold",
        label = "mix_rule"
      )
    ) +
      ggplot2::geom_errorbarh(size = line_size) +
      ggplot2::geom_point(size = point_size) +
      ggplot2::geom_vline(xintercept = 0, alpha = .25, linetype = "dotted",
                          size = line_size) +
      ggplot2::labs(x = "ATE", y = "Fold", color = "") +
      ggplot2::ggtitle(title) +
      ggplot2::theme_classic() +
      ggplot2::theme(text = text_theme, axis.text = axis_text_theme,
                     legend.position = "none") +
      ggplot2::geom_text(size = 4, hjust = hjust, vjust = 0,
                         colour = "#3C3C3C", nudge_x = -0, nudge_y = 0.0) +
      viridis::scale_color_viridis(discrete = TRUE) +
      hrbrthemes::theme_ipsum()

    plot_list[[i]] <- plot
  }

  names(plot_list) <- names(v_intxn_results)

  return(plot_list)
}
