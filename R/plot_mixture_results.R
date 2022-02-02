#' Create dot-whisker plots for the mixture results found
#'
#' @param v_intxn_results Table of the interaction results for a mixture
#'
#' @import ggplot2
#'
#' @export

plot_mixture_results <- function(v_intxn_results, move_right) {
  plot_list <- list()
  for (i in seq(v_intxn_results)) {
    intxn_results <- v_intxn_results[[i]]
    intxn_results$`Mixture ATE` <- round(as.numeric(intxn_results$`Mixture ATE`), 3)
    intxn_results$`Lower CI` <- round(as.numeric(intxn_results$`Lower CI`), 3)
    intxn_results$`Upper CI` <- round(as.numeric(intxn_results$`Upper CI`), 3)

    intxn_results$fold <- rownames(intxn_results)
    title <- paste("Interaction for:", unique(intxn_results$Variables))
    text_size <- 12
    line_size <- 2
    point_size <- 3
    plot_width <- 20
    plot_height <- 6
    text_theme <- ggplot2::element_text(size = text_size, color = "black")
    axis_text_theme <- ggplot2::element_text(size = text_size, color = "black")

    plot <- ggplot2::ggplot(
      intxn_results,
      ggplot2::aes_string(
        x = "`Mixture ATE`", y = "fold",
        xmin = "`Lower CI`", xmax = "`Upper CI`",
        color = "fold", label="`Mixture Interaction Rules`"
      )
    ) +
      ggplot2::geom_errorbarh(size = line_size) +
      ggplot2::geom_point(size = point_size) +
      ggplot2::geom_vline(xintercept = 0, alpha = .25, linetype = "dotted", size = line_size) +
      ggplot2::labs(x = "ATE", y = "Fold", color = "") +
      ggplot2::ggtitle(title) +
      ggplot2::theme_classic() +
      ggplot2::theme(text = text_theme, axis.text = axis_text_theme, legend.position = "none") +
      ggplot2::geom_text(size=4, hjust=move_right, vjust=0, colour = "#3C3C3C", nudge_x = -0, nudge_y = 0.0)

    # ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45))

    plot_list[[i]] <- plot
  }

  names(plot_list) <- names(v_intxn_results)

  return(plot_list)
}
