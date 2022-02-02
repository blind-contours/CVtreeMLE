#' Create dot-whisker plots for the margina results found
#'
#' @param v_marginal_results Table of the marginal results for mixture components
#' @param mix_comps Character list of mixture variables
#'
#' @import ggplot2
#'
#' @export

plot_marginal_results <- function(v_marginal_results, mix_comps, hjust) {
  plot_list <- list()

  for (var in mix_comps) {
    marg_results <- v_marginal_results[stringr::str_detect(names(v_marginal_results), var)]
    if (length(marg_results) != 0) {
      marg_data <- do.call(rbind, marg_results)
      marg_data$`Marginal ATE` <- round(as.numeric(marg_data$`Marginal ATE`), 3)
      marg_data$`Lower CI` <- round(as.numeric(marg_data$`Lower CI`), 3)
      marg_data$`Upper CI` <- round(as.numeric(marg_data$`Upper CI`), 3)

      title <- paste("Marginal Results for ", var,".", " Pooled Ref is: ", unique(marg_data$Reference[marg_data$fold== "Pooled"]), sep = "")
      text_size <- 12
      line_size <- 2
      point_size <- 3
      plot_width <- 10
      plot_height <- 8
      text_theme <- ggplot2::element_text(size = text_size, color = "black")
      axis_text_theme <- ggplot2::element_text(size = text_size, color = "black")

      marg_data$Type <- factor(marg_data$comparison,
        levels = unique(marg_data$comparison)
        # labels = marg_data$Comparison[marg_data$fold == "Pooled"]
      )

      plot <- ggplot2::ggplot(
        marg_data,
        ggplot2::aes_string(
          x = "`Marginal ATE`", y = "fold", color = "fold",
          xmin = "`Lower CI`", xmax = "`Upper CI`", label="`Comparison`"
        )
      ) +
        ggplot2::facet_wrap(~Type) +
        ggplot2::geom_errorbarh(size = line_size) +
        ggplot2::geom_point(size = point_size) +
        ggplot2::geom_vline(xintercept = 0, alpha = .25, linetype = "dotted", size = line_size) +
        ggplot2::labs(x = "ATE", y = "Fold", color = "") +
        ggplot2::ggtitle(title) +
        ggplot2::theme_classic() +
        ggplot2::theme(text = text_theme, axis.text = axis_text_theme, legend.position = "none") +
        ggplot2::geom_text(size=5, hjust=hjust, vjust=0, colour = "#3C3C3C", nudge_x = 0, nudge_y = 0.0)

      # ggplot2::scale_x_continuous(breaks = round(seq(min(marg_data$`Marginal ATE`), max(marg_data$`Marginal ATE`), by =(max(marg_data$`Marginal ATE`) - min(marg_data$`Marginal ATE`)) / 5), 2))
      # ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45))

      plot_list[[var]] <- plot
    }
  }

  return(plot_list)
}
