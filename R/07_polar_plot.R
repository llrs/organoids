library("readxl")
library("dplyr")
library("tidyr")
library("ggplot2")
library("patchwork")

excels <- list.files("data/IPA/", full.names = TRUE)
# Information about the sheets present
sheets <- lapply(excels, excel_sheets)

# raw data ###
raw <- lapply(excels, read_xls, skip = 1, .name_repair = "universal")
names(raw) <- sub(".xls", "", basename(excels))
add0 <- function(x) {
  y <- x[1, , drop = FALSE]
  y$p.value.of.overlap <- max(x$p.value.of.overlap)
  rbind(y, x)
}



polar_graph <- function(i, raw, title = NULL) {
  title <- ifelse(is.null(title), names(raw)[i], title )
  data_plot <- raw[[i]] %>%
    mutate(p.value.of.overlap = sub(",", ".", p.value.of.overlap),
           Expr.Fold.Change = sub(",", ".", Expr.Fold.Change)) %>%
    mutate(p.value.of.overlap = as.numeric(p.value.of.overlap),
           Expr.Fold.Change = as.numeric(Expr.Fold.Change),
           group = title,
    ) %>%
    filter(!is.na(Predicted.Activation.State))


  if (NROW(data_plot) != 0)  {
    data_plot %>%
      group_by(Predicted.Activation.State) %>%
      arrange(p.value.of.overlap) %>%
      slice_head(n = 10) %>%
      ungroup() %>%
      arrange(p.value.of.overlap) %>%
      mutate(Upstream.Regulator = forcats::fct_reorder(Upstream.Regulator,
                                                       p.value.of.overlap)) %>%
      ggplot(aes(Upstream.Regulator, -log10(p.value.of.overlap))) +
      geom_polygon(aes(group = group), fill = "transparent", col = "black") +
      geom_point(aes(col = Predicted.Activation.State,
                     shape = Predicted.Activation.State), size = 3) +
      coord_polar(start = - pi / 20, clip = "off")  + # From https://r-graphics.org/recipe-axes-polar#discussion-74
      ylim(0, NA) +
      theme_minimal() +
      labs(y = "-log10(p.value)", x = "Upstream regulator", title = title,
           col = "State", shape = "State") +
      theme(panel.grid = element_line(colour = "gray"),
            axis.text.x = element_text(hjust = 0, vjust = 0, face = "bold"),
            axis.title.y = element_text(angle = 0))

  } else {
    NULL

  }
}

l <- lapply(seq_along(raw), polar_graph, raw)
k <- sapply(l, is.null)
l2 <- l[!k]

l2[[1]] + l2[[2]] +
l2[[3]] + l2[[4]] +
l2[[5]] + l2[[6]] +
l2[[7]]  + l2[[8]] +
l2[[9]] + l2[[10]] +
l2[[11]] + l2[[12]] +
l2[[13]] + l2[[14]] + plot_layout(byrow = TRUE, guides='collect')
