#' Minimal publication-figure ggplot2 theme
#'
#' No title/subtitle/caption chrome, no gridlines, no gratuitous background --
#' just axis titles, tick labels, and the legend. Intended so plot functions
#' never bake in a title (callers add one via `patchwork`/`cowplot`/knitr
#' chunk captions if they want one at all).
#'
#' @param base_size Base font size.
#' @return A ggplot2 theme.
#' @export
mp_theme_pub <- function(base_size = 11) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank(),
      legend.key.size = ggplot2::unit(0.8, "lines"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold")
    )
}
