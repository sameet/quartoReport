#' Convert a ggplot object to a plotly object
#'
#' This function takes a ggplot object and converts it into an interactive
#' plotly object using plotly::ggplotly().
#'
#' @param gg_plot A ggplot object.
#' @return A plotly object.
#' @export
#' @examples
#' \dontrun{
#' library(ggplot2)
#' p_static <- ggplot(data = mtcars, aes(x = mpg, y = wt)) + geom_point()
#' p_interactive <- to_interactive_plot(p_static)
#' }
to_interactive_plot <- function(gg_plot) {
  return(plotly::ggplotly(gg_plot))
}
