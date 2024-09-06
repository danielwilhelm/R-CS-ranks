#' @describeIn plotranking Plot \code{csranks} output
#' 
#' @param x An \code{csranks} object, likely produced by \code{\link{csranks}}.
#' @param ... Other arguments, passed to \code{plotranking}.
#' @export

plot.csranks <- function(x, ...){
  plotranking(x$rank, x$L, x$U, ...)
}

#' Plot ranking with confidence sets
#' 
#' Display ranks together with their confidence set bounds.
#'
#' @param ranks vector of ranks
#' @param L vector of lower bounds of confidence sets for the ranks
#' @param U vector of lower bounds of confidence sets for the ranks
#' @param popnames vector containing names of the populations whose ranks are in \code{ranks}. If \code{popnames=NULL} (default), then populations are automatically numbered.
#' @param title character string containing the main title of the graph. \code{title=NULL} (default) means no title.
#' @param subtitle character string containing the subtitle of the graph. \code{subtitle=NULL} (default) means no subtitle.
#' @param caption character string containing the caption of the graph. \code{caption=NULL} (default) means no caption.
#' @param colorbins integer indicating the number of quantile bins into which populations are grouped and color-coded. Value has to lie between 1 (default) and the number of populations.
#' @param horizontal logical. Should be the bars displayed horizontally, or vertically?
#' @return A ggplot plot displaying confidence sets.

#' @examples
#' x <- seq(1, 3, length = 10)
#' V <- diag(rep(0.04, 10))
#' CS <- csranks(x, V)
#' grid::current.viewport()
#' plot(CS)
#' # Equivalent: 
#' plotranking(CS$rank, CS$L, CS$U)
#' 
#' # plotranking returns a ggplot object. It can be customized further:
#' library(ggplot2)
#' pl <- plot(CS)
#' pl + xlab("position in ranking") + ylab("population label") + theme_gray()
#' 
#' # horizontal = FALSE uses ggplot2::coord_flip underneath. The x and y axes swap places.
#' pl <- plot(CS, horizontal = FALSE)
#' pl + xlab("position in ranking") + # Note, that xlab refers to vertical axis now
#'   ylab("population label") + theme_gray()
#' @export
#' @import ggplot2
#' @importFrom stats reorder
plotranking <- function(ranks, L, U, popnames = NULL, title = NULL, subtitle = NULL,
                        caption = NULL, colorbins = 1, horizontal = TRUE) {
  # initializations
  check_plotranking_args(ranks, L, U, popnames, title, subtitle,
                         caption, colorbins, horizontal)
  cc <- scales::seq_gradient_pal("#66a182", "#d1495b", "Lab")(seq(0, 1, length.out = 4))
  p <- length(ranks)
  max_rank <- max(U)
  min_rank <- min(L)
  if (is.null(popnames)) popnames <- 1:p
  dat <- data.frame(ranks = ranks, L = L, U = U, popnames = popnames)
  coltitle <- ifelse(colorbins == 4, "quartile", "quantile bins")

  # plot
  if (colorbins > 1) {
    pl <- ggplot(data = dat, aes(
      y = reorder(popnames, ranks),
      x = ranks, color = createbins(ranks, colorbins)
    ))
  } else {
    pl <- ggplot(data = dat, aes(
      y = reorder(popnames, ranks),
      x = ranks
    ))
  }

  # To understand geom_errorbar behavior, it's better to imagine example with barplots and errorbars
  # The absolute width of bars depends on amount of categories. The more categories, the narrower bars
  # `width` in geom_errorbar is width relative to width of associated bar
  # So again, the more categories (or populations in our case), the narrower errorbar
  
  desired_errorbar_width <- .0125 # This is a value, which looked nice for PISA example
  errorbar_width <- desired_errorbar_width * p # Inside the geom_errorbar, it will be divided back to desired value

  pl <- pl +
    theme_bw() +
    geom_point()
  
  if(nrow(dat) >= 50){ # Chosen with trial and error
    pl <- pl + geom_errorbar(aes(xmin = L, xmax = U))
  } else {
    pl <- pl + geom_errorbar(aes(xmin = L, xmax = U),
                  width = errorbar_width,
                  linewidth = 1, position = position_dodge(0.1)) + 
      scale_x_continuous(limits = c(min_rank, max_rank), 
                         breaks = unique(c(min_rank, seq(ceiling(min_rank / 5)*5, max_rank, by = 5), max_rank)), 
                         labels = unique(c(min_rank, seq(ceiling(min_rank / 5)*5, max_rank, by = 5), max_rank)))
  }
    
  if (!horizontal) {
    pl <- pl + coord_flip()
  }
  if (colorbins > 1) {
    pl <- pl +
      labs(title = title, subtitle = subtitle, caption = caption, color = coltitle)
  } else {
    pl <- pl +
      labs(title = title, subtitle = subtitle, caption = caption)
  }
  
  pl <- pl + ylab(NULL) + xlab("rank")

  if (colorbins > 1) {
    pl <- pl +
      theme(
        legend.position = "inside",
        legend.position.inside = c(.95, .05),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right"
      ) +
      scale_color_manual(values = cc)
  }

  return(pl)
}
