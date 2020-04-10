#' Plot ranking with confidence sets
#'
#' @param ranks vector of ranks
#' @param L vector of lower bounds of confidence sets for the ranks
#' @param U vector of lower bounds of confidence sets for the ranks
#' @param popnames vector containing names of the populations whose ranks are in \code{ranks}. If \code{popnames=NULL} (default), then populations are automatically numbered.
#' @param title character string containing the main title of the graph. \code{title=NULL} (default) means no title.
#' @param subtitle character string containing the subtitle of the graph. \code{subtitle=NULL} (default) means no subtitle.
#' @param caption character string containing the caption of the graph. \code{caption=NULL} (default) means no caption.

#' @return ggplot

#' @examples
#' x <- seq(1,3,length=10)
#' sd <- rep(0.2,10)
#' ranks <- xrank(x)
#' CS <- csranks(x, sd)
#' grid::current.viewport()
#' plotranking(ranks, CS$L, CS$U)

#' @export
#' @import ggplot2
#' @importFrom stats reorder
plotranking <- function(ranks, L, U, popnames=NULL, title=NULL, subtitle=NULL, caption=NULL) {

	# initializations
	cc <- scales::seq_gradient_pal("#66a182", "#d1495b", "Lab")(seq(0,1,length.out=4))
	p <- length(ranks)
	if (is.null(popnames)) popnames <- 1:p
	dat <- data.frame(ranks=ranks, L=L, U=U, popnames=popnames)

	# plot
	pl <- ggplot(data=dat, aes(x=reorder(popnames, ranks), y=ranks, color=createquartiles(ranks))) + 
	  theme_bw() + 
	  geom_point() +
	  geom_errorbar(aes(ymin=L, ymax=U), width=.5, size=1, position=position_dodge(0.1)) +
	  coord_flip() +
	  labs(title = title, subtitle = subtitle, caption=caption, color="quartile") +
	  scale_y_continuous(name="rank", limits=c(1, p), breaks=unique(c(1,seq(5,p,by=5),p)),labels=unique(c(1,seq(5,p,by=5),p))) +
	  xlab(NULL) +
	  theme(legend.position = c(.95, .05),
	    legend.justification = c("right", "bottom"),
	    legend.box.just = "right") +
	  scale_color_manual(values=cc)

	 return(pl)
}