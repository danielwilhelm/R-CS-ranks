#' Plot ranking with confidence sets
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
#' x <- seq(1,3,length=10)
#' sd <- rep(0.2,10)
#' ranks <- xrank(x)
#' CS <- csranks(x, sd)
#' grid::current.viewport()
#' plotranking(ranks, CS$L, CS$U)

#' @export
#' @import ggplot2
#' @importFrom stats reorder
plotranking <- function(ranks, L, U, popnames=NULL, title=NULL, subtitle=NULL, 
                        caption=NULL, colorbins=1, horizontal=TRUE) {

	# initializations
	cc <- scales::seq_gradient_pal("#66a182", "#d1495b", "Lab")(seq(0,1,length.out=4))
	p <- length(ranks)
	stopifnot(colorbins>=1 & colorbins<=p)
	if (is.null(popnames)) popnames <- 1:p
	dat <- data.frame(ranks=ranks, L=L, U=U, popnames=popnames)
	coltitle <- ifelse(colorbins==4, "quartile", "quantile bins")

	# plot
	if (colorbins>1) {
		pl <- ggplot(data=dat, aes(x=reorder(popnames, ranks),
		                           y=ranks, color=createbins(ranks,colorbins)))	
	} else {
		pl <- ggplot(data=dat, aes(x=reorder(popnames, ranks), 
		                           y=ranks))	
	}
	
	# To understand geom_errorbar behavior, it's better to imagine example with barplots and errorbars
	# The absolute width of bars depends on amount of categories. The more categories, the narrower bars
	# `width` in geom_errorbar is width relative to width of associated bar
	# So again, the more categories (or populations in our case), the narrower errorbar
	desired_errorbar_width <- .0125 # This is a value, which looked nice for PISA example
	errorbar_width <- desired_errorbar_width * p # Inside the geom_errorbar, it will be divided back to desired value
	
	pl = pl +
	  	theme_bw() + 
	  	geom_point() +
	  	geom_errorbar(aes(ymin=L, ymax=U), width=errorbar_width, 
	  	              size=1, position=position_dodge(0.1))
  if(horizontal)
    pl <- pl + coord_flip()
 	if (colorbins>1) {
 		pl = pl +
 			labs(title = title, subtitle = subtitle, caption=caption, color=coltitle)
 	} else {
 		pl = pl +
 			labs(title = title, subtitle = subtitle, caption=caption)
 	}

 	pl = pl +
	  	scale_y_continuous(name="rank", limits=c(1, p), breaks=unique(c(1,seq(5,p,by=5),p)),labels=unique(c(1,seq(5,p,by=5),p))) +
	  	xlab(NULL)

	if (colorbins>1) {
		pl = pl + 
			theme(legend.position = c(.95, .05),
		    legend.justification = c("right", "bottom"),
		    legend.box.just = "right") +
		  	scale_color_manual(values=cc)
	}
	  
	 return(pl)
}
