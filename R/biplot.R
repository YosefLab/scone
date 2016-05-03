## When porting the package to S4, we can make this a biplot method

#' Another implementation of the biplot function
#'
#' Function to plot a biplot with no point labels and with points color-coded
#' according to a certain quantitative variable, for instance the rank of the
#' normalization performance or the expression of a certain gene.
#'
#' This function implements the biplot only for \code{\link[stats]{prcomp}}
#' objects. Eventually, we will turn this into an S4 method.
#'
#' @param x the result of a call to \code{\link[stats]{prcomp}}.
#' @param y the rank value that should be used to color the points.
#' @param choices which principal components to plot. Only 2D plots are
#'   possible for now. Default to first two PCs.
#' @param expand numeric value used to adjust the spread of the arrows relative
#'   to the points.
#' @param ... passed to plot.
#'
#' @importFrom grDevices colorRampPalette
#' @export
#'
#' @examples
#' mat <- matrix(rnorm(1000), ncol=10)
#' colnames(mat) <- paste("X", 1:ncol(mat), sep="")
#'
#' pc <- prcomp(mat)
#'
#' biplot_colored(pc, rank(pc$x[,1]))
#'
biplot_colored <- function(x, y, choices=1:2, expand=1, ...) {

  lam <- x$sdev[choices]
  n <- NROW(x$x)
  lam <- lam * sqrt(n)

  xx <- t(t(x$x[, choices])/lam)
  yy <- t(t(x$rotation[, choices]) * lam)

  ratio <- max(range(yy)/range(xx))/expand

  cols <- rev(colorRampPalette(c("black","navyblue","mediumblue","dodgerblue3","aquamarine4","green4","yellowgreen","yellow"))(length(y)))[y]
  plot(xx, pch=19, col=cols, ...)

  labs <- rownames(yy)

  text(yy/ratio, labels=labs, col=2)
  arrows(0, 0, yy[, 1] * 0.8/ratio, yy[, 2] * 0.8/ratio, col = 2, length = 0.1)

  invisible(xx)
}
