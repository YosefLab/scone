#' Function for biplotting with no point labels and with
#' points color-coded according to a quantitative variable.
#' For example: the rank of normalization performance.
#' 
#' This function implements biplot for \code{\link[stats]{prcomp}} objects.
#' 
#' @param x \code{\link[stats]{prcomp}} object.
#' @param y numeric. Quantitative values used to color the points. If rank is 
#'   FALSE, all values must be positive integers and less than or equal to the 
#'   length of y.
#' @param rank logical. If TRUE (default) y will be transformed by the rank() 
#'   function
#' @param ties_method character. ties.method used by the rank() function
#' @param choices numeric. 2 principal components to plot. Default to first two
#'   PCs.
#' @param expand numeric. value used to adjust the spread of the arrows
#'   relative to the points.
#' @param ... arguments passed to plot.
#'   
#' @importFrom grDevices colorRampPalette
#' @export
#' 
#' @return Invisibly returns scaled point coordinates used in plot.
#'   
#' @examples
#' mat <- matrix(rnorm(1000), ncol=10)
#' colnames(mat) <- paste("X", 1:ncol(mat), sep="")
#' 
#' pc <- prcomp(mat)
#' 
#' biplot_color(pc, rank(pc$x[,1]))
#' 
biplot_color <- function(x, y, rank = TRUE, 
                         ties_method = c("max", "min", 
                                         "first", "last", "random"),
                         choices = 1:2, expand = 1, ...) {

  if(rank){
    
    ties_method <- match.arg(ties_method)
    y = rank(y,ties.method = ties_method)
    
  }else{
    
    if(any(abs(y - round(y)) > .Machine$double.eps^0.5)){
      stop("ranks must be integer")
    }else{y = as.integer(y)}
    
    if(any(y <= 0)){
      stop("ranks must be positive")
    }
    
    if(any(y > length(y))){
      stop("ranks must be less than or equal to total number of elements")
    }
    
  }
  
  lam <- x$sdev[choices]
  n <- NROW(x$x)
  lam <- lam * sqrt(n)

  xx <- t(t(x$x[, choices])/lam)
  yy <- t(t(x$rotation[, choices]) * lam)

  ratio <- max(range(yy)/range(xx))/expand

  cols <- rev(colorRampPalette(c("black","navyblue","mediumblue",
                                 "dodgerblue3","aquamarine4","green4",
                                 "yellowgreen","yellow"))(length(y)))[y]
  plot(xx, pch=19, col=cols, ...)

  labs <- rownames(yy)

  text(yy/ratio, labels=labs, col=2)
  arrows(0, 0, yy[, 1] * 0.8/ratio, yy[, 2] * 0.8/ratio, col = 2, length = 0.1)

  invisible(xx)
}
