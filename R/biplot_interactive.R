#' Interactive biplot
#'
#' This is a wrapper around \code{\link{biplot_colored}}, which creates a shiny
#' gadget to allow the user to select specific points in the graph.
#'
#' @details Since this is based on the shiny gadget feature, it will not work in
#'   static documents, such as vignettes or markdown / knitr documents.
#'   See \code{biplot_colored} for more details on the internals.
#'
#' @param data a data.frame containing the data to be plotted.
#' @param scores a numeric vector used to color the points.
#' @param ... passed to \code{\link{biplot_colored}}.
#'
#' @importFrom miniUI gadgetTitleBar miniContentPanel miniPage gadgetTitleBar
#' @importFrom shiny plotOutput renderPlot observeEvent brushedPoints runGadget verbatimTextOutput stopApp renderText
#'
#' @export
#'
#' @examples
#' \dontrun{
#' mat <- matrix(rnorm(1000), ncol=10)
#' colnames(mat) <- paste("X", 1:ncol(mat), sep="")
#'
#' biplot_interactive(mat, mat[,1])
#' }
biplot_interactive <- function(data, scores, ...) {

  data <- as.data.frame(data)
  scores <- as.numeric(scores)

  ui <- miniPage(
    gadgetTitleBar("Drag to select points"),
    miniContentPanel(
      # The brush="brush" argument means we can listen for
      # brush events on the plot using input$brush.
      plotOutput("plot1", height = "80%", brush = "plot_brush"),
      verbatimTextOutput("info")
    )
  )

  server <- function(input, output, session) {

    # Compute PCA
    pc_obj <- prcomp(data, center = TRUE, scale = FALSE)
    bp_obj <- biplot_colored(pc_obj, y = scores)

    # Render the plot
    output$plot1 <- renderPlot({
      # Biplot
      biplot_colored(pc_obj, y = scores, ...)
    })

    data_out <- cbind(data, bp_obj)

    output$info <- renderText({
      xy_range_str <- function(e) {
        if(is.null(e)) return("NULL\n")
        idx <- which(bp_obj[,1] >= e$xmin & bp_obj[,1] <= e$xmax &
                       bp_obj[,2] >= e$ymin & bp_obj[,2] <= e$ymax)
        paste0(rownames(data)[idx], collapse = "\n")
      }
      xy_range_str(input$plot_brush)
    })

    # Handle the Done button being pressed.
    observeEvent(input$done, {
      # Return the brushed points. See ?shiny::brushedPoints.
      stopApp(brushedPoints(data_out, input$plot_brush, xvar="PC1", yvar="PC2"))
    })
  }

  runGadget(ui, server)
}
