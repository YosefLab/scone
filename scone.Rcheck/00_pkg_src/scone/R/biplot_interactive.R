#' Interactive biplot
#'
#' This is a wrapper around \code{\link{biplot_color}}, creating a shiny gadget
#' to allow the user to select specific points in the graph.
#'
#' @details Since this is based on the shiny gadget feature, it will not work
#'   in static documents, such as vignettes or markdown / knitr documents. See
#'   \code{biplot_color} for more details on the internals.
#'
#' @param x a \code{\link{SconeExperiment}} object.
#' @param ... passed to \code{\link{biplot_color}}.
#'
#' @export
#'
#' @return A \code{\link{SconeExperiment}} object representing
#'   selected methods.
#'
#' @examples
#' mat <- matrix(rpois(1000, lambda = 5), ncol=10)
#' colnames(mat) <- paste("X", 1:ncol(mat), sep="")
#' obj <- SconeExperiment(mat)
#' res <- scone(obj, scaling=list(none=identity,
#'    uq=UQ_FN, deseq=DESEQ_FN,  fq=FQT_FN),
#' evaluate=TRUE, k_ruv=0, k_qc=0, eval_kclust=2,
#'    bpparam = BiocParallel::SerialParam())
#' \dontrun{
#' biplot_interactive(res)
#' }
#'
biplot_interactive <- function(x, ...) {

  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("shiny package needed for biplot_interactive")
  }

  if (!requireNamespace("miniUI", quietly = TRUE)) {
    stop("miniUI package needed for biplot_interactive")
  }

  data <- as.data.frame(apply(t(get_scores(x)),1,rank))
  scores <- get_score_ranks(x)

  ui <- miniUI::miniPage(
    miniUI::gadgetTitleBar("Drag to select points"),
    miniUI::miniContentPanel(
      # The brush="brush" argument means we can listen for
      # brush events on the plot using input$brush.
      shiny::plotOutput("plot1", height = "80%", brush = "plot_brush"),
      shiny::verbatimTextOutput("info")
    )
  )

  server <- function(input, output, session) {

    # Compute PCA
    pc_obj <- prcomp(data, center = TRUE, scale = FALSE)
    bp_obj <- biplot_color(pc_obj, y = scores)

    # Render the plot
    output$plot1 <- shiny::renderPlot({
      # Biplot
      biplot_color(pc_obj, y = scores, ...)
    })

    data_out <- cbind(data, bp_obj)

    output$info <- shiny::renderText({
      xy_range_str <- function(e) {
        if(is.null(e)) return("NULL\n")
        idx <- which(bp_obj[,1] >= e$xmin & bp_obj[,1] <= e$xmax &
                       bp_obj[,2] >= e$ymin & bp_obj[,2] <= e$ymax)
        paste0(rownames(data)[idx], collapse = "\n")
      }
      xy_range_str(input$plot_brush)
    })

    # Handle the Done button being pressed.
    shiny::observeEvent(input$done, {
      # Return the brushed points. See ?shiny::brushedPoints.
      names <- rownames(shiny::brushedPoints(data_out,
                                      input$plot_brush,
                                      xvar="PC1",
                                      yvar="PC2"))
      out <- select_methods(x, names)
      shiny::stopApp(invisible(out))
    })
  }

  shiny::runGadget(ui, server)
}
