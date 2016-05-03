#' Positive and negative control genes
#'
#' Sets of positive and negative control genes, useful for input in
#' \code{\link{scone}}.
#'
#' These gene sets can be used as negative or positive controls, either for RUV
#' normalization or for evaluation and ranking of the normalization strategies.
#'
#' @details The datasets are in the form of \code{data.frame}, with at least one
#'   column containing the gene symbols and optionally a second column
#'   containing additional information (such as cortical layer or cell cycle
#'   phase).
#'
#' @details Note that the gene symbols follow the mouse conventions (i.e.,
#'   capitalized) or the human conventions (i.e, all upper-case), based on the
#'   original pubblication. One can use the \code{\link[base]{toupper}},
#'   \code{\link[base]{tolower}}, and \code{\link[tools]{toTitleCase}} functions
#'   to convert the symbols.
#'
#' @details The genes in \code{cortical_markers} are from Figure 3 of Molyneaux
#'   et al. (2007). The genes in \code{housekeeping} are from Eisenberg and
#'   Levanon (2003) and in \code{housekeeping_revised} are from Eisenberg and
#'   Levanon (2013). The genes in \code{cellcycle_genes} are adapted from
#'   Kowalczyk et al. (2015).
#'
#' @references Molyneaux, B.J., Arlotta, P., Menezes, J.R. and Macklis, J.D..
#'   Neuronal subtype specification in the cerebral cortex. Nature Reviews
#'   Neuroscience, 2007, 8(6):427-437.
#' @references Eisenberg E, Levanon EY. Human housekeeping genes are compact.
#'   Trends in Genetics, 2003, 19(7):362-5.
#' @references Eisenberg E, Levanon EY. Human housekeeping genes, revisited.
#'   Trends in Genetics, 2013, 29(10):569-74.
#' @references Kowalczyk, M.S., Tirosh, I., Heckl, D., Rao, T.N., Dixit, A.,
#'   Haas, B.J., Schneider, R.K., Wagers, A.J., Ebert, B.L. and Regev, A.
#'   Single-cell RNA-seq reveals changes in cell cycle and differentiation
#'   programs upon aging of hematopoietic stem cells. Genome research, 2015.
#'
#' @name control_genes
#'
#' @docType data
#' @aliases cortical_markers housekeeping housekeeping_revised cellcycle_genes
#'
#' @examples
#' data(housekeeping)
#' data(housekeeping_revised)
#' data(cellcycle_genes)
#' data(cortical_markers)
NULL
