#' Data: Positive and Negative Control Genes
#' 
#' Sets of "positive" and "negative" control genes, useful arguments for 
#' \code{\link{scone}}.
#' 
#' These gene sets can be used as negative or positive controls, either for RUV
#' factor normalization or for evaluation and ranking of the normalization
#' workflows.
#' 
#' @details Gene set datasets are in the form of \code{data.frame}, with the
#'   first column containing the gene symbols and an (optional) second column
#'   containing additional information (such as cortical layer or cell cycle
#'   phase).
#'   
#' @details Note that the gene symbols follow the mouse conventions (i.e. 
#'   capitalized) or the human conventions (i.e, all upper-case), based on the
#'   original publication. One can use the \code{\link[base]{toupper}}, 
#'   \code{\link[base]{tolower}}, and \code{\link[tools]{toTitleCase}} 
#'   functions to alter symbol conventions.
#'   
#' @details Mouse gene symbols in \code{cortical_markers} are transcribed from
#'   Figure 3 of Molyneaux et al. (2007): "laminar-specific expression of 66
#'   genes within the neocortex."
#'   
#' @details Human gene symbols in \code{housekeeping} are derived from the list
#'   of "housekeeping" genes from the cDNA microarray analysis of Eisenberg
#'   and Levanon (2003): "[HK genes] belong to the class of genes that are
#'   EXPRESSED in all tissues." "... from 47 different human tissues and cell
#'   lines."
#'   
#' @details Human gene symbols in \code{housekeeping_revised} from Eisenberg
#'   and Levanon (2013): "This list provided ... is based on analysis of
#'   next-generation sequencing (RNA-seq) data. At least one variant of these
#'   genes is expressed in all tissues uniformly... The RefSeq transcript
#'   according to which we deemed the gene 'housekeeping' is given."
#'   Housekeeping exons satisfy "(i) expression observed in all tissues; (ii)
#'   low variance over tissues: standard-deviation [log2(RPKM)]<1; and (iii) no
#'   exceptional expression in any single tissue; that is, no log-expression
#'   value differed from the averaged log2(RPKM) by two (fourfold) or more."
#'   "We define a housekeeping gene as a gene for which at least one RefSeq
#'   transcript has more than half of its exons meeting the previous criteria
#'   (thus being housekeeping exons)."
#'   
#' @details Human gene symbols in \code{cellcycle_genes} from Macosko et al.
#'   (2015) and represent a set of genes marking G1/S, S, G2/M, M, and M/G1 
#'   phases.
#'   
#' @references Molyneaux, B.J., Arlotta, P., Menezes, J.R. and Macklis, J.D.. 
#'   Neuronal subtype specification in the cerebral cortex. Nature Reviews 
#'   Neuroscience, 2007, 8(6):427-437.
#' @references Eisenberg E, Levanon EY. Human housekeeping genes are compact. 
#'   Trends in Genetics, 2003, 19(7):362-5.
#' @references Eisenberg E, Levanon EY. Human housekeeping genes, revisited. 
#'   Trends in Genetics, 2013, 29(10):569-74.
#' @references Macosko, E. Z., et al. Highly parallel genome-wide expression 
#'   profiling of individual cells using nanoliter droplets. Cell, 2015, 
#'   161.5:1202-1214.
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
