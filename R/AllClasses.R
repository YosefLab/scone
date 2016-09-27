#' Class SconeExperiment
#'
#' @description Objects of this class store,  at the minimum, the original gene
#'   expression matrix and a set of covariates (aka sample metadata) useful for
#'   running \code{\link{scone}}. These include, the quality control (QC)
#'   metrics, batch information, and biological classes of interest (if
#'   available).
#'
#' @description The typical way of creating \code{SconeExperiment} objects is
#'   via a call to the \code{\link{sconeExperiment}} function or to the
#'   \code{\link{scone}} function. If the object is a result to a
#'   \code{\link{scone}} call, it will contain the results, e.g., the
#'   performance metrics, scores, and the normalization schemes compared. (See
#'   Slots for a full list).
#'
#' @description This object extends the
#'   \code{\linkS4class{SummarizedExperiment}} class.
#'
#' @details The QC matrix, biological class, and batch information are stored as
#'   elements of the `colData` of the object.
#' @details The positive and negative control genes are stored as elements of
#'   the `rowData` of the object.
#'
#' @import methods
#' @import SummarizedExperiment
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#'
#' @name SconeExperiment-class
#' @import methods
#' @aliases SconeExperiment
#'
#' @export
#'
#' @slot which_qc integer. Index of columns of `colData` that contain the QC
#'   metrics.
#' @slot which_bio integer. Index of the column of `colData` that contains the
#'   biological classes information (it must be a factor).
#' @slot which_batch integer. Index of the column of `colData` that contains the
#'   batch information (it must be a factor).
#' @slot which_negconruv integer. Index of the column of `rowData` that contains
#'   a logical vector indicating which genes to use as negative controls to
#'   infer the factors of unwanted variation in RUV.
#' @slot which_negconeval integer. Index of the column of `rowData` that
#'   contains a logical vector indicating which genes to use as negative
#'   controls to evaluate the performance of the normalizations.
#' @slot which_poscon integer. Index of the column of `rowData` that contains a
#'   logical vector indicating which genes to use as positive controls to
#'   evaluate the performance of the normalizations.
#' @slot hdf5_pointer character. A string specifying to which file to write /
#'   read the normalized data.
#' @slot scone_metrics matrix. Matrix containing the "raw" performance metrics.
#'   See \code{\link{scone}} for a description of each metric.
#' @slot scone_scores matrix. Matrix containing the performance scores
#'   (transformed metrics). See \code{\link{scone}} for a discussion on the
#'   difference between scores and metrics.
#' @slot scone_parmas data.frame. A data frame containing the normalization
#'   schemes applied to the data and compared.
#' @slot design_mats list. A list of design matrices, one for each normalization
#'   scheme.
#' @slot scone_run character Whether \code{\link{scone}} and in which mode
#'   ("no", "in_memory", "hdf5").
setClass(
  Class = "SconeExperiment",
  contains = "SummarizedExperiment",
  slots = list(
    which_qc = "integer",
    which_bio = "integer",
    which_batch = "integer",
    which_negconruv = "integer",
    which_negconeval = "integer",
    which_poscon = "integer",
    hdf5_pointer = "character",
    scone_metrics = "matrix",
    scone_scores = "matrix",
    scone_params = "data.frame",
    design_mats = "list",
    scone_run = "character"
  )
)

setValidity("SconeExperiment", function(object) {

  ## check that the indeces are not out of bounds
  if(!all(which_qc %in% seq_len(NCOL(colData(object))))) {
    return("`which_qc` index out of bounds.")
  }

  if(!all(which_bio %in% seq_len(NCOL(colData(object))))) {
    return("`which_bio` index out of bounds.")
  }

  if(!all(which_batch %in% seq_len(NCOL(colData(object))))) {
    return("`which_batch` index out of bounds.")
  }

  if(!all(which_negconruv %in% seq_len(NCOL(rowData(object))))) {
    return("`which_negconruv` index out of bounds.")
  }

  if(!all(which_negconeval %in% seq_len(NCOL(rowData(object))))) {
    return("`which_negconeval` index out of bounds.")
  }

  if(!all(which_poscon %in% seq_len(NCOL(rowData(object))))) {
    return("`which_poscon` index out of bounds.")
  }

  ## check that which_bio and which_batch are of length 1
  if(length(which_bio) > 1) {
    return("Only one `bio` variable can be specified.")
  }
  if(length(which_batch) > 1) {
    return("Only one `batch` variable can be specified.")
  }

  ## check that which_poscon and which_negcon are of length 1
  if(length(which_negconruv) > 1) {
    return("Only one set of negative controls for RUV can be specified.")
  }
  if(length(which_negconeval) > 1) {
    return("Only one set of negative controls for evaluation can be specified.")
  }
  if(length(which_poscon) > 1) {
    return("Only one set of positive controls can be specified.")
  }

  ## check that all QC columns are numeric
  if(any(lapply(colData(object)[,which_qc], class) != "numeric")) {
    return("Only numeric QC metrics are allowed.")
  }

  ## check that bio is a factor
  if(!is.factor(colData(object)[,which_bio])) {
    return("`bio` must be a factor.")
  }

  ## check that batch is a factor
  if(!is.factor(colData(object)[,which_batch])) {
    return("`batch` must be a factor.")
  }

  ## check that poscon and negcon are logical
  if(!is.logical(rowData(object)[,which_negconruv])) {
    return("`negconruv` must be a logical vector.")
  }
  if(!is.logical(rowData(object)[,which_negconeval])) {
    return("`negconeval` must be a logical vector.")
  }
  if(!is.logical(rowData(object)[,which_poscon])) {
    return("`poscon` must be a logical vector.")
  }

  ## check if hdf5 file exist (depends on scone_run)
  if(scone_run == "hdf5" && !file.exists(hdf5_pointer)) {
    return(paste0("File ", hdf5_pointer, " not found."))
  }
  if(scone_run == "no" && file.exists(hdf5_pointer)) {
    return(paste0("File ", hdf5_pointer, " exists. Please specify a new file."))
  }
})
