#' Class SconeExperiment
#'
#' @description Objects of this class store, at minimum, a gene expression
#'   matrix, sample metadata, and feature data (e.g. gene sets) useful for
#'   running \code{\link{scone}}.
#'
#' @description A \code{SconeExperiment} objects is created via a call to the
#'   \code{\link{SconeExperiment}} constructor or the main \code{\link{scone}}
#'   function. If the object results from a \code{\link{scone}} call, it will
#'   contain normalization performance metrics, scores, and all data required
#'   to recover individual normalized expression matrices. (See Slots for a
#'   full list).
#'
#' @description This object extends the
#'   \code{\linkS4class{SummarizedExperiment}} class.
#'
#' @details The quality control (QC) metrics, biological class, and batch
#'   information are stored as elements of the `colData` of the object.
#' @details Positive and negative control gene sets are stored as elements of
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
#'   biological class information (it must be a factor).
#' @slot which_batch integer. Index of the column of `colData` that contains
#'   the batch information (it must be a factor).
#' @slot which_negconruv integer. Index of the column of `rowData` that
#'   contains a logical vector indicating which genes are negative controls for
#'   RUV normalization.
#' @slot which_negconeval integer. Index of the column of `rowData` that
#'   contains a logical vector indicating which genes are negative controls for
#'   evaluating the normalization performance.
#' @slot which_poscon integer. Index of the column of `rowData` that contains a
#'   logical vector indicating which genes are positive controls for evaluating
#'   the normalization performance.
#' @slot hdf5_pointer character. A string specifying to which file to write /
#'   read the normalized expression data.
#' @slot imputation_fn list of functions used by \code{\link{scone}} in the
#'   imputation step.
#' @slot scaling_fn list of functions used by \code{\link{scone}} in the
#'   scaling step.
#' @slot scone_metrics matrix. Matrix containing the \code{\link{scone}}
#'   performance metrics.
#' @slot scone_scores matrix. Matrix containing the \code{\link{scone}}
#'   performance scores (transformed metrics) and aggregate score.
#' @slot scone_params data.frame. A data frame containing specifications for
#'   normalization schemes applied by \code{\link{scone}}.
#' @slot scone_run character. Whether \code{\link{scone}} was run and in which
#'   mode ("no", "in_memory", "hdf5").
#' @slot is_log logical. Are the expression data in log- scale? (currently not
#'   supported by \code{\link{scone}})
#' @slot nested logical. Is batch nested within bio? (Automatically set by
#'   \code{\link{scone}}).
#' @slot rezero logical. Did \code{\link{scone}} restore zero values following
#'   the scaling step? TRUE if \code{\link{scone}} was run with
#'   \code{zero="preadjust"} or \code{zero="strong"}.
#' @slot fixzero logical. Did \code{\link{scone}} fix zero values following the
#'   adjustment step? TRUE if \code{\link{scone}} was run with
#'   \code{zero="postadjust"} or \code{zero="strong"}.
#' @slot impute_args list. Arguments passed to all imputation functions.
#'
#' @seealso \code{\link{get_normalized}}, \code{\link{get_params}},
#' \code{\link{get_batch}}, \code{\link{get_bio}}, \code{\link{get_design}},
#' \code{\link{get_negconeval}}, \code{\link{get_negconruv}},
#' \code{\link{get_poscon}}, \code{\link{get_qc}},
#' \code{\link{get_scores}}, and \code{\link{get_score_ranks}}
#' to access internal fields, \code{\link{select_methods}} for subsetting
#' by method, and \code{\link{scone}} for running scone workflows.
#'
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
    imputation_fn = "list",
    scaling_fn = "list",
    scone_metrics = "matrix",
    scone_scores = "matrix",
    scone_params = "data.frame",
    scone_run = "character",
    is_log = "logical",
    nested = "logical",
    rezero = "logical",
    fixzero = "logical",
    impute_args = "list"
  )
)

setValidity("SconeExperiment", function(object) {
  ## check that the indeces are not out of bounds
  if (!all(object@which_qc %in% seq_len(NCOL(colData(object))))) {
    return("`which_qc` index out of bounds.")
  }
  
  if (!all(object@which_bio %in% seq_len(NCOL(colData(object))))) {
    return("`which_bio` index out of bounds.")
  }
  
  if (!all(object@which_batch %in% seq_len(NCOL(colData(object))))) {
    return("`which_batch` index out of bounds.")
  }
  
  if (!all(object@which_negconruv %in% seq_len(NCOL(rowData(object))))) {
    return("`which_negconruv` index out of bounds.")
  }
  
  if (!all(object@which_negconeval %in% seq_len(NCOL(rowData(object))))) {
    return("`which_negconeval` index out of bounds.")
  }
  
  if (!all(object@which_poscon %in% seq_len(NCOL(rowData(object))))) {
    return("`which_poscon` index out of bounds.")
  }
  
  ## check that which_bio and which_batch are of length 1
  if (length(object@which_bio) > 1) {
    return("Only one `bio` variable can be specified.")
  }
  if (length(object@which_batch) > 1) {
    return("Only one `batch` variable can be specified.")
  }
  
  ## check that which_poscon and which_negcon are of length 1
  if (length(object@which_negconruv) > 1) {
    return("Only one set of negative controls for RUV can be specified.")
  }
  if (length(object@which_negconeval) > 1) {
    return(paste0(
      "Only one set of negative controls ",
      "for evaluation can be specified."
    ))
  }
  if (length(object@which_poscon) > 1) {
    return("Only one set of positive controls can be specified.")
  }
  
  ## check that all QC columns are numeric
  if (length(object@which_qc) > 0) {
    if (any(lapply(get_qc(object), class) != "numeric")) {
      return("Only numeric QC metrics are allowed.")
    }
  }
  
  ## check that bio is a factor
  if (length(object@which_bio) > 0) {
    if (!is.factor(get_bio(object))) {
      return("`bio` must be a factor.")
    }
  }
  
  ## check that batch is a factor
  if (length(object@which_batch) > 0) {
    if (!is.factor(get_batch(object))) {
      return("`batch` must be a factor.")
    }
  }
  
  ## check that poscon and negcon are logical
  if (length(object@which_negconruv) > 0) {
    if (!is.logical(get_negconruv(object))) {
      return("`negconruv` must be a logical vector.")
    }
  }
  if (length(object@which_negconeval) > 0) {
    if (!is.logical(get_negconeval(object))) {
      return("`negconeval` must be a logical vector.")
    }
  }
  if (length(object@which_poscon) > 0) {
    if (!is.logical(get_poscon(object))) {
      return("`poscon` must be a logical vector.")
    }
  }
  
  ## check if hdf5 file exist (depends on scone_run)
  if (length(object@hdf5_pointer) > 0) {
    if (object@scone_run == "hdf5" &&
        !file.exists(object@hdf5_pointer)) {
      return(paste0("File ", object@hdf5_pointer, " not found."))
    }
    if (object@scone_run == "no" &&
        file.exists(object@hdf5_pointer)) {
      return(paste0(
        "File ",
        object@hdf5_pointer,
        " exists. Please specify a new file."
      ))
    }
  }
  
  return(TRUE)
  
})

## Constructor

#' @rdname SconeExperiment-class
#'
#' @description The constructor \code{SconeExperiment} creates an object of the
#'   class \code{SconeExperiment}.
#'
#' @param object Either a matrix or a \code{\link{SummarizedExperiment}}
#'   containing the raw gene expression.
#' @param ... see specific S4 methods for additional arguments.
#'
#' @export
#'
#' @examples
#' set.seed(42)
#' nrows <- 200
#' ncols <- 6
#' counts <- matrix(rpois(nrows * ncols, lambda=10), nrows)
#' rowdata <- data.frame(poscon=c(rep(TRUE, 10), rep(FALSE, nrows-10)))
#' coldata <- data.frame(bio=gl(2, 3))
#' se <- SummarizedExperiment(assays=SimpleList(counts=counts),
#'                           rowData=rowdata, colData=coldata)
#'
#' scone1 <- SconeExperiment(assay(se), bio=coldata$bio, poscon=rowdata$poscon)
#'
#' scone2 <- SconeExperiment(se, which_bio=1L, which_poscon=1L)
#'
#'
setGeneric(
  name = "SconeExperiment",
  def = function(object, ...) {
    standardGeneric("SconeExperiment")
  }
)

#' @rdname SconeExperiment-class
#'
#' @param which_qc integer. Index of columns of `colData` that contain the QC
#'   metrics.
#' @param which_bio integer. Index of the column of `colData` that contains the
#'   biological class information (it must be a factor).
#' @param which_batch integer. Index of the column of `colData` that contains
#'   the batch information (it must be a factor).
#' @param which_negconruv integer. Index of the column of `rowData` that
#'   contains a logical vector indicating which genes are negative controls for
#'   RUV normalization.
#' @param which_negconeval integer. Index of the column of `rowData` that
#'   contains a logical vector indicating which genes are negative controls for
#'   evaluating the normalization performance.
#' @param which_poscon integer. Index of the column of `rowData` that contains
#'   a logical vector indicating which genes are positive controls for
#'   evaluating the normalization performance.
#' @param is_log logical. Are the expression data in log- scale? (currently not
#'   supported by \code{\link{scone}})
#'
#' @export
#'
setMethod(
  f = "SconeExperiment",
  signature = signature("SummarizedExperiment"),
  definition = function(object,
                        which_qc = integer(),
                        which_bio = integer(),
                        which_batch = integer(),
                        which_negconruv = integer(),
                        which_negconeval = integer(),
                        which_poscon = integer(),
                        is_log = FALSE) {
    out <- new(
      "SconeExperiment",
      object,
      which_qc = which_qc,
      which_bio = which_bio,
      which_batch = which_batch,
      which_negconruv = which_negconruv,
      which_negconeval = which_negconeval,
      which_poscon = which_poscon,
      hdf5_pointer = character(),
      imputation_fn = list(),
      scaling_fn = list(),
      scone_metrics = matrix(),
      scone_scores = matrix(),
      scone_params = data.frame(),
      scone_run = "no",
      is_log = is_log,
      nested = FALSE,
      rezero = FALSE,
      fixzero = FALSE,
      impute_args = list()
    )
    
    validObject(out)
    return(out)
  }
)


#' @rdname SconeExperiment-class
#'
#' @param qc matrix. numeric matrix of QC metrics (samples-by-metrics).
#' @param bio factor. biological class of interest.
#' @param batch factor. batch information.
#' @param negcon_ruv logical. vector indicating which genes are negative
#'   controls for RUV.
#' @param negcon_eval logical. vector indicating which genes are negative
#'   controls for evaluation.
#' @param poscon logical. vector indicating which genes are positive controls.
#'
#' @export
#'
#' @return A \code{\link{SconeExperiment}} object.
#'
setMethod(
  f = "SconeExperiment",
  signature = signature("matrix"),
  definition = function(object,
                        qc,
                        bio,
                        batch,
                        negcon_ruv = NULL,
                        negcon_eval = negcon_ruv,
                        poscon = NULL,
                        is_log = FALSE) {
    which_qc <- which_bio <- which_batch <- integer()
    which_negconruv <- which_negconeval <- which_poscon <- integer()
    
    coldata <- as.data.frame(matrix(nrow = NCOL(object), ncol = 0))
    
    if (!missing(qc)) {
      coldata <- cbind(coldata, qc)
      which_qc <- seq_len(NCOL(coldata))
    }
    
    if (!missing(bio)) {
      coldata <- cbind(coldata, bio)
      which_bio <- NCOL(coldata)
    }
    
    if (!missing(batch)) {
      coldata <- cbind(coldata, batch)
      which_batch <- NCOL(coldata)
    }
    
    rowdata <- as.data.frame(matrix(nrow = NROW(object), ncol = 0))
    
    if (!is.null(negcon_ruv)) {
      rowdata <- cbind(rowdata, negcon_ruv)
      which_negconruv <- NCOL(rowdata)
    }
    
    if (!is.null(negcon_eval)) {
      rowdata <- cbind(rowdata, negcon_eval)
      which_negconeval <- NCOL(rowdata)
    }
    
    if (!is.null(poscon)) {
      rowdata <- cbind(rowdata, poscon)
      which_poscon <- NCOL(rowdata)
    }
    
    se <-
      SummarizedExperiment(object, rowData = rowdata, colData = coldata)
    SconeExperiment(
      se,
      which_qc,
      which_bio,
      which_batch,
      which_negconruv,
      which_negconeval,
      which_poscon,
      is_log
    )
  }
)

