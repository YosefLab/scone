#' Class SconeExperiment
#' 
#' @description Objects of this class store, at minimum, a gene expression 
#'   matrix and a set of covariates (sample metadata) useful for running 
#'   \code{\link{scone}}. These include, the quality control (QC) metrics,
#'   batch information, and biological classes of interest (if available).
#'   
#' @description The typical way of creating \code{SconeExperiment} objects is 
#'   via a call to the \code{\link{SconeExperiment}} function or to the 
#'   \code{\link{scone}} function. If the object is a result to a 
#'   \code{\link{scone}} call, it will contain the results, e.g., the 
#'   performance metrics, scores, and normalization workflow comparisons. (See 
#'   Slots for a full list).
#'   
#' @description This object extends the 
#'   \code{\linkS4class{SummarizedExperiment}} class.
#'   
#' @details The QC matrix, biological class, and batch information are 
#'   stored as elements of the `colData` of the object.
#' @details The positive and negative control genes are stored as 
#'   elements of the `rowData` of the object.
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
#' @slot which_qc integer. Index of columns of `colData` that contain the
#'   QC metrics.
#' @slot which_bio integer. Index of the column of `colData` that contains
#'   the biological classes information (it must be a factor).
#' @slot which_batch integer. Index of the column of `colData`
#'   that contains the batch information (it must be a factor).
#' @slot which_negconruv integer. Index of the column of `rowData` that
#'   contains a logical vector indicating which genes to use as negative
#'   controls to infer the factors of unwanted variation in RUV.
#' @slot which_negconeval integer. Index of the column of `rowData` that 
#'   contains a logical vector indicating which genes to use as negative 
#'   controls to evaluate the performance of the normalizations.
#' @slot which_poscon integer. Index of the column of `rowData` that
#'   contains a logical vector indicating which genes to use as positive
#'   controls to evaluate the performance of the normalizations.
#' @slot hdf5_pointer character. A string specifying to which 
#'   file to write / read the normalized data.
#' @slot imputation_fn list of functions used by scone for 
#'   the imputation step.
#' @slot scaling_fn list of functions used by scone for the scaling step.
#' @slot scone_metrics matrix. Matrix containing the "raw" 
#'   performance metrics. See \code{\link{scone}} for a 
#'   description of each metric.
#' @slot scone_scores matrix. Matrix containing the performance scores 
#'   (transformed metrics). See \code{\link{scone}} for a discussion on the 
#'   difference between scores and metrics.
#' @slot scone_params data.frame. A data frame containing
#'   the normalization schemes applied to the data and compared.
#' @slot scone_run character. Whether \code{\link{scone}} was 
#'   run and in which mode ("no", "in_memory", "hdf5").
#' @slot is_log logical. Are the expression data in log scale?
#' @slot nested logical. Is batch nested within bio? 
#'   (Automatically set by \code{\link{scone}}).
#' @slot rezero logical. TRUE if \code{\link{scone}} was run with 
#'   \code{rezero=TRUE}.
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
    impute_args = "list"
  )
)

setValidity("SconeExperiment", function(object) {

  ## check that the indeces are not out of bounds
  if(!all(object@which_qc %in% seq_len(NCOL(colData(object))))) {
    return("`which_qc` index out of bounds.")
  }

  if(!all(object@which_bio %in% seq_len(NCOL(colData(object))))) {
    return("`which_bio` index out of bounds.")
  }

  if(!all(object@which_batch %in% seq_len(NCOL(colData(object))))) {
    return("`which_batch` index out of bounds.")
  }

  if(!all(object@which_negconruv %in% seq_len(NCOL(rowData(object))))) {
    return("`which_negconruv` index out of bounds.")
  }

  if(!all(object@which_negconeval %in% seq_len(NCOL(rowData(object))))) {
    return("`which_negconeval` index out of bounds.")
  }

  if(!all(object@which_poscon %in% seq_len(NCOL(rowData(object))))) {
    return("`which_poscon` index out of bounds.")
  }

  ## check that which_bio and which_batch are of length 1
  if(length(object@which_bio) > 1) {
    return("Only one `bio` variable can be specified.")
  }
  if(length(object@which_batch) > 1) {
    return("Only one `batch` variable can be specified.")
  }

  ## check that which_poscon and which_negcon are of length 1
  if(length(object@which_negconruv) > 1) {
    return("Only one set of negative controls for RUV can be specified.")
  }
  if(length(object@which_negconeval) > 1) {
    return(paste0("Only one set of negative controls ",
                  "for evaluation can be specified."))
  }
  if(length(object@which_poscon) > 1) {
    return("Only one set of positive controls can be specified.")
  }

  ## check that all QC columns are numeric
  if(length(object@which_qc) > 0) {
    if(any(lapply(get_qc(object), class) != "numeric")) {
      return("Only numeric QC metrics are allowed.")
    }
  }

  ## check that bio is a factor
  if(length(object@which_bio) > 0) {
    if(!is.factor(get_bio(object))) {
      return("`bio` must be a factor.")
    }
  }

  ## check that batch is a factor
  if(length(object@which_batch) > 0) {
    if(!is.factor(get_batch(object))) {
      return("`batch` must be a factor.")
    }
  }

  ## check that poscon and negcon are logical
  if(length(object@which_negconruv) > 0) {
    if(!is.logical(get_negconruv(object))) {
      return("`negconruv` must be a logical vector.")
    }
  }
  if(length(object@which_negconeval) > 0) {
    if(!is.logical(get_negconeval(object))) {
      return("`negconeval` must be a logical vector.")
    }
  }
  if(length(object@which_poscon) > 0) {
    if(!is.logical(get_poscon(object))) {
      return("`poscon` must be a logical vector.")
    }
  }

  ## check if hdf5 file exist (depends on scone_run)
  if(length(object@hdf5_pointer) > 0) {
    if(object@scone_run == "hdf5" && !file.exists(object@hdf5_pointer)) {
      return(paste0("File ", object@hdf5_pointer, " not found."))
    }
    if(object@scone_run == "no" && file.exists(object@hdf5_pointer)) {
      return(paste0("File ", object@hdf5_pointer,
                    " exists. Please specify a new file."))
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
#' @param which_qc index that specifies which columns of `colData` 
#'   correspond to QC measures.
#' @param which_bio index that specifies which column of `colData`
#'   corresponds to `bio`.
#' @param which_batch index that specifies which column of `colData`
#'   corresponds to `batch`.
#' @param which_negconruv index that specifies which column of `rowData`
#'   has information on negative controls for RUV.
#' @param which_negconeval index that specifies which column of `rowData`
#'   has information on negative controls for evaluation.
#' @param which_poscon index that specifies which column of `rowData` has 
#'   information on positive controls.
#' @param is_log are the expression data in log scale?
#' @export
#' 
setMethod(
  f = "SconeExperiment",
  signature = signature("SummarizedExperiment"),
  definition = function(object, which_qc=integer(), which_bio=integer(),
                        which_batch=integer(),
                        which_negconruv=integer(), which_negconeval=integer(),
                        which_poscon=integer(), is_log=FALSE) {

    out <- new("SconeExperiment",
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
               impute_args = list()
               )

    validObject(out)
    return(out)
  }
)


#' @rdname SconeExperiment-class
#'   
#' @param qc numeric matrix with the QC measures.
#' @param bio factor with the biological class of interest.
#' @param batch factor with the batch information.
#' @param negcon_ruv a logical vector indicating which genes to use as negative
#'   controls for RUV.
#' @param negcon_eval a logical vector indicating which genes to use as 
#'   negative controls for evaluation.
#' @param poscon a logical vector indicating which genes to use as positive
#'   controls.
#'   
#' @export
#' 
#' @return A \code{\link{SconeExperiment}} object.
#'   
setMethod(
  f = "SconeExperiment",
  signature = signature("matrix"),
  definition = function(object, qc, bio, batch,
                        negcon_ruv=NULL, negcon_eval=negcon_ruv,
                        poscon=NULL, is_log=FALSE) {

    which_qc <- which_bio <- which_batch <- integer()
    which_negconruv <- which_negconeval <- which_poscon <- integer()

    coldata <- as.data.frame(matrix(nrow=NCOL(object), ncol=0))

    if(!missing(qc)) {
      coldata <- cbind(coldata, qc)
      which_qc <- seq_len(NCOL(coldata))
    }

    if(!missing(bio)) {
      coldata <- cbind(coldata, bio)
      which_bio <- NCOL(coldata)
    }

    if(!missing(batch)) {
      coldata <- cbind(coldata, batch)
      which_batch <- NCOL(coldata)
    }

    rowdata <- as.data.frame(matrix(nrow=NROW(object), ncol=0))

    if(!is.null(negcon_ruv)) {
      rowdata <- cbind(rowdata, negcon_ruv)
      which_negconruv <- NCOL(rowdata)
    }

    if(!is.null(negcon_eval)) {
      rowdata <- cbind(rowdata, negcon_eval)
      which_negconeval <- NCOL(rowdata)
    }

    if(!is.null(poscon)) {
      rowdata <- cbind(rowdata, poscon)
      which_poscon <- NCOL(rowdata)
    }

    se <- SummarizedExperiment(object, rowData=rowdata, colData=coldata)
    SconeExperiment(se,  which_qc, which_bio, which_batch,
                    which_negconruv, which_negconeval, which_poscon, is_log)
  }
)

