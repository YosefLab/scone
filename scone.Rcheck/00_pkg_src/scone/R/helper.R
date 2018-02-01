#' @rdname get_params
#'   
#' @param x an object of class \code{\link{SconeExperiment}}.
#'   
#' @return A data.frame containing workflow parameters for each scone workflow.
#'   
#' @export
#' 
setMethod(
  f = "get_params",
  signature = signature(x = "SconeExperiment"),
  definition = function(x) {
    return(x@scone_params)
  })

#' @rdname get_scores
#'   
#' @param x an object of class \code{\link{SconeExperiment}}.
#'   
#' @return \code{get_scores} returns a matrix with all (non-missing) scone 
#'   scores, ordered by average score rank.
#'   
#' @export
#' 
setMethod(
  f = "get_scores",
  signature = signature(x = "SconeExperiment"),
  definition = function(x) {
      scores <- t(na.omit(t(x@scone_scores[,-NCOL(x@scone_scores)])))
      return(scores)
})

#' @rdname get_scores
#'   
#' @return \code{get_score_ranks} returns a vector of average score ranks.
#'   
#' @export
#' 
setMethod(
  f = "get_score_ranks",
  signature = signature(x = "SconeExperiment"),
  definition = function(x) {
    return(x@scone_scores[,NCOL(x@scone_scores)])
  })


#' @rdname get_negconruv
#'   
#' @param x an object of class \code{\link{SconeExperiment}}.
#'   
#' @return NULL or a logical vector.
#'   
#' @return For \code{get_negconruv} the returned vector indicates which genes
#'   are negative controls to be used for RUV.
#'   
#' @export
#' 
setMethod(
  f = "get_negconruv",
  signature = signature(x = "SconeExperiment"),
  definition = function(x) {
    if(length(x@which_negconruv) == 0) {
      return(NULL)
    } else {
      return(rowData(x)[,x@which_negconruv])
    }
  }
)

#' @rdname get_negconruv
#'   
#' @return For \code{get_negconeval} the returned vector indicates which genes
#'   are negative controls to be used for evaluation.
#'
#' @export
#' 
setMethod(
  f = "get_negconeval",
  signature = signature(x = "SconeExperiment"),
  definition = function(x) {
    if(length(x@which_negconeval) == 0) {
      return(NULL)
    } else {
      return(rowData(x)[,x@which_negconeval])
    }
  }
)

#' @rdname get_negconruv
#'   
#' @return For \code{get_poscon} the returned vector indicates which genes are 
#'   positive controls to be used for evaluation.
#'
#' @export
#'
setMethod(
  f = "get_poscon",
  signature = signature(x = "SconeExperiment"),
  definition = function(x) {
    if(length(x@which_poscon) == 0) {
      return(NULL)
    } else {
      return(rowData(x)[,x@which_poscon])
    }
  }
)

#' @rdname get_qc
#'   
#' @param x an object of class \code{\link{SconeExperiment}}.
#'   
#' @return NULL or the quality control (QC) metric matrix.
#'
#' @export
#'
setMethod(
  f = "get_qc",
  signature = signature(x = "SconeExperiment"),
  definition = function(x) {
    if(length(x@which_qc) == 0) {
      return(NULL)
    } else {
      retval <- as.matrix(colData(x)[, x@which_qc, drop=FALSE])
      return(retval)
    }
  }
)

#' @rdname get_bio
#'   
#' @param x an object of class \code{\link{SconeExperiment}}.
#'   
#' @return NULL or a factor containing bio or batch covariate.
#'
#' @export
#'
setMethod(
  f = "get_bio",
  signature = signature(x = "SconeExperiment"),
  definition = function(x) {
    if(length(x@which_bio) == 0) {
      return(NULL)
    } else {
      return(colData(x)[,x@which_bio])
    }
  }
)

#' @rdname get_bio
#'
#' @export
#'
setMethod(
  f = "get_batch",
  signature = signature(x = "SconeExperiment"),
  definition = function(x) {
    if(length(x@which_batch) == 0) {
      return(NULL)
    } else {
      return(colData(x)[,x@which_batch])
    }
  }
)

#' Parse rows
#' 
#' This function is used internally in scone to parse the variables used to
#' generate the design matrices.
#' 
#' @param pars character. A vector of parameters corresponding to a row of
#'   workflow parameters.
#' @param bio factor. The biological covariate.
#' @param batch factor. The batch covariate.
#' @param ruv_factors list. A list containing the factors of unwanted variation
#'   (RUVg) for all upstream workflows.
#' @param qc matrix. The principal components of the QC metric matrix.
#'   
#' @return A list with the variables to be passed to make_design.
#'   
#' @keywords internal
#'   
.parse_row <- function(pars, bio, batch, ruv_factors, qc) {
  
  # Define upstream workflow: imputation x scaling
  sc_name <- paste(pars[1:2], collapse="_")

  W <- out_bio <- out_batch <- NULL
  
  if(pars[3]!="no_uv") {
    parsed <- strsplit(as.character(pars[3]), "=")[[1]]
    if(grepl("ruv", parsed[1])) {
      W <- ruv_factors[[sc_name]][,seq_len(as.numeric(parsed[2]))]
    } else {
      W <- qc[,seq_len(as.numeric(parsed[2]))]
    }
  }

  if(pars[4]=="bio") {
    out_bio <- bio
  }

  if(pars[5]=="batch") {
    out_batch <- batch
  }

  return(list(sc_name=sc_name, W=W, bio=out_bio, batch=out_batch))
}

#' Make a Design Matrix
#' 
#' This function builds a design matrix for the Adjustment Normalization Step,
#' in which covariates are two (possibly nested) categorical factors and one or
#' more continuous variables.
#' 
#' @details If nested=TRUE a nested design is used, i.e. the batch variable is
#'   assumed to be nested within the bio variable. Here, nested means that each
#'   batch is composed of samples from only *one* level of bio, while each
#'   level of bio may contain multiple batches.
#'   
#' @export
#' 
#' @param bio factor. The biological covariate.
#' @param batch factor. The batch covariate.
#' @param W numeric. Either a vector or matrix containing one or more
#'   continuous covariates (e.g. RUVg factors).
#' @param nested logical. Whether or not to consider a nested design
#'   (see details).
#'   
#' @return The design matrix.
#'   
#' @examples 
#' 
#' bio = as.factor(rep(c(1,2),each = 2))
#' batch = as.factor(rep(c(1,2),2))
#' design_mat = make_design(bio,batch, W = NULL)
#' 
make_design <- function(bio, batch, W, nested=FALSE) {
  
  if(nested & (is.null(bio) | is.null(batch))) {
    stop("Nested design can be used only if both batch and bio are specified.")
  }
  
  if(!is.null(bio)) {
    if(class(bio)!="factor") {
      stop("bio must be a factor.")
    }
  }
  
  if(!is.null(batch)){
    if(class(batch)!="factor") {
      stop("batch must be a factor.")
    }
  }

  f <- "~ 1"
  if(!is.null(bio)) {
    f <- paste(f, "bio", sep="+")
  }
  if(!is.null(batch)) {
    f <- paste(f, "batch", sep="+")
  }
  if(!is.null(W)) {
    f <- paste(f, "W", sep="+")
  }

  if(is.null(bio) & is.null(batch) & is.null(W)) {
    return(NULL)
  } else if (!is.null(bio) & !is.null(batch) & nested) {
    
    n_vec <- tapply(batch, bio, function(x) nlevels(droplevels(x)))

    mat = matrix(0,nrow = sum(n_vec),ncol = sum(n_vec - 1))
    xi = 1
    yi = 1
    for(i in 1:length(n_vec)){
      if(n_vec[i] > 1){
        cs = contr.sum(n_vec[i])
        dd = dim(cs)
        mat[xi:(xi + dd[1] - 1),yi:(yi + dd[2] - 1)] = cs
        xi = xi + dd[1]
        yi = yi + dd[2]
      }else{
        xi = xi + 1
      }
    }

    return(model.matrix(as.formula(f), 
                        contrasts=list(bio=contr.sum, batch=mat)))
  } else {
    return(model.matrix(as.formula(f)))
  }
}

#' Linear Adjustment Normalization
#' 
#' Given a matrix with log expression values and a design matrix, this function
#' fits a linear model and removes the effects of the batch factor
#' as well as of the linear variables encoded in W.
#' 
#' @details The function assumes that the columns of the design matrix
#'   corresponding to the variable for which expression needs to be adjusted,
#'   start with either the word "batch" or the letter "W" (case sensitive). Any
#'   other covariate (including the intercept) is kept.
#'   
#' @importFrom limma lmFit
#' @export
#' 
#' @param log_expr matrix. The log gene expression (genes in row, samples in
#'   columns).
#' @param design_mat matrix. The design matrix (usually the result of
#'   make_design).
#' @param batch factor. A factor with the batch information, identifying batch
#'   effect to be removed.
#' @param weights matrix. A matrix of weights.
#' @return The corrected log gene expression.
#'   
#' @examples
#' 
#' set.seed(141)
#' bio = as.factor(rep(c(1,2),each = 2))
#' batch = as.factor(rep(c(1,2),2))
#' design_mat = make_design(bio,batch, W = NULL)
#' 
#' log_expr = matrix(rnorm(20),ncol = 4)
#' adjusted_log_expr = lm_adjust(log_expr = log_expr,
#'   design_mat = design_mat,
#'   batch = batch)
#' 
lm_adjust <- function(log_expr, design_mat, batch=NULL, weights=NULL) {
  lm_object <- lmFit(log_expr, design = design_mat, weights = weights)

  uvind <- grep("^W", colnames(design_mat))
  bind <- grep("^batch", colnames(design_mat))

  if(length(uvind)) {
    uv_term <- t(design_mat[,uvind] %*% t(lm_object$coefficients[,uvind]))
  } else {
    uv_term <- 0
  }

  if(length(bind)) {
    if(is.character(attr(design_mat,"contrasts")$batch)) {
      contr <- get(attr(design_mat,"contrasts")$batch)(nlevels(batch))
    } else {
      contr <- attr(design_mat,"contrasts")$batch
    }
    batch_term <- t(contr %*% t(lm_object$coefficients[,bind]))[,batch]
  } else {
    batch_term <- 0
  }

  return(log_expr - batch_term  - uv_term)
}
