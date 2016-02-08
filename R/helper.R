#' Internal Function
#' 
#' This function is used internally in scone to parse the variables used to generate the design matrices.
#' 
#' @param pars character. A vector of parameters corresponding to a row of params.
#' @param bio factor. The biological factor of interest.
#' @param batch factor. The known batch effects.
#' @param ruv_factors list. A list containing the factors of unwanted variation.
#' @param qc matrix. The principal components of the QC metrics.
#' 
#' @return A list with the variables to be passed to make_design.
parse_row <- function(pars, bio, batch, ruv_factors, qc) {
  sc_name <- paste(pars[1:2], collapse="_")
  
  W <- out_bio <- out_batch <- NULL
  
  if(pars[3]!="no_uv") {
    parsed <- strsplit(as.character(pars[3]), "=")[[1]]
    if(grepl("ruv", parsed[1])) {
      W <- ruv_factors[[sc_name]][,1:as.numeric(parsed[2])]
    } else {
      W <- qc[,1:as.numeric(parsed[2])]
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

#' Function to make a design matrix
#' 
#' This function is useful to create a design matrix, when the covariates are two (possibly nested) factors
#' and one or more continuous variables.
#' 
#' @details If nested=TRUE a nested design is used, i.e., the batch variable is assumed to be nested within
#' the bio variable. Here, nested means that each batch is made of observations from only one level of bio,
#' while each level of bio may contain multiple batches.
#'  
#' @export
#'  
#' @param bio factor. The biological factor of interest.
#' @param batch factor. The known batch effects.
#' @param W numeric. Either a vector or matrix containing one or more continuous covariates (e.g. RUV factors).
#' @param nested logical. Whether or not to consider a nested design (see details).
#' 
#' @return The design matrix.
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
    
    return(model.matrix(as.formula(f), contrasts=list(bio=contr.sum, batch=mat)))    
  } else {
    return(model.matrix(as.formula(f)))    
  }
}

#' Function to perform linear batch effect correction
#' 
#' Given a matrix with log expression values and a design matrix, this function fits a linear model
#' and removes the effects of the batch factor as well as of the linear variables encoded in W.
#' 
#' @details The function assumes that the columns of the design matrix corresponding to the variable
#' for which expression needs to be adjusted, start with either the word "batch" or the letter "W" (case sensitive).
#' Any other covariate (including the intercept) is kept.
#'  
#' @importFrom limma lmFit
#' @export
#' 
#' @param log_expr matrix. The log gene expression (genes in row, samples in columns).
#' @param design_mat matrix. The design matrix (usually the result of make_design).
#' @param batch factor. A factor with the batch information.
#' @param weights matrix. A matrix of weights.
#' @return The corrected log gene expression.
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
  
  log_norm <- log_expr - batch_term  - uv_term
}
