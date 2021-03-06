% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/helper.R
\name{get_scores}
\alias{get_scores}
\alias{get_scores,SconeExperiment-method}
\alias{get_score_ranks}
\alias{get_score_ranks,SconeExperiment-method}
\title{Extract scone scores}
\usage{
get_scores(x)

get_score_ranks(x)

\S4method{get_scores}{SconeExperiment}(x)

\S4method{get_score_ranks}{SconeExperiment}(x)
}
\arguments{
\item{x}{an object of class \code{\link{SconeExperiment}}.}
}
\value{
\code{get_scores} returns a matrix with all (non-missing) scone 
  scores, ordered by average score rank.

\code{get_score_ranks} returns a vector of average score ranks.
}
\description{
Extract scone scores
}
\examples{
set.seed(42)
mat <- matrix(rpois(500, lambda = 5), ncol=10)
colnames(mat) <- paste("X", 1:ncol(mat), sep="")
obj <- SconeExperiment(mat)
res <- scone(obj, scaling=list(none=identity, uq=UQ_FN),
           evaluate=TRUE, k_ruv=0, k_qc=0,
           eval_kclust=2, bpparam = BiocParallel::SerialParam())
scores = get_scores(res)
score_ranks = get_score_ranks(res)

}
