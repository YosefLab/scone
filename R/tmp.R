library(brainUtils)
library(EDASeq)
library(clusterCells)
library(matrixStats)

setMethod(
  f = "plotRLE",
  signature = signature(x="matrix"),
  definition = function(x,...) {
    y <- log(x+1)
    median <- apply(y, 1, median)
    rle <- apply(y, 2, function(x) x - median)
    
    boxplot(rle, ...)
    abline(h=0, lty=2)
    invisible(rle)
  }
)

load("~/git/brainAnalysis/brain/data/Nov15/SCONE/NOIMPUTE_FQ_NOWEIGHT_BIO_BATCH_HK_1.Rdata")
norm <- scone_out$exp_mat
norm[norm<0] <- 0
log_norm <- log1p(norm)

qc <- read.table("~/git/brainAnalysis/brain/data/Nov15/SCONE/filtered_qc.txt")

type <- read.table("~/git/brainAnalysis/brain/data/Nov15/SCONE/cell_type.txt", stringsAsFactors = FALSE)[,2]
type <- as.factor(type)

batch <- read.table("~/git/brainAnalysis/brain/data/Nov15/SCONE/batch.txt", stringsAsFactors = FALSE)[,2]
batch <- factor(batch, levels=c(unlist(tapply(batch, type, function(x) names(table(x))))))

names(batch) <- names(type) <- colnames(norm)
batch <- sort(batch)
type <- type[names(batch)]
norm <- norm[,names(batch)]
qc <- qc[names(batch),]

data(housekeeping)
negcon <- intersect(housekeeping, rownames(norm))

markers <- read.table("~/git/brainAnalysis/brain/data/macklis_markers.txt", as.is=TRUE)
poscon <- intersect(markers[,1], rownames(norm))

pca <- prcomp(t(log1p(norm)), center=TRUE, scale=TRUE)
pca_pc <- prcomp(t(log1p(norm[poscon,])), center=TRUE, scale=TRUE)
pca_nc <- prcomp(t(log1p(norm[negcon,])), center=TRUE, scale=TRUE)

vars <- rowVars(log_norm)
names(vars) <- rownames(log_norm)
mv_var <- names(sort(vars)[1:100])

cvs <- rowSds(log_norm)/rowMeans(log_norm)
mv_cv <- names(sort(cvs)[1:100])

pca_cv <- prcomp(t(log1p(norm[mv_cv,])), center=TRUE, scale=TRUE)
pca_var <- prcomp(t(log1p(norm[mv_var,])), center=TRUE, scale=TRUE)

cors <- sapply(1:10, function(i) cor(pca$x[,i], qc))
colnames(cors) <- paste("PC", 1:NCOL(cors), sep="")
rownames(cors) <- colnames(qc)

choose_color <- function(input) {
  if(input$color_code=="type") {
    cols <- cc[type]
  } else {
    cols <- cc[batch]
  }
  return(cols)
}

rle <- plotRLE(norm, outline=FALSE)

sub <- lapply(2:10, function(kk) subsampleClustering(pca$x[,1:50], k=kk))

pam <- lapply(2:10, function(kk) clusterAll(pca$x[,1:50], subsample=FALSE, 
                  sequential=FALSE, clusterFunction = "pam",
                  clusterDArgs = list(k=kk)))

qc_pca <- prcomp(qc, center=TRUE, scale=TRUE)

source("~/git/YosefCode/packages/RNASeq/OFBIT/SCONE/estimateFNR.R")
data(housekeeping)
hk <- intersect(housekeeping, rownames(norm))
fnr <- estimateFNR(norm, bulk_model=TRUE, is.expressed=which(rownames(norm) %in% hk))
w <- (norm > 0) + (1 - fnr$Z)

logistic <- binomial()$linkinv

fitted <- t(fnr$W %*% fnr$Alpha)
save.image("tmp.rda")