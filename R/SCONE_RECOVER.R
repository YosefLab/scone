out_dir = "/data/yosef2/HIV/dat/CMV_HIV_150wk/SCONE_out_rerun/data"
flist = paste(out_dir,list.files(out_dir),sep = "/")

# Scoring
scores = NULL
snames = NULL
for(fil in flist){
  load(fil)
  nom = scone_out$method
  sc = scone_out$evaluation
  if(!is.na(sc)[1]){
    print(nom)
    snames = c(snames,nom)
    scores = rbind(scores,sc$scores)
  }else{
    # print(nom)
  }
}
rownames(scores) = snames
escores = t(t(scores)*c(-1,-1,1,1,1,-1,1))
colnames(escores) = paste(colnames(scores),c("_INV","")[1.5 + c(-1,-1,1,1,1,-1,1)/2],sep = "")

scores = t(na.omit(t(scores)))
escores = t(na.omit(t(escores)))

# Score Heatmaps
require(gplots)
heatmap.2(t(apply(scores ,2,rank)[order(-apply((apply(escores ,2,rank)),1,min)),][1:20,]),density.info = 'n',trace = 'n',key.xlab = "Score",key.title = NA,
          Colv = NA,col = colorRampPalette(c("purple","black","yellow"))(100),
          margins = c(15,10),cexRow = .75,cexCol = .5)
title(main = "Top 20 SCONE",cex.lab=0.5)

my_rank = (1 + max(rank(apply((apply(escores ,2,rank)),1,min))) - rank(apply((apply(escores ,2,rank)),1,min)))

names(which(my_rank == 1))
load("/data/yosef2/HIV/dat/CMV_HIV_150wk/SCONE_out_rerun/data/NOIMPUTE_NONE.Rdata")
is_zero = scone_out$exp_mat == 0
load("/data/yosef2/HIV/dat/CMV_HIV_150wk/SCONE_out_rerun/data/NOIMPUTE_FQ_NOWEIGHT_NOBIO_NOBATCH_HK_1.Rdata")
myma = scone_out$exp_mat
hk = read.table("~/YosefCode/packages/RNASeq/summary/EXAMPLE/reference_files/house_keeping_human_names.txt")$V1
myma[is_zero] = 0
fnr_out = estimateFNR(myma,bulk_model = T,is.expressed = rownames(myma) %in% as.character(hk) )
w = 0 + !is_zero
w[is_zero] = 1-fnr_out$Z[is_zero]
wpc = wPCA(log(scone_out$exp_mat + 1),w,filt = T)
plot(wpc$x[,1:2], col = cols)

rownames(wpc$rotation) = rownames(scone_out$exp_mat)[wpc$to.pass]
sort(-abs(wpc$rotation[,2]))[1:20]

load("/data/yosef2/HIV/dat/CMV_HIV_150wk/SCONE_out_rerun/evaluate_material.Rdata")
pchs = 16
cols = rainbow(length(levels(evaluate_material$batch)))[evaluate_material$batch]
plot(scone_out$evaluation$pc_val[,c(1,2)], xlab = "wPC1", ylab = "wPC2",col = cols,pch = pchs,main = "Top SCONE (Filtered): T Cells")
legend("bottomleft",legend = (levels(evaluate_material$batch)),pch = 16, col = rainbow(length(levels(evaluate_material$batch))), cex = .5)
outy = log(scone_out$exp_mat+1)
cols = rainbow(2)[scone_out$evaluation$clusters]
plot(scone_out$evaluation$pc_val[,c(1,3)], xlab = "wPC1", ylab = "wPC2",col = cols,pch = pchs,main = "Top SCONE (Filtered): T Cells")

write.table(outy,file = "/data/yosef2/HIV/dat/CMV_HIV_150wk/SCONE_out_rerun/NOIMPUTE_FQ_NOWEIGHT_NOBIO_NOBATCH_HK_1_log_counts.txt",quote = F,sep = "\t",col.names = T,row.names = T)
  evaluate_material$
evaluate_material$dim_eval
plot(evaluate_material$DE_scores[,c(1,2)], xlab = "wPC1", ylab = "wPC2",col = cols,pch = pchs,main = "Top SCONE (Filtered): T Cells")


my_rank = (1 + max(rank(apply((apply(escores ,2,rank)),1,min))) - rank(apply((apply(escores ,2,rank)),1,min)))
write.table(cbind(scores,my_rank),file = "/data/yosef2/HIV/dat/CMV_HIV_150wk/SCONE_out_rerun/scores.txt",quote = F,sep = "\t")


