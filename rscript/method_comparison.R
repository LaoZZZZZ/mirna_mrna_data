# Methods comparison
setwd("/Users/lukez/Documents/mirna_mrna/original_data/")
require(miRLAB)
library(PMA);

all_mirna_et_no = read.csv("processed_mirna.csv", header = T, stringsAsFactors = F)
all_mrna_et_no = read.csv("processed_mrna.csv", header = T, stringsAsFactors = F)

mi_condition = all_mirna_et_no[1,-1]
mr_condition = all_mrna_et_no[1,-1]

mi_data = all_mirna_et_no[-1,-1]
rownames(mi_data) = all_mirna_et_no[-1,1]
colnames(mi_data) = colnames(all_mirna_et_no)[-1]
write.csv(mi_data[,mi_condition == 'ET'],
          "mirna_et.csv",
          col.names = T,
          row.names = T,
          quote = F)

write.csv(mi_data[,mi_condition == 'NO'],
          "mirna_no.csv",
          col.names = T,
          row.names = T,
          quote = F)


mr_data = all_mrna_et_no[-1,-1]
rownames(mr_data) =all_mrna_et_no[-1,1]
colnames(mr_data) = colnames(all_mrna_et_no)[-1]
write.csv(mr_data[,mr_condition == 'ET'],
          "mrna_et.csv",
          col.names = T,
          row.names = T,
          quote = F)

write.csv(mr_data[,mr_condition == 'NO'],
          "mrna_no.csv",
          col.names = T,
          row.names = T,
          quote = F)


top_mirna = 200
top_mrna = 400
padj_mirna = 0.01
padj_mrna = 0.01
differientially_expressed_limma = DiffExpAnalysis("mirna_et.csv",
                                                  "mirna_no.csv",
                                                  "mrna_et.csv",
                                                  "mrna_no.csv",
                                                  top_mirna,
                                                  top_mrna,
                                                  padj_mirna,
                                                  padj_mrna)


# all genes in the SCCA result belongs to the differentially expressed mrna subset.
erya_genes = c("HSD17B12", "GLA", "MMP1", "PKIG" , "SERPINI1", "CAV2", "WASF1",
               "NME4", "TIMP1", "TGFB1I1")



sum(erya_genes %in% colnames(differientially_expressed_limma))

sum(erya_genes %in% colnames(mrna.all)[out.un$v != 0])
write.csv(differientially_expressed_limma,
          "differentially_expressed_mirna_mrna.csv",
          col.names = F,
          row.names = F,
          quote = F)

write.csv(differientially_expressed_limma[,1:61],
          "differentially_expressed_mirna.csv",
          col.names = T,
          row.names = T, 
          quote = F)

write.csv(differientially_expressed_limma[,62:80],
          "differentially_expressed_mrna.csv",
          col.names = T,
          row.names = T, 
          quote = F)
processed_data = cbind(mirna.all,mrna.all)


write.csv(processed_data,
          "differential_mirna_mrna_erya.csv",
          col.names = T,
          row.names = F,
          quote = F)



# Sparse canonical partial correlation.
# This function takes 
scca <- function(diff_expressed_mirna_file,
                 diff_expressed_mrna_file,
                 k = 1) {
  mirna_data = read.csv(diff_expressed_mirna_file, header = T,
                        stringsAsFactors = F)
  mirna.all = apply(mirna_data[,-1],2,as.numeric)
  rownames(mirna.all) = mirna_data[,1]
  
  
  mrna_data = read.csv(diff_expressed_mrna_file, header = T,
                       stringsAsFactors = F)
  mrna.all = apply(mrna_data[,-1],2,as.numeric)
  rownames(mrna.all) = mrna_data[,1]
  
  out.un = CCA(x=mirna.all, z=mrna.all, typex="standard", typez="standard", K=k, 
               xnames = colnames(mirna.all), znames = colnames(mrna.all),
               penaltyx = 0.2, penaltyz = 0.4)
  #print(out.un, verbose=T)
  return(out.un)
}

fillValue <- function(u_array, v_array, mat) {
    
  index = as.matrix(expand.grid(which(u_array != 0), which(v_array !=0)), ncol= 2)
  for (i in 1:nrow(index)) {
    mat[index[i,1],index[i,2]] = 1
  }
}

buildSCCAMatrix <- function(scca_result) {
  num_mirna = nrow(scca_result$u)
  num_components = ncol(scca_result$u)
  
  num_mrna = dim(scca_result$v)[1]
  # Invalid input
  if (ncol(scca_result$v) != num_components) {
    print("Invalid input!\n")
    sys.on.exit()
  }
  result = matrix(0, nrow = num_mirna, ncol = num_mrna)
  for (i in 1:num_components) {
    
    index = as.matrix(expand.grid(which(scca_result$u[,i] != 0), which(scca_result$v[,i] !=0)), ncol= 2)
    for (i in 1:nrow(index)) {
      result[index[i,1],index[i,2]] = 1
    }
    #fillValue(scca_result$u[,i], scca_result$v[,i], result)
  }
  result
}



dev.off()

cause = 1:61
effect = 62:80


input_file = 'differentially_expressed_mirna_mrna.csv'
# Association study
pearson = Pearson(input_file, cause, effect)
spearman = Spearman(input_file, cause, effect)
kendall = Kendall(input_file, cause, effect)
dcov = Dcov(input_file, cause, effect)
mi = MI(input_file, cause, effect)
hoeffding = Hoeffding(input_file, cause, effect)
rdc = RDC(input_file, cause, effect)

# regression study
# Lasso 
lasso = Lasso(input_file, cause, effect)
# Elastic-net regression 
elastic = Elastic(input_file, cause, effect)


# Bayesian network

ida = IDA(input_file, cause, effect)

# Other method
zscore = Zscore(input_file, cause, effect)
promise = ProMISe(input_file, cause, effect)


all = list(pearson,spearman, kendall,dcov, hoeffding, mi, ida, rdc, lasso, elastic,
           zscore, promise)

allresultsTop100 = experiment(all, topk = 20, "ground_truth.csv",2.0)

ValidateAll(ida, topk = 100, "Groundtruth", 1.0)


### plot the heatmap in the same order

reorderAndPlot <- function(data, rowIndex, colIndex, percent = 0.95, ...) {
  

  data = data[rowIndex,colIndex]
  data = gateMatrix(data,percent)
  
  heatmap.3(data, dendrogram = 'none', 
            Rowv = F, Colv = F,...)
}


sca_result = scca("differentially_expressed_mirna.csv",
                  "differentially_expressed_mrna.csv",
                  3)
mat = t(buildSCCAMatrix(sca_result))
colnames(mat) = colnames(differientially_expressed_limma)[cause]
rownames(mat) = colnames(differientially_expressed_limma)[effect]

dev.off()
cluster_index = heatmap.2(mat, main = 'scca')

pdf("Comparison between SCCA and other methods.pdf")

reorderAndPlot(mat, cluster_index$rowInd, cluster_index$colInd,
               main = 'Sparse Canonical correlation Analysis')

reorderAndPlot(pearson, cluster_index$rowInd, cluster_index$colInd,
               main = "Pearson Correlation")

reorderAndPlot(spearman, cluster_index$rowInd, cluster_index$colInd,
               main = "Spearman Correlation")

reorderAndPlot(kendall, cluster_index$rowInd, cluster_index$colInd,
               main = "Kendall Rank Correlation")

reorderAndPlot(hoeffding, cluster_index$rowInd, cluster_index$colInd,
               main = "Hoeffding's D measure")


reorderAndPlot(dcov, cluster_index$rowInd, cluster_index$colInd,
               main = "Distance correlation")


reorderAndPlot(rdc, cluster_index$rowInd, cluster_index$colInd,
               main = "Randomised Dependence Coefficient")


reorderAndPlot(lasso, cluster_index$rowInd, cluster_index$colInd,
               main = "Lasso regression")


reorderAndPlot(elastic, cluster_index$rowInd, cluster_index$colInd,
               main = "Elastic-net regression")



reorderAndPlot(ida, cluster_index$rowInd, cluster_index$colInd,
               main = "Intervention calculus when the DAG is Absent")


reorderAndPlot(zscore, cluster_index$rowInd, cluster_index$colInd,
               main = "Z-score")


dev.off()


# Compare result between methods other than SCCA

## Take top 5-10 percent of values in each result to compare with 

getPercentage <- function(data, percent) {
  d = quantile(abs(as.numeric(data)), probs = percent)
  return(as.numeric(d))
}

gateMatrix <- function(data, percentage) {
  
  threshold = getPercentage(data, percentage)
  cp = data
  cp[abs(data) < threshold] = 0
  return(cp)
}


library(GMD)

pdf("ida_result.pdf", width = 10, height = 10)
heatmap.3(ida)
dev.off()


pdf("rdc_result.pdf", width = 10, height = 10, dendrogram = "none")
heatmap.3(rdc)
dev.off()

pdf("kendall_result.pdf", width = 10, height = 10)
heatmap.3(kendall)
dev.off()

pdf("pearson_result.pdf", width = 10, height = 10)
heatmap.3(pearson)
dev.off()

pdf("lasso_result.pdf", width = 10, height = 10)
heatmap.3(lasso)
dev.off()


pdf("dcov_result.pdf", width = 10, height = 10)
heatmap.3(dcov)
dev.off()


pdf("hoeffding_result.pdf", width = 10, height = 10)
heatmap.3(hoeffding)
dev.off()



pdf("zscore_result.pdf", width = 10, height = 10)
heatmap.3(zscore)
dev.off()


pdf("promise_result.pdf", width = 10, height = 10)
heatmap.3(promise)
dev.off()


pdf("mi_result.pdf", width = 10, height = 10)
heatmap.3(mi)
dev.off()
# save truncated ida heat map.
pdf("ida_90_result.pdf")
d = gateMatrix(ida,0.90)
heatmap.3(d)
dev.off()


pdf("kendal_95_result.pdf")
d = gateMatrix(kendall,0.95)
heatmap.3(d)
dev.off()


pdf("lasso_95_result.pdf")
d = gateMatrix(lasso, 0.95)
heatmap.3(d)
dev.off()

pdf("spearman_95_result.pdf")
d = gateMatrix(spearman, 0.95)
heatmap.3(d)
dev.off()


pdf("mi_95_result.pdf")
d = gateMatrix(mi, 0.95)
heatmap.3(d)
dev.off()

pdf("rdc_95_result.pdf")
d = gateMatrix(rdc, 0.95)
heatmap.3(d)
dev.off()


pdf("zscore_95_result.pdf")
d = gateMatrix(zscore, 0.95)
heatmap.3(d)
dev.off()



dev.off()
cluster_index = heatmap.2(gateMatrix(ida,.95), main = 'ida')

pdf("Comparison between IDA and other methods.pdf")

reorderAndPlot(ida, cluster_index$rowInd, cluster_index$colInd,
               main = "Intervention calculus when the DAG is Absent")

reorderAndPlot(pearson, cluster_index$rowInd, cluster_index$colInd,
               main = "Pearson Correlation")

reorderAndPlot(spearman, cluster_index$rowInd, cluster_index$colInd,
               main = "Spearman Correlation")

reorderAndPlot(kendall, cluster_index$rowInd, cluster_index$colInd,
               main = "Kendall Rank Correlation")

reorderAndPlot(hoeffding, cluster_index$rowInd, cluster_index$colInd,
               main = "Hoeffding's D measure")


reorderAndPlot(dcov, cluster_index$rowInd, cluster_index$colInd,
               main = "Distance correlation")


reorderAndPlot(rdc, cluster_index$rowInd, cluster_index$colInd,
               main = "Randomised Dependence Coefficient")


reorderAndPlot(lasso, cluster_index$rowInd, cluster_index$colInd,
               main = "Lasso regression")


reorderAndPlot(elastic, cluster_index$rowInd, cluster_index$colInd,
               main = "Elastic-net regression")



reorderAndPlot(mat, cluster_index$rowInd, cluster_index$colInd,
               main ='Sparse Canonical correlation Analysis' )


reorderAndPlot(zscore, cluster_index$rowInd, cluster_index$colInd,
               main = "Z-score")


dev.off()

