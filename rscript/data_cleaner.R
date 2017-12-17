# miRNA and mRNA regulatory relationship study on ET samples.
require(miRLAB)
setwd("/Users/lukez/Documents/mirna_mrna/original_data/")
dataset=system.file("extdata", "ToyEMT.csv", package="miRLAB")
impdata=ImputeNormData(dataset, 0.1)
dataset = system.file("extdata", "EMT.csv", pacage = "miRLAB")
data = read.csv("EMT.csv", header = F)
ps=Pearson(dataset, cause=1:3, effect=4:18)

dataset=system.file("extdata", "ToyEMT.csv", package="miRLAB")
toydata=read.csv(dataset, 0.1)




# load the dat and reshape the data

mirna_data = read.csv("miRNA_compiled_data.csv", header = F, stringsAsFactors = F)

condition = mirna_data[3,-1]
colnames(condition) = NULL
rownames(condition) = NULL
mi_sample_names = mirna_data[1,-1]


reshaped_mirna_data = t(mirna_data[-c(1,2,3), -1])

rownames(reshaped_mirna_data) = NULL
colnames(reshaped_mirna_data) = NULL

reshaped_mirna_data = data.frame(cbind(as.character(mi_sample_names),t(condition[1,]),reshaped_mirna_data))
colnames(reshaped_mirna_data) = c("sample_name","condition", as.character(mirna_data[-c(1,2,3),1]))

mi_replicates = reshaped_mirna_data[grepl("_",reshaped_mirna_data$sample_name),]

non_replicates = reshaped_mirna_data[!(reshaped_mirna_data$sample_name %in% mi_replicates$sample_name),]
non_replicates[,1] = as.character(non_replicates[,1])
non_replicates[,-c(1,2)] = apply(non_replicates[,-c(1,2)],1, as.numeric)
non_replicates[,c(1,2)] = apply(non_replicates[,c(1,2)], 2, as.character)


mi_replicates[,1] = as.character(mi_replicates[,1])
mi_replicates[,3:941] = apply(mi_replicates[,3:941],1, as.numeric)

# use the mean value of the repcliates to represent the valud for this sample

rep_sample_name = unique(unlist(strsplit(mi_replicates[,1], '_'))[seq(1,30,by = 2)])

for (r_name in rep_sample_name) {
  
  rep = mi_replicates[grepl(r_name,mi_replicates[,1]),]
  mean_value = apply(rep[,-c(1,2)],2,mean)
  #print(mean_value)
  cond = rep[,2][1]
  print(cond)
  #print(cond)
  non_replicates = rbind(non_replicates, c(r_name,as.character(cond),mean_value))
  
}


reshaped_mirna_data = non_replicates

mrna_data = read.csv("Platelet_95samples_newest.csv", header = F,
                     stringsAsFactors = F)

ET_mrna_data = mrna_data[,mrna_data[3,-1] != 'RT']

matched_mrna_data = ET_mrna_data[,ET_mrna_data[1,-1] %in% reshaped_mirna_data[,1]]


write.csv(matched_mrna_data[-c(1,2),-1], 
          "matched_mrna_data.csv", 
          col.names = F, row.names = F,
          quote = F)
matched_dataset = system.file("extdata", "/Users/lukez/Documents/mirna_mrna/original_data/matched_mrna_data.csv", package="miRLAB")
imputed_matched_mrna_data = ImputeNormData("/Users/lukez/Documents/mirna_mrna/original_data/matched_mrna_data.csv", r = 0.2)

mr_sample_names = mrna_data[1,-1]
conditions = mrna_data[3,-1]

reshaped_mrna_data = t(apply(mrna_data[-c(1,2,3), -1],2, as.numeric))

reshaped_mrna_data = cbind(as.character(mr_sample_names),
                           t(conditions),
                           reshaped_mrna_data)


colnames(reshaped_mrna_data) = c("sample_name", "condition", as.character(mrna_data[-c(1,2,3), 1]))

# merge mirna and mrna data into one data frame
overlaped_sample = merge(reshaped_mirna_data, reshaped_mrna_data, by = "sample_name")

# The paired ET and control samples
ET_paired_samples = overlaped_sample[overlaped_sample$condition.x != 'RT',]

ET_paired_samples = ET_paired_samples[,-942]
ET_paired_samples[,-c(1,2)] = apply(ET_paired_samples[,-c(1,2)],2,as.numeric)

### handle missing value ####
### count the percentile of each mirna present in
## each group###

library(data.table)

countZero <- function(ve) {
  return(sum(as.numeric(na.omit(ve)) == 0))
}

countNA <- function(ve) {
  return(sum(is.na(ve)))
}

mirna_range = 3:941
mrna_range = 942:1447


############  filter out extreme values #####################
#data.log <- log2(data.filt);  #choose to use logged data or not
# No entry is filtred out...
#cutoff <- rowMeans(ET_paired_samples[,mirna_range], na.rm = TRUE) - 3*apply(ET_paired_samples[,mirna_range], 1, sd, na.rm=T);
#ET_paired_samples[,mirna_range][ET_paired_samples[,mirna_range] < min(cutoff)] <- NA;

#write.table(data.log,file = "test.csv",sep = ',');


##### Keep miRNA with !=A > 70% in at least one group ####
##### that is, filter out those have A > 30% in all three groups ####
mi_part = ET_paired_samples[,mirna_range]

countZero <- function(ve, label) {
  ET_count = round(sum(!is.na(ve[label == 'ET']) & ve[label == 'ET'] == 0) / sum(label == 'ET'),4)
  NO_count = round(sum(!is.na(ve[label == 'NO']) & ve[label == 'NO'] == 0) / sum(label == 'NO'),4)
  return(c(ET_count, NO_count))
}

missing_percentagle = apply(mi_part, 2, countZero, label = ET_paired_samples[,2])
missing_cutoff = 0.3
removed_mirna = apply(missing_percentagle,2,max) < 0.3
mi_part = mi_part[!removed_mirna] # nothing removed


############  quantile normalization ####################
mi.normalized <- normalizeBetweenArrays(t(mi_part),method = "quantile");


############ knn imputation ########################   
filt_fun<-function(x) round(sum(is.na(x))/length(x),digits=2);
filt_ind <- apply(mi.normalized,1,filt_fun);
mi.normalized <- mi.normalized[!(filt_ind > 0.5),];
imputed <- impute.knn(mi.normalized, k = 10);

# use the log2 scale
mirna_part <- log2(imputed$data); 




# preprocess mrna part
mrna_part = t(ET_paired_samples[,mrna_range])
mrna_name = colnames(ET_paired_samples)[mrna_range]
######   handle missing value at mrna part ##############
na_percentage = 0.5
na_count = apply(mrna_part, 1, countNA) / nrow(ET_paired_samples)



# use the log2 scale of the mrna expression value
mrna_part = log2(mrna_part[na_count < na_percentage,])
mrna_name = mrna_name[na_count < na_percentage]

mrna_part_normalized <- normalizeBetweenArrays(mrna_part,method = "quantile");
mrna_part_normalized_imputed <- impute.knn(mrna_part_normalized, k = 10);

mrna_part_normalized_imputed <- mrna_part_normalized_imputed$data;

################## Convert custom ids into gene names #################

combined_processed_data = cbind(mi_part, t(mrna_part_normalized_imputed))
colnames(combined_processed_data) = c(colnames(ET_paired_samples)[mirna_range], mrna_name)
rownames(combined_processed_data) = ET_paired_samples[,1]

combined_processed_data = cbind(ET_paired_samples[,c(1,2)],combined_processed_data)
write.csv(combined_processed_data[,-c(1,2)],
          "processed_mirna_mrna.csv",
          col.names = F,
          row.names = F,
          quote = F)
cause = 1:ncol(mi_part) + 2

effect = 1:length(mrna_name) + nrow(mi_part) + 2

hc = hclust(dist(combined_processed_data[,cause]),method = 'complete')
plot(hc,labels = combined_processed_data[,2])


#find differentially expressed 

library(miRLAB)
top_mirna = 200
top_mrna = 6000

padj_mirna = 0.5
padj_mrna = 0.05
write.csv(t(combined_processed_data[combined_processed_data[,2] == 'ET',cause]),
          "mirna_et.csv",
          col.names = T,
          row.names = T)
write.csv(t(combined_processed_data[combined_processed_data[,2] == 'NO',cause]),
          "mirna_NO.csv",
          col.names = T,
          row.names = T)

write.csv(t(combined_processed_data[combined_processed_data[,2] == 'ET',effect]),
          "mrna_et.csv",
          col.names = T,
          row.names = T)

write.csv(t(combined_processed_data[combined_processed_data[,2] == 'NO',effect]),
          "mrna_NO.csv",
          col.names = T,
          row.names = T)

# USE the miRLAB to find the differentially expressed microRNA
# NO microRNA is differentially expressed under the significant level 0.05
diff_expressed_mirlab = DiffExpAnalysis("mirna_et.csv",
                                 "mirna_no.csv",
                                 "mrna_et.csv",
                                 "mrna_NO.csv",
                                 top_mirna,
                                 top_mrna,
                                 padj_mirna,
                                 padj_mrna)

# USE the samr package to find the differentially expressed microRNA

library(samr)
mi_sam_data = list(x = t(as.matrix(combined_processed_data[,cause])),
                   y = as.numeric(factor(combined_processed_data[,2])),
                   logged2=FALSE,
                   geneid = as.character(length(cause)),
                   gene.names = colnames(combined_processed_data[,cause]))
samr.obj<-samr(mi_sam_data,resp.type="Two class unpaired", nperms=100)
delta.table <- samr.compute.delta.table(samr.obj, min.foldchange= 2)

siggenes.table <- samr.compute.siggenes.table(samr.obj, del=0, data, delta.table,all.genes=F)

down_regulated_mirna = combined_processed_data[,cause][,as.numeric(siggenes.table$genes.lo[,3])]
up_regulated_mirna = combined_processed_data[,cause][,as.numeric(siggenes.table$genes.up[,3])]
fc_2_down = down_regulated_mirna[,which(as.numeric(siggenes.table$genes.lo[,"Fold Change"]) <= 0.5)]
fc_2_up = up_regulated_mirna[,which(as.numeric(siggenes.table$genes.up[,"Fold Change"]) >= 2)]
genes_name_down_fc_2 = colnames(combined_processed_data[,cause])[which(as.numeric(siggenes.table$genes.lo[,"Fold Change"]) <= 0.5)]
genes_name_up_fc_2 = colnames(combined_processed_data[,cause])[which(as.numeric(siggenes.table$genes.up[,"Fold Change"]) >= 2)]
colnames(fc_2_down) = genes_name_down_fc_2
colnames(fc_2_up) = genes_name_up_fc_2
differentially_expressed_mirnas_SAM = cbind(fc_2_down, fc_2_up)

# Another way to call samr package

mi.cl <- as.numeric(factor(combined_processed_data[,2])) - 1;
mi.cl <- as.numeric(!mi.cl);
sam.out <- sam(t(combined_processed_data[,cause]), mi.cl, rand = 820, gene.names = colnames(combined_processed_data[,cause]));

differentially_expressed_mirnas_sam = combined_processed_data[,cause][sam.out@p.value < 0.95 ]

print(sam.out, seq(0.6, 0.9, 0.01));
# print(sam.out, seq(0.6, 0.65, 0.0005));

glist <- summary(sam.out,0.618);
sam.sig <- glist@mat.sig[glist@mat.sig$q.value < 0.05 & (glist@mat.sig$R.fold > 2 | glist@mat.sig$R.fold < (1/2)),];
sigs <- rownames(sam.sig);
# NO: 1; ET: 0 #
#mi.cl <- as.numeric(!mi.cl);
sam.out <- SAM(t(combined_processed_data[,cause]), 
                mi.cl,
                rand = 820, 
                gene.names = colnames(combined_processed_data[,cause]),
                resp.type = "Two class unpaired",
                testStatistic = "standard",
                regression.method = 'standard');
print(sam.out)
delta.table <- samr.compute.delta.table(sam.out, min.foldchange=0.1,nvals=200)

sig_genes = samr.compute.siggenes.table(sam.out)

# differential mirna from erya's script
sig_diff_mirna_erya = sigs

mirna_expression_erya = combined_processed_data[,colnames(combined_processed_data) %in% sig_diff_mirna_erya]
boxplot(mirna_expression_erya[,3] ~ combined_processed_data[,2])
whole_differentially_data = cbind(mirna_expression_erya, combined_processed_data[,effect])
write.csv(whole_differentially_data,
          "differential_mirna_mrna_erya.csv",
          col.names = T,
          quote = F)
#
# Association study

cause = 1:11
effect = 12:ncol(whole_differentially_data) 

pearson = Pearson("differential_mirna_mrna_erya.csv", cause, effect)
spearman = Spearman("differential_mirna_mrna_erya.csv", cause, effect)
kendall = Kendall("differential_mirna_mrna_erya.csv", cause, effect)
dcov = Dcov("differential_mirna_mrna_erya.csv", cause, effect)
mi = MI("differential_mirna_mrna_erya.csv", cause, effect)
hoeffding = Hoeffding("differential_mirna_mrna_erya.csv", cause, effect)

# 


