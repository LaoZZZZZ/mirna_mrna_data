###################### Data Process ########################

##### miRNA: Agilent, human (after log2); PMA: p-present; m-marginal~p; a-absent; ####
##### mRNA: custom array (before log2) and affymetrix ####

setwd('/Users/lukez/Documents/mirna_mrna/original_data/');
library(siggenes);
library(limma);
library(impute);
library(multtest);
library(gplots);

############  Read in microRNA samples ###############
num1 <- read.table('miRNA_compiled_data.csv', sep = ",", na.strings ="",header=T,stringsAsFactors=F);
num2 <- read.table('miRNA_compiled_PMA.csv', sep = ",", na.strings ="",header=T,stringsAsFactors=F);


# combined duplicates

dat <- as.matrix(num1[-c(1,2),-1]);
data <- apply(dat,2,as.numeric);
# figure out the condition
unique_samples = gsub('_1','',colnames(num1)[-1][!grepl('_2',colnames(num1)[-1])])
mi.conditions = t(data.frame(as.character(num1[2,-1][!grepl('_2',colnames(num1)[-1])])))
colnames(mi.conditions) = unique_samples


#mi.group <- as.character(colnames(num1)[-1][grepl('_2',colnames(num1)[-1])]);

mi.names <- as.character(num1[-c(1,2),1]);
samples <- as.character(colnames(num1)[-1]);
colnames(data) <- samples;
rownames(data) <- mi.names;

pma.calls <- num2[-c(1,2),-1];
colnames(pma.calls) <- samples;
rownames(pma.calls) <- mi.names;

### handle technical replicates ######
combinePMA <- function(data) {
  if ('P' %in% data) {
    return('P')
  } else {
    return('A')
  }
}
combinePMASamples <- function(samples_rep) {
  return(unlist(apply(samples_rep,1,combinePMA)))
  
}
combineReplicates <- function(data, func) {
  
  replicate_samples = data[,grepl('_',colnames(data))]
  non_replicates = data[,!(colnames(data) %in% colnames(replicate_samples))]
  print(colnames(replicate_samples))
  rep_sample_name = unique(unlist(strsplit(colnames(replicate_samples),'_')))
  for (r_name in rep_sample_name) {
    if (grepl('THR', r_name)) {
      
      rep = replicate_samples[,grepl(r_name,colnames(replicate_samples))]
      print(c(r_name,colnames(rep)))
      rep_values = unlist(apply(rep,1,func))
      original_colname = colnames(non_replicates)
      non_replicates = cbind(non_replicates, rep_values)
      colnames(non_replicates) = c(original_colname, r_name)
    }
  }
  return(non_replicates)
}

dedup_data = combineReplicates(data, mean)
dedup_pma.calls = combineReplicates(pma.calls,combinePMA)
mi.group = as.vector(mi.conditions[,colnames(dedup_data)])

data = dedup_data
pma.calls = dedup_pma.calls
############  To filter the data  #################

##### Keep miRNA with !=A > 70% in at least one group ####
##### that is, filter out those have A > 30% in all three groups ####

filt_fun<-function(weight, group){
  grpid <- unique(group);
  buf <- list();
  index <- c(0,0);
  for (k in 1:length(grpid)) {
    buf[[k]] <- weight[group == grpid[k]];
    if (sum(buf[[k]]!='A') > 0.7 * length(buf[[k]])) {
      index[1] <- index[1] + 1;}
    else if (sum(buf[[k]]=='A') > 0.7 * length(buf[[k]])) {
      index[2] <- index[2] + 1;}
  }
  return(index); 
}
buffer <- apply(pma.calls,1,filt_fun,group = mi.group);
filt_ind <- t(buffer);
data.pass <- data[filt_ind[,1] > 0,];
rm(buffer,dat,num1,num2);

############  To find miRNAs only appear on particular one or two groups ###########
# data.interest <- data[(filt_ind[,1] > 0) && (filt_ind[,2] > 0)];
# print(paste('There are',length(data.interest), 'parculiar miRNAs identified'));

############  count absent value percentiles in each sample ####################  
count = rep(-1,dim(data.pass)[2]);
for (i in 1:dim(data.pass)[2]) {
  count[i] = sum(data.pass[,i] == 0)/dim(data.pass)[1];
  count[i] = round(count[i],4) * 100;
};
percentile <- cbind(colnames(data.pass),count);

############  Filter out bad samples ####################
data.filt <- data.pass[,as.numeric(percentile[,2]) < 40];
group.filt <- mi.group[as.numeric(percentile[,2]) < 40];
samples <- colnames(data.filt)
# THR 127 is deleted too (as well as THR30, THR31, THR528, THR126, THR536, THR011, THR037) 

############  Deal with THR118 duplicates ##########################  
# THR118 <- rowMeans(cbind(data.filt[,samples == "THR118_1"],data.filt[,samples == "THR118_2"]));
# buffer <- data.filt[,!(samples == "THR118_1") & (!(samples == "THR118_2"))];
# mi.group <- group.filt[!(samples == "THR118_1") & (!(samples == "THR118_2"))];
# mi.group <- c(mi.group,"ET");
# data.filt <- cbind(buffer,THR118);
# samples <- colnames(data.filt);

############  samples names ####################
# data.filt[data.filt == 0] <- NA;
# colnames(data.filt) <- gsub('_2',replacement = "", colnames(data.filt));
# colnames(data.filt) <- gsub('R0',replacement = 'R',colnames(data.filt));
# colnames(data.filt) <- gsub('NO0',replacement = 'NO',colnames(data.filt));

############  filter out extreme values #####################
#data.log <- log2(data.filt);  #choose to use logged data or not
cutoff <- colMeans(data.filt, na.rm = TRUE) - 3*apply(data.filt, 2, sd, na.rm=T);
data.filt[data.filt < min(cutoff)] <- NA;
which(data.filt < min(cutoff))
#write.table(data.log,file = "test.csv",sep = ',');

############ quantile normalization ##################
mi.normalized <- normalizeBetweenArrays(data.filt,method = "quantile");
rownames(mi.normalized) <- rownames(data.filt);
colnames(mi.normalized) <- colnames(data.filt);

############ knn imputation ########################   
filt_fun<-function(x) round(sum(is.na(x))/length(x),digits=2);
filt_ind <- apply(mi.normalized,1,filt_fun);
mi.normalized <- mi.normalized[!(filt_ind > 0.5),];
### hsa-miR-135a* was deleted #### 
imputed <- impute.knn(mi.normalized, k = 10);
midata.all <- imputed$data; 
mi.samples <- colnames(midata.all);
mi.names <- rownames(midata.all);
mi.group <- mi.group

##########################################################################################
############ Read in custom array mRNA data ##################################################
##########################################################################################
custom.data <- read.csv("Platelet_95samples_newest.csv",sep=",",
                          na.strings ="",header = T,stringsAsFactors = F);
custom.dat <- data.matrix(custom.data[-c(1,2),-1]);
custom.group <- as.character(custom.data[2,-1]);
custom.gnames <- as.character(custom.data[-c(1,2),1]);
custom.samples <- colnames(custom.data)[-1];
row.names(custom.dat) <- custom.gnames;
#colnames(custom.dat) <- gsub('R0',replacement = 'R',colnames(custom.dat));
#colnames(custom.dat) <- gsub('NO0',replacement = 'NO',colnames(custom.dat));

############ Filtering by 50% missing ############ 
filt_fun<-function(x) round(sum(is.na(x))/length(x),digits=2);
ids <- apply(custom.dat,1,filt_fun);
# use the original scale
custom.gdata <- log2(custom.dat[ids < 0.5,]);
#custom.gdata <- custom.dat[ids < 0.5,];

custom.gnames <- as.character(custom.gnames[ids < 0.5]);
rm(custom.data,ids);

############ quantile normalization ##################
custom.normalized <- normalizeBetweenArrays(custom.gdata,method = "quantile");
rownames(custom.normalized) <- rownames(custom.gdata);
colnames(custom.normalized) <- colnames(custom.gdata);

############ Knn Imputation ############ 
##library(impute);
custom.impute <- impute.knn(custom.normalized, k = 10);
custom.all <- custom.impute$data;
custom.samples <- colnames(custom.all)
custom.gnames <- rownames(custom.all)

################## Convert custom ids into gene names #################
library("hgu133a.db");
custom.ids <- custom.gnames;
symbols <- hgu133aSYMBOL[custom.ids];
convert.tab <- toTable(symbols);

# custom.data <- mdata2[custom.ids %in% convert.tab[,1],];
dup.ind <- duplicated(convert.tab[,2]);
dup.list <- as.vector(unique(convert.tab[dup.ind,2]));
# dup.index <- convert.tab[,2] %in% dup.list;

l <- dup.list; 
n <- length(l);
n.sp <- length(custom.samples)
data_dup <- matrix(0,n,n.sp); colnames(data_dup) = custom.samples; 
r <- as.vector(rep(0,n))

for (i in 1:n) {
  sym_temp <- l[i];
  id_temp <- convert.tab[convert.tab[,2]==sym_temp,1];
  data_temp <- custom.all[id_temp,];
  data_dup[i,] = colMeans(data_temp);
  r[i] = sym_temp; }
rownames(data_dup) <- r;

sym <- convert.tab[!dup.ind,];
dup_temp <- sym[,2] %in% dup.list;
uni.list <- sym[!dup_temp,];
data_uni <- custom.all[uni.list[,1],];
rownames(data_uni) <- uni.list[,2]

data_sym <- rbind(data_dup,data_uni)
symbol <- rownames (data_sym);
# write.csv(data_sym,file = 'C:/Users/HEY/Desktop/data_sym.csv',sep = " ", row.names=T, col.names = T);###

#########################################################################################
############ Select paired samples in miRNA and mRNA ##################################################
##########################################################################################

mi.etno <- mi.group!= "RT";  custom.etno <- custom.group!="RT";
mi.samp.com <- mi.samples[mi.etno]; cust.samp.com <- custom.samples[custom.etno]
paired <- mi.samp.com %in% cust.samp.com;
paired.samples <- mi.samp.com[paired];
paired.group <- mi.group[mi.etno][paired];
gp <- paired.group; samp <- paired.samples;
n.et <- sum(gp=="ET");  n.no <- sum(gp=="NO");

custom.etno <- data_sym[,paired.samples];
midata.etno <- midata.all[,paired.samples]


# output peprocessed mrna.
write.csv(rbind(mi.conditions[,colnames(custom.etno)],custom.etno), "processed_mrna.csv",
          col.names = T,
          row.names = T,
          quote = F)
# output preprocessed mirna.
write.csv(rbind(mi.conditions[,colnames(midata.etno)],midata.etno), 
          "processed_mirna.csv",
          col.names = T,
          row.names = T,
          quote = F
          )


# use limma to identify differentially expressed mirna
# Need to install miRLAB
require(miRLAB)

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
# Set the parameters for DE analysis
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

write.csv(differientially_expressed_limma, "de_mirna_mrna_et_no.csv",
          col.names = F, row.names = F)


write.csv(differientially_expressed_limma, 'de_mirna_mrna_et_no.csv',
          col.names = T, row.names = T)



### For 22 sig miRNA using SAM & FC ### 
# mi.etno.cm <- mi.data[,paired.samples]

############ 2 groups comparison SAM Analysis ###############

######################## miRNA ######################## 

# NO: 0; ET: 1 #
# mi.cl <- as.numeric(factor(gp)) - 1;
# mi.cl <- as.numeric(!mi.cl);
# sam.out <- sam(midata.etno, mi.cl, rand = 820, gene.names = rownames(midata.etno));
# 
# fdr_rate = 0.95
# sam.out;
# print(sam.out, seq(0.6, 0.9, 0.01));
# # print(sam.out, seq(0.6, 0.65, 0.0005));
# 
# glist <- summary(sam.out,0.618);
# # sam.sig <- glist@mat.sig[glist@mat.sig$q.value < 0.05 & (glist@mat.sig$R.fold > 1.5 | glist@mat.sig$R.fold < (1/1.5)),];
# sam.sig <- glist@mat.sig[glist@mat.sig$q.value < 0.05 & (glist@mat.sig$R.fold > 2 | glist@mat.sig$R.fold < (1/2)),];
# sigs <- rownames(sam.sig);
# 
# gene_list_fc = glist@mat.sig[ glist@mat.sig$R.fold > 2 | glist@mat.sig$R.fold < (1/2),];
# 
# # Find differentially expressed mRNA
# sam.out_mrna <- sam(custom.etno, mi.cl, rand = 820, gene.names = rownames(custom.etno));
# 
# fdr_rate = 0.95
# print(sam.out_mrna, seq(0.6, 0.9, 0.01));
# # print(sam.out, seq(0.6, 0.65, 0.0005));
# 
# glist <- summary(sam.out_mrna,0.618);
# # sam.sig <- glist@mat.sig[glist@mat.sig$q.value < 0.05 & (glist@mat.sig$R.fold > 1.5 | glist@mat.sig$R.fold < (1/1.5)),];
# sam.sig <- glist@mat.sig[glist@mat.sig$q.value < 0.05 & (glist@mat.sig$R.fold > 2 | glist@mat.sig$R.fold < (1/2)),];
# sigs <- rownames(sam.sig);


