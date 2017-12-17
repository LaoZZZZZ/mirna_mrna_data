# Data analysis
setwd("/Users/lukez/Documents/mirna_mrna/original_data/")
# Load the data

de_mi_mr = read.csv('de_mirna_mrna_et_no.csv', header = T, stringsAsFactors = F)

samples = de_mi_mr[,1]
data = de_mi_mr[,-1]

mirnas = 1:61
mrnas = 62:80

de_mirna = data[,mirnas]
de_mrna = data[,mrnas]


mirna_file = 'de_mirna_et_no.csv'
mrna_file = 'de_mrna_et_no.csv'
write.csv(de_mirna, mirna_file, col.names = T, quote = F)
write.csv(de_mrna, mrna_file, col.names = T, quote = F)


# Pearsong correlation analysis

library(psych)
pearson_cor = corr.test(data, method = 'pearson')
sig_p = 0.01

correlation_pearson = pearson_cor$r
correlation_pearson[pearson_cor$p > sig_p] = 0



# Sparse Canonical component analysis
library(PMA)
scca <- function(diff_expressed_mirna_file,
                 diff_expressed_mrna_file,
                 px,pz,k = 1) {
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
               penaltyx = px, penaltyz = pz)
  #print(out.un, verbose=T)
  return(out.un)
}

penaltyx = 0.3
penaltyz = 0.7
scca.result = scca(mirna_file,mrna_file, penaltyx, penaltyz)
print(scca.result, verbose = T)

# Sparse partial correlation analysis
library(space)

# set the parameters
p = ncol(data)
n = nrow(data)

data = apply(data,2,as.numeric)
alpha=1
l1=1/sqrt(n)*qnorm(1-alpha/(2*p^2))
iter=3

#### Joint method with no weight
par_result=space.joint(data, lam1=l1*n*1.56, lam2=0, iter=iter)
par_corr = par_result$ParCor

# Sparse Bayesian network
