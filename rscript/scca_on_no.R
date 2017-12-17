# this script is designed to find regulatory network using sCCA
setwd("/Users/lukez/Documents/mirna_mrna/original_data/")

if (TRUE) {
mrna_data = read.csv("processed_mrna.csv", header = T, stringsAsFactors = F)
mrna.all = apply(t(mrna_data[-1,-1]), 2, as.numeric)
colnames(mrna.all) = mrna_data[-1,1]
rownames(mrna.all) = colnames(mrna_data)[-1]

mirna_data = read.csv("processed_mirna.csv", header = T, stringsAsFactors = F)

mirna.all = apply(t(mirna_data[-1,-1]), 2, as.numeric)
colnames(mirna.all) = mirna_data[-1,1]
rownames(mirna.all) = colnames(mirna_data)[-1]
} else {
  
  mrna_data = read.csv("differentially_expressed_mrna.csv", header = T, stringsAsFactors = F)
  mrna.all = apply(mrna_data[,-1],2,as.numeric)
  rownames(mrna.all) = mrna_data[,1]
  
  mirna_data = read.csv("differentially_expressed_mirna.csv", header = T, stringsAsFactors = F)
  mirna.all = apply(mirna_data[,-1],2,as.numeric)
  rownames(mirna.all) = mirna_data[,1]
}

c = mi.conditions[, rownames(mirna.all)]


no_mi =mirna.all[c == 'NO',]
no_mr = mrna.all[c == 'NO',]


et_mi = mirna.all[c != 'NO',]
et_mr = mrna.all[c != 'NO',]
library(GMD)
library(CCA)
scca <- function(mirna,
                 mrna,
                 k = 1,
                 ...) {
  
  out.un = CCA(x=mirna, z=mrna, typex="standard", typez="standard", K=k, 
               xnames = colnames(mirna.all), znames = colnames(mrna.all),
               ...)
  #print(out.un, verbose=T)
  return(out.un)
}

num_comp = 3

u_pen = 0.1
v_pen = 0.1
components = scca(et_mi, et_mr,
                  k = num_comp,
                  penaltyx = u_pen,
                  penaltyz = v_pen)
print(components$cors)

print(components, verbose = T)


u = components$u
v = components$v

mirna_cca_feature = mirna.all %*% u
colnames(mirna_cca_feature) = paste("cca_component", seq(1,num_comp), sep = '_')
mrna_cca_feature = mrna.all %*% v
colnames(mrna_cca_feature) = paste("cca_component", seq(1,num_comp), sep = '_')


# use t test to find the different network between ET and NO

t_test <- function(x, y) {
  test <- t.test(x,y, paired=TRUE)
  out <- data.frame(stat = test$statistic,
                    df   = test$parameter,
                    pval = test$p.value,
                    conl = test$conf.int[1],
                    conh = test$conf.int[2]
  )
  return(out)
}

sapply(seq(num_comp), function(x) t_test(mirna_cca_feature[c == 'NO', x],
                                         )
k = 1
x_lim = c(min(mirna_cca_feature[,k]),max(mirna_cca_feature[,k]))
y_lim = c(min(mrna_cca_feature[,k]),max(mrna_cca_feature[,k]))
plot(mirna_cca_feature[c == 'NO',k], mrna_cca_feature[c =='NO',k], pch = 19, col = 'blue',
     xlab = 'miRNA',
     ylab = 'mRNA',
     xlim = x_lim,
     ylim = y_lim)

points(mirna_cca_feature[c != 'NO',k], mrna_cca_feature[c !='NO',k],
       pch = 19, col = 'red')

mirna_cca_feature_et = 
