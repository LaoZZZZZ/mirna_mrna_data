library(PMA);
library(corrplot);

## midata.etno; custom.etno; misig.data; gp=("NO","ET")

# ### sparse CCA for ET group ###
# # mirna.et <- t(midata.etno)[gp=="ET",];
# mirna.et <- t(misig.data)[gp=="ET",];
# mrna.et <- t(custom.etno)[gp=="ET",];
# out.a = CCA(mirna.et, mrna.et, typex="standard", typez="standard", K=5)
# print(out.a)
# 
# ### sparse CCA for NO group ###
# # mirna.no <- t(midata.etno)[gp=="NO",];
# mirna.no <- t(misig.data)[gp=="NO",];
# mrna.no <- t(custom.etno)[gp=="NO",];
# out.b = CCA(mirna.no, mrna.no, typex="standard", typez="standard", K=5)
# print(out.b)

setwd("/Users/lukez/Documents/mirna_mrna/original_data/")
# Load the differentially expressed mirna.
mirna_data = read.csv("differentially_expressed_mirna.csv", header = T,
                      stringsAsFactors = F)
mirna.all = apply(mirna_data[,-1],2,as.numeric)
rownames(mirna.all) = mirna_data[,1]


mrna_data = read.csv("differentially_expressed_mrna.csv", header = T,
                     stringsAsFactors = F)
mrna.all = apply(mrna_data[,-1],2,as.numeric)
rownames(mrna.all) = mrna_data[,1]
# mrna_data = read.csv("processed_mrna.csv", header = T, stringsAsFactors = F)
# mrna.all = apply(t(mrna_data[-1,-1]), 2, as.numeric)
# colnames(mrna.all) = mrna_data[-1,1]
# rownames(mrna.all) = colnames(mrna_data)[-1]



out.un = CCA(x=mirna.all, z=mrna.all, 
             typex="standard", typez="standard", 
             K=1, outcome = 'multiclass',
            xnames = colnames(mirna.all), znames = colnames(mrna.all),
            penaltyx = 0.25, penaltyz = 0.4)#, y = mi.conditions[,rownames(mirna.all)])
print(out.un, verbose=T)
# sink("C:/Users/HEY/Desktop/unsupervised.txt");
# out.un;print(out.un, verbose=T);
# sink();
gp = mi.conditions[,rownames(mirna.all)]

# Plot of canonical variables by group #
cv.et <- cbind(mirna.all[gp == "ET",] %*% out.un$u, mrna.all[gp=="ET",] %*% out.un$v);
cv.no <- cbind(mirna.all[gp == "NO",] %*% out.un$u, mrna.all[gp=="NO",] %*% out.un$v);
colnames(cv.et) <- c("miRNA","mRNA"); colnames(cv.no) <- colnames(cv.et);
r <- range(c(cv.et,cv.no));

# Boxplot #
tiff(filename = '/Users/lukez/Documents/mirna_mrna/erya_method.tif')
boxplot(cv.et, at = 0.5:1.5-0.2, boxwex=0.30, col="lightblue", ylim=r, xaxt="n", main="Canonical Variables by Groups")
boxplot(cv.no, add=T, at = 0.5:1.5+0.2, boxwex=0.30, col="lightgreen",axes=F,xaxt="n")
axis(1, at=0.5:1.5, labels=c("miRNA","mRNA"), tick=F)
legend(x="left",y="top", inset=0,c("ET","NO"), fill=c("lightblue", "lightgreen"))
dev.off();

# Scatter plot #
tiff(filename = '/Users/lukez/Documents/mirna_mrna/scatter_unsuper.tif')
plot(cv.et[,1] ~ cv.et[,2], pch=6, col="red", xlab = 'miRNA', ylab = 'mRNA')
points(cv.no[,1] ~ cv.no[,2], pch = 8, col = 'blue')
#points(cv.no,pch=8, col="blue")
legend("topleft",c("ET","NO"), col=c("red", "blue"), pch=c(6,8))
dev.off();
 
# # Permutation for p-value #
# pvalue.un <- CCA.permute(x=mirna.all, z=mrna.all, typex="standard", typez="standard", nperms=25)
# print(pvalue.un);
# plot(pvalue.un);
# out.perm = CCA(x=mirna.all, z=mrna.all, typex="standard", typez="standard", K=1, penaltyx=pvalue.un$bestpenaltyx, 
#                penaltyz=pvalue.un$bestpenaltyz, xnames = colnames(mirna.all), znames = colnames(mrna.all))


### sparse CCA for both groups, supervised for outcome prediction ###

#mirna.all <- t(misig.data);
#mrna.all <- t(custom.etno);
y.all <- as.numeric(factor(gp))

out.c = CCA(mirna.all, mrna.all, typex="standard", typez="standard", K=1, outcome="multiclass", y = y.all,
            xnames = colnames(mirna.all), znames = colnames(mrna.all), penaltyx = 0.25, penaltyz = 0.5)
print(out.c, verbose=T)
# sink("C:/Users/HEY/Desktop/supervised.txt");
# out.c;print(out.c, verbose=T);
# sink();

# Plot of canonical variables by group #
cv.et <- cbind(mirna.all[gp=="ET",] %*% out.c$u, mrna.all[gp=="ET",] %*% out.c$v);
cv.no <- cbind(mirna.all[gp=="NO",] %*% out.c$u, mrna.all[gp=="NO",] %*% out.c$v);
colnames(cv.et) <- c("miRNA","mRNA"); colnames(cv.no) <- colnames(cv.et);
r <- range(c(cv.et,cv.no));

# Boxplot #
tiff(filename = '/Users/lukez/Documents/mirna_mrna/boxplot_super.tif')
boxplot(cv.et, at = 0.5:1.5-0.2, boxwex=0.30, col="lightblue", ylim=r, xaxt="n", main="Canonical Variables by Groups")
boxplot(cv.no, add=T, at = 0.5:1.5+0.2, boxwex=0.30, col="lightgreen",axes=F,xaxt="n")
axis(1, at=0.5:1.5, labels=c("miRNA","mRNA"), tick=F)
legend(x="topleft", inset=0,c("ET","NO"), fill=c("lightblue", "lightgreen"))
dev.off();

# Scatter plot #
tiff(filename = '/Users/lukez/Documents/mirna_mrna//scatter_super.tif')
plot(cv.et, pch=6, col="red", ylim=c(r[1]-5, r[2]+5), xlim=c(3,r[2]))
points(cv.no,pch=8, col="blue")
legend("topleft",inset=0,c("ET","NO"), col=c("red", "blue"), pch=c(6,8))
dev.off();

# # Permutation for p-value #
pvalue.super <- CCA.permute(x=mirna.all, z=mrna.all, typex="standard", typez="standard", nperms=50, 
                         y = y.all, penaltyxs=0.3, penaltyzs=0.3)
print(pvalue.super, verbose = T);
sink("/Users/lukez/Documents/mirna_mrna/supervised.txt");
out.c;print(out.c, verbose=T); pvalue.super;
sink();

# Analysis with first K of supervised result #
# Data #
name1 = colnames(mirna.all); name2 = colnames(mrna.all)
mi = name1[which(out.c$u!=0)]
m = name2[which(out.c$v!=0)]
mirna = mirna.all[,mi]
mrna = mrna.all[,m]
main = cbind(mirna,mrna)

# Correlation matrix plot #
cor.all <- cor(main)
corrplot(cor.all, type="upper", tl.col="black", tl.cex=0.6, tl.srt=45, method="circle",
         addCoef.col="black", addCoefasPercent = TRUE,)

# Partial Correlation by Kith #
p <- data.frame(); 
l = dim(main)[1];
group = as.numeric(gp=="ET");
for (i in 1:(l-1)) {
  for (j in (i+1):l) {
    y1 = main[i,]; y2 = main[j,]; y3 = t(main[-c(i,j),]);
    resid1 <- y1-predict(lm(y1~y3));
    resid2 <- y2-predict(lm(y2~y3));
    p1 <- summary(lm(resid1~resid2 + resid2 * group))$coefficients[16];
    p2 <- summary(lm(resid2~resid1 + resid1 * group))$coefficients[16];
    p = rbind(p, cbind(rownames(main)[i],rownames(main)[j],mean(p1,p2)))
  }
}
names(p) <- c("node1","node2","p.value");
p.m <- as.matrix(p)
sig.gp <- p.m[as.numeric(p.m[,3])<0.05,]

# Partial Correlation by group #
main.et <- main[,gp=="ET"]
main.no <- main[,gp=="NO"]
parcorr <- data.frame(); 
pcor.value <- matrix(0,l,l); rownames(pcor.value) <- rownames(main); colnames(pcor.value) <- rownames(pcor.value)
l = dim(main)[1];
group = as.numeric(gp=="ET");
for (i in 1:(l-1)) {
  for (j in (i+1):l) {
    y1 = main.et[i,]; y2 = main.et[j,]; y3 = t(main.et[-c(i,j),]);
    resid1 <- y1-predict(lm(y1~y3));
    resid2 <- y2-predict(lm(y2~y3));
    y1.no = main.no[i,]; y2.no = main.no[j,]; y3.no = t(main.no[-c(i,j),]);
    resid1.no <- y1.no-predict(lm(y1.no~y3.no));
    resid2.no <- y2.no-predict(lm(y2.no~y3.no));
    parcorr = rbind(parcorr, cbind(rownames(main)[i],rownames(main)[j],
                             cor(resid1,resid2), cor(resid1.no,resid2.no)))
    pcor.value[i,j] <- cor(resid1,resid2); pcor.value [j,i] <- cor(resid1.no,resid2.no);
  }
}
names(parcorr) <- c("node1","node2","Parcorr_ET","Parcorr_NO");
corrplot(pcor.value) 

# scca <- list();
# n = out.c$K;
# # l = order(out.c$cors, decreasing=T);
# name1 = rownames(misig.data); name2 = rownames(custom.etno)
# # name1 = rownames(midata.etno); name2 = rownames(custom.etno)
# 
# for (i in 1:n) {
#   mi = name1[which(out.c$u[,i]!=0)]
#   m = name2[which(out.c$v[,i]!=0)]
#   n.mi <- length(mi); n.m <- length(m);
#   scca[[i]] <- list(mi, m, n.mi, n.m, out.c$cors[i])
#   names(scca[[i]]) <- c("miRNA", "mRNA", "# of miRNA", "# of mRNA", "CCA");
# }
# 
# # sink("C:/Users/HEY/Desktop/scca_all.txt");
# sink("C:/Users/HEY/Desktop/scca_sig.txt");
# scca;
# sink();

### sparse CCA for both groups, with permutation for best tuning parameter ###

# perm.out <- CCA.permute(mirna.all, mrna.all, typex="standard", typez="standard", nperms=7)
# print(perm.out)
# out.d = CCA(mirna.all, mrna.all, typex="standard", typez="standard", K=5, outcome="multiclass", y = y.all, 
#             penaltyx=perm.out$bestpenaltyx, penaltyz=perm.out$bestpenaltyz)
# print(out.d)
# 
# scca.perm <- list();
# n = out.d$K;
# name1 = rownames(misig.data); name2 = rownames(custom.etno)
# for (i in 1:n) {
#   mi = name1[which(out.d$u[,i]!=0)]
#   m = name2[which(out.d$v[,i]!=0)]
#   n.mi <- length(mi); n.m <- length(m);
#   scca.perm[[i]] <- list(mi, m, n.mi, n.m, out.d$cors[i])
#   names(scca.perm[[i]]) <- c("miRNA", "mRNA", "# of miRNA", "# of mRNA", "CCA");
# }
# 
# sink("C:/Users/HEY/Desktop/scca_sig_perm.txt");
# scca.perm;
# sink();
