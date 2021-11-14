# Celniker&graveley data re-analysis.
# by ezeq Perez

################################################################################--- 
## working with the develomental stage data (graveley et al 2011).





#### libraries and wd ----

library("dplyr")
library("edgeR")
setwd("D:/Documents/Uncoma/_Pasantía Dto FM/00-Circadian clock project/RNA-seq analysis")
load("N&G-graveley-anal_rawcts.Rdata")



#### data loading ----

tt <- read.table("data2/gene_rpkm_report_fb_2021_02.tsv", quote = "", sep = "\t", header = T)
GOI <- read.table("goi.tsv", sep = "\t", header = T)


#### data filtering and management ----

tt$counts = as.numeric(tt$X478) * as.numeric(tt$X299) * 1000 # calculo de conteos de reads crudos. RPKM -> crudos.
raw.counts = tt[ ,c(3, 7, 13)]
names(raw.counts) = c("genes", "exp", "counts")
expi = c("mE_mRNA_A_MateM_1d_head", "mE_mRNA_A_MateM_4d_head", "mE_mRNA_A_MateM_20d_head", 
         "mE_mRNA_A_MateF_1d_head", "mE_mRNA_A_MateF_4d_head", "mE_mRNA_A_MateF_20d_head",
         "mE_mRNA_A_VirF_1d_head", "mE_mRNA_A_VirF_4d_head", "mE_mRNA_A_VirF_20d_head")
lista = list(mode=list, length=length(expi))
for (i in 1:length(expi)){
  lista[[i]] = raw.counts[which(raw.counts$exp == expi[i]), ]
}
names(lista) = expi
length(lista)

rawcts = full_join(lista[[1]], lista[[2]], by = "genes")
for (i in 3:length(expi)){
  rawcts = full_join(rawcts, lista[[i]], by = "genes")
}
rawcts = rawcts[ , seq(-18, -2, by = 2)]
colnames(rawcts) = c("genes", expi)
rawcts[is.na(rawcts)] = 0 # Na generados los transforma en 0.

#cheking all is right. 
head(rawcts)
length(rawcts$genes)

# tabla de reads crudos !!!!!
write.table(rawcts, file = "data-graveley/graveley-raw_counts.tsv", quote = F, row.names = F, col.names = T, na = "0", sep = "\t")





##------------=------------------=

#### EdgeR pipeline ----

      ## MALES VS MATED FEMALES

# bulding an object DGElist
samples<-c("Male", "Male", "MFemale", "MFemale")      # first, it is needed a vector pointing out the treatment which each sample belongs to.
mcounts<-as.matrix(rawcts[ , c(2,3,5,6)])             # second, it is needed a matrix with the raw counts, not a dataframe!!
rownames(mcounts)<-rawcts[,1]
head(mcounts)
class(mcounts)
lcounts<-DGEList(mcounts, group=samples) # construction of the object.
head(lcounts)

##### Filtering ---

dim(lcounts$counts)                        # dimension of my counts
filtb_counts <- rowSums(cpm(lcounts)>1)>= 3 # vector with boolean values. Filter all transcripts which have less than 1 cpm in at least 3 samples
filt_counts <- lcounts[filtb_counts,]      # selecting of the rows that are "TRUE" according with the boolean vector.
head(filt_counts)

##### Normalization ---
norm_counts <- calcNormFactors(filt_counts)
head(norm_counts)

##### Dispersion estimate ---

# Estimate common dispersion
norm_counts <- estimateCommonDisp(norm_counts)

# Estimate trended and tagwise dispersion:
norm_counts <- estimateTrendedDisp(norm_counts) # tiro error.
norm_counts <- estimateTagwiseDisp(norm_counts)
head(norm_counts) # cheking what I got.

##### DE analysis ---

de <- exactTest(norm_counts, dispersion="tagwise", rejection.region="doubletail")
head(de) # cheking

# FDR correction
fdr_de <- topTags(de, n = nrow(de))
head(fdr_de) # cheking
write.table(fdr_de, "data2/graveley-DE_analysis-MaHa-egdeR.tsv", sep = "\t", row.names = T, col.names = T, quote = F)

# subset of significant fold change 
fctable = fdr_de$table # whole table
fdr05_de <- subset(fdr_de$table, FDR < 0.05)
fdr10_de <- subset(fdr_de$table, FDR < 0.10)
head(fdr05_de)
head(fdr10_de)
write.table(fdr05_de, "data2/graveley-sig05_DE-MaHa-egdeR.tsv", sep = "\t", row.names = T, col.names = T, quote = F)
write.table(fdr10_de, "data2/graveley-sig10_DE-MaHa-egdeR.tsv", sep = "\t", row.names = T, col.names = T, quote = F)

##  Check Pvalue distribution
hist(de$table$PValue, breaks=20)

##  Lets make a MA plot or Smear plot to see DE genes vs rest
plotSmear(de, main="plotSmear")
plotSmear(de, de.tags=rownames(fdr05_de), main="Machos vs Hembras apareadas")
plotSmear(de, de.tags=rownames(fdr10_de), main="Machos vs Hembras apareadas")

# plots propios
goi <- GOI$symbol      # lista de etiquetas goi
goi10 = filter(fdr10_de, rownames(fdr10_de)%in%goi)
goi10
goi05 = filter(fdr05_de, rownames(fdr05_de)%in%goi)
goi05
fcgoi = filter(fctable, rownames(fctable)%in%goi)
head(fcgoi)

attach(fctable)
op = par(cex.axis = 2, cex.lab = 2, cex.main = 2, mar = c(5, 5, 5, 1))
plot(logCPM, logFC, cex=1.9, col='grey85', main = "Machos vs Hembras apareadas", 
     ylab = 'LogFC', xlab = 'Average logCPM')
abline(h=0)
for (i in 1:length(fcgoi$logFC)){
  points(fcgoi$logCPM[i], fcgoi$logFC[i], col="black", pch=20, cex = 2)
}
for (i in 1:length(fdr10_de$logFC)){
  points(fdr10_de$logCPM[i], fdr10_de$logFC[i], col="pink", pch=20, cex = 2)
}
for (i in 1:length(fdr05_de$logFC)){
  points(fdr05_de$logCPM[i], fdr05_de$logFC[i], col="red", pch=20, cex = 2)
}
for (i in 1:length(goi10$logFC)){
  text(goi10$logCPM[i], goi10$logFC[i], labels = rownames(goi10)[i], col = "black", pos = 1, cex = 1.8)
}
legend(x = "bottomright", legend = c('GOIs', 'FC p-value < .10', ' FC p-value < .05'), pt.cex = 2,
       cex = 1.2, bty = 'n', pch = 20, y.intersp = .5, x.intersp = .5, col = c('black', 'pink', 'red'))
par(op)
detach(fctable)




      ## MALES VS VIRGIN FEMALES


# bulding an object DGElist
samples<-c("Male", "Male", "VFemale", "VFemale")      # first, it is needed a vector pointing out the treatment which each sample belongs to.
mcounts<-as.matrix(rawcts[ , c(2,3,8,9)])             # second, it is needed a matrix with the raw counts, not a dataframe!!
rownames(mcounts)<-rawcts[,1]
head(mcounts)
class(mcounts)
lcounts<-DGEList(mcounts, group=samples) # construction of the object.
head(lcounts)

##### Filtering ---

dim(lcounts$counts)                        # dimension of my counts
filtb_counts <- rowSums(cpm(lcounts)>1)>=3 # vector with boolean values. Filter all transcripts which have less than 1 cpm in at least 3 samples
filt_counts <- lcounts[filtb_counts,]      # selecting of the rows that are "TRUE" according with the boolean vector.
head(filt_counts)

##### Normalization ---
norm_counts <- calcNormFactors(filt_counts)
head(norm_counts)

##### Dispersion estimate ---

# Estimate common dispersion
norm_counts <- estimateCommonDisp(norm_counts)

# Estimate trended and tagwise dispersion:
norm_counts <- estimateTrendedDisp(norm_counts) # tiro error.
norm_counts <- estimateTagwiseDisp(norm_counts)
head(norm_counts) # cheking what I got.

##### DE analysis ---

de <- exactTest(norm_counts, dispersion="tagwise", rejection.region="doubletail")
head(de) # cheking

# FDR correction
fdr_de <- topTags(de, n = nrow(de))
head(fdr_de) # cheking
write.table(fdr_de, "data2/graveley-DE_analysis-MaHv-egdeR.tsv", sep = "\t", row.names = T, col.names = T, quote = F)

# subset of significant fold change 
fctable = fdr_de$table # whole table
fdr05_de <- subset(fdr_de$table, FDR < 0.05)
fdr10_de <- subset(fdr_de$table, FDR < 0.10)
head(fdr05_de)
head(fdr10_de)
write.table(fdr05_de, "data2/graveley-sig05_DE-MaHv-egdeR.tsv", sep = "\t", row.names = T, col.names = T, quote = F)
write.table(fdr10_de, "data2/graveley-sig10_DE-MaHv-egdeR.tsv", sep = "\t", row.names = T, col.names = T, quote = F)

##  Check Pvalue distribution
hist(de$table$PValue, breaks=20)

##  Lets make a MA plot or Smear plot to see DE genes vs rest
plotSmear(de, main="plotSmear")
plotSmear(de, de.tags=rownames(fdr05_de), main="Machos vs Hembras virgenes")
plotSmear(de, de.tags=rownames(fdr10_de), main="Machos vs Hembras virgenes")

# plots propios
goi <- GOI$symbol      # lista de etiquetas goi
goi10 = filter(fdr10_de, rownames(fdr10_de)%in%goi)
goi10
goi05 = filter(fdr05_de, rownames(fdr05_de)%in%goi)
goi05
fcgoi = filter(fctable, rownames(fctable)%in%goi)
head(fcgoi)

attach(fctable)
op = par(cex.axis = 2, cex.lab = 2, cex.main = 2, mar = c(5, 5, 5, 1))
plot(logCPM, logFC, cex=1.9, col='grey85', main = "Machos vs Hembras virgenes", 
     ylab = 'LogFC', xlab = 'Average logCPM')
abline(h=0)
for (i in 1:length(fcgoi$logFC)){
  points(fcgoi$logCPM[i], fcgoi$logFC[i], col="black", pch=20, cex = 2)
}
for (i in 1:length(fdr10_de$logFC)){
  points(fdr10_de$logCPM[i], fdr10_de$logFC[i], col="pink", pch=20, cex = 2)
}
for (i in 1:length(fdr05_de$logFC)){
  points(fdr05_de$logCPM[i], fdr05_de$logFC[i], col="red", pch=20, cex = 2)
}
for (i in 1:length(goi10$logFC)){
  text(goi10$logCPM[i], goi10$logFC[i], labels = rownames(goi10)[i], col = "black", pos = 1, cex = 1.8)
}
legend(x = "bottomright", legend = c('GOIs', 'FC p-value < .10', ' FC p-value < .05'), pt.cex = 2,
       cex = 1.2, bty = 'n', pch = 20, y.intersp = .5, x.intersp = .5, col = c('black', 'pink', 'red'))
par(op)
detach(fctable)



## MATED FEMALES VS VIRGIN FEMALES


# bulding an object DGElist
samples<-c("MFemale", "MFemale", "VFemale", "VFemale")      # first, it is needed a vector pointing out the treatment which each sample belongs to.
mcounts<-as.matrix(rawcts[ , c(5,6,8,9)])             # second, it is needed a matrix with the raw counts, not a dataframe!!
rownames(mcounts)<-rawcts[,1]
head(mcounts)
class(mcounts)
lcounts<-DGEList(mcounts, group=samples) # construction of the object.
head(lcounts)

##### Filtering ---

dim(lcounts$counts)                        # dimension of my counts
filtb_counts <- rowSums(cpm(lcounts)>1)>=3 # vector with boolean values. Filter all transcripts which have less than 1 cpm in at least 3 samples
filt_counts <- lcounts[filtb_counts,]      # selecting of the rows that are "TRUE" according with the boolean vector.
head(filt_counts)

##### Normalization ---
norm_counts <- calcNormFactors(filt_counts)
head(norm_counts)

##### Dispersion estimate ---

# Estimate common dispersion
norm_counts <- estimateCommonDisp(norm_counts)

# Estimate trended and tagwise dispersion:
norm_counts <- estimateTrendedDisp(norm_counts) # tiro error.
norm_counts <- estimateTagwiseDisp(norm_counts)
head(norm_counts) # cheking what I got.

##### DE analysis ---

de <- exactTest(norm_counts, dispersion="tagwise", rejection.region="doubletail")
head(de) # cheking

# FDR correction
fdr_de <- topTags(de, n = nrow(de))
head(fdr_de) # cheking
write.table(fdr_de, "data2/graveley-DE_analysis-HaHv-egdeR.tsv", sep = "\t", row.names = T, col.names = T, quote = F)

# subset of significant fold change 
fctable = fdr_de$table # whole table
fdr05_de <- subset(fdr_de$table, FDR < 0.05)
fdr10_de <- subset(fdr_de$table, FDR < 0.10)
head(fdr05_de)
head(fdr10_de)
write.table(fdr05_de, "data2/graveley-sig05_DE-HaHv-egdeR.tsv", sep = "\t", row.names = T, col.names = T, quote = F)
write.table(fdr10_de, "data2/graveley-sig10_DE-HaHv-egdeR.tsv", sep = "\t", row.names = T, col.names = T, quote = F)

##  Check Pvalue distribution
hist(de$table$PValue, breaks=20)

##  Lets make a MA plot or Smear plot to see DE genes vs rest
plotSmear(de, main="plotSmear")
plotSmear(de, de.tags=rownames(fdr05_de), main="Hembras apareadeas vs Hembras virgenes")
plotSmear(de, de.tags=rownames(fdr10_de), main="Hembras apareadeas vs Hembras virgenes")

# plots propios
goi <- GOI$symbol      # lista de etiquetas goi
goi10 = filter(fdr10_de, rownames(fdr10_de)%in%goi)
goi10
goi05 = filter(fdr05_de, rownames(fdr05_de)%in%goi)
goi05
fcgoi = filter(fctable, rownames(fctable)%in%goi)
head(fcgoi)

attach(fctable)
op = par(cex.axis = 2, cex.lab = 2, cex.main = 2, mar = c(5, 5, 5, 1))
plot(logCPM, logFC, cex=1.9, col='grey85', main = "Hembras apareadeas vs Hembras virgenes", 
     ylab = 'LogFC', xlab = 'Average logCPM')
abline(h=0)
for (i in 1:length(fcgoi$logFC)){
  points(fcgoi$logCPM[i], fcgoi$logFC[i], col="black", pch=20, cex = 2)
}
for (i in 1:length(fdr10_de$logFC)){
  points(fdr10_de$logCPM[i], fdr10_de$logFC[i], col="pink", pch=20, cex = 2)
}
for (i in 1:length(fdr05_de$logFC)){
  points(fdr05_de$logCPM[i], fdr05_de$logFC[i], col="red", pch=20, cex = 2)
}
for (i in 1:length(goi10$logFC)){
  text(goi10$logCPM[i], goi10$logFC[i], labels = rownames(goi10)[i], col = "black", pos = 1, cex = 1.8)
}
legend(x = "bottomright", legend = c('GOIs', 'FC p-value < .10', ' FC p-value < .05'), pt.cex = 2,
       cex = 1.2, bty = 'n', pch = 20, y.intersp = .5, x.intersp = .5, col = c('black', 'pink', 'red'))
par(op)
detach(fctable)







#------------------------------------------------------------------------------------------------------------------
save.image(file = "N&G-graveley-anal_rawcts.Rdata")
