# Celniker&graveley data re-analysis.
# by ezeq Perez

################################################################################--- 
## working with the develomental stage data (graveley et al 2011).





#### libraries and wd ----

library("reshape")
library("ggplot2")
library("car")
library("lsmeans")
#library("contrast")

options(contrasts=c("contr.sum","contr.poly"))  ## needed for correct conditional tests

setwd("D:/Documents/Uncoma/_Pasantía Dto FM/00-Circadian clock project/RNA-seq analysis")
load("graveley_var.Rdata")





#### Preparing the data ----

# reading tables
t <- read.table("data2/gene_rpkm_matrix_fb_2021_01.tsv", quote = "", sep = "\t", header = T)
GOI <- read.table("goi.tsv", sep = "\t", header = T)

# Looking for genes of interest in the data set.
goi <- GOI$symbol # lista de etiquetas.

gidx <- numeric(length = length(goi)) # goi indexs en la base de datos.
for (i in 1:length(goi)){  # file number of the GOIs.
  gidx[i]<-which(t$gene_symbol==goi[i])
}

# extrayendo las columnas de interes, transformacion log10.
ma4d = log(t$mE_mRNA_A_MateM_4d_head_.FBlc0000216, 10) # macho apareado 4 dias.
ha4d = log(t$mE_mRNA_A_MateF_4d_head_.FBlc0000213.,10) # hembra apareada 4 dias.
hv4d = log(t$mE_mRNA_A_VirF_4d_head_.FBlc0000211., 10) # hembra virgen 4 dias.
ma1d = log(t$mE_mRNA_A_MateM_1d_head_.FBlc0000209.,10) # macho apareado 1 dia.
ha1d = log(t$mE_mRNA_A_MateF_1d_head_.FBlc0000207.,10) # hembra apareada 1 dia.
hv1d = log(t$mE_mRNA_A_VirF_1d_head_.FBlc0000210., 10) # hembra virgen 1 dia.
ma20d= log(t$mE_mRNA_A_MateM_20d_head_.FBlc0000214.,10)# macho apareado 20 dias.
ha20d= log(t$mE_mRNA_A_MateF_20d_head_.FBlc0000212.,10)# hembra apareada 20 dias.
hv20d= log(t$mE_mRNA_A_VirF_20d_head_.FBlc0000231., 10)# hembra virgen 20 dias.

gt = data.frame(10^(ma1d), 10^(ma4d), 10^(ha1d), 10^(ha4d), t$gene_symbol) # data frame con las variables extraidas destransformadas.
names(gt) = c("ma1d", "ma4d", "ha1d", "ha4d", "genes")
head(gt)

gt_filt = gt[gidx, ] # filtrado de dataframe por genes de interes
names(gt_filt) = c("male", "male", "female", "female", "genes")
gt_filt

genet = melt(gt_filt, id.vars = 'genes')
names(genet) = c("genes", "sexo", "expresion")
genet





#### exploring the data ----

# plots exploratorios de DE

# Machos
op = par(mfrow = c(2, 2), cex.lab = 1.7)
plot(ma1d, ma4d) # machos 1 vs 4 dias.
abline(0, 1, col='red')
for (i in 1:length(goi)){
  points(ma1d[gidx[i]], ma4d[gidx[i]], col=50+i, pch=20, cex = 1.2)
}
for (i in 1:length(goi)){
  text(ma1d[gidx[i]], ma4d[gidx[i]], labels = t$gene_symbol[gidx[i]], col = 50+i, pos = 4, cex = 1.5)
}
plot(ma1d, ma20d) # machos 1 vs 20 dias.
abline(0, 1, col='red')
for (i in 1:length(goi)){
  points(ma1d[gidx[i]], ma20d[gidx[i]], col=50+i, pch=20, cex = 1.2)
}
for (i in 1:length(goi)){
  text(ma1d[gidx[i]], ma20d[gidx[i]], labels = t$gene_symbol[gidx[i]], col = 50+i, pos = 4, cex = 1.5)
}
plot(ma4d, ma20d) # machos 4 vs 20 dias.
abline(0, 1, col='red')
for (i in 1:length(goi)){
  points(ma4d[gidx[i]], ma20d[gidx[i]], col=50+i, pch=20, cex = 1.2)
}
for (i in 1:length(goi)){
  text(ma4d[gidx[i]], ma20d[gidx[i]], labels = t$gene_symbol[gidx[i]], col = 50+i, pos = 4, cex = 1.5)
}
par(op)

# hembras apareadas.
op = par(mfrow = c(2, 2), cex.lab = 1.7)
plot(ha1d, ha4d) # hembras ap. 1 vs 4 dias.
abline(0, 1, col='red')
for (i in 1:length(goi)){
  points(ha1d[gidx[i]], ha4d[gidx[i]], col=50+i, pch=20, cex = 1.2)
}
for (i in 1:length(goi)){
  text(ha1d[gidx[i]], ha4d[gidx[i]], labels = t$gene_symbol[gidx[i]], col = 50+i, pos = 4, cex = 1.5)
}
plot(ha1d, ha20d) # hembras ap. 1 vs 20 dias.
abline(0, 1, col='red')
for (i in 1:length(goi)){
  points(ha1d[gidx[i]], ha20d[gidx[i]], col=50+i, pch=20, cex = 1.2)
}
for (i in 1:length(goi)){
  text(ha1d[gidx[i]], ha20d[gidx[i]], labels = t$gene_symbol[gidx[i]], col = 50+i, pos = 4, cex = 1.5)
}
plot(ha4d, ha20d) # hembras ap. 4 vs 20 dias.
abline(0, 1, col='red')
for (i in 1:length(goi)){
  points(ha4d[gidx[i]], ha20d[gidx[i]], col=50+i, pch=20, cex = 1.2)
}
for (i in 1:length(goi)){
  text(ha4d[gidx[i]], ha20d[gidx[i]], labels = t$gene_symbol[gidx[i]], col = 50+i, pos = 4, cex = 1.5)
}
par(op)

# hembras virgenes.
op = par(mfrow = c(2, 2), cex.lab = 1.7)
plot(hv1d, hv4d) # hembras virg. 1 vs 4 dias
abline(0, 1, col='red')
for (i in 1:length(goi)){
  points(hv1d[gidx[i]], hv4d[gidx[i]], col=50+i, pch=20, cex = 1.2)
}
for (i in 1:length(goi)){
  text(hv1d[gidx[i]], hv4d[gidx[i]], labels = t$gene_symbol[gidx[i]], col = 50+i, pos = 4, cex = 1.5)
}
plot(hv1d, hv20d) # hembras virg. 1 vs 20 dias
abline(0, 1, col='red')
for (i in 1:length(goi)){
  points(hv1d[gidx[i]], hv20d[gidx[i]], col=50+i, pch=20, cex = 1.2)
}
for (i in 1:length(goi)){
  text(hv1d[gidx[i]], hv20d[gidx[i]], labels = t$gene_symbol[gidx[i]], col = 50+i, pos = 4, cex = 1.5)
}
plot(hv4d, hv20d) # hembras virg. 4 vs 20 dias
abline(0, 1, col='red')
for (i in 1:length(goi)){
  points(hv4d[gidx[i]], hv20d[gidx[i]], col=50+i, pch=20, cex = 1.2)
}
for (i in 1:length(goi)){
  text(hv4d[gidx[i]], hv20d[gidx[i]], labels = t$gene_symbol[gidx[i]], col = 50+i, pos = 4, cex = 1.5)
}
par(op)





#### Graficos definitivos ----


                                              ## log-log plots

# macho ap vs hembra ap
plot_1 = plot(ma4d, ha4d, xlab="Macho apareado", ylab="Hembra apareada", cex=0.5, col='gray91', main = "Expression in Log10(RPKM)") 
abline(0, 1)
# Coloring the GOIs 
for (i in 1:length(goi)){
  points(ma4d[gidx[i]], ha4d[gidx[i]], col=50+i, pch=20, cex = 1.2)
}
for (i in 1:length(goi)){
  text(ma4d[gidx[i]], ha4d[gidx[i]], labels = t$gene_symbol[gidx[i]], col = 50+i, pos = 4, cex = 0.5)
}
# legend(x = "bottomright", legend = goi, fill = seq(51, 51+length(goi),1), 
      # title = "GOI",
      # cex = 0.7)


# hembra virgen vs hembra apareada
plot_2 = plot(hv4d, ha4d, xlab = 'hembra virgen', ylab = 'hembra apareada', cex = 0.5, col = 'gray91') 
abline(0, 1)
# Coloring the GOIs 
for (i in 1:length(goi)){
  points(hv4d[gidx[i]],ha4d[gidx[i]], col=50+i, pch=20, cex = 1.2)
}
for (i in 1:length(goi)){
  text(hv4d[gidx[i]], ha4d[gidx[i]], labels = t$gene_symbol[gidx[i]], col = 50+i, pos = 4, cex = 0.5)
}
# legend(x = "bottomright", legend = goi, fill = seq(51, 51+length(goi),1), 
       #title = "GOI",
       #cex = 0.7)


# macho ap. vs hembra virgen
plot_3 = plot(ma4d, hv4d, xlab = 'macho apareado', ylab = 'hembra virgen', cex = 0.5, col = 'gray91') # Macho ap vs hembra virg.
abline(0, 1)
# Coloring the GOIs 
for (i in 1:length(goi)){
  points(ma4d[gidx[i]],hv4d[gidx[i]], col=50+i, pch=20, cex = 1.2)
}
for (i in 1:length(goi)){
  text(ma4d[gidx[i]], hv4d[gidx[i]], labels = t$gene_symbol[gidx[i]], col = 50+i, pos = 4, cex = 0.5)
}
#legend(x = "bottomright", legend = goi, fill = seq(51, 51+length(goi),1), 
       #title = "GOI",
       #cex = 0.7)


                                        ## histogramas 
par(mfrow=c(2,1))

# machos ap. vs hembras ap.
h_ma4d = hist(ma4d, br=20, main = 'machos apareados')
h_ha4d = hist(ha4d, br=20, main = 'hembras apareadas')
# lineas de densidad hembras vs machos
lines(h_ma4d$mids, h_ma4d$counts, col='red')
lines(h_ha4d$mids, h_ha4d$counts, col='blue')

# hembras virg. vs hembras ap.
h_hv4d = hist(hv4d,br=20, main = 'hembras virgenes')
h_ha4d = hist(ha4d,br=20, main = 'hembras apareadas')
# lineas de densidad hembras virg. vs hembras ap.
lines(h_hv4d$mids, h_hv4d$counts, col='red')
lines(h_ha4d$mids, h_ha4d$counts, col='blue')

# machos ap. vs hembras virg.
h_ma4d = hist(ma4d, br=20, main = 'machos apareados')
h_hv4d = hist(hv4d,br=20, main = 'hembras virgenes')
# lineas de densidad hembras virg. vs hembras ap.
lines(h_ma4d$mids, h_ma4d$counts, col='red')
lines(h_hv4d$mids, h_hv4d$counts, col='blue')

par(mfrow=c(1,1))


                                      ## MAPlots

# machos vs hembras
avrg = rowMeans(cbind(ma4d,ha4d)) # promedio expresion de las dos condiciones
l2r = log2(ma4d/ha4d) # diferencia expresion entre condiciones

plot(avrg,l2r,cex=0.5,col='grey91', main = 'MAplot: Ma/Ha', ylab = 'Log2 fold change', xlab = 'Mean of normalizated counts')
abline(h=0)
# Coloring the GOIs
for (i in 1:length(goi)){
  points(avrg[gidx[i]], l2r[gidx[i]], col=100+i, pch=20)
}
for (i in 1:length(goi)){
  text(avrg[gidx[i]], l2r[gidx[i]], labels = t$gene_symbol[gidx[i]], col = 100+i, pos = 4, cex = 0.45)
}
 #legend(x = "bottomright", legend = goi, fill = seq(51, 51+length(goi),1), 
       #title = "GOI",
      #cex = 0.7)


# Hembras virg. vs hembras ap.
avrg_2 = rowMeans(cbind(hv4d,ha4d)) # promedio expresion de las dos condiciones
l2r_2 = log2(hv4d/ha4d) # diferencia expresion entre condiciones

plot(avrg_2, l2r_2, cex=0.5, col='grey91', main = 'MAplot: Hv/Ha', ylab = 'Log2 fold change', xlab = 'Mean of normalizated counts')
abline(h=0)
# Coloring the GOIs
for (i in 1:length(goi)){
  points(avrg_2[gidx[i]], l2r_2[gidx[i]], col=100+i, pch=20)
}
for (i in 1:length(goi)){
  text(avrg_2[gidx[i]], l2r_2[gidx[i]], labels = t$gene_symbol[gidx[i]], col = 100+i, pos = 4, cex = 0.45)
}
#legend(x = "bottomright", legend = goi, fill = seq(51, 51+length(goi),1), 
       #title = "GOI",
       #cex = 0.7)


# machos vs hembras virg.
avrg_3 = rowMeans(cbind(ma4d,hv4d)) # promedio expresion de las dos condiciones
l2r_3 = log2(ma4d/hv4d) # diferencia expresion entre condiciones

plot(avrg_3, l2r_3, cex=0.5, col='grey91', main = 'MAplot: Ma/Hv', ylab = 'Log2 fold change', xlab = 'Mean of normalizated counts')
abline(h=0)
# Coloring the GOIs
for (i in 1:length(goi)){
  points(avrg_3[gidx[i]], l2r_3[gidx[i]], col=100+i, pch=20)
}
for (i in 1:length(goi)){
  text(avrg_3[gidx[i]], l2r_3[gidx[i]], labels = t$gene_symbol[gidx[i]], col = 100+i, pos = 4, cex = 0.45)
}
#legend(x = "bottomright", legend = goi, fill = seq(51, 51+length(goi),1), 
       #title = "GOI",
       #cex = 0.7)



                                      ## Barplots

# gen de interes
gen_idx = gidx[which(goi=='Pdf')] # cambiar segun gen buscado

# sin transformacion
gen = as.numeric(t[gen_idx, 5:128]) # Expresion del gen para todas las condiciones y tejidos.
names(gen) = colnames(t)[5:128]
gen_filt = gen[gen>0] # filtro condiciones y tejidos que no tengan conteos.

# escala logaritmica
lgen = log(as.numeric(t[gen_idx, 5:128])+1, 10) # goi numero 1
names(lgen) = colnames(t)[5:128]
lgen_filt = lgen[lgen>0]

# barplot con expresion del gen escogido en cada tejido y estadio, sin transformacion ni filtrado
par(mar=c(12,5,5,2))
barplot(gen, names=names(gen), las=2, cex.names=0.3)

# barplot con expresion del gen escogido en cada tejido y estadio, sin transformacion y filtrados los genes con 0 expresion
barplot(gen_filt, names=names(gen_filt), las=2, cex.names=0.4)

# barplot con expresion del gen escogido en cada tejido y estadio, escala RPKM logaritmica sin filtrado.
barplot(lgen, names=names(lgen), las=2, cex.names=0.4)

# barplot con expresion del gen escogido en cada tejido y estadio, escala RPKM logaritmica y filtrados los genes con 0 expresion
barplot(lgen_filt, names=names(lgen_filt), las=2, cex.names=0.5)
par(mar=c(2.5,2.5,1,1))



                                      ## comparacion entre transcriptomas en gral

Ma4d_goi = ma4d[gidx] # solo me quedo con los genes de interes.
Ha4d_goi = ha4d[gidx]
Hv4d_goi = hv4d[gidx]
tabla_log = cbind(Ma4d_goi, Ha4d_goi, Hv4d_goi)
tabla_log = melt(tabla_log)
tabla_log = tabla_log[ , 2:3]
class(tabla_log)
colnames(tabla_log) = c('mosca', 'Log10(RPKM)')
levels(tabla_log$mosca)
head(tabla_log) # checking things.

ma4d_nl = t$mE_mRNA_A_MateM_4d_head_.FBlc0000216 # columnnas de interes sin transformar (no logaritmicas)
ha4d_nl = t$mE_mRNA_A_MateF_4d_head_.FBlc0000213.
hv4d_nl = t$mE_mRNA_A_VirF_4d_head_.FBlc0000211.

Ma4d_goi2 = ma4d_nl[gidx[-12]] # quito el valor atipico de Yp1
Ha4d_goi2 = ha4d_nl[gidx[-12]]
Hv4d_goi2 = hv4d_nl[gidx[-12]]
tabla = cbind(Ma4d_goi2, ha4d_nl, Hv4d_goi2)
tabla = melt(tabla)
tabla = tabla[ , 2:3]
class(tabla)
colnames(tabla) = c('mosca', 'Log(RPKM)')
levels(tabla$mosca)
head(tabla)


# plots
ggplot(tabla_log, aes(x=mosca, y=`Log10(RPKM)`, color=mosca, shape=mosca))+
  geom_point(size=2.5, position=position_jitterdodge(0.3)) +
  stat_summary(aes(x=mosca, y=`Log10(RPKM)`), fun = 'mean', geom = 'crossbar', color = 'black', position = position_dodge(0.8), 
               width = 0.45) +
  stat_summary(aes(x=mosca, y=`Log10(RPKM)`), geom = 'errorbar', color = 'black', position = position_dodge(0.75), width = 0.3) +
  ylab("Expresion Log(RPKM)") + 
  xlab("Mosca") + 
  theme_classic()


ggplot(tabla, aes(x=mosca, y=`Log(RPKM)`, color=mosca, shape=mosca))+
  geom_point(size=2.5, position=position_jitterdodge(0.3)) +
  stat_summary(aes(x=mosca, y=`Log(RPKM)`), fun = 'mean', geom = 'crossbar', color = 'black', position = position_dodge(0.8), 
               width = 0.45) +
  stat_summary(aes(x=mosca, y=`Log(RPKM)`), geom = 'errorbar', color = 'black', position = position_dodge(0.75), width = 0.3) +
  ylab("Expresion (RPKM)") + 
  xlab("Mosca") + 
  theme_classic()





#### Analisis de expresion diferencial ----

# analisis propio con los genes de interes.
mod1 = glm(expresion ~ sexo * genes, data = genet, family = poisson) # modelo lineal generalizado con distribucion de poisson
summary(mod1)
1-pchisq(mod1$deviance, mod1$df.residual)
Anova(mod1, type=3)
meds = lsmeans(mod1, ~sexo*genes) # medias del modelo.
meds
#contrast(meds, method = tukey.emmc(goi), interaction = T, adjust = "bonferroni")
#contrast(mod1, a = list(sexo = levels(genet$sexo), genes = levels(as.factor(genet$genes))))

save.image(file = "graveley_var.Rdata")
##----------------------------------------------------------------------------##