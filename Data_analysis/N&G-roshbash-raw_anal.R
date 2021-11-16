# Roshbash data re-analysis.
# by Guzkiel


################################################################################--- 
## working with the single-cell data.

# libraries, variables and wd

setwd("D:/Documents/Uncoma/_Pasantía Dto FM/00-Circadian clock project/RNA-seq analysis")
load("N&G-roshbash-raw_anal.RData")

##### data loading ----

## RNA-seq data

files<-list.files("data-roshbash/raw/")
head(files)

lista<-vector(mode = 'list', length=length(files))
for (i in 1:length(files)){                             # reading tables.
lista[[i]]<-read.files(files[i], header = T, sep = ',')
}

length(lista)

for (i in 1:length(files)){                             # Checking cells in each sample.
print(ncol(lista[[i]]))
}

# Cheking the table
for(i in 1:length(lista)){
  print(colnames(lista[[i]])[1])
}


## molecules of interest list

GOI<-read.table("goi.tsv", sep = "\t", header = T)      # No lo uso mas adelante pero podría ser de utilidad.
class(GOI)
head(GOI)







#### data filtering and storing ----






## filtering by cells of interest (Pdf_only)

pdf_list<-vector(mode = "list", length = length(lista))
for (i in 1:length(lista)){
  pdf_list[[i]]<-lista[[i]][ , lista[[i]][which(lista[[i]]$X=="Pdf"), ]!=0] # df[allrows, df[rowIndex, allcolumns]!=0].
  # En este caso me con todas las células que expresan al menos 1 umi de pdf.
}
length(pdf_list)
names(pdf_list)<-files

## Filtered data storing
for (i in 1:length(lista)){
  write.table(pdf_list[[i]], file = paste("data-roshbash/filtered/pdf-only_", files[i], sep=""), sep = ",", 
              row.names = FALSE, col.names = TRUE, quote = FALSE)
}

## checking things
for (i in 1:length(pdf_list)){                             # Checking cells in each filtered sample.
  print(ncol(pdf_list[[i]]))
}








# filtering by roshbash criteria

ngenes <- numeric(length(lista))
ncells <- numeric(length(lista))
for (i in 1:length(lista)){
  ngenes[i] = nrow(lista[[i]])-1
  ncells[i] = ncol(lista[[i]])-1
}

ngenes
ncells
sum(ncells)

cont_0 <- function(vec){                                     # función que cuenta la cantidad de NO ceros en un vector.
  cont = length(which(vec>0))
  return(cont)
}



# transforming data frame into matrix
lista_ma<-vector(mode = "list", length = length(lista))
for(i in 1:length(lista)){
  lista_ma[[i]]<-as.matrix[]
}



list_filt<-vector(mode = "list", length = length(lista))

## (1) fewer than 1000 or more than 6000 detected genes (where each gene had to have at least one UMI aligned)
for (i in 1:length(lista)){
  list_filt[[i]]<-lista[[i]][ , append(apply(as.matrix(lista[[i]][ , -1]), MARGIN = 2, cont_0) > 1000, TRUE, after = 0)]  
  # En este caso me quedo con todas las células que expresan mas de 1000 genes
}
for (i in 1:length(list_filt)){
  list_filt[[i]]<-list_filt[[i]][ , append(apply(as.matrix(list_filt[[i]][ , -1]), MARGIN = 2, cont_0) < 6000, TRUE, after = 0)] 
  # En este caso me quedo con todas las células que expresan menos de 6000 genes
}

# cheking
length(list_filt)
fil_cells = numeric(length(list_filt))
for (i in 1:84){
  fil_cells[i] = ncol(list_filt[[i]])
}
fil_cells
sum(fil_cells)

# Cheking the table
for(i in 1:length(list_filt)){
  print(colnames(list_filt[[i]])[1])
}




## (2) fewer than 6000 or more than 75000 total UMI;
for (i in 1:length(lista)){
  list_filt[[i]]<-list_filt[[i]][ , append(apply(as.matrix(list_filt[[i]][ , -1]), MARGIN = 2, sum) > 6000, TRUE, after = 0)] 
  # En este caso me quedo con todas las células que tienen mas de 6000 umis
}

for (i in 1:length(list_filt)){
  list_filt[[i]]<-list_filt[[i]][ , append(apply(as.matrix(list_filt[[i]][ , -1]), MARGIN = 2, sum) < 75000, TRUE, after = 0)] 
  # En este caso me quedo con todas las células que tienen menos de 75000 umis
}

# cheking
length(list_filt)
for (i in 1:84){
  fil_cells[i] = ncol(list_filt[[i]])
}
fil_cells
sum(fil_cells)

# Cheking the table
for(i in 1:length(list_filt)){
  print(colnames(list_filt[[i]])[1])
}





## (3) gene expression entropy smaller than 5.5, where entropy was defined as -nUMI * ln(nUMI) for genes with nUMI >0,
##  where nUMI was a number of UMI in a cell.
entropy <- function(vec){                       # función que calcula la entropía
  vec  = vec[vec > 0]
  prop = vec/sum(vec)
  res  = sum(-prop*log(prop))
  return(res)
}

for (i in 1:length(list_filt)){
  list_filt[[i]]<-list_filt[[i]][ , append(apply(as.matrix(list_filt[[i]][ , -1]), MARGIN = 2, entropy) > 5.5, TRUE, after = 0)]
}


# cheking
length(list_filt)
for (i in 1:84){
  fil_cells[i] = ncol(list_filt[[i]])
}
fil_cells
sum(fil_cells)-84

# Cheking the table
for(i in 1:length(list_filt)){
  print(colnames(list_filt[[i]])[1])
}








## Filtered data storing
for (i in 1:length(list_filt)){
  write.table(list_filt[[i]], file = paste("data-roshbash/filtered-by-quality/Qfilt_", files[i], sep=""), sep = ",", 
              row.names = FALSE, col.names = TRUE, quote = FALSE)
}






save.image(file = "N&G-roshbash-raw_anal.RData")
##----------------------------------------------------------------------------##