##########################################
######BayesPrism for SimBu data###########
##########################################
#
#
# ========= #
# Clean Env #
# ========= #
rm(list = ls())
#############################################
####Libraries and BayesPrism installation####
#library("devtools")
#install_github("Danko-Lab/BayesPrism/BayesPrism")
#
#
#
#Bioconductor
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")
#
#
#############################################
library(BayesPrism)
library(tidyverse)
#############################################
#############################################
#in and out paths
cpath <- "C:/Users/moren/PycharmProjects/TFM/R/simBU"
ipath <- paste0(cpath,'/output_SimBu/')
opath <- paste0(cpath,'/results_BP_SimBu/')
#
#
############################################
############################################
#
#Load csv matrix count
#
data <- read.csv(paste0(ipath,'pseudobulk_even-sampling_countmat.csv'))
colnames(data)[1] <- "ID"
bk.data <- data[,-1]
#bk.data.prueba <- bk.data[1:2,]
#
#Load sc data and celltype info (same lenght)
#
sc.dat <-read.csv("C:/Users/moren/PycharmProjects/TFM/R/BayesPrism/snRNA-minimal_counts_major_celltypes.csv")
cell.types.labels <- read.csv("C:/Users/moren/PycharmProjects/TFM/R/BayesPrism/snRNA-minimal_metadata_major_celltypes.csv")
#
#
#Tables by cell type (main and subtypes)
#
sort(table(cell.types.labels$major_cell_type))
sort(table(cell.types.labels$subtype_cell_type))
#
#Define cell type level
#
major.cell.types <- cell.types.labels$major_cell_type
#
#
#
#Correlation plot between cell types
#
plot.cor.phi (input=sc.dat, 
              input.labels=major.cell.types, 
              title="cell type correlation",
              #specify pdf.prefix if need to output to pdf
              #pdf.prefix="gbm.cor.ct",
              cexRow=0.5, cexCol=0.5,
)
#
#clean up chr X/Y and mitochondrial genes
sc.dat.filtered <- cleanup.genes (input=sc.dat,
                                  input.type="count.matrix",
                                  species="hs", 
                                  gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
                                  exp.cells=5)
#
#plot gene filtered
plot.bulk.vs.sc (sc.input = sc.dat.filtered,
                 bulk.input = bk.data
                 #pdf.prefix="gbm.bk.vs.sc" specify pdf.prefix if need to output to pdf
)
#
#filtered by protein coding
sc.dat.filtered.pc <-  select.gene.type (sc.dat.filtered,
                                         gene.type = "protein_coding")
#
#
###########################################
###########################################
#Setting BayesPrism parameters
#
myPrism <- new.prism(
  reference=sc.dat.filtered.pc, 
  mixture=bk.data,
  input.type="count.matrix", 
  cell.type.labels = major.cell.types, 
  cell.state.labels = major.cell.types,
  key=NULL,
  outlier.cut=0.01,
  outlier.fraction=0.1,
)
#
#
#Run BayesPrism
#
bp.res <- run.prism(prism = myPrism, n.cores=2)
#
###########################################
###########################################
#Get cell type fraction
theta <- get.fraction (bp=bp.res,
                       which.theta="final",
                       state.or.type="type")
#
theta.df <- as.data.frame(theta)
ID <- data[,1]
row.names(theta.df) <- ID
#ID as column not index

theta.df <- cbind(ID = rownames(theta.df), theta.df)
rownames(theta.df) <- NULL

write.csv(theta.df, 'theta.csv', row.names=T)
print(paste0('Theta file saved'))
#
#
#
#Get expression csv
#
cell.types <-levels(as.factor(major.cell.types))
#
for (types in cell.types){

  types.AD <- get.exp (bp=bp.res,
                      state.or.type='type',
                      cell.name=types) 
  types.df <- as.data.frame(types.AD)
  

  write.csv(types.df, paste0(opath,types,'.csv'), row.names=T)
  print(paste0(types,' file saved'))
}
######################################################
######################################################



