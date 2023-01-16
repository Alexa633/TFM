# ========= #
# Clean Env #
# ========= #
rm(list = ls())
#############################################
####Libraries and BayesPrism installation####
library("devtools")
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
cpath <- getwd()
ipath <- paste0(cpath,'/input/ALL_FTD_S2_FINAL_CUT_CUT_FINAL.csv')
opath_exp <- paste0(cpath,'/output_exp/')
opath_theta <- paste0(cpath,'/output_theta/')
#
#
############################################
############################################
#
#Load csv matrix count
#
data <- read.csv(ipath, sep= '\t', header=T)
#traspose, subjects in row-wise
x1 <- t(data)
str(x1)
bk.data <- as.data.frame(x1)
str(bk.data)
#gene as columns
names(bk.data) <- bk.data[1,]
bk.data <- bk.data[-1,]
#matrix-count as numeric
bk.data <- data.frame(lapply(bk.data, function(a) as.numeric(as.character(a))))
str(bk.data)

#Add ID column
ID <- names(data)[-1]
row.names(bk.data) <- ID
bk.data <- cbind(ID = rownames(bk.data), bk.data)
rownames(bk.data) <- NULL
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
row.names(theta.df) <- ID
#ID as column not index

theta.df <- cbind(ID = rownames(theta.df), theta.df)
rownames(theta.df) <- NULL

write.csv(theta.df, paste0(opath_theta,'theta_ELA_DFT.csv'), row.names=T)
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
  row.names(types.df) <- ID
  types.df <- cbind(ID = rownames(types.df), types.df)
  rownames(types.df) <- NULL
  

  write.csv(types.df, paste0(opath_exp,types,'_ELA_DFT.csv'), row.names=T)
  print(paste0('Saved: ',types))
}
######################################################
######################################################



