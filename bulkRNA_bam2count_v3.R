rm(list = ls())

# Libraries
# ==========
library(Rsubread)

# ======================
# Description:
#   1) Load a .txt file with a list of bam_files to read
#   2) Obtain counts from each .bam file individually.
#   3) Save output counts in a .txt file
# ========================

# ==========
# Defaults
# ==========
print("Reading files..")
listfile <- '/media/neuroimagen2/MU-data-cluster/Users/vmontalb/student/2022-Alexandra-UOC/data/list.list'
rawdatapath <-'/media/neuroimagen2/MU-data-cluster/Users/vmontalb/student/2022-Alexandra-UOC/data/'

ingenome <- '/media/neuroimagen2/MU-data-cluster/Users/vmontalb/student/2022-Alexandra-UOC/utils/gencode.v14.annotation.gtf'

# ==========
# Load Files
# ==========
inlist <- read.delim(listfile, header=FALSE, stringsAsFactors=FALSE)

# ==========
# Obtain gene counts per file
# ==========
allsubj <- inlist[[1]]

print("Computing gene read-counts...")
for (csubj in allsubj)
{
  print(paste0("Working with file: ",csubj))
  
  # Set path
  dir.create(paste0(rawdatapath,csubj,"/feaureCounts"), showWarnings = FALSE)
  csubjin <- paste0(rawdatapath,csubj,"/",csubj,'.bam')      # 
  csubjoutcount <- paste0(rawdatapath,csubj,"/feaureCounts/",csubj,'_countReads.txt')      # 
  csubjoutsumm <- paste0(rawdatapath,csubj,"/feaureCounts/",csubj,'_countReads.summary')      # 
  
  # Run cmd count
  ccount <- featureCounts(csubjin, isGTFAnnotationFile=TRUE, annot.ext=ingenome, GTF.attrType='gene_name', 
                          isPairedEnd=TRUE, strandSpecific=2, countMultiMappingReads=TRUE)
  
  # Save data
  write.table(ccount$counts, csubjoutcount, sep=" ", row.names=TRUE)
  write.table(ccount$stat, csubjoutsumm, sep=" ", row.names=TRUE)
  
}