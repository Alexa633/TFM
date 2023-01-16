# ========= #
# Clean Env #
# ========= #
rm(list = ls())

# ========= #
# Libraries #
# ========= #

library("devtools")
devtools::install_github("omnideconv/SimBu")

library("SimBu")
library("Matrix")


# ================= #
# Load RNA seq data #
# ================= #
cpath <- "C:/Users/moren/PycharmProjects/TFM/R/simBU"
ipath <- paste0(cpath,"/input")
opath <- paste0(cpath,"/input")
# isnrna <- paste0(ipath,"/snRNA-countmat.csv")
# isntype <- paste0(ipath,"/snRNA-celltype.csv")
isnrna <- paste0(ipath,"/snRNA-minimal_counts_major_celltypes.csv")
isntype <- paste0(ipath,"/snRNA-minimal_metadata_major_celltypes.csv")

print("Loading snRNA data..")
snrnadf <- as.matrix(read.csv(isnrna, sep=',', row.names = 1, header = TRUE))
print("    Convert to sparse Matrix object")
snrnadf <- as(t(snrnadf),"Matrix")
print("    Example rowname and colname (debug mode)")
print(row.names(snrnadf)[1])
print(colnames(snrnadf)[1])
print("    Type input matrix")
print(class(snrnadf))

sntype <- read.csv(isntype, sep=',', row.names = 1, header = TRUE)
sntype$cell_type <- sntype$major_cell_type
sntype$ID <- row.names(sntype)
print("   snRNA shape: ")
dim(snrnadf)
print("   snRNA major cell type: ")
dim(sntype)
sort(table(sntype$cell_type))

# ================================= #
# Generate objects to run SimBu     #
# ================================= #
print("Generate objects to run SimBu..")
# - All cells from snRNA Mathys
print("   All Mathys cells")
dsall <- SimBu::dataset(annotation = sntype,
                     count_matrix = snrnadf,
                     tpm_matrix = NULL,
                     filter_genes = FALSE,
                     name = "mathys",
                     variance_cutoff = 0,
                     scale_tpm = FALSE)

# ========== #
# Run SimBu  #
# ========== #
print("Run simBU using all data from Mathys with even percentage cell type")
sim_dsall_even<- SimBu::simulate_bulk(data = dsall,
                                    scenario = "even", 
                                    scaling_factor = "read_number", 
                                    ncells=1000, 
                                    nsamples = 10, 
                                    BPPARAM = BiocParallel::MulticoreParam(workers = 4),  
                                    run_parallel = TRUE)


print("Run simBU with extra weight to Microglia")
sim_dsall_Mic_weighted<- SimBu::simulate_bulk(data = dsall,
                                      scenario = "weighted", 
                                      scaling_factor = "read_number", 
                                      weighted_cell_type = "Mic",
                                      weighted_amount = 0.5,
                                      ncells=1000, 
                                      nsamples = 10, 
                                      BPPARAM = BiocParallel::MulticoreParam(workers = 4),  
                                      run_parallel = TRUE)

print("Run simBU excluding Exc")
sim_dsall_Ex_blacklist<- SimBu::simulate_bulk(data = dsall,
                                      scenario = "even", 
                                      scaling_factor = "read_number", 
                                      blacklist = c("Ex"),
                                      ncells=1000, 
                                      nsamples = 10, 
                                      BPPARAM = BiocParallel::MulticoreParam(workers = 4),  
                                      run_parallel = TRUE)


# ======================= #
# Export pseudo-bulk RNA  #
# ======================= #
print("Exporting simulated datasets")
sim_dsall_even_count  <- SummarizedExperiment::assays(sim_dsall_even$bulk)[["bulk_counts"]]
sim_dsall_even_fraction  <- sim_dsall_even$cell_fractions

sim_dsall_highMic_count <- SummarizedExperiment::assays(sim_dsall_Mic_weighted$bulk)[["bulk_counts"]]
sim_dsall_highMic_fraction <- sim_dsall_Mic_weighted$cell_fractions

sim_dsall_blacklistEx_count <- SummarizedExperiment::assays(sim_dsall_Ex_blacklist$bulk)[["bulk_counts"]]
sim_dsall_blacklistEx_fraction <- sim_dsall_Ex_blacklist$cell_fractions

ofilecount <- paste0(opath,"/pseudobulk_even-sampling_countmat.csv")
ofilefraction<- paste0(opath,"/pseudobulk_even-sampling_fraction.csv")
write.csv(t(as.matrix(sim_dsall_even_count)), ofilecount, row.names=TRUE, )
write.csv(sim_dsall_even_fraction, ofilefraction, row.names=TRUE, )

ofilecount <- paste0(opath,"/pseudobulk_highMic_countmat.csv")
ofilefraction<- paste0(opath,"/pseudobulk_highMic_fraction.csv")
write.csv(t(as.matrix(sim_dsall_highMic_count)), ofilecount, row.names=TRUE, )
write.csv(sim_dsall_highMic_fraction, ofilefraction, row.names=TRUE, )

ofilecount <- paste0(opath,"/pseudobulk_blacklistEx_countmat.csv")
ofilefraction<- paste0(opath,"/pseudobulk_blacklistEx_fraction.csv")
write.csv(t(as.matrix(sim_dsall_blacklistEx_count)), ofilecount, row.names=TRUE, )
write.csv(sim_dsall_blacklistEx_fraction, ofilefraction, row.names=TRUE, )


