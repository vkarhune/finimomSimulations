# rm(list=ls())

args <- commandArgs(trailingOnly = T)

reps <- as.integer(args[1])
seedstart <- as.integer(args[2])

clump <- FALSE

.libPaths("/projappl/minmanni/project_rpackages_421/")

library(bigsnpr)
library(corpcor)
library(finimom) # remotes::install_github("vkarhune/finimom")

set.seed(123) # the first seed is for genotype imputation

bedfile <- "data/synthetic_filtered.bed"

rds <- snp_readBed(bedfile, backingfile = tempfile())

obj <- snp_attach(rds)

# str(obj)

geno <- obj$genotypes

get_geno_mat <- function(geno){
 geno[seq_len(dim(geno)[1]), seq_len(dim(geno)[2])]
}

# check MAF
mafs <- apply(get_geno_mat(geno), 2, function(x) sum(x)/(2*length(x)))

geno <- geno[,mafs >= 0.01]

eafs <- mafs[mafs >= 0.01]

variants <- obj$map$marker.ID[mafs >= 0.01]

###

### reference data
bedfile_ref <- "data/synthetic_ref_filtered.bed"
rds_ref <- snp_readBed(bedfile_ref, backingfile = tempfile())
obj_ref <- snp_attach(rds_ref)
# str(obj)
geno_ref <- obj_ref$genotypes
geno_ref <- geno_ref[,obj_ref[["map"]]$marker.ID %in% variants]
# dim(geno)
# geno[1:5,1:5]
# impute missing
# geno_imputed <- snp_fastImputeSimple(geno, method = "random")
mafs_ref <- apply(get_geno_mat(geno_ref), 2, function(x) sum(x)/(2*length(x)))
#eafs <- mafs[mafs >= 0.01]
#variants <- obj$map$marker.ID[mafs >= 0.01]



### 

n <- nrow(geno)
p <- ncol(geno)

windows <- c(500, 1000, 2000, 5000, 10000)

sapply(seq_along(windows), function(i){
# i <- 1
 windowsize <- windows[i]
 genostart <- ifelse(i == 1, 1, cumsum(windows)[i - 1] + 1)
 genostop <- cumsum(windows)[i]
 
 geno_window <- geno[,genostart:genostop]

 # save LD file
 cat(sprintf("Calculating LD matrix for %i variants\n", windowsize))
 prc <- proc.time()
 LDmat <- cor(geno_window)
 cat(sprintf("LD matrix done in %.2f seconds\n", (proc.time() - prc)[[3]]))
 LD <- list(LDmat, eafs[genostart:genostop], variants[genostart:genostop])

 LDfile <- paste0("data/LDmatsynthetic_", paste0(c(windowsize, n, clump), collapse = "_"), ".Rds")

 cat(sprintf("Number of observations: %i\n", n))
 cat(sprintf("Number of variants: %i\n", windowsize))

 saveRDS(LD, file = LDfile)
 
 ### and the reference data
 geno_ref_window <- geno_ref[,genostart:genostop]
 cat(sprintf("Calculating reference LD matrix for %i variants\n", windowsize))
 prc <- proc.time()
 LDmat_ref <- cor(geno_ref_window)
 cat(sprintf("LD matrix done in %.2f seconds\n", (proc.time() - prc)[[3]]))
 LDref <- list(LDmat_ref, mafs_ref[genostart:genostop], variants[genostart:genostop])
 LDreffile <- paste0("data/refLDmatsynthetic_", paste0(c(windowsize, n, clump), collapse = "_"), ".Rds")
 saveRDS(LDref, file = LDreffile)
 ###

 
 # simulate phenotype
 allnumcausals <- c(1, 2, 5)
 allR2 <- c(0.001, 0.015)

 d_pars <- expand.grid(allnumcausals, allR2)

 d_pars <- do.call("rbind", replicate(reps, d_pars, simplify = FALSE))
 names(d_pars) <- c("numcausals", "R2")

 d_pars <- d_pars[with(d_pars, order(numcausals, R2)),]

 d_pars$seed <- seq_len(nrow(d_pars))

 d_pars <- d_pars[,c("seed", "numcausals", "R2")]

 d_parslist <- split(d_pars, d_pars$seed)



invisible(lapply(d_parslist, function(x){
# x <- d_parslist[[177]]

seed <- x[["seed"]] + seedstart + - 1
numcausals <- x[["numcausals"]]
R2 <- x[["R2"]]

set.seed(seed)

causals <- select_causals(numcausals = numcausals, LDmat = LDmat,
 mincorr = 0, maxcorr = 0.95, minanycorr = 0)

pheno <- simulate_phenotype_data(X = geno_window, N = n, p = ncol(geno_window),
 seed = seed, causals = causals, R2 = R2, minpower = 0.8, alpha = 0.05)

summarystats <- generate_sumstats(phenotype = pheno[[1]], X = geno_window)

out <- list(summarystats, causals, pheno[[2]][causals], LDmat[causals,causals])

outname <- paste0("data/simdatasynthetic", paste0(c(seed, paste0(windowsize, "variants"), n, clump, numcausals, R2), collapse = "_"), ".Rds")

saveRDS(out, outname)

}))

cat(sprintf("%i\n", windowsize))

})

Sys.Date()

sessionInfo()
