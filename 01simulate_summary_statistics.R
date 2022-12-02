# rm(list=ls())

args <- commandArgs(trailingOnly = T)

generegion <- args[1] # remember to use here the unpruned genotype data
reps <- as.integer(args[2])
clump <- as.logical(args[3])
seedstart <- as.integer(args[4])
simulate_geno <- as.logical(args[5])
if(length(args) == 6){ simN <- as.integer(args[6]) }

if(0){
 generegion <- "UMPS"
 reps <- 100
 clump <- FALSE
 simulate_geno <- FALSE
}

.libPaths("/home/vkarhune/rpackages/")

library(bigsnpr)
library(finimom) # remotes::install_github("vkarhune/finimom")

set.seed(123) # the first seed is for genotype imputation

bedfile <- paste0("data/", generegion, "_filtered.bed")

rds <- snp_readBed(bedfile, backingfile = tempfile())

obj <- snp_attach(rds)

# str(obj)

geno <- obj$genotypes

# dim(geno)
# geno[1:5,1:5]

# impute missing
geno_imputed <- snp_fastImputeSimple(geno, method = "random")

get_geno_mat <- function(geno){
 geno[seq_len(dim(geno)[1]), seq_len(dim(geno)[2])]
}

# check MAF
mafs <- apply(get_geno_mat(geno_imputed), 2, function(x) sum(x)/(2*length(x)))

# summary(mafs)

geno_imputed <- geno_imputed[,mafs >= 0.01]

eafs <- mafs[mafs >= 0.01]

variants <- obj$map$marker.ID[mafs >= 0.01]

# not used
if(simulate_geno){

}


# not used
if(clump){

}

# 

n <- nrow(geno_imputed)
p <- ncol(geno_imputed)

LDmat <- cor(geno_imputed)

LD <- list(LDmat, eafs, variants)

LDfile <- paste0("data/LDmat_", paste0(c(generegion, n, clump), collapse = "_"), ".Rds")
if(simulate_geno) { LDfile <- paste0("data/LDmat_", paste0(c(generegion, paste0(n, "sim"), clump), collapse = "_"), ".Rds") }

saveRDS(LD, file = LDfile)

cat(sprintf("Number of observations: %i\n", n))
cat(sprintf("Number of variants: %i\n", p))

allnumcausals <- c(1, 2, 5)
allR2 <- c(0.015, 0.03)

d_pars <- expand.grid(allnumcausals, allR2)

d_pars <- do.call("rbind", replicate(reps, d_pars, simplify = FALSE))
names(d_pars) <- c("numcausals", "R2")

d_pars <- d_pars[with(d_pars, order(numcausals, R2)),]

d_pars$seed <- seq_len(nrow(d_pars))

d_pars <- d_pars[,c("seed", "numcausals", "R2")]

d_parslist <- split(d_pars, d_pars$seed)



invisible(lapply(d_parslist, function(x){
# x <- d_parslist[[177]]

seed <- x[["seed"]] + seedstart - 1
numcausals <- x[["numcausals"]]
R2 <- x[["R2"]]

set.seed(seed)

causals <- select_causals(numcausals = numcausals, LDmat = LDmat,
 mincorr = 0, maxcorr = 0.95, minanycorr = 0)

pheno <- simulate_phenotype_data(X = geno_imputed, N = n, p = p,
 seed = seed, causals = causals, R2 = R2, minpower = 0.8, alpha = 0.05)

summarystats <- generate_sumstats(phenotype = pheno[[1]], X = geno_imputed)

out <- list(summarystats, causals, pheno[[2]][causals], LDmat[causals,causals])

outname <- paste0("data/simdata", paste0(c(seed, generegion, n, clump, numcausals, R2), collapse = "_"), ".Rds")

saveRDS(out, outname)

}))

Sys.Date()

sessionInfo()
